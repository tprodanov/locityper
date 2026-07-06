#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

readonly SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]:-$0}")"
readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]:-$0}")"

function help_message {
  cat <<HELP
Usage: $SCRIPT_NAME \\
    (-a assemblies.agc | -g assemblies_dir) \\
    (-t targets.fa | -b targets.bed -r reference.fa) \\
    -o directory [args]
then:  $SCRIPT_NAME -b targets.bed --combine directory -o final_directory

│ Maps target sequences to assembly genomes and extracts corresponding subregions.
│ Multiple instances of this script can be run in parallel on the same output directory,
│ as it creates "*.lock" files for haplotypes in progress and ".ok" files for the finished haplotypes.

Input/output arguments:
    -a, --agc         FILE  Input AGC file.
    -g, --genomes     DIR   Directory with various genome assemblies (.fa[.gz]).
                            Mutually exclusive with -a/--agc.
    -n, --names       FILE  Optional: replace genome names (first column) with another name (second column).
                            In case of -g/--genomes, first column should match file basename without extension.
                            Use "skip" in the second column to skip this samples.
    -t, --targets     FILE  FASTA file with target sequences. Name lines should not contain spaces.
    -b, --targets-bed FILE  Instead of -t/--targets, provide BED file with reference locus coordinates.
    -r, --reference   FILE  Provide reference genome (must be indexed) to extract reference haplotypes.
                            Only relevant together with -b/--target-bed.
    -o, --output      DIR   Output directory.
        --combine     DIR   After all genomes were processed (potentially in multiple parallel instances),
                            combine fasta files from DIR into a new final directory.
                            Can be specified multiple times to provide multiple input directories.
                            If -b is provided, copies it with new fifth column to the output directory.

Filter arguments:
    -d, --distance    INT   Merge PAF entries if distance is smaller than INT [${distance}].
    -l, --min-len     NUM   Consider regions longer than NUM * target length [${min_len}].
    -s, --min-simil   NUM   Consider regions with approximate similarity over NUM [${min_simil}].
    -N, --count       INT   Extract top N sequences (by length) for each haplotype and target [${count}].

Minimap2 arguments:
    -x, --preset      STR   Minimap2 preset [${preset}].
    -@, --threads     INT   Minimap2 threads [${threads}].
    -S, --secondary   INT   Number of secondary alignments [${secondary}].
    -p, --score-ratio NUM   Secondary-to-primary score ratio [${score_ratio}].
    --  ARGUMENTS           If necessary, provide additional minimap2 arguments.

Other arguments:
    -h, --help              Print this help and exit.
HELP
}

# [TODO] Separate provide minimap arguments.

function setup_colors {
    readonly RED="\e[31m"
    readonly ENDCOLOR="\e[0m"
}

function msg {
    echo -e "$*" >&2
}

function err {
    msg "${RED}[ERROR]${ENDCOLOR} $*"
}

function panic {
    err "$1"
    exit "${2-1}" # Return 1 by default.
}

function parse_params {
    combine=()

    min_len=0.5
    min_simil=0.3
    distance=5000
    count=1

    preset=asm20
    threads=3
    secondary=50
    score_ratio=0.5

    long1=help,agc:,genomes:,names:,targets:,targets-bed:,reference:,output,combine:
    long2=distance:,min-len:,min-simil:,count:
    long3=preset:,threads:,secondary:,score-ratio:
    ARGS="$(getopt -o a:g:n:t:b:r:o:d:l:s:N:x:@:p:h --long "${long1},${long2},${long3}" \
        --name "$SCRIPT_NAME" -- "$@")"
    eval set -- "$ARGS"
    while :; do
        case "$1" in
            -a | --agc )
                agc_file="$2"; shift 2 ;;
            -g | --genomes )
                genomes_dir="$2"; shift 2 ;;
            -t | --targets )
                targets_fa="$2"; shift 2 ;;
            -b | --targets-bed )
                targets_bed="$2"; shift 2 ;;
            -r | --reference )
                reference="$2"; shift 2 ;;
            -n | --names )
                names_file="$2"; shift 2 ;;
            -o | --output )
                output="$2"; shift 2 ;;
            --combine )
                combine+=( "$2" ); shift 2 ;;
            -d | --distance )
                distance="$2"; shift 2 ;;
            -l | --min-len )
                min_len="$2"; shift 2 ;;
            -s | --min-simil )
                min_simil="$2"; shift 2 ;;
            -N | --count )
                count="$2"; shift 2 ;;
            -x | --preset )
                preset="$2"; shift 2 ;;
            -@ | --threads )
                threads="$2"; shift 2 ;;
            -S | --secondary )
                secondary="$2"; shift 2 ;;
            -p | --score-ratio )
                score_ratio="$2"; shift 2 ;;
            -h | --help)
                help_message; exit 0;
                ;;
            -- ) shift; break ;;
            * )  panic "Unexpected argument $1" ;;
        esac
    done

    minimap2_args=( -c -x "$preset" -t "$threads" -N "$secondary" -p "$score_ratio" )
    if [[ ${#@} -ne 0 ]]; then
        minimap2_args+=( "$@" )
    fi

    [[ ! -z "${output-}" ]]   || panic "Missing required parameter -o/--output"
    [[ -z "${targets_fa-}" ]] && have_targets_fa=n || have_targets_fa=y
    [[ -z "${targets_bed-}" ]] && have_targets_bed=n || have_targets_bed=y

    if [[ ${combine[@]} ]]; then
        [[ "$have_targets_bed" = y ]] || panic "--combine requires -b"
    else
        [[ $have_targets_fa = y || $have_targets_bed = y ]] || panic "Either -t or -b is required"
        if [[ $have_targets_fa = y && $have_targets_bed = y ]]; then
            msg "Both -t and -b are provided, will use FASTA from -t"
            have_targets_bed=n
        fi

        [[ $have_targets_fa = n || -f "${targets_fa-}" ]] || panic "Targets file ${targets_fa-} not found"
        [[ $have_targets_bed = n || ( -f "${targets_bed-}" && -f "${reference-}" ) ]] \
            || panic "Targets BED file ${targets_bed-} or reference FASTA ${reference-} not found or not provided"

        [[ -z "${agc_file-}" ]] && have_agc=n || have_agc=y
        [[ -z "${genomes_dir-}" ]] && have_genomes=n || have_genomes=y
        [[ $have_agc != $have_genomes ]] || panic "Require either -a or -g, but not both"
    fi
}

function combine_files {
    [[ ${combine[@]} ]] || return 0

    local lock_file="${output}/lock"
    ( set -C; 2>/dev/null > "$lock_file" ) || \
        panic "Output directory is locked. --combine is not supposed to work in parallel"
    trap 'rm -f "${lock_file}"; exit 1' INT TERM ERR

    rm -f "${output}/"*.fa.gz
    find "${combine[@]}" -name "*.fa.gz" | \
        awk -F/ 'BEGIN {OFS=FS} {
            sample = $(NF - 1)
            basename = $NF
            if (seen[(sample basename)]++) {
                if (!seen[sample]++) {
                    printf("Ignoring second occurance of sample \"%s\"\n", sample) > "/dev/stderr";
                }
            } else {
                print;
            }
        }' | while read filename; do
            cat "$filename" >> "${output}/$(basename "$filename")"
        done

    cat "${targets_bed}" | while read chrom start end target extra; do
        if [[ -f "${output}/${target}.fa.gz" ]]; then
            # Output BED file (fifth column = FASTA with locus haplotypes).
            echo -e "${chrom}\t${start}\t${end}\t${target}\t${target}.fa.gz"
        fi
    done > "${output}/targets.bed"

    rm -f "${lock_file}"
    trap - INT TERM ERR
    exit 0
}

function load_names {
    [[ ! -z "${names_file-}" ]] || return 0
    while read name upd_name; do
        names["$name"]="$upd_name"
    done < "$names_file"
}

function prepare_targets {
    [[ $have_targets_bed = y ]] || return 0

    targets_fa="${output}/__targets__.fa"
    local targets_tmp="${targets_fa}.tmp"
    local targets_lock="${targets_fa}.lock"

    msg "Extracting target sequences"
    # Try to take lock 10 times.
    for i in {1..10}; do
        if ! ( set -C; 2>/dev/null > "$targets_lock" ); then
            # If lock exists, wait 20 seconds until it releases.
            timeout 20 sh -c "while [[ -f \"${targets_lock}\" ]]; do sleep 0.1; done"
            [[ $? -eq 0 ]] || panic "Lock file ${targets_lock} exists for too long, perhaps relevant process was killed"
            [[ ! -f "${targets_fa}" ]] || return
        else
            trap 'rm -f "${targets_lock}"; exit 1' INT TERM ERR
            cat "${targets_bed}" | while read chrom start end name extra; do
                samtools faidx "${reference}" "${chrom}:$((start+1))-${end}" | \
                    seqtk seq -U -l 120 | \
                    sed "1c>$name"
            done > "${targets_tmp}"
            mv "${targets_tmp}" "${targets_fa}"
            rm -f "${targets_lock}"
            trap - INT TERM ERR
            return
        fi
    done
    panic "Most likely, reference target sequences could not be extracted in other instances of this program"
}

function process_genome {
    local arg="$1"
    local genome_name
    [[ $have_agc = y ]] && genome_name="$arg" || genome_name="$(basename "${arg%.fa*}")"

    local short_name
    # :- if unset or empty, use $genome_name
    short_name="${names["$genome_name"]:-"$genome_name"}"
    [[ "$short_name" != skip ]] || return

    local prefix="${output}/${short_name}"
    local ok_file="${prefix}.ok"
    local lock_file="${prefix}.lock"
    [[ ! -f "$ok_file" ]] || return
    ( set -C; 2>/dev/null > "$lock_file" ) || return
    trap 'rm -f "${lock_file}"; exit 1' INT TERM ERR

    # ===== START ======
    msg "Processing $short_name"
    mkdir -p "$prefix"

    local genome_fasta
    if [[ $have_agc = y ]]; then
        msg "    Extracting genome sequence"
        genome_fasta="${prefix}.fa"
        agc getset "$agc_file" "$genome_name" > "${genome_fasta}"
        samtools faidx "${genome_fasta}"
    else
        genome_fasta="$arg"
    fi

    local paf_filename="${prefix}.paf.gz"
    if [[ ! -f "$paf_filename" ]]; then
        msg "    Mapping targets to assembly"
        minimap2 "${minimap2_args[@]}" "$genome_fasta" "$targets_fa" 2> /dev/null | \
            gzip > "${paf_filename}.tmp" \
            && mv "${paf_filename}"{.tmp,}
    fi

    msg "    Extracting subsequences"
    "$SCRIPT_DIR/inner/merge_hits.py" "$paf_filename" -g "${short_name}" -o "${prefix}.bed.gz" -d "$distance" \
        2> "${prefix}.warnings.csv"
    # At this point, $prefix.bed.gz will have columns
    # chrom, start, end, strand (+/-), target name, length fraction, similarity.

    rm -f "${prefix}/"*.fa{,.gz}
    # Take $count first entries,
    #     convert into: region, strand faidx argument ("-i" or ""), suffix ("" or "-<INDEX>").
    # Then, fetch regions from the current assembly.
    zcat "${prefix}.bed.gz" | \
        awk -F$'\t' -v name="$short_name" -v min_len="$min_len" -v min_simil="$min_simil" -v max_count="$count" \
            -v warnings_file="${prefix}.warnings.csv" \
            'BEGIN{OFS=";"} $6 >= min_len && $7 >= min_simil {
                target = $5;
                ind = ++target_count[target];
                if (ind <= max_count) {
                    region = $1 ":" ($2+1) "-" $3;
                    strand_arg = $4 == "+" ? "" : "-i";
                    suffix = ind == 1 ? "" : ("-" ind);
                    print target, region, strand_arg, suffix
                }
            } END {
                for (target in target_count) {
                    if (target_count[target] > max_count) {
                        printf("%s\t%s\tFound too many hits (%d)\n", target, name, target_count[target]) >> warnings_file
                    }
                }
            }' | \
        while IFS=";" read target region strand_arg suffix; do
            samtools faidx "$genome_fasta" "$region" $strand_arg | \
                seqtk seq -U -l 120 | \
                sed "1c>${short_name}${suffix} ${region}" >> "${prefix}/${target}.fa"
        done
    # This way we don't have to worry about zero matches.
    find "${prefix}" -name "*.fa" -exec gzip {} ';'

    [[ $have_agc = n ]] || rm "${genome_fasta}"{,.fai}
    # ===== END ======

    touch "${ok_file}"
    rm -f "${lock_file}"
    trap - INT TERM ERR
}

setup_colors
parse_params "$@"
combine_files
mkdir -p "$output"

declare -A names
load_names

prepare_targets
# zcat -f opens plain files as well. sed -n does not print by default.
readarray -t target_names < <(zcat -f "$targets_fa" | sed -n 's/>//p' | sort -u)

if [[ $have_agc = y ]]; then
    agc listset "$agc_file" | while read genome; do
        process_genome "$genome"
    done
else
    for filename in "${genomes_dir}"/*.fa{,sta}{,.gz}; do
        process_genome "$filename"
    done
fi
