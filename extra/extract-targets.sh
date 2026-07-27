#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

readonly SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]:-$0}")"
readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]:-$0}")"

function help_message {
  cat <<HELP
Usage: $SCRIPT_NAME \\
    -i assemblies [-i assemblies2 ...] [-n aliases.txt] \\
    -c targets.bed -r reference.fa [-t targets.fa] \\
    -o directory [args]

│ Maps target sequences to assembly genomes and extracts corresponding subregions.
│ Multiple instances of this script can be run in parallel on the same output directory,
│ as it creates "*.lock" files for haplotypes in progress and ".ok" files for the finished haplotypes.

Input/output arguments:
    -i, --input   FILE|DIR  Input assemblies: either AGC file, FA[.gz] file, or directory with FA[.gz] files.
                            Can be repeated multiple times. Repeated assemblies will be skipped.
                            All extensions must be in lowercase.
    -n, --names       FILE  Optional: replace genome names (first column) with another name (second column).
                            In case of the FASTA files, first column should match basename without extension.
                            Use "skip" in the second column to skip this assembly.
    -c, --coordinates FILE  BED file with reference locus coordinates.
    -r, --reference   FILE  Reference genome.
    -t, --targets     FILE  Optional: instead of fetching target sequences from the reference,
                            directly provide them in this FASTA file. Names must match -c/--coordinates.
    -o, --output      DIR   Output directory. Will contain subdirectories
                            "by_assembly" with local haplotypes extracted from each assembly,
                            and "panels" with combined reference panels for each locus.

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
    -T, --timeout     INT   Number of minutes the script will wait for other instances to finish
                            in order to start combining reference panels [${timeout}].
                            If timeout is exceeded, this instance will stop.
    -h, --help              Print this help and exit.
HELP
}

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

function parse_args {
    input=()

    min_len=0.5
    min_simil=0.3
    distance=5000
    count=1

    preset=asm20
    threads=1
    secondary=50
    score_ratio=0.5
    timeout=10

    long1=help,input:,names:,coordinates:,reference:,targets:,output:
    long2=distance:,min-len:,min-simil:,count:
    long3=preset:,threads:,secondary:,score-ratio:,timeout:
    ARGS="$(getopt -o i:n:c:r:t:o:d:l:s:N:x:@:p:T:h --long "${long1},${long2},${long3}" \
        --name "$SCRIPT_NAME" -- "$@")"
    eval set -- "$ARGS"
    while :; do
        case "$1" in
            -i | --input )
                input+=( "$2" ); shift 2 ;;
            -n | --names )
                names_file="$2"; shift 2 ;;
            -c | --coord | --coordinates )
                targets_bed="$2"; shift 2 ;;
            -r | --reference )
                reference="$2"; shift 2 ;;
            -t | --targets )
                targets_fa="$2"; shift 2 ;;
            -o | --output )
                output="$2"; shift 2 ;;
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
            -T | --timeout )
                timeout="$2"; shift 2 ;;
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

    [[ ${input[@]:+${input[@]}} ]] || panic "Missing -i/--input"
    [[ -n "${targets_bed-}" ]] || panic "Missing required parameter -c/--coordinates"
    [[ -n "${reference-}" ]] || panic "Missing required parameter -r/--reference"
    [[ -n "${output-}" ]] || panic "Missing required parameter -o/--output"
    [[ -n "${targets_fa-}" ]] && have_targets_fa=y || have_targets_fa=n

    [[ $have_targets_fa = n || -f "${targets_fa-}" ]] || panic "Targets file ${targets_fa-} not found"
    [[ -f "${targets_bed-}" && -f "${reference-}" ]] \
        || panic "Targets BED file ${targets_bed-} or reference FASTA ${reference-} not found"
}

function load_names {
    [[ -n "${names_file-}" ]] || return 0
    while read name upd_name; do
        names["$name"]="$upd_name"
    done < "$names_file"
}

function prepare_targets {
    [[ $have_targets_fa = n ]] || return 0
    targets_fa="${output}/__targets__.fa"
    [[ ! -f "${targets_fa}" ]] || return 0

    local targets_tmp="${targets_fa}.tmp"
    local targets_lock="${targets_fa}.lock"

    msg "Extracting target sequences"
    # Try to take lock 10 times.
    for i in {1..10}; do
        if ! ( set -C; 2>/dev/null > "$targets_lock" ); then
            # If lock exists, wait 30 seconds until it is released.
            timeout 30 sh -c "while [[ -f \"${targets_lock}\" ]]; do sleep 0.1; done"
            [[ $? -eq 0 ]] || panic "Lock file ${targets_lock} exists for too long, perhaps relevant process was killed"
            [[ ! -f "${targets_fa}" ]] || return 0
        else
            trap 'rm -f "${targets_lock}"; exit 1' INT TERM ERR EXIT
            cat "${targets_bed}" | while read chrom start end name extra; do
                samtools faidx "${reference}" "${chrom}:$((start+1))-${end}" | \
                    seqtk seq -U -l 120 | \
                    sed "1c>$name"
            done > "${targets_tmp}"
            mv "${targets_tmp}" "${targets_fa}"
            rm -f "${targets_lock}"
            trap - INT TERM ERR EXIT
            return 0
        fi
    done
    panic "Most likely, reference target sequences could not be extracted in other instances of this program"
}

function process_assembly {
    local arg="$1"
    local agc_file="${2-}"
    local genome_name
    [[ -n "$agc_file" ]] && genome_name="$arg" || genome_name="$(basename "${arg%.fa*}")"

    local short_name
    # :- if unset or empty, use $genome_name
    short_name="${names["$genome_name"]:-"$genome_name"}"
    [[ "$short_name" != skip ]] || return 0

    local prefix="${output1}/${short_name}"
    local ok_file="${prefix}.ok"
    local lock_file="${prefix}.lock"
    # Needed to account for all unfinished threads.
    local todo_file="${prefix}.todo"
    [[ ! -f "$ok_file" ]] || return 0
    ( set -C; 2>/dev/null > "$lock_file" ) || return 0
    trap 'rm -f "${lock_file}"; exit 1' INT TERM ERR EXIT
    touch "${todo_file}"

    # ===== START ======
    msg "Processing $short_name"
    mkdir -p "$prefix"

    local genome_fasta
    if [[ -n "$agc_file" ]]; then
        msg "    Extracting assembly from AGC file"
        genome_fasta="${prefix}.fa"
        agc getset "$agc_file" "$genome_name" > "${genome_fasta}"
        samtools faidx "${genome_fasta}"
    else
        genome_fasta="$arg"
    fi

    local paf_filename="${prefix}.paf.gz"
    if [[ ! -f "$paf_filename" ]]; then
        msg "    Mapping targets to the assembly"
        minimap2 "${minimap2_args[@]}" "$genome_fasta" "$targets_fa" 2> /dev/null | \
            gzip > "${paf_filename}.tmp" \
            && mv "${paf_filename}"{.tmp,}
    fi

    msg "    Extracting target subsequences"
    "$SCRIPT_DIR/inner/merge_hits.py" "$paf_filename" \
        -b "$targets_bed" -g "$short_name" \
        -d "$distance" -l "$min_len" -s "$min_simil" \
        -o "${prefix}.bed.gz" -c "${prefix}.copy_num.csv.gz"
        2> "${prefix}.warnings.csv"
    # At this point, $prefix.bed.gz will have columns
    # chrom, start, end, strand (+/-), target name, length fraction, similarity.

    rm -f "${prefix}/"*.fa{,.gz}
    # Take $count first entries,
    #     convert into: region, strand faidx argument ("-i" or ""), suffix ("" or "-<INDEX>").
    # Then, fetch regions from the current assembly.
    zcat "${prefix}.bed.gz" | \
        awk -F$'\t' -v name="$short_name" -v max_count="$count" \
            -v cn_csv="${prefix}.copy_num.csv" \
            'BEGIN{OFS=";"} {
                target = $5;
                ind = ++target_cn[target];
                if (ind <= max_count) {
                    region = $1 ":" ($2+1) "-" $3;
                    strand_arg = $4 == "+" ? "" : "-i";
                    suffix = ind == 1 ? "" : ("-" ind);
                    print target, region, strand_arg, suffix
                }
            }' | \
        while IFS=";" read target region strand_arg suffix; do
            samtools faidx "$genome_fasta" "$region" $strand_arg | \
                seqtk seq -U -l 120 | \
                sed "1c>${short_name}${suffix} ${region}" >> "${prefix}/${target}.fa"
        done
    # This way we don't have to worry about zero matches.
    find "${prefix}" -name "*.fa" -exec gzip {} ';'

    [[ -z "$agc_file" ]] || rm "${genome_fasta}"{,.fai}
    # ===== END ======

    touch "${ok_file}"
    rm -f "${lock_file}" "${todo_file}"
    trap - INT TERM ERR EXIT
}

function process_assemblies {
    for curr in "${input[@]}"; do
        if [[ -d "$curr" ]]; then
            for filename in "${curr}"/*.fa{,sta}{,.gz}; do
                process_assembly "$filename"
            done
        elif [[ "$curr" =~ \.fa(sta)?(\.gz)?$ ]]; then
            process_assembly "$curr"
        elif [[ "$curr" =~ \.agc$ ]]; then
            agc listset "$curr" | while read assembly; do
                process_assembly "$assembly" "$curr"
            done
        else
            panic "Unexpected value --input $curr, expected either a directory, `
                `FASTA or AGC file (with lowercase extension)."
        fi
    done
}

function combine_locus {
    local target="$1"

    local prefix="${output2}/${target}"
    local ok_file="${prefix}.ok"
    local lock_file="${prefix}.lock"
    local todo_file="${prefix}.todo"
    if [[ -f "$ok_file" ]]; then
        if [[ "$ok_file" -ot "$latest_ok_file" ]]; then
            rm "$ok_file"
        else
            return 0
        fi
    fi
    ( set -C; 2>/dev/null > "$lock_file" ) || return 0
    trap 'rm -f "${lock_file}"; exit 1' INT TERM ERR EXIT
    touch "${todo_file}"

    # ===== START ======
    msg "    Combining haplotypes for ${target}"
    # Use `find` so that we don't exceed the max number of arguments.
    find "$output1" -mindepth 2 -maxdepth 2 -name "${target}.fa.gz" | \
        xargs -P1 -n 50 cat > "${output2}/${target}.fa.gz"
    # ===== END ======

    touch "${ok_file}"
    rm -f "${lock_file}" "${todo_file}"
    trap - INT TERM ERR EXIT
}

function combine_panels {
    local timeout_sec
    timeout_sec=$((timeout * 60))
    local ready=n
    for i in $(seq 0 "$timeout_sec"); do
        local unfinished
        unfinished="$(find "${output1}" -mindepth 1 -maxdepth 1 -name "*.todo" -print -quit)"
        if [[ -z "$unfinished" ]]; then
            ready=y
            break
        else
            if [[ "$i" = 0 ]]; then
                msg "Waiting $timeout minutes for unfinished jobs (such as $unfinished)"
            fi
            sleep 1
        fi
    done
    [[ "$ready" == y ]] || \
        panic "Cannot combine reference panels: exceeded timeout (-T) and there are unfinished jobs"

    msg "Combining reference panels"
    latest_ok_file="$(find "${output1}" -mindepth 1 -maxdepth 1 -name "*.ok" -printf "%T@\t%p\n" | sort -k1,1gr | \
        head -n1 | cut -f2)"
    cut -f4 "$targets_bed" | while read target; do
        combine_locus "$target"
    done

    # Test will return false if the files don't exist
    if [[ "${output2}/targets.bed" -nt "$latest_ok_file" && "${output2}/warnings.csv" -nt "$latest_ok_file" ]]; then
        return 0
    fi
    local lock_file="${output2}/targets.lock"
    ( set -C; 2>/dev/null > "$lock_file" ) || return 0
    trap 'rm -f "${lock_file}"; exit 1' INT TERM ERR EXIT

    cut -f-4 "$targets_bed" | awk -F$'\t' 'BEGIN{OFS=FS} { print $0, ($4 ".fa.gz") }' > "${output2}/targets.bed.tmp"
    mv "${output2}/targets.bed"{.tmp,}
    find "$output1" -mindepth 1 -maxdepth 1 -name "*.warnings.csv" | \
        xargs -P1 -n 50 cat | sort > "${output2}/warnings.csv.tmp"
    mv "${output2}/warnings.csv"{.tmp,}
    find "$output1" -mindepth 1 -maxdepth 1 -name "*.copy_num.csv.gz" | \
        xargs -P1 -n 50 cat | sort > "${output2}/copy_num.csv.gz.tmp"
    mv "${output2}/copy_num.csv.gz"{.tmp,}

    rm -f "${lock_file}"
    trap - INT TERM ERR EXIT
}

function check_completion {
    if [[ -f "${output2}/targets.bed"
            && -z "$(find "$output1" "$output2" -mindepth 1 -maxdepth 1 -name "*.todo" -print -quit)" ]]; then
        msg "All jobs completed"
        touch "${output2}/ok"
    else
        msg "Some jobs incomplete, please wait for other instances to complete, `
            `run additional instances, or check for error messages"
    fi
}

setup_colors
parse_args "$@"
output1="$output/by_assembly"
output2="$output/panels"
mkdir -p "$output1" "$output2"

declare -A names
load_names

prepare_targets
# zcat -f opens plain files as well. sed -n does not print by default.
readarray -t target_names < <(zcat -f "$targets_fa" | sed -n 's/>//p' | sort -u)
process_assemblies
combine_panels
check_completion
