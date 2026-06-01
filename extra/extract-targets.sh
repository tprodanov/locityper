#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]:-$0}")"

function help_message {
  cat <<HELP
Usage: $SCRIPT_NAME (-a FILE | -g DIR) -t FILE -o DIR [-d INT] [args] [-- minimap-args]

Maps target sequences to assembly genomes and extracts corresponding subregions.
Multiple instances of this script can be run in parallel on the same output directory,
as it creates "*.lock" files for haplotypes in progress and ".ok" files for the finished haplotypes.

Available options:
    -a, --agc      FILE  Input AGC file.
    -g, --genomes  DIR   Directory with various genome assemblies (.fa[.gz]).
                         Mutually exclusive with -a/--agc.
    -n, --names    FILE  Optional: replace genome names (first column) with another name (second column).
                         In case of -g/--genomes, first column should match file basename without extension.
    -t, --targets  FILE  FASTA file with target sequences. Name lines should not contain spaces.
    -o, --output   DIR   Output directory.
    -d, --distance INT   Merge PAF entries if distance is smaller than INT [${distance}].
    -f, --min-frac NUM   Output regions longer than NUM * target length [${min_frac}].
    -N, --count    INT   Extract top N sequences (by length) for each haplotype and target [${count}].
    -h, --help           Print this help and exit.

Provide minimap2 arguments after --
    Default arguments are "-cx asm20 -t 3 -N 10 -p 0.5"
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

function parse_params {
    min_frac=0.7
    distance=5000
    count=1
    names_file=

    ARGS="$(getopt -o a:g:n:t:o:d:f:N:h \
        --long agc:,genomes:,names:,targets:,output:,distance:,min-frac:,count:,help \
        --name "$SCRIPT_NAME" -- "$@")"
    eval set -- "$ARGS"
    while :; do
        case "$1" in
            -a | --agc )
                agc_file="$2"; shift 2 ;;
            -g | --genomes )
                genomes_dir="$2"; shift 2 ;;
            -t | --targets )
                targets_file="$2"; shift 2 ;;
            -n | --names )
                names_file="$2"; shift 2 ;;
            -o | --output )
                output="$2"; shift 2 ;;
            -d | --distance )
                distance="$2"; shift 2 ;;
            -f | --min-frac )
                min_frac="$2"; shift 2 ;;
            -N | --count )
                count="$2"; shift 2 ;;
            -h | --help)
                help_message; exit 0;
                ;;
            -- ) shift; break ;;
            * )  panic "Unexpected argument $1" ;;
        esac
    done

    if [[ ${#@} -ne 0 ]]; then
        minimap2_args=( "$@" )
    else
        minimap2_args=( -cx asm20 -t 3 -N 10 -p 0.5 )
    fi

    [[ ! -z "${targets_file-}" ]] || panic "Missing required parameter -t/--targets"
    [[ -f "${targets_file}" ]] || panic "Targets file ${targets_file} not found"
    [[ ! -z "${output-}" ]]   || panic "Missing required parameter -o/--output"

    [[ -z "${agc_file-}" ]] && have_agc=n || have_agc=y
    [[ -z "${genomes_dir-}" ]] && have_genomes=n || have_genomes=y
    [[ $have_agc != $have_genomes ]] || panic "Require either -a or -g, but not both"
}

function load_names {
    [[ ! -z "$names_file" ]] || return 0
    while read name upd_name; do
        names["$name"]="$upd_name"
    done < "$names_file"
}

function process_genome {
    local arg="$1"
    local genome_name
    [[ $have_agc = y ]] && genome_name="$arg" || genome_name="$(basename "${arg%.fa*}")"

    local short_name
    # :- if unset or empty, use $genome_name
    short_name="${names["$genome_name"]:-"$genome_name"}"

    local prefix="${output}/${short_name}"
    local ok_file="${prefix}.ok"
    local lock_file="${prefix}.lock"
    if [[ -f "$ok_file" ]]; then
        return
    fi
    if ! ( set -C; 2>/dev/null > "$lock_file" ); then
        return
    fi
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
        minimap2 "${minimap2_args[@]}" "$genome_fasta" "$targets_file" 2> /dev/null | \
            gzip > "${paf_filename}.tmp" \
            && mv "${paf_filename}"{.tmp,}
    fi

    msg "    Extracting subsequences"
    for target in "${target_names[@]}"; do
        # * Take PAF for given target and convert to BED file
        #       columns: chrom, start, end, strand (±1 * length), target length
        # * Sort and merge, sum fourth column
        # * Convert columns: chrom, start, end, strand (+/-), target name, length / target_length.
        zcat "${paf_filename}" | \
            awk -F$'\t' -v target="$target" \
                'BEGIN{OFS=FS} $1 == target { print $6, $8, $9, ($5 == "+" ? 1 : -1) * ($4 - $3), $2 }' | \
            sort -k1,1V -k2,2n | \
            bedtools merge -d "$distance" -c 4,5 -o sum,distinct 2> /dev/null | \
            awk -F$'\t' -v target="$target" \
                'BEGIN{OFS=FS} { print $1, $2, $3, ($4 >= 0 ? "+" : "-"), target, ($3-$2) / $5 }' |
                sort -k6,6gr > "${prefix}/${target}.bed"

        # Take $count first entries,
        #     convert into: region, strand faidx argument ("-i" or ""), suffix ("" or "-<INDEX>").
        # Then, fetch regions from the current assembly.
        awk -F$'\t' -v count="$count" -v min_frac="$min_frac" \
            'BEGIN{OFS=";"} NR <= count && $6 >= min_frac {
                print ($1 ":" ($2+1) "-" $3), $4, $3 - $2
            }' "${prefix}/${target}.bed" | \
            while IFS=";" read region strand_arg suffix; do
                samtools faidx "$genome_fasta" "$region" $strand_arg | \
                    sed "1c>${short_name}${suffix} ${region}"
            done | gzip > "${prefix}/${target}.fa.gz"
    done

    [[ $have_agc = n ]] || rm "${genome_fasta}"{,.fai}

    cat "${prefix}"/*.bed | gzip > "${prefix}.bed.gz"
    rm "${prefix}"/*.bed
    # ===== END ======

    touch "${ok_file}"
    rm -f "${lock_file}"
    trap - INT TERM ERR
}

setup_colors
parse_params "$@"
declare -A names
load_names

# zcat -f opens plain files as well. sed -n does not print by default.
readarray -t target_names < <(zcat -f "$targets_file" | sed -n 's/>//p' | sort -u)

mkdir -p "$output"
if [[ $have_agc = y ]]; then
    agc listset "$agc_file" | while read genome; do
        process_genome "$genome"
    done
else
    for filename in "${genomes_dir}"/*.fa{,sta}{,.gz}; do
        process_genome "$filename"
    done
fi
