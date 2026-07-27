#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]:-$0}")"

function help_message {
  cat <<HELP
Usage: $SCRIPT_NAME dir [args]

│ Extract reads, recruited by Locityper genotyping.

Input/output arguments:
        input     DIR  Locityper genotyping directory.
    -d, --edit    INT  Output reads with this maximum edit distance [$edit].
    -f, --frac    NUM  Instead of edit distance, filter by edit distance / read length.
    -b, --both    y|n  Should both (y) or only one (n) of the reads pass the filter [$both].
    -z, --gzip    y|n  Should the output files be gzipped? [$compress].
    -o, --output  DIR  Optional: output directory.
                       If not specified, read files are placed in the same output directory.

Other arguments:
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
    edit=4
    frac=
    both=y
    compress=y
    output=

    ARGS="$(getopt -o d:f:b:o:h --long "edit:,frac:,both:,output:,help" --name "$SCRIPT_NAME" -- "$@")"
    eval set -- "$ARGS"
    while :; do
        case "$1" in
            -d | --edit )
                edit="$2"; shift 2 ;;
            -f | --frac )
                frac="$2"; shift 2 ;;
            -b | --both )
                both="$2"; shift 2 ;;
            -z | --gzip )
                compress="$2"; shift 2 ;;
            -o | --output )
                output="$2"; shift 2 ;;
            -h | --help)
                 help_message; exit 0 ;;
            -- ) shift; break ;;
            * )  panic "Unexpected argument $1" ;;
        esac
    done

    [[ $# -eq 1 ]] || panic "Missing/too many positional arguments"
    input="$1"

    [[ "$both" = y || "$both" = n ]] || panic '--both must be either "y" or "n"'
    [[ "$compress" = y || "$compress" = n ]] || panic '--compress must be either "y" or "n"'
}

function process_dir {
    local in_dir
    in_dir="$1"

    local prefix
    if [[ -n "$output" ]]; then
        prefix="${output}/$(basename "$in_dir").reads"
    else
        prefix="${in_dir}/reads"
    fi

    echo "Extracting reads $in_dir -> ${prefix}*.fq"
    samtools view -h -F 2304 "${in_dir}/aln.bam" | \
    awk -F "\t" -v edit_thresh="$edit" -v frac_thresh="$frac" -v both="$both" 'BEGIN {
            OFS = FS;
            both = (both == "y");
            UNDEF = 2147483647;
        } {
            if ($1 ~ /^@/) { print }
            else {
                edit = UNDEF;
                for (i = 12; i <= NF; ++i) {
                    if ($i ~ /^NM/) {
                        edit = int(substr($i, 6));
                        break;
                    }
                }
                keep_this = 0;
                if ($6 != "*" && edit == UNDEF) {
                    printf("Read mapping for %s does not have edit distance (NM)\n", $1) > "/dev/stderr";
                } else if (frac_thresh != "" && $10 == "*") {
                    printf("Sequence missing for read %s\n", $1) > "/dev/stderr";
                } else {
                    keep_this = frac_thresh == "" ? (edit <= edit_thresh) : (edit / length($10) <= frac_thresh);
                }

                if (first != "" && first_name != $1) {
                    if (keep_first) {
                        print first;
                    }
                    first = "";
                }

                if (first == "") {
                    first_name = $1;
                    keep_first = keep_this;
                    $2 = or($2, 128);
                    first = $0;
                } else {
                    $2 = or($2, 64);
                    if ((keep_this && keep_first) || (!both && keep_this) || (!both && keep_first)) {
                        print first;
                        print $0;
                    }
                    first = "";
                }
            }
        } END {
            if (first != "" && keep_first) {
                print first;
            }
        }' | \
        samtools view -bo "${prefix}.keep.bam"
    samtools fastq -1 "${prefix}1.fq" -2 "${prefix}2.fq" "${prefix}.keep.bam"
    rm "${prefix}.keep.bam"

    [[ -s "${prefix}2.fq" ]] || rm "${prefix}2.fq"

    if [[ "$compress" = y ]]; then
        gzip -f "${prefix}1.fq"
        [[ ! -f "${prefix}2.fq" ]] || gzip -f "${prefix}2.fq"
    fi
}

setup_colors
parse_args "$@"
[[ -z "$output" ]] || mkdir -p "$output"
find "$input/loci" -mindepth 1 -maxdepth 1 -type d | while read d; do
    process_dir "$d"
done
