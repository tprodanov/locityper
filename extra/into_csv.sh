#!/bin/sh

set -euo pipefail

if [ $# = 0 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    2>&1 echo "Usage:  into_csv.sh DIR OUT.csv.gz"
    2>&1 echo "    where DIR/SAMPLE contains Locityper results for SAMPLE"
    2>&1 echo "    all found DIR/*/loci/*/res.json.gz will be taken"
    exit 0
fi

input_dir="$1"
output_csv="$(realpath "$2")"

cd "$input_dir"
(echo -e "sample\tlocus\tgenotype";
find -name res.json.gz | \
    parallel -t -P8 \
        --rpl '{firstdirname} s:^\./::; s:/.*::' \
        --rpl '{lastdirname} $Global::use{"File::Basename"} ||=
            eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
        echo -ne {firstdirname}'"\t"'{lastdirname} ';' \
        zgrep -m1 genotype {} | \
    sed 's/    "genotype": "/\t/; s/",//') | \
    gzip > "$output_csv"
