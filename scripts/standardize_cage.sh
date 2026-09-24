#!/usr/bin/env bash
# Convert a bedtools-intersect CAGE GFF (feature GFF + CAGE cluster columns)
# into a keyed TSV that DECLTR reads with the keyed_tsv preset.
#
# Input columns: 1-9 feature GFF, 10 cluster start, 11 cluster end, 12 strand,
# 13 dominant TSS, 14 cluster attributes carrying Shape=...
# Output: key, CAGE_Start, CAGE_End, CAGE_Str, CAGE_dTSS, CAGE_Shape,
# first occurrence per feature key.
#
# Usage: standardize_cage.sh in.formatted.gff out.tsv

set -euo pipefail
[ $# -eq 2 ] || { echo "Usage: $0 in.gff out.tsv" >&2; exit 1; }

awk -F'\t' 'BEGIN { OFS = "\t"; print "key", "CAGE_Start", "CAGE_End", "CAGE_Str", "CAGE_dTSS", "CAGE_Shape" }
/^#/ { next }
NF >= 14 {
    id = ""; if (match($9, /(^|;)ID=[^;]+/)) { id = substr($9, RSTART, RLENGTH); sub(/^;?ID=/, "", id) }
    if (id == "" || (id in seen)) next
    seen[id] = 1
    shape = "NA"; if (match($14, /Shape=[^;]+/)) { shape = substr($14, RSTART + 6, RLENGTH - 6) }
    print id, $10, $11, $12, $13, shape
}' "$1" > "$2"
echo "Wrote $(($(wc -l < "$2") - 1)) rows to $2" >&2
