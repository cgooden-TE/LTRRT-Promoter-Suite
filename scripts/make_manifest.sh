#!/usr/bin/env bash
# Draft a DECLTR manifest from input directories and files.
#
# Emits one row per sample-file pairing with sample_id, platform, preset, role,
# and scope guessed from file names. Tissue and replicate are guessed only where
# the name makes it obvious; review and edit every row before running DECLTR.
#
# Usage:
#   make_manifest.sh -o manifest.tsv [chip=DIR] [umr=FILE|DIR] [isoforms=DIR]
#                    [motif=DIR] [cage=DIR|FILE] [illumina=FILE]
#
# Paths in the manifest are written as given (absolute paths recommended).

set -euo pipefail

out=""
declare -a specs=()
while [ $# -gt 0 ]; do
    case "$1" in
        -o) out="$2"; shift 2 ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        *=*) specs+=("$1"); shift ;;
        *) echo "Unrecognised argument: $1" >&2; exit 1 ;;
    esac
done
[ -n "$out" ] || { echo "-o manifest.tsv is required" >&2; exit 1; }
[ ${#specs[@]} -gt 0 ] || { echo "No inputs given" >&2; exit 1; }

header="sample_id	platform	role	preset	path	column	tissue	replicate	scope	options"
echo "$header" > "$out"

# Append one manifest row; empty fields are allowed.
row() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$out"; }

# Expand a directory into its regular files, or pass a file through.
expand() {
    if [ -d "$1" ]; then find "$1" -maxdepth 1 -type f | sort; else echo "$1"; fi
}

for spec in "${specs[@]}"; do
    kind="${spec%%=*}"
    target="${spec#*=}"
    case "$kind" in
        chip)
            for f in $(expand "$target"); do
                b=$(basename "$f")
                case "$b" in *UMR*) continue ;; esac
                sid=${b%_peaks.intsct.gff}; sid=${sid%.intsct.gff}; sid=${sid%.gff}
                tissue=${sid%%_*}
                rep=$(echo "$sid" | grep -oE '\.[0-9]+$' | tr -d '.' || true)
                row "$sid" chip score intersect_gff_chip "$f" "" "$tissue" "$rep" both ""
            done ;;
        umr)
            for f in $(expand "$target"); do
                b=$(basename "$f")
                sid=${b%_intsct.gff}; sid=${sid%.intsct.gff}; sid=${sid%.gff}
                row "$sid" umr score intersect_gff_umr "$f" "" "" "" both ""
            done ;;
        isoforms)
            for f in $(expand "$target"); do
                b=$(basename "$f")
                case "$b" in *antisense*) continue ;; esac
                platform=pacbio; sid=""; tissue=""; rep=""
                if echo "$b" | grep -qiE '^ONT'; then
                    platform=ont
                    sid=$(echo "$b" | grep -oE '^ONT(_CT[0-9]+|_CTMerge|PB)' || true)
                    rep=$(echo "$sid" | grep -oE 'CT[0-9]+' || true)
                elif echo "$b" | grep -qE '^PB_'; then
                    sid=$(echo "$b" | grep -oE '^PB_[A-Za-z0-9]+' || true)
                    tissue=${sid#PB_}
                fi
                [ -n "$sid" ] || sid=${b%%.*}
                case "$b" in
                    *_Gene_TSS.tsv)          row "$sid" $platform score tss_summary_v1 "$f" "" "$tissue" "$rep" gene "" ;;
                    *_LTR_TSS.tsv)           row "$sid" $platform score tss_summary_v1 "$f" "" "$tissue" "$rep" te "" ;;
                    *_sense.tss_summary.tsv) row "$sid" $platform score tss_summary_v2 "$f" "" "$tissue" "$rep" both "" ;;
                    *_sense.isoforms.tsv)    row "$sid" $platform extra isoforms_v2 "$f" "" "$tissue" "$rep" both "" ;;
                    *isoforms.tsv)           row "$sid" $platform extra isoforms_v1 "$f" "" "$tissue" "$rep" te "" ;;
                    *) echo "Skipping unrecognised isoform file: $b" >&2 ;;
                esac
            done ;;
        motif)
            for f in $(expand "$target"); do
                b=$(basename "$f")
                sid=${b%%_*}
                scope=te; case "$b" in *Gene*|*gene*) scope=gene ;; esac
                row "$sid" motif extra motif_tsv "$f" "" "" "" $scope ""
            done ;;
        cage)
            for f in $(expand "$target"); do
                b=$(basename "$f")
                sid=${b%%_*}
                row "$sid" cage extra keyed_tsv "$f" "" "$sid" "" both "keep_cols=*"
            done ;;
        illumina)
            f="$target"
            head -1 "$f" | tr '\t' '\n' | tr -d '\r' \
              | grep -vxE 'Chr|Source|Name|Start|End|Score|Strand|Phase|Attributes|ID' \
              | while read -r col; do
                    row "$col" illumina score count_matrix "$f" "$col" "" "" both ""
                done ;;
        *) echo "Unknown input kind: $kind" >&2; exit 1 ;;
    esac
done

n=$(($(wc -l < "$out") - 1))
echo "Wrote $n rows to $out. Fill in tissue and replicate before running DECLTR." >&2
