#!/usr/bin/env bash
# Checks the consensus TSS call and the --write-state / --merge-sheet path on simulated nested loci:
#   1. unit tests of tss_consensus.py
#   2. consensus and mode runs credit identical reads (only TSS columns may differ)
#   3. merging one sample's state reproduces that sample's own run exactly
#   4. merging two halves of a BAM matches one run over both halves (rows compared sorted)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="$HERE/output"
ISO="$ROOT/IsoClassifier.py"
rm -rf "$OUT" && mkdir -p "$OUT" && cd "$OUT"
status=0

python3 "$HERE/test_tss_consensus.py" || status=1

python3 "$ROOT/test/nested_em/sim_nested.py" --outdir sim --hard > /dev/null
SIM="--gff sim/sim.gff --gene-gff sim/sim_genes.gff3 --include-partial"

run() {  # run <outdir> [args]: one IsoClassifier run with the standard output names
    local o=$1; shift
    mkdir -p "$o"
    python3 "$ISO" $SIM --output "$o/iso" --tss_out "$o/tss.tsv" --gene_out "$o/gene.tsv" "$@" \
        > "$o/log.txt" 2>&1 || { echo "RUN FAILED: $o"; tail -5 "$o/log.txt"; exit 1; }
}

# Summary outputs a merge writes; per-read and EM tables stay with each sample.
SUMMARIES="tss.tsv gene.tsv iso_sense.isoforms.tsv iso_antisense.isoforms.tsv
  iso_10site.tss_summary.tsv iso_sense.tss_summary.tsv iso_antisense.tss_summary.tsv
  iso_10site.cleavage_summary.tsv iso_10site.gene_summary.tsv iso_primary_tss_density.tsv
  iso_secondary_tss_density.tsv iso_gene_primary_density.tsv iso_gene_secondary_density.tsv
  iso_tss_consensus.tsv iso_tss_by_library.tsv"

# --- 2. consensus vs mode: same reads credited ---
run mode --bam sim/sim.bam --assign em --tss-method mode
run cons --bam sim/sim.bam --assign em --write-state
python3 - <<'EOF' || status=1
import pandas as pd, sys
ok = True
tss_cols = {'TSS_Peak', 'TSS_Called', 'TSS_Peak_Frac'}
for o in ('sense', 'antisense'):
    a = pd.read_csv(f'mode/iso_{o}.isoforms.tsv', sep='\t').sort_values(['Feature', 'Orientation'])
    b = pd.read_csv(f'cons/iso_{o}.isoforms.tsv', sep='\t').sort_values(['Feature', 'Orientation'])
    if a.empty and b.empty:
        print(f"ok   consensus vs mode, {o}: no records")
        continue
    cols = [c for c in a.columns if c not in tss_cols]
    same = a[cols].reset_index(drop=True).equals(b[cols].reset_index(drop=True))
    moved = (a['TSS_Peak'].values != b['TSS_Peak'].values).sum()
    print(f"{'ok  ' if same else 'FAIL'} consensus vs mode, {o}: counts identical; "
          f"{moved}/{len(a)} TSS_Peak moved")
    ok &= same
sys.exit(0 if ok else 1)
EOF

# --- 3. merge of one state == that sample's own run ---
printf "Source\tLibrary\tPlatform\ncons/iso\t\t\n" > one.tsv
run merge1 --merge-sheet one.tsv
for f in $SUMMARIES; do
    if cmp -s "cons/$f" "merge1/$f"; then echo "ok   merge-one: $f"
    else echo "FAIL merge-one: $f"; status=1; fi
done

# --- 4. two halves merged == one run over both halves ---
samtools view -h sim/sim.bam | awk -F'\t' '/^@/ || NR % 2' | samtools view -b -o half1.bam -
samtools view -h sim/sim.bam | awk -F'\t' '/^@/ || !(NR % 2)' | samtools view -b -o half2.bam -
samtools index half1.bam && samtools index half2.bam
run h1 --bam half1.bam --write-state
run h2 --bam half2.bam --write-state
run both --bam half1.bam half2.bam
printf "Source\tLibrary\tPlatform\nh1/iso\thalf1\tont\nh2/iso\thalf2\tont\n" > two.tsv
run merge2 --merge-sheet two.tsv
for f in $SUMMARIES; do
    if cmp -s <(sort "both/$f") <(sort "merge2/$f"); then echo "ok   merge-two vs direct: $f"
    else echo "FAIL merge-two vs direct: $f"; status=1; fi
done

# --- guard rails ---
if python3 "$ISO" --gff sim/sim.gff --bam sim/sim.bam --merge-sheet two.tsv --output x \
        --tss_out x --gene_out x > /dev/null 2>&1; then
    echo "FAIL --bam with --merge-sheet accepted"; status=1
else echo "ok   --bam with --merge-sheet refused"; fi

exit $status
