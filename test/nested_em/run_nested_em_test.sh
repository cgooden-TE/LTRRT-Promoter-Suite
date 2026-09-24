#!/usr/bin/env bash
# Benchmarks --assign em against the innermost rule on simulated nested loci, and checks that the
# default mode is unchanged relative to the unpatched script. The optional argument is either a git
# ref (default HEAD, i.e. run before committing) or a path to a saved copy of the unpatched script.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="$HERE/output"
BASE_REF="${1:-HEAD}"
[[ -f "$BASE_REF" ]] && BASE_REF="$(cd "$(dirname "$BASE_REF")" && pwd)/$(basename "$BASE_REF")"
mkdir -p "$OUT"
cd "$OUT"

run() {  # run <script> <simdir> <outdir> <mode> [extra args]
    local script=$1 sim=$2 out=$3 mode=$4; shift 4
    mkdir -p "$out"
    python3 "$script" --gff "$sim/sim.gff" --gene-gff "$sim/sim_genes.gff3" --bam "$sim/sim.bam" \
        --output "$out/iso" --tss_out "$out/tss.tsv" --gene_out "$out/gene.tsv" \
        --include-partial --assign "$mode" "$@" > "$out/log.txt" 2>&1
}

python3 "$HERE/sim_nested.py" --outdir sim_clean > /dev/null
python3 "$HERE/sim_nested.py" --outdir sim_hard --hard > /dev/null

# --- regression: default mode vs unpatched script ---
if [[ -f "$BASE_REF" ]]; then
    cp "$BASE_REF" IsoClassifier_base.py
else
    git -C "$ROOT" show "$BASE_REF:IsoClassifier.py" > IsoClassifier_base.py
fi
# A baseline from after the EM commit imports nested_em from its own directory; innermost-mode output
# does not depend on its contents, so the current copy serves.
cp "$ROOT/nested_em.py" "$ROOT/tss_consensus.py" .
mkdir -p base && python3 IsoClassifier_base.py --gff sim_clean/sim.gff --gene-gff sim_clean/sim_genes.gff3 \
    --bam sim_clean/sim.bam --output base/iso --tss_out base/tss.tsv --gene_out base/gene.tsv \
    --include-partial > base/log.txt 2>&1
# The pre-patch script used the modal TSS; --tss-method mode reproduces it for the byte comparison.
run "$ROOT/IsoClassifier.py" sim_clean default innermost --tss-method mode
status=0
for f in iso_sense.isoforms.tsv iso_antisense.isoforms.tsv tss.tsv gene.tsv \
         iso_10site.cleavage_summary.tsv iso_10site.tss_summary.tsv \
         iso_primary_tss_density.tsv iso_gene_primary_density.tsv; do
    if cmp -s "base/$f" "default/$f"; then echo "regression OK: $f"; else echo "REGRESSION: $f"; status=1; fi
done
python3 - <<'EOF' || status=1
import pandas as pd, sys
ok = True
for suf in ['_3p_softclip_per_read_spliced.tsv', '_3p_softclip_per_read_nonspliced.tsv',
            '_te_exon_stats_per_read.tsv', '_gene_exon_stats_per_read.tsv']:
    a = pd.read_csv('base/iso' + suf, sep='\t', dtype=str)
    b = pd.read_csv('default/iso' + suf, sep='\t', dtype=str)
    same = a.equals(b[a.columns])
    print(f"regression {'OK' if same else 'FAILED'}: per-read {suf} (original columns)")
    ok &= same
sys.exit(0 if ok else 1)
EOF

# --- benchmarks ---
run "$ROOT/IsoClassifier.py" sim_clean clean_em em
run "$ROOT/IsoClassifier.py" sim_hard hard_inner innermost
run "$ROOT/IsoClassifier.py" sim_hard hard_em em
run "$ROOT/IsoClassifier.py" sim_hard hard_em_t4 em --threads 4
if cmp -s hard_em/iso_em_units.tsv hard_em_t4/iso_em_units.tsv; then
    echo "determinism OK: --threads 4 matches --threads 1"
else
    echo "DETERMINISM FAILED"; status=1
fi
python3 "$HERE/eval_assign.py" default/iso=clean_innermost clean_em/iso=clean_em \
    hard_inner/iso=hard_innermost hard_em/iso=hard_em
exit $status
