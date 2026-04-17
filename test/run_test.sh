#!/usr/bin/env bash
# =============================================================================
# LTRRT-Promoter-Suite Test Runner
# =============================================================================
# Tests the full pipeline: AccuMap -> IsoClassifier -> WindowScrubber -> Query_WSDB
#
# This script uses a small test dataset (200 ONT reads, 1 Mb reference region
# from B73 chr1) to verify that all tools run correctly after installation.
#
# Prerequisites:
#   conda activate LTRPromSuite_pipeline
#
# Usage:
#   cd LTRRT-Promoter-Suite
#   bash test/run_test.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE_DIR="$(dirname "$SCRIPT_DIR")"
OUT_DIR="$SCRIPT_DIR/output"

mkdir -p "$OUT_DIR"

echo "============================================="
echo " LTRRT-Promoter-Suite Test"
echo "============================================="
echo "Test data directory: $SCRIPT_DIR"
echo "Output directory:    $OUT_DIR"
echo ""

# -------------------------------------------------------
# Step 1: AccuMap (PyChopper -> Cutadapt -> Minimap2)
# -------------------------------------------------------
echo "---------------------------------------------"
echo " Step 1: AccuMap"
echo "---------------------------------------------"
python3 "$SUITE_DIR/AccuMap.py" \
    --fq "$SCRIPT_DIR/test_reads.fastq" \
    --sample "$OUT_DIR/test_sample" \
    --ref "$SCRIPT_DIR/test_reference.fa" \
    --run_pyc \
    --run_cut \
    --run_map \
    --pyc_threads 2 \
    --cut_threads 2 \
    --map_threads 2

# Verify key output exists
if [ ! -f "$OUT_DIR/test_sample.STtagged.sorted.bam" ]; then
    echo "ERROR: AccuMap did not produce expected BAM output."
    exit 1
fi
echo "AccuMap completed successfully."
echo ""

# -------------------------------------------------------
# Step 2: IsoClassifier
# -------------------------------------------------------
echo "---------------------------------------------"
echo " Step 2: IsoClassifier"
echo "---------------------------------------------"
python3 "$SUITE_DIR/IsoClassifier.py" \
    --gff "$SCRIPT_DIR/test_annotations.gff" \
    --bam "$OUT_DIR/test_sample.STtagged.sorted.bam" \
    --min_mapq 30 \
    --threads 1 \
    --output "$OUT_DIR/test_ltr_isoforms" \
    --tss_out "$OUT_DIR/test_ltr_tss_summary.tsv" \
    --gene_out "$OUT_DIR/test_gene_summary.tsv" \
    --genome-fasta "$SCRIPT_DIR/test_reference.fa"

# Verify key outputs exist
if [ ! -f "$OUT_DIR/test_ltr_isoforms.tsv" ]; then
    echo "ERROR: IsoClassifier did not produce expected isoform output."
    exit 1
fi
# Verify U3/promoter extraction outputs (gene outputs should always exist;
# LTR outputs depend on whether elements had classified reads)
for f in test_ltr_isoforms_gene_2kb_proms.bed \
         test_ltr_isoforms_gene_2kb_proms.fa \
         test_ltr_isoforms_gene_dummy_u3.fa \
         test_ltr_isoforms_u3_seqs.fa \
         test_ltr_isoforms_ltr_seqs.fa; do
    if [ ! -f "$OUT_DIR/$f" ]; then
        echo "WARNING: U3/promoter output $f not found (may be empty if no reads classified)."
    fi
done
echo "IsoClassifier completed successfully."
echo ""

# -------------------------------------------------------
# Step 3: WindowScrubber
# -------------------------------------------------------
# WindowScrubber requires LTR/U3 FASTA and TSS data.
# Use the pre-built test files (from real data) rather than
# IsoClassifier output, since the small read set may not
# produce sufficient LTR classifications.
echo "---------------------------------------------"
echo " Step 3: WindowScrubber"
echo "---------------------------------------------"
python3 "$SUITE_DIR/WindowScrubber.py" \
    -l "$SCRIPT_DIR/test_ltr_seqs.fa" \
    -u3 "$SCRIPT_DIR/test_u3_seqs.fa" \
    -t "$SCRIPT_DIR/test_tss_summary.tsv" \
    -db "$OUT_DIR/test_motif_hits.db"

if [ ! -f "$OUT_DIR/test_motif_hits.db" ]; then
    echo "ERROR: WindowScrubber did not produce expected database."
    exit 1
fi
echo "WindowScrubber completed successfully."
echo ""

# -------------------------------------------------------
# Step 4: Query_WSDB
# -------------------------------------------------------
echo "---------------------------------------------"
echo " Step 4: Query_WSDB"
echo "---------------------------------------------"

echo "  4a. Database statistics:"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$OUT_DIR/test_motif_hits.db" \
    --stats

echo ""
echo "  4b. Best TATA per element:"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$OUT_DIR/test_motif_hits.db" \
    --best-per-element TATA \
    --order-by score \
    -o "$OUT_DIR/test_best_tata.tsv"

echo ""
echo "  4c. Best hits for all motifs (wide format):"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$OUT_DIR/test_motif_hits.db" \
    --best-all-motifs \
    -o "$OUT_DIR/test_best_all_motifs.tsv"

echo ""
echo "  4d. CA dinucleotide summary:"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$OUT_DIR/test_motif_hits.db" \
    --ca-summary \
    -o "$OUT_DIR/test_ca_runs.tsv"

echo "Query_WSDB completed successfully."
echo ""

# -------------------------------------------------------
# Summary
# -------------------------------------------------------
echo "============================================="
echo " All tests passed!"
echo "============================================="
echo ""
echo "Output files:"
ls -lh "$OUT_DIR/"
