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
if [ ! -f "$OUT_DIR/test_ltr_isoforms_sense.isoforms.tsv" ]; then
    echo "ERROR: IsoClassifier did not produce expected isoform output."
    exit 1
fi
# Verify U3/promoter extraction outputs (gene outputs should always exist;
# LTR outputs depend on whether elements had classified reads)
for f in test_ltr_isoforms_gene_2kb_proms.bed \
         test_ltr_isoforms_gene_2kb_proms.fa \
         test_ltr_isoforms_gene_dummy_u3.fa \
         test_ltr_isoforms_u3_seqs_sense.fa \
         test_ltr_isoforms_ltr_seqs_sense.fa; do
    if [ ! -f "$OUT_DIR/$f" ]; then
        echo "WARNING: U3/promoter output $f not found (may be empty if no reads classified)."
    fi
done
echo "IsoClassifier completed successfully."
echo ""

# -------------------------------------------------------
# Step 3: WindowScrubber
# -------------------------------------------------------
# One run covers every element class. LTR_structural elements whose TSS sits in an
# LTR use the U3 schema (LTR bounds from ltr_seqs_sense.fa); genes, DNA TEs, LINE/SINE
# and LTR fragments use a fixed genome window around the TSS.
echo "---------------------------------------------"
echo " Step 3: WindowScrubber"
echo "---------------------------------------------"
python3 "$SUITE_DIR/WindowScrubber.py" \
    -t "$OUT_DIR/test_ltr_isoforms_sense.tss_summary.tsv" \
    -g "$SCRIPT_DIR/test_reference.fa" \
    -l "$OUT_DIR/test_ltr_isoforms_ltr_seqs_sense.fa" \
    -db "$OUT_DIR/test_motif_hits.db"

if [ ! -s "$OUT_DIR/test_motif_hits.db" ]; then
    echo "ERROR: WindowScrubber did not produce a motif database."
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
DB="$OUT_DIR/test_motif_hits.db"

echo "  4a. Database statistics:"
python3 "$SUITE_DIR/Query_WSDB.py" -db "$DB" --stats

echo ""
echo "  4b. Best TATA per element:"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$DB" \
    --best-per-element TATA \
    --order-by score \
    -o "$OUT_DIR/test_best_tata.tsv"

echo ""
echo "  4c. Best hits for all motifs (wide format), per schema:"
for MODE in U3 TSS_window; do
    python3 "$SUITE_DIR/Query_WSDB.py" \
        -db "$DB" \
        --best-all-motifs \
        --anchor-mode "$MODE" \
        -o "$OUT_DIR/test_${MODE}_best_all_motifs.tsv"
done

echo ""
echo "  4d. Best TATA for genes only:"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$DB" \
    --best-per-element TATA \
    --class Gene \
    -o "$OUT_DIR/test_gene_best_tata.tsv"

echo ""
echo "  4e. CA dinucleotide summary:"
python3 "$SUITE_DIR/Query_WSDB.py" \
    -db "$DB" \
    --ca-summary \
    -o "$OUT_DIR/test_ca_runs.tsv"

echo "Query_WSDB completed successfully."
echo ""

# -------------------------------------------------------
# Step 5: DECLTR (needs the DECLTR-env R environment)
# -------------------------------------------------------
# Runs on the IsoClassifier v2 outputs with a two-row manifest and no
# Illumina/ChIP/UMR inputs, so every platform is exercised as optional.
# Set DECLTR_RSCRIPT to the Rscript of the DECLTR-env if it is not on PATH.
echo "---------------------------------------------"
echo " Step 5: DECLTR"
echo "---------------------------------------------"
RSCRIPT="${DECLTR_RSCRIPT:-Rscript}"
if "$RSCRIPT" -e 'suppressMessages({library(qs); library(segmented); library(yaml); library(data.table)})' >/dev/null 2>&1; then
    "$RSCRIPT" "$SUITE_DIR/DECLTR.r" \
        --manifest "$SCRIPT_DIR/decltr_manifest.tsv" \
        --config "$SCRIPT_DIR/decltr_config.yml" \
        --validate
    "$RSCRIPT" "$SUITE_DIR/DECLTR.r" \
        --manifest "$SCRIPT_DIR/decltr_manifest.tsv" \
        --config "$SCRIPT_DIR/decltr_config.yml"
    if [ ! -s "$OUT_DIR/decltr_test_labels.tsv" ] || [ "$(wc -l < "$OUT_DIR/decltr_test_labels.tsv")" -lt 2 ]; then
        echo "ERROR: DECLTR did not produce a labels table."
        exit 1
    fi
    echo "DECLTR completed successfully."
else
    echo "WARNING: R with qs/segmented/yaml/data.table not found; skipping DECLTR (set DECLTR_RSCRIPT)."
fi
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
