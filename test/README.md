# Test Dataset

Minimal test data for verifying installation of the LTRRT-Promoter-Suite pipeline.

## Data Description

All files are derived from a 1 Mb region of maize B73 chromosome 1 (original coordinates 30,500,000–31,500,000, offset to 1–1,000,001 in the test files):

| File | Description |
|------|-------------|
| `test_reads.fastq` | 212 ONT cDNA reads (untrimmed, raw); includes reads aligning to two structural LTR-RTs so IsoClassifier emits non-empty U3/LTR FASTA outputs |
| `test_reference.fa` | 1 Mb reference genome subset (B73 chr1) |
| `test_annotations.gff` | TE and gene annotations (658 features including 30 structural LTR-RTs, 25 genes) |
| `test_ltr_seqs.fa` | Full LTR sequences for 2 elements with expression data |
| `test_u3_seqs.fa` | U3 promoter region sequences for the same 2 elements |
| `test_tss_summary.tsv` | TSS positions for the 2 expressed LTR-RTs |

## Running the Test

```bash
conda activate LTRPromSuite_pipeline
cd LTRRT-Promoter-Suite
bash test/run_test.sh
```

The script runs all four tools in order:
1. **AccuMap** — PyChopper + Cutadapt + Minimap2 on the test reads
2. **IsoClassifier** — Read classification + U3/promoter sequence extraction against the test annotations and reference
3. **WindowScrubber** — Motif scanning on the pre-built LTR/U3 sequences
4. **Query_WSDB** — Example queries against the motif database

Output is written to `test/output/`. The test should complete in under 5 minutes.

## Cleaning Up

```bash
rm -rf test/output/
```

---

## Inspecting Published Results

`test/` also contains two result objects from the full LTR-RT analysis in maize. These are independent of the test run and require the DECLTR R environment.

**Activate the environment:**

```bash
conda activate DECLTR-env
R
```

**DECLTR results** (`DECLTR_results.qs`) — labeled output data frame from the DECLTR filtering pipeline:

```r
library(qs)
decltr <- qread("test/DECLTR_results.qs")
head(decltr)
colnames(decltr)
```

**WGCNA results** (`WGCNA_results.rds`) — named list with network and module data:

```r
wgcna <- readRDS("test/WGCNA_results.rds")
names(wgcna)
# $bwnet         — blockwiseModules output object
# $moduleColors  — per-gene module color assignments
# $mes           — module eigengenes (sample × module matrix)
# $traits        — trait/metadata matrix used for correlations
# $datExpr       — VST-normalized expression matrix (samples × genes)
# $softpwr       — soft-thresholding power selected for network construction
```

Run these commands from the `LTRRT-Promoter-Suite/` repo root, or adjust paths accordingly.
