# LTRRT-Promoter-Suite
[![DOI](https://zenodo.org/badge/1154798632.svg)](https://doi.org/10.5281/zenodo.19634858)

A suite of tools for characterizing the U3 promoter regions and transcription events of long terminal repeat retrotransposons (LTR-RTs) in maize.

---

## Table of Contents

- [Installation](#installation)
- [AccuMap](#accumap)
- [IsoClassifier](#isoclassifier)
- [WindowScrubber](#windowscrubber)
- [Query_WSDB](#query_wsdb)
- [DECLTR](#decltr)
- [Developer Notes](#developer-notes)

---

## Installation

Clone the repository and build the conda environments:

```bash
git clone https://github.com/cgooden-TE/LTRRT-Promoter-Suite
cd LTRRT-Promoter-Suite

# Environment for AccuMap, IsoClassifier, and WindowScrubber
mamba env create -f LTRPromSuite_pipeline.yml
mamba activate LTRPromSuite_pipeline

# Environment for DECLTR (and WGCNA)
mamba env create -f DECLTR_env.yml
mamba activate DECLTR-env
```

---

## AccuMap

Read directionalization, trimming, and alignment pipeline for long-read transcriptomics. Wraps PyChopper, Cutadapt, and Minimap2 into a single command, annotating output BAM files with strand-of-origin tags for downstream isoform analysis. Can be run for non-ONT reads by omitting the 'pyc' arguments.

### Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `--fq` | Yes | Input FASTQ of untrimmed, demultiplexed reads |
| `--sample` | Yes | Sample name prefix used for all output filenames |
| `--ref` | Yes (with `--run_map`) | Reference genome FASTA |
| `--run_pyc` | No | Run PyChopper for ONT primer removal |
| `--run_cut` | No | Run Cutadapt for homopolymer trimming |
| `--run_map` | No | Run Minimap2 splice-aware alignment |
| `--kit` | No | ONT sequencing kit (default: `PCB114`). Set to `none` for non-ONT data |
| `--map_preset` | No | Minimap2 preset: `splice` (default, ONT), `splice:hq` (PacBio), or `none` |
| `--pyc_threads` | No | PyChopper threads (default: 8) |
| `--cut_threads` | No | Cutadapt threads (default: 16) |
| `--map_threads` | No | Minimap2 threads (default: 24) |
| `--map_sec` | No | Allow secondary alignments (default: `no`) |
| `--map_gap` | No | Max intron length for Minimap2 (default: 5000) |

### Outputs

All files are prefixed with the `--sample` name:

| File | Description |
|------|-------------|
| `<sample>.pychopped.fastq` | Reads after PyChopper primer removal |
| `<sample>.cutadapt.fastq` | Reads after Cutadapt homopolymer trimming |
| `<sample>.minimap2.sorted.bam` | Sorted, indexed BAM from Minimap2 |
| `<sample>.strandtags.tsv` | Read name to PyChopper strand orientation mapping |
| `<sample>.STtagged.sorted.bam` | Sorted BAM annotated with PyChopper strand tags (ST) |
| `<sample>.STtagged.bed` | BED file with strand-resolved alignments |
| `<sample>.pychopper.log` | PyChopper log |
| `<sample>.pychopper.report.pdf` | PyChopper QC report |
| `<sample>.pychopper.rescued.fastq` | PyChopper rescued reads |
| `<sample>.cutadapt.log` | Cutadapt log |
| `<sample>.minimap2.log` | Minimap2 log |

### Example

```bash
# Full pipeline (ONT reads)
python3 AccuMap.py \
    --fq raw_reads.fastq \
    --sample CT_1 \
    --ref Zm-B73-REFERENCE-NAM-5.0.fa \
    --run_pyc \
    --run_cut \
    --run_map

# PacBio reads (skip PyChopper, use splice:hq preset)
python3 AccuMap.py \
    --fq hifi_reads.fastq \
    --sample PB_Ear_1 \
    --ref Zm-B73-REFERENCE-NAM-5.0.fa \
    --kit none \
    --map_preset splice:hq \
    --run_cut \
    --run_map
```

---

## IsoClassifier

Classifies reads overlapping structurally intact LTR retrotransposons into isoform categories (left-LTR(5'), right-LTR(3'), spanning(coding), read-through-5', read-through-3') and summarizes gene-level read counts and TSS positions. Uses interval trees to identify and exclude reads originating from nested transposable elements. Parallelizes classification by chromosome.

### Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `--gff` | Yes | GFF annotation with LTR retrotransposons (including left/right LTRs), genes, and nested features |
| `--bam` | Yes | One or more BAM files (space-separated) |
| `--output` | Yes | Output prefix for LTR isoform TSV and per-read files |
| `--tss_out` | Yes | Output path for LTR isoform TSS summary |
| `--gene_out` | Yes | Output path for gene read count and TSS summary |
| `--min_mapq` | No | Minimum mapping quality filter (default: 30) |
| `--threads` | No | Parallel processes for per-chromosome classification (default: 1) |
| `--genome-fasta` | No | Genome FASTA for U3/promoter sequence extraction (enables `u3_seq_extraction` after classification) |

### Outputs

All LTR output files use the `--output` prefix:

| File | Description |
|------|-------------|
| `<output>.tsv` | Per-element LTR isoform summary with read counts, percentages, mean lengths, splice counts, and unique junctions per category |
| `<tss_out>` | LTR TSS summary: top isoform, strand, total reads, and top 2 TSS positions per element |
| `<gene_out>` | Gene summary: total reads and top 2 TSS positions per gene |
| `<output>_primary_tss_density.tsv` | LTR primary TSS read density histogram |
| `<output>_secondary_tss_density.tsv` | LTR secondary TSS read density histogram |
| `<output>_gene_primary_density.tsv` | Gene primary TSS read density histogram |
| `<output>_gene_secondary_density.tsv` | Gene secondary TSS read density histogram |
| `<output>_ltr_3p_softclip_per_read_spliced.tsv` | Per-read 3' soft-clip data for spliced LTR reads |
| `<output>_ltr_3p_softclip_per_read_nonspliced.tsv` | Per-read 3' soft-clip data for non-spliced LTR reads |
| `<output>_ltr_exon_stats_per_read.tsv` | Per-read exon/intron structure for LTR-contained spliced reads |
| `<output>_gene_exon_stats_per_read.tsv` | Per-read exon/intron structure for gene spliced reads |
| `<output>_u3_seqs.fa` | U3 region sequences for expressed LTR elements (requires `--genome-fasta`) |
| `<output>_ltr_seqs.fa` | Full LTR sequences for expressed LTR elements (requires `--genome-fasta`) |
| `<output>_gene_2kb_proms.bed` | BED of ±1000 bp promoter windows around gene TSS (requires `--genome-fasta`) |
| `<output>_gene_2kb_proms.fa` | Strand-aware FASTA of gene promoter regions (requires `--genome-fasta`) |
| `<output>_gene_dummy_u3.fa` | Dummy U3 FASTA entries for upstream 1 kb of gene TSS (requires `--genome-fasta`) |

### Example

```bash
python3 IsoClassifier.py \
    --gff Zm-B73_TE_and_gene_annotations.gff \
    --bam CT_1.STtagged.sorted.bam CT_2.STtagged.sorted.bam CT_3.STtagged.sorted.bam \
    --min_mapq 30 \
    --threads 8 \
    --output CT_LTR_isoforms \
    --tss_out CT_LTR_TSS_summary.tsv \
    --gene_out CT_Gene_summary.tsv \
    --genome-fasta Zm-B73-REFERENCE-NAM-5.0.fa
```

---

## WindowScrubber

Motif scanner for LTR U3 promoter regions. Identifies core promoter elements (TATA box, CCAAT box, Y Patch, Inr, DPE, secondary TATA) and TA-rich regions using PWM log-odds scoring with FIMO-style p-values, and detects statistically significant CA dinucleotide runs. All hits are stored in a SQLite database for flexible downstream querying.

### Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `-l` / `--ltr-fasta` | Yes | FASTA of full LTR sequences (headers: `>Feature\|chr:start-end(strand)`) |
| `-u3` / `--u3-fasta` | Yes | FASTA of U3 region sequences (same header format) |
| `-t` / `--tss-summary` | Yes | TSS summary TSV from IsoClassifier (Feature, TSS1 columns) |
| `-db` / `--database` | Yes | Output SQLite database filename |
| `--pwm-json` | No | Custom PWM definitions in JSON format |
| `--tata_mismatch` | No | Allowed mismatches for TATA box (default: 0) |
| `--ccaat_mismatch` | No | Allowed mismatches for CCAAT box (default: 0) |
| `--ypatch_mismatch` | No | Allowed mismatches for Y Patch (default: 0) |
| `--inr_mismatch` | No | Allowed mismatches for Inr (default: 0) |
| `--window-size` | No | Sliding window size for TA-rich detection (default: 10) |
| `--ta-threshold` | No | TA fraction threshold for TA-rich windows (default: 0.75) |
| `--sig-alpha` | No | Per-motif significance level (default: 1e-4) |
| `--ca-min-length` | No | Minimum CA dinucleotide run length (default: 10) |
| `--ca-alpha` | No | Significance level for CA runs (default: 0.05) |
| `--tata-max-dist` | No | Max upstream distance for primary TATA from TSS (default: 100) |
| `--inr-max-dist` | No | Max offset around TSS for Inr motif (default: 10) |
| `--bg-{A,C,G,T}` | No | Background nucleotide frequencies (default: 0.25 each) |
| `--bg-n` | No | Background k-mers for null distribution (default: 200000) |

### Outputs

| Output | Description |
|--------|-------------|
| `<database>.db` | SQLite database containing five tables: |
| &emsp; `elements` | Feature metadata, coordinates, TSS position, and U3/LTR sequence lengths |
| &emsp; `ta_regions` | TA-rich regions per element (relative and absolute coordinates) |
| &emsp; `ca_runs` | All significant CA dinucleotide runs with p-values |
| &emsp; `motif_hits` | All motif matches (TATA, CCAAT, Y Patch, Inr, DPE, secondary TATA) with scores, p-values, and distances to TSS |
| &emsp; `thresholds` | Per-motif score thresholds and significance cutoffs |

### Example

```bash
python3 WindowScrubber.py \
    -l CT_full_ltr_sequences.fa \
    -u3 CT_u3_regions.fa \
    -t CT_LTR_TSS_summary.tsv \
    -db CT_motif_hits.db \
    --tata_mismatch 1 \
    --ccaat_mismatch 0 \
    --sig-alpha 1e-4
```

---

## Query_WSDB

Query interface for the WindowScrubber SQLite database. Retrieves, filters, and summarizes motif hits, CA dinucleotide runs, and promoter element composition across LTR retrotransposons. Supports flexible ranking, distance-based filtering, threshold enforcement, and multi-motif presence queries.

### Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `-db` / `--database` | Yes | SQLite database file produced by WindowScrubber |

**Query mode** (mutually exclusive, one required):

| Flag | Description |
|------|-------------|
| `--best-per-element MOTIF` | Best hit per element for a specific motif type |
| `--best-all-motifs` | Best hit for every motif type per element (wide format, one row per element) |
| `--all-hits` | All hits with optional filtering |
| `--require-motifs MOTIF1,MOTIF2,...` | Elements containing all specified motif types |
| `--stats` | Print database summary statistics to stdout |
| `--ca-summary` | Export CA dinucleotide run data |

**Filtering options:**

| Flag | Description |
|------|-------------|
| `--motif-type` | Filter by motif type (TATA, CCAAT, YPATCH, INR7, DPE7, SEC_TATA) |
| `--min-score` / `--max-score` | PWM score range |
| `--min-dist` / `--max-dist` | Absolute distance-to-TSS range |
| `--in-ta-region` | Only hits within TA-rich regions |
| `--features` | Comma-separated list of feature IDs to include |
| `--p-max` | Maximum p-value threshold |
| `--range` | Single distance range as `lo:hi` (relative to TSS) |
| `--ranges` | Per-motif distance ranges (e.g., `TATA:-100:0,CCAAT:-460:-140,INR7:-10:10`) |
| `--apply-threshold` | Require score >= stored threshold from the database |
| `--threshold-method` | Threshold method to apply (default: `empirical_p`) |
| `--threshold-param` | Specific parameter string (e.g., `alpha=1e-4;bg=iid;n=200000`) |
| `--order-by` | Column to rank by for `--best-per-element`: `score`, `dist_to_tss`, or `y_count` (default: `score`) |
| `--ascending` | Pick lowest value instead of highest when ranking |
| `--ca-min-length` | Minimum CA run length |
| `--ca-all` | Include non-significant CA runs |
| `-o` / `--output` | Output TSV path. If omitted, prints first 10 results to stdout |

### Outputs

All query modes write tab-separated (TSV) files (or print to stdout when `-o` is omitted):

| Query Mode | Output Description |
|------------|--------------------|
| `--best-per-element` | One row per element with the top-ranked hit for the requested motif (feature, motif_type, coordinates, dist_to_tss, score, p_value, sequence, in_ta_region, y_count) |
| `--best-all-motifs` | One row per element in wide format with columns for each motif type (`{MOTIF}_start_abs`, `{MOTIF}_end_abs`, `Dist_TSS_to_{MOTIF}`; `-1` if absent) |
| `--all-hits` | One row per hit with full metadata (feature, motif_type, coordinates, score, p_value, sequence, chrom, strand, tss_abs) |
| `--require-motifs` | One row per qualifying element with a `found_motifs` column listing all motifs present |
| `--ca-summary` | One row per CA run (feature, coordinates, length, p_value, is_significant, chrom, strand) |
| `--stats` | Summary printed to stdout: element count, per-motif hit counts, score statistics, TA-region counts, and CA run counts |

### Examples

```bash
# Database summary statistics
python3 Query_WSDB.py -db CT_motif_hits.db --stats

# Best TATA hit per element ranked by score
python3 Query_WSDB.py -db CT_motif_hits.db \
    --best-per-element TATA \
    --order-by score \
    -o CT_best_tata.tsv

# Best hit for every motif type per element (wide format)
python3 Query_WSDB.py -db CT_motif_hits.db \
    --best-all-motifs \
    -o CT_best_all_motifs.tsv

# All TATA hits within TA-rich regions
python3 Query_WSDB.py -db CT_motif_hits.db \
    --all-hits \
    --motif-type TATA \
    --in-ta-region \
    -o CT_tata_in_ta.tsv

# All hits with motif-specific distance ranges and threshold enforcement
python3 Query_WSDB.py -db CT_motif_hits.db \
    --all-hits \
    --ranges "TATA:-100:0,CCAAT:-460:-140,INR7:-10:10,DPE7:0:50" \
    --apply-threshold \
    -o CT_filtered_hits.tsv

# Elements that contain TATA, CCAAT, and INR7 motifs
python3 Query_WSDB.py -db CT_motif_hits.db \
    --require-motifs TATA,CCAAT,INR7 \
    -o CT_complete_promoters.tsv

# Export significant CA dinucleotide runs
python3 Query_WSDB.py -db CT_motif_hits.db \
    --ca-summary \
    -o CT_ca_runs.tsv
```

---

## DECLTR

R-based multi-omic integration and activity classification pipeline for LTR retrotransposons and genes. Merges expression data from Illumina, PacBio, and ONT platforms with ChIP-seq (H3K4me3, H3K27ac), DNA methylation (UMR), CAGE, and promoter motif annotations. Uses segmented regression to estimate per-platform expression thresholds, applies sigmoid-transformed scoring across tissue groups, and assigns activity labels (Constitutive, Facultative, Tissue-Specific, Developmental, Vegetative, or Inactive).

> **Note:** DECLTR is provided for transparency. The current implementation is tuned to specific data types and EDTA naming conventions used in this analysis. Reuse will require adapting the hardcoded file paths and column-naming logic for your data.

### Inputs

DECLTR reads from directory structures rather than command-line arguments. The following paths are set at the top of the script and must be edited before running:

| Variable | Description |
|----------|-------------|
| `data_dir` | Directory of ChIP-seq and UMR intersect GFF files (bedtools intersect outputs) |
| `isoform_dir` | Directory of IsoClassifier output TSVs (isoform counts and TSS summaries) |
| `motif_dir` | Directory of WindowScrubber/Query_WSDB motif summary TSVs |
| `combined_ref` | GFF annotation file with LTR retrotransposons and genes |
| `illumina_path` | TSV of Illumina short-read expression counts per feature |

### Outputs

| File | Description |
|------|-------------|
| `Segmented_Breakpoints.csv` | Per-platform (Illumina, PacBio, ONT) segmented-regression expression thresholds on log1p scale |
| `*.qs` | Serialized R object (qs format) of the fully labeled data frame with all merged annotations, activity scores, and classification labels |

### Example

```bash
# Edit paths at top of DECLTR.r, then:
Rscript DECLTR.r
```

---

## Developer Notes

**AccuMap:** Developed primarily for ONT long reads. Requires that reads have not been demultiplexed or trimmed, only base-called, for optimal usage with PyChopper. For non-ONT data, start at the `--run_cut` step as long as reads still contain homopolymers (poly(A) or poly(T)).

**IsoClassifier:** Best used with a GFF containing annotations for all transposons and genes to enable nested-element removal. The GFF *must* contain records for the left and right long terminal repeats (lLTR and rLTR as named by EDTA) of each LTR-RT, in addition to the full `LTR_retrotransposon` annotation for each structurally intact locus. When `--genome-fasta` is provided, IsoClassifier runs U3/promoter extraction internally (superseding the standalone `U3_Seq_Extractor.py`).

**WindowScrubber:** Uses IsoClassifier outputs to define search windows for each motif context. External TSS data (e.g., from CAGE) can substitute for IsoClassifier TSS calls, but the U3/LTR FASTA files (generated by IsoClassifier with `--genome-fasta`) are still required and the TSS file must match IsoClassifier's column format (`Feature`, `TSS1`).

**DECLTR:** Provided for development transparency. The methods are tuned to the data types and EDTA naming conventions used in this analysis. Reuse requires re-tooling the hardcoded paths and parsing logic for your data until the generalized version is released.

**SpliceJunTest:** This is a file written for testing canonical splice junctions in our data. It is published here since it is mentioned in the associated manuscript but should not be considered part of the official LTRRT Promoter Suite tools. 
