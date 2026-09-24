# LTRRT-Promoter-Suite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19634858.svg)](https://doi.org/10.5281/zenodo.19634858)

Tools developed to characterize the U3 promoter regions and transcription events of long terminal repeat retrotransposons.

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

Classifies long reads into isoform categories at every annotated transposable element and gene, and
summarizes read counts, TSS positions, and 3'-end (cleavage) positions. Reads are assigned by
5' end: of every feature whose span contains a read's TSS, the smallest one wins, so a read
starting inside a TE nested in a gene is credited to the TE and flagged as nested in the gene.
Sense and antisense transcripts are both retained and written to separate isoform tables.
Parallelizes classification over 1 Mb genomic chunks.

### Element classes

| Class | Source rows | TSS assigned from | Isoform vocabulary |
|-------|-------------|-------------------|--------------------|
| `LTR_structural` | `*LTR_retrotransposon` with `Method=structural` and both LTRs resolved | the LTR the read initiates in | `ltr5_contained`, `ltr3_contained`, `spanning`, `readout_5ltr`, `readout_3ltr`, `readout_internal`, `spliced_ltr5`, `spliced_ltr3`, `spliced_spanning`, `partial` |
| `LTR_fragment` | `*LTR_retrotransposon` with `Method=homology`, plus structural elements missing an LTR | read coordinate in the element body | `full_length`, `spliced`, `readout`, `partial` |
| `TIR` | `*_TIR_transposon` (CACTA, hAT, Mutator, PIF/Harbinger, Tc1/Mariner) | element body | same |
| `Helitron` | `helitron` | element body | same |
| `LINE`, `SINE` | `*LINE*`, `*SINE*` | element body | same |
| `Gene` | `gene` | element body, with full-length judged on the exon union | same |

Annotation bookkeeping (`repeat_region`, `target_site_duplication`) and non-transposon repeats
(`knob`, `centromeric_repeat`, `subtelomere`, `rDNA_intergenic_spacer_element`, `low_complexity`)
are never loaded as elements. Use `--te-classes` to narrow the set further.

Categories are named in the **read's** frame: the LTR a read initiates in is its 5' LTR whether or
not the read agrees with the element's annotated strand. That is the whole of the sense/antisense
inversion — an antisense read on a `+` element is evaluated against the element's right LTR.

One consequence, since it is the only way the category arises: a read that starts in the element's
annotated **3'** LTR and ends in its annotated **5'** LTR is not a backwards transcript, it is an
antisense one. In its own frame it ran 5' → 3', so it is scored `spanning` and appears in the
antisense table. That traversal is geometrically impossible for a sense read.

The three spliced categories are kept distinct because they are different transcripts:

| Category | Meaning |
|---|---|
| `spliced_ltr5`, `spliced_ltr3` | Contained within the LTR the read initiated in, carrying a real intron. Typically the dominant phenotype — so `ltr5_contained` / `ltr3_contained` mean contained *and unspliced* |
| `spliced_spanning` | Spans both LTRs but skipped most of the internal domain (coverage below `--spliced-cov-frac`) |

A spanning read always initiates in its 5' LTR, so folding it into `spliced_ltr5` would make the
pair asymmetric — the demotion could never yield `spliced_ltr3`, and the contained phenotype could
not be counted on its own. Giving it `spliced_spanning` keeps both readable.

Read-out is likewise split by where the transcript started, so the three are never mixed in one
column: `readout_5ltr` and `readout_3ltr` initiate in an LTR, `readout_internal` initiates in the
internal domain (not an LTR-driven transcript at all). Non-LTR classes use plain `readout`.

### Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `--gff` | Yes | GFF annotation with TE features (EDTA-style), `long_terminal_repeat` records, and genes |
| `--bam` | Yes* | One or more BAM files (space-separated). *Not used with `--merge-sheet` |
| `--output` | Yes | Output prefix for isoform tables, summaries, and per-read files |
| `--tss_out` | Yes | Output path for the combined isoform TSS summary |
| `--gene_out` | Yes | Output path for the gene read count and TSS summary |
| `--gene-gff` | No | GFF3 carrying `exon` records (e.g. the Ensembl gene annotation). Exons are keyed to genes by **gene ID**, so seqid naming need not match `--gff`, but coordinates must be on the same assembly. Without it, gene `full_length` falls back to gene-span coverage, which almost never fires for spliced transcripts (the median maize exon union is ~54% of the gene span) |
| `--canonical-exons-only` | No | Build each gene's exon union from the `Ensembl_canonical` transcript only (default: union across all transcripts) |
| `--te-classes` | No | Comma-separated classes to classify (default: all six above) |
| `--min_mapq` | No | Minimum mapping quality filter (default: 30) |
| `--threads` | No | Parallel worker processes (default: 1) |
| `--min-intron-len` | No | Minimum CIGAR `N` length counted as a real intron (default: 69) |
| `--full-length-frac` | No | Covered fraction at or above which a contained read is `full_length` (default: 0.8) |
| `--spliced-cov-frac` | No | Covered fraction below which a spliced read is `spliced` rather than `full_length`/`spanning` (default: 0.5) |
| `--tss-window` | No | Width in bp of the TSS agreement window, centred on the modal 5' end (default: 3) |
| `--tss-min-frac` | No | Fraction of a feature's reads that must fall in that window for `TSS_Called=1` (default: 0.5). **Reported only — never discards reads** |
| `--progress-every` | No | Print one progress line per N reads assessed, aggregated across workers, with percent complete and an ETA derived from the BAM indexes (default: 1000000) |
| `--include-partial` | No | Include the `partial` category in the isoform tables. Partial reads are always written to the per-read tables regardless |
| `--trust-st-tag` | No | Trust the PyChopper/AccuMap `ST` tag as a genomic strand call (off by default; see the `infer_read_strand` docstring) |
| `--genome-fasta` | No | Genome FASTA for U3/promoter sequence extraction |
| `--assign` | No | `innermost` (default): a read whose 5' end lies in several features goes to the smallest. `em`: the read is shared among those features by expectation-maximization using TSS peaks, 3' ends, read length and splicing, and credited to the most probable one. Recommended for nested TEs and TEs in genes. |
| `--em-tss-halfwidth` | No | Half-width (bp) used to call TSS peaks and pool 5' ends (default 8; about 3 is reasonable for PacBio HiFi) |
| `--em-tss-trunc-ratio` | No | A TSS peak must hold at least this many times the 5' ends that read-through transcription plus RT drop-off would leave in its window (default 20; 0 disables the test). Stops an RT stall site inside a transcribed feature from being sold as a TSS to whichever unit happens to have that position in its promoter zone |
| `--em-eps-evidence-k` | No | How much evidence (certain reads + called TSS peaks) a feature needs before it may claim reads that started at a TSS it has no peak for (default 5; 0 restores a flat allowance for every feature regardless of evidence) |
| `--em-alpha-prior` | No | Dirichlet prior on a feature's abundance, so reads it has already absorbed do not make it more attractive to the next read (default 0.5) |
| `--em-trunc-min-sample` | No | Reads a BAM needs before it gets its own read-length survival curve instead of the pooled one (default 2000). Libraries that differ in RT processivity should not share a curve |
| `--em-exclusive-nested-promoters` | No | Give a nested feature exclusive claim to TSS peaks inside it even where it overlaps the host's own promoter (a gene's 5' end, an element's LTRs). Off by default: a short annotation across a gene's TSS would otherwise take that gene's promoter outright, leaving the gene with none. Features nested in a host's *body* always keep exclusive promoters |
| `--em-antisense-prior` | No | Prior odds that a read arose from antisense transcription rather than sense at the same locus (default 0.25; 1.0 removes the prior). Antisense is real but rarer, and nothing else in the model says so |
| `--em-junc-init-frac` | No | Share of initiation the model must give a read before it may use a feature's learned read-out junctions (default 0.5). Lower retains more read-out; higher is stricter about truncated reads |
| `--em-readout-bp` | No | Length of the read-out region past a feature's 3' end (default 10000) |
| `--em-gene-promoter-bp` | No | Window at a gene's annotated 5' end treated as promoter, in addition to the first exon (default 300) |
| `--em-max-iter`, `--em-tol` | No | EM stopping rules (default 500 iterations; stop when no feature's expected count moves by more than 0.001 reads) |
| `--tss-method` | No | `consensus` (default): library-balanced TSS call (below). `mode`: the original modal 5' end of all pooled reads, to reproduce results from earlier versions |
| `--platform` | No | Platform of every `--bam` not in `--library-sheet`: `ont` (default), `pacbio_ccs` or `pacbio_clr`. Sets the TSS half-width |
| `--tss-halfwidths` | No | Platform half-widths in bp (default `ont=8,pacbio_ccs=3,pacbio_clr=8`); add new platforms the same way |
| `--tss-lib-k` | No | Library weight `n/(n+k)` for a feature with `n` reads in that library (default 10: a library with 10 reads for the feature counts half) |
| `--tss-min-lib-reads` | No | Reads a library needs for a feature to count toward `Libraries_Qualifying`/`Samples_Supporting` (default 5) |
| `--tss-evidence` | No | External TSS evidence as `name=path[,path]`, repeatable, **highest priority first** (e.g. `smar2c2=...` then `cage=shoot.gff3,root.gff3`). Formats by extension: `.gff`/`.gff3` (CAGE clusters, dominant TSS in column 8, score in column 6), `.bed`, or a TSS table with `seq`/`TSS`/`strand`/`nTAGs` columns (TSRexplorer TSS set) |
| `--tss-evidence-window` | No | Evidence counted within ± this many bp of a candidate site (default 10) |
| `--tss-evidence-min` | No | Summed evidence (tags or CAGE score) a track needs at a candidate to choose TSS1 (default 3) |
| `--tss-evidence-min-reads` | No | Reads that must start at a position before evidence may make it TSS1 (default 2) |
| `--library-sheet` | No | TSV with `Source` (BAM path or basename), `Library`, `Platform`. BAMs sharing a `Library` (runs, SMRT cells, barcodes of one biosample) are summed before the TSS call |
| `--write-state` | No | Also write `<output>_state_{features,5p,3p}.tsv.gz`, the per-sample counts that `--merge-sheet` combines |
| `--merge-sheet` | No | Merge mode: TSV with `Source` (a `--write-state` run's `--output` prefix) and optional `Library`/`Platform` overrides. Reads no BAMs |
| `--merge-allow-mismatch` | No | Merge states whose annotation or counting arguments differ, with a warning instead of an error |

### Outputs

| File | Description |
|------|-------------|
| `<output>_sense.isoforms.tsv` | Per-feature isoform table for reads matching the element strand |
| `<output>_antisense.isoforms.tsv` | Same, for reads on the opposite strand |
| `<tss_out>` | Isoform TSS summary, all orientations, with an `Orientation` column |
| `<output>_10site.tss_summary.tsv` | As above with the top 10 TSS positions |
| `<output>_{sense,antisense}.tss_summary.tsv` | Single-orientation copies, with `Chrom`. **Feed these to WindowScrubber**, which keys its TSS lookup on `Feature` and needs one row per feature |
| `<output>_10site.cleavage_summary.tsv` | 3'-end distribution per feature and orientation |
| `<gene_out>`, `<output>_10site.gene_summary.tsv` | Gene read counts and TSS positions, with an `Orientation` column |
| `<output>_{primary,secondary}_tss_density.tsv` | TE TSS read-density histograms (read-frame distances) |
| `<output>_gene_{primary,secondary}_density.tsv` | Gene TSS read-density histograms |
| `<output>_3p_softclip_per_read_{spliced,nonspliced}.tsv` | Per-read 3' soft-clip data, with `Class`, `Orientation`, `Nested_In`, `Terminates_In`. The read's 5' end is `TSS`; `Aln_Reverse` is the BAM reverse flag, distinct from `Orientation` |
| `<output>_te_exon_stats_per_read.tsv` | Per-read exon/intron structure for spliced TE reads |
| `<output>_gene_exon_stats_per_read.tsv` | Per-read exon/intron structure for spliced gene reads |
| `<output>_u3_seqs_{sense,antisense}.fa` | U3 region sequences, structural LTR-RTs only (requires `--genome-fasta`) |
| `<output>_ltr_seqs_{sense,antisense}.fa` | Full LTR sequences, structural LTR-RTs only |
| `<output>_gene_2kb_proms.{bed,fa}` | ±1000 bp promoter windows around gene TSS (sense only). Deprecated: WindowScrubber now cuts gene windows from the genome |
| `<output>_gene_dummy_u3.fa` | Dummy U3 entries for the upstream 1 kb of each gene TSS. Deprecated, no longer used by WindowScrubber |
| `<output>_em_reads.tsv` | (`--assign em`) Every read whose 5' end lies in more than one feature: innermost feature, assigned feature, posterior, probability that the 5' end is a genuine TSS, and all candidates with their posteriors |
| `<output>_em_units.tsv` | (`--assign em`) Per feature/orientation in contested loci: certain reads, reads the innermost rule would give, reads credited by EM (most probable), expected reads (fractional; use for quantification), expected TSS-initiated reads, called TSS peaks, EM TSS peak and poly(A) peak, 3'-end, splicing and intronic fractions |
| `<output>_em_read_survival.tsv` | (`--assign em`) Probability of a read extending a given length from its 3' end (Kaplan–Meier), i.e. RT processivity and degradation. With several BAMs it gains a `Sample` column and carries one curve per BAM plus the pooled one |
| `<output>_tss_consensus.tsv` | Per feature/orientation: `TSS1`/`TSS2` with their scores, `Libraries_Qualifying`, `Samples_Supporting` (libraries whose own peak lies within their half-width of `TSS1`), `Support_Frac`, `Low_Agreement` (fewer than half agree; a candidate tissue-specific promoter), `TSS_Source` (`reads`, or the evidence track that chose `TSS1`), `Evidence_Support`, and `Reads_TSS1` (what the reads alone called) |
| `<output>_tss_by_library.tsv` | Each library's reads, weight and own 5'-end peak per feature/orientation, and whether that peak agrees with `TSS1` |
| `<output>_state_{features,5p,3p}.tsv.gz` | (`--write-state`) Category counts and junctions, 5' ends per sample, 3' ends; the input to `--merge-sheet` |

Isoform-table columns beyond the per-category counts:

| Column | Meaning |
|--------|---------|
| `Class` | Element class from the table above |
| `Orientation` | `sense` or `antisense`, read strand versus element strand |
| `Strand`, `Strand_Source` | Element strand, and whether it was `annotated` or `inferred` from read support (features annotated `.` or `?`) |
| `Nested_In`, `Nested_In_Class` | Features this one sits inside, and their classes (`NA` when not nested) |
| `TSS_Called`, `TSS_Peak`, `TSS_Peak_Frac` | The TSS call (`TSS1`, below), the library-weighted fraction of the feature's reads within `--tss-window` of it, and whether that clears `--tss-min-frac` |

Per category, four columns at most: `<cat>_reads`, `mean_len_<cat>`, `n_spliced_<cat>`,
`unique_juncts_<cat>`. The last two are omitted where the category's own membership test already
fixes them, rather than written out as constants:

| Categories | Omitted | Why |
|---|---|---|
| `spliced_ltr5`, `spliced_ltr3`, `spliced_spanning`, `spliced` | `n_spliced_` | membership requires an intron, so it always equals `<cat>_reads` |
| `ltr5_contained`, `ltr3_contained` | `n_spliced_`, `unique_juncts_` | membership requires *no* intron, so both are always 0 |

The counter prefix is `n_spliced_`, not `spliced_`, because several categories are themselves
named `spliced_*` — `spliced_spanning` as a counter prefix collided with the `spliced_spanning`
category's own columns.

Columns that were pure restatements of others are not written. Recover them with:

| Wanted | From |
|---|---|
| `<cat>_pct` | `<cat>_reads / Total_Reads * 100` |
| `TSS_Peak_Count` | `TSS_Peak_Frac * Total_Reads` |
| `Nested` | `Nested_In != "NA"` |

`TSS1` in the TSS summary and `TSS_Peak` in the isoform table are the same quantity — both pool 5'
ends across every category. (`TSS1` previously ranked only the top category, which made the two
files disagree.)

**How the TSS is called.** One call per feature/orientation feeds every output that reports or
uses a TSS (isoform tables, TSS and gene summaries, densities, U3/LTR and promoter FASTAs), so they
always agree. For each library `s` with `n_s` reads for the feature:

- `c̃_s(x)` = the library's 5' ends within ±`h_s` of `x`, where `h_s` is its platform's half-width
  (ONT 8 bp, PacBio CCS 3 bp, PacBio CLR 8 bp);
- `p_s(x) = c̃_s(x) / n_s`, the share of the feature's reads starting near `x`, so library depth
  cancels;
- `w_s = n_s / (n_s + k)`, so a library with few reads for the feature counts little and a deep
  library counts no more than any other well-covered one;
- `Score(x) = Σ w_s p_s(x) / Σ w_s`, evaluated at observed 5' ends.

`TSS1` is the top-scoring position (ties: more raw reads, then more upstream). `TSS2`… are the next
best, each outside ±max(`h_s`) of those already chosen, so `TSS2` is a second promoter rather than
`TSS1`±1. `Count{i}` stays the raw number of reads starting exactly at `TSS{i}`. With one library
this is a smoothed modal 5' end; `--tss-method mode` restores the unsmoothed pooled mode.

**External TSS evidence.** Long-read 5' ends usually include a pile at the true TSS but cannot say
which pile it is: 5' truncation and RT stalls build piles too. With `--tss-evidence`, tracks are
tried in priority order; the first with at least `--tss-evidence-min` support near any position where
`--tss-evidence-min-reads` reads start picks, among those positions, the one it supports most. With
no track qualifying, the reads-only call stands. `TSS1` therefore always sits where reads in these
libraries start — evidence chooses among those sites, it never adds one — and `TSS_Source` says what
decided it.

On chr1 of one ONT cDNA library (`OG_ONT_CT1`), scored on genes with Smar2C2 support. Smar2C2 tags
were split at random into two halves, one used as evidence and the other held out as truth, and CAGE
is scored independently wherever it was not itself the evidence:

| TSS1 from | within 10 bp of held-out Smar2C2 TSS | within 10 bp of CAGE dominant TSS |
|---|---|---|
| reads, `--tss-method mode` | 37.5% | 30.9% |
| reads, consensus | 38.2% | 31.9% |
| reads + Smar2C2 | **52.4%** | **39.8%** |
| reads + CAGE | **42.3%** | (not independent) |

With one library, consensus and mode are statistically indistinguishable; the gain comes from the
evidence. Both external sets are shoot (CAGE also root), so a tissue-specific TSS in other tissues
can be passed over for a shoot one; `TSS_Source` and `Reads_TSS1` make every such choice visible.

```bash
python3 IsoClassifier.py ... \
    --tss-evidence smar2c2=B73_ForAdap_TSSset_NoScaf.txt \
    --tss-evidence cage=dTSS_shoot_B73-AGPv5.gff3,dTSS_root_B73-AGPv5.gff3
```

### Running samples separately and merging

Holding every BAM in one run is memory-bound. Run each sample on its own with `--write-state`, then
merge. The merge reads only the state files and writes the same summary files as a normal run:

```bash
# 1. one run per BAM
python3 IsoClassifier.py --gff TEs_LTRs_Genes_UpdStr.gff --gene-gff genes.gff3 \
    --bam S1.STtagged.sorted.bam --platform ont --write-state --assign em \
    --output per_sample/S1 --tss_out per_sample/S1.tss.tsv --gene_out per_sample/S1.gene.tsv

# 2. merge sheet (tab-separated); runs of one biosample share a Library
#    Source            Library      Platform
#    per_sample/S1     S1           ont
#    per_sample/R1     biosampleA   pacbio_clr
#    per_sample/R2     biosampleA   pacbio_clr

# 3. merge
python3 IsoClassifier.py --gff TEs_LTRs_Genes_UpdStr.gff --gene-gff genes.gff3 \
    --merge-sheet merge_sheet.tsv --genome-fasta genome.fa \
    --output merged/all --tss_out merged/all.tss.tsv --gene_out merged/all.gene.tsv
```

- The merge stops if the states were made with a different annotation or different counting
  arguments (`--assign`, `--min_mapq`, `--te-classes`, ...), since their counts would not be the
  same quantity.
- Unknown-strand features are resolved on the summed support of all samples.
- Each sample's reads are counted as that sample assigned them; EM is not re-run, and per-read and
  EM tables stay with the per-sample runs.
- Merging one sample reproduces that sample's own summary files exactly, and merging two halves of
  a BAM reproduces one run over both (`test/tss_consensus/run_tss_consensus_test.sh`).

`Terminates_In` in the per-read tables names the innermost feature containing a read's 3' end when
it differs from the assigned feature. Together with `Nested_In` it isolates chimerism candidates:
transcripts that start in a gene and terminate inside a nested TE, or start in a TE and run out
through a gene.

`_em_units.tsv` also carries `Promoter_Anchor` (`observed` when a called TSS peak lies in the
feature's promoter zone, `annotated` when the zone rests on the annotated 5' end alone) and
`TSS_Offset` (bp from the annotated 5' end to the EM TSS peak, positive downstream). On chr1 of an
ONT cDNA library 337 of 1012 gene sense units were `annotated`, and 62 genes had the observed peak
more than 500 bp from the annotation -- those promoter zones, and anything built on them, are the
least trustworthy rows in the table.

Per-read tables gain `Assign_Method` (`unique`/`innermost`/`em`), `Assign_Posterior`,
`TSS_Init_Prob`, `Innermost_Feature` and `EM_Candidates`.

Interpreting `--assign em`:
- Use `EM_Reads` (fractional) for expression and DE, and the per-read assignment for read-level
  work.
- Posteriors between about 0.4 and 0.6 mean the read cannot be resolved.
- EM can only choose among annotated features. A read whose 5' end falls outside every feature it
  truly belongs to cannot be recovered by any 5'-end-based rule. An example is a TE-initiated
  transcript whose read lost its 5' end and begins in a downstream gene exon.
- `test/nested_em/preflight_5p_truncation.py` measures how often each library loses 5' ends, on
  genes that overlap no other feature; run it with the same `--gff`/`--gene-gff`/`--bam` to judge
  whether `--assign em` matters for your data.
- Treat `antisense` assignments in contested loci with suspicion. Where a read's splice junctions
  settle the question independently (`test/nested_em/junction_audit.py`), antisense calls on
  chr1 of an ONT cDNA library were almost always wrong: the reads were 5'-truncated sense reads of
  an overlapping gene. `--em-tss-trunc-ratio`, `--em-eps-evidence-k` and `--em-alpha-prior`
  between them cut those calls by more than half, but they have not eliminated them, so a unit
  whose `Unique_Reads` is 0 or 1 and whose `EM_Reads` is large still deserves a look in a browser
  before it is believed.

Auditing a run against splice junctions:

```bash
python3 test/nested_em/junction_audit.py \
    --em-reads out_em_reads.tsv --bam sample.STtagged.sorted.bam \
    --gene-gff genes_with_exons.gff3 --out junction_audit.tsv
```

A read whose junctions match the annotated intron chain of exactly one candidate, and of no other,
was transcribed from that candidate whatever its 5' end suggests. The EM never looks at junction
identity, so this is independent evidence; the script scores `--assign em` and the innermost rule
on the same reads.

### Example

```bash
python3 IsoClassifier.py \
    --gff TEs_LTRs_Genes_UpdStr.gff \
    --gene-gff Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.gff3 \
    --bam CT_1.STtagged.sorted.bam CT_2.STtagged.sorted.bam CT_3.STtagged.sorted.bam \
    --min_mapq 30 \
    --threads 8 \
    --output CT_isoforms \
    --tss_out CT_TSS_summary.tsv \
    --gene_out CT_Gene_summary.tsv \
    --genome-fasta B73.PLATINUM.pseudomolecules-v1.fasta
```

Loading the full maize TE + gene annotation (~1.5 M GFF rows) yields ~1.22 M elements in ~1 minute
at ~2.3 GB resident, shared across workers by `fork`.

---

## WindowScrubber

TSS-anchored promoter motif scanner for genes and every TE class. Identifies core promoter elements (TATA box, CCAAT box, Y Patch, Inr, DPE, secondary TATA) and TA-rich regions using PWM log-odds scoring with FIMO-style p-values, and detects statistically significant CA dinucleotide runs. All hits are stored in a SQLite database for flexible downstream querying.

Every element is anchored on IsoClassifier's `TSS1` and scanned with the same TSS-relative motif bands. The element class decides which of two schemas sets the scanned sequence and the TA-rich search zone:

| Schema (`anchor_mode`) | Elements | Scanned sequence | TA-rich zone |
|---|---|---|---|
| `U3` | `LTR_structural` whose TSS lies in the TSS-bearing LTR (`-l` given) | That LTR, in the transcript's frame | U3: upstream LTR edge to TSS |
| `TSS_window` | `Gene`, `TIR`, `Helitron`, `LINE`, `SINE`, `LTR_fragment`, and structural LTRs without a usable U3 | Genome, TSS − `--flank-up` to TSS + `--flank-down` | Fixed core band, TSS − `--core-up` to TSS + `--core-down` |

No U3 is imputed for non-LTR elements. U3 lengths vary while the core band is fixed, so compare `in_ta_region` / TA-rich rates within one `anchor_mode` only. All sequence is fetched from the genome; the LTR FASTA only supplies coordinates.

### Inputs

| Flag | Required | Description |
|------|----------|-------------|
| `-t` / `--tss-summary` | Yes | IsoClassifier `<output>_{sense,antisense}.tss_summary.tsv` (`Feature`, `Class`, `Orientation`, `Strand`, `Chrom`, 0-based `TSS1`) |
| `-g` / `--genome-fasta` | No | faidx-indexed genome (default: B73v5 PLATINUM) |
| `-l` / `--ltr-fasta` | No | IsoClassifier `<output>_ltr_seqs_<orient>.fa`; enables the U3 schema. Without it every element uses `TSS_window` |
| `-i` / `--isoforms` | No | IsoClassifier `<output>_<orient>.isoforms.tsv`; supplies `Chrom` for TSS summaries written before that column existed |
| `-u3` / `--u3-fasta` | No | Deprecated and ignored; U3 is derived from the LTR bounds and `TSS1` |
| `--orientation` | No | `sense` (default) or `antisense` rows of the TSS summary |
| `--classes` | No | Comma-separated classes to scan (default: all) |
| `--flank-up` / `--flank-down` | No | `TSS_window` span around the TSS (default: 500 / 200) |
| `--core-up` / `--core-down` | No | `TSS_window` TA-rich zone around the TSS (default: 100 / 0) |
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
| `<database>.db` | SQLite database containing six tables: |
| &emsp; `elements` | Feature, `element_class`, `orientation`, `anchor_mode`, scanned span (`seq_start`/`seq_end`), `u3_start`/`u3_end` (NULL outside U3 mode), TA-search zone, TSS, read support, `seq_truncated` (window clipped at a contig end) |
| &emsp; `ta_regions` | TA-rich regions per element (relative and absolute coordinates) |
| &emsp; `ca_runs` | All significant CA dinucleotide runs with p-values |
| &emsp; `motif_hits` | All motif matches (TATA, CCAAT, Y Patch, Inr, DPE, secondary TATA) with scores, p-values, and distances to TSS |
| &emsp; `thresholds` | Per-motif score thresholds and significance cutoffs |
| &emsp; `run_params` | Inputs and schema geometry (flanks, core band) used to build the database |

### Example

```bash
python3 WindowScrubber.py \
    -t CT_sense.tss_summary.tsv \
    -g B73.PLATINUM.pseudomolecules-v1.fasta \
    -l CT_ltr_seqs_sense.fa \
    -db CT_motif_hits.db \
    --tata_mismatch 1 \
    --ccaat_mismatch 0 \
    --sig-alpha 1e-4
```

---

## Query_WSDB

Query interface for the WindowScrubber SQLite database. Retrieves, filters, and summarizes motif hits, CA dinucleotide runs, and promoter element composition across genes and TEs. Supports flexible ranking, distance-based filtering, threshold enforcement, and multi-motif presence queries.

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
| `--class` | Comma-separated element classes (e.g. `Gene,TIR`); applies to every query mode |
| `--anchor-mode` | `U3` or `TSS_window`; applies to every query mode |
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

R-based multi-omic integration and activity classification for LTR retrotransposons and genes. Merges expression evidence from any number of Illumina, PacBio, and ONT samples with ChIP-seq peaks, DNA methylation (UMR), CAGE, and promoter-motif tables, estimates per-sample expression thresholds by segmented regression, fuses evidence per tissue group with sigmoid scoring, and assigns activity labels (Constitutive, Facultative, Tissue-Specific, Developmental, Vegetative, Repressed, Background, Silent).

Inputs are declared in a **manifest** (one row per sample-file pairing) and a **config** (reference GFF, tissue groups, scoring parameters). Nothing about file names or column names is assumed by the code.

### Usage

```bash
mamba activate DECLTR-env
# 1. Draft a manifest from your input directories, then edit tissue/replicate by hand
scripts/make_manifest.sh -o samples.tsv \
    chip=ChIP_intersects/ umr=UMR_intersect.gff isoforms=IsoClassifier_out/ \
    motif=motifs/ illumina=merged_counts.tsv
# 2. Check paths, columns, and key resolution without running anything
Rscript DECLTR.r --manifest samples.tsv --config decltr_config.yml --validate
# 3. Run
Rscript DECLTR.r --manifest samples.tsv --config decltr_config.yml --out results/run1 --threads 8
```

Worked examples: `configs/decltr_manifest.example.tsv` and `configs/decltr_config.example.yml` reproduce the maize B73 analysis; `test/decltr_manifest.tsv` and `test/decltr_config.yml` run on the test dataset with a single ONT sample.

### Manifest

Tab-separated, one row per sample-file pairing. Paths may be absolute or relative to the manifest's directory.

| Column | Meaning |
|---|---|
| `sample_id` | Sample label used in output column names (e.g. `PB_Ear`, `ONT_CT1`, `Ears_H3K27ac.1`) |
| `platform` | `illumina`, `pacbio`, `ont`, `chip`, `umr`, `cage`, `motif` (free text; expression platforms are the ones listed under `scoring.weights` in the config) |
| `role` | `score`: the file's value column feeds the scoring assay for this sample. `extra`: columns are carried into the output only |
| `preset` | Reader preset (table below) |
| `path` | Input file |
| `column` | For matrix files (`count_matrix`): which column holds this sample |
| `tissue` | Tissue label used to build groups |
| `replicate` | Replicate label; samples sharing platform, tissue and replicate form one threshold unit |
| `scope` | `gene`, `te`, or `both`: which features the file describes (labels extras; several rows of one sample are summed into one assay column) |
| `options` | `key=value;key=value` overrides of the preset: `key_col`, `key_regex`, `value_col`, `agg` (`max`/`sum`/`first`), `fill`, `filter` (`col==value`), `keep_cols` (comma list or `*`) |

Every input is read by one generic routine: a delimited table with a key column (optionally a regex on it), a value column, an aggregation for repeated keys, a fill value for features absent from the file, and pass-through columns. Keys are resolved against both the feature `ID` and its `Parent` in the reference GFF, so IsoClassifier tables keyed on `repeat_region_N` and Illumina tables keyed on `LTRRT_N` both attach to the same locus.

| Preset | Reads | Key | Value (score rows) | Aggregation / fill |
|---|---|---|---|---|
| `intersect_gff_chip` | bedtools intersect of the annotation with peaks | col 9 `ID=` | col 14 (peak signal) | max / 0 |
| `intersect_gff_umr` | bedtools intersect with UMR calls | col 9 `ID=` | col 16 (mean methylation %) | first / 100 |
| `count_matrix` | one TSV with GFF columns, `Attributes`, and one count column per sample | `Attributes` `ID=` | column named in `column` | sum / 0 |
| `tss_summary_v1` | pre-v2 IsoClassifier `*_LTR_TSS.tsv` / `*_Gene_TSS.tsv` | `Feature` or `Gene` | `Total_Reads` | first / 0; keeps TSS1..Count2 |
| `tss_summary_v2` | current IsoClassifier `*_sense.tss_summary.tsv` | `Feature` | `Total_Reads` | filters `Orientation==sense`; keeps Class, Top_Isoform, TSS calls |
| `isoforms_v1` | pre-v2 `*_isoforms.tsv` (extras only) | `attrs` `Parent=` | – | keeps the category count columns |
| `isoforms_v2` | current `*_sense.isoforms.tsv` (extras only) | `Feature` | – | keeps every category column |
| `motif_tsv` | Query_WSDB `--best-all-motifs` output (extras only) | `feature` | – | keeps all columns; `*_present` become logical |
| `keyed_tsv` | any table with the key in column 1 | col 1 | col 2 | first / 0; use `keep_cols=*` for extras |

Anything else is best converted to a keyed TSV with a few lines of awk rather than taught to the reader; `scripts/standardize_cage.sh` does this for CAGE cluster GFFs.

### Config

YAML. Only `reference_gff` is required; everything else has defaults. `configs/decltr_config.example.yml` is annotated.

| Key | Meaning |
|---|---|
| `reference_gff`, `feature_types` | Annotation and regexes on GFF column 3 selecting the features to classify |
| `drop_contigs_regex`, `drop_attributes` | Contigs to exclude (e.g. scaffolds); GFF attributes not carried to the output |
| `groups` | Tissue groups for scoring: `Name: {platform: [members]}`, a member being a tissue or `{tissue:, replicate:, sample_id:}`. Omit to get one group per tissue per platform |
| `aliases` | Group -> label-level tissue (groups sharing an alias are collapsed by max before labeling) |
| `dev_groups`, `veg_groups` | Label-level tissues treated as developmental / vegetative |
| `scoring.weights`, `scoring.s` | Per-platform fusion weights (listing a platform makes it an expression platform) and sigmoid scale |
| `labels.*` | Activity thresholds (`active_thr`, `weak_thr`, `silent_thr`, `repress_thr`, `const_frac`, `dom_margin`, ...) |
| `thresholds.*` | Segmented-regression search range, ChIP range, UMR minimum signal, seed, `pass_by` (`unit` = pass if any replicate unit clears its threshold, default `total`) |
| `assay_columns` | `per_sample` (default): the score rows of one sample are summed into one column. `per_row`: one column per manifest row |
| `legacy_collapse_duplicate_loci` | Reproduce the pre-refactor merging of features sharing coordinates. Regression use only |
| `output_prefix` | Default output prefix (overridden by `--out`) |

### How the score is computed

1. **Threshold units.** Samples sharing platform, tissue and replicate are summed; a segmented regression of loci-remaining against count threshold gives each unit's breakpoint (`*_breakpoints.csv`).
2. **Candidate filter.** A feature is scored if its log1p total on any expression platform clears that platform's breakpoint (per unit for platforms with `pass_by: unit`), or any ChIP column clears its breakpoint, or its unmethylated signal is at least `umr_min_signal`.
3. **Group evidence.** Per group and platform, the row-median of log1p member columns minus log1p of the group threshold (median of member unit breakpoints) goes through a sigmoid; platforms are fused by weighted mean.
4. **Labels.** Groups are collapsed through `aliases`; breadth, dominance margin, and developmental vs vegetative means assign the activity label. Chromatin support turns weak-but-open loci into `Repressed`.

### Outputs

| File | Description |
|---|---|
| `<prefix>.qs` | Full data frame: annotation, one assay column per sample (`platform.sample_id`), extras (`platform.sample_id[.scope].field`), `Passed_Platforms`, and the label columns |
| `<prefix>_labels.tsv` | ID, coordinates, `Passed_Platforms`, `Activity`, `top_tissue`, `top_score`, `margin`, `breadth_active`, `dev_score`, `veg_score` |
| `<prefix>_breakpoints.csv` | Segmented-regression thresholds per unit, raw and log1p |
| `<prefix>_manifest.tsv`, `<prefix>_config.yml` | Copies of the inputs used |

### Tests

**Unit tests** cover the input layer and the data model. Fixtures are generated into a temp directory, so no test data is needed and the suite runs in seconds:

```bash
conda activate DECLTR-env
bash test/run_unit_tests.sh          # or: DECLTR_RSCRIPT=/path/to/Rscript bash test/run_unit_tests.sh
```

Results are written to `test/unit_results/`: `unit_test_results.txt` names every test and its expectation count, `unit_test_summary.csv` is the same as a table. The runner exits non-zero on any failure.

| Area | What is covered |
|---|---|
| Option parsing | `key=value` pairs, comma lists, numeric coercion, values containing `=` and `;`, malformed input |
| Column resolution | By name, by index, candidate lists, error messages when absent or out of range |
| Row filtering | `==` and `!=`, NA rows, unparseable expressions |
| Reference loading | `feature_types` selection, attribute expansion, `drop_attributes`, `Parent` handling, legacy collapse |
| Key resolution | ID and Parent both resolving to one locus; unresolvable keys counted rather than dropped |
| Every reader preset | Aggregation (`max`, `first`, `sum`), fill values, key regexes, `filter`, `keep_cols`, extras naming and typing |
| Assay construction | `per_sample` versus `per_row`, and the rules for combining a sample's score rows |
| Groups and output | `member_columns`, `resolve_groups` explicit and default, `filter_features`, wide-output naming and fill |

**Regression tests** compare a full run against a previous one. `test/compare_decltr_baseline.R baseline.qs baseline_breakpoints.csv <prefix>` checks the label columns and breakpoints; `test/diff_decltr_labels.R old_labels.tsv new_labels.tsv` prints a label-transition table. With `legacy_collapse_duplicate_loci: true` and `assay_columns: per_row` the refactored code reproduces the pre-refactor results exactly.

**Pipeline test.** `bash test/run_test.sh` runs the whole suite end to end on the bundled test dataset, ending with DECLTR on the IsoClassifier outputs. Set `DECLTR_RSCRIPT` if the DECLTR environment's `Rscript` is not on `PATH`; the DECLTR step is skipped with a warning when it is unavailable.

### Changes from the pre-refactor script

- Long-read group evidence was computed as the row-median across the gene table column and the TE table column of the same sample. One of the two is always zero for any locus, so PacBio evidence was halved and the ONT leaf group used half of the smallest replicate. `assay_columns: per_sample` sums the two tables first. On the B73 dataset this changes 10,833 of 80,143 scored labels, almost all toward broader activity (4,831 Vegetative and 3,169 Facultative loci become Constitutive; 2,333 Vegetative become Facultative); pass flags and breakpoints are unchanged.
- Features sharing chromosome, type, start, end and strand (about 195 loci in B73) were merged into one row whose ID matched nothing. They are now kept as separate features; because the segmented breakpoints sit close to integer boundaries, this alone moved several Illumina thresholds by one count.
- Motif `*_present` flags were parsed with `as.logical("1")`, which is `NA` in R, so every locus was `FALSE`. They are now read as `== "1"`.
- The candidate set no longer depends on Illumina being present; any platform may be absent.
- Score rows of one sample are combined according to their fill value: counts (fill 0) are summed, while a non-zero fill such as the UMR sentinel of 100 means "no data", so those rows must cover disjoint features and are refused if they overlap.
- `key_regex` overrides in the `options` column could not contain a semicolon, which made the documented `ID=([^;]+)` form unusable. Options now split only at a `;` that begins the next `key=`.

---

## Developer Notes

**AccuMap:** Developed primarily for ONT long reads. Requires that reads have not been demultiplexed or trimmed, only base-called, for optimal usage with PyChopper. For non-ONT data, start at the `--run_cut` step as long as reads still contain homopolymers (poly(A) or poly(T)).

**IsoClassifier:** Best used with a GFF containing annotations for all transposons and genes to enable nested-element removal. The GFF *must* contain records for the left and right long terminal repeats (lLTR and rLTR as named by EDTA) of each LTR-RT, in addition to the full `LTR_retrotransposon` annotation for each structurally intact locus. When `--genome-fasta` is provided, IsoClassifier runs U3/promoter extraction internally (superseding the standalone `U3_Seq_Extractor.py`).

**WindowScrubber:** Uses IsoClassifier outputs to define search windows for each motif context. External TSS data (e.g., from CAGE) can substitute for IsoClassifier TSS calls if the file matches IsoClassifier's column format (`Feature`, `Class`, `Strand`, `Chrom`, 0-based `TSS1`). The LTR FASTA (IsoClassifier with `--genome-fasta`) is needed only for the U3 schema. Databases built before the class-aware schema reverse-complemented already-oriented FASTA records, so minus-strand elements in them were scanned on the wrong strand; rebuild them.

**DECLTR:** Feature keys are resolved through both `ID` and `Parent`, so the reference GFF must carry `Parent=repeat_region_N` on structural LTR-RT records (EDTA convention) for IsoClassifier and WindowScrubber tables to attach. Segmented-regression breakpoints are floored to integers and several sit within a few hundredths of a boundary; small changes in the feature set can move a threshold by one count, so keep the same `feature_types` and reference when comparing runs.

**SpliceJunTest:** This is a file written for testing canonical splice junctions in our data. It is published here since it is mentioned in the associated manuscript but should not be considered part of the official LTRRT Promoter Suite tools. 
