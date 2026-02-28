# LTRRT-Promoter-Suite
Tools developed to characterize the U3 promoter regions and transcription events of long terminal repeat retrotransposons in maize.

### Developer Notes
**AccuMap:** AccuMap was developed primarily for use with ONT long-reads. It requires that reads not have been demultiplexed or trimemed, only base-called, for optimal usage conditions in PyChopper. If not using ONT sequencing data, the pipeline can be started at the --run_cut step as long as reads have not been previously trimmed and still contain homoploymers (poly(A) or poly(T)). \
**IsoClassifier:** IsoClassifier is best used with a feature reference containing annotations for all transposons and genes for your organism of interest in order to remove nested features. It *must* contain records for the long terminal repeats (lTLR and rLTR as named by EDTA) of each LTR-RT in addition to the full annotation for each structurally intact LTR-RT locus. \
**WindowScrubber:** WindowScrubber uses the outputs of IsoClassifier to set its search window for each motif context. If you have your own TSS identifications (for example from CAGE), these can be used, but the U3/LTR sequence files from IsoClassifier are still required and the file formatting / name conventions must match IsoClassifier's TSS output to find the primary TSS coordinate. \
**DECLTR:** DECLTR is made publicly available for purposes of development transparency, but it should be noted that the methods implemented are tuned to the specific data types and naming conventions utilized in this analysis. Reuse of this method will require re-tooling the functions for your data until the more broadly-applicable version is released.


### Usage
First, git clone this repository and build the environments included as .yml files. 
```
git clone https://github.com/cgooden-TE/LTRRT-Promoter-Suite
cd LTRRT-Promoter-Suite

# Env used for AccuMap, IsoClassifier, and WindowScrubber
mamba env create -f Nanopore_pipeline_final.yml

# Env used for WGCNA and DECLTR
mamba env create -f R-DESeq-env.yml

# Activate the environment
mamba activate <env>
```

#### AccuMap
**Minimum inputs:**\
    --fq, Input FASTQ of untrimmed, demultiplexed reads\
    --pyc, PyChopper output FASTQ filename\
    --sample, Sample name prefix for output files (ex. Control_1)\
    --ref, Reference genome. Ensure chromosome nomenclature is consistent.\

**Optional inputs:**\
    --run_pyc, Runs PyChopper with default settings\
    --run_cut, Runs Cutadapt with default settings\
    --run_map, Runs Minimap2 with default settings\
    --kit, for ONT reads, default="PCB114"\

**Outputs (using provided sample name prefix):**\
    - cutadapt.fastq : FASTQ after Cutadapt trimming\
    - cutadapt.log : Log file from Cutadapt run\
    - pychopped.fastq : Reads FASTQ after PyChopper primer removal\
    - pychopper.log : Log file from PyChopper run\
    - pychopper.report.pdf : PDF report from PyChopper\
    - pychopper.rescued.fastq : FASTQ of rescued reads from PyChopper\
    - minimap2.log : Log file from Minimap2 run\
    - minimap2.sam : SAM file from Minimap2 alignment\
    - strandtags.tsv : TSV file of read names and PyChopper strand orientations\
    - STtagged.bam : BAM file with PyChopper strand tags (ST) added\
    - STtagged.bed : BED file conversion of BAM with strand information\
    - STtagged.sorted.bam : Sorted BAM file with PyChopper strand tags\


#### IsoClassifier
**Minimum inputs:**\
    --gff, GFF reference with LTR elements, terminal repeated, and genes\
    --bam, One or more BAM files with mapped reads preferentially from AccuMap\
    --output, Prefix for default files\
    --tss_out, File name for the file of primary and secondary TSS information for LTR-RTs\
    --gene_out, File name for the file of primary and secondary TSS information for genes\

**Optional inputs:**\
    --length_mode

**Outputs (using provided sample name prefix):**\
    - TSS summary for LTR-RT loci\
    - TSS summary for gene loci\
    - LTR isoforms with read counts, lengths, splicing, and junctions per locus, isoforms.tsv\
    - Top 10 predicted TSS sites and their read contributions for LTR-RTs, 10site.ltr_tss_summary.tsv\
    - Top 10 predicted TSS sites and their read contributions for genes, 10site.gene_tss_summary.tsv\
    - Cleavage summary: Cleavage peak, fraction of reads contributing to that peak, top 10 predicted poly(A)-associated cleavage sites, their read contributions for LTR-RTs, 10site.ltr_cleavage_summary.tsv\
    - Read contribution to -10:+10 nt positions centered on predicted primary gene TSS, gene_primary_tss_density.tsv\
    - Read contribution to -10:+10 nt positions centered on predicted secondary gene TSS, gene_secondary_tss_density.tsv\
    - Read contribution to -10:+10 nt positions centered on predicted primary LTR-RT TSS, ltr_primary_tss_density.tsv\
    - Read contribution to -10:+10 nt positions centered on predicted secondary LTR-RT TSS, ltr_secondary_tss_density.tsv\
    - Count and lengths for introns and exons found in genes per read ID linked to feature and isoform category, gene_exon_stats_per_read.tsv\
    - Count and lengths for introns and exons found in LTR-RTs per read ID linked to feature and isoform category, ltr_exon_stats_per_read.tsv\
    - Summary of soft clipping found at the 3'-end of spliced reads aligned to LTR-RTs, ltr_3p_softclip_per_read_spliced.tsv\
    - Summary of soft clipping found at the 3'-end of nonspliced reads aligned to LTR-RTs, ltr_3p_softclip_per_read_nonspliced.tsv\


#### WindowScrubber

#### DECLTR

#### WGCNA
