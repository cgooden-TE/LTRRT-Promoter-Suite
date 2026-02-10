## Revised script, v06/17/25 — fixes Isoform, TSS, and Motif merges
## NOTE that this script, because it is based on outputs of previous steps of the LTR pipeline, 
##  is designed explicitly to study intact (structural) LTRs. The renaming and merging is
##  built on the assumption that the input files are based EDTA gff annotation naming convention for 
##  LTRs and their repeat regions. 

#----------------------------------------------
#  Load libraries
#----------------------------------------------
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(DESeq2)
library(BiocParallel)
library(openxlsx)
library(qs)
register(SerialParam())

#----------------------------------------------
#  Paths and file lists
#----------------------------------------------
data_dir <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/DECLTR/GenomicData"
isoform_dir <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/DECLTR/IsoformData"
motif_dir <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/DECLTR/MotifData"
combined_ref <- "/home/caleb/data/genome_and_annotations/TE_B73_UpdatedStrands_wGenes.gff"

data_fps <- list.files(data_dir, full.names = TRUE)
isoform_fps <- list.files(isoform_dir, pattern = ".tsv$", full.names = TRUE)
motif_fps <- list.files(motif_dir, pattern = ".tsv$", full.names = TRUE)

#----------------------------------------------
#  Reset merged_df to start fresh
#----------------------------------------------
if (exists("merged_df")) rm(merged_df)

#----------------------------------------------
#  Read base GFF annotation into DataFrame
#----------------------------------------------
te_df <- read_tsv(combined_ref,
  comment = "#", col_names = FALSE,
  col_types = cols(.default = "c")
) %>%
  tidyr::separate_rows(X9, sep = ";") %>%
  tidyr::separate(X9, into = c("key", "value"), sep = "=", fill = "right") %>%
  dplyr::select(X1, X3, X4, X5, X7, key, value) %>%
  tidyr::pivot_wider(
    names_from = key, 
    values_from = value,
    values_fn = list(value = toString)
  ) %>%
  dplyr::rename(
    Chr   = X1,
    Type  = X3,
    Start = X4,
    End   = X5,
    Strand = X7
  ) %>%
  filter(
    Type == "gene" |
      str_detect(Type, "LTR_retrotransposon")
  ) %>%
  relocate(ID) %>%
  dplyr::select(-Type)

#----------------------------------------------
#  Merge overlap/vector data
#----------------------------------------------

# Record the IDs
vectors_df <- data.frame(ID = te_df$ID, stringsAsFactors = FALSE)

for (fp in data_fps) {
  print(paste("Processing file:", basename(fp)))
  # strip off suffix to get name like "Ears_H3K27ac.1" or "LTRRTs_Genes_UMR"
  nm <- sub("_peaks\\.intsct\\.gff$|_intsct\\.gff$", "", basename(fp))

  # read in the file
  df <- read.table(fp,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE
  )

  # extract TE IDs from the 9th column
  ids <- sapply(
    strsplit(df[, 9], ";"),
    function(x) sub(".*ID=([^;]+).*", "\\1", x[1])
  )

  if (grepl("H3K", nm)) {
    # ChIP‐seq intersects: take the maximum of column 14 (read count)
    tmp <- data.frame(
      ID = ids,
      val = df[, 14],
      stringsAsFactors = FALSE
    )
    agg <- tmp %>%
      group_by(ID) %>%
      summarize(val = max(val, na.rm = TRUE), .groups = "drop")
    default_val <- 0
  } else if (grepl("UMR", nm)) {
    # UMR intersects: take just the first column 16 value (mean meth)
    tmp <- data.frame(
      ID = ids,
      val = df[, 16],
      stringsAsFactors = FALSE
    )
    agg <- tmp[!duplicated(tmp$ID), , drop = FALSE]
    default_val <- 100
  } else {
    # skip anything else
    next
  }

  # Rename and merge, fill NAs with file‐specific default
  colnames(agg)[2] <- nm
  vectors_df <- left_join(vectors_df, agg, by = "ID")
  vectors_df[[nm]][is.na(vectors_df[[nm]])] <- default_val
}

# Attach new columns back to te_df
merged_df <- left_join(te_df, vectors_df, by = "ID")

#----------------------------------------------
#  Merge TSS & totals (Gene & LTR)
#----------------------------------------------
tss_files <- list.files(isoform_dir, pattern = "TSS", full.names = TRUE)
for (fp in tss_files) {
  print(paste("Processing file:", basename(fp)))
  label <- if (grepl("CT", basename(fp), ignore.case = TRUE)) {
    sub("(?i)[._-]?TSS.*$", "", basename(fp), perl = TRUE)
  } else if (grepl("Cold", basename(fp), ignore.case = TRUE)) {
    "Cold"
  } else if (grepl("^PB", basename(fp), ignore.case = TRUE)) {
    # for files like PB_Ear_TSS.tsv, this gives "PB_Ear"
    sub("(?i)[._-]?TSS.*$", "", basename(fp), perl = TRUE)
  } else if (grepl("ONTPB", basename(fp), ignore.case = TRUE)) {
	  label <- "ONTPB"
  } else {
    stop("unrecognized classifier file: ", fp)
  }
  df_tss <- readr::read_tsv(fp, col_types = cols(.default = "c"))
  if (grepl("Gene", basename(fp), ignore.case = TRUE)) {
    df_tss <- df_tss %>%
      dplyr::rename(ID = Gene) %>%
      dplyr::select(ID, Total_Reads, TSS1, Count1, TSS2, Count2) %>%
      dplyr::mutate(
        Total_Reads = as.numeric(Total_Reads),
        TSS1        = as.integer(TSS1),
        Count1      = as.numeric(Count1),
        TSS2        = as.integer(TSS2),
        Count2      = as.numeric(Count2)
      ) %>%
      dplyr::rename_with(~ paste0("gene_", .x, "_", label), -ID)
    # join genes by ID → ID
    merged_df <- left_join(merged_df, df_tss, by = "ID")
  } else if (grepl("LTR", basename(fp), ignore.case = TRUE)) {
    df_tss <- df_tss %>%
      dplyr::rename(ID = Feature) %>%
      dplyr::select(ID, Total_Reads, TSS1, Count1, TSS2, Count2) %>%
      dplyr::mutate(
        Total_Reads = as.numeric(Total_Reads),
        TSS1        = as.integer(TSS1),
        Count1      = as.numeric(Count1),
        TSS2        = as.integer(TSS2),
        Count2      = as.numeric(Count2)
      ) %>%
      dplyr::rename_with(~ paste0("LTR_", .x, "_", label), -ID)
    merged_df <- dplyr::left_join(merged_df, df_tss,
      by = c("Parent" = "ID")
    )
  } else {
    next
}
}
# fill NAs post-TSS merge
merged_df <- merged_df %>%
  dplyr::mutate(
    dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)),
    dplyr::across(where(is.logical), ~ tidyr::replace_na(., FALSE))
  )

#----------------------------------------------
#  Merge classifier isoform files (counts)
#----------------------------------------------
classifier_files <- list.files(isoform_dir, pattern = "isoforms.tsv$", full.names = TRUE)
for (fp in classifier_files) {
  # pull the filename into a variable
  base <- basename(fp)
  print(paste("Processing file:", base))
  # decide label
  if (grepl("CT", base, ignore.case = TRUE)) {
    label <- "CT"
  } else if (grepl("Cold", base, ignore.case = TRUE)) {
    label <- "Cold"
  } else if (grepl("^PB", base, ignore.case = TRUE)) {
    # for files like PB_Ear_isoforms.tsv, this gives "PB_Ear"
    label <- sub("^(PB_[^_]+).*", "\\1", base)
  } else if (grepl("ONTPB", base, ignore.case = TRUE)) {
    label <- "ONTPB"
  } else {
    stop("unrecognized classifier file: ", fp)
  }

  # Check if file is Gene or LTR
  is_gene <- grepl("Gene", base, ignore.case = TRUE)

  if (is_gene) {
    # Gene input: standard GFF, only one read count (col 10)
    df_counts <- readr::read_tsv(fp, col_types = cols(.default = "c")) %>%
      dplyr::mutate(
        ID = sub(".*ID=([^;]+).*", "\\1", attrs), # column 9 = attributes
        gene_count = as.numeric(int_count) # column 10 = read count
      ) %>%
      dplyr::select(ID, gene_count) %>%
      dplyr::rename_with(~ paste0(.x, "_", label), -ID)
    merged_df <- dplyr::left_join(merged_df, df_counts, by = "ID")
  } else {
    df_counts <- readr::read_tsv(fp, col_types = cols(.default = "c")) %>%
      dplyr::mutate(
        # extract the Parent ID as *character*
        ID = sub(".*Parent=([^;]+).*", "\\1", attrs),
        # convert only the count columns to numeric
        across(c(
          total_reads, ltr_left_reads, spliced_ltr_left, ltr_right_reads, spliced_ltr_right,
          spanning_reads, spliced_spanning, ro5_reads, spliced_ro5, ro3_reads, spliced_ro3
        ), as.numeric)
      ) %>%
      dplyr::select(
        ID, total_reads, ltr_left_reads, spliced_ltr_left, ltr_right_reads, spliced_ltr_right,
        spanning_reads, spliced_spanning, ro5_reads, spliced_ro5, ro3_reads, spliced_ro3
      ) %>%
      dplyr::rename_with(~ paste0(.x, "_", label), -ID)

    merged_df <- dplyr::left_join(
      merged_df,
      df_counts,
      by = c("Parent" = "ID")
    )
  }
}

# fill any NAs
merged_df <- merged_df %>%
  mutate(
    across(where(is.numeric), ~ replace_na(., 0)),
    across(where(is.logical), ~ replace_na(., FALSE))
  )

#----------------------------------------------
#  Merge ONT-only motif data (both LTRs and genes)
#----------------------------------------------
ont_files <- motif_fps[grepl("ONT", basename(motif_fps))]
for (fp in ont_files) {
  print(paste("Processing file:", basename(fp)))
  label <- ifelse(grepl("CT", basename(fp), ignore.case = TRUE), "CT", "Merge")

  # 1) read in, ensure an ID column
  df_m <- readr::read_tsv(fp, col_types = cols(.default = "c"))
  if ("Feature" %in% names(df_m)) {
    names(df_m)[names(df_m) == "Feature"] <- "ID"
  } else if ("feature" %in% names(df_m)) {
    names(df_m)[names(df_m) == "feature"] <- "ID" 
  } else if (!"ID" %in% names(df_m)) {
    stop("No 'Feature' or 'ID' column in motif file: ", fp)
  }

  # 2) coerce types
  df_m <- df_m %>%
    mutate(
      across(contains("present"), as.logical),
      across(matches("^(Dist_|TSS_abs)"), as.integer)
    )

  # 3) rename all motif columns to avoid collisions
  #    - for LTRs prefix with "LTR_motif_",
  #    - for genes "gene_motif_"
  is_gene <- grepl("Gene", basename(fp), ignore.case = TRUE)
  prefix <- if (is_gene) "gene_motif_" else "LTR_motif_"
  df_m <- df_m %>% rename_with(~ paste0(prefix, .x, "_", label), -ID)

  # 4) join back into merged_df
  if (is_gene) {
    # join on ID == ID
    merged_df <- left_join(merged_df, df_m, by = "ID")
  } else {
    # join on Parent == ID for LTRs
    merged_df <- left_join(merged_df, df_m, by = c("Parent" = "ID"))
  }
}

# 5) finally fill in NAs
merged_df <- merged_df %>%
  mutate(
    across(where(is.numeric), ~ replace_na(., 0)),
    across(where(is.logical), ~ replace_na(., FALSE))
  )

#----------------------------------------------
#  Add CAGE data
#----------------------------------------------
cage_dir <- "/home/caleb/data/PaperWritingReruns/GenBio_25/CAGE"
cage_files <- list.files(cage_dir, pattern = "formatted.gff$", full.names = TRUE)

attach_cage_from_gff <- function(df, gff_path, id_col = "ID", col_prefix = NULL) {
  g <- read.table(gff_path,
    sep = "\t", header = FALSE,
    quote = "", comment.char = "#", fill = TRUE,
    stringsAsFactors = FALSE, colClasses = "character"
  )
  stopifnot(ncol(g) >= 14)

  id9 <- sub(".*\\bID=([^;]+).*", "\\1", g[[9]])

  as_int_safe <- function(x) suppressWarnings(as.integer(x))
  cage_start <- as_int_safe(g[[10]])
  cage_end <- as_int_safe(g[[11]])
  cage_str <- g[[12]]
  cage_dTSS <- as_int_safe(g[[13]])

  shape_field <- g[[14]]
  m <- regexpr("Shape=([^;]+)", shape_field, perl = TRUE)
  cage_shape <- ifelse(m > 0, sub("Shape=([^;]+).*", "\\1", regmatches(shape_field, m)), NA)

  cage_df <- data.frame(
    ID = id9,
    CAGE_Start = cage_start,
    CAGE_End = cage_end,
    CAGE_Str = cage_str,
    CAGE_dTSS = cage_dTSS,
    CAGE_Shape = cage_shape,
    stringsAsFactors = FALSE
  )

  # ---- NEW: prefix the added columns when requested ----
  if (!is.null(col_prefix) && nzchar(col_prefix)) {
    nc <- names(cage_df) != "ID"
    names(cage_df)[nc] <- paste0(col_prefix, "_", names(cage_df)[nc])
  }
  # ------------------------------------------------------

  cage_df <- cage_df[!duplicated(cage_df$ID), ]

  names(cage_df)[names(cage_df) == "ID"] <- id_col
  merge(df, cage_df, by = id_col, all.x = TRUE, sort = FALSE)
}

for (fp in cage_files) {
  base <- basename(fp)
  message("Processing file: ", base)

  # first word before the first underscore, e.g. "Root" from "Root_dTSS_..."
  prefix <- sub("_.*$", "", base)

  merged_df <- attach_cage_from_gff(merged_df, fp, id_col = "ID", col_prefix = prefix)
}

#----------------------------------------------
#  Quick Cleanup
#----------------------------------------------
merged_df <- merged_df %>% dplyr::select(
  -motif, -tsd, -TSD, -TIR,
  -biotype, -logic_name, -Sequence_ontology
)
merged_df <- merged_df %>%
  dplyr::select(
    -dplyr::contains("classifier_source"),
    -dplyr::contains("source_file")
  )
#----------------------------------------------
#  Filter and save
#----------------------------------------------
filtered_merged_df <- merged_df %>%
  filter(!str_detect(Chr, regex("scaf", ignore_case = TRUE)))

qsave(filtered_merged_df, "020426_samtools-Clip_Exon_Cleave_df.qs",
  preset = "balanced"
)

#----------------------------------------------
#  DESeq2 and Plotting 
#----------------------------------------------
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(ggrepel)
library(BiocParallel)
register(SerialParam())

## THIS SECTION SHOULD BE USED WITH A COUNTS TABLE, EITHER BY
## FEATURE-COUNTS OR MADE USING THE TSS_FINDER
## Dependencies: dplyr, DESeq2
# Load required libraries

# Read in counts table
counts <- read.delim("/home/caleb/data/XingliSeq/ReprocComplete/FeatureCounts/All_Samples_count_matrix_wGenes.tsv",
 header = TRUE, row.names = 1
)

# Make a condition vector — adjust to match actual condition setup
condition <- factor(c("control", "control", "control", "cold", "cold", "cold"),
  levels = c("control", "cold")
)
col_data <- data.frame(row.names = colnames(counts), condition)

# Split into genes and TEs
gene_counts <- counts[grepl("^Zm", rownames(counts)), ]
all_counts <- counts

# Create DESeqDataSet using genes only
dds_genes <- DESeqDataSetFromMatrix(countData = gene_counts, colData = col_data, design = ~condition)

# Estimate size factors based on gene expression
dds_genes <- estimateSizeFactors(dds_genes)
size_factors <- sizeFactors(dds_genes)

# Apply size factors to full dataset (genes + TEs)
dds_all <- DESeqDataSetFromMatrix(countData = all_counts, colData = col_data, design = ~condition)
sizeFactors(dds_all) <- size_factors

# Filter: remove rows with low counts (e.g., row sum < 10)
keep <- rowSums(counts(dds_all, normalized = TRUE)) >= 10
dds_all_filtered <- dds_all[keep, ]

# Run DESeq on the full dataset using gene-based normalization
dds_all_filtered <- DESeq(dds_all_filtered)

# Get results
res <- results(dds_all_filtered)
res_order <- res[order(res$padj, na.last = NA), ]
deseq_results_df <- data.frame(
  Name = rownames(res_order),
  Log2FoldChange = res_order$log2FoldChange,
  PValue = res_order$pvalue,
  PAdj = res_order$padj
)
saveRDS(deseq_results_df, file = "deseq2_counts-table_results.rds")

## THIS SECTION SHOULD BE USED WITH MAPPED BAM FILES AS INPUT
# Use the combined ref defined above to do DESeq on TEs and genes
# Make into GRanges object, add column for names to use later
deseq_gff <- read.table(combined_ref, header = FALSE, stringsAsFactors = FALSE)
colnames(deseq_gff) <- c(
  "seqid", "source", "sequence_ontology",
  "start", "end", "score", "strand", "phase", "attributes"
)

deseq_gff$Name <- sapply(strsplit(deseq_gff$attributes, ";"), function(attr) {
  # Extract TE family in "Name=", gene returns NA
  name_attr <- grep("Name=", attr, value = TRUE)
  if (length(name_attr) > 0) {
    sub("Name=", "", name_attr)
  } else {
    NA  # return NA if no Name attribute is found
  }
})
deseq_gff$Class <- sapply(strsplit(deseq_gff$attributes, ";"), function(attr) {
  # Look for TE superfamily in "Classification=", gene NA
  class_attr <- grep("Classification=", attr, value = TRUE)
  if (length(class_attr) > 0) {
    sub("Classification=", "", class_attr)
  } else {
    NA # return NA if no Class attribute is found
  }
})
deseq_gff$ID <- sapply(strsplit(deseq_gff$attributes, ";"), function(attr) {
  # Look for superfamily in "ID="
  id_attr <- grep("ID=", attr, value = TRUE)
  if (length(id_attr) > 0) {
    sub("ID=", "", id_attr)
  } else {
    NA # return NA if no Class attribute is found
  }
})

## Filter TE reference by TEs which are retained in TE_OverlapDF
filtered_deseq_gff <- deseq_gff[deseq_gff$ID %in% rownames(te_gene_anno_df), ]
filtered_deseq_gff$strand[filtered_deseq_gff$strand == "?"] <- "*"
te_grob <- makeGRangesFromDataFrame(filtered_deseq_gff,
                                    seqnames = "seqid",
                                    start.field = "start",
                                    end.field = "end",
                                    strand = "strand",
                                    keep.extra.columns = TRUE)
names(te_grob) <- filtered_deseq_gff$ID

# Paths to bam alignments, using long-read cDNA reads
bampath <- "/home/caleb/data/XingliSeq/MatrixReady/BAMData/"
bamfiles <- list.files(
  path = bampath,
  full.names = TRUE
)

# Count the reads
se <- summarizeOverlaps(features = te_grob, reads = bamfiles, mode = "Union",
                        singleEnd = TRUE, ignore.strand = TRUE,
                        fragments = FALSE)

condition <- factor(c(rep("Cold", 3), rep("Control", 3)))
condition <- relevel(condition, ref = "Control")
colData(se) <- DataFrame(condition = condition)
deseqobj <- DESeqDataSet(se, design = ~ condition)

# Run DESeq2 with its object and order results by p-value
deseqobj <- DESeq(deseqobj)
res <- results(deseqobj)
res_order <- res[order(res$padj, na.last = NA), ]
deseq_results_df <- data.frame(
  Name = rownames(res_order),
  Log2FoldChange = res_order$log2FoldChange,
  PValue = res_order$pvalue,
  PAdj = res_order$padj
)
saveRDS(deseq_results_df, file = "deseq2_results.rds")



## THIS SECTION INCLUDES SEVERAL PRE-MADE PLOTTING FUNCTIONS AND CAN BE ADJUSTED 
##  TO COLOR CERTAIN ELEMENTS OF INTEREST 
# Add a column for significance (optional)
deseq_results_df$Significant <- with(
  deseq_results_df, ifelse(PAdj < 0.05 & abs(Log2FoldChange) > 1, "Yes", "No")
)

# Create the plot
deseq_volcano <- ggplot(
  deseq_results_df, aes(x = Log2FoldChange, y = -log10(PAdj))) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  labs(
    title = "Volcano Plot of DESeq2 Results",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted P-value)"
  ) +
  theme_minimal()
ggsave("volcano_plot.png",
  plot = deseq_volcano, width = 8, height = 6, dpi = 300
)

## Creat a Volcano Plot for Significantly DE LTRs
# Create a filtered set
ltrrt_sig <- deseq_results_df %>%
  filter(PAdj < 0.05, grepl("^LTRRT", Name))

# Add a column to the full dataset for LTRRT significance
deseq_results_df$SigLTR <- ifelse(
  deseq_results_df$Name %in% ltrrt_sig$Name, "Significant LTR", "Other"
)
# Separate data for plotting layers
df_other <- deseq_results_df %>% filter(SigLTR == "Other")
df_ltrrt <- deseq_results_df %>% filter(SigLTR == "Significant LTR")

# Plot with layers
ltr_volcano <- ggplot() +
  geom_point(
    data = df_other, aes(x = Log2FoldChange, y = -log10(PAdj)),
    color = "gray", alpha = 0.5, size = 1
  ) +
  geom_point(
    data = df_ltrrt, aes(x = Log2FoldChange, y = -log10(PAdj)),
    color = "red", size = 1, alpha = 0.9
  ) +
  labs(
    title = "Volcano Plot: Highlighting Significant LTRRT Elements",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Legend"
  ) +
  theme_minimal()

# Save the plot
ggsave("LTR_Sig_Cand_volcano_plot.png",
  plot = ltr_volcano, width = 8, height = 6, dpi = 300
)

## Save excel files for DE elements; separate sig LTRs
# Create workbook
library(openxlsx)
deseq_wb <- createWorkbook()
# Add each data frame to its own sheet
addWorksheet(deseq_wb, "All Results")
writeData(deseq_wb, "All Results", deseq_results_df)
addWorksheet(deseq_wb, "LTRRT Significant")
writeData(deseq_wb, "LTRRT Significant", ltrrt_sig)
# Save the workbook
saveWorkbook(deseq_wb, file = "deseq2_results_sheets.xlsx", overwrite = TRUE)

################################################################################
# Running DESeq for genes hit by LTR-originated reads
# Run and edit as necessary
################################################################################

hitgenes_path <- " "
hitgenes_list <- list.files(path = hitgenes_path, full.names = TRUE)

# Function to count the occurrences of each transcript ID and store it as DF
read_and_counttrxpt <- function(file) {
  exon_info <- readLines(file)
  exon_counts <- as.data.frame(table(exon_info))
  colnames(exon_counts) <- c("exon", "count")
  exon_counts
}
# Read all the count files for when LTR-origin read ran into an exon of gene
exon_counts_list <- lapply(hitgenes_list, read_and_counttrxpt)
# Merge
exon_counts_matrix <- Reduce(function(x, y) {
  full_join(x, y, by = "exon")
}, exon_counts_list)
# Replace missing values with 0
exon_counts_matrix[is.na(exon_counts_matrix)] <- 0

sample_names <- gsub(".*/|\\.gff$", "", hitgenes_list)
colnames(exon_counts_matrix)[-1] <- sample_names

# Set row names to exon transcript IDs and remove the "exon" column
rownames(exon_counts_matrix) <- exon_counts_matrix$exon
exon_counts_matrix <- exon_counts_matrix[-1]

hitgene_cond <- c(rep("control", 3), rep("cold_stress", 3))
hitegene_data <- data.frame(
  condition = factor(hitgene_cond),
  row.names = colnames(exon_counts_matrix)
)
hitegene_data$condition <- relevel(hitegene_data$condition, ref = "control")
# Create DESeq object, run pipeline, convert results to DF
hitgene_dds <- DESeqDataSetFromMatrix(countData = exon_counts_matrix,
                                      colData = hitegene_data,
                                      design = ~ condition)
hitgene_dds <- DESeq(hitgene_dds)
hitgene_res <- results(hitgene_dds)
hitgene_res <- as.data.frame(hitgene_res)
hitgene_res$gene <- rownames(hitgene_res)

################################################################################
# MANIPULATING DATA FOR LTRs AND TOP EXPRESSION
# Run and edit as necessary
################################################################################
res_order_df <- as.data.frame(res_order)
