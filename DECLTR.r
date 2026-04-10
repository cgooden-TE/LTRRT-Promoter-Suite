#!/usr/bin/env Rscript

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


#----------------------------------------------
#  Paths and file lists
#----------------------------------------------
data_dir <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/DECLTR/GenomicData"
isoform_dir <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/DECLTR/IsoformData"
motif_dir <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/DECLTR/MotifData"
combined_ref <- "/home/caleb/data/genome_and_annotations/TE_B73_UpdatedStrands_wGenes.gff"
illumina_path <- "/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/NAM_Reproc/Intersects/B73NAM_MergedCounts_StrucLTRsGenes.tsv"

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
  if (grepl("CT_[0-9]+", base, ignore.case = TRUE)) {
    label <- sub("(?i)(CT_[0-9]+).*", "\\1", base, perl = TRUE)
  } else if (grepl("CTMerge", base, ignore.case = TRUE)) {
    label <- "CTMerge"
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
  label <- ifelse(grepl("ONTPB", basename(fp), ignore.case = TRUE), "ONTPB", "ONTPB")

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
#  Add B73 NAM data
#----------------------------------------------
read_illumina_counts <- function(path, id_key = "ID", count_prefix = NULL) {
  df <- readr::read_tsv(path, col_types = readr::cols(.default = "c"), comment = "")

  # normalize header
  names(df) <- trimws(gsub("\uFEFF", "", names(df)))

  if (!"Attributes" %in% names(df)) stop("Illumina file must have an 'Attributes' column: ", path)

  # pull ID from Attributes
  df$ID <- sub(paste0(".*\\b", id_key, "=([^;]+).*"), "\\1", df$Attributes)
  df$ID <- trimws(df$ID)

  # detect count columns: keep everything that is NOT GFF-ish and NOT Attributes/ID
  gff_cols <- intersect(names(df), c("Chr", "Source", "Name", "Start", "End", "Score", "Strand", "Phase"))
  drop_cols <- unique(c(gff_cols, "Attributes", "ID"))

  count_cols <- setdiff(names(df), drop_cols)

  if (length(count_cols) == 0) stop("No count columns detected in Illumina file: ", path)

  # convert counts to numeric safely
  df <- df |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(count_cols),
        ~ suppressWarnings(as.numeric(.x))
      )
    )

  # if you want to prefix to avoid collisions with existing columns, do it here
  if (!is.null(count_prefix) && nzchar(count_prefix)) {
    df <- df |> dplyr::rename_with(~ paste0(count_prefix, "_", .x), dplyr::all_of(count_cols))
    count_cols <- paste0(count_prefix, "_", count_cols)
  }

  # aggregate duplicates (very common if file has multiple rows per ID)
  df_out <- df |>
    dplyr::select(ID, dplyr::all_of(count_cols)) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  df_out
}
message("Processing Illumina counts: ", basename(illumina_path))

illumina_df <- read_illumina_counts(illumina_path, id_key = "ID", count_prefix = NULL)

merged_df <- dplyr::left_join(merged_df, illumina_df, by = "ID")

# fill any NAs introduced by join (counts -> 0)
merged_df <- merged_df |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(setdiff(names(illumina_df), "ID")), ~ tidyr::replace_na(.x, 0))
  )


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
#  Filter, Pass to Thresholding
#----------------------------------------------
filtered_merged_df <- merged_df %>%
  filter(!str_detect(Chr, regex("scaf", ignore_case = TRUE)))

decltr_res_df <- filtered_merged_df |>
  dplyr::mutate(Classification = coalesce(Classification, "Gene"))

decltr_res_df$ID <- trimws(as.character(decltr_res_df$ID))

# Coalesce helper (treat missing counts as 0 for all downstream sums)
coalesce0 <- function(x) dplyr::coalesce(x, 0)

# Expand a semicolon-separated key=value attribute column into wide columns
expand_attr_column <- function(df, attr_col, prefix = NULL) {
  if (!attr_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in data frame.", attr_col))
  }

  base_cols <- setdiff(names(df), attr_col)

  out <- df |>
    mutate(
      .rid  = row_number(),
      .attr = .data[[attr_col]]
    ) |>
    separate_rows(.attr, sep = ";\\s*") |>
    separate(.attr, into = c("key", "value"), sep = "=", fill = "right", extra = "merge") |>
    mutate(key = str_trim(key), value = str_trim(value)) |>
    filter(!is.na(key), key != "") |>
    pivot_wider(
      id_cols     = c(.rid, all_of(base_cols)),
      names_from  = key,
      values_from = value,
      values_fn   = ~ paste(unique(na.omit(.x)), collapse = ";")
    )

  if (!is.null(prefix) && nzchar(prefix)) {
    new_cols <- setdiff(names(out), c(".rid", base_cols))
    names(out)[match(new_cols, names(out))] <- paste0(prefix, new_cols)
  }

  out |> select(-.rid)
}

# Logistic on log-delta. delta=0 => evidence 0.5.
sigmoid01_delta <- function(delta, s = 0.35) {
  1 / (1 + exp(-delta / s))
}

# Median of thresholds for a group of Illumina columns
# (unified: no feature_type distinction)
group_thr_illumina <- function(cols, THRESH_ILLUMINA) {
  thr_vec <- THRESH_ILLUMINA[cols]
  thr <- suppressWarnings(median(as.numeric(thr_vec), na.rm = TRUE))
  if (!is.finite(thr)) thr <- 1
  thr
}

# PB: read from THRESH_PB by tissue name
# (unified: no feature_type distinction)
group_thr_pb <- function(tissue, THRESH_PB) {
  thr <- THRESH_PB[[tissue]]
  thr <- as.numeric(thr)
  if (!is.finite(thr)) thr <- 1
  thr
}

# ONT: one unified threshold
group_thr_ont <- function(THRESH_ONT_CT) {
  thr <- as.numeric(THRESH_ONT_CT)
  if (!is.finite(thr)) thr <- 1
  thr
}

# =========================
# 2) Column discovery
# =========================
nm <- names(decltr_res_df)

## Illumina expression columns (counts)
illumina_cols <- grep("^B73NAM_", nm, value = TRUE)

## All PB_* columns (for numeric casting)
pb_all_cols <- grep("PB_", nm, value = TRUE)

## All per-locus PB expression columns (total reads)
pb_expr_cols_gene <- grep("^gene_Total_Reads_PB_", nm, value = TRUE)
pb_expr_cols_ltr  <- grep("^LTR_Total_Reads_PB_",  nm, value = TRUE)
pb_expr_cols_all  <- c(pb_expr_cols_gene, pb_expr_cols_ltr)

## Helper: PB expression columns per tissue (for grouping)
pb_expr_cols_by_tissue <- function(tissue) {
  grep(paste0("Total_Reads_PB_", tissue, "_"), pb_expr_cols_all, value = TRUE)
}

## Tissues present in the PB expression columns
pb_tissues_expr <- unique(sub(".*Total_Reads_PB_([^_]+)_.*", "\\1", pb_expr_cols_all))
pb_tissues_expr <- sort(pb_tissues_expr)

## ONT CT merged columns, per feature type
ont_ct_merge_cols_gene <- grep("^gene_Total_Reads_ONT_CTMerge", nm, value = TRUE)
ont_ct_merge_cols_ltr  <- grep("^LTR_Total_Reads_ONT_CTMerge",  nm, value = TRUE)
ont_ct_merge_cols      <- c(ont_ct_merge_cols_gene, ont_ct_merge_cols_ltr)

## ONT Cold columns (if any)
cold_cols <- grep("Cold", nm, value = TRUE)

## ChIP & UMR
chip_cols <- grep("H3K4me3|H3K27ac", nm, value = TRUE)
umr_col   <- grep("_UMR$", nm, value = TRUE)[1]
if (length(umr_col) == 0) stop("Could not find a column ending with '_UMR'.")

## Ensure relevant numeric columns are numeric
to_numeric <- unique(c(
  illumina_cols,
  pb_all_cols,
  ont_ct_merge_cols,
  cold_cols,
  chip_cols,
  umr_col
))
for (cc in to_numeric) {
  if (cc %in% names(decltr_res_df)) {
    suppressWarnings({
      decltr_res_df[[cc]] <- as.numeric(decltr_res_df[[cc]])
    })
  }
}

# =========================
# 2.5) Regression Helper
# =========================
estimate_seg_threshold <- function(counts, max_k = 20L, default = 1L) {
  counts <- coalesce0(counts)
  max_count <- max(counts, na.rm = TRUE)
  if (!is.finite(max_count) || max_count < 1) return(default)

  max_k <- min(max_k, max_count)

  seg_df <- data.frame(
    threshold      = seq_len(max_k),
    loci_remaining = vapply(seq_len(max_k),
                            function(t) sum(counts >= t, na.rm = TRUE),
                            integer(1))
  )

  lm0 <- lm(loci_remaining ~ threshold, data = seg_df)
  seg <- try(
    segmented::segmented(lm0, seg.Z = ~ threshold, psi = list(threshold = 5)),
    silent = TRUE
  )

  if (inherits(seg, "try-error")) return(default)

  est <- seg$psi[1, "Est."]
  thr <- floor(est)
  thr <- max(1L, min(thr, max_k))
  thr
}

# =========================
# 3) ONT segmented regression — UNIFIED across all loci
# =========================

# All ONT merged counts (genes + LTRs together)
ont_all_counts <- if (length(ont_ct_merge_cols)) {
  coalesce0(rowSums(decltr_res_df[, ont_ct_merge_cols, drop = FALSE], na.rm = TRUE))
} else {
  rep(0, nrow(decltr_res_df))
}

ONT_THRESH <- estimate_seg_threshold(ont_all_counts, default = 1L)

THRESH_ONT_CT <- ONT_THRESH   # single numeric value

## TEMP
ONT_THRESH_COLD <- 1

# =========================
# 4) PacBio segmented regression per tissue — UNIFIED across all loci
# =========================

PB_THRESH <- setNames(rep(NA_real_, length(pb_tissues_expr)), pb_tissues_expr)

for (tissue in pb_tissues_expr) {
  # Collect both gene and LTR columns for this tissue
  gene_col <- grep(paste0("^gene_Total_Reads_PB_", tissue, "_"), nm, value = TRUE)
  ltr_col  <- grep(paste0("^LTR_Total_Reads_PB_",  tissue, "_"), nm, value = TRUE)

  all_tissue_cols <- c(gene_col, ltr_col)

  if (length(all_tissue_cols)) {
    # Sum across all columns for this tissue (or just use first col if only one)
    if (length(all_tissue_cols) == 1L) {
      counts_all <- coalesce0(decltr_res_df[[all_tissue_cols]])
    } else {
      counts_all <- coalesce0(rowSums(
        decltr_res_df[, all_tissue_cols, drop = FALSE], na.rm = TRUE
      ))
    }
    PB_THRESH[[tissue]] <- estimate_seg_threshold(counts_all, default = 1L)
  } else {
    PB_THRESH[[tissue]] <- 1L
  }
}

THRESH_PB <- PB_THRESH   # named vector: tissue -> single threshold

# =========================
# 5) Illumina segmented regression per sample — UNIFIED across all loci
# =========================

ILL_THRESH <- setNames(rep(NA_real_, length(illumina_cols)), illumina_cols)

for (cc in illumina_cols) {
  v <- coalesce0(decltr_res_df[[cc]])
  ILL_THRESH[[cc]] <- estimate_seg_threshold(v, default = 1L)
}

THRESH_ILLUMINA <- ILL_THRESH   # named vector: column -> single threshold

# ===============================================
# 5b) ChIP segmented regression, UMR thresholding
# ===============================================

CHIP_THRESH <- setNames(rep(NA_real_, length(chip_cols)), chip_cols)
for (cc in chip_cols) {
  v <- coalesce0(decltr_res_df[[cc]])
  CHIP_THRESH[[cc]] <- estimate_seg_threshold(v, max_k = 50L, default = 1L)
}

present_chip <- decltr_res_df[, chip_cols, drop = FALSE]
for (cc in chip_cols) {
  present_chip[[cc]] <- coalesce0(decltr_res_df[[cc]]) >= CHIP_THRESH[[cc]]
}

umr_signal <- 1 - (replace(decltr_res_df[[umr_col]], decltr_res_df[[umr_col]] == 100, NA) / 100)
present_umr <- is.finite(umr_signal) & (umr_signal >= 0.1)

# =========================
# Tissue groups (Illumina + PB/ONT)
# =========================

b73_tissue <- function(pattern) {
  grep(paste0("^B73.*_", pattern, "_"), nm, value = TRUE, ignore.case = TRUE)
}

v11_base_cols   <- grep("^B73.*Base_",   nm, value = TRUE, ignore.case = TRUE)
v11_middle_cols <- grep("^B73.*Middle_", nm, value = TRUE, ignore.case = TRUE)
v11_tip_cols    <- grep("^B73.*Tip_",    nm, value = TRUE, ignore.case = TRUE)
shoot_cols      <- grep("^B73.*_Shoot_", nm, value = TRUE, ignore.case = TRUE)

groups <- list(
  Embryo     = list(illumina = b73_tissue("Embryo"),    pb = pb_expr_cols_by_tissue("Embryo")),
  Endosperm  = list(illumina = b73_tissue("Endosperm"), pb = pb_expr_cols_by_tissue("Endosperm")),
  Root       = list(illumina = b73_tissue("Root"),      pb = pb_expr_cols_by_tissue("Root")),
  Shoot      = list(illumina = shoot_cols),
  V11_Base   = list(illumina = v11_base_cols),
  V11_Middle = list(illumina = v11_middle_cols),
  V11_Tip    = list(illumina = v11_tip_cols),
  Ear        = list(illumina = b73_tissue("Ear"),    pb = pb_expr_cols_by_tissue("Ear")),
  Tassel     = list(illumina = b73_tissue("Tassel"), pb = pb_expr_cols_by_tissue("Tassel")),
  Anther     = list(illumina = b73_tissue("Anther")),
  Pollen     = list(pb = pb_expr_cols_by_tissue("Pollen")),
  Leaf       = list(
    illumina = c(v11_base_cols, v11_middle_cols, v11_tip_cols),
    ont      = ont_ct_merge_cols
  )
)

# ===============================================
# 6) Score-based candidate filter — UNIFIED (no gene/LTR split)
# ===============================================

# ---- per-view totals (fast) ----
ill_tot <- if (length(illumina_cols)) {
  X <- as.matrix(decltr_res_df[, illumina_cols, drop = FALSE])
  storage.mode(X) <- "numeric"
  X[!is.finite(X)] <- 0
  rowSums(X)
} else {
  rep(0, nrow(decltr_res_df))
}
pb_tot <- if (length(pb_expr_cols_all)) {
  X <- as.matrix(decltr_res_df[, pb_expr_cols_all, drop = FALSE])
  storage.mode(X) <- "numeric"
  X[!is.finite(X)] <- 0
  rowSums(X)
} else {
  rep(0, nrow(decltr_res_df))
}
ont_tot <- if (length(ont_ct_merge_cols)) {
  X <- as.matrix(decltr_res_df[, ont_ct_merge_cols, drop = FALSE])
  storage.mode(X) <- "numeric"
  X[!is.finite(X)] <- 0
  rowSums(X)
} else {
  rep(0, nrow(decltr_res_df))
}
ill_tot_log  <- log1p(ill_tot)
pb_tot_log   <- log1p(pb_tot)
ont_tot_log  <- log1p(ont_tot)

# ---- UNIFIED summary thresholds (all loci together) ----
ILL_SUM_THR_log  <- estimate_seg_threshold(ill_tot_log,  default = 1L)
PB_SUM_THR_log   <- estimate_seg_threshold(pb_tot_log,   default = 1L)
ONT_SUM_THR_log  <- estimate_seg_threshold(ont_tot_log,  default = 1L)

segmented_breakpoints <- data.frame(
  Platform         = c("Illumina", "PacBio", "ONT"),
  Breakpoint_log1p = c(ILL_SUM_THR_log, PB_SUM_THR_log, ONT_SUM_THR_log)
)
segmented_breakpoints$Breakpoint_raw_scale <- exp(segmented_breakpoints$Breakpoint_log1p) - 1
segmented_breakpoints
write.csv(segmented_breakpoints, "Segmented_Breakpoints.csv", row.names = FALSE)

# ---- UNIFIED pass flags (same threshold for genes and LTRs) ----
pass_ill_sum <- ill_tot_log >= ILL_SUM_THR_log
pass_pb_sum  <- pb_tot_log  >= PB_SUM_THR_log
pass_ont_sum <- ont_tot_log >= ONT_SUM_THR_log

# Note: Classification column is retained in the data for interpretability in
# the output, but plays no role in filtering or threshold decisions.

decltr_res_df$Passed_Platforms <- paste0(
  ifelse(pass_ill_sum, "Illumina;", ""),
  ifelse(pass_pb_sum,  "PacBio;",   ""),
  ifelse(pass_ont_sum, "ONT;",      "")
)
decltr_res_df$Passed_Platforms[decltr_res_df$Passed_Platforms == ""] <- NA

# ---- chromatin keep ----
keep_any_chip <- if (length(chip_cols)) {
  rowSums(as.matrix(present_chip[, chip_cols, drop = FALSE]), na.rm = TRUE) > 0
} else {
  rep(FALSE, nrow(decltr_res_df))
}

umr_vec <- as.numeric(decltr_res_df[[umr_col]])
umr_vec[umr_vec == 100] <- NA
umr_signal       <- 1 - (umr_vec / 100)
umr_signal       <- pmin(pmax(umr_signal, 0), 1)
keep_any_umr     <- is.finite(umr_signal) & (umr_signal >= 0.1)

# ---- final candidate set (unified filter) ----
keep_for_scoring <- pass_ill_sum | pass_pb_sum | pass_ont_sum | keep_any_chip | keep_any_umr

# Single unified subset — genes and LTRs filtered identically
decltr_sub <- decltr_res_df[keep_for_scoring, , drop = FALSE]

nrow(decltr_sub)

compute_chrom_support <- function(df, chip_cols, present_chip_df, umr_col) {
  chip_pass_frac <- if (length(chip_cols)) {
    M <- as.matrix(present_chip_df[, chip_cols, drop = FALSE])
    storage.mode(M) <- "logical"
    rowMeans(M, na.rm = TRUE)
  } else {
    rep(0, nrow(df))
  }

  umr_vec    <- as.numeric(df[[umr_col]])
  umr_vec[umr_vec == 100] <- NA
  umr_signal <- 1 - (umr_vec / 100)
  umr_signal <- pmin(pmax(umr_signal, 0), 1)
  umr_signal[!is.finite(umr_signal)] <- NA_real_

  chrom_support <- pmax(chip_pass_frac, umr_signal, na.rm = TRUE)
  chrom_support[!is.finite(chrom_support)] <- 0
  chrom_support
}

present_chip_sub <- if (length(chip_cols)) {
  out <- decltr_sub[, chip_cols, drop = FALSE]
  for (cc in chip_cols) out[[cc]] <- coalesce0(decltr_sub[[cc]]) >= CHIP_THRESH[[cc]]
  out
} else {
  NULL
}
chrom_support <- compute_chrom_support(decltr_sub, chip_cols, present_chip_sub, umr_col)
names(chrom_support) <- decltr_sub$ID

log1p_safe <- function(x) log1p(coalesce0(as.numeric(x)))

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  mu  <- rowMeans(mat, na.rm = TRUE)
  sd  <- matrixStats::rowSds(mat, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sd, "/")
}

dominance_margin <- function(v_named) {
  v <- sort(as.numeric(v_named), decreasing = TRUE)
  if (length(v) < 2) return(Inf)
  v[1] - v[2]
}

# ------------------------------------------------------------
# BUILD OMICS BLOCKS (raw -> log1p -> row-z) from decltr_res_df
# ------------------------------------------------------------
build_omics_matrices <- function(decltr_res_df,
                                 groups,
                                 illumina_cols,
                                 pb_expr_cols_all,
                                 ont_ct_merge_cols,
                                 chip_cols,
                                 umr_col) {
  df <- decltr_res_df %>%
    mutate(ID = trimws(as.character(ID))) %>%
    distinct(ID, .keep_all = TRUE)

  loci <- df$ID

  summarize_group_median <- function(df, cols) {
    cols <- intersect(cols, names(df))
    if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
    x <- as.matrix(df[, cols, drop = FALSE])
    matrixStats::rowMedians(x, na.rm = TRUE)
  }

  # ---- Illumina group summaries ----
  ill_group_names <- names(groups)[vapply(groups, function(g) !is.null(g$illumina), logical(1))]
  ill_mat <- sapply(ill_group_names, function(gname) {
    cols <- groups[[gname]]$illumina
    tmp  <- df
    if (length(cols)) tmp[, cols] <- lapply(tmp[, cols, drop = FALSE], log1p_safe)
    summarize_group_median(tmp, cols)
  })
  ill_mat <- t(ill_mat)
  rownames(ill_mat) <- ill_group_names
  colnames(ill_mat) <- loci

  # ---- PacBio group summaries ----
  pb_mat <- NULL
  pb_group_names <- names(groups)[vapply(groups, function(g) !is.null(g$pb), logical(1))]
  if (length(pb_group_names)) {
    pb_mat <- sapply(pb_group_names, function(gname) {
      cols <- groups[[gname]]$pb
      tmp  <- df
      if (length(cols)) tmp[, cols] <- lapply(tmp[, cols, drop = FALSE], log1p_safe)
      summarize_group_median(tmp, cols)
    })
    pb_mat <- t(pb_mat)
    rownames(pb_mat) <- pb_group_names
    colnames(pb_mat) <- loci
  }

  # ---- ONT (leaf only) ----
  ont_mat <- NULL
  if (length(ont_ct_merge_cols) > 0) {
    tmp <- df
    tmp[, ont_ct_merge_cols] <- lapply(tmp[, ont_ct_merge_cols, drop = FALSE], log1p_safe)
    ont_leaf <- summarize_group_median(tmp, ont_ct_merge_cols)
    ont_mat  <- matrix(ont_leaf, nrow = 1, dimnames = list("ONT_Leaf", loci))
  }

  # ---- ChIP ----
  chip_mat <- NULL
  if (length(chip_cols) > 0) {
    tmp <- df
    tmp[, chip_cols] <- lapply(tmp[, chip_cols, drop = FALSE], log1p_safe)
    chip_mat <- t(as.matrix(tmp[, chip_cols, drop = FALSE]))
    rownames(chip_mat) <- chip_cols
    colnames(chip_mat) <- loci
  }

  # ---- UMR ----
  umr_vec    <- as.numeric(df[[umr_col]])
  umr_vec[umr_vec == 100] <- NA
  umr_signal <- 1 - (umr_vec / 100)
  umr_signal <- pmin(pmax(umr_signal, 0), 1)
  umr_mat    <- matrix(umr_signal, nrow = 1, dimnames = list("UMR_leaf_like", loci))

  omics_log <- list(Illumina = ill_mat, PacBio = pb_mat, ONT = ont_mat, ChIP = chip_mat, UMR = umr_mat)
  omics_log <- omics_log[!vapply(omics_log, is.null, logical(1))]

  list(loci = loci, omics_log = omics_log)
}

# ---- ChIP evidence per locus ----
chip_pass_frac <- if (length(chip_cols)) {
  M <- as.matrix(present_chip[, chip_cols, drop = FALSE])
  storage.mode(M) <- "logical"
  rowMeans(M, na.rm = TRUE)
} else {
  rep(0, nrow(decltr_res_df))
}
umr_ev <- umr_signal
umr_ev[!is.finite(umr_ev)] <- NA_real_

# ------------------------------------------------------------
# SCORE MODEL (per tissue group)
# ------------------------------------------------------------

sigmoid01 <- function(z, z0 = 0.0, s = 1.0) {
  1 / (1 + exp(-(z - z0) / s))
}

# Fuse evidence per tissue group across platforms (Illumina/PB/ONT)
compute_activity_scores_log <- function(
    omics_log,
    groups,
    THRESH_ILLUMINA,
    THRESH_PB,
    THRESH_ONT_CT,
    w = list(Illumina = 1.0, PacBio = 1.0, ONT = 0.7),
    s = 0.25) {
  loci         <- colnames(omics_log$Illumina)
  tissue_names <- names(groups)

  butter <- matrix(NA_real_,
    nrow = length(tissue_names), ncol = length(loci),
    dimnames = list(tissue_names, loci)
  )

  for (t in tissue_names) {
    # --- Illumina evidence ---
    ill_ev <- rep(NA_real_, length(loci))
    if (!is.null(omics_log$Illumina) && t %in% rownames(omics_log$Illumina) && !is.null(groups[[t]]$illumina)) {
      thr   <- group_thr_illumina(groups[[t]]$illumina, THRESH_ILLUMINA)
      mid   <- log1p(thr)
      delta <- as.numeric(omics_log$Illumina[t, ]) - mid
      ill_ev <- sigmoid01_delta(delta, s = s)
    }

    # --- PacBio evidence ---
    pb_ev <- rep(NA_real_, length(loci))
    if (!is.null(omics_log$PacBio) && t %in% rownames(omics_log$PacBio)) {
      thr   <- group_thr_pb(t, THRESH_PB)
      mid   <- log1p(thr)
      delta <- as.numeric(omics_log$PacBio[t, ]) - mid
      pb_ev <- sigmoid01_delta(delta, s = s)
    }

    # --- ONT evidence (Leaf only) ---
    ont_ev <- rep(NA_real_, length(loci))
    if (!is.null(omics_log$ONT) && "ONT_Leaf" %in% rownames(omics_log$ONT) && t %in% c("Leaf")) {
      thr   <- group_thr_ont(THRESH_ONT_CT)
      mid   <- log1p(thr)
      delta <- as.numeric(omics_log$ONT["ONT_Leaf", ]) - mid
      ont_ev <- sigmoid01_delta(delta, s = s)
    }

    # --- fuse (weighted mean, NA-safe per locus) ---
    num <- rep(0, length(loci))
    den <- rep(0, length(loci))

    ok <- is.finite(ill_ev);  num[ok] <- num[ok] + w$Illumina * ill_ev[ok]; den[ok] <- den[ok] + w$Illumina
    ok <- is.finite(pb_ev);   num[ok] <- num[ok] + w$PacBio   * pb_ev[ok];  den[ok] <- den[ok] + w$PacBio
    ok <- is.finite(ont_ev);  num[ok] <- num[ok] + w$ONT      * ont_ev[ok]; den[ok] <- den[ok] + w$ONT

    out <- rep(NA_real_, length(loci))
    ok2 <- den > 0
    out[ok2] <- num[ok2] / den[ok2]
    butter[t, ] <- out
  }

  butter
}

# ------------------------------------------------------------
# RULE LABELS from tissue activity matrix
# ------------------------------------------------------------
label_loci_from_activity <- function(
  butter,
  chrom_support    = NULL,
  dev_groups       = c("Embryo", "Endosperm", "Anther", "Pollen"),
  veg_groups       = c("Leaf", "Root", "Ear", "Tassel"),
  dom_margin       = 0.12,
  min_facultative  = 2L
) {
  active_thr    <- 0.75
  weak_thr      <- 0.45
  silent_thr    <- 0.30
  repress_thr   <- 0.35
  const_frac    <- 0.50
  dev_veg_gap   <- 0.12
  veg_cap       <- 0.55
  veg_dev_gap   <- 0.12
  dev_cap       <- 0.55

  loci        <- colnames(butter)
  tissues_raw <- rownames(butter)

  leaf_alias <- function(x) {
    x <- as.character(x)
    x[x %in% c("V11_Base", "V11_Middle", "V11_Tip")] <- "Leaf"
    x
  }
  tissues_alias <- leaf_alias(tissues_raw)
  alias_levels  <- unique(tissues_alias)

  butter_alias <- sapply(alias_levels, function(tt) {
    rows <- which(tissues_alias == tt)
    if (length(rows) == 1L) return(butter[rows, ])
    apply(butter[rows, , drop = FALSE], 2, max, na.rm = TRUE)
  })
  butter_alias <- t(butter_alias)
  rownames(butter_alias) <- alias_levels
  colnames(butter_alias) <- loci

  dev_rows  <- intersect(dev_groups, rownames(butter_alias))
  veg_rows  <- intersect(veg_groups, rownames(butter_alias))
  dev_score <- if (length(dev_rows)) colMeans(butter_alias[dev_rows, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, length(loci))
  veg_score <- if (length(veg_rows)) colMeans(butter_alias[veg_rows, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, length(loci))

  top_tissue <- apply(butter_alias, 2, function(x) {
    if (all(!is.finite(x))) return(NA_character_)
    rownames(butter_alias)[which.max(x)]
  })
  top_score <- apply(butter_alias, 2, function(x) if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE))

  margin <- apply(butter_alias, 2, function(x) {
    x <- sort(x[is.finite(x)], decreasing = TRUE)
    if (length(x) < 2) return(Inf)
    x[1] - x[2]
  })

  breadth_active <- apply(butter_alias, 2, function(x) sum(x >= active_thr, na.rm = TRUE))
  n_tis   <- nrow(butter_alias)
  const_n <- max(3L, floor(n_tis * const_frac))

  label <- rep("Background", length(loci))
  label[is.finite(top_score) & top_score < silent_thr] <- "Silent"
  idx_bg <- is.finite(top_score) & (top_score >= silent_thr) & (top_score < weak_thr)
  label[idx_bg] <- "Background"

  is_active <- is.finite(top_score) & (top_score >= active_thr)

  if (!is.null(chrom_support)) {
    cs  <- chrom_support[loci]
    cs[!is.finite(cs)] <- 0
    idx_rep <- (!is_active) & is.finite(top_score) & (top_score >= silent_thr) & (cs >= repress_thr)
    label[idx_rep] <- "Repressed"
  }

  idx_const  <- is_active & (breadth_active >= const_n)
  idx_single <- is_active & (breadth_active == 1L)
  label[idx_const] <- "Constitutive"

  idx_not_const_active <- is_active & !idx_const

  idx_dev <- idx_not_const_active &
    is.finite(dev_score) & is.finite(veg_score) &
    ((dev_score - veg_score) >= dev_veg_gap) &
    (veg_score <= veg_cap)

  idx_veg <- idx_not_const_active &
    is.finite(dev_score) & is.finite(veg_score) &
    ((veg_score - dev_score) >= veg_dev_gap) &
    (dev_score <= dev_cap)

  label[idx_dev] <- "Developmental"
  label[idx_veg] <- "Vegetative"

  idx_fac <- is_active &
    (breadth_active >= min_facultative) &
    (breadth_active < const_n) &
    !(idx_dev | idx_veg)
  label[idx_fac] <- "Facultative"

  idx_ts <- idx_single & (breadth_active == 1L) & (margin >= dom_margin)
  label[idx_ts] <- paste0("Tissue-Specific:", top_tissue[idx_ts])

  data.frame(
    ID             = loci,
    Activity       = label,
    top_tissue     = top_tissue,
    top_score      = top_score,
    margin         = margin,
    breadth_active = breadth_active,
    dev_score      = dev_score,
    veg_score      = veg_score,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# BUILD OMICS (on unified filtered subset) -> activity scores -> labels
# ------------------------------------------------------------

mats <- build_omics_matrices(
  decltr_res_df     = decltr_sub,
  groups            = groups,
  illumina_cols     = illumina_cols,
  pb_expr_cols_all  = pb_expr_cols_all,
  ont_ct_merge_cols = ont_ct_merge_cols,
  chip_cols         = chip_cols,
  umr_col           = umr_col
)

butter <- compute_activity_scores_log(
  omics_log       = mats$omics_log,
  groups          = groups,
  THRESH_ILLUMINA = THRESH_ILLUMINA,
  THRESH_PB       = THRESH_PB,
  THRESH_ONT_CT   = THRESH_ONT_CT,
  s               = 0.25
)

labels <- label_loci_from_activity(
  butter,
  chrom_support = chrom_support,
  dev_groups    = c("Embryo", "Endosperm", "Anther", "Pollen"),
  veg_groups    = c("Leaf", "Root", "Ear", "Tassel")
)

decltr_labeled <- decltr_res_df %>%
  dplyr::left_join(labels, by = "ID")

qsave(decltr_labeled, "0310_samtools-Clip_Exon_Cleave_NAM_NewScoring_GeneLTR_Unif_PycFix.qs",
  preset = "balanced"
)
