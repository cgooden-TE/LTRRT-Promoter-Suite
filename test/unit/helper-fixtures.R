# Fixtures for the DECLTR unit tests.
#
# Every fixture is written from this file rather than checked in, so the expected
# values in a test sit next to the input that produces them. Files land in a
# per-session temp directory that R removes on exit.

# Source the DECLTR library. DECLTR_ROOT is set by run_unit_tests.sh; the fallback
# covers running test_dir() from test/unit by hand.
decltr_root <- Sys.getenv("DECLTR_ROOT", unset = normalizePath("../..", mustWork = FALSE))
for (f in c("io.R", "model.R", "thresholds.R", "scoring.R", "labels.R", "validate.R"))
  source(file.path(decltr_root, "decltr", f))

FIXTURE_DIR <- file.path(tempdir(), "decltr_fixtures")
dir.create(FIXTURE_DIR, showWarnings = FALSE, recursive = TRUE)

# Write tab-separated lines built from character vectors, one row per element.
write_tsv_lines <- function(name, rows, header = NULL) {
  path <- file.path(FIXTURE_DIR, name)
  writeLines(c(header, rows), path)
  path
}

tsv <- function(...) paste(..., sep = "\t")

# --- reference annotation -------------------------------------------------
# 2 genes, 1 structural LTR-RT (keyed on Parent downstream), 1 homology fragment,
# 1 helitron (excluded by feature_types), 1 scaffold gene, and 2 genes sharing
# identical coordinates (the pre-refactor collapse case).
fx_reference_gff <- function() {
  write_tsv_lines("reference.gff", c(
    tsv("chr1", "EDTA", "gene", "100", "200", ".", "+", ".",
        "ID=Zm00001eb000010;biotype=protein_coding;logic_name=cshl_gene"),
    tsv("chr1", "EDTA", "gene", "300", "400", ".", "-", ".",
        "ID=Zm00001eb000020;biotype=protein_coding"),
    tsv("chr1", "EDTA", "Gypsy_LTR_retrotransposon", "1000", "2000", ".", "+", ".",
        "ID=LTRRT_1;Parent=repeat_region_1;Classification=LTR/Gypsy;Method=structural;tsd=AAAAA"),
    tsv("chr1", "EDTA", "Copia_LTR_retrotransposon", "3000", "3500", ".", "-", ".",
        "ID=TE_homo_1;Classification=LTR/Copia;Method=homology"),
    tsv("chr1", "EDTA", "helitron", "4000", "4500", ".", "+", ".",
        "ID=TE_homo_2;Classification=DNA/Helitron"),
    tsv("scaf_9", "EDTA", "gene", "10", "20", ".", "+", ".",
        "ID=Zm00001eb999999;biotype=protein_coding"),
    tsv("chr1", "EDTA", "gene", "5000", "5100", ".", "+", ".", "ID=Zm00001ebDUPA"),
    tsv("chr1", "EDTA", "gene", "5000", "5100", ".", "+", ".", "ID=Zm00001ebDUPB")
  ))
}

# Features as the readers see them, without the legacy collapse.
fx_features <- function() {
  read_reference_gff(fx_reference_gff(), c("^gene$", "LTR_retrotransposon"),
                     drop_attributes = c("tsd", "biotype", "logic_name"))
}

# --- ChIP intersect (headerless, key col 9, value col 14) -----------------
# Zm00001eb000010 appears three times so agg=max is distinguishable from first/sum.
fx_chip_gff <- function() {
  row <- function(id, signal) tsv("chr1", "EDTA", "gene", "100", "200", ".", "+", ".",
                                  paste0("ID=", id), "chr1", "90", "210", "150", signal, "peak_1")
  write_tsv_lines("chip.intsct.gff", c(
    row("Zm00001eb000010", "7"),
    row("Zm00001eb000010", "42"),
    row("Zm00001eb000010", "13"),
    row("LTRRT_1", "5")
  ))
}

# --- UMR intersect (headerless, key col 9, value col 16) ------------------
# Zm00001eb000010 appears twice so agg=first is distinguishable from max.
fx_umr_gff <- function() {
  row <- function(id, meth) tsv("chr1", "EDTA", "gene", "100", "200", ".", "+", ".",
                                paste0("ID=", id), "chr1", "90", "210", "UMR", "35", "36", meth, "6.25")
  write_tsv_lines("umr.intsct.gff", c(
    row("Zm00001eb000010", "4.5"),
    row("Zm00001eb000010", "88.0"),
    row("LTRRT_1", "95.0")
  ))
}

# --- Illumina count matrix (key from Attributes, one column per sample) ---
# Zm00001eb000020 appears on two rows so agg=sum is distinguishable from first.
fx_count_matrix <- function() {
  hdr <- tsv("Chr", "Source", "Name", "Start", "End", "Score", "Strand", "Phase", "Attributes",
             "SampleA", "SampleB")
  row <- function(id, a, b) tsv("chr1", "EDTA", "x", "1", "2", ".", "+", ".",
                                paste0("ID=", id, ";Name=n"), a, b)
  write_tsv_lines("counts.tsv", c(
    row("Zm00001eb000010", "10", "1"),
    row("Zm00001eb000020", "3", "2"),
    row("Zm00001eb000020", "4", "5"),
    row("LTRRT_1", "6", "0")
  ), header = hdr)
}

# --- IsoClassifier v1 TSS summaries --------------------------------------
# Gene table keys on 'Gene'; LTR table keys on 'Feature' and carries the Parent id.
fx_tss_v1_gene <- function() {
  write_tsv_lines("sample_Gene_TSS.tsv", c(
    tsv("Zm00001eb000010", "25", "150", "20", "151", "5"),
    tsv("Zm00001eb000020", "8", "350", "8", "NA", "0")
  ), header = tsv("Gene", "Total_Reads", "TSS1", "Count1", "TSS2", "Count2"))
}

fx_tss_v1_ltr <- function() {
  write_tsv_lines("sample_LTR_TSS.tsv", c(
    tsv("repeat_region_1", "ro3", "+", "12", "1500", "9", "1502", "3")
  ), header = tsv("Feature", "Top_Isoform", "Top_Isoform_Strand", "Total_Reads",
                  "TSS1", "Count1", "TSS2", "Count2"))
}

# --- IsoClassifier v2 (Class column, sense/antisense in one file) ---------
fx_tss_v2 <- function() {
  write_tsv_lines("sample_sense.tss_summary.tsv", c(
    tsv("Zm00001eb000010", "Gene", "sense", "spliced", "+", "25", "1", "0.4000", "150", "20", "151", "5"),
    tsv("repeat_region_1", "LTR_structural", "sense", "ltr5_contained", "+", "12", "1", "0.7500", "1500", "9", "1502", "3"),
    tsv("Zm00001eb000020", "Gene", "antisense", "readout", "-", "99", "0", "0.1000", "350", "9", "NA", "0")
  ), header = tsv("Feature", "Class", "Orientation", "Top_Isoform", "Strand", "Total_Reads",
                  "TSS_Called", "TSS_Peak_Frac", "TSS1", "Count1", "TSS2", "Count2"))
}

fx_isoforms_v2 <- function() {
  write_tsv_lines("sample_sense.isoforms.tsv", c(
    tsv("Zm00001eb000010", "chr1", "100", "200", "Gene", "n", "+", "annotated", "sense",
        "25", "NA", "NA", "7", "631.4", "ID=Zm00001eb000010"),
    tsv("repeat_region_1", "chr1", "1000", "2000", "LTR_structural", "n", "+", "annotated", "sense",
        "12", "NA", "NA", "4", "500.0", "ID=LTRRT_1;Parent=repeat_region_1")
  ), header = tsv("Feature", "Chrom", "Start", "End", "Class", "Name", "Strand", "Strand_Source",
                  "Orientation", "Total_Reads", "Nested_In", "Nested_In_Class",
                  "spliced_reads", "mean_len_spliced", "Attrs"))
}

# --- IsoClassifier v1 isoform table (key from attrs Parent=) --------------
fx_isoforms_v1 <- function() {
  write_tsv_lines("sample_isoforms.tsv", c(
    tsv("chr1", "elem", "1000", "2000",
        "ID=LTRRT_1;Parent=repeat_region_1;Classification=LTR/Gypsy",
        "12", "4", "2", "3", "1", "0", "0", "0", "0", "5", "2")
  ), header = tsv("chrom", "name", "start", "end", "attrs", "total_reads",
                  "ltr_left_reads", "spliced_ltr_left", "ltr_right_reads", "spliced_ltr_right",
                  "spanning_reads", "spliced_spanning", "ro5_reads", "spliced_ro5",
                  "ro3_reads", "spliced_ro3"))
}

# --- Motif table (the *_present flags are "1"/"0", not TRUE/FALSE) --------
fx_motif <- function() {
  write_tsv_lines("sample_motifs.tsv", c(
    tsv("Zm00001eb000010", "150", "1", "120", "124", "30", "-1", "-1", "-1"),
    tsv("Zm00001eb000020", "350", "0", "-1", "-1", "-1", "340", "347", "10")
  ), header = tsv("feature", "tss_abs", "TA_rich_present", "CCAAT_start_abs", "CCAAT_end_abs",
                  "Dist_TSS_to_CCAAT", "TATA_start_abs", "TATA_end_abs", "Dist_TSS_to_TATA"))
}

# --- Generic keyed TSV (the escape hatch preset) --------------------------
fx_keyed <- function() {
  write_tsv_lines("keyed.tsv", c(
    tsv("Zm00001eb000010", "11", "Sharp"),
    tsv("LTRRT_1", "22", "Broad")
  ), header = tsv("key", "CAGE_dTSS", "CAGE_Shape"))
}

# --- Manifest construction ------------------------------------------------
# Build a manifest data frame in memory and hand it through read_manifest so the
# tests exercise the same validation and row_id derivation the CLI uses.
make_manifest <- function(...) {
  rows <- list(...)
  cols <- c("sample_id", "platform", "role", "preset", "path", "column",
            "tissue", "replicate", "scope", "options")
  df <- do.call(rbind, lapply(rows, function(r) {
    r <- utils::modifyList(list(sample_id = NA, platform = NA, role = "score", preset = NA,
                                path = NA, column = NA, tissue = NA, replicate = NA,
                                scope = "both", options = NA), r)
    as.data.frame(r[cols], stringsAsFactors = FALSE)
  }))
  path <- file.path(FIXTURE_DIR, "manifest.tsv")
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  read_manifest(path)
}

# A single manifest row, ready for read_keyed_table().
manifest_row <- function(...) make_manifest(list(...))[1, , drop = FALSE]

# Keymap over the fixture reference.
fx_keymap <- function() build_keymap(fx_features())
