# Unit tests for the DECLTR input layer: option parsing, column resolution,
# reference loading, key resolution, and every reader preset.

# ---------------------------------------------------------------------------
# Option string parsing
# ---------------------------------------------------------------------------

test_that("parse_options splits key=value pairs and types the values", {
  o <- parse_options("value_col=14;agg=max;fill=0")
  expect_equal(o$value_col, 14)          # numeric-looking values become numeric
  expect_equal(o$agg, "max")
  expect_equal(o$fill, 0)
})

test_that("parse_options turns comma lists into vectors", {
  o <- parse_options("keep_cols=TSS1,Count1,TSS2")
  expect_equal(o$keep_cols, c("TSS1", "Count1", "TSS2"))
})

test_that("parse_options keeps '=' inside a value", {
  o <- parse_options("key_regex=ID=([^;]+)")
  expect_equal(o$key_regex, "ID=([^;]+)")
})

test_that("parse_options returns an empty list for blank input", {
  expect_equal(parse_options(NA_character_), list())
  expect_equal(parse_options("   "), list())
})

test_that("parse_options rejects a malformed pair", {
  expect_error(parse_options("novalue"), "Malformed option")
})

# ---------------------------------------------------------------------------
# Column resolution
# ---------------------------------------------------------------------------

test_that("pick_col resolves by name, by index, and by candidate list", {
  dt <- data.table::data.table(a = 1, b = 2, c = 3)
  expect_equal(pick_col(dt, "b", "value"), "b")
  expect_equal(pick_col(dt, 3, "value"), "c")
  expect_equal(pick_col(dt, c("zz", "c", "a"), "key"), "c")   # first that exists wins
  expect_null(pick_col(dt, NULL, "value"))
})

test_that("pick_col errors informatively when the column is absent", {
  dt <- data.table::data.table(a = 1, b = 2)
  expect_error(pick_col(dt, "NoSuchColumn", "value"),
               "value column 'NoSuchColumn' not found")
  expect_error(pick_col(dt, 9, "value"), "index 9 exceeds 2 columns")
})

# ---------------------------------------------------------------------------
# Row filtering
# ---------------------------------------------------------------------------

test_that("apply_filter keeps == matches and drops NA rows", {
  dt <- data.table::data.table(Orientation = c("sense", "antisense", NA), v = 1:3)
  expect_equal(apply_filter(dt, "Orientation==sense")$v, 1L)
  expect_equal(apply_filter(dt, "Orientation!=sense")$v, 2L)
})

test_that("apply_filter is a no-op without an expression", {
  dt <- data.table::data.table(a = 1:3)
  expect_equal(nrow(apply_filter(dt, NULL)), 3L)
})

test_that("apply_filter rejects an unparseable expression", {
  dt <- data.table::data.table(a = 1:3)
  expect_error(apply_filter(dt, "a > 1"), "Cannot parse filter")
})

# ---------------------------------------------------------------------------
# Reference GFF and key resolution
# ---------------------------------------------------------------------------

test_that("read_reference_gff keeps only the requested feature types", {
  f <- fx_features()
  expect_setequal(f$ID, c("Zm00001eb000010", "Zm00001eb000020", "LTRRT_1", "TE_homo_1",
                          "Zm00001eb999999", "Zm00001ebDUPA", "Zm00001ebDUPB"))
  expect_false("TE_homo_2" %in% f$ID)          # helitron excluded
})

test_that("read_reference_gff expands attributes and honours drop_attributes", {
  f <- fx_features()
  expect_equal(f$Classification[f$ID == "LTRRT_1"], "LTR/Gypsy")
  expect_equal(f$Parent[f$ID == "LTRRT_1"], "repeat_region_1")
  expect_true(is.na(f$Parent[f$ID == "Zm00001eb000010"]))
  expect_false(any(c("tsd", "biotype", "logic_name") %in% names(f)))
})

test_that("read_reference_gff keeps coordinate-sharing features apart by default", {
  f <- fx_features()
  expect_true(all(c("Zm00001ebDUPA", "Zm00001ebDUPB") %in% f$ID))
})

test_that("legacy_collapse merges coordinate-sharing features into one row", {
  f <- read_reference_gff(fx_reference_gff(), c("^gene$", "LTR_retrotransposon"),
                          legacy_collapse = TRUE)
  expect_false("Zm00001ebDUPA" %in% f$ID)
  expect_true("Zm00001ebDUPA, Zm00001ebDUPB" %in% f$ID)
})

test_that("build_keymap resolves both the feature ID and its Parent", {
  km <- fx_keymap()
  expect_equal(km$ID[km$feat_key == "LTRRT_1"], "LTRRT_1")
  expect_equal(km$ID[km$feat_key == "repeat_region_1"], "LTRRT_1")
  expect_equal(km$ID[km$feat_key == "Zm00001eb000010"], "Zm00001eb000010")
})

# ---------------------------------------------------------------------------
# Preset: intersect_gff_chip
# ---------------------------------------------------------------------------

test_that("intersect_gff_chip reads column 14 and takes the maximum per feature", {
  row <- manifest_row(sample_id = "ChIP1", platform = "chip", role = "score",
                      preset = "intersect_gff_chip", path = fx_chip_gff())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000010"], 42)   # not 7, not 62
  expect_equal(res$score$value[res$score$ID == "LTRRT_1"], 5)
  expect_null(res$extras)
})

test_that("intersect_gff_chip defaults to a fill of 0", {
  row <- manifest_row(sample_id = "ChIP1", platform = "chip", role = "score",
                      preset = "intersect_gff_chip", path = fx_chip_gff())
  expect_equal(as.numeric(reader_spec(row)$fill), 0)
})

# ---------------------------------------------------------------------------
# Preset: intersect_gff_umr
# ---------------------------------------------------------------------------

test_that("intersect_gff_umr reads column 16 and takes the first value per feature", {
  row <- manifest_row(sample_id = "UMR", platform = "umr", role = "score",
                      preset = "intersect_gff_umr", path = fx_umr_gff())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000010"], 4.5)  # not 88.0
})

test_that("intersect_gff_umr defaults to a fill of 100 (no data, not zero methylation)", {
  row <- manifest_row(sample_id = "UMR", platform = "umr", role = "score",
                      preset = "intersect_gff_umr", path = fx_umr_gff())
  expect_equal(as.numeric(reader_spec(row)$fill), 100)
})

# ---------------------------------------------------------------------------
# Preset: count_matrix
# ---------------------------------------------------------------------------

test_that("count_matrix selects the sample's column and sums duplicate rows", {
  row <- manifest_row(sample_id = "SampleA", platform = "illumina", role = "score",
                      preset = "count_matrix", path = fx_count_matrix(), column = "SampleA")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000010"], 10)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000020"], 7)   # 3 + 4
})

test_that("count_matrix reads a different sample from the same file", {
  row <- manifest_row(sample_id = "SampleB", platform = "illumina", role = "score",
                      preset = "count_matrix", path = fx_count_matrix(), column = "SampleB")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000020"], 7)   # 2 + 5
  expect_equal(res$score$value[res$score$ID == "LTRRT_1"], 0)
})

test_that("count_matrix errors when the named column is absent", {
  row <- manifest_row(sample_id = "Nope", platform = "illumina", role = "score",
                      preset = "count_matrix", path = fx_count_matrix(), column = "NoSuchColumn")
  expect_error(read_keyed_table(row, fx_keymap(), verbose = FALSE),
               "value column 'NoSuchColumn' not found")
})

# ---------------------------------------------------------------------------
# Preset: tss_summary_v1
# ---------------------------------------------------------------------------

test_that("tss_summary_v1 keys a gene table on 'Gene' and scores Total_Reads", {
  row <- manifest_row(sample_id = "PB_Ear", platform = "pacbio", role = "score",
                      preset = "tss_summary_v1", path = fx_tss_v1_gene(), scope = "gene")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$key_col, "Gene")
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000010"], 25)
})

test_that("tss_summary_v1 resolves an LTR table keyed on the Parent id", {
  row <- manifest_row(sample_id = "PB_Ear", platform = "pacbio", role = "score",
                      preset = "tss_summary_v1", path = fx_tss_v1_ltr(), scope = "te")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$key_col, "Feature")
  expect_equal(res$score$ID, "LTRRT_1")          # repeat_region_1 resolved to the feature ID
  expect_equal(res$score$value, 12)
})

test_that("tss_summary_v1 carries the TSS columns through as extras", {
  row <- manifest_row(sample_id = "PB_Ear", platform = "pacbio", role = "score",
                      preset = "tss_summary_v1", path = fx_tss_v1_gene(), scope = "gene")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(names(res$extras),
               c("ID", paste("pacbio.PB_Ear.gene", c("TSS1", "Count1", "TSS2", "Count2"), sep = ".")))
  expect_equal(res$extras[["pacbio.PB_Ear.gene.Count1"]][res$extras$ID == "Zm00001eb000010"], 20)
})

# ---------------------------------------------------------------------------
# Preset: tss_summary_v2 (current IsoClassifier)
# ---------------------------------------------------------------------------

test_that("tss_summary_v2 keeps only sense rows", {
  row <- manifest_row(sample_id = "ONT1", platform = "ont", role = "score",
                      preset = "tss_summary_v2", path = fx_tss_v2())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_setequal(res$score$ID, c("Zm00001eb000010", "LTRRT_1"))
  expect_false("Zm00001eb000020" %in% res$score$ID)     # the antisense row
})

test_that("tss_summary_v2 resolves genes and structural LTR-RTs from one file", {
  row <- manifest_row(sample_id = "ONT1", platform = "ont", role = "score",
                      preset = "tss_summary_v2", path = fx_tss_v2())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000010"], 25)
  expect_equal(res$score$value[res$score$ID == "LTRRT_1"], 12)
})

test_that("tss_summary_v2 carries Class through as an extra", {
  row <- manifest_row(sample_id = "ONT1", platform = "ont", role = "score",
                      preset = "tss_summary_v2", path = fx_tss_v2())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$extras[["ont.ONT1.Class"]][res$extras$ID == "LTRRT_1"], "LTR_structural")
})

# ---------------------------------------------------------------------------
# Preset: isoforms_v1 / isoforms_v2 (extras only)
# ---------------------------------------------------------------------------

test_that("isoforms_v1 keys on Parent= inside the attrs column", {
  row <- manifest_row(sample_id = "ONT1", platform = "ont", role = "extra",
                      preset = "isoforms_v1", path = fx_isoforms_v1(), scope = "te")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_null(res$score)
  expect_equal(res$extras$ID, "LTRRT_1")
  expect_equal(res$extras[["ont.ONT1.te.ltr_left_reads"]], 4)
})

test_that("isoforms_v2 keeps every column but the structural ones", {
  row <- manifest_row(sample_id = "ONT1", platform = "ont", role = "extra",
                      preset = "isoforms_v2", path = fx_isoforms_v2())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  kept <- sub("^ont\\.ONT1\\.", "", setdiff(names(res$extras), "ID"))
  expect_true(all(c("Class", "Total_Reads", "spliced_reads") %in% kept))
  expect_false(any(c("Chrom", "Start", "End", "Orientation", "Attrs") %in% kept))
})

# ---------------------------------------------------------------------------
# Preset: motif_tsv
# ---------------------------------------------------------------------------

test_that("motif_tsv reads *_present flags written as 1/0 as logicals", {
  row <- manifest_row(sample_id = "M1", platform = "motif", role = "extra",
                      preset = "motif_tsv", path = fx_motif(), scope = "gene")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  present <- res$extras[["motif.M1.gene.TA_rich_present"]]
  expect_type(present, "logical")
  expect_equal(present[res$extras$ID == "Zm00001eb000010"], TRUE)
  expect_equal(present[res$extras$ID == "Zm00001eb000020"], FALSE)
})

test_that("motif_tsv keeps numeric motif distances numeric", {
  row <- manifest_row(sample_id = "M1", platform = "motif", role = "extra",
                      preset = "motif_tsv", path = fx_motif(), scope = "gene")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$extras[["motif.M1.gene.Dist_TSS_to_CCAAT"]][res$extras$ID == "Zm00001eb000010"], 30)
})

# ---------------------------------------------------------------------------
# Preset: keyed_tsv and option overrides
# ---------------------------------------------------------------------------

test_that("keyed_tsv scores column 2 keyed on column 1", {
  row <- manifest_row(sample_id = "CAGE", platform = "cage", role = "score",
                      preset = "keyed_tsv", path = fx_keyed())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "LTRRT_1"], 22)
})

test_that("keyed_tsv with keep_cols=* returns every non-key column as extras", {
  row <- manifest_row(sample_id = "CAGE", platform = "cage", role = "extra",
                      preset = "keyed_tsv", path = fx_keyed(), options = "keep_cols=*")
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_setequal(setdiff(names(res$extras), "ID"),
                  c("cage.CAGE.CAGE_dTSS", "cage.CAGE.CAGE_Shape"))
})

test_that("manifest options override the preset defaults", {
  row <- manifest_row(sample_id = "ChIP1", platform = "chip", role = "score",
                      preset = "intersect_gff_chip", path = fx_chip_gff(),
                      options = "agg=sum;fill=3")
  spec <- reader_spec(row)
  expect_equal(spec$agg, "sum")
  expect_equal(spec$fill, 3)
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$score$value[res$score$ID == "Zm00001eb000010"], 62)   # 7 + 42 + 13
})

test_that("keep_cols naming an absent column is an error, not a silent drop", {
  row <- manifest_row(sample_id = "M1", platform = "motif", role = "extra",
                      preset = "motif_tsv", path = fx_motif(),
                      options = "keep_cols=tss_abs,NotAColumn")
  expect_error(read_keyed_table(row, fx_keymap(), verbose = FALSE),
               "keep_cols not found")
})

test_that("read_keyed_table reports how many keys resolved to a feature", {
  row <- manifest_row(sample_id = "ChIP1", platform = "chip", role = "score",
                      preset = "intersect_gff_chip", path = fx_chip_gff())
  res <- read_keyed_table(row, fx_keymap(), verbose = FALSE)
  expect_equal(res$n_keys, 2L)          # Zm00001eb000010 and LTRRT_1
  expect_equal(res$n_resolved, 2L)
})

test_that("read_keyed_table errors on a missing file", {
  row <- manifest_row(sample_id = "X", platform = "chip", role = "score",
                      preset = "intersect_gff_chip", path = file.path(FIXTURE_DIR, "absent.gff"))
  expect_error(read_keyed_table(row, fx_keymap(), verbose = FALSE), "File not found")
})

# ---------------------------------------------------------------------------
# Manifest validation
# ---------------------------------------------------------------------------

test_that("read_manifest derives a row_id that includes the scope when it narrows the file", {
  m <- make_manifest(
    list(sample_id = "PB_Ear", platform = "pacbio", preset = "tss_summary_v1",
         path = fx_tss_v1_gene(), scope = "gene"),
    list(sample_id = "PB_Ear", platform = "pacbio", preset = "tss_summary_v1",
         path = fx_tss_v1_ltr(), scope = "te"),
    list(sample_id = "ONT1", platform = "ont", preset = "tss_summary_v2",
         path = fx_tss_v2(), scope = "both"))
  expect_equal(m$row_id, c("pacbio.PB_Ear.gene", "pacbio.PB_Ear.te", "ont.ONT1"))
})

test_that("read_manifest rejects an unknown preset, role, or scope", {
  expect_error(make_manifest(list(sample_id = "A", platform = "ont", preset = "nope",
                                  path = fx_tss_v2())), "Unknown preset")
  expect_error(make_manifest(list(sample_id = "A", platform = "ont", role = "maybe",
                                  preset = "tss_summary_v2", path = fx_tss_v2())),
               "role must be")
  expect_error(make_manifest(list(sample_id = "A", platform = "ont", preset = "tss_summary_v2",
                                  path = fx_tss_v2(), scope = "exon")), "scope must be")
})

test_that("read_manifest rejects two score rows that would claim the same assay column", {
  expect_error(make_manifest(
    list(sample_id = "ONT1", platform = "ont", preset = "tss_summary_v2", path = fx_tss_v2()),
    list(sample_id = "ONT1", platform = "ont", preset = "tss_summary_v2", path = fx_tss_v2())),
    "Duplicate platform/sample/scope among score rows")
})

test_that("read_manifest allows extras rows to share a row_id with a score row", {
  m <- make_manifest(
    list(sample_id = "ONT1", platform = "ont", role = "score", preset = "tss_summary_v2",
         path = fx_tss_v2()),
    list(sample_id = "ONT1", platform = "ont", role = "extra", preset = "isoforms_v2",
         path = fx_isoforms_v2()))
  expect_equal(nrow(m), 2L)
})
