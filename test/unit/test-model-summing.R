# Unit tests for the DECLTR data model: how manifest rows become assay columns,
# how the score rows of one sample are combined, and how groups and output are built.

# A manifest covering one PacBio sample split across a gene table and an LTR table,
# which is the shape that made the pre-refactor group median span two columns.
split_sample_manifest <- function() {
  make_manifest(
    list(sample_id = "PB_Ear", platform = "pacbio", role = "score", preset = "tss_summary_v1",
         path = fx_tss_v1_gene(), tissue = "Ear", replicate = "1", scope = "gene"),
    list(sample_id = "PB_Ear", platform = "pacbio", role = "score", preset = "tss_summary_v1",
         path = fx_tss_v1_ltr(), tissue = "Ear", replicate = "1", scope = "te"))
}

# ---------------------------------------------------------------------------
# assay_columns: per_sample vs per_row
# ---------------------------------------------------------------------------

test_that("per_sample merges the gene and LTR tables of one sample into one column", {
  m <- build_model(split_sample_manifest(), fx_features(), verbose = FALSE,
                   assay_columns = "per_sample")
  expect_equal(colnames(m$assays$pacbio), "pacbio.PB_Ear")
  v <- m$assays$pacbio[, "pacbio.PB_Ear"]
  expect_equal(unname(v["Zm00001eb000010"]), 25)   # from the gene table
  expect_equal(unname(v["LTRRT_1"]), 12)           # from the LTR table
  expect_equal(unname(v["TE_homo_1"]), 0)          # in neither, filled
})

test_that("per_row keeps one column per manifest row (pre-refactor behaviour)", {
  m <- build_model(split_sample_manifest(), fx_features(), verbose = FALSE,
                   assay_columns = "per_row")
  expect_setequal(colnames(m$assays$pacbio), c("pacbio.PB_Ear.gene", "pacbio.PB_Ear.te"))
  # Each locus carries its count in one column and the fill in the other. Taking a
  # row median across the pair is what halved long-read evidence before the fix.
  expect_equal(unname(m$assays$pacbio["LTRRT_1", "pacbio.PB_Ear.gene"]), 0)
  expect_equal(unname(m$assays$pacbio["LTRRT_1", "pacbio.PB_Ear.te"]), 12)
  expect_equal(median(m$assays$pacbio["LTRRT_1", ]), 6)
})

test_that("build_model rejects an unknown assay_columns setting", {
  expect_error(build_model(split_sample_manifest(), fx_features(), verbose = FALSE,
                           assay_columns = "per_thing"),
               "assay_columns must be")
})

# ---------------------------------------------------------------------------
# Combining rules by fill value
# ---------------------------------------------------------------------------

test_that("zero-fill score rows of one sample are summed", {
  # Two count files for one sample: absent means zero, so overlap is a real sum.
  mf <- make_manifest(
    list(sample_id = "SampleA", platform = "illumina", role = "score", preset = "count_matrix",
         path = fx_count_matrix(), column = "SampleA", tissue = "Ear", scope = "gene"),
    list(sample_id = "SampleA", platform = "illumina", role = "score", preset = "count_matrix",
         path = fx_count_matrix(), column = "SampleB", tissue = "Ear", scope = "te"))
  m <- build_model(mf, fx_features(), verbose = FALSE, assay_columns = "per_sample")
  v <- m$assays$illumina[, "illumina.SampleA"]
  expect_equal(unname(v["Zm00001eb000010"]), 11)   # 10 from SampleA + 1 from SampleB
  expect_equal(unname(v["Zm00001eb000020"]), 14)   # (3+4) + (2+5)
})

test_that("non-zero-fill rows covering disjoint features take the non-fill value", {
  # UMR fills 100 for "no data". A second row covering other features must not
  # be summed, or an unmethylated locus would end up above 100.
  umr_other <- write_tsv_lines("umr_other.gff", c(
    tsv("chr1", "EDTA", "gene", "300", "400", ".", "-", ".", "ID=Zm00001eb000020",
        "chr1", "290", "410", "UMR", "35", "36", "12.5", "6.25")))
  mf <- make_manifest(
    list(sample_id = "UMR", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = fx_umr_gff(), scope = "gene"),
    list(sample_id = "UMR", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = umr_other, scope = "te"))
  m <- build_model(mf, fx_features(), verbose = FALSE, assay_columns = "per_sample")
  v <- m$assays$umr[, "umr.UMR"]
  expect_equal(unname(v["Zm00001eb000010"]), 4.5)    # only in the first file
  expect_equal(unname(v["Zm00001eb000020"]), 12.5)   # only in the second, not 112.5
  expect_equal(unname(v["TE_homo_1"]), 100)          # in neither, stays at the fill
})

test_that("non-zero-fill rows overlapping on a feature are refused", {
  # Both files carry Zm00001eb000010, so there is no correct way to combine them.
  mf <- make_manifest(
    list(sample_id = "UMR", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = fx_umr_gff(), scope = "gene"),
    list(sample_id = "UMR", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = fx_umr_gff(), scope = "te"))
  expect_error(build_model(mf, fx_features(), verbose = FALSE, assay_columns = "per_sample"),
               "must cover disjoint features")
})

test_that("score rows of one sample declaring different fills are refused", {
  mf <- make_manifest(
    list(sample_id = "Mixed", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = fx_umr_gff(), scope = "gene"),
    list(sample_id = "Mixed", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = fx_umr_gff(), scope = "te", options = "fill=0"))
  expect_error(build_model(mf, fx_features(), verbose = FALSE, assay_columns = "per_sample"),
               "different fill values")
})

test_that("features absent from every input keep the preset fill", {
  mf <- make_manifest(
    list(sample_id = "UMR", platform = "umr", role = "score", preset = "intersect_gff_umr",
         path = fx_umr_gff(), scope = "both"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  expect_equal(unname(m$assays$umr["Zm00001eb000020", "umr.UMR"]), 100)
})

# ---------------------------------------------------------------------------
# Threshold units and row sums
# ---------------------------------------------------------------------------

test_that("unit_label groups a sample by platform, tissue and replicate", {
  s <- data.frame(platform = c("ont", "ont", "illumina"),
                  tissue = c("Leaf", "Leaf", NA),
                  replicate = c("CT1", "CT2", NA),
                  sample_id = c("ONT_CT1", "ONT_CT2", "IllA"),
                  stringsAsFactors = FALSE)
  expect_equal(unit_label(s), c("ont.Leaf_CT1", "ont.Leaf_CT2", "illumina.IllA"))
})

test_that("assay_rowsums totals only the named columns", {
  # Columns are c1 = (1, 2), c2 = (NA, 4), c3 = (Inf, 6).
  mat <- matrix(c(1, 2, NA, 4, Inf, 6), nrow = 2,
                dimnames = list(c("a", "b"), c("c1", "c2", "c3")))
  expect_equal(unname(assay_rowsums(mat, c("c1", "c2"))), c(1, 6))
})

test_that("assay_rowsums treats NA and Inf as zero", {
  mat <- matrix(c(1, 2, NA, 4, Inf, 6), nrow = 2,
                dimnames = list(c("a", "b"), c("c1", "c2", "c3")))
  expect_equal(unname(assay_rowsums(mat, colnames(mat))), c(1, 12))
})

# ---------------------------------------------------------------------------
# Feature filtering stays consistent across blocks
# ---------------------------------------------------------------------------

test_that("filter_features drops the same rows from features and every assay", {
  mf <- make_manifest(
    list(sample_id = "SampleA", platform = "illumina", role = "score", preset = "count_matrix",
         path = fx_count_matrix(), column = "SampleA", tissue = "Ear"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  expect_true("Zm00001eb999999" %in% m$features$ID)
  f <- filter_features(m, "scaf")
  expect_false("Zm00001eb999999" %in% f$features$ID)
  expect_equal(nrow(f$assays$illumina), nrow(f$features))
  expect_equal(rownames(f$assays$illumina), f$features$ID)
})

test_that("filter_features is a no-op without a regex", {
  mf <- make_manifest(
    list(sample_id = "SampleA", platform = "illumina", role = "score", preset = "count_matrix",
         path = fx_count_matrix(), column = "SampleA", tissue = "Ear"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  expect_equal(nrow(filter_features(m, NULL)$features), nrow(m$features))
})

# ---------------------------------------------------------------------------
# Group resolution
# ---------------------------------------------------------------------------

group_model <- function() {
  mf <- make_manifest(
    list(sample_id = "ONT_CT1", platform = "ont", role = "score", preset = "tss_summary_v2",
         path = fx_tss_v2(), tissue = "Leaf", replicate = "CT1"),
    list(sample_id = "ONT_CT2", platform = "ont", role = "score", preset = "tss_summary_v2",
         path = fx_tss_v2(), tissue = "Leaf", replicate = "CT2"),
    list(sample_id = "IllEar", platform = "illumina", role = "score", preset = "count_matrix",
         path = fx_count_matrix(), column = "SampleA", tissue = "Ear", replicate = "1"))
  build_model(mf, fx_features(), verbose = FALSE)
}

test_that("member_columns matches by tissue and by a tissue/replicate pair", {
  m <- group_model()
  expect_setequal(member_columns(m$samples, "ont", list("Leaf")),
                  c("ont.ONT_CT1", "ont.ONT_CT2"))
  expect_equal(member_columns(m$samples, "ont", list(list(tissue = "Leaf", replicate = "CT2"))),
               "ont.ONT_CT2")
  expect_equal(member_columns(m$samples, "ont", list(list(sample_id = "ONT_CT1"))),
               "ont.ONT_CT1")
})

test_that("resolve_groups follows an explicit config", {
  m <- group_model()
  cfg <- list(scoring = list(weights = list(ont = 0.7, illumina = 1)),
              groups = list(Leaf = list(ont = list("Leaf")),
                            Ear = list(illumina = list("Ear"))))
  g <- resolve_groups(m, cfg)
  expect_setequal(names(g), c("Leaf", "Ear"))
  expect_setequal(g$Leaf$ont, c("ont.ONT_CT1", "ont.ONT_CT2"))
  expect_equal(g$Ear$illumina, "illumina.IllEar")
})

test_that("resolve_groups falls back to one group per tissue per platform", {
  m <- group_model()
  cfg <- list(scoring = list(weights = list(ont = 0.7, illumina = 1)), groups = list())
  g <- resolve_groups(m, cfg)
  expect_setequal(names(g), c("Leaf", "Ear"))
  expect_setequal(g$Leaf$ont, c("ont.ONT_CT1", "ont.ONT_CT2"))
})

test_that("resolve_groups ignores platforms with no scoring weight", {
  m <- group_model()
  cfg <- list(scoring = list(weights = list(ont = 0.7)), groups = list())
  g <- resolve_groups(m, cfg)
  expect_null(g$Ear)              # illumina carries no weight, so it forms no group
  expect_setequal(names(g), "Leaf")
})

# ---------------------------------------------------------------------------
# Wide output
# ---------------------------------------------------------------------------

test_that("model_to_wide emits one column per assay column plus the extras", {
  mf <- make_manifest(
    list(sample_id = "PB_Ear", platform = "pacbio", role = "score", preset = "tss_summary_v1",
         path = fx_tss_v1_gene(), tissue = "Ear", scope = "gene"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  w <- model_to_wide(m)
  expect_true("pacbio.PB_Ear" %in% names(w))
  expect_true("pacbio.PB_Ear.gene.TSS1" %in% names(w))
  expect_equal(nrow(w), nrow(m$features))
})

test_that("model_to_wide fills missing extras with 0 and FALSE when asked", {
  mf <- make_manifest(
    list(sample_id = "M1", platform = "motif", role = "extra", preset = "motif_tsv",
         path = fx_motif(), scope = "gene"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  w <- model_to_wide(m, extras_fill = TRUE)
  i <- which(w$ID == "LTRRT_1")               # absent from the motif table
  expect_equal(w[["motif.M1.gene.tss_abs"]][i], 0)
  expect_equal(w[["motif.M1.gene.TA_rich_present"]][i], FALSE)
})

test_that("model_to_wide leaves missing extras as NA when fill is off", {
  mf <- make_manifest(
    list(sample_id = "M1", platform = "motif", role = "extra", preset = "motif_tsv",
         path = fx_motif(), scope = "gene"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  w <- model_to_wide(m, extras_fill = FALSE)
  expect_true(is.na(w[["motif.M1.gene.tss_abs"]][w$ID == "LTRRT_1"]))
})

test_that("model_to_wide disambiguates fields two extras blocks share", {
  # The v2 TSS summary and the v2 isoform table both carry Class for one sample.
  mf <- make_manifest(
    list(sample_id = "ONT1", platform = "ont", role = "score", preset = "tss_summary_v2",
         path = fx_tss_v2(), tissue = "Leaf"),
    list(sample_id = "ONT1", platform = "ont", role = "extra", preset = "isoforms_v2",
         path = fx_isoforms_v2(), tissue = "Leaf"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  w <- model_to_wide(m)
  expect_true("ont.ONT1.Class" %in% names(w))
  expect_true("ont.ONT1.isoforms_v2.Class" %in% names(w))
  expect_equal(anyDuplicated(names(w)), 0L)
})

# ---------------------------------------------------------------------------
# Read statistics
# ---------------------------------------------------------------------------

test_that("build_model records the key resolution rate per file", {
  m <- build_model(split_sample_manifest(), fx_features(), verbose = FALSE)
  expect_equal(nrow(m$read_stats), 2L)
  expect_true(all(m$read_stats$frac_resolved == 1))
})

test_that("unresolvable keys are counted, not silently dropped", {
  stray <- write_tsv_lines("stray_Gene_TSS.tsv", c(
    tsv("NotAFeature", "5", "1", "5", "NA", "0"),
    tsv("Zm00001eb000010", "7", "1", "7", "NA", "0")
  ), header = tsv("Gene", "Total_Reads", "TSS1", "Count1", "TSS2", "Count2"))
  mf <- make_manifest(
    list(sample_id = "S", platform = "ont", role = "score", preset = "tss_summary_v1",
         path = stray, tissue = "Leaf", scope = "gene"))
  m <- build_model(mf, fx_features(), verbose = FALSE)
  expect_equal(m$read_stats$n_keys, 2L)
  expect_equal(m$read_stats$n_resolved, 1L)
  expect_equal(m$read_stats$frac_resolved, 0.5)
})
