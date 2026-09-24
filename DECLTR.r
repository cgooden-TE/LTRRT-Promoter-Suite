#!/usr/bin/env Rscript
#
# Multi-omic integration and activity classification for transposable elements and genes.
#
# Merges expression evidence from any number of Illumina, PacBio, and ONT samples with
# ChIP-seq peaks, DNA methylation (UMR), CAGE, and promoter-motif tables onto one reference
# annotation, estimates a per-sample expression threshold by segmented regression, fuses the
# evidence per tissue group, and assigns an activity label to every locus that clears at
# least one filter. Genes and TEs are scored by the same rules; the element class is carried
# through for interpretation but never used to branch.
#
# Nothing about file names or column names is assumed by the code. Which files exist, what
# each one contains, and which sample it belongs to are declared in a manifest; the biology
# (which samples form a tissue, which tissues are developmental) lives in a config.
#
# Feature model
# -------------
# Features come from the reference GFF, filtered by `feature_types` regexes on column 3.
# Input files may key on either a feature's ID or its Parent, because IsoClassifier and
# WindowScrubber report structural LTR-RTs under the EDTA Parent (repeat_region_N) while the
# annotation and read-count tables use the element ID (LTRRT_N). Both resolve to one locus.
#
# Scoring
# -------
# 1. Threshold units. Score rows sharing platform, tissue and replicate are summed, and a
#    segmented regression of loci-remaining against count threshold gives that unit's
#    breakpoint. Reported in <prefix>_breakpoints.csv.
# 2. Candidate filter. A locus is scored if its log1p total on any expression platform
#    clears that platform's breakpoint, or any ChIP column clears its own, or its
#    unmethylated signal reaches `umr_min_signal`. Everything else is left unlabelled.
# 3. Group evidence. Per tissue group and platform, the row median of log1p member columns
#    minus log1p of the group threshold passes through a sigmoid; platforms are then fused
#    by a weighted mean, skipping platforms with no data at that locus.
# 4. Labels. Groups are collapsed through `aliases`, then breadth of activity, dominance
#    margin, and the developmental-versus-vegetative means assign the label. Chromatin
#    support promotes weak-but-open loci to Repressed.
#
# Minimum inputs:
#     --manifest, TSV declaring one row per sample-file pairing (columns below)
#     --config, YAML declaring the reference GFF, tissue groups, and scoring parameters
#     --out, Output prefix (optional; falls back to output_prefix in the config)
#     --threads, data.table threads, default 1
#     --validate, Check inputs and exit without scoring
#
# Manifest columns (all required, blank where not applicable):
#     sample_id  : sample label used in output column names
#     platform   : illumina, pacbio, ont, chip, umr, cage, motif (free text; the expression
#                  platforms are whichever are listed under scoring.weights in the config)
#     role       : score (feeds the assay matrix) or extra (carried to the output only)
#     preset     : reader preset, one of intersect_gff_chip, intersect_gff_umr,
#                  count_matrix, tss_summary_v1, tss_summary_v2, isoforms_v1, isoforms_v2,
#                  motif_tsv, keyed_tsv
#     path       : input file, absolute or relative to the manifest
#     column     : for count_matrix, which column holds this sample
#     tissue     : tissue label used to build groups
#     replicate  : replicate label; platform + tissue + replicate defines a threshold unit
#     scope      : gene, te, or both, the features this file describes
#     options    : key=value;... overrides of the preset (key_col, key_regex, value_col,
#                  agg, fill, filter, keep_cols)
#
# Outputs (using the provided prefix):
#     - .qs : full data frame, annotation plus one column per sample, extras, and labels
#     - _labels.tsv : ID, coordinates, Passed_Platforms, Activity, and the score columns
#     - _breakpoints.csv : segmented-regression thresholds per unit, raw and log1p
#     - _manifest.tsv : copy of the manifest used, for provenance
#     - _config.yml : copy of the config used, for provenance
#
# Dependencies: data.table, yaml, qs, segmented, matrixStats (DECLTR-env)
#
# Usage (validate first, then run):
#     Rscript DECLTR.r \
#         --manifest configs/decltr_manifest.example.tsv \
#         --config configs/decltr_config.example.yml \
#         --validate
#
#     Rscript DECLTR.r \
#         --manifest configs/decltr_manifest.example.tsv \
#         --config configs/decltr_config.example.yml \
#         --out results/b73_run1 \
#         --threads 8
#
# Draft a manifest for a new dataset with scripts/make_manifest.sh, then fill in the
# tissue and replicate columns by hand. Unit tests: bash test/run_unit_tests.sh

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parse_args <- function(args) {
  out <- list(manifest = NULL, config = NULL, out = NULL, validate = FALSE, threads = 1L)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    take <- function() { if (i + 1 > length(args)) stop("Missing value for ", a); i <<- i + 1; args[i] }
    switch(a,
      "--manifest" = out$manifest <- take(),
      "--config"   = out$config <- take(),
      "--out"      = out$out <- take(),
      "--threads"  = out$threads <- as.integer(take()),
      "--validate" = out$validate <- TRUE,
      "-h" =, "--help" = { print_header(); quit(status = 0) },
      stop("Unknown argument: ", a))
    i <- i + 1
  }
  if (is.null(out$manifest) || is.null(out$config)) stop("--manifest and --config are required")
  out
}

# Print the comment block at the top of this file as the help text, so the two
# can never drift apart.
print_header <- function() {
  lines <- readLines(script_path())
  lines <- lines[-1]                                   # drop the shebang
  lines <- lines[seq_len(which(!startsWith(lines, "#"))[1] - 1)]
  cat(sub("^#[ ]?", "", lines), sep = "\n")
  cat("\n")
}

script_path <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(f[1]) else normalizePath("DECLTR.r")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
lib_dir <- file.path(dirname(script_path()), "decltr")
for (f in c("io.R", "model.R", "thresholds.R", "scoring.R", "labels.R", "validate.R"))
  source(file.path(lib_dir, f))
suppressPackageStartupMessages(library(qs))
data.table::setDTthreads(opt$threads)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
cfg <- read_config(opt$config)
manifest <- read_manifest(opt$manifest)
prefix <- if (!is.null(opt$out)) opt$out else cfg$output_prefix
if (is.null(prefix)) stop("Give --out or set output_prefix in the config")
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)

feats <- read_reference_gff(cfg$reference_gff, cfg$feature_types, cfg$drop_attributes,
                            legacy_collapse = isTRUE(cfg$legacy_collapse_duplicate_loci))

if (opt$validate) {
  v <- validate_inputs(manifest, cfg, feats)
  print_validation(v)
  quit(status = if (v$ok) 0 else 1)
}

model <- build_model(manifest, feats, assay_columns = cfg$assay_columns)
print(model)
model <- filter_features(model, cfg$drop_contigs_regex)
if ("Classification" %in% names(model$features)) {
  cl <- model$features$Classification
  model$features$Classification <- ifelse(is.na(cl), "Gene", cl)
}

# ---------------------------------------------------------------------------
# Thresholds and candidate filter
# ---------------------------------------------------------------------------
message("Estimating thresholds")
thr <- compute_thresholds(model, cfg)
keep <- Reduce(`|`, thr$pass) | thr$keep_chip | thr$keep_umr
idx <- which(keep)
message(length(idx), " of ", nrow(model$features), " features pass at least one platform or chromatin filter")

write.csv(thr$breakpoints, paste0(prefix, "_breakpoints.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Scores and labels
# ---------------------------------------------------------------------------
message("Scoring ", length(idx), " candidate features across ", length(thr$groups), " groups")
mats <- build_omics_matrices(model, thr$groups, idx, thr$expr_platforms)
butter <- compute_activity_scores_log(mats$omics_log, thr$groups, thr$group_thr,
                                      cfg$scoring$weights, s = cfg$scoring$s, loci = mats$loci)
chrom_support <- compute_chrom_support(thr, idx)
names(chrom_support) <- mats$loci
labels <- label_loci_from_activity(butter, chrom_support, aliases = cfg$aliases,
                                   dev_groups = cfg$dev_groups, veg_groups = cfg$veg_groups,
                                   lab = cfg$labels)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
message("Writing output to ", prefix, ".*")
wide <- model_to_wide(model, extras_fill = isTRUE(cfg$extras_fill))
wide$Passed_Platforms <- passed_platforms_string(thr)
wide <- merge(wide, labels, by = "ID", all.x = TRUE, sort = FALSE)

qsave(wide, paste0(prefix, ".qs"), preset = "balanced")
label_cols <- c("ID", "Chr", "Start", "End", "Strand", "Classification", "Passed_Platforms", names(labels)[-1])
label_cols <- intersect(label_cols, names(wide))
data.table::fwrite(wide[, label_cols], paste0(prefix, "_labels.tsv"), sep = "\t", na = "NA")
file.copy(opt$manifest, paste0(prefix, "_manifest.tsv"), overwrite = TRUE)
file.copy(opt$config, paste0(prefix, "_config.yml"), overwrite = TRUE)
print(table(wide$Activity, useNA = "ifany"))
message("Done.")
