# DECLTR input layer: manifest, config, reference GFF, and one generic keyed-table reader.
#
# Every input file is treated as a delimited table with a key column (optionally a regex
# applied to it), an optional value column that feeds a scoring assay, an aggregation rule
# for repeated keys, and pass-through columns kept as extras. Presets supply defaults per
# known format; the manifest's `options` column overrides them.

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

# key_col / value_col may be a column name or a 1-based index (for headerless files).
# keep_cols: character vector of extras columns, "*" for every non-key column, or NULL.
PRESETS <- list(
  intersect_gff_chip = list(header = FALSE, key_col = 9, key_regex = "ID=([^;]+)",
                            value_col = 14, agg = "max", fill = 0, keep_cols = NULL),
  intersect_gff_umr  = list(header = FALSE, key_col = 9, key_regex = "ID=([^;]+)",
                            value_col = 16, agg = "first", fill = 100, keep_cols = NULL),
  count_matrix       = list(header = TRUE, key_col = "Attributes", key_regex = "ID=([^;]+)",
                            value_col = NULL, agg = "sum", fill = 0, keep_cols = NULL),
  tss_summary_v1     = list(header = TRUE, key_col = c("Feature", "Gene"), key_regex = NULL,
                            value_col = "Total_Reads", agg = "first", fill = 0,
                            keep_cols = c("TSS1", "Count1", "TSS2", "Count2")),
  tss_summary_v2     = list(header = TRUE, key_col = "Feature", key_regex = NULL,
                            value_col = "Total_Reads", agg = "first", fill = 0,
                            filter = "Orientation==sense",
                            keep_cols = c("Class", "Top_Isoform", "TSS_Called", "TSS_Peak_Frac",
                                          "TSS1", "Count1", "TSS2", "Count2")),
  isoforms_v1        = list(header = TRUE, key_col = "attrs", key_regex = "Parent=([^;]+)",
                            value_col = NULL, agg = "first", fill = 0,
                            keep_cols = c("total_reads", "ltr_left_reads", "spliced_ltr_left",
                                          "ltr_right_reads", "spliced_ltr_right", "spanning_reads",
                                          "spliced_spanning", "ro5_reads", "spliced_ro5",
                                          "ro3_reads", "spliced_ro3")),
  isoforms_v2        = list(header = TRUE, key_col = "Feature", key_regex = NULL,
                            value_col = NULL, agg = "first", fill = 0,
                            filter = "Orientation==sense", keep_cols = "*",
                            drop_cols = c("Chrom", "Start", "End", "Name", "Strand",
                                          "Strand_Source", "Orientation", "Attrs")),
  motif_tsv          = list(header = TRUE, key_col = c("feature", "Feature", "ID"), key_regex = NULL,
                            value_col = NULL, agg = "first", fill = 0, keep_cols = "*"),
  keyed_tsv          = list(header = TRUE, key_col = 1, key_regex = NULL,
                            value_col = 2, agg = "first", fill = 0, keep_cols = NULL)
)

# ---------------------------------------------------------------------------
# Manifest and config
# ---------------------------------------------------------------------------

MANIFEST_COLS <- c("sample_id", "platform", "role", "preset", "path", "column",
                   "tissue", "replicate", "scope", "options")

# Parse "a=1;b=x,y" into a named list; comma-separated values become vectors.
#
# Splitting only on a ';' that starts the next 'key=' lets a value keep its own
# semicolons, which key_regex values need: ID=([^;]+) is one option, not two.
parse_options <- function(s) {
  if (is.na(s) || !nzchar(trimws(s))) return(list())
  parts <- strsplit(s, ";(?=[A-Za-z_][A-Za-z0-9_.]*=)", perl = TRUE)[[1]]
  out <- list()
  for (p in parts) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) stop("Malformed option '", p, "' (expected key=value)")
    key <- trimws(kv[1]); val <- trimws(paste(kv[-1], collapse = "="))
    if (grepl(",", val, fixed = TRUE)) val <- trimws(strsplit(val, ",", fixed = TRUE)[[1]])
    num <- suppressWarnings(as.numeric(val))
    out[[key]] <- if (length(val) == 1 && !is.na(num)) num else val
  }
  out
}

# Resolve a possibly relative path against a base directory.
resolve_path <- function(p, base_dir) {
  if (is.na(p) || !nzchar(p)) return(p)
  if (grepl("^(/|~)", p)) return(path.expand(p))
  normalizePath(file.path(base_dir, p), mustWork = FALSE)
}

read_manifest <- function(path) {
  m <- fread(path, sep = "\t", colClasses = "character", na.strings = c("", "NA"),
             fill = TRUE, quote = "")
  missing <- setdiff(MANIFEST_COLS, names(m))
  if (length(missing)) stop("Manifest is missing columns: ", paste(missing, collapse = ", "))
  m <- as.data.frame(m, stringsAsFactors = FALSE)[, MANIFEST_COLS]
  m$path <- vapply(m$path, resolve_path, character(1), base_dir = dirname(normalizePath(path)))
  m$role[is.na(m$role)] <- "extra"
  m$scope[is.na(m$scope)] <- "both"
  m$platform <- tolower(m$platform)
  bad <- !m$preset %in% names(PRESETS)
  if (any(bad)) stop("Unknown preset(s): ", paste(unique(m$preset[bad]), collapse = ", "))
  bad <- !m$role %in% c("score", "extra")
  if (any(bad)) stop("role must be 'score' or 'extra' (rows ", paste(which(bad), collapse = ","), ")")
  bad <- !m$scope %in% c("gene", "te", "both")
  if (any(bad)) stop("scope must be gene, te, or both (rows ", paste(which(bad), collapse = ","), ")")
  # A unique label per row: sample plus scope when the scope narrows the file.
  m$row_id <- ifelse(m$scope == "both", m$sample_id, paste(m$sample_id, m$scope, sep = "."))
  m$row_id <- paste(m$platform, m$row_id, sep = ".")
  # Score rows become assay columns, so their labels must be unique; extras may share one.
  sc <- m$row_id[m$role == "score"]
  dup <- duplicated(sc)
  if (any(dup)) stop("Duplicate platform/sample/scope among score rows: ",
                     paste(unique(sc[dup]), collapse = ", "))
  m
}

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  base_dir <- dirname(normalizePath(path))
  if (is.null(cfg$reference_gff)) stop("config needs reference_gff")
  cfg$reference_gff <- resolve_path(cfg$reference_gff, base_dir)
  if (!is.null(cfg$output_prefix)) cfg$output_prefix <- resolve_path(cfg$output_prefix, base_dir)
  defaults <- list(
    feature_types = c("^gene$", "LTR_retrotransposon"),
    drop_contigs_regex = "scaf",
    drop_attributes = character(0),
    groups = list(), aliases = list(),
    dev_groups = character(0), veg_groups = character(0),
    scoring = list(weights = list(illumina = 1, pacbio = 1, ont = 0.7), s = 0.25),
    labels = list(active_thr = 0.75, weak_thr = 0.45, silent_thr = 0.30, repress_thr = 0.35,
                  const_frac = 0.50, dom_margin = 0.12, min_facultative = 2L,
                  dev_veg_gap = 0.12, veg_cap = 0.55, veg_dev_gap = 0.12, dev_cap = 0.55),
    thresholds = list(default_max_k = 20L, chip_max_k = 50L, umr_min_signal = 0.1, seed = 1L,
                      pass_by = list(ont = "unit")),
    extras_fill = TRUE,
    assay_columns = "per_sample",
    legacy_collapse_duplicate_loci = FALSE
  )
  for (k in names(defaults)) {
    if (is.null(cfg[[k]])) cfg[[k]] <- defaults[[k]]
    else if (is.list(defaults[[k]]) && is.list(cfg[[k]]))
      for (kk in names(defaults[[k]])) if (is.null(cfg[[k]][[kk]])) cfg[[k]][[kk]] <- defaults[[k]][[kk]]
  }
  cfg
}

# ---------------------------------------------------------------------------
# Reference GFF -> features table
# ---------------------------------------------------------------------------

# Vectorised attribute expansion: every key present in the kept rows becomes a column.
# Repeated keys within one row are joined with ", " (mirrors the previous behaviour).
expand_gff_attributes <- function(attrs) {
  n <- length(attrs)
  pieces <- strsplit(attrs, ";", fixed = TRUE)
  lens <- lengths(pieces)
  long <- data.table(row = rep.int(seq_len(n), lens), kv = unlist(pieces, use.names = FALSE))
  long <- long[nzchar(kv)]
  long[, key := sub("=.*$", "", kv)]
  long[, value := sub("^[^=]*=?", "", kv)]
  long[, key := trimws(key)]
  wide <- dcast(long, row ~ key, value.var = "value",
                fun.aggregate = function(x) paste(x, collapse = ", "), fill = NA_character_)
  out <- data.table(row = seq_len(n))
  out <- merge(out, wide, by = "row", all.x = TRUE, sort = TRUE)
  out[, row := NULL]
  as.data.frame(out, stringsAsFactors = FALSE)
}

read_reference_gff <- function(path, feature_types, drop_attributes = character(0),
                               legacy_collapse = FALSE) {
  message("Reading reference GFF: ", path)
  g <- fread(path, sep = "\t", header = FALSE, quote = "", fill = TRUE,
             colClasses = "character", col.names = paste0("V", 1:9), showProgress = FALSE)
  g <- g[!startsWith(V1, "#")]
  keep <- Reduce(`|`, lapply(feature_types, function(p) grepl(p, g$V3)))
  g <- g[keep]
  message("  ", nrow(g), " features of type(s) ", paste(feature_types, collapse = "|"))
  attrs <- expand_gff_attributes(g$V9)
  if (!"ID" %in% names(attrs)) stop("Reference GFF features have no ID attribute")
  feats <- data.frame(ID = attrs$ID, Chr = g$V1, Type = g$V3, Start = g$V4, End = g$V5,
                      Strand = g$V7, stringsAsFactors = FALSE)
  other <- setdiff(names(attrs), c("ID", drop_attributes))
  feats <- cbind(feats, attrs[, other, drop = FALSE])
  if (!"Parent" %in% names(feats)) feats$Parent <- NA_character_
  if (legacy_collapse) feats <- collapse_duplicate_loci(feats)
  feats$ID <- trimws(feats$ID)
  if (anyDuplicated(feats$ID)) {
    d <- sum(duplicated(feats$ID))
    warning(d, " duplicate feature IDs in reference; keeping the first of each")
    feats <- feats[!duplicated(feats$ID), , drop = FALSE]
  }
  rownames(feats) <- NULL
  feats
}

# Reproduce the pre-refactor behaviour where features sharing Chr/Type/Start/End/Strand
# were merged into one row with attributes joined by ", " (so their IDs matched nothing).
# Only for regression against results produced before the refactor.
collapse_duplicate_loci <- function(feats) {
  dt <- as.data.table(feats)
  key <- c("Chr", "Type", "Start", "End", "Strand")
  n_before <- nrow(dt)
  dt <- dt[, lapply(.SD, function(x) if (length(x) == 1) x else paste(x, collapse = ", ")), by = key]
  setcolorder(dt, c("ID", key))
  message("  legacy collapse: ", n_before - nrow(dt), " duplicate-locus rows merged")
  as.data.frame(dt, stringsAsFactors = FALSE)
}

# Lookup from any key an input file may use (feature ID or its Parent) to the feature ID.
build_keymap <- function(feats) {
  km <- data.table(feat_key = feats$ID, ID = feats$ID)
  has_parent <- !is.na(feats$Parent) & nzchar(feats$Parent)
  if (any(has_parent))
    km <- rbind(km, data.table(feat_key = trimws(feats$Parent[has_parent]), ID = feats$ID[has_parent]))
  unique(km)
}

# ---------------------------------------------------------------------------
# Generic keyed-table reader
# ---------------------------------------------------------------------------

# Merge preset defaults, the manifest `column` field, and the options string.
reader_spec <- function(row) {
  spec <- PRESETS[[row$preset]]
  opts <- parse_options(row$options)
  for (k in names(opts)) spec[[k]] <- opts[[k]]
  if (!is.na(row$column) && nzchar(row$column)) spec$value_col <- row$column
  if (identical(spec$header, FALSE) && is.character(spec$key_col))
    spec$key_col <- as.integer(spec$key_col)
  spec
}

# Pick a column by name (first that exists among candidates) or 1-based index.
pick_col <- function(dt, col, what) {
  if (is.null(col)) return(NULL)
  if (is.numeric(col)) {
    if (col > ncol(dt)) stop(what, " column index ", col, " exceeds ", ncol(dt), " columns")
    return(names(dt)[col])
  }
  hit <- col[col %in% names(dt)]
  if (!length(hit)) stop(what, " column '", paste(col, collapse = "|"), "' not found; columns are: ",
                         paste(head(names(dt), 40), collapse = ", "))
  hit[1]
}

# Apply a simple "col==value" / "col!=value" filter expression.
apply_filter <- function(dt, expr) {
  if (is.null(expr) || !nzchar(expr)) return(dt)
  m <- regmatches(expr, regexec("^\\s*([^=!]+?)\\s*(==|!=)\\s*(.+?)\\s*$", expr))[[1]]
  if (length(m) != 4) stop("Cannot parse filter '", expr, "' (use col==value or col!=value)")
  col <- pick_col(dt, m[2], "filter")
  keep <- if (m[3] == "==") dt[[col]] == m[4] else dt[[col]] != m[4]
  dt[keep %in% TRUE]
}

# Read one manifest row. Returns list(score = data.table(ID, value) or NULL,
# extras = data.table(ID, ...) or NULL, n_keys, n_resolved).
read_keyed_table <- function(row, keymap, verbose = TRUE) {
  spec <- reader_spec(row)
  if (!file.exists(row$path)) stop("File not found: ", row$path)
  if (verbose) message("Reading ", row$row_id, ": ", basename(row$path))
  dt <- fread(row$path, sep = "\t", header = spec$header, quote = "", fill = TRUE,
              colClasses = "character", na.strings = c("", "NA"), showProgress = FALSE)
  dt <- apply_filter(dt, spec$filter)

  # Resolve every column reference by name now: later merges reorder columns, so
  # positional indexes would silently point at the wrong column afterwards.
  key_col <- pick_col(dt, spec$key_col, "key")
  vcol <- if (row$role == "score") pick_col(dt, spec$value_col, "value") else NULL
  keep <- spec$keep_cols
  if (!is.null(keep) && length(keep)) {
    if (identical(keep, "*")) {
      keep <- setdiff(names(dt), c(key_col, vcol, spec$drop_cols))
    } else {
      absent <- setdiff(keep, names(dt))
      if (length(absent)) stop("keep_cols not found in ", basename(row$path), ": ",
                               paste(absent, collapse = ", "))
    }
  }
  keys <- dt[[key_col]]
  if (!is.null(spec$key_regex)) {
    hit <- regexpr(spec$key_regex, keys, perl = TRUE)
    cap <- attr(hit, "capture.start")
    caplen <- attr(hit, "capture.length")
    keys <- ifelse(hit > 0, substr(keys, cap[, 1], cap[, 1] + caplen[, 1] - 1), NA_character_)
  }
  dt[, .key := trimws(keys)]
  dt <- dt[!is.na(.key) & nzchar(.key)]
  n_keys <- uniqueN(dt$.key)

  # Resolve keys to reference IDs (a key may map to several features, e.g. shared Parent).
  dt <- merge(dt, keymap, by.x = ".key", by.y = "feat_key", allow.cartesian = TRUE, sort = FALSE)
  n_resolved <- uniqueN(dt$.key)

  score <- NULL
  if (row$role == "score") {
    dt[, .value := suppressWarnings(as.numeric(get(vcol)))]
    score <- switch(spec$agg,
      max   = dt[, .(value = suppressWarnings(max(.value, na.rm = TRUE))), by = ID],
      sum   = dt[, .(value = sum(.value, na.rm = TRUE)), by = ID],
      first = dt[!duplicated(ID), .(ID, value = .value)],
      stop("Unknown agg '", spec$agg, "'"))
    score[!is.finite(value), value := NA_real_]
  }

  extras <- NULL
  if (!is.null(keep) && length(keep)) {
    extras <- dt[!duplicated(ID), c("ID", keep), with = FALSE]
    extras <- coerce_extras(extras)
    setnames(extras, keep, paste(row$row_id, keep, sep = "."))
  }

  list(score = score, extras = extras, n_keys = n_keys, n_resolved = n_resolved,
       n_rows = nrow(dt), key_col = key_col)
}

# Type extras: *_present -> logical from "1"/"0"/TRUE/FALSE, all-numeric columns -> numeric.
coerce_extras <- function(dt) {
  for (cn in setdiff(names(dt), "ID")) {
    x <- dt[[cn]]
    if (grepl("present", cn, ignore.case = TRUE)) {
      set(dt, j = cn, value = x %in% c("1", "TRUE", "True", "true"))
    } else {
      num <- suppressWarnings(as.numeric(x))
      if (all(is.na(num) == is.na(x))) set(dt, j = cn, value = num)
    }
  }
  dt
}
