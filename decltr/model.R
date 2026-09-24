# DECLTR data model: features + samples + per-platform assay matrices + extras.
#
# assays[[platform]] is a numeric matrix (features x columns). In the current
# version each manifest score row is one column, named by its row_id
# (platform.sample_id[.scope]); this preserves the previous behaviour where gene
# and TE counts of one long-read sample were separate columns.

suppressPackageStartupMessages(library(data.table))

# Build the model from a manifest and features table.
#
# assay_columns = "per_sample": every score row of one sample_id (e.g. the gene and TE
# tables of one long-read library) is summed into a single column, so a group median
# never spans a gene column and a TE column of the same sample. "per_row" keeps one
# column per manifest row, reproducing the pre-refactor behaviour.
build_model <- function(manifest, feats, verbose = TRUE, assay_columns = "per_sample") {
  if (!assay_columns %in% c("per_sample", "per_row"))
    stop("assay_columns must be 'per_sample' or 'per_row'")
  keymap <- build_keymap(feats)
  id_index <- setNames(seq_len(nrow(feats)), feats$ID)
  assays <- list()
  extras <- list()
  fills <- list()
  stats <- vector("list", nrow(manifest))
  manifest$assay_col <- NA_character_

  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    res <- read_keyed_table(row, keymap, verbose = verbose)
    stats[[i]] <- data.frame(row_id = row$row_id, path = basename(row$path),
                             n_keys = res$n_keys, n_resolved = res$n_resolved,
                             frac_resolved = if (res$n_keys) res$n_resolved / res$n_keys else NA_real_,
                             stringsAsFactors = FALSE)
    if (!is.null(res$score)) {
      spec <- reader_spec(row)
      v <- rep(as.numeric(spec$fill), nrow(feats))
      idx <- id_index[res$score$ID]
      v[idx] <- res$score$value
      v[is.na(v)] <- as.numeric(spec$fill)
      p <- row$platform
      cname <- if (assay_columns == "per_sample") paste(p, row$sample_id, sep = ".") else row$row_id
      manifest$assay_col[i] <- cname
      if (is.null(assays[[p]])) {
        assays[[p]] <- matrix(v, ncol = 1, dimnames = list(feats$ID, cname))
      } else if (cname %in% colnames(assays[[p]])) {
        fill_v <- as.numeric(spec$fill)
        if (!identical(fills[[cname]], fill_v))
          stop("Score rows of sample '", row$sample_id, "' declare different fill values (",
               fills[[cname]], " and ", fill_v, "); they cannot share one assay column")
        prev <- assays[[p]][, cname]
        if (fill_v == 0) {
          # Counts: absent means zero, so combining rows of one sample is a plain sum.
          assays[[p]][, cname] <- prev + v
        } else {
          # Non-zero fill marks "no data" rather than zero, so rows must cover disjoint
          # features; summing would add the sentinel into the result.
          overlap <- (prev != fill_v) & (v != fill_v)
          if (any(overlap))
            stop("Score rows of sample '", row$sample_id, "' both carry data at ", sum(overlap),
                 " feature(s) (e.g. ", paste(head(feats$ID[overlap], 3), collapse = ", "),
                 "). With a non-zero fill (", fill_v, ") they must cover disjoint features; ",
                 "give them different sample_id values or set assay_columns: per_row")
          prev[v != fill_v] <- v[v != fill_v]
          assays[[p]][, cname] <- prev
        }
      } else {
        assays[[p]] <- cbind(assays[[p]], v)
        colnames(assays[[p]])[ncol(assays[[p]])] <- cname
      }
      fills[[cname]] <- as.numeric(spec$fill)
    }
    if (!is.null(res$extras)) extras[[paste(row$row_id, row$preset, sep = "|")]] <- res$extras
  }

  samples <- manifest[manifest$role == "score", , drop = FALSE]
  samples$unit <- unit_label(samples)
  # One sample row per assay column; per_row keeps every manifest row.
  samples <- samples[!duplicated(samples$assay_col), , drop = FALSE]
  samples$row_id <- samples$assay_col
  structure(list(features = feats, samples = samples, assays = assays, extras = extras,
                 fills = fills, read_stats = do.call(rbind, stats)),
            class = "decltr_model")
}

# Threshold unit: one (platform, tissue, replicate) combination. Columns of a
# unit are summed before threshold estimation.
unit_label <- function(samples) {
  tissue <- ifelse(is.na(samples$tissue), "", samples$tissue)
  rep <- ifelse(is.na(samples$replicate), "", samples$replicate)
  u <- ifelse(nzchar(tissue) & nzchar(rep), paste(tissue, rep, sep = "_"),
              ifelse(nzchar(tissue), tissue, ifelse(nzchar(rep), rep, samples$sample_id)))
  paste(samples$platform, u, sep = ".")
}

print.decltr_model <- function(x, ...) {
  cat("DECLTR model:", nrow(x$features), "features\n")
  for (p in names(x$assays))
    cat(sprintf("  assay %-9s %d column(s): %s\n", p, ncol(x$assays[[p]]),
                paste(head(colnames(x$assays[[p]]), 6), collapse = ", ")))
  cat("  extras blocks:", length(x$extras), "\n")
  invisible(x)
}

# Drop features whose contig matches a regex (e.g. scaffolds) from every block.
filter_features <- function(model, drop_regex) {
  if (is.null(drop_regex) || !nzchar(drop_regex)) return(model)
  keep <- !grepl(drop_regex, model$features$Chr, ignore.case = TRUE)
  message("Dropping ", sum(!keep), " features on contigs matching '", drop_regex, "'")
  model$features <- model$features[keep, , drop = FALSE]
  for (p in names(model$assays)) model$assays[[p]] <- model$assays[[p]][keep, , drop = FALSE]
  model
}

# Column names of a platform's assay for samples matching a group member spec.
# A member is either a tissue name or a list with any of tissue / replicate / sample_id.
member_columns <- function(samples, platform, members) {
  s <- samples[samples$platform == platform, , drop = FALSE]
  if (!nrow(s)) return(character(0))
  hit <- rep(FALSE, nrow(s))
  for (m in members) {
    if (is.character(m)) m <- list(tissue = m)
    ok <- rep(TRUE, nrow(s))
    for (f in intersect(names(m), c("tissue", "replicate", "sample_id")))
      ok <- ok & (s[[f]] %in% m[[f]])
    hit <- hit | ok
  }
  s$row_id[hit]
}

# Expand config groups into platform -> columns per group. When the config has no
# explicit groups, one group per distinct tissue per platform is created.
resolve_groups <- function(model, cfg) {
  samples <- model$samples
  expr_platforms <- intersect(names(cfg$scoring$weights), names(model$assays))
  groups <- list()
  if (length(cfg$groups)) {
    for (g in names(cfg$groups)) {
      spec <- cfg$groups[[g]]
      groups[[g]] <- list()
      for (p in intersect(names(spec), expr_platforms))
        groups[[g]][[p]] <- member_columns(samples, p, spec[[p]])
      groups[[g]] <- groups[[g]][lengths(groups[[g]]) > 0]
    }
  } else {
    for (p in expr_platforms) {
      s <- samples[samples$platform == p & !is.na(samples$tissue), , drop = FALSE]
      for (t in unique(s$tissue)) groups[[t]][[p]] <- s$row_id[s$tissue == t]
    }
  }
  groups <- groups[lengths(groups) > 0]
  if (!length(groups)) stop("No tissue groups could be resolved from the manifest and config")
  groups
}

# Wide data.frame of features + assay columns + extras, for output.
model_to_wide <- function(model, extras_fill = TRUE) {
  out <- as.data.table(model$features)
  for (p in names(model$assays)) {
    m <- model$assays[[p]]
    for (j in seq_len(ncol(m))) set(out, j = colnames(m)[j], value = unname(m[, j]))
  }
  # Two blocks of one sample may share field names (e.g. Class in both the v2 TSS
  # summary and isoform table); the later block's colliding columns get the preset inserted.
  used <- names(out)
  for (b in names(model$extras)) {
    e <- model$extras[[b]]
    idx <- match(out$ID, e$ID)
    row_id <- sub("\\|.*$", "", b); preset <- sub("^.*\\|", "", b)
    for (cn0 in setdiff(names(e), "ID")) {
      cn <- cn0
      if (cn %in% used) cn <- paste(row_id, preset, sub(paste0("^", row_id, "\\."), "", cn0), sep = ".")
      used <- c(used, cn)
      v <- e[[cn0]][idx]
      if (extras_fill) {
        if (is.numeric(v)) v[is.na(v)] <- 0
        if (is.logical(v)) v[is.na(v)] <- FALSE
      }
      set(out, j = cn, value = v)
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
