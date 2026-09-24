# Per-tissue-group activity scores: log1p group medians -> sigmoid evidence -> weighted fusion.

sigmoid01_delta <- function(delta, s = 0.35) {
  1 / (1 + exp(-delta / s))
}

log1p_safe <- function(x) log1p(coalesce0(as.numeric(x)))

# One matrix per expression platform: rows = groups with that platform, cols = loci,
# values = row-median of log1p member columns.
build_omics_matrices <- function(model, groups, idx, expr_platforms) {
  loci <- model$features$ID[idx]
  omics <- list()
  for (p in expr_platforms) {
    gnames <- names(groups)[vapply(groups, function(g) !is.null(g[[p]]), logical(1))]
    if (!length(gnames)) next
    mat <- sapply(gnames, function(g) {
      cols <- groups[[g]][[p]]
      x <- model$assays[[p]][idx, cols, drop = FALSE]
      x <- log1p_safe(x); dim(x) <- c(length(idx), length(cols))
      matrixStats::rowMedians(x, na.rm = TRUE)
    })
    mat <- matrix(mat, ncol = length(gnames))
    mat <- t(mat)
    dimnames(mat) <- list(gnames, loci)
    omics[[p]] <- mat
  }
  list(loci = loci, omics_log = omics)
}

# Fuse per-platform evidence per group (weighted mean, NA-safe per locus).
compute_activity_scores_log <- function(omics_log, groups, group_thr, weights, s = 0.25, loci) {
  tissue_names <- names(groups)
  butter <- matrix(NA_real_, nrow = length(tissue_names), ncol = length(loci),
                   dimnames = list(tissue_names, loci))
  for (t in tissue_names) {
    num <- rep(0, length(loci))
    den <- rep(0, length(loci))
    for (p in names(omics_log)) {
      if (!t %in% rownames(omics_log[[p]]) || is.null(groups[[t]][[p]])) next
      thr <- group_thr[[t]][[p]]
      mid <- log1p(thr)
      delta <- as.numeric(omics_log[[p]][t, ]) - mid
      ev <- sigmoid01_delta(delta, s = s)
      w <- as.numeric(weights[[p]])
      ok <- is.finite(ev)
      num[ok] <- num[ok] + w * ev[ok]
      den[ok] <- den[ok] + w
    }
    out <- rep(NA_real_, length(loci))
    ok2 <- den > 0
    out[ok2] <- num[ok2] / den[ok2]
    butter[t, ] <- out
  }
  butter
}
