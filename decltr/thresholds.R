# Threshold estimation (segmented regression) and platform pass flags.

PLATFORM_LABELS <- c(ont = "ONT", pacbio = "PacBio", illumina = "Illumina")
platform_label <- function(p) ifelse(p %in% names(PLATFORM_LABELS), PLATFORM_LABELS[p], p)

coalesce0 <- function(x) { x[is.na(x)] <- 0; x }

# Breakpoint of loci-remaining vs threshold; unchanged from the original implementation.
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

# Row sums over a set of assay columns, with non-finite values treated as 0.
assay_rowsums <- function(mat, cols) {
  X <- mat[, cols, drop = FALSE]
  X[!is.finite(X)] <- 0
  rowSums(X)
}

# Compute unit thresholds, group thresholds, platform pass flags, and chromatin support.
compute_thresholds <- function(model, cfg) {
  set.seed(cfg$thresholds$seed)
  samples <- model$samples
  assays <- model$assays
  n <- nrow(model$features)
  expr_platforms <- intersect(names(cfg$scoring$weights), names(assays))
  max_k <- as.integer(cfg$thresholds$default_max_k)

  # --- unit thresholds: one (platform, tissue, replicate) unit = summed columns ---
  units <- list()
  for (p in expr_platforms) {
    s <- samples[samples$platform == p, , drop = FALSE]
    for (u in unique(s$unit)) {
      cols <- s$row_id[s$unit == u]
      sids <- unique(s$sample_id[s$unit == u])
      counts <- assay_rowsums(assays[[p]], cols)
      units[[u]] <- list(platform = p, cols = cols,
                         label = if (length(sids) == 1) sids else sub("^[^.]+\\.", "", u),
                         thr = estimate_seg_threshold(counts, max_k = max_k, default = 1L))
    }
  }
  unit_thr <- vapply(units, function(x) as.numeric(x$thr), numeric(1))

  # --- group thresholds: median of member unit thresholds ---
  groups <- resolve_groups(model, cfg)
  group_thr <- list()
  for (g in names(groups)) {
    group_thr[[g]] <- list()
    for (p in names(groups[[g]])) {
      us <- unique(samples$unit[samples$row_id %in% groups[[g]][[p]]])
      thr <- suppressWarnings(median(unit_thr[us], na.rm = TRUE))
      if (!is.finite(thr)) thr <- 1
      group_thr[[g]][[p]] <- thr
    }
  }

  # --- platform pass flags on log1p totals ---
  pass <- list()
  sum_thr <- list()
  for (p in expr_platforms) {
    mode <- cfg$thresholds$pass_by[[p]]
    if (is.null(mode)) mode <- "total"
    if (mode == "unit") {
      us <- names(units)[vapply(units, function(x) x$platform == p, logical(1))]
      flags <- lapply(us, function(u) {
        tot <- log1p(assay_rowsums(assays[[p]], units[[u]]$cols))
        thr <- estimate_seg_threshold(tot, max_k = max_k, default = 1L)
        sum_thr[[u]] <<- thr
        tot >= thr
      })
      pass[[p]] <- Reduce(`|`, flags)
    } else {
      tot <- log1p(assay_rowsums(assays[[p]], colnames(assays[[p]])))
      thr <- estimate_seg_threshold(tot, max_k = max_k, default = 1L)
      sum_thr[[p]] <- thr
      pass[[p]] <- tot >= thr
    }
  }

  # --- ChIP: per-column threshold, present if at or above ---
  chip_present <- NULL
  keep_chip <- rep(FALSE, n)
  chip_thr <- numeric(0)
  if (!is.null(assays$chip)) {
    m <- assays$chip
    chip_present <- matrix(FALSE, nrow = n, ncol = ncol(m), dimnames = dimnames(m))
    for (cc in colnames(m)) {
      chip_thr[[cc]] <- estimate_seg_threshold(m[, cc], max_k = as.integer(cfg$thresholds$chip_max_k), default = 1L)
      chip_present[, cc] <- coalesce0(m[, cc]) >= chip_thr[[cc]]
    }
    keep_chip <- rowSums(chip_present) > 0
  }

  # --- UMR: percent methylation -> unmethylated signal in [0, 1]; 100 (fill) means no data ---
  umr_signal <- rep(NA_real_, n)
  keep_umr <- rep(FALSE, n)
  if (!is.null(assays$umr)) {
    if (ncol(assays$umr) > 1) warning("Several UMR columns supplied; using the first: ", colnames(assays$umr)[1])
    v <- as.numeric(assays$umr[, 1])
    v[v == 100] <- NA
    umr_signal <- pmin(pmax(1 - (v / 100), 0), 1)
    keep_umr <- is.finite(umr_signal) & (umr_signal >= cfg$thresholds$umr_min_signal)
  }

  # --- breakpoint report (raw-count thresholds per unit) ---
  order_p <- c(intersect(names(PLATFORM_LABELS), expr_platforms),
               setdiff(expr_platforms, names(PLATFORM_LABELS)))
  bp <- do.call(rbind, lapply(order_p, function(p) {
    us <- names(units)[vapply(units, function(x) x$platform == p, logical(1))]
    if (!length(us)) return(NULL)
    data.frame(Platform = unname(platform_label(p)),
               Sample = vapply(us, function(u) units[[u]]$label, character(1)),
               Breakpoint_raw = unname(unit_thr[us]), stringsAsFactors = FALSE)
  }))
  if (!is.null(bp)) {
    bp$Breakpoint_log1p <- log1p(bp$Breakpoint_raw)
    rownames(bp) <- NULL
  }

  list(units = units, unit_thr = unit_thr, groups = groups, group_thr = group_thr,
       pass = pass, sum_thr = sum_thr, chip_thr = chip_thr, chip_present = chip_present,
       keep_chip = keep_chip, umr_signal = umr_signal, keep_umr = keep_umr,
       breakpoints = bp, expr_platforms = expr_platforms)
}

# "Illumina;PacBio;" style string of platforms each feature passed, NA when none.
passed_platforms_string <- function(thr) {
  n <- length(thr$keep_chip)
  order_p <- c(intersect(c("illumina", "pacbio", "ont"), thr$expr_platforms),
               setdiff(thr$expr_platforms, c("illumina", "pacbio", "ont")))
  s <- rep("", n)
  for (p in order_p) s <- paste0(s, ifelse(thr$pass[[p]], paste0(platform_label(p), ";"), ""))
  s[s == ""] <- NA
  s
}

# Fraction of ChIP columns present, or unmethylated signal, whichever is larger.
compute_chrom_support <- function(thr, idx) {
  chip_pass_frac <- if (!is.null(thr$chip_present)) rowMeans(thr$chip_present[idx, , drop = FALSE]) else rep(0, length(idx))
  umr <- thr$umr_signal[idx]
  umr[!is.finite(umr)] <- NA_real_
  cs <- pmax(chip_pass_frac, umr, na.rm = TRUE)
  cs[!is.finite(cs)] <- 0
  cs
}
