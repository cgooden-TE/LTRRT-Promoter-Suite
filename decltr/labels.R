# Rule-based activity labels from the group x locus activity matrix.

label_loci_from_activity <- function(butter, chrom_support = NULL, aliases = list(),
                                     dev_groups = character(0), veg_groups = character(0),
                                     lab = list()) {
  active_thr    <- lab$active_thr
  weak_thr      <- lab$weak_thr
  silent_thr    <- lab$silent_thr
  repress_thr   <- lab$repress_thr
  const_frac    <- lab$const_frac
  dev_veg_gap   <- lab$dev_veg_gap
  veg_cap       <- lab$veg_cap
  veg_dev_gap   <- lab$veg_dev_gap
  dev_cap       <- lab$dev_cap
  dom_margin    <- lab$dom_margin
  min_facultative <- as.integer(lab$min_facultative)

  loci        <- colnames(butter)
  tissues_raw <- rownames(butter)

  # Collapse aliased groups (e.g. leaf sections and ONT replicates) to one label-level tissue.
  tissues_alias <- tissues_raw
  hit <- tissues_alias %in% names(aliases)
  tissues_alias[hit] <- unlist(aliases[tissues_alias[hit]], use.names = FALSE)
  alias_levels <- unique(tissues_alias)

  butter_alias <- sapply(alias_levels, function(tt) {
    rows <- which(tissues_alias == tt)
    if (length(rows) == 1L) return(butter[rows, ])
    apply(butter[rows, , drop = FALSE], 2, max, na.rm = TRUE)
  })
  butter_alias <- matrix(butter_alias, ncol = length(alias_levels))
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
