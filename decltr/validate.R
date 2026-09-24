# Manifest validation: files, columns, key resolution, and group membership.

validate_inputs <- function(manifest, cfg, feats) {
  problems <- character(0)
  missing <- !file.exists(manifest$path)
  if (any(missing))
    problems <- c(problems, paste0("row ", which(missing), " (", manifest$row_id[missing],
                                   "): file not found ", manifest$path[missing]))
  if (length(problems)) return(list(ok = FALSE, problems = problems, stats = NULL, model = NULL))

  model <- tryCatch(build_model(manifest, feats, verbose = TRUE, assay_columns = cfg$assay_columns),
                    error = function(e) { problems <<- c(problems, conditionMessage(e)); NULL })
  if (is.null(model)) return(list(ok = FALSE, problems = problems, stats = NULL, model = NULL))

  st <- model$read_stats
  # Intersect files legitimately carry many feature types outside the reference set, so a
  # low match rate is only fatal when essentially nothing resolves.
  none <- st$frac_resolved < 0.01 | is.na(st$frac_resolved)
  if (any(none))
    problems <- c(problems, paste0(st$row_id[none], ": only ",
                                   round(100 * st$frac_resolved[none]), "% of keys match a reference ID"))
  low <- !none & st$frac_resolved < 0.5
  if (any(low))
    message("Note: under half of the keys resolve for ", paste(st$row_id[low], collapse = ", "),
            " (expected for intersect files spanning other feature types)")

  groups <- tryCatch(resolve_groups(model, cfg),
                     error = function(e) { problems <<- c(problems, conditionMessage(e)); NULL })
  if (!is.null(groups)) {
    gs <- do.call(rbind, lapply(names(groups), function(g)
      data.frame(group = g, platform = names(groups[[g]]),
                 columns = vapply(groups[[g]], function(x) paste(x, collapse = ","), character(1)),
                 stringsAsFactors = FALSE)))
    rownames(gs) <- NULL
    attr(st, "groups") <- gs
    for (a in c("dev_groups", "veg_groups")) {
      absent <- setdiff(cfg[[a]], names(groups))
      if (length(absent)) problems <- c(problems, paste0(a, " names not in groups: ", paste(absent, collapse = ", ")))
    }
  }
  list(ok = !length(problems), problems = problems, stats = st, model = model)
}

print_validation <- function(v) {
  if (!is.null(v$stats)) {
    cat("\nPer-file key resolution:\n")
    st <- v$stats
    st$frac_resolved <- round(st$frac_resolved, 3)
    print(st, row.names = FALSE, right = FALSE)
    gs <- attr(v$stats, "groups")
    if (!is.null(gs)) { cat("\nResolved groups:\n"); print(gs, row.names = FALSE, right = FALSE) }
  }
  if (length(v$problems)) {
    cat("\nPROBLEMS:\n"); cat(paste0("  - ", v$problems), sep = "\n")
  } else cat("\nManifest OK.\n")
}
