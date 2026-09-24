#!/usr/bin/env Rscript
# Regression check: compare a DECLTR run against a baseline qs and breakpoints CSV.
#
# Usage: Rscript test/compare_decltr_baseline.R baseline.qs baseline_breakpoints.csv new_prefix
#
# Compares the label and score columns by feature ID, and the breakpoint values by
# platform and order. Extras column names are not compared (they changed by design).

suppressPackageStartupMessages({ library(qs); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: compare_decltr_baseline.R baseline.qs baseline_breakpoints.csv new_prefix")

base <- as.data.table(qread(args[1]))
new  <- as.data.table(qread(paste0(args[3], ".qs")))
bp_base <- fread(args[2])
bp_new  <- fread(paste0(args[3], "_breakpoints.csv"))

cat("Baseline rows:", nrow(base), " new rows:", nrow(new), "\n")
cat("IDs only in baseline:", sum(!base$ID %in% new$ID), " (collapsed 'a, b' IDs:",
    sum(grepl(", ", base$ID)), ")\n")
cat("IDs only in new:", sum(!new$ID %in% base$ID), "\n")

cols <- c("Activity", "top_tissue", "top_score", "margin", "breadth_active",
          "dev_score", "veg_score", "Passed_Platforms")
j <- merge(base[, c("ID", cols), with = FALSE], new[, c("ID", cols), with = FALSE],
           by = "ID", suffixes = c(".base", ".new"))
cat("Joined rows:", nrow(j), "\n\n")

ok_all <- TRUE
for (cc in cols) {
  a <- j[[paste0(cc, ".base")]]; b <- j[[paste0(cc, ".new")]]
  if (is.numeric(a)) {
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & (a == b | abs(a - b) <= 1e-9 | (is.infinite(a) & is.infinite(b))))
  } else {
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  }
  n_diff <- sum(!same)
  cat(sprintf("%-18s differing: %d\n", cc, n_diff))
  if (n_diff) {
    ok_all <- FALSE
    if (!is.numeric(a)) print(head(sort(table(paste(a[!same], "->", b[!same])), decreasing = TRUE), 10))
    else print(head(j[!same, c("ID", paste0(cc, ".base"), paste0(cc, ".new")), with = FALSE], 5))
  }
}

cat("\nBreakpoints (baseline vs new, by platform and order):\n")
for (p in unique(bp_base$Platform)) {
  a <- bp_base[Platform == p]; b <- bp_new[Platform == p]
  cat(sprintf("  %-9s n=%d/%d  raw equal: %s\n", p, nrow(a), nrow(b),
              if (nrow(a) == nrow(b)) all(a$Breakpoint_raw == b$Breakpoint_raw) else "n mismatch"))
  if (nrow(a) == nrow(b) && !all(a$Breakpoint_raw == b$Breakpoint_raw)) {
    ok_all <- FALSE
    print(data.table(base = a$Sample, base_raw = a$Breakpoint_raw, new = b$Sample, new_raw = b$Breakpoint_raw))
  }
}
cat("\nRESULT:", if (ok_all) "IDENTICAL on compared columns" else "DIFFERENCES FOUND", "\n")
quit(status = if (ok_all) 0 else 1)
