#!/usr/bin/env Rscript
# Label-transition table between two DECLTR runs (e.g. before and after a scoring change).
# Usage: Rscript test/diff_decltr_labels.R old_labels.tsv new_labels.tsv
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2) stop("Usage: diff_decltr_labels.R old_labels.tsv new_labels.tsv")
o <- fread(a[1]); n <- fread(a[2])
j <- merge(o[, .(ID, old = Activity, old_pass = Passed_Platforms)],
           n[, .(ID, new = Activity, new_pass = Passed_Platforms)], by = "ID")
j[is.na(old), old := "Not scored"]; j[is.na(new), new := "Not scored"]
cat("Features compared:", nrow(j), " changed label:", sum(j$old != j$new), "\n\n")
tr <- j[old != new, .N, by = .(old, new)][order(-N)]
print(tr, nrows = 60)
cat("\nPassed_Platforms changes:\n")
pp <- j[!identical(old_pass, new_pass)][is.na(old_pass) != is.na(new_pass) | (old_pass != new_pass) %in% TRUE, .N, by = .(old_pass, new_pass)][order(-N)]
print(pp, nrows = 30)
