#!/usr/bin/env bash
# =============================================================================
# DECLTR unit tests
# =============================================================================
# Exercises the input layer (option parsing, column resolution, reference
# loading, key resolution, every reader preset) and the data model (assay column
# construction, the rules for combining a sample's score rows, group resolution,
# wide output). Fixtures are generated in a temp directory, so no test data files
# are needed and nothing outside test/unit_results/ is written.
#
# Prerequisites:
#   conda activate DECLTR-env
#   (or set DECLTR_RSCRIPT to an Rscript with testthat, data.table, yaml, qs)
#
# Usage:
#   bash test/run_unit_tests.sh
#
# Results:
#   test/unit_results/unit_test_results.txt   full test log
#   test/unit_results/unit_test_summary.csv   one row per test file
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE_DIR="$(dirname "$SCRIPT_DIR")"
RESULTS_DIR="$SCRIPT_DIR/unit_results"
RSCRIPT="${DECLTR_RSCRIPT:-Rscript}"

mkdir -p "$RESULTS_DIR"

echo "============================================="
echo " DECLTR Unit Tests"
echo "============================================="
echo "Test directory:   $SCRIPT_DIR/unit"
echo "Results directory: $RESULTS_DIR"
echo ""

if ! "$RSCRIPT" -e 'suppressMessages(library(testthat))' >/dev/null 2>&1; then
    echo "ERROR: testthat not available to '$RSCRIPT'."
    echo "       Activate DECLTR-env, or set DECLTR_RSCRIPT to an Rscript that has it."
    exit 1
fi

DECLTR_ROOT="$SUITE_DIR" "$RSCRIPT" - "$SCRIPT_DIR" "$RESULTS_DIR" <<'EOF'
suppressMessages(library(testthat))

args <- commandArgs(trailingOnly = TRUE)
test_dir_path <- file.path(args[1], "unit")
results_dir <- args[2]
txt_path <- file.path(results_dir, "unit_test_results.txt")
csv_path <- file.path(results_dir, "unit_test_summary.csv")

res <- test_dir(test_dir_path, reporter = "summary", stop_on_failure = FALSE)

df <- as.data.frame(res)
keep <- intersect(c("file", "test", "nb", "failed", "skipped", "error", "warning"), names(df))
write.csv(df[, keep, drop = FALSE], csv_path, row.names = FALSE)

n_tests <- nrow(df)
n_pass <- sum(df$nb, na.rm = TRUE)
n_failed <- sum(df$failed, na.rm = TRUE)
n_error <- sum(df$error, na.rm = TRUE)
n_skip <- sum(df$skipped, na.rm = TRUE)

# Written from the results rather than the console stream, so the log names every
# test rather than recording a row of progress dots.
status_of <- function(r) {
  if (isTRUE(r$error)) "ERROR" else if (r$failed > 0) "FAIL"
  else if (isTRUE(r$skipped)) "SKIP" else "ok"
}

lines <- c("DECLTR unit test results",
           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           paste("R:", R.version.string),
           paste("testthat:", as.character(utils::packageVersion("testthat"))),
           "")
for (f in unique(df$file)) {
  sub <- df[df$file == f, , drop = FALSE]
  lines <- c(lines, paste0(f, "  (", sum(sub$nb, na.rm = TRUE), " expectations in ",
                           nrow(sub), " tests)"), strrep("-", 78))
  for (i in seq_len(nrow(sub)))
    lines <- c(lines, sprintf("  %-5s %-64s %2d", status_of(sub[i, ]), sub$test[i], sub$nb[i]))
  lines <- c(lines, "")
}
lines <- c(lines, strrep("=", 78),
           sprintf("test blocks: %d   expectations: %d   failed: %d   errors: %d   skipped: %d",
                   n_tests, n_pass, n_failed, n_error, n_skip),
           sprintf("RESULT: %s", if (n_failed > 0 || n_error > 0) "FAILURES PRESENT" else "ALL PASSED"))
writeLines(lines, txt_path)

cat("\n---------------------------------------------\n")
cat(sprintf("test blocks: %d   expectations passed: %d   failed: %d   errors: %d   skipped: %d\n",
            n_tests, n_pass, n_failed, n_error, n_skip))
cat(sprintf("log:     %s\n", txt_path))
cat(sprintf("summary: %s\n", csv_path))

quit(status = if (n_failed > 0 || n_error > 0) 1 else 0)
EOF

echo ""
echo "============================================="
echo " Unit tests passed"
echo "============================================="
