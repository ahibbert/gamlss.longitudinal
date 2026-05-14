#!/usr/bin/env Rscript
setwd("c:\\Users\\Aydin\\OneDrive - The University of Sydney (Students)\\gamlss.longitudinal")

results <- read.csv('results/coverage_suite/margin_copula_fit_results.csv')

cat('=== COVERAGE TEST RESULTS ===\n\n')
cat('Total combos:', nrow(results), '\n')
cat('Successful:', sum(results$success), '\n')
cat('Failed:', sum(!results$success), '\n')
cat('Success rate:', round(100*sum(results$success)/nrow(results), 1), '%\n\n')

# Copula summary
cat('=== RESULTS BY COPULA ===\n')
by_copula <- table(results$copula, results$success)
print(by_copula)

cat('\n=== TOP ERROR PATTERNS ===\n')
failed <- results[!results$success, ]
if (nrow(failed) > 0) {
  cat('Failed combos:', nrow(failed), '\n\n')
  error_types <- sort(table(failed$error), decreasing=TRUE)
  for (i in seq_len(min(10, length(error_types)))) {
    cat(sprintf("%d) %s (%d combos)\n", i, names(error_types)[i], error_types[i]))
  }
}

# Success metrics
cat('\n=== SUCCESS METRICS ===\n')
success_results <- results[results$success, ]
if (nrow(success_results) > 0) {
  cat('Mean absolute error:', 
      round(mean(success_results$mean_abs_error, na.rm=TRUE), 6), '\n')
  cat('Mean relative error:', 
      round(mean(success_results$mean_abs_rel_error, na.rm=TRUE), 6), '\n')
  cat('Max relative error:', 
      round(max(success_results$max_abs_rel_error, na.rm=TRUE), 6), '\n')
}
