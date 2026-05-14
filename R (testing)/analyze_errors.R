res <- read.csv('results/coverage_suite/margin_copula_fit_results.csv')
cat('Total rows:', nrow(res), '\n')
cat('Successful fits:', sum(res$success), '\n')
cat('Failed fits:', sum(!res$success), '\n\n')

# Categorize error types
errors <- res[!res$success, 'error']
error_types <- table(errors)
cat('Error types (sorted by frequency) - Top 15:\n')
top_errors <- sort(error_types, decreasing=TRUE)[1:15]
print(top_errors)

cat('\n\nError type labels and their frequencies:\n')
for (i in seq_along(top_errors)) {
  cat(sprintf('%2d. %s: %d occurrences\n', i, substr(names(top_errors)[i], 1, 70), top_errors[i]))
}

# Identify specific problematic combinations by error type
cat('\n\nExamples of major error categories:\n')
failed_rows <- which(!res$success)

# Group by copula
cat('\n--- Failures by Copula ---\n')
by_copula <- table(res[failed_rows, 'copula'])
print(by_copula)

# Specific patterns
cat('\n--- Specific Examples ---\n')
for (i in head(failed_rows, 8)) {
  cat(sprintf('%s + %s: %s\n', 
    res$margin[i], res$copula[i], substr(res$error[i], 1, 90)))
}
