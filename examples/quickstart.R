# Run after installing scRefAudit. The first computation may set up Python.
library(scRefAudit)
x <- demo_audit_data()
fit <- audit_fit(x$cost, x$labels, epsilon = x$epsilon, cell_ids = x$cell_ids)
print(fit)

# 1% actual change in reference composition for both strategies.
checks <- stress_test(fit, budget = 16L, radius = 0.01, metric = "tv", seed = 37L)
print(checks$comparison)
print(review_cells(fit, top_n = 10L))

# Save ordinary R results, not the fit containing a live Python pointer.
saveRDS(checks, "scRefAudit_example_results.rds")
write.csv(review_cells(fit), "scRefAudit_review_cells.csv", row.names = FALSE)
