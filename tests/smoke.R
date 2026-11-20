library(scRefAudit)

if (identical(Sys.getenv("SCREFAUDIT_RUN_R_TESTS"), "true")) {
  x <- demo_audit_data()
  fit <- audit_fit(x$cost, x$labels, epsilon = x$epsilon)
  stopifnot(inherits(fit, "scRefAudit_fit"), length(fit$classes) == 4L)
  ps <- rbind(fit$baseline_proportions, .9 * fit$baseline_proportions + .1 * c(.1,.2,.3,.4))
  a <- audit_compositions(fit, ps, order = 2L, tol = 1e-10,
                          batch_size = 1L, return_probabilities = TRUE)
  stopifnot(identical(dim(a$probabilities), c(96L, 2L, 4L)),
            max(a$marginal_residual) < 1.01e-10,
            identical(colnames(a$baseline_probabilities), fit$classes))
  # Named columns in a different order must be realigned correctly.
  b <- audit_compositions(fit, ps[,4:1,drop=FALSE], order=2L, tol=1e-10,
                          return_probabilities=TRUE)
  stopifnot(max(abs(a$probabilities-b$probabilities)) < 1e-7)
  p <- propose_scenarios(fit, anchors = c(1L, 3L), metric = "tv", radius = .01)
  p1 <- propose_scenarios(fit, anchors = 1L, metric = "tv", radius = .01)
  stopifnot(identical(p$anchor_index, c(1L,3L)), identical(p1$anchor_index,1L),
            max(abs(p$actual_tv-.01)) < 1e-12)
  r <- review_cells(fit, top_n=10L)
  stopifnot(nrow(r)==10L, all(r$cell_index>=1L), all(diff(r$linearized_radius)>=0))
  s <- stress_test(fit, budget=5L, metric="tv", radius=.01, seed=37L,
                   return_probabilities=TRUE)
  stopifnot(nrow(s$comparison)==2L, identical(dim(s$audit$probabilities),c(96L,10L,4L)),
            identical(s$directed$anchor_index, s$random$anchor_index))
  tmp <- tempfile(fileext=".rds"); saveRDS(s,tmp); stopifnot(is.list(readRDS(tmp))); unlink(tmp)
  cat("scRefAudit R interface smoke tests passed.\n")
} else {
  cat("R numerical tests not requested; set SCREFAUDIT_RUN_R_TESTS=true with a working Python backend.\n")
}
