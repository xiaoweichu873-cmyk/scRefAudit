.scRefAudit_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c("numpy>=1.26", "scipy>=1.11"),
                         python_version = ">=3.10")
}

.backend <- function() {
  if (is.null(.scRefAudit_env$module)) {
    path <- system.file("python", package = "scRefAudit", mustWork = TRUE)
    .scRefAudit_env$module <- reticulate::import_from_path(
      "screfaudit", path = path, convert = TRUE)
  }
  .scRefAudit_env$module
}

.scalar_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  as.integer(x)
}

.check_fit <- function(fit) {
  if (!inherits(fit, "scRefAudit_fit")) {
    stop("fit must be returned by audit_fit().", call. = FALSE)
  }
  if (reticulate::py_is_null_xptr(fit$engine)) {
    stop("This fit belongs to an earlier R session. Re-run audit_fit(); save results, not the Python fit pointer.",
         call. = FALSE)
  }
  invisible(fit)
}

.columns_frame <- function(x) {
  as.data.frame(lapply(x, as.vector), stringsAsFactors = FALSE,
                optional = TRUE, check.names = FALSE)
}

.decorate_audit <- function(x, fit) {
  x$cells <- .columns_frame(x$cells)
  dimnames(x$proportions) <- list(NULL, fit$classes)
  dimnames(x$baseline_probabilities) <- list(fit$cell_ids, fit$classes)
  if (!is.null(x$probabilities)) {
    dimnames(x$probabilities) <- list(fit$cell_ids, NULL, fit$classes)
  }
  x$diagnostics$fallback_indices <- as.integer(x$diagnostics$fallback_indices) + 1L
  x
}

.decorate_proposals <- function(x, fit) {
  x$anchor_index <- as.integer(x$anchor_index) + 1L
  dimnames(x$proportions) <- list(NULL, fit$classes)
  dimnames(x$directions) <- list(NULL, fit$classes)
  x
}

#' Fit the fixed-cost reference annotation model
#' @export
audit_fit <- function(cost, labels, epsilon = 0.05, a = NULL, b = NULL,
                      cell_ids = rownames(cost), base_tol = 1e-11) {
  if (!is.matrix(cost) || !is.numeric(cost) || any(!is.finite(cost)) ||
      nrow(cost) < 1L || ncol(cost) < 1L || any(cost < 0)) {
    stop("cost must be a finite nonnegative numeric matrix.", call. = FALSE)
  }
  if (length(labels) != ncol(cost) || anyNA(labels)) {
    stop("labels must align with cost columns and contain no NA.", call. = FALSE)
  }
  if (!is.null(cell_ids) && anyNA(cell_ids)) {
    stop("cell_ids cannot contain NA.", call. = FALSE)
  }
  storage.mode(cost) <- "double"
  engine <- .backend()$ReferenceAudit(
    cost, as.list(as.character(labels)), epsilon = epsilon,
    a = if (is.null(a)) NULL else as.list(as.numeric(a)),
    b = if (is.null(b)) NULL else as.list(as.numeric(b)),
    cell_ids = if (is.null(cell_ids)) NULL else as.list(as.character(cell_ids)),
    base_tol = base_tol)
  info <- engine$describe()
  info$classes <- as.character(info$classes)
  info$cell_ids <- as.character(info$cell_ids)
  info$baseline_proportions <- stats::setNames(as.numeric(info$baseline_proportions), info$classes)
  structure(c(list(engine = engine), info), class = "scRefAudit_fit")
}

#' Solve specified reference composition scenarios
#' @export
audit_compositions <- function(fit, proportions, order = 0L,
                               corrector = c("chord", "aa2"), tol = 1e-7,
                               batch_size = 32L, return_probabilities = FALSE) {
  .check_fit(fit)
  corrector <- match.arg(corrector)
  if (is.numeric(proportions) && is.null(dim(proportions))) {
    nm <- names(proportions)
    proportions <- matrix(proportions, nrow = 1L)
    colnames(proportions) <- nm
  }
  if (!is.matrix(proportions) || !is.numeric(proportions)) {
    stop("proportions must be a numeric vector or matrix.", call. = FALSE)
  }
  nm <- colnames(proportions)
  if (!is.null(nm)) {
    if (anyDuplicated(nm) || length(nm) != length(fit$classes) ||
        !setequal(nm, fit$classes)) {
      stop("proportions column names must exactly match fit$classes.", call. = FALSE)
    }
    proportions <- proportions[, match(fit$classes, nm), drop = FALSE]
  }
  storage.mode(proportions) <- "double"
  result <- fit$engine$audit_compositions(
    proportions, order = .scalar_integer(order, "order", 0L),
    corrector = corrector, tol = tol,
    batch_size = .scalar_integer(batch_size, "batch_size"),
    return_probabilities = isTRUE(return_probabilities))
  .decorate_audit(result, fit)
}

#' Generate derivative-guided or random composition scenarios
#' @export
propose_scenarios <- function(fit, budget = 32L, radius = 0.1,
                              metric = c("log_weight", "tv"),
                              method = c("directed", "random"), seed = 1L,
                              anchors = NULL) {
  .check_fit(fit)
  metric <- match.arg(metric)
  method <- match.arg(method)
  if (!is.null(anchors)) {
    if (!is.numeric(anchors) || !length(anchors) || any(!is.finite(anchors)) ||
        any(anchors != floor(anchors)) || any(anchors < 1L) ||
        any(anchors > fit$n_query) || anyDuplicated(anchors)) {
      stop("anchors must be unique one-based query row indices.", call. = FALSE)
    }
    anchors <- as.list(as.integer(anchors) - 1L)
  }
  result <- fit$engine$propose_scenarios(
    budget = .scalar_integer(budget, "budget"), radius = radius,
    metric = metric, method = method, seed = .scalar_integer(seed, "seed", 0L),
    anchors = anchors)
  .decorate_proposals(result, fit)
}

#' Rank cells for review using the local linearized radius
#' @export
review_cells <- function(fit, top_n = NULL) {
  .check_fit(fit)
  if (!is.null(top_n)) top_n <- .scalar_integer(top_n, "top_n")
  result <- .columns_frame(fit$engine$review_cells(top_n = top_n))
  result$cell_index <- as.integer(result$cell_index) + 1L
  result
}

#' Compare directed and random checks at a common scenario budget
#' @export
stress_test <- function(fit, budget = 32L, radius = 0.01,
                         metric = c("tv", "log_weight"), seed = 1L,
                         order = 1L, tol = 1e-7, batch_size = 32L,
                         return_probabilities = FALSE) {
  .check_fit(fit)
  metric <- match.arg(metric)
  result <- fit$engine$stress_test(
    budget = .scalar_integer(budget, "budget"), radius = radius,
    metric = metric, seed = .scalar_integer(seed, "seed", 0L),
    order = .scalar_integer(order, "order", 0L), tol = tol,
    batch_size = .scalar_integer(batch_size, "batch_size"),
    return_probabilities = isTRUE(return_probabilities))
  result$comparison <- .columns_frame(result$comparison)
  result$anchors <- .columns_frame(result$anchors)
  result$anchors$cell_index <- as.integer(result$anchors$cell_index) + 1L
  result$directed <- .decorate_proposals(result$directed, fit)
  result$random <- .decorate_proposals(result$random, fit)
  result$audit <- .decorate_audit(result$audit, fit)
  result
}

#' Small deterministic example
#' @export
demo_audit_data <- function() {
  x <- .backend()$demo_data()
  x$labels <- as.character(x$labels)
  x$cell_ids <- as.character(x$cell_ids)
  rownames(x$cost) <- x$cell_ids
  x
}

#' Inspect the Python environment used by this package
#' @export
backend_info <- function() {
  module <- .backend()
  list(scRefAudit = as.character(utils::packageVersion("scRefAudit")),
       python_package = module$`__version__`, python = reticulate::py_config())
}

#' @export
print.scRefAudit_fit <- function(x, ...) {
  cat("scRefAudit fit\n")
  cat("  Query cells:", x$n_query, " Reference cells:", x$n_reference, "\n")
  cat("  Classes:", paste(x$classes, collapse = ", "), "\n")
  cat("  Baseline marginal residual:", format(x$baseline_residual), "\n")
  invisible(x)
}
