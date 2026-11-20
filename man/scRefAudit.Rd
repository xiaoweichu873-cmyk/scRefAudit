\name{scRefAudit}
\alias{scRefAudit}
\alias{audit_fit}
\alias{audit_compositions}
\alias{propose_scenarios}
\alias{review_cells}
\alias{stress_test}
\alias{demo_audit_data}
\alias{backend_info}
\alias{print.scRefAudit_fit}
\title{Reference Composition Sensitivity Analysis for Single-Cell Annotation}
\description{Fit a fixed-cost balanced entropic optimal transport model,
reuse curvature across composition scenarios, generate directed checks,
and prioritize cells for annotation review.}
\usage{
audit_fit(cost, labels, epsilon = 0.05, a = NULL, b = NULL,
          cell_ids = rownames(cost), base_tol = 1e-11)

audit_compositions(fit, proportions, order = 0L,
                   corrector = c("chord", "aa2"), tol = 1e-7,
                   batch_size = 32L, return_probabilities = FALSE)

propose_scenarios(fit, budget = 32L, radius = 0.1,
                  metric = c("log_weight", "tv"),
                  method = c("directed", "random"), seed = 1L,
                  anchors = NULL)

review_cells(fit, top_n = NULL)

stress_test(fit, budget = 32L, radius = 0.01,
            metric = c("tv", "log_weight"), seed = 1L,
            order = 1L, tol = 1e-7, batch_size = 32L,
            return_probabilities = FALSE)

demo_audit_data()
backend_info()
\method{print}{scRefAudit_fit}(x, ...)
}
\arguments{
\item{cost}{Finite nonnegative numeric matrix, query cells by reference cells.
Representation and distance scaling must be chosen before fitting.}
\item{labels}{Reference labels aligned with columns of \code{cost}. At least two classes.}
\item{epsilon}{Positive entropic regularization in cost units.}
\item{a, b}{Optional positive query and reference marginals summing to one. Default uniform.}
\item{cell_ids}{Unique query identifiers. Defaults to row names, or generated identifiers.}
\item{base_tol}{Baseline marginal residual tolerance.}
\item{fit}{An object returned by \code{audit_fit}.}
\item{proportions}{Positive class-proportion vector or scenario-by-class matrix;
each row sums to one. Named columns are aligned to \code{fit$classes}.}
\item{order}{Zero, one or two: reuse baseline, linear or quadratic initialization.}
\item{corrector}{Shared-curvature chord correction or batched accelerated Sinkhorn.}
\item{tol}{Final marginal residual tolerance.}
\item{batch_size}{Maximum scenarios solved per chunk.}
\item{return_probabilities}{Retain the full query-by-scenario-by-class array.}
\item{budget}{Number of scenarios per method, capped by query count.}
\item{radius}{Positive perturbation size. With \code{metric="tv"}, 0.01 means one percent.}
\item{metric}{Centered log-weight Euclidean radius or matched actual total variation.}
\item{method}{Derivative-guided or random proposal.}
\item{seed}{Nonnegative integer random seed.}
\item{anchors}{Optional unique one-based query indices. If supplied, they determine the scenario count.}
\item{top_n}{Maximum cells returned by the review ranking; default all cells.}
\item{x}{A fitted object.}
\item{...}{Additional print arguments, currently unused.}
}
\value{
\code{audit_fit} returns a \code{scRefAudit_fit} with class order, baseline
proportions, query identifiers, solver diagnostics and a live Python engine.
\code{audit_compositions} returns cell-level summaries, corrected marginal
residuals, iteration counts, fallback diagnostics and optionally probabilities.
\code{propose_scenarios} returns class proportions, directions, anchor identifiers,
and actual total variation. These are proposals until numerically corrected.
\code{review_cells} returns a data frame ordered by local linearized radius.
\code{stress_test} returns directed/random comparisons, paired anchor hits,
proposals and corrected audit results. Its \code{detected_fraction} measures
the union across queries; \code{anchor_hit_fraction} measures the assigned
scenario for each anchor. \code{demo_audit_data} returns a synthetic example.
\code{backend_info} returns software and Python environment information.
}
\details{
The R interface calls the bundled Python kernel using \pkg{reticulate}.
NumPy and SciPy dependencies are declared on package load; computation may
require an initial download when no compatible Python environment exists.

Cell indices in R are one-based. Sensitivity describes annotation dependence
on reference weights. The local radius is an approximation, not a certified
exact boundary or a biological error probability. An unreachable total
variation request raises an error. Positive marginals and a representable
positive kernel are required. Kernel underflow is rejected.

Reference cells, query cells, representation, costs and epsilon stay fixed
within an audit. Dense kernel/curvature storage is proportional to
\eqn{nm + m^2}, where \eqn{n} is query count and \eqn{m} reference count.
Changing the cells or representation requires a new fit.

Save ordinary result lists with \code{saveRDS}; recreate fit objects in a new
R session because their Python pointers cannot be restored from RDS.
}
\examples{
\dontrun{
library(scRefAudit)
x <- demo_audit_data()
fit <- audit_fit(x$cost, x$labels, epsilon = x$epsilon)
out <- stress_test(fit, budget = 16L, metric = "tv", radius = 0.01)
out$comparison
review_cells(fit, top_n = 10L)
}
}
