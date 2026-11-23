# scRefAudit

**Single-cell reference-composition audit**

scRefAudit checks how balanced optimal-transport cell annotations change when
the proportions of reference cell types change. It provides a Python package
and an R interface for the same numerical implementation.

- **Repeated computation:** reuse a baseline curvature factorization across
  composition scenarios, with optional Taylor initialization and numerical
  correction.
- **Targeted checks:** use composition derivatives to propose a change likely
  to alter a particular cell's assignment.
- **Review priorities:** rank cells by a local linearized boundary radius.

The numerical kernel is copied unchanged from the manuscript code. Its hashes
are recorded in `validation/kernel_provenance.json`. This software release does
not change the paper's results or rerun its performance benchmark.

## Exploratory unbalanced OT extension — 23 November 2026

A separate [KL-unbalanced OT addendum](manuscript_data/unbalanced_ot_20261123/README.md)
provides the derivation, all 96 benchmark configurations, numerical checks, and
exact code and prepared inputs. The revised method was faster in 57/96 settings,
with a median paired audit ratio of **1.34×** (**1.26×** including common
initialization). The median was **2.26×** for local perturbations and **0.98×**
for larger perturbations; all speedups and slowdowns are retained.

The initial five-control comparison and the subsequent fallback retest used the
same development inputs. This is an experimental solver extension with a
first-order predictor, separate from the current balanced R/Python interface.
See the addendum for timing definitions, scope and reproduction commands.

## Installation

### R / RStudio

Install `reticulate` first, then the local source archive supplied with this release:

```r
install.packages("reticulate")
install.packages("scRefAudit_0.1.0.tar.gz", repos = NULL, type = "source")
library(scRefAudit)
```

If working from the unpacked repository directory, run `R CMD INSTALL .` in a
terminal instead. The R interface requires R >= 4.1 and reticulate >= 1.41.
Reticulate declares NumPy and SciPy dependencies and can provision Python at
the first computation. This step needs internet access if dependencies are
not already available. To use an existing compatible environment, set
`RETICULATE_PYTHON` to that interpreter **before** initializing Python in R.

This is a development release, not a CRAN or Bioconductor release. See
`validation/VALIDATION.md` for the checks actually performed and R test status.

### Python

From the repository directory:

```sh
python -m pip install .
python examples/quickstart.py
```

Alternatively install the supplied wheel. Python >= 3.10, NumPy >= 1.26 and
SciPy >= 1.11 are required. A PyPI release has not been published.

## Quick start in R

```r
library(scRefAudit)
x <- demo_audit_data()           # synthetic example, not study data
fit <- audit_fit(x$cost, x$labels, epsilon = x$epsilon)

# Compare 16 directed and 16 random checks, each moving 1% of class mass.
checks <- stress_test(fit, budget = 16L, radius = 0.01, metric = "tv", seed = 37L)
checks$comparison
review_cells(fit, top_n = 10L)

# Supply your own reference compositions, one row per scenario.
p0 <- fit$baseline_proportions
p1 <- 0.9 * p0 + 0.1 * rep(1 / length(p0), length(p0))
result <- audit_compositions(fit, rbind(p0, p1), order = 2L)
head(result$cells)
```

Named proportion columns are aligned to `fit$classes`; unnamed columns must
already follow that order. Values must be positive and rows must sum to one.
R indices, including anchors and fallback indices, are **one-based**.

## Quick start in Python

```python
from screfaudit import ReferenceAudit, demo_data

fit = ReferenceAudit(**demo_data())
checks = fit.stress_test(budget=16, radius=0.01, metric="tv", seed=37)
print(checks["comparison"])
print(fit.review_cells(top_n=10))
```

Python indices are **zero-based**. Probability arrays, when requested with
`return_probabilities=True`, have axes `(query cell, scenario, class)`.
`CompositionAudit` remains available for compatibility with the manuscript API.

## Using real single-cell data

Supply a **fixed query-by-reference cost matrix** and reference annotations.
Choose normalization, features, representation and distance scaling before
fitting. Keep query cells, reference cells, costs and epsilon fixed within an
audit. Reference weights are varied within the existing classes; zero weights,
missing reference classes and unbalanced OT are outside this implementation.

The package does not load raw counts or refit an atlas automatically. It can
be inserted after a reference-based representation and distance matrix have
been calculated. In R, `result$cells` is an ordinary data frame with `cell_id`;
it can be joined to a Seurat metadata table by that identifier. See
`docs/R_WORKFLOW.md` for a minimal example.

## Reading the results

| Output | Meaning |
|---|---|
| `switch_frequency` | Fraction of evaluated scenarios changing a cell's baseline label |
| `max_probability_change` | Largest class-score change from baseline |
| `probability_span` | Largest class-score range across evaluated scenarios |
| `linearized_radius` | Local approximation to a score-boundary distance in centered log-weight L2 space |
| `detected_fraction` | Fraction of all queries changed by at least one scenario |
| `anchor_hit_fraction` | Fraction of selected anchors changed by their own paired scenario |
| `actual_tv` | Half the sum of absolute changes in reference class proportions |

All finite proposed scenarios are corrected before labels are reported.
Fallback use and final marginal residuals remain in the returned diagnostics.
Different `log_weight` directions can produce different actual TV changes.
Use `metric="tv"` for comparisons at matched actual composition change.
An unreachable TV request raises an error rather than being silently reduced.

Sensitivity is a dependence on reference weights. It is not a probability of
misannotation, an absolute cell count, or experimental validation of cell identity.
The linearized radius is not an exact certified robustness radius.

## Computational scope

This is a dense CPU implementation. Kernel/curvature storage is
`O(n_query*n_reference + n_reference^2)`. Scenario batching limits temporary
arrays; returning full probabilities additionally stores
`n_query*n_scenarios*n_classes` values. Derivative-based ranking has further
class-Jacobian storage. Use a suitable reference subset for large atlases.
Small easy workloads can favor ordinary accelerated Sinkhorn. The new R and
Python wrapper has not been assigned the manuscript's 2.55-fold speedup.

## Tests and provenance

```sh
python -m unittest discover -s tests/python -v
SCREFAUDIT_RUN_R_TESTS=true Rscript tests/smoke.R
```

Tests cover an independently implemented log-domain Sinkhorn reference,
finite-difference derivatives, matched TV, batching, intended-anchor endpoints,
fallback, and invalid input. `validation/` records checks performed for this
release. The accompanying GitHub workflow tests Python and the R interface on
Linux and macOS when it is run; its results are not preclaimed here.

## Manuscript and citation

Preprint: *Efficient reference composition sensitivity analysis for single-cell
annotation*, DOI [10.21203/rs.3.rs-11637549/v1](https://doi.org/10.21203/rs.3.rs-11637549/v1).
The newer submission manuscript contains additional controls and analyses.

The three supplementary result workbooks are available in
[`manuscript_data/`](manuscript_data/README.md), with download links and SHA-256
checksums. They contain manuscript tables, figure source values, timing records
and numerical diagnostics. The workbooks are unchanged from the submission
files; no analyses were rerun for this data deposit.

This repository is not the complete experimental archive: full processed
benchmarking inputs, all scenario predictions and comparator build materials
are not supplied by these result workbooks. No permanent data/software DOI has
yet been assigned to this deposit.

## License and development

MIT; see `LICENSE.md`. This uses the license proposed in the earlier software
bundle; authors should confirm the release metadata before public deposition.
NumPy, SciPy and reticulate are external dependencies under their own licenses.
SSNS/SPLR source code is not bundled here. Software packaging and interface
development were assisted by OpenAI ChatGPT/Codex.

Chinese instructions: [README_zh.md](README_zh.md).
