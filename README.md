# scRefAudit

**Find which single-cell labels to review when reference proportions are uncertain.**

In reference-based single-cell annotation, the same cell can get a different
label if one cell type is given more weight in the reference, even when the
cells and their expression profiles stay fixed. scRefAudit checks this
dependence for balanced optimal-transport (OT) annotation before the labels
are used in downstream comparisons.

- **Prioritize review.** `review_cells()` ranks cells by the estimated change
  in reference proportions needed to switch their label.
- **Check a cell of interest.** `stress_test()` proposes targeted changes in
  reference proportions and checks whether the selected cell's label switches.
- **Compare specific reference compositions.** `audit_compositions()` reports
  which cells switch and how their class scores change across scenarios.

In eight configurations from two immune-cell datasets, the top 100 cells
ranked by scRefAudit contained an average of **70.63 composition-sensitive
cells**, versus **58.25** when ranked by the gap between the two highest label
scores, under the same independent set of sampled reference compositions. A
label switch is a reason to review the assignment, not proof it is wrong.

Start with a fixed query-by-reference cost matrix and reference cell labels.
The R and Python interfaces use the same numerical implementation.
[Quick start in R](#quick-start-in-r) · [Quick start in Python](#quick-start-in-python)
· [中文说明](README_zh.md)

## Installation

### R / RStudio

[Download the source ZIP](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/refs/heads/main.zip) and unzip it. In RStudio, replace the example path with the path to your unzipped `scRefAudit-main` folder:

```r
install.packages("reticulate")
install.packages("/your/path/scRefAudit-main", repos = NULL, type = "source")
library(scRefAudit)
```

Alternatively, run `R CMD INSTALL .` in a terminal opened in the unzipped
folder. The R interface requires R >= 4.1 and reticulate >= 1.41.
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

Python >= 3.10, NumPy >= 1.26 and
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

## Acceleration and computational scope

For repeated checks, scRefAudit reuses a baseline curvature factorization
and can use Taylor predictions to initialize new solves. In the manuscript's
single-threaded CPU benchmark, the research implementation was faster than
Anderson-accelerated Sinkhorn (AA2) in 10 of 11 completed paired configurations,
with a median speedup of 2.55×. The numerical kernel is unchanged from the
manuscript code; hashes are in `validation/kernel_provenance.json`.

This is a dense CPU implementation. Kernel/curvature storage is
`O(n_query*n_reference + n_reference^2)`. Scenario batching limits temporary
arrays; returning full probabilities additionally stores
`n_query*n_scenarios*n_classes` values. Derivative-based ranking has further
class-Jacobian storage. Use a suitable reference subset for large atlases.
Small easy workloads can favor ordinary accelerated Sinkhorn. The new R and
Python wrapper has not been assigned the manuscript's 2.55-fold speedup.

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

The [fixed supporting-data archive](https://github.com/xiaoweichu873-cmyk/scRefAudit/tree/859f41bc4c9752bb02580d557bf2afe2d79b8e58) includes the three result workbooks, 78 prepared tissue inputs for Figure 7, source records, plotting scripts and the exploratory KL supplement. [Download the full archive](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/859f41bc4c9752bb02580d557bf2afe2d79b8e58.zip).

Preprint: *Efficient reference composition sensitivity analysis for single-cell
annotation*, DOI [10.21203/rs.3.rs-11637549/v1](https://doi.org/10.21203/rs.3.rs-11637549/v1).
The newer submission manuscript contains additional controls and analyses.

The three supplementary result workbooks are available in
[`manuscript_data/`](manuscript_data/README.md), with download links and SHA-256
checksums. They contain manuscript tables, figure source values, timing records
and numerical diagnostics. The workbooks are unchanged from the submission
files; no analyses were rerun for this data deposit.

The [balanced-OT reproducibility download](manuscript_data/balanced_ot_v5/README.md) supplies the exact prepared inputs, research drivers, comparator build sources and expected outputs for Figure 2B / Table 1 and Figure 5. It starts at the prepared embeddings. The complete experimental record for all other analyses is outside the scope of that download.

## License and development

MIT; see `LICENSE.md`. This uses the license proposed in the earlier software
bundle; authors should confirm the release metadata before public deposition.
NumPy, SciPy and reticulate are external dependencies under their own licenses.
The balanced benchmark archive retains the SSNS/SPLR adapter source and its upstream license. Software packaging and interface
development were assisted by OpenAI ChatGPT/Codex.

Chinese instructions: [README_zh.md](README_zh.md).
