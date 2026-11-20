# Validation status — 2026-11-21

This is a development software release. Tests check the new package interfaces
and their agreement with the existing numerical implementation; they do not
constitute a new manuscript benchmark or a new biological validation cohort.

## Completed

- All five files in the original `sopc_ot` numerical package are byte-for-byte
  unchanged. SHA-256 records are in `kernel_provenance.json`.
- Python wheel built and tested after installation into an isolated environment
  using the already installed NumPy/SciPy dependencies.
- Nine Python acceptance tests pass: independent log-domain Sinkhorn reference;
  three initialization orders and two correctors; finite-difference derivative;
  chunked probability and summary equivalence; matched-TV generation and seeds;
  ranking and intended-anchor endpoints; invalid input; tied baseline scores
  and single-query inputs; forced fallback.
- The runnable synthetic example completes and its expected output is included.
- Three archived real inputs were replayed: HCC myeloid (500 queries), lung
  epithelium (236 queries), and lung mesenchyme (40 queries), each with 384
  references. Baseline class scores matched the archived scores exactly in
  this environment. Directed vectors agreed to floating-point precision.
  All generated 1% TV scenarios met the 1e-7 marginal criterion. The three
  checks retained 0, 2 and 14 fallback scenarios respectively.

## R interface status

R/Rscript is not installed in the packaging environment. Therefore **R package
installation, reticulate conversion, R examples and R CMD check were not run**.
The R archive is assembled from the package source; it is not represented as
an R CMD build/check-certified release. `tests/smoke.R` and the GitHub Actions
workflow provide the next execution checks, including one-based indices,
named-class alignment, dimensions and result serialization on Linux/macOS.

The package has not been submitted to CRAN, Bioconductor or PyPI. No remote
GitHub repository or DOI was created during this packaging task.

## Re-run

```sh
python -m pip install .
python validation/validate_release.py
python validation/replay_real_inputs.py /path/to/matched_tv_20261118/prepared
SCREFAUDIT_RUN_R_TESTS=true Rscript tests/smoke.R
```

The original prepared real-data inputs are in the separate study archive and
are not duplicated in this software repository.
