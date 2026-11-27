# scRefAudit balanced optimal transport experiments

This download accompanies the Technical Note *Efficient reference composition sensitivity analysis for single-cell annotation*. It contains the prepared input matrices, frozen cell-level results, the research implementation and the two benchmark protocols for **Figure 2B and Table 1** (review yield) and **Figure 5A–B** (timing). It uses the original balanced-transport research drivers. The public scRefAudit R and Python interfaces are released separately at [the project repository](https://github.com/xiaoweichu873-cmyk/scRefAudit). 中文说明见 [README_zh.md](README_zh.md).

## Start here

Use CPython 3.12 in a Linux x86-64 environment for all three checks. The original Figure 5 timing worker uses Linux process affinity and the bundled CPython 3.12 AVX2/OpenMP extension. Its source, patch, upstream license and build metadata are in `timing/source/regot_followup/` and `timing/results/followup/exact_build_environment.json`. The cross-platform table calculation and the Figure 2B solver do not load the extension.

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-reviewer.txt
sha256sum -c SHA256SUMS
python scripts/public_tables.py --out output/tables
python scripts/review_cells.py --mode frozen --out output/frozen_review.csv
python scripts/review_cells.py --mode solve --compare-archive --out output/fresh_review.csv
python scripts/timing_trial.py --out output/timing_trial
```

## What each command tests

| Command | Inputs | Expected output |
|---|---|---|
| `public_tables.py` | Public supplementary workbooks in `workbooks/`. Every Figure 5 ratio is calculated from individual trial records. | `figure5_ratios.csv`, `table1.csv`, `figure5_from_workbooks.pdf` and `figure2B_from_workbooks.pdf`. Figure 5A shows 10 wins in 11 complete pairs and a median 2.55-fold speedup (0.60–8.10). Figure 5B shows 8 wins in 12, median 1.10-fold. |
| `review_cells.py --mode frozen` | Baseline probabilities and scores plus archived random-scenario switch masks under `frozen_cells/range/`. | Four top-100 counts for each of eight configurations. Compare with `expected/figure2B_from_frozen.csv`. This step reranks previously computed outputs. |
| `review_cells.py --mode solve` | Four donor-separated tasks in `timing/data/prepared/{hao_f0,hao_f1,kang_f0,kang_f1}/inputs.npz`. Starts from the prepared embeddings, fits the baseline, evaluates derivatives and solves new compositions. | The 32 top-100 counts match `expected/table1_counts.csv`. The rounded means for confidence, entropy, top-two margin and linearized radius are 50.25, 35.88, 58.25 and 70.63. With `--compare-archive`, all 768,000 scenario labels are checked after fitting. The maximum marginal residual in the archived replay was below 10⁻⁷. |
| `timing_trial.py` | Prepared embeddings and archived mixtures; two configurations, three methods, three independent worker runs. Additional reference sizes for the full grid are included under `timing/data/prepared/`. | 18 trial records in `timings.csv`, worker logs and probabilities. The archived independent check found median AA2/proposed ratios of 0.86 for Kang fold 1 (1,200 references, ε=0.10, 16 scenarios) and 4.78 for Hao P1 (1,200 references, ε=0.05, 64 scenarios). These illustrative timings depend on the hardware. |

For the full five-arm timing grid, run `python scripts/timing_trial.py --grid --out output/full_timing`. Use a fresh output directory. The historical grid contained 180 planned trials: 152 completed, 10 timed out and 18 were skipped after a previous failure. All 36 proposed-method trials completed. The median 2.55-fold value is the historical research-benchmark result. The fresh trial in this release covers two configurations. Time for rebuilding common geometry and the common baseline is recorded separately.

Each review configuration uses 1,000 query cells, 300 or 1,200 reference cells, and 96 random compositions across seeds 20261106, 20261107 and 20261108, at centered log-weight radius 0.10. A cell counts once if its top label changes in any scenario. A change identifies reference-composition sensitivity. The review benchmark does not adjudicate cell identity. The raw GEO counts to PCA preprocessing step was not rerun in this check.

## Versions and provenance

- The fixed workbook release is commit `abdbacfcfe5868e08c4e87e4f87691225b31d003`. The source workbook sheet `Table_4` corresponds to **Table 1** in the Technical Note. Older filenames, original code and workbook indexes use **Figure 7B / Table 4** for the review counts and **Figure 4** for the timing comparison.
- The research implementation in the two archives identifies itself as `sopc_ot 0.2.0.dev0`. The archived Figure 5 driver `timing/run_followup.py` has SHA-256 prefix `3b065dac975b`. The archived range driver `review/range_extension/scripts/run_range.py` has prefix `b39fb79002b8`. The original derivative/score and random-only ranking scripts have prefixes `1fdf88559a40` and `566f67f8f568`, respectively. The protocol and amendment, software environment, raw archived output hashes and exact figures tested are recorded in `expected/verification.json`; `SHA256SUMS` covers the files shipped here.
- The included replay scripts use v5 display numbers. The four 1,200-reference starting archives are shared between the two analyses and stored once. To run the broader original range study using its original path convention, run `python scripts/restore_archived_layout.py` first. Restoration checks each file against the author's original digest.
- Original expression datasets: Kang [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583) and Hao [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378). Source provenance records and task-specific metadata accompany the computed embeddings. The original sources retain their terms. scRefAudit software uses the project MIT license. Third-party RegOT code keeps its own MIT notice under `timing/source/regot_followup/LICENSE`.
