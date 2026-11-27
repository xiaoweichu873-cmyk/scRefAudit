"""Run fresh isolated balanced-OT timing trials with the authors' workers."""
import argparse
import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--grid", action="store_true", help="Run all 12 configurations and all five arms (180 scheduled trials).")
args = parser.parse_args()
out = args.out.resolve()
if out.exists() and any(out.iterdir()):
    parser.error("Use a new, empty output directory to avoid loading saved results.")
out.mkdir(parents=True, exist_ok=True)
os.environ["SOPC_FOLLOWUP_OUT"] = str(out)
for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
sys.path.insert(0, str(ROOT / "timing"))
import run_followup as bench
from threadpoolctl import threadpool_limits

configs = bench.PROTOCOL["main_configs"] + bench.PROTOCOL["scale_configs"]
arms = bench.PROTOCOL["methods"]
if not args.grid:
    selected = {("kang_f0", 1200, 0.1, 16), ("hao_f0", 1200, 0.05, 64)}
    arms = [m for m in arms if m in ("reuse_AA2", "reuse_chord", "quadratic_chord")]
columns = ["task", "m", "eps", "B", "method", "rep", "status", "total_s", "max_residual", "fallback_count"]
with threadpool_limits(limits=1), (out / "timings.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for i, config in enumerate(configs):
        if not args.grid and tuple(config[k] for k in ("task", "m", "eps", "B")) not in selected:
            continue
        bench.prepare(config)
        for rep in range(bench.PROTOCOL["repeats"]):
            shift = (i + rep) % len(arms)
            for method in arms[shift:] + arms[:shift]:
                row = bench.monitor(config, method, rep)
                writer.writerow(row)
                handle.flush()
                if not args.grid and row["status"] != "ok":
                    raise RuntimeError(f"{bench.key(config)} {method} {rep}: {row}")
