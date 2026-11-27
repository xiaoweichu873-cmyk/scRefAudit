"""Recompute the four top-100 counts in Table 1; Figure 2B uses three of them."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import softmax

ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / "review"
FROZEN = ROOT / "frozen_cells/range"
PROTOCOL = json.loads((REVIEW / "range_extension/protocol.json").read_text())
TASKS = sorted(PROTOCOL["tasks"])
SIZES = [size for size, taxonomy in PROTOCOL["configurations"] if taxonomy == "native"]
SEEDS = PROTOCOL["seeds"]


def original_validation_module():
    spec = importlib.util.spec_from_file_location("author_validation", REVIEW / "scripts/02_validation.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def geometry(task):
    """Use the same random-number stream and full-matrix scale as run_range.py."""
    with np.load(ROOT / "timing/data/prepared" / task / "inputs.npz") as source:
        yr = source["yr"]
        if np.issubdtype(yr.dtype, np.integer):
            yr = source["classes"][yr]
        rng = np.random.default_rng(20261106)
        query = np.sort(rng.choice(len(source["yq"]), 1000, replace=False))
        full = cdist(source["zq"], source["zr"], "sqeuclidean")
        cost = full[query] / np.median(full)
    core = []
    for label in np.unique(yr):
        ids = np.flatnonzero(yr == label)
        core.extend(rng.choice(ids, min(2, len(ids)), replace=False))
    pool = np.setdiff1d(np.arange(len(yr)), core)
    small = np.sort(np.r_[core, rng.choice(pool, 300 - len(core), replace=False)])
    return cost, yr, query, small


def summarize(prob, scores, target):
    gap = np.diff(np.sort(prob, axis=1)[:, -2:], axis=1).ravel()
    columns = {"confidence": scores[:, 0], "entropy": scores[:, 1],
               "margin": -gap, "linear_radius": scores[:, 6]}
    return {name: int(target[np.argsort(-values, kind="stable")[:100]].sum())
            for name, values in columns.items()}


def from_frozen(folder):
    with np.load(folder / "baseline.npz") as b:
        prob, scores = b["probabilities"], b["scores"]
    target = np.zeros(len(prob), dtype=bool)
    residual = 0.0
    for seed in SEEDS:
        with np.load(folder / f"{seed}_r0.1_random.npz") as saved:
            target |= saved["switched"].any(axis=1)
            residual = max(residual, float(saved["marginal_residual"].max()))
    return summarize(prob, scores, target), int(target.sum()), residual


def from_prepared(module, cost, yr, query, refs, folder, check):
    selected_cost = cost[:, refs]
    model = module.CompositionAudit(selected_cost, yr[refs], epsilon=0.05)
    model.run(model.p0[None, :], order=1, corrector="chord", tol=1e-7)
    p0 = model.base_probabilities
    initial_label = p0.argmax(axis=1)
    jacobian = module.log_jacobian(model)
    scores, _, _, _ = module.scores(p0, jacobian, selected_cost, np.repeat(p0[:, None, :], 2, axis=1))
    if check:
        with np.load(folder / "baseline.npz") as saved:
            np.testing.assert_array_equal(query, saved["q_index"])
            np.testing.assert_array_equal(refs, saved["r_index"])
            np.testing.assert_array_equal(initial_label, saved["probabilities"].argmax(axis=1))
    target = np.zeros(len(p0), dtype=bool)
    residual = 0.0
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        rng.choice(len(p0), 32, replace=False)  # Anchors were sampled before random directions.
        directions = rng.normal(size=(32, len(model.classes)))
        directions -= directions.mean(axis=1, keepdims=True)
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        composition = softmax(np.log(model.p0) + 0.1 * directions, axis=1)
        fit = model.run(composition, order=1, corrector="chord", tol=1e-7)
        labels = fit["probabilities"].argmax(axis=2)
        changed = labels != initial_label[:, None]
        target |= changed.any(axis=1)
        residual = max(residual, float(fit["marginal_residual"].max()))
        if check:
            with np.load(folder / f"{seed}_r0.1_random.npz") as saved:
                np.testing.assert_allclose(composition, saved["proportions"], atol=1e-13, rtol=1e-13)
                np.testing.assert_array_equal(labels, saved["scenario_labels"])
                np.testing.assert_array_equal(changed, saved["switched"])
    return summarize(p0, scores, target), int(target.sum()), residual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["frozen", "solve"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--compare-archive", action="store_true", help="In solve mode, compare after each fit.")
    args = parser.parse_args()
    module = original_validation_module() if args.mode == "solve" else None
    args.out.parent.mkdir(exist_ok=True, parents=True)
    cols = ["task", "references", "confidence", "entropy", "margin", "linear_radius", "sampled_sensitive", "max_residual"]
    with args.out.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=cols)
        writer.writeheader()
        for task in TASKS:
            if module:
                cost, yr, query, small = geometry(task)
            for size in SIZES:
                folder = FROZEN / f"{task}_m{size}_native"
                if module:
                    refs = small if size == 300 else np.arange(len(yr))
                    counts, total, residual = from_prepared(module, cost, yr, query, refs, folder, args.compare_archive)
                else:
                    counts, total, residual = from_frozen(folder)
                if not np.isfinite(residual) or residual > 1e-7:
                    raise ValueError(f"Unmet marginal tolerance for {task}, {size}: {residual}")
                writer.writerow({"task": task, "references": size, "sampled_sensitive": total, "max_residual": residual, **counts})
                out.flush()
                print(task, size, counts, f"max residual {residual:.3g}", flush=True)

if __name__ == "__main__":
    main()
