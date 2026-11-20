"""Replay representative archived inputs; original data are supplied separately.

Usage: python validation/replay_real_inputs.py /path/to/matched_tv_20261118/prepared
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from screfaudit import ReferenceAudit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    args = parser.parse_args()
    rows = []
    for group in ["HCC_Myeloid", "IPF_Epithelial", "IPF_Mesenchymal"]:
        p = sorted((args.prepared / group).rglob("*.npz"))[0]
        with np.load(p, allow_pickle=False) as z:
            fit = ReferenceAudit(z["cost"], z["reference_labels"], epsilon=.05,
                                  b=z["b"], cell_ids=z["barcodes"])
            diff = float(np.max(np.abs(fit.baseline_probabilities - z["archived_baseline_probabilities"])))
            proposal = fit.propose_scenarios(anchors=z["anchors"], metric="tv", radius=.01)
            cosine = np.einsum("ij,ij->i", proposal["directions"], z["directed_directions"])
            out = fit.audit_compositions(proposal["proportions"], order=1, batch_size=16)
            row = dict(group=group, source=str(p.relative_to(args.prepared)),
                       source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                       n_query=fit.n_query, n_reference=fit.n_reference, n_classes=len(fit.classes),
                       max_baseline_difference=diff,
                       baseline_label_disagreements=int(np.sum(fit.baseline_probabilities.argmax(1) !=
                                                              z["archived_baseline_probabilities"].argmax(1))),
                       minimum_archived_direction_cosine=float(cosine.min()),
                       max_tv_error=float(np.max(np.abs(proposal["actual_tv"] - .01))),
                       max_residual=float(out["marginal_residual"].max()),
                       fallback_count=out["diagnostics"]["fallback_count"])
        assert diff < 1e-6 and row["baseline_label_disagreements"] == 0
        assert cosine.min() > 1 - 1e-5
        assert row["max_tv_error"] < 1e-12 and row["max_residual"] <= 1.01e-7
        rows.append(row)
    result = dict(purpose="Three archived real-input compatibility checks; not a new efficacy or timing study",
                  checks=rows)
    path = Path(__file__).parent / "real_input_replay.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
