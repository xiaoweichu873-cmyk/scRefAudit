"""Run from the installed Python package: python examples/quickstart.py."""
import json
from screfaudit import ReferenceAudit, demo_data

data = demo_data()
fit = ReferenceAudit(**data)
result = fit.stress_test(budget=16, radius=0.01, metric="tv", seed=37)
review = fit.review_cells(top_n=10)
summary = dict(
    software="scRefAudit 0.1.0", example="synthetic, not a manuscript benchmark",
    classes=fit.classes.tolist(),
    comparison={k: v.tolist() for k, v in result["comparison"].items()},
    first_ten_review_ids=review["cell_id"].tolist(),
    max_marginal_residual=float(result["audit"]["marginal_residual"].max()),
    fallback_count=int(result["audit"]["diagnostics"]["fallback_count"]))
print(json.dumps(summary, indent=2))
