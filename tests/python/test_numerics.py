"""Numerical acceptance tests using an independent log-domain reference."""
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
from scipy.special import logsumexp, softmax
from screfaudit import ReferenceAudit, demo_data


def independent_sinkhorn(cost, a, b, epsilon, tol=1e-13):
    """Plain log-domain alternating scaling, independent of the bundled core."""
    logK = -cost / epsilon
    logv = np.zeros(len(b))
    for iteration in range(100000):
        logu = np.log(a) - logsumexp(logK + logv[None, :], axis=1)
        logv = np.log(b) - logsumexp(logK + logu[:, None], axis=0)
        if iteration % 5 == 0:
            G = np.exp(logu[:, None] + logK + logv[None, :])
            residual = max(np.abs(G.sum(axis=1) - a).sum(),
                           np.abs(G.sum(axis=0) - b).sum())
            if residual <= tol:
                return G
    raise RuntimeError("Independent log-domain solver did not converge")


class NumericalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.example = demo_data()
        cls.fit = ReferenceAudit(**cls.example)
        cls.rng = np.random.default_rng(194)
        cls.ps = 0.9 * cls.fit.baseline_proportions + 0.1 * cls.rng.dirichlet(np.ones(4), 9)

    def test_01_preserved_numerical_source(self):
        root = Path(__file__).resolve().parents[2]
        recorded = json.loads((root / "validation/kernel_provenance.json").read_text())
        for item in recorded["files"]:
            self.assertEqual(hashlib.sha256((root / item["file"]).read_bytes()).hexdigest(),
                             item["sha256"])

    def test_02_independent_solver_all_orders(self):
        fit = self.fit
        a = fit._fit.a
        refs = []
        for p in self.ps[:3]:
            G = independent_sinkhorn(self.example["cost"], a, fit._fit.W @ p, fit.epsilon)
            refs.append(G @ fit._fit.Y / a[:, None])
        refs = np.stack(refs, axis=1)
        for order in (0, 1, 2):
            for corrector in ("chord", "aa2"):
                out = fit.audit_compositions(self.ps[:3], order=order, corrector=corrector,
                                             tol=1e-10, return_probabilities=True)
                np.testing.assert_allclose(out["probabilities"], refs, rtol=0, atol=2e-7)
                np.testing.assert_array_equal(out["probabilities"].argmax(2), refs.argmax(2))
                self.assertLessEqual(out["marginal_residual"].max(), 1.01e-10)

    def test_03_probability_derivative(self):
        fit = self.fit
        J = fit.composition_jacobian()
        h = 1e-4
        d = np.array([0.3, -0.1, 0.2, -0.4]); d /= np.linalg.norm(d)
        P = []
        for sign in (-1, 1):
            p = softmax(np.log(fit.baseline_proportions) + sign * h * d)
            G = independent_sinkhorn(self.example["cost"], fit._fit.a,
                                     fit._fit.W @ p, fit.epsilon)
            P.append(G @ fit._fit.Y / fit._fit.a[:, None])
        predicted = np.einsum("ijk,j->ik", J, d)
        np.testing.assert_allclose((P[1] - P[0]) / (2*h), predicted, rtol=1e-5, atol=1e-7)
        np.testing.assert_allclose(J.sum(axis=1), 0, atol=1e-12)
        np.testing.assert_allclose(J.sum(axis=2), 0, atol=1e-9)

    def test_04_chunking_and_summary(self):
        x = self.fit.audit_compositions(self.ps, order=2, batch_size=2,
                                       tol=1e-10, return_probabilities=True)
        y = self.fit.audit_compositions(self.ps, order=2, batch_size=99,
                                       tol=1e-10, return_probabilities=True)
        np.testing.assert_allclose(x["probabilities"], y["probabilities"], atol=1e-8, rtol=0)
        base = self.fit.baseline_probabilities
        hits = (x["probabilities"].argmax(2) != base.argmax(1)[:, None]).sum(1)
        np.testing.assert_array_equal(x["cells"]["switch_count"], hits)
        expected_span = (x["probabilities"].max(1) - x["probabilities"].min(1)).max(1)
        np.testing.assert_allclose(x["cells"]["probability_span"], expected_span)
        z = self.fit.audit_compositions(self.ps, batch_size=2)
        self.assertIsNone(z["probabilities"])

    def test_05_tv_matching_and_reproducibility(self):
        for metric, radius in [("tv", .01), ("tv", .02), ("log_weight", .1)]:
            proposals = []
            for method in ("directed", "random"):
                p = self.fit.propose_scenarios(budget=11, metric=metric, radius=radius,
                                               method=method, seed=18)
                again = self.fit.propose_scenarios(budget=11, metric=metric, radius=radius,
                                                   method=method, seed=18)
                np.testing.assert_array_equal(p["proportions"], again["proportions"])
                np.testing.assert_allclose(p["directions"].sum(1), 0, atol=1e-14)
                np.testing.assert_allclose(np.linalg.norm(p["directions"], axis=1), 1, atol=1e-14)
                if metric == "tv":
                    np.testing.assert_allclose(p["actual_tv"], radius, atol=1e-13, rtol=0)
                proposals.append(p)
            np.testing.assert_array_equal(proposals[0]["anchor_index"], proposals[1]["anchor_index"])

    def test_06_rank_direction_and_endpoints(self):
        r = self.fit.review_cells()
        self.assertTrue(np.all(np.diff(r["linearized_radius"]) >= 0))
        winner, competitor, radius, direction, valid = self.fit._boundaries()
        J = self.fit.composition_jacobian();P = self.fit.baseline_probabilities
        row = np.arange(self.fit.n_query)
        g = J[row, :, winner] - J[row, :, competitor]
        slopes = np.einsum("ij,ij->i", g, direction)
        self.assertTrue(np.all(slopes[valid] < 0))
        np.testing.assert_allclose(P[row, winner] - P[row, competitor] + radius * slopes, 0, atol=1e-12)
        result = self.fit.stress_test(budget=7, seed=29, return_probabilities=True, tol=1e-10)
        labels = result["audit"]["probabilities"].argmax(2)
        changed = labels != winner[:, None]
        np.testing.assert_array_equal(result["comparison"]["detected_cells"],
                                      [changed[:, :7].any(1).sum(), changed[:, 7:].any(1).sum()])
        anchors = result["anchors"]["cell_index"]
        np.testing.assert_array_equal(result["anchors"]["directed_hit"], changed[anchors, np.arange(7)])

    def test_07_invalid_inputs_and_unreachable_tv(self):
        for ps in [np.zeros((2, 4)), np.full((2, 4), np.nan), np.ones((2, 4))]:
            with self.assertRaises(ValueError): self.fit.audit_compositions(ps)
        for kwargs in [dict(budget=1.2), dict(anchors=[0, 0]), dict(anchors=[-1]),
                       dict(metric="tv", radius=1), dict(radius=0), dict(seed=-1)]:
            with self.assertRaises(ValueError): self.fit.propose_scenarios(**kwargs)
        with self.assertRaises(ValueError):
            self.fit.propose_scenarios(metric="tv", radius=.99999)
        with self.assertRaises(ValueError):
            ReferenceAudit(np.ones((3, 4)), ["A", "A", "A", "A"])
        with self.assertRaises(ValueError):
            ReferenceAudit(np.ones((3, 4)), ["A", None, "B", "B"])
        with self.assertRaises(ValueError):
            ReferenceAudit(np.ones((3, 4)), ["A", "A", "B", "B"], cell_ids=["x", "x", "z"])

    def test_08_tied_scores_and_single_query(self):
        flat = ReferenceAudit(np.ones((3, 4)), ["A", "A", "B", "B"], epsilon=.5)
        # A constant expression distance still responds to marginal composition:
        # tied baseline scores give a zero boundary radius, not a zero Jacobian.
        np.testing.assert_allclose(flat.review_cells()["linearized_radius"], 0, atol=1e-12)
        directed = flat.propose_scenarios(budget=1, metric="tv", radius=.01)
        self.assertEqual(directed["proportions"].shape, (1, 2))
        random = flat.propose_scenarios(method="random", budget=1, metric="tv", radius=.01)
        self.assertEqual(random["proportions"].shape, (1, 2))
        one = ReferenceAudit(np.array([[0.0, .3, .4, .1]]), ["A", "A", "B", "B"], epsilon=.5)
        out = one.audit_compositions([.51, .49], return_probabilities=True)
        self.assertEqual(out["probabilities"].shape, (1, 1, 2))

    def test_09_forced_fallback(self):
        from sopc_ot.chord import chord_solve, factor_target
        fit = self.fit._fit
        U,V,it,res,info = chord_solve(fit.K, fit.a, fit.W @ self.ps.T,
                                     np.tile(fit.v[:,None], (1,len(self.ps))),
                                     factor_target(fit.G,fit.a,fit.b), maxiter=0, tol=1e-10)
        self.assertEqual(info["fallback_count"], len(self.ps))
        self.assertLessEqual(res.max(), 1.01e-10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
