"""Documented user interface around the unchanged manuscript numerical kernel.

Arrays of probabilities use (query cell, scenario, class) order. Python cell
indices are zero-based; the R interface converts them to one-based indices.
"""
from numbers import Integral
import time

import numpy as np
from scipy.special import softmax
from scipy.spatial.distance import cdist
from sopc_ot import CompositionAudit


def _integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return int(value)


def _strings(values, name, n, unique=False):
    x = np.asarray(values, dtype=object)
    if x.shape != (n,):
        raise ValueError(f"{name} must have length {n}")
    for value in x:
        if value is None or str(value).strip() == "":
            raise ValueError(f"{name} cannot contain missing or empty values")
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            raise ValueError(f"{name} cannot contain missing values")
    x = x.astype(str)
    if unique and len(np.unique(x)) != n:
        raise ValueError(f"{name} must be unique")
    return x


def demo_data():
    """Return a small deterministic synthetic example, not study data."""
    rng = np.random.default_rng(37)
    query = rng.normal(size=(96, 4))
    reference = rng.normal(size=(80, 4))
    cost = cdist(query, reference, "sqeuclidean")
    cost /= np.median(cost)
    return dict(cost=cost, labels=np.array([f"type_{j % 4 + 1}" for j in range(80)]),
                cell_ids=np.array([f"cell_{j + 1}" for j in range(96)]), epsilon=0.2)


class ReferenceAudit:
    """Fit a fixed-cost balanced OT model and examine reference-class weights.

    Parameters
    ----------
    cost : array (n_query, n_reference)
        Finite nonnegative distances in an already chosen representation.
        No normalization, PCA or integration is performed by this constructor.
    labels : sequence of length n_reference
        Reference annotations; at least two distinct classes are needed.
    epsilon : positive float
        Entropic regularization in the same units as ``cost``.
    a, b : optional positive vectors summing to one
        Query and reference marginals. Uniform if omitted.
    cell_ids : optional unique strings
        Identifiers in cost-row order.
    base_tol : positive float
        Baseline marginal tolerance, default 1e-11.

    Notes
    -----
    Class order is lexicographic and is available in ``classes``. The kernel
    deliberately rejects numerical underflow; rescale distances or use a
    different stabilized solver for that regime. Dense kernel and curvature
    storage scale as O(n_query*n_reference + n_reference**2).
    """

    def __init__(self, cost, labels, epsilon=0.05, a=None, b=None,
                 cell_ids=None, base_tol=1e-11):
        cost = np.array(cost, dtype=float, copy=True, order="C")
        if cost.ndim != 2 or not cost.size or not np.isfinite(cost).all():
            raise ValueError("cost must be a finite nonempty matrix")
        if np.any(cost < 0):
            raise ValueError("cost must contain nonnegative distances")
        self.n_query, self.n_reference = cost.shape
        labels = _strings(labels, "labels", self.n_reference)
        if len(np.unique(labels)) < 2:
            raise ValueError("reference-composition analysis requires at least two classes")
        if cell_ids is None:
            cell_ids = [f"cell_{j + 1}" for j in range(self.n_query)]
        self._cell_ids = _strings(cell_ids, "cell_ids", self.n_query, unique=True)
        a = None if a is None else np.array(a, dtype=float, copy=True)
        b = None if b is None else np.array(b, dtype=float, copy=True)
        self._fit = CompositionAudit(cost, labels, epsilon=float(epsilon),
                                     a=a, b=b, base_tol=float(base_tol))
        self.epsilon = float(epsilon)
        self._jacobian = None
        self._boundary_cache = None

    @property
    def classes(self):
        return self._fit.classes.copy()

    @property
    def cell_ids(self):
        return self._cell_ids.copy()

    @property
    def baseline_proportions(self):
        return self._fit.p0.copy()

    @property
    def baseline_probabilities(self):
        return self._fit.base_probabilities.copy()

    def describe(self):
        return dict(n_query=self.n_query, n_reference=self.n_reference,
                    n_classes=len(self.classes), classes=self.classes,
                    baseline_proportions=self.baseline_proportions,
                    cell_ids=self.cell_ids, epsilon=self.epsilon,
                    baseline_residual=float(self._fit.base_residual),
                    baseline_solver=self._fit.base_solver)

    def composition_jacobian(self):
        """Return dP/dz for p(z)=softmax(log(p0)+z), shape (n,q,q).

        Axis 1 is the changed log-composition coordinate, axis 2 the output
        class. Centering makes the constant-log-weight null direction explicit.
        """
        if self._jacobian is None:
            p0 = self._fit.p0
            Q = np.diag(p0) - np.outer(p0, p0)
            H = self._fit.local_sensitivity()
            J = np.einsum("ijk,jl->ilk", H, Q)
            J -= J.mean(axis=1, keepdims=True)
            self._jacobian = J
        return self._jacobian.copy()

    def _boundaries(self):
        if self._boundary_cache is not None:
            return self._boundary_cache
        P = self._fit.base_probabilities
        J = self.composition_jacobian()
        row = np.arange(self.n_query)
        winner = P.argmax(axis=1)
        gaps = P[row, winner, None] - P
        gradients = J[row, :, winner][:, :, None] - J
        lengths = np.linalg.norm(gradients, axis=1)
        radii = np.full_like(gaps, np.inf)
        np.divide(np.maximum(gaps, 0), lengths, out=radii, where=lengths > 1e-14)
        radii[row, winner] = np.inf
        competitor = radii.argmin(axis=1)
        radius = radii[row, competitor]
        valid = np.isfinite(radius)
        directions = np.zeros((self.n_query, len(self.classes)))
        directions[valid] = (-gradients[row[valid], :, competitor[valid]] /
                             lengths[row[valid], competitor[valid], None])
        directions -= directions.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(directions, axis=1)
        valid &= norms > 1e-14
        directions[valid] /= norms[valid, None]
        radius[~valid] = np.inf
        self._boundary_cache = (winner, competitor, radius, directions, valid)
        return self._boundary_cache

    def review_cells(self, top_n=None):
        """Rank cells by a local linearized boundary radius, smallest first.

        This radius is an approximation in centered log-weight L2 space, not a
        certified distance to an exact boundary and not biological error risk.
        Ties retain input cell order. Nonresponsive cells have infinite radius.
        """
        if top_n is None:
            top_n = self.n_query
        top_n = min(_integer(top_n, "top_n"), self.n_query)
        winner, competitor, radius, _, valid = self._boundaries()
        order = np.argsort(radius, kind="stable")[:top_n]
        P = self._fit.base_probabilities
        sorted_P = np.sort(P, axis=1)
        competitors = self.classes[competitor].astype(object)
        competitors[~valid] = ""
        return dict(rank=np.arange(1, len(order) + 1), cell_index=order,
                    cell_id=self.cell_ids[order], baseline_label=self.classes[winner[order]],
                    competitor=competitors[order], linearized_radius=radius[order],
                    confidence=P.max(axis=1)[order],
                    winning_margin=(sorted_P[:, -1] - sorted_P[:, -2])[order])

    def _proportions(self, proportions):
        ps = np.asarray(proportions, dtype=float)
        if ps.ndim == 1:
            ps = ps[None, :]
        if ps.ndim != 2 or ps.shape[1] != len(self.classes) or not len(ps):
            raise ValueError("proportions must have one column per fit.classes")
        if not np.isfinite(ps).all() or np.any(ps <= 0):
            raise ValueError("every class proportion must be finite and strictly positive")
        if not np.allclose(ps.sum(axis=1), 1, rtol=0, atol=1e-12):
            raise ValueError("each composition must sum to one")
        return ps.copy()

    def audit_compositions(self, proportions, order=0, corrector="chord", tol=1e-7,
                           batch_size=32, return_probabilities=False):
        """Correct all scenarios, retain fallbacks, and summarize label changes.

        Chunks limit temporary scenario storage. The cached core correction
        system is reused across chunks. Set return_probabilities=True to also
        retain the full (query,scenario,class) array.
        """
        ps = self._proportions(proportions)
        batch_size = _integer(batch_size, "batch_size")
        order = _integer(order, "order", minimum=0)
        if order not in (0, 1, 2):
            raise ValueError("order must be 0, 1 or 2")
        if corrector not in ("chord", "aa2"):
            raise ValueError("corrector must be 'chord' or 'aa2'")
        n, B, q = self.n_query, len(ps), len(self.classes)
        baseline = self._fit.base_probabilities
        base_index = baseline.argmax(axis=1)
        hits = np.zeros(n, dtype=int)
        max_change = np.zeros(n)
        lower, upper = np.full((n, q), np.inf), np.full((n, q), -np.inf)
        saved = np.empty((n, B, q)) if return_probabilities else None
        residuals, iterations = np.empty(B), np.empty(B, dtype=int)
        fallback_indices, evaluations, backtracks = [], [], []
        predictor_fallback_count = 0
        setup_seconds = 0.0
        started = time.perf_counter()
        for start in range(0, B, batch_size):
            stop = min(start + batch_size, B)
            result = self._fit.run(ps[start:stop], order=order, corrector=corrector, tol=float(tol))
            P = result["probabilities"]
            hits += (P.argmax(axis=2) != base_index[:, None]).sum(axis=1)
            max_change = np.maximum(max_change, np.abs(P - baseline[:, None, :]).max(axis=(1, 2)))
            lower = np.minimum(lower, P.min(axis=1))
            upper = np.maximum(upper, P.max(axis=1))
            if saved is not None:
                saved[:, start:stop] = P
            residuals[start:stop] = result["marginal_residual"]
            iterations[start:stop] = result["iterations"]
            d = result["diagnostics"]
            fallback_indices.extend(start + int(i) for i in d.get("fallback_indices", []))
            evaluations.extend(d.get("evaluations", []))
            backtracks.extend(d.get("backtracks", []))
            predictor_fallback_count += int(d.get("predictor_fallback_count", 0))
            setup_seconds += float(d["setup_seconds"])
        return dict(classes=self.classes, proportions=ps, probabilities=saved,
                    baseline_probabilities=baseline.copy(), marginal_residual=residuals,
                    iterations=iterations,
                    cells=dict(cell_id=self.cell_ids, baseline_label=self.classes[base_index],
                               switch_count=hits, switch_frequency=hits / B,
                               max_probability_change=max_change,
                               probability_span=(upper - lower).max(axis=1)),
                    diagnostics=dict(n_scenarios=B, fallback_count=len(fallback_indices),
                                     fallback_indices=np.array(fallback_indices, dtype=int),
                                     predictor_fallback_count=predictor_fallback_count,
                                     evaluations=np.array(evaluations), backtracks=np.array(backtracks),
                                     setup_seconds=setup_seconds,
                                     total_seconds=time.perf_counter() - started,
                                     baseline_seconds=float(self._fit.base_seconds)))

    def propose_scenarios(self, budget=32, radius=0.1, metric="log_weight",
                          method="directed", seed=1, anchors=None):
        """Generate directions; proposals alone do not establish label changes.

        metric='log_weight' sets a centered-log-weight L2 radius.
        metric='tv' sets actual total variation (0.01 means 1%). TV matching
        scales the existing direction; it does not optimize under a TV norm.
        A requested TV unreachable along a direction raises an error.
        """
        budget = _integer(budget, "budget")
        seed = _integer(seed, "seed", minimum=0)
        radius = float(radius)
        if not np.isfinite(radius) or radius <= 0:
            raise ValueError("radius must be finite and positive")
        if metric not in ("log_weight", "tv") or method not in ("directed", "random"):
            raise ValueError("use metric='log_weight'/'tv' and method='directed'/'random'")
        if metric == "tv" and radius >= 1:
            raise ValueError("TV must be less than one; use 0.01 for 1 percent")
        rng = np.random.default_rng(seed)
        if anchors is None:
            anchors = rng.choice(self.n_query, min(budget, self.n_query), replace=False)
        else:
            anchors = np.asarray(anchors)
            if (anchors.ndim != 1 or not len(anchors) or anchors.dtype.kind not in "iu"
                    or np.any(anchors < 0) or np.any(anchors >= self.n_query)
                    or len(np.unique(anchors)) != len(anchors)):
                raise ValueError("anchors must be unique zero-based integer query indices")
        if method == "directed":
            _, _, _, directions, valid = self._boundaries()
            if not valid[anchors].all():
                bad = self.cell_ids[anchors[~valid[anchors]]].tolist()
                raise ValueError(f"No nonzero local boundary direction for anchors {bad}; use random proposals")
            directions = directions[anchors].copy()
        else:
            directions = rng.normal(size=(len(anchors), len(self.classes)))
            directions -= directions.mean(axis=1, keepdims=True)
            norms = np.linalg.norm(directions, axis=1)
            if np.any(norms <= 1e-14):
                raise RuntimeError("Degenerate random direction; try another seed")
            directions /= norms[:, None]
        logp = np.log(self._fit.p0)
        if metric == "log_weight":
            scales = np.full(len(anchors), radius)
            ps = softmax(logp + scales[:, None] * directions, axis=1)
        else:
            ps = np.empty_like(directions)
            scales = np.empty(len(directions))
            for i, direction in enumerate(directions):
                max_mask = direction == direction.max()
                attainable = 1 - self._fit.p0[max_mask].sum()
                if radius >= attainable - 1e-12:
                    raise ValueError(f"Requested TV is unreachable along proposal {i}; reduce radius")
                low, high = 0.0, 1.0
                for _ in range(60):
                    p = softmax(logp + high * direction)
                    if 0.5 * np.abs(p - self._fit.p0).sum() >= radius:
                        break
                    high *= 2
                else:
                    raise ValueError(f"Could not bracket TV for proposal {i}")
                for _ in range(80):
                    mid = (low + high) / 2
                    p = softmax(logp + mid * direction)
                    if 0.5 * np.abs(p - self._fit.p0).sum() < radius:
                        low = mid
                    else:
                        high = mid
                scales[i] = (low + high) / 2
                ps[i] = softmax(logp + scales[i] * direction)
        if np.any(ps <= 0):
            raise ValueError("Proposed weights underflowed to zero; reduce radius")
        return dict(classes=self.classes, proportions=ps, directions=directions,
                    anchor_index=anchors.copy(), anchor_id=self.cell_ids[anchors],
                    log_weight_scale=scales, actual_tv=0.5 * np.abs(ps - self._fit.p0).sum(axis=1),
                    metric=metric, radius=radius, method=method, seed=seed)

    def stress_test(self, budget=32, radius=0.01, metric="tv", seed=1,
                    order=1, tol=1e-7, batch_size=32, return_probabilities=False):
        """Compare directed/random proposals at a common budget and radius.

        Returns both all-cell union detection and paired intended-anchor hits.
        If metric is log_weight, inspect actual_tv: the two methods need not
        induce equally large changes in actual class proportions.
        """
        directed = self.propose_scenarios(budget, radius, metric, "directed", seed)
        random = self.propose_scenarios(budget, radius, metric, "random", seed)
        if not np.array_equal(directed["anchor_index"], random["anchor_index"]):
            raise RuntimeError("Internal anchor pairing mismatch")
        B = len(directed["proportions"])
        result = self.audit_compositions(np.vstack([directed["proportions"], random["proportions"]]),
                                         order=order, tol=tol, batch_size=batch_size,
                                         return_probabilities=True)
        base = self._fit.base_probabilities.argmax(axis=1)
        switched = result["probabilities"].argmax(axis=2) != base[:, None]
        anchors = directed["anchor_index"]
        anchor_directed = switched[anchors, np.arange(B)]
        anchor_random = switched[anchors, B + np.arange(B)]
        comparison = dict(method=np.array(["directed", "random"]),
                          n_scenarios=np.array([B, B]),
                          detected_cells=np.array([switched[:, :B].any(axis=1).sum(),
                                                   switched[:, B:].any(axis=1).sum()]),
                          detected_fraction=np.array([switched[:, :B].any(axis=1).mean(),
                                                      switched[:, B:].any(axis=1).mean()]),
                          anchor_hits=np.array([anchor_directed.sum(), anchor_random.sum()]),
                          anchor_hit_fraction=np.array([anchor_directed.mean(), anchor_random.mean()]),
                          median_actual_tv=np.array([np.median(directed["actual_tv"]),
                                                     np.median(random["actual_tv"])]))
        if not return_probabilities:
            result["probabilities"] = None
        return dict(comparison=comparison, directed=directed, random=random,
                    anchors=dict(cell_index=anchors, cell_id=self.cell_ids[anchors],
                                 directed_hit=anchor_directed, random_hit=anchor_random), audit=result)
