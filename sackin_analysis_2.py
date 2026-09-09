#!/usr/bin/env python3
"""Reproducible Sackin-index computations for unlabelled, non-plane full k-ary trees.

This module consolidates the two R scripts used for the computational section of
"The Sackin Index on Full k-ary Trees".  It provides:

* extremal Sackin values and admissibility checks;
* dynamic programming for Sackin-value multiplicities;
* weighted maximum-likelihood Beta fits on the conditional interior support;
* exact discrete-versus-continuous Kolmogorov distances;
* binary and higher-k tables used in the manuscript;
* the empirical log-log regressions for the binary fitted Beta parameters;
* publication-ready PDF figures and CSV/LaTeX exports.

Multiplicity coefficients use Python arbitrary-precision integers, with exact
polynomial convolution and checked integer division. Probabilities, moments,
Beta fits and CDF comparisons are evaluated in floating point. Exact counting
can be substantially slower than the former floating-point implementation.
The fft_threshold arguments and --fft-threshold option are retained for API
compatibility but are ignored: FFT convolution is never used.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, betainc, digamma

DEFAULT_FFT_THRESHOLD = 1_000_000_000


# -----------------------------------------------------------------------------
# Structural formulas
# -----------------------------------------------------------------------------

def is_admissible(n: int, k: int) -> bool:
    """Return True iff a full k-ary rooted tree with n leaves can exist."""
    return (
        isinstance(n, (int, np.integer))
        and isinstance(k, (int, np.integer))
        and n >= 1
        and k >= 2
        and (n - 1) % (k - 1) == 0
    )


def ceil_log_k(n: int, k: int) -> Tuple[int, int]:
    """Return (d, k**d), where d = ceil(log_k(n)), using integer arithmetic."""
    if n < 1 or k < 2:
        raise ValueError("Require n >= 1 and k >= 2.")
    d = 0
    kd = 1
    while kd < n:
        kd *= k
        d += 1
    return d, kd


def sackin_min(n: int, k: int = 2) -> int:
    """Minimum Sackin Index on T_{n,k}."""
    if not is_admissible(n, k):
        raise ValueError(
            f"No full {k}-ary rooted tree with n={n} leaves exists; "
            f"admissible n satisfy n = 1 mod {k - 1}."
        )
    if n == 1:
        return 0
    d, kd = ceil_log_k(n, k)
    value = n * d - (kd - n) // (k - 1)
    return int(value)


def sackin_max(n: int, k: int = 2) -> int:
    """Maximum Sackin Index on T_{n,k}."""
    if not is_admissible(n, k):
        raise ValueError(
            f"No full {k}-ary rooted tree with n={n} leaves exists; "
            f"admissible n satisfy n = 1 mod {k - 1}."
        )
    if n == 1:
        return 0
    num = (n - 1) * (n + k)
    den = 2 * (k - 1)
    if num % den != 0:
        raise ArithmeticError("Extremal formula unexpectedly produced a non-integer.")
    return num // den


def admissible_n_values(n_max: int, k: int, include_one: bool = False) -> List[int]:
    """Admissible leaf counts up to n_max for a fixed branching degree k."""
    if k < 2:
        raise ValueError("k must be at least 2.")
    start = 1 if include_one else k
    if n_max < start:
        return []
    return list(range(start, n_max + 1, k - 1))


# -----------------------------------------------------------------------------
# Polynomial utilities
# -----------------------------------------------------------------------------

def integer_coefficients(values: np.ndarray) -> np.ndarray:
    """Validate coefficients and promote fixed-width integers to Python ints."""
    values = np.asarray(values, dtype=object)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("Coefficient arrays must be nonempty and one-dimensional.")
    result = np.empty(values.size, dtype=object)
    for i, value in enumerate(values):
        if not isinstance(value, (int, np.integer)):
            raise TypeError("Coefficients must be integers, not floating-point values.")
        if value < 0:
            raise ValueError("Counting coefficients must be nonnegative.")
        result[i] = int(value)
    return result


def poly_convolution(
    a: np.ndarray,
    b: np.ndarray,
    fft_threshold: int = DEFAULT_FFT_THRESHOLD,
) -> np.ndarray:
    """Exact integer convolution; fft_threshold is ignored for compatibility."""
    return np.convolve(integer_coefficients(a), integer_coefficients(b))


def exact_divide(coefficients: np.ndarray, divisor: int) -> np.ndarray:
    """Divide integer coefficients, rejecting any nonzero remainder."""
    if not isinstance(divisor, (int, np.integer)) or divisor <= 0:
        raise ValueError("The divisor must be a positive integer.")
    coefficients = integer_coefficients(coefficients)
    result = np.empty(coefficients.size, dtype=object)
    for i, value in enumerate(coefficients):
        quotient, remainder = divmod(value, int(divisor))
        if remainder:
            raise ArithmeticError(f"Coefficient {i} is not divisible by {divisor}.")
        result[i] = quotient
    return result


def relative_probabilities(counts: np.ndarray) -> np.ndarray:
    """Convert exact counts to float ratios without first converting the total."""
    total = sum(counts, 0)
    if total <= 0:
        raise ValueError("The total count must be positive.")
    return np.array([count / total for count in counts], dtype=float)


def nondecreasing_k_tuples(total: int, k: int) -> Iterator[Tuple[int, ...]]:
    """Yield q_1 <= ... <= q_k, q_i >= 0, with sum q_i = total."""
    if total < 0 or k < 1:
        return

    prefix: List[int] = []

    def rec(remaining: int, parts_left: int, lower: int) -> Iterator[Tuple[int, ...]]:
        if parts_left == 1:
            if remaining >= lower:
                yield tuple(prefix + [remaining])
            return

        upper = remaining // parts_left
        for x in range(lower, upper + 1):
            prefix.append(x)
            yield from rec(remaining - x, parts_left - 1, x)
            prefix.pop()

    yield from rec(int(total), int(k), 0)


def multiset_sackin_poly(
    counts: np.ndarray,
    r: int,
    fft_threshold: int = DEFAULT_FFT_THRESHOLD,
) -> np.ndarray:
    """Sackin polynomial for an unordered multiset of r equal-size subtrees.

    If counts[q] is the number of tree shapes with Sackin Index minS + q,
    the output coefficients count unordered multisets of r such shapes by the
    total relative Sackin increment.  The complete-homogeneous recurrence

        m h_m = sum_{i=1}^m p_i h_{m-i}

    is used, where p_i substitutes z -> z^i. 

    See https://math.berkeley.edu/~corteel/MATH249/macdonald.pdf#page=34
    equation (2.11), page 23
    """
    counts = integer_coefficients(counts)
    r = int(r)
    if r < 0:
        raise ValueError("r must be non-negative.")
    if r == 0:
        return np.array([1], dtype=object)
    if r == 1:
        return counts.copy()

    max_rel = counts.size - 1
    h: List[np.ndarray] = [np.array([1], dtype=object)]

    for m in range(1, r + 1):
        acc = np.zeros(m * max_rel + 1, dtype=object)
        for i in range(1, m + 1):
            p_i = np.zeros(i * max_rel + 1, dtype=object)
            p_i[::i] = counts
            term = poly_convolution(p_i, h[m - i], fft_threshold)
            if term.size < acc.size:
                term = np.pad(term, (0, acc.size - term.size))
            acc += term[: acc.size]

        h.append(exact_divide(acc, m))

    return h[r]


# -----------------------------------------------------------------------------
# Dynamic programming states
# -----------------------------------------------------------------------------

@dataclass
class DPState:
    min_sackin: int
    max_sackin: int
    counts: np.ndarray


class BinarySackinDP:
    """Incremental dynamic program for full binary tree shapes."""

    def __init__(self, fft_threshold: int = DEFAULT_FFT_THRESHOLD) -> None:
        self.fft_threshold = int(fft_threshold)
        self.dp: Dict[int, DPState] = {
            1: DPState(0, 0, np.array([1], dtype=object))
        }
        self.max_built = 1

    def _extend_to(self, n: int) -> None:
        if n <= self.max_built:
            return

        for N in range(self.max_built + 1, n + 1):
            min_n = sackin_min(N, 2)
            max_n = sackin_max(N, 2)
            current = np.zeros(max_n - min_n + 1, dtype=object)

            for i in range(1, N // 2 + 1):
                j = N - i
                state_i = self.dp[i]
                state_j = self.dp[j]

                if i < j:
                    conv = poly_convolution(
                        state_i.counts, state_j.counts, self.fft_threshold
                    )
                    first_s = N + state_i.min_sackin + state_j.min_sackin
                else:
                    ordered = poly_convolution(
                        state_i.counts, state_i.counts, self.fft_threshold
                    )
                    diagonal = np.zeros_like(ordered)
                    diagonal[::2] = state_i.counts
                    conv = exact_divide(ordered + diagonal, 2)
                    first_s = N + 2 * state_i.min_sackin

                start = first_s - min_n
                current[start : start + conv.size] += conv

            self.dp[N] = DPState(min_n, max_n, current)

        self.max_built = n

    def frequencies(self, n: int) -> pd.DataFrame:
        if not isinstance(n, (int, np.integer)) or n < 1:
            raise ValueError("n must be a positive integer.")
        n = int(n)
        self._extend_to(n)
        state = self.dp[n]
        values = np.arange(state.min_sackin, state.max_sackin + 1, dtype=int)
        return pd.DataFrame(
            {
                "Sackin": values,
                "absolute_frequency": state.counts.copy(),
                "relative_frequency": relative_probabilities(state.counts),
            }
        )


class KarySackinDP:
    """Incremental dynamic program for unlabelled, non-plane full k-ary trees."""

    def __init__(self, k: int, fft_threshold: int = DEFAULT_FFT_THRESHOLD) -> None:
        if not isinstance(k, (int, np.integer)) or k < 2:
            raise ValueError("k must be an integer >= 2.")
        self.k = int(k)
        self.fft_threshold = int(fft_threshold)
        self.dp: Dict[int, DPState] = {
            1: DPState(0, 0, np.array([1], dtype=object))
        }
        self.max_p_built = 0
        self.multiset_cache: Dict[Tuple[int, int], np.ndarray] = {}

    def _extend_to_p(self, p_target: int) -> None:
        if p_target <= self.max_p_built:
            return

        k = self.k
        for p in range(self.max_p_built + 1, p_target + 1):
            N = 1 + p * (k - 1)
            min_n = sackin_min(N, k)
            max_n = sackin_max(N, k)
            current = np.zeros(max_n - min_n + 1, dtype=object)

            # Q = (N-k)/(k-1) = p-1.
            for q_vec in nondecreasing_k_tuples(p - 1, k):
                j_vec = tuple(1 + (k - 1) * q for q in q_vec)
                size_table = Counter(j_vec)

                subtree_poly = np.array([1], dtype=object)
                min_subtree_sum = 0

                for j, r in sorted(size_table.items()):
                    child = self.dp[j]
                    key = (j, r)
                    group_poly = self.multiset_cache.get(key)
                    if group_poly is None:
                        group_poly = multiset_sackin_poly(
                            child.counts, r, self.fft_threshold
                        )
                        self.multiset_cache[key] = group_poly

                    subtree_poly = poly_convolution(
                        subtree_poly, group_poly, self.fft_threshold
                    )
                    min_subtree_sum += r * child.min_sackin

                first_s = N + min_subtree_sum
                start = first_s - min_n
                current[start : start + subtree_poly.size] += subtree_poly

            self.dp[N] = DPState(min_n, max_n, current)

        self.max_p_built = p_target

    def frequencies(self, n: int) -> pd.DataFrame:
        if not isinstance(n, (int, np.integer)) or n < 1:
            raise ValueError("n must be a positive integer.")
        n = int(n)
        if not is_admissible(n, self.k):
            raise ValueError(
                f"No full {self.k}-ary rooted tree with n={n} leaves exists; "
                f"admissible n satisfy n = 1 mod {self.k - 1}."
            )

        if n == 1:
            return pd.DataFrame(
                {
                    "Sackin": [0],
                    "absolute_frequency": np.array([1], dtype=object),
                    "relative_frequency": [1.0],
                }
            )

        p = (n - 1) // (self.k - 1)
        self._extend_to_p(p)
        state = self.dp[n]

        values = np.arange(
            state.min_sackin, state.max_sackin + 1, self.k - 1, dtype=int
        )
        indices = values - state.min_sackin
        counts = state.counts[indices]

        return pd.DataFrame(
            {
                "Sackin": values,
                "absolute_frequency": counts.copy(),
                "relative_frequency": relative_probabilities(counts),
            }
        )


class SackinAnalysis:
    """Cache-aware facade providing frequency distributions for any fixed k."""

    def __init__(self, fft_threshold: int = DEFAULT_FFT_THRESHOLD) -> None:
        self.fft_threshold = int(fft_threshold)
        self.binary = BinarySackinDP(fft_threshold)
        self.kary: Dict[int, KarySackinDP] = {}

    def frequencies(self, n: int, k: int) -> pd.DataFrame:
        if k == 2:
            return self.binary.frequencies(n)
        if k not in self.kary:
            self.kary[k] = KarySackinDP(k, self.fft_threshold)
        return self.kary[k].frequencies(n)


# -----------------------------------------------------------------------------
# Beta fit and diagnostics
# -----------------------------------------------------------------------------

@dataclass
class BetaFitResult:
    n: int
    k: int
    internal_nodes: int
    alpha: float
    beta: float
    mu_exact: float
    mu_beta: float
    error_mu: float
    var_exact: float
    var_beta: float
    error_var: float
    D: float
    D_full: float
    x_max: float
    side_max: str
    extreme_mass: float
    optimizer_success: bool
    optimizer_message: str


def _beta_logpdf(x: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    return (
        (alpha - 1.0) * np.log(x)
        + (beta - 1.0) * np.log1p(-x)
        - betaln(alpha, beta)
    )


def beta_fit(
    n: int,
    k: int,
    analysis: SackinAnalysis | None = None,
    plot_cdf: bool = False,
    cdf_path: str | Path | None = None,
) -> BetaFitResult:
    """Fit Beta(alpha,beta) to the conditional interior Sackin distribution.

    The reported exact mean, variance, moment errors, MLE and Kolmogorov
    distance D all use the conditional distribution given 0 < X < 1.
    D_full and extreme_mass are separate full-distribution diagnostics.
    """
    if analysis is None:
        analysis = SackinAnalysis()

    # extract non-zero cases
    freq = analysis.frequencies(int(n), int(k)).copy()
    freq = freq[freq["absolute_frequency"] > 0].reset_index(drop=True)

    min_s = int(freq["Sackin"].min())
    max_s = int(freq["Sackin"].max())
    span = max_s - min_s
    if span == 0:
        raise ValueError(
            f"The Sackin Index is constant for n={n}, k={k}; a Beta fit is undefined."
        )

    x_full = (freq["Sackin"].to_numpy(dtype=float) - min_s) / span
    p_full = freq["relative_frequency"].to_numpy(dtype=float)

    interior = (x_full > 0.0) & (x_full < 1.0)
    x = x_full[interior]
    w_raw = p_full[interior]
    extreme_mass = float(np.sum(p_full[~interior]))

    if x.size < 2 or float(np.sum(w_raw)) <= 0.0:
        raise ValueError(
            f"Not enough interior Sackin values for n={n}, k={k} to fit a Beta distribution."
        )
    # force area 1
    w = w_raw / np.sum(w_raw)
    # compute mean and var
    mu_cond = float(np.sum(w * x))
    var_cond = float(np.sum(w * (x - mu_cond) ** 2))
    # All table moments and errors describe the same conditional fit target.
    mu_exact = mu_cond
    var_exact = var_cond
    lam = mu_cond * (1.0 - mu_cond) / var_cond - 1.0 if var_cond > 0 else -1.0

    if not np.isfinite(lam) or lam <= 0.0:
        alpha0 = beta0 = 2.0
    else:
        alpha0 = mu_cond * lam
        beta0 = (1.0 - mu_cond) * lam

    mean_log_x = float(np.sum(w * np.log(x)))
    mean_log_1mx = float(np.sum(w * np.log1p(-x)))

    def objective(log_par: np.ndarray) -> float:
        alpha, beta = np.exp(log_par)
        value = (
            -(alpha - 1.0) * mean_log_x
            -(beta - 1.0) * mean_log_1mx
            + betaln(alpha, beta)
        )
        return float(value) if np.isfinite(value) else 1e100

    def gradient(log_par: np.ndarray) -> np.ndarray:
        alpha, beta = np.exp(log_par)
        common = digamma(alpha + beta)
        d_alpha = -mean_log_x + digamma(alpha) - common
        d_beta = -mean_log_1mx + digamma(beta) - common
        # Chain rule for eta=(log alpha, log beta).
        return np.array([alpha * d_alpha, beta * d_beta], dtype=float)

    opt = minimize(
        objective,
        x0=np.log([alpha0, beta0]),
        jac=gradient,
        method="BFGS",
        options={"gtol": 1e-8, "maxiter": 10_000},
    )
    alpha_hat, beta_hat = np.exp(opt.x)
    alpha_hat = float(alpha_hat)
    beta_hat = float(beta_hat)

    mu_beta = alpha_hat / (alpha_hat + beta_hat)
    var_beta = (
        alpha_hat
        * beta_hat
        / ((alpha_hat + beta_hat) ** 2 * (alpha_hat + beta_hat + 1.0))
    )

    F_right = np.cumsum(w)
    F_left = np.concatenate(([0.0], F_right[:-1]))
    F_beta = betainc(alpha_hat, beta_hat, x)
    diff_right = np.abs(F_right - F_beta)
    diff_left = np.abs(F_left - F_beta)

    d_right = float(np.max(diff_right))
    d_left = float(np.max(diff_left))
    if d_right >= d_left:
        idx = int(np.argmax(diff_right))
        D = d_right
        side = "right"
    else:
        idx = int(np.argmax(diff_left))
        D = d_left
        side = "left"
    x_max = float(x[idx])

    F_full_right = np.cumsum(p_full)
    F_full_left = np.concatenate(([0.0], F_full_right[:-1]))
    F_beta_full = betainc(alpha_hat, beta_hat, x_full)
    D_full = float(
        max(
            np.max(np.abs(F_full_right - F_beta_full)),
            np.max(np.abs(F_full_left - F_beta_full)),
        )
    )

    if plot_cdf:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        xx = np.linspace(0.0, 1.0, 2000)
        ax.plot(xx, betainc(alpha_hat, beta_hat, xx), color="red", label="Beta CDF")
        ax.step(x, F_right, where="post", color="black", label="Exact conditional CDF")
        ax.set_xlabel("Normalized Sackin Index")
        ax.set_ylabel("CDF")
        ax.set_title(f"Full {k}-ary trees: n={n}, D={D:.5f}")
        ax.legend(frameon=False)
        fig.tight_layout()
        if cdf_path is not None:
            fig.savefig(Path(cdf_path), bbox_inches="tight")
            plt.close(fig)

    return BetaFitResult(
        n=int(n),
        k=int(k),
        internal_nodes=(int(n) - 1) // (int(k) - 1),
        alpha=alpha_hat,
        beta=beta_hat,
        mu_exact=mu_exact,
        mu_beta=float(mu_beta),
        error_mu=abs(mu_exact - mu_beta),
        var_exact=var_exact,
        var_beta=float(var_beta),
        error_var=abs(var_exact - var_beta),
        D=D,
        D_full=D_full,
        x_max=x_max,
        side_max=side,
        extreme_mass=extreme_mass,
        optimizer_success=bool(opt.success),
        optimizer_message=str(opt.message),
    )


# -----------------------------------------------------------------------------
# Tables, regressions and figures
# -----------------------------------------------------------------------------

def binary_table(
    analysis: SackinAnalysis,
    ns: Iterable[int] = range(10, 101, 5),
) -> pd.DataFrame:
    rows = []
    for n in ns:
        fit = beta_fit(int(n), 2, analysis)
        rows.append(
            {
                "n": fit.n,
                "alpha": fit.alpha,
                "beta": fit.beta,
                "mu_exact": fit.mu_exact,
                "mu_beta": fit.mu_beta,
                "error_mu": fit.error_mu,
                "var_exact": fit.var_exact,
                "var_beta": fit.var_beta,
                "error_var": fit.error_var,
                "D": fit.D,
            }
        )
    return pd.DataFrame(rows)


def higher_k_table(
    analysis: SackinAnalysis,
    p: int = 20,
    ks: Iterable[int] = range(3, 11),
) -> pd.DataFrame:
    rows = []
    for k in ks:
        n = 1 + (int(k) - 1) * int(p)
        fit = beta_fit(n, int(k), analysis)
        rows.append(
            {
                "p": int(p),
                "k": fit.k,
                "n": fit.n,
                "alpha": fit.alpha,
                "beta": fit.beta,
                "D": fit.D,
                "mu_exact": fit.mu_exact,
                "mu_beta": fit.mu_beta,
                "error_mu": fit.error_mu,
                "var_exact": fit.var_exact,
                "var_beta": fit.var_beta,
                "error_var": fit.error_var,
            }
        )
    return pd.DataFrame(rows)


def higher_k_grid(
    analysis: SackinAnalysis,
    ps: Sequence[int] = (10, 20, 30, 40),
    ks: Sequence[int] = tuple(range(3, 11)),
) -> pd.DataFrame:
    """Compute the complete (p,k) grid efficiently, reusing the DP for each k."""
    rows = []
    for k in ks:
        for p in sorted(ps):
            n = 1 + (int(k) - 1) * int(p)
            fit = beta_fit(n, int(k), analysis)
            rows.append(
                {
                    "p": int(p),
                    "k": fit.k,
                    "n": fit.n,
                    "alpha": fit.alpha,
                    "beta": fit.beta,
                    "D": fit.D,
                    "mu_exact": fit.mu_exact,
                    "mu_beta": fit.mu_beta,
                    "error_mu": fit.error_mu,
                    "var_exact": fit.var_exact,
                    "var_beta": fit.var_beta,
                    "error_var": fit.error_var,
                }
            )
    return pd.DataFrame(rows).sort_values(["p", "k"]).reset_index(drop=True)


def max_error_table(grid: pd.DataFrame) -> pd.DataFrame:
    return (
        grid.groupby("p", as_index=False)
        .agg(
            max_D=("D", "max"),
            max_error_mu=("error_mu", "max"),
            max_error_var=("error_var", "max"),
        )
        .sort_values("p")
        .reset_index(drop=True)
    )


def power_law_regressions(binary: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    """Fit log(parameter) = intercept + exponent * log(n)."""
    x = np.log(binary["n"].to_numpy(dtype=float))
    out: Dict[str, Dict[str, float]] = {}
    for name in ("alpha", "beta"):
        y = np.log(binary[name].to_numpy(dtype=float))
        exponent, intercept = np.polyfit(x, y, 1)
        y_hat = intercept + exponent * x
        ss_res = float(np.sum((y - y_hat) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r2 = 1.0 - ss_res / ss_tot
        out[name] = {
            "intercept": float(intercept),
            "exponent": float(exponent),
            "coefficient": float(math.exp(intercept)),
            "r_squared": float(r2),
        }
    return out


def plot_distribution_beta(
    n: int,
    k: int,
    analysis: SackinAnalysis,
    output_path: str | Path,
    title: str | None = None,
) -> BetaFitResult:
    """Plot the exact density-scaled distribution and its fitted Beta density.

    The exact heights use the *complete* probabilities P_{n,k}(s) divided by
    Delta=(k-1)/(M-m), as defined in the manuscript.  The Beta fit itself is
    conditional on 0 < X < 1.
    """
    freq = analysis.frequencies(n, k)
    fit = beta_fit(n, k, analysis)

    min_s = int(freq["Sackin"].min())
    max_s = int(freq["Sackin"].max())
    span = max_s - min_s
    if span <= 0:
        raise ValueError("A non-degenerate Sackin range is required for plotting.")

    x = (freq["Sackin"].to_numpy(dtype=float) - min_s) / span
    p = freq["relative_frequency"].to_numpy(dtype=float)
    delta = (k - 1) / span
    heights = p / delta

    xx = np.linspace(0.0001, 0.9999, 2000)
    yy = np.exp(_beta_logpdf(xx, fit.alpha, fit.beta))

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(x, heights, linewidth=2.0, color="black", label="Exact Sackin distribution")
    ax.plot(
        xx,
        yy,
        linewidth=2.0,
        color="red",
        label=rf"Beta $B({fit.alpha:.3f},{fit.beta:.3f})$",
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("Normalized Sackin Index")
    ax.set_ylabel("Density")
    if title is None:
        title = f"Sackin Index distribution for n = {n}" if k == 2 else f"Full {k}-ary trees, n = {n}"
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(Path(output_path), bbox_inches="tight")
    plt.close(fig)
    return fit


def _fmt_sci(x: float, digits: int = 2) -> str:
    if x == 0:
        return "0"
    exponent = int(math.floor(math.log10(abs(x))))
    mantissa = x / (10 ** exponent)
    return rf"${mantissa:.{digits}f}\times10^{{{exponent}}}$"


def write_latex_tables(
    binary: pd.DataFrame,
    higher: pd.DataFrame,
    max_errors: pd.DataFrame,
    outdir: Path,
) -> None:
    """Write LaTeX tables; all moments and errors use 0 < X < 1."""
    lines = [
        r"% All moments and errors refer to the conditional distribution 0 < X < 1.",
        r"\begin{tabular}{@{}rrrrrrrrrr@{}}",
        r"\hline",
        r"$n$ & $\widehat\alpha_n$ & $\widehat\beta_n$ & $\mu_n$ & $\mu_n^B$ & $e_{\mu,n}$ & $\sigma_n^2$ & $(\sigma_n^B)^2$ & $e_{\sigma^2,n}$ & $D_n$ \\",
        r"\hline",
    ]
    for _, r in binary.iterrows():
        lines.append(
            f"{int(r.n)} & {r.alpha:.6f} & {r.beta:.6f} & {r.mu_exact:.6f} & "
            f"{r.mu_beta:.6f} & {_fmt_sci(r.error_mu)} & {r.var_exact:.6f} & "
            f"{r.var_beta:.6f} & {_fmt_sci(r.error_var)} & {r.D:.6f} \\\\" 
        )
    lines += [r"\hline", r"\end{tabular}"]
    (outdir / "table_binary.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    lines = [
        r"% All moments and errors refer to the conditional distribution 0 < X < 1.",
        r"\begin{tabular}{@{}rrrrrrrrrrr@{}}",
        r"\hline",
        r"$k$ & $n$ & $\widehat\alpha$ & $\widehat\beta$ & $D_{p,k}$ & $\mu$ & $\mu^B$ & $e_\mu$ & $\sigma^2$ & $(\sigma^B)^2$ & $e_{\sigma^2}$ \\",
        r"\hline",
    ]
    for _, r in higher.iterrows():
        lines.append(
            f"{int(r.k)} & {int(r.n)} & {r.alpha:.4f} & {r.beta:.4f} & {r.D:.5f} & "
            f"{r.mu_exact:.5f} & {r.mu_beta:.5f} & {_fmt_sci(r.error_mu)} & "
            f"{r.var_exact:.5f} & {r.var_beta:.5f} & {_fmt_sci(r.error_var)} \\\\" 
        )
    lines += [r"\hline", r"\end{tabular}"]
    (outdir / "table_higher_k_p20.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    lines = [
        r"% All moments and errors refer to the conditional distribution 0 < X < 1.",
        r"\begin{tabular}{rrrr}",
        r"\hline",
        r"$p$ & $\max D_{p,k}$ & $\max e_{\mu,p,k}$ & $\max e_{\sigma^2,p,k}$ \\",
        r"\hline",
    ]
    for _, r in max_errors.iterrows():
        lines.append(
            f"{int(r.p)} & {r.max_D:.5f} & {_fmt_sci(r.max_error_mu)} & "
            f"{_fmt_sci(r.max_error_var)} \\\\" 
        )
    lines += [r"\hline", r"\end{tabular}"]
    (outdir / "table_max_errors.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def reproduce_paper(
    output_dir: str | Path,
    fft_threshold: int = DEFAULT_FFT_THRESHOLD,
) -> None:
    """Generate the tables, regressions and figures used in the manuscript."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    analysis = SackinAnalysis(fft_threshold=fft_threshold)

    print("[1/5] Binary table n=10,15,...,100")
    binary = binary_table(analysis)
    binary.to_csv(outdir / "binary_beta_table.csv", index=False)

    regressions = power_law_regressions(binary)
    (outdir / "binary_power_law_regressions.json").write_text(
        json.dumps(regressions, indent=2) + "\n", encoding="utf-8"
    )

    print("[2/5] Binary figures n=10,15,20,25")
    for n in (10, 15, 20, 25):
        plot_distribution_beta(n, 2, analysis, outdir / f"Sackin{n}.pdf")

    print("[3/5] Higher-k table p=20, k=3,...,10")
    higher = higher_k_table(analysis, p=20, ks=range(3, 11))
    higher.to_csv(outdir / "higher_k_p20_table.csv", index=False)

    print("[4/5] Higher-k figures p=20, k=3,4,5,6")
    for k in (3, 4, 5, 6):
        n = 1 + (k - 1) * 20
        plot_distribution_beta(n, k, analysis, outdir / f"Sackin{n}_{k}.pdf")

    print("[5/5] Higher-k grid p in {10,20,30,40}, k=3,...,10")
    grid = higher_k_grid(analysis)
    grid.to_csv(outdir / "higher_k_grid.csv", index=False)
    maxima = max_error_table(grid)
    maxima.to_csv(outdir / "higher_k_max_errors.csv", index=False)

    write_latex_tables(binary, higher, maxima, outdir)
    print(f"Done. Results written to: {outdir.resolve()}")


# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sackin Index multiplicities and Beta approximations for full k-ary tree shapes."
    )
    parser.add_argument(
        "--fft-threshold",
        type=int,
        default=DEFAULT_FFT_THRESHOLD,
        help="Ignored compatibility option; counting always uses exact integer convolution.",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    p_freq = sub.add_parser("frequencies", help="Compute the exact multiplicity distribution.")
    p_freq.add_argument("--n", type=int, required=True)
    p_freq.add_argument("--k", type=int, default=2)
    p_freq.add_argument("--output", type=Path)

    p_fit = sub.add_parser("fit", help="Fit a Beta distribution and report diagnostics.")
    p_fit.add_argument("--n", type=int, required=True)
    p_fit.add_argument("--k", type=int, default=2)
    p_fit.add_argument("--cdf", type=Path, help="Optional PDF/PNG path for a CDF comparison plot.")

    p_plot = sub.add_parser("plot", help="Plot density-scaled exact distribution and Beta fit.")
    p_plot.add_argument("--n", type=int, required=True)
    p_plot.add_argument("--k", type=int, default=2)
    p_plot.add_argument("--output", type=Path, required=True)

    p_binary = sub.add_parser("binary-table", help="Generate the binary table n=10,15,...,100.")
    p_binary.add_argument("--output", type=Path, default=Path("binary_beta_table.csv"))

    p_higher = sub.add_parser("higher-table", help="Generate a fixed-p higher-k table.")
    p_higher.add_argument("--p", type=int, default=20)
    p_higher.add_argument("--k-min", type=int, default=3)
    p_higher.add_argument("--k-max", type=int, default=10)
    p_higher.add_argument("--output", type=Path, default=Path("higher_k_table.csv"))

    p_grid = sub.add_parser("grid", help="Generate the p x k grid and maximum-error summary.")
    p_grid.add_argument("--ps", type=int, nargs="+", default=[10, 20, 30, 40])
    p_grid.add_argument("--k-min", type=int, default=3)
    p_grid.add_argument("--k-max", type=int, default=10)
    p_grid.add_argument("--output", type=Path, default=Path("higher_k_grid.csv"))
    p_grid.add_argument("--max-output", type=Path, default=Path("higher_k_max_errors.csv"))

    p_paper = sub.add_parser("paper", help="Reproduce all computational outputs used in the paper.")
    p_paper.add_argument("--output-dir", type=Path, default=Path("results"))

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    analysis = SackinAnalysis(fft_threshold=args.fft_threshold)

    if args.command == "frequencies":
        df = analysis.frequencies(args.n, args.k)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(args.output, index=False)
        else:
            print(df.to_string(index=False))

    elif args.command == "fit":
        result = beta_fit(
            args.n,
            args.k,
            analysis,
            plot_cdf=args.cdf is not None,
            cdf_path=args.cdf,
        )
        print(json.dumps(asdict(result), indent=2))

    elif args.command == "plot":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result = plot_distribution_beta(args.n, args.k, analysis, args.output)
        print(json.dumps(asdict(result), indent=2))

    elif args.command == "binary-table":
        df = binary_table(analysis)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False)
        print(df.to_string(index=False))

    elif args.command == "higher-table":
        df = higher_k_table(analysis, args.p, range(args.k_min, args.k_max + 1))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False)
        print(df.to_string(index=False))

    elif args.command == "grid":
        grid = higher_k_grid(
            analysis,
            ps=tuple(args.ps),
            ks=tuple(range(args.k_min, args.k_max + 1)),
        )
        maxima = max_error_table(grid)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.max_output.parent.mkdir(parents=True, exist_ok=True)
        grid.to_csv(args.output, index=False)
        maxima.to_csv(args.max_output, index=False)
        print(maxima.to_string(index=False))

    elif args.command == "paper":
        reproduce_paper(args.output_dir, fft_threshold=args.fft_threshold)


if __name__ == "__main__":
    main()
