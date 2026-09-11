"""Independent total-count checks for unordered full k-ary tree shapes.

Based on Andres (2010), Ars Combinatoria 94, 465-469, Theorem 1:
https://combinatorialpress.com/article/ars/Volume%20094/volume-94-paper-45.pdf

Uses the equivalent multiset factor comb(t+r-1, r), not Sackin polynomials.
The returned sequence is indexed by INTERNAL nodes p; leaves n=1+(k-1)*p.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path
import sys


def andres_sequence(k: int, p_max: int) -> list[int]:
    """Return totals for p=0,...,p_max, where p is the internal-node count."""
    if not isinstance(k, int) or k < 2:
        raise ValueError("k must be an integer >= 2")
    if not isinstance(p_max, int) or p_max < 0:
        raise ValueError("p_max must be an integer >= 0")
    totals = [1]  # The single leaf, p=0.
    for p in range(1, p_max + 1):
        target = p - 1  # Internal nodes available below the root.
        limit = min(k, target)
        # dp[c][s]: multisets with c non-leaf children, s internal nodes.
        dp = [[0] * (target + 1) for _ in range(limit + 1)]
        dp[0][0] = 1
        for q in range(1, target + 1):
            t = totals[q]  # Available shapes with q internal nodes.
            updated = [row.copy() for row in dp]  # Choose zero of this size.
            for c in range(limit + 1):
                for s, ways in enumerate(dp[c]):
                    if not ways:
                        continue
                    max_r = min(limit - c, (target - s) // q)
                    for r in range(1, max_r + 1):
                        updated[c + r][s + r * q] += (
                            ways * math.comb(t + r - 1, r)
                        )
            dp = updated
        # Fill the remaining root-child positions with identical leaves.
        totals.append(sum(dp[c][target] for c in range(limit + 1)))
    return totals


def load_analysis(path: Path):
    """Load the explicit script, avoiding a same-named installed package."""
    spec = importlib.util.spec_from_file_location('_andres_sackin_target', path)
    if spec is None or spec.loader is None:
        raise ImportError(f'Cannot load analysis script: {path}')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.SackinAnalysis()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ks', type=int, nargs='+', default=[2, 3, 4, 5, 6])
    parser.add_argument('--p-max', type=int, default=15,
                        help='Maximum internal-node count (default: 15).')
    parser.add_argument('--analysis', type=Path,
                        default=Path(__file__).with_name('sackin_analysis.py'),
                        help='Analysis script to check (default: sibling sackin_analysis.py).')
    parser.add_argument('--counts-only', action='store_true',
                        help='Print sequences without importing or running the analysis.')
    args = parser.parse_args()
    if args.p_max < 0 or any(k < 2 for k in args.ks):
        parser.error('Require p-max >= 0 and every k >= 2.')
    analysis = None if args.counts_only else load_analysis(args.analysis)
    checked = 0
    for k in dict.fromkeys(args.ks):
        totals = andres_sequence(k, args.p_max)
        if args.counts_only:
            print(f'k={k}, p=0..{args.p_max}: {totals}', flush=True)
            continue
        for p, expected in enumerate(totals):
            n = 1 + (k - 1) * p
            freq = analysis.frequencies(n, k)
            actual = sum(int(v) for v in freq['absolute_frequency'])
            if actual != expected:
                raise AssertionError(
                    f'n={n}, k={k}, p={p}: Sackin total={actual}, '
                    f'Andres-based total={expected}'
                )
            checked += 1
        print(f'PASS k={k}: p=0..{args.p_max}', flush=True)
    if not args.counts_only:
        print(f'All {checked} Andres-based total-count comparisons passed.')


if __name__ == '__main__':
    main()