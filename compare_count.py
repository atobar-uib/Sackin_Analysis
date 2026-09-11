import argparse

from check_andres import andres_sequence
from sackin_analysis import SackinAnalysis


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare Sackin multiplicity totals with Andres-based counts."
    )
    parser.add_argument(
        "--ks", type=int, nargs="+", default=list(range(2, 11)),
        help="Branching degrees to check (default: 2 through 10).",
    )
    parser.add_argument(
        "--p-max", type=int, default=15,
        help="Maximum internal-node count, inclusive (default: 15).",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print both counts for every case.",
    )
    args = parser.parse_args()

    if any(k < 2 for k in args.ks):
        parser.error("Every k must be at least 2.")
    if args.p_max < 0:
        parser.error("--p-max must be nonnegative.")

    analysis = SackinAnalysis()
    checked = 0

    for k in dict.fromkeys(args.ks):
        reference = andres_sequence(k, args.p_max)

        for p, expected in enumerate(reference):
            n = 1 + (k - 1) * p
            frequencies = analysis.frequencies(n, k)
            actual = sum(int(v) for v in frequencies["absolute_frequency"])

            if actual != expected:
                raise AssertionError(
                    f"k={k}, n={n}, p={p}: "
                    f"our count={actual}, Andres count={expected}"
                )

            if args.verbose:
                print(
                    f"k={k}, n={n}, p={p}: "
                    f"ours={actual}, Andres={expected} — PASS"
                )

            checked += 1

        print(f"PASS: k={k}, internal nodes p=0..{args.p_max}")

    print(f"All {checked} comparisons passed.")


if __name__ == "__main__":
    main()