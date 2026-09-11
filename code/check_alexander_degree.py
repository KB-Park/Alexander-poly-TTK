"""Verify the Alexander breadth formula against the manuscript polynomial."""

from __future__ import annotations

import argparse
from math import gcd

from check_fiberedness_end_families import alexander_from_formula


def predicted_breadth(p: int, q: int, r: int, s: int):
    if s >= 0:
        return (p - 1) * (q - 1) + s * r * (r - 1)
    return (-s) * r * (r - 1) - (p - 1) * (q - 1) - 2 * (r - 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-max", type=int, default=15)
    parser.add_argument("--s-min", type=int, default=-4)
    parser.add_argument("--s-max", type=int, default=4)
    args = parser.parse_args()

    checked = 0
    for p in range(3, args.p_max + 1):
        for q in range(2, p):
            if gcd(p, q) != 1:
                continue
            for r in range(p + 1, p + q):
                for s in range(args.s_min, args.s_max + 1):
                    coefficients = alexander_from_formula(p, q, r, s)
                    actual = len(coefficients) - 1
                    predicted = predicted_breadth(p, q, r, s)
                    if actual != predicted:
                        raise AssertionError(
                            f"breadth mismatch for {(p, q, r, s)}: "
                            f"{actual} != {predicted}"
                        )
                    checked += 1
    print(f"verified {checked} parameter quadruples")


if __name__ == "__main__":
    main()
