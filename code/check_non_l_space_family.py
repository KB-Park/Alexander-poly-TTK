"""Check the coefficient obstruction for the non-L-space theorem."""

from __future__ import annotations

import argparse

from check_fiberedness_end_families import alexander_from_formula


def check(k: int):
    p = 3 * k + 1
    q = p - 1
    r = p + 2
    coefficients = alexander_from_formula(p, q, r, -1)
    obstructing_coefficient = coefficients[p + 6]
    breadth = len(coefficients) - 1
    if obstructing_coefficient != -3:
        raise AssertionError(
            f"unexpected coefficient for k={k}: {obstructing_coefficient}"
        )
    if breadth != 4 * p - 2:
        raise AssertionError(f"unexpected breadth for k={k}: {breadth}")
    return {
        "k": k,
        "parameters": (p, q, r, -1),
        "coefficient_index": p + 6,
        "coefficient": obstructing_coefficient,
        "alexander_breadth": breadth,
        "l_space_coefficient_obstruction": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k-max", type=int, default=20)
    args = parser.parse_args()
    for k in range(2, args.k_max + 1):
        print(check(k))


if __name__ == "__main__":
    main()
