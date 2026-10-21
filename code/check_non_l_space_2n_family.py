#!/usr/bin/env python3
"""Check the family T(2n+1,2n;3n,s).

For negative s=-u, n >= 3, the normalized Alexander polynomial has
coefficient 2 at exponent 4(3u-2)n.  Positive twists have two further
uniform coefficient obstructions; see ``obstructing_term`` below.
"""

from __future__ import annotations

import argparse
from typing import Optional

from alexander_polynomial import alexander_coefficients


def obstructing_term(n: int, s: int) -> Optional[tuple[int, int]]:
    """Return a proved coefficient obstruction, when one applies."""
    if s < 0 and n >= 3:
        u = -s
        return 4 * (3 * u - 2) * n, 2
    if s >= 2 and n >= 4:
        return (15 * s + 9) * n + 3, -2
    if s == 1 and n >= 5:
        return 29 * n + 4, -2
    return None


def expected_breadth(n: int, s: int) -> int:
    """Specialize the paper's breadth formula to this family."""
    if s >= 0:
        return (9 * s + 4) * n * n - (3 * s + 2) * n
    u = -s
    return (9 * u - 4) * n * n - (3 * u + 4) * n + 2


def check(n: int, s: int = -1) -> dict[str, object]:
    """Return, and internally verify, the data for one pair (n,s)."""
    if n < 1:
        raise ValueError("n must be positive")

    parameters = (2 * n + 1, 2 * n, 3 * n, s)
    coefficients = alexander_coefficients(*parameters)
    breadth = len(coefficients) - 1
    predicted_breadth = expected_breadth(n, s)
    if breadth != predicted_breadth:
        raise AssertionError(
            f"unexpected breadth for (n,s)={(n, s)}: "
            f"{breadth} != {predicted_breadth}"
        )

    obstruction = obstructing_term(n, s)
    first_obstructing_index = next(
        (index for index, value in enumerate(coefficients) if abs(value) > 1),
        None,
    )
    coefficient_index = obstruction[0] if obstruction else None
    coefficient = coefficients[coefficient_index] if obstruction else None
    if obstruction and coefficient != obstruction[1]:
        raise AssertionError(
            f"unexpected coefficient for (n,s)={(n, s)}: "
            f"[t^{coefficient_index}] Delta = {coefficient}"
        )
    if obstruction and first_obstructing_index != coefficient_index:
        raise AssertionError(
            f"unexpected first obstructing index for (n,s)={(n, s)}: "
            f"{first_obstructing_index} != {coefficient_index}"
        )

    return {
        "n": n,
        "s": s,
        "parameters": parameters,
        "alexander_breadth": breadth,
        "coefficient_index": coefficient_index,
        "coefficient": coefficient,
        "first_obstructing_index": first_obstructing_index,
        "maximum_coefficient_absolute_value": max(map(abs, coefficients)),
        "l_space_coefficient_obstruction": obstruction is not None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-min", type=int, default=1)
    parser.add_argument("--n-max", type=int, default=100)
    parser.add_argument("--s-min", type=int, default=-1)
    parser.add_argument("--s-max", type=int, default=-1)
    args = parser.parse_args()
    if args.n_min < 1 or args.n_max < args.n_min:
        parser.error("require 1 <= n-min <= n-max")
    if args.s_max < args.s_min:
        parser.error("require s-min <= s-max")

    for s in range(args.s_min, args.s_max + 1):
        for n in range(args.n_min, args.n_max + 1):
            print(check(n, s))


if __name__ == "__main__":
    main()
