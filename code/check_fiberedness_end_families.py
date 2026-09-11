"""Knot Floer checks for two rigorously braided subfamilies with r > p > q.

The two end cases are r = p + 1 and r = p + q - 1.  In the first case we
positively stabilize the standard p-strand torus braid once.  In the second
case we positively destabilize once from the (p+q)-strand Lorenz shuffle
braid.  Appending s full twists on all r strands gives the required twisted
torus knot diagram in these two cases.

For every example, the Euler characteristic of HFK is checked against the
Alexander polynomial formula in the manuscript before the reported
fiberedness value is accepted.
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import gcd

from spherogram import Link


def _counter(exponents):
    return Counter(exponents)


def _convolve(left, right):
    result = Counter()
    for left_exp, left_coeff in left.items():
        for right_exp, right_coeff in right.items():
            result[left_exp + right_exp] += left_coeff * right_coeff
    return result


def alexander_from_formula(p: int, q: int, r: int, s: int):
    """Return the normalized coefficient tuple from the main theorem."""
    j_star = (-r * pow(p, -1, q)) % q
    s_values = [(-j * p) % q for j in range(q)]
    lengths_l = [(r - 1 - value) // q + 1 for value in s_values]

    # The congruence alone represents i_*=p by 0 at the boundary r=q.
    # Keep the integer representative used in the unified manuscript.
    numerator_i = p * j_star + r
    if numerator_i % q:
        raise ArithmeticError("the unified index i_star is not integral")
    i_star = p + 1 - numerator_i // q
    t_values = [(-i * q) % p for i in range(p)]
    lengths_m = [(r - 1 - value) // p + 1 for value in t_values]

    cumulative_l = []
    total = 0
    for value in lengths_l:
        total += value
        cumulative_l.append(total)

    cumulative_m = []
    total = 0
    for value in lengths_m:
        total += value
        cumulative_m.append(total)

    lambda_star = (lengths_l[j_star] - 1) * s + 1
    e_a = (
        r * s * (cumulative_l[j_star - 1] if j_star else 0)
        + p * j_star
        + r * lambda_star
    )
    e_c = r * r * s + p * q

    x_a = _counter(
        r * s * cumulative_l[m] + p * m for m in range(j_star)
    )
    x_c = _counter(
        r * s * cumulative_l[m] + p * m for m in range(q)
    )
    z_a = _counter(
        e_a
        - q * h
        - r * s * sum(lengths_m[p - h : p])
        for h in range(p - i_star + 1)
    )
    z_c = _counter(
        e_c - r * s * cumulative_m[m] - q * m for m in range(p)
    )

    numerator = _convolve(z_a, x_c)
    for exponent, coefficient in _convolve(x_a, z_c).items():
        numerator[exponent] -= coefficient
    numerator = Counter(
        {exponent: coefficient for exponent, coefficient in numerator.items() if coefficient}
    )

    # Divide by 1+t+...+t^(r-1).  Laurent shifts do not affect this division.
    quotient = Counter()
    remainder = numerator.copy()
    while remainder:
        top = max(remainder)
        coefficient = remainder[top]
        shift = top - r + 1
        quotient[shift] += coefficient
        for offset in range(r):
            exponent = shift + offset
            remainder[exponent] -= coefficient
            if remainder[exponent] == 0:
                del remainder[exponent]

    lowest = min(quotient)
    coefficients = tuple(
        quotient[lowest + offset]
        for offset in range(max(quotient) - lowest + 1)
    )
    if coefficients[-1] < 0:
        coefficients = tuple(-coefficient for coefficient in coefficients)
    return coefficients


def full_twist_word(r: int, s: int):
    if s > 0:
        return list(range(1, r)) * (r * s)
    if s < 0:
        return list(range(-1, -r, -1)) * (-r * s)
    return []


def lower_end_braid(p: int, q: int):
    """Base braid for r=p+1: one positive stabilization."""
    return list(range(1, p)) * q + [p]


def upper_end_braid(p: int, q: int):
    """Base braid for r=p+q-1: one Lorenz-braid destabilization."""
    n = p + q
    shuffle = []
    for j in range(q):
        shuffle.extend(range(p + j, j, -1))
    highest_generator = n - 1
    position = shuffle.index(highest_generator)
    # Cyclically move the unique highest generator to the end and destabilize.
    return shuffle[position + 1 :] + shuffle[:position]


def euler_characteristic(hfk_data):
    coefficients = Counter()
    for (alexander, maslov), rank in hfk_data["ranks"].items():
        coefficients[alexander] += (-1) ** maslov * rank
    lowest = min(coefficients)
    result = tuple(
        round(coefficients[lowest + offset])
        for offset in range(max(coefficients) - lowest + 1)
    )
    if result[-1] < 0:
        result = tuple(-coefficient for coefficient in result)
    return result


def check_example(p: int, q: int, endpoint: str, s: int):
    if endpoint == "lower":
        r = p + 1
        base = lower_end_braid(p, q)
    elif endpoint == "upper":
        r = p + q - 1
        base = upper_end_braid(p, q)
    else:
        raise ValueError("endpoint must be 'lower' or 'upper'")

    link = Link(braid_closure=base + full_twist_word(r, s))
    hfk = link.knot_floer_homology()
    formula = alexander_from_formula(p, q, r, s)
    euler = euler_characteristic(hfk)
    if euler != formula:
        raise AssertionError(
            f"Alexander mismatch for {(p, q, r, s)}: {euler} != {formula}"
        )
    return {
        "parameters": (p, q, r, s),
        "endpoint": endpoint,
        "crossings": len(link.crossings),
        "fibered": hfk["fibered"],
        "seifert_genus": hfk["seifert_genus"],
        "top_hfk_rank": sum(
            rank
            for (alexander, _), rank in hfk["ranks"].items()
            if alexander == hfk["seifert_genus"]
        ),
        "total_hfk_rank": hfk["total_rank"],
    }


def main():
    parser = argparse.ArgumentParser()
    # p=6 includes the substantially slower upper case (6,5,10,-2).
    parser.add_argument("--p-max", type=int, default=5)
    parser.add_argument(
        "--endpoints",
        choices=("lower", "upper"),
        nargs="+",
        default=("lower", "upper"),
    )
    parser.add_argument("--twists", type=int, nargs="+", default=[-2, -1, 1])
    args = parser.parse_args()

    seen = set()
    for p in range(3, args.p_max + 1):
        for q in range(2, p):
            if gcd(p, q) != 1:
                continue
            for endpoint in args.endpoints:
                for s in args.twists:
                    r = p + 1 if endpoint == "lower" else p + q - 1
                    parameters = (p, q, r, s)
                    if parameters in seen:
                        continue
                    seen.add(parameters)
                    print(check_example(p, q, endpoint, s), flush=True)


if __name__ == "__main__":
    main()
