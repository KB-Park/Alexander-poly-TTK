#!/usr/bin/env python3
"""Compute the Alexander polynomial of a twisted torus knot.

The input is T(p,q;r,s), where p and q are nonzero relatively prime
integers, 1 < r < |p|+|q|, and s is any integer.  The output is normalized
up to the usual multiplication by +/-t^k: the lowest exponent is zero and
the leading coefficient is positive.

This module uses only the Python standard library.
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import gcd


def normalize_parameters(p: int, q: int, r: int, s: int) -> tuple[int, int, int, int]:
    """Reduce the parameters to P > Q > 0 using the paper's symmetries.

    Simultaneous sign reversal of p and q preserves the knot.  If p and q
    have opposite signs, mirroring makes them have the same sign and changes
    s to -s; the Alexander polynomial is unchanged by mirroring.  Exchanging
    p and q then gives P > Q > 0.
    """
    for name, value in (("p", p), ("q", q), ("r", r), ("s", s)):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer")
    if p == 0 or q == 0:
        raise ValueError("p and q must be nonzero")
    if gcd(abs(p), abs(q)) != 1:
        raise ValueError("p and q must be relatively prime")
    if not (1 < r < abs(p) + abs(q)):
        raise ValueError("require 1 < r < |p|+|q|")

    p_abs, q_abs = abs(p), abs(q)
    if p_abs == q_abs:
        # Relative primeness forces p_abs=q_abs=1, for which no admissible
        # integer r exists; keep the error explicit for direct callers.
        raise ValueError("require |p| and |q| to be distinct")
    P, Q = max(p_abs, q_abs), min(p_abs, q_abs)
    S = s if p * q > 0 else -s
    return P, Q, r, S


def _convolve(left: Counter[int], right: Counter[int]) -> Counter[int]:
    result: Counter[int] = Counter()
    for left_exp, left_coefficient in left.items():
        for right_exp, right_coefficient in right.items():
            result[left_exp + right_exp] += left_coefficient * right_coefficient
    return Counter({exponent: value for exponent, value in result.items() if value})


def _divide_by_geometric_sum(
    numerator: Counter[int], length: int
) -> tuple[int, ...]:
    """Divide exactly by 1+t+...+t^(length-1), allowing Laurent shifts."""
    if not numerator:
        raise ArithmeticError("Alexander numerator unexpectedly vanished")

    low, high = min(numerator), max(numerator)
    quotient: list[int] = []
    previous_window = 0
    for exponent in range(low, high + 1):
        coefficient = numerator.get(exponent, 0) - previous_window
        quotient.append(coefficient)
        previous_window += coefficient
        if len(quotient) >= length:
            previous_window -= quotient[-length]

    if previous_window != 0:
        raise ArithmeticError("nonzero remainder in Alexander division")
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    while len(quotient) > 1 and quotient[0] == 0:
        quotient.pop(0)
    if quotient[-1] < 0:
        quotient = [-coefficient for coefficient in quotient]
    return tuple(quotient)


def _torus_alexander_coefficients(p: int, q: int) -> tuple[int, ...]:
    """Normalized dense coefficients of Delta_{T(p,q)}(t)."""
    numerator: Counter[int] = Counter()
    for exponent, coefficient in (
        (0, 1),
        (1, -1),
        (p * q, -1),
        (p * q + 1, 1),
    ):
        numerator[exponent] += coefficient
        if numerator[exponent] == 0:
            del numerator[exponent]

    denominator = _convolve(
        Counter({0: 1, p: -1}), Counter({0: 1, q: -1})
    )
    degree = (p - 1) * (q - 1)
    remainder = [numerator.get(index, 0) for index in range(p * q + 2)]
    quotient: list[int] = []
    for exponent in range(degree + 1):
        coefficient = remainder[exponent]
        quotient.append(coefficient)
        for offset, value in denominator.items():
            remainder[exponent + offset] -= coefficient * value
    if any(remainder):
        raise ArithmeticError("torus Alexander polynomial division was not exact")
    if quotient[-1] < 0:
        quotient = [-coefficient for coefficient in quotient]
    return tuple(quotient)


def _positive_formula(p: int, q: int, r: int, s: int) -> tuple[int, ...]:
    """Apply the unified block formula for p > q > 0."""
    if s == 0:
        return _torus_alexander_coefficients(p, q)

    # Python defines the unique residue modulo 1 to be zero, so the same
    # formulas also cover q=1.
    j_star = (-r * pow(p, -1, q)) % q
    s_values = [(-j * p) % q for j in range(q)]
    lengths_l = [(r - 1 - value) // q + 1 for value in s_values]

    numerator_i = p * j_star + r
    if numerator_i % q:
        raise ArithmeticError("the index i_star is not integral")
    i_star = p + 1 - numerator_i // q
    if not (1 <= i_star <= p):
        raise ArithmeticError("the index i_star is outside 1,...,p")

    t_values = [(-i * q) % p for i in range(p)]
    lengths_m = [(r - 1 - value) // p + 1 for value in t_values]

    cumulative_l: list[int] = []
    total = 0
    for value in lengths_l:
        total += value
        cumulative_l.append(total)

    cumulative_m: list[int] = []
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

    x_a = Counter(
        r * s * cumulative_l[index] + p * index
        for index in range(j_star)
    )
    x_c = Counter(
        r * s * cumulative_l[index] + p * index
        for index in range(q)
    )

    suffix_m = [0] * (p + 1)
    for index in range(p - 1, -1, -1):
        suffix_m[index] = suffix_m[index + 1] + lengths_m[index]
    z_a = Counter(
        e_a - q * h - r * s * suffix_m[p - h]
        for h in range(p - i_star + 1)
    )
    z_c = Counter(
        e_c - r * s * cumulative_m[index] - q * index
        for index in range(p)
    )

    numerator = _convolve(z_a, x_c)
    for exponent, coefficient in _convolve(x_a, z_c).items():
        numerator[exponent] -= coefficient
        if numerator[exponent] == 0:
            del numerator[exponent]

    coefficients = _divide_by_geometric_sum(numerator, r)
    if coefficients != coefficients[::-1]:
        raise ArithmeticError("the computed Alexander polynomial is not symmetric")
    return coefficients


def alexander_coefficients(p: int, q: int, r: int, s: int) -> tuple[int, ...]:
    """Return normalized dense coefficients in ascending exponent order.

    If the return value is (a_0,...,a_d), it represents

        a_0 + a_1*t + ... + a_d*t^d.
    """
    P, Q, R, S = normalize_parameters(p, q, r, s)
    return _positive_formula(P, Q, R, S)


def format_polynomial(coefficients: tuple[int, ...], variable: str = "t") -> str:
    """Format a dense coefficient tuple as a descending polynomial."""
    if not variable or any(character.isspace() for character in variable):
        raise ValueError("variable must be a nonempty string without whitespace")

    terms: list[tuple[str, str]] = []
    for exponent in range(len(coefficients) - 1, -1, -1):
        coefficient = coefficients[exponent]
        if coefficient == 0:
            continue
        magnitude = abs(coefficient)
        if exponent == 0:
            body = str(magnitude)
        elif exponent == 1:
            body = variable if magnitude == 1 else f"{magnitude}*{variable}"
        else:
            power = f"{variable}^{exponent}"
            body = power if magnitude == 1 else f"{magnitude}*{power}"
        terms.append(("+" if coefficient > 0 else "-", body))

    if not terms:
        return "0"
    first_sign, first_body = terms[0]
    result = ("" if first_sign == "+" else "-") + first_body
    for sign, body in terms[1:]:
        result += f" {sign} {body}"
    return result


def alexander_polynomial(
    p: int, q: int, r: int, s: int, variable: str = "t"
) -> str:
    """Return the normalized Alexander polynomial as a readable string."""
    return format_polynomial(alexander_coefficients(p, q, r, s), variable)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p", type=int)
    parser.add_argument("q", type=int)
    parser.add_argument("r", type=int)
    parser.add_argument("s", type=int)
    parser.add_argument(
        "--coefficients",
        action="store_true",
        help="also print the dense coefficient tuple",
    )
    args = parser.parse_args()

    normalized = normalize_parameters(args.p, args.q, args.r, args.s)
    coefficients = alexander_coefficients(args.p, args.q, args.r, args.s)
    print(f"normalized parameters: {normalized}")
    print(f"Delta(t) = {format_polynomial(coefficients)}")
    if args.coefficients:
        print(f"coefficients (ascending): {coefficients}")


if __name__ == "__main__":
    main()
