"""Tests for the standalone general Alexander-polynomial function."""

from __future__ import annotations

from math import gcd
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from alexander_polynomial import (  # noqa: E402
    alexander_coefficients,
    alexander_polynomial,
    normalize_parameters,
)
from check_fiberedness_end_families import (  # noqa: E402
    alexander_from_formula as reference_formula,
    euler_characteristic,
    full_twist_word,
)
from spherogram import Link  # noqa: E402


class GeneralAlexanderPolynomialTests(unittest.TestCase):
    def test_normalization(self):
        self.assertEqual(normalize_parameters(3, 5, 4, 2), (5, 3, 4, 2))
        self.assertEqual(normalize_parameters(-3, -5, 4, 2), (5, 3, 4, 2))
        self.assertEqual(normalize_parameters(3, -5, 4, 2), (5, 3, 4, -2))
        self.assertEqual(normalize_parameters(-3, 5, 4, 2), (5, 3, 4, -2))

    def test_elementary_symmetries(self):
        for p in range(3, 9):
            for q in range(1, p):
                if gcd(p, q) != 1:
                    continue
                for r in range(2, p + q):
                    for s in (-2, -1, 0, 1, 2):
                        expected = alexander_coefficients(p, q, r, s)
                        self.assertEqual(expected, alexander_coefficients(q, p, r, s))
                        self.assertEqual(expected, alexander_coefficients(-p, -q, r, s))
                        self.assertEqual(expected, alexander_coefficients(p, -q, r, -s))
                        self.assertEqual(expected, alexander_coefficients(-p, q, r, -s))

    def test_matches_manuscript_implementation(self):
        for p in range(3, 12):
            for q in range(1, p):
                if gcd(p, q) != 1:
                    continue
                for r in range(2, p + q):
                    for s in (-3, -1, 0, 1, 3):
                        self.assertEqual(
                            alexander_coefficients(p, q, r, s),
                            reference_formula(p, q, r, s),
                            (p, q, r, s),
                        )

    def test_zero_twist_is_independent_of_r(self):
        for p, q in ((5, 2), (7, 3), (8, 5), (9, 1)):
            values = {
                alexander_coefficients(p, q, r, 0)
                for r in range(2, p + q)
            }
            self.assertEqual(len(values), 1)

    def test_q_equal_one_is_covered(self):
        for p in range(2, 12):
            for r in range(2, p + 1):
                for s in (-3, -1, 0, 1, 3):
                    coefficients = alexander_coefficients(p, 1, r, s)
                    self.assertEqual(coefficients, coefficients[::-1])

    def test_q_equal_one_against_hfk_euler_characteristic(self):
        for p in range(3, 7):
            torus_unknot_braid = list(range(1, p))
            for r in range(2, p + 1):
                for s in (-1, 1):
                    link = Link(
                        braid_closure=torus_unknot_braid + full_twist_word(r, s)
                    )
                    link.simplify("pickup")
                    actual = (
                        euler_characteristic(link.knot_floer_homology())
                        if link.crossings
                        else (1,)
                    )
                    self.assertEqual(
                        actual,
                        alexander_coefficients(p, 1, r, s),
                        (p, 1, r, s),
                    )

    def test_readable_output(self):
        self.assertEqual(alexander_polynomial(3, 2, 2, 0), "t^2 - t + 1")
        self.assertEqual(
            alexander_polynomial(-3, -2, 2, 0),
            alexander_polynomial(3, 2, 2, 0),
        )

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            alexander_coefficients(0, 3, 2, 1)
        with self.assertRaises(ValueError):
            alexander_coefficients(6, 4, 3, 1)
        with self.assertRaises(ValueError):
            alexander_coefficients(5, 2, 0, 1)
        with self.assertRaises(ValueError):
            alexander_coefficients(5, 2, 8, 1)
        with self.assertRaises(TypeError):
            alexander_coefficients(5, 2, 3, 1.0)


if __name__ == "__main__":
    unittest.main()
