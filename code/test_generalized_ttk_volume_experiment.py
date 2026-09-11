#!/usr/bin/env python3
"""Regression tests for generalized_ttk_volume_experiment.py."""

from __future__ import annotations

import argparse
from math import gcd
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from check_fiberedness_end_families import (  # noqa: E402
    alexander_from_formula as manuscript_reference,
    euler_characteristic,
)
from generalized_ttk_volume_experiment import (  # noqa: E402
    Parameters,
    alexander_coefficients,
    generalized_braid_word,
    known_torus_type,
    parameter_range,
)


class AlexanderTests(unittest.TestCase):
    def test_matches_existing_manuscript_implementation(self):
        for p in range(3, 9):
            for q in range(1, p):
                if gcd(p, q) != 1:
                    continue
                for r in range(2, p + q):
                    for s in (1, 2):
                        item = Parameters(p, q, r, s)
                        self.assertEqual(
                            alexander_coefficients(item),
                            manuscript_reference(p, q, r, s),
                            item,
                        )

    def test_r_equal_q_is_the_expected_torus_knot(self):
        self.assertEqual(
            alexander_coefficients(Parameters(3, 2, 2, 1)),
            (1, -1, 1, -1, 1),  # T(2,5)
        )

    def test_r_one_is_the_torus_knot(self):
        self.assertEqual(alexander_coefficients(Parameters(5, 1, 1, 3)), (1,))
        self.assertEqual(
            alexander_coefficients(Parameters(5, 2, 1, 3)),
            (1, -1, 1, -1, 1),
        )

    def test_previous_paper_range_has_7450_cases(self):
        args = argparse.Namespace(
            p_min=3,
            p_max=20,
            q_min=1,
            q_max=None,
            r_min=None,
            r_max=None,
            s_min=1,
            s_max=5,
            region="old",
        )
        self.assertEqual(sum(1 for _ in parameter_range(args)), 7450)

    def test_elementary_torus_boundaries(self):
        self.assertEqual(known_torus_type(Parameters(7, 4, 1, 2)), "T(7,4)")
        self.assertEqual(known_torus_type(Parameters(7, 4, 4, 2)), "T(4,15)")
        self.assertEqual(known_torus_type(Parameters(7, 4, 7, 2)), "T(7,18)")
        self.assertEqual(known_torus_type(Parameters(7, 4, 8, 2)), "")


class BraidTests(unittest.TestCase):
    def test_generalized_braid_is_a_knot_and_matches_formula(self):
        try:
            import snappy  # noqa: F401
            from spherogram import Link
        except ImportError as error:
            self.skipTest(str(error))

        item = Parameters(4, 3, 5, 1)
        word, strands = generalized_braid_word(item)
        self.assertEqual(strands, 7)
        link = Link(braid_closure=word)
        self.assertEqual(len(link.link_components), 1)
        link.simplify("basic")
        link.simplify("pickup")
        computed = euler_characteristic(link.knot_floer_homology())
        self.assertEqual(computed, alexander_coefficients(item))

    def test_negative_twist_braid_matches_formula(self):
        try:
            import snappy  # noqa: F401
            from spherogram import Link
        except ImportError as error:
            self.skipTest(str(error))

        item = Parameters(4, 3, 5, -1)
        word, strands = generalized_braid_word(item)
        self.assertEqual(strands, 7)
        self.assertIn(-1, word)
        link = Link(braid_closure=word)
        self.assertEqual(len(link.link_components), 1)
        link.simplify("basic")
        link.simplify("pickup")
        computed = euler_characteristic(link.knot_floer_homology())
        self.assertEqual(computed, alexander_coefficients(item))


if __name__ == "__main__":
    unittest.main()
