"""Endpoint regressions against independent braid/HFK and torus calculations."""
from math import gcd
import unittest
from alexander_polynomial import alexander_coefficients, _torus_alexander_coefficients
from generalized_ttk_volume_experiment import (
    Parameters, alexander_coefficients as experiment_coefficients,
    generalized_braid_word, known_torus_type, parameter_range, parse_args,
)

class EndpointTests(unittest.TestCase):
    def test_endpoint_algebra_and_symmetries(self):
        for p in range(2, 10):
            for q in range(1, p):
                if gcd(p, q) != 1:
                    continue
                for s in (-3, -1, 0, 1, 3):
                    r = p + q
                    result = alexander_coefficients(p, q, r, s)
                    self.assertEqual(result, result[::-1])
                    self.assertEqual(abs(sum(result)), 1)
                    self.assertEqual(result, experiment_coefficients(Parameters(p, q, r, s)))
                    self.assertEqual(result, alexander_coefficients(-q, p, r, -s))
                    if s == 0:
                        self.assertEqual(result, _torus_alexander_coefficients(p, q))

    def test_endpoint_q_one_and_trivial_p_q(self):
        for p in range(1, 8):
            for s in (-2, -1, 0, 1, 2):
                r = p + 1
                self.assertEqual(alexander_coefficients(p, 1, r, s),
                    _torus_alexander_coefficients(r, abs(r * s + 1)))

    def test_torus_boundaries_for_both_twist_signs(self):
        for s in (-2, -1, 0, 1, 2):
            for r, a, b in ((1, 7, 4), (4, 4, 7 + 4*s), (7, 7, 4 + 7*s)):
                expected = _torus_alexander_coefficients(a, abs(b))
                self.assertEqual(alexander_coefficients(7, 4, r, s), expected)
                self.assertEqual(experiment_coefficients(Parameters(7, 4, r, s)), expected)

    def test_scan_opt_in_and_signed_twists(self):
        options = ['--p-min', '3', '--p-max', '3', '--q-min', '2',
                   '--s-min', '-1', '--s-max', '1']
        historical = list(parameter_range(parse_args(options)))
        extended = list(parameter_range(parse_args(options + ['--include-endpoint'])))
        self.assertTrue(all(item.r < item.p + item.q for item in historical))
        self.assertEqual(len(extended) - len(historical), 3)
        self.assertIn(Parameters(3, 2, 5, -1), extended)
        self.assertEqual(known_torus_type(Parameters(3, 2, 5, 0)), 'T(3,2)')
        self.assertEqual(known_torus_type(Parameters(3, 2, 5, 3)), '')

    def test_endpoint_matches_braid_hfk(self):
        from spherogram import Link
        from check_fiberedness_end_families import euler_characteristic
        for p, q, s in ((2, 1, -1), (3, 2, -1), (3, 2, 1), (4, 3, 1)):
            item = Parameters(p, q, p + q, s)
            word, strands = generalized_braid_word(item)
            self.assertEqual(strands, p + q)
            link = Link(braid_closure=word)
            self.assertEqual(len(link.link_components), 1)
            link.simplify('pickup')
            self.assertEqual(euler_characteristic(link.knot_floer_homology()),
                             experiment_coefficients(item), item)

if __name__ == '__main__':
    unittest.main()
