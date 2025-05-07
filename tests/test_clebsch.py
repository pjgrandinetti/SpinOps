import unittest

from sympy import N, S
from sympy.physics.wigner import clebsch_gordan as sympy_cg

from spinOps import \
    clebsch  # Uses the Cython wrapper that internally converts 2*j


def cg_ref(j1, m1, j2, m2, j, m):
    # Call sympy_cg with (j1, j2, j, m1, m2, m) ordering
    return float(N(sympy_cg(S(j1), S(j2), S(j), S(m1), S(m2), S(m)), 15))


class TestClebsch(unittest.TestCase):
    def test_valid_clebsch(self):
        result = clebsch(j1=1, m1=1, j2=1, m2=-1, j=1, m=0)
        expected = cg_ref(1, 1, 1, -1, 1, 0)
        self.assertAlmostEqual(abs(result), abs(expected), places=6)

        result = clebsch(j1=1, m1=0, j2=1, m2=0, j=1, m=0)
        expected = cg_ref(1, 0, 1, 0, 1, 0)
        self.assertAlmostEqual(abs(result), abs(expected), places=6)

    def test_edge_cases(self):
        result = clebsch(j1=0.5, m1=0.5, j2=0.5, m2=-0.5, j=0, m=0)
        expected = cg_ref(0.5, 0.5, 0.5, -0.5, 0, 0)
        self.assertAlmostEqual(abs(result), abs(expected), places=6)

        result = clebsch(j1=0.5, m1=-0.5, j2=0.5, m2=0.5, j=0, m=0)
        expected = cg_ref(0.5, -0.5, 0.5, 0.5, 0, 0)
        self.assertAlmostEqual(abs(result), abs(expected), places=6)

    def test_selection_rule_violation(self):
        """
        These combinations violate CG selection rules and must return zero.
        """
        self.assertAlmostEqual(
            clebsch(j1=1, m1=1, j2=1, m2=-1, j=1, m=1), 0.0, places=10
        )  # m1 + m2 ≠ m
        self.assertAlmostEqual(
            clebsch(j1=1, m1=1, j2=1, m2=-1, j=3, m=0), 0.0, places=10
        )  # |j1 - j2| > j

    def test_symmetry_identity(self):
        j1, m1 = 1.0, 0.0
        j2, m2 = 0.5, -0.5
        j, m = 1.5, -0.5
        c1 = clebsch(j1=j1, m1=m1, j2=j2, m2=m2, j=j, m=m)
        c2 = clebsch(j1=j2, m1=m2, j2=j1, m2=m1, j=j, m=m)
        phase = (-1) ** int(2 * (j1 + j2 - j))
        self.assertAlmostEqual(c1, phase * c2, places=6)

    def test_consistency_with_sympy(self):
        j1 = j2 = 0.5
        for m1 in [-0.5, 0.5]:
            for m2 in [-0.5, 0.5]:
                for j in [0.0, 1.0]:
                    m = m1 + m2
                    if abs(m) <= j:
                        result = clebsch(j1=j1, m1=m1, j2=j2, m2=m2, j=j, m=m)
                        expected = cg_ref(j1, m1, j2, m2, j, m)
                        self.assertAlmostEqual(abs(result), abs(expected), places=6)
