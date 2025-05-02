import unittest
import numpy as np
from sympy.physics.wigner import clebsch_gordan as sympy_cg
from sympy.physics.wigner import wigner_d_small as sympy_wigner_d
from spinOps import create_single_spin_Ix, clebsch, wigner_d


class TestClebsch(unittest.TestCase):
    def test_valid_clebsch(self):
        # Test valid Clebsch-Gordan coefficients
        result = clebsch(j1=1, m1=1, j2=1, m2=-1, j=1, m=0)
        expected = float(sympy_cg(1, 1, 1, 1, -1, 0))  # Sympy reference
        self.assertAlmostEqual(result, expected, places=6)

        result = clebsch(j1=1, m1=0, j2=1, m2=0, j=1, m=0)
        expected = float(sympy_cg(1, 1, 1, 0, 0, 0))  # Sympy reference
        self.assertAlmostEqual(result, expected, places=6)

    def test_invalid_magnetic_quantum_numbers(self):
        # Test invalid magnetic quantum numbers (m1 + m2 != m)
        with self.assertRaises(ValueError):
            clebsch(1, 1, 1, -1, 1, 1)

    def test_invalid_total_angular_momentum(self):
        # Test invalid total angular momentum (|j1 - j2| > j or j > j1 + j2)
        with self.assertRaises(ValueError):
            clebsch(j1=1, m1=1, j2=1, m2=-1, j=3, m=0)

    def test_edge_cases(self):
        # Test edge cases
        result = clebsch(j1=0.5, m1=0.5, j2=0.5, m2=-0.5, j=1, m=0)
        expected = float(sympy_cg(0.5, 0.5, 1, 0.5, -0.5, 0))  # Sympy reference
        self.assertAlmostEqual(result, expected, places=6)

        result = clebsch(j1=0.5, m1=-0.5, j2=0.5, m2=0.5, j=1, m=0)
        expected = float(sympy_cg(0.5, 0.5, 1, -0.5, 0.5, 0))  # Sympy reference
        self.assertAlmostEqual(result, expected, places=6)



class TestWignerD(unittest.TestCase):
    # Test valid reduced Wigner d-matrix elements
    def test_identity_at_zero_angle(self):
        # d^l_{m1,m2}(0) should equal delta_{m1,m2}
        for l in [0, 1, 2, 3]:
            for m1 in range(-l, l+1):
                for m2 in range(-l, l+1):
                    result = wigner_d(l=l, m1=m1, m2=m2, beta=0.0)
                    expected = 1.0 if m1 == m2 else 0.0
                    self.assertAlmostEqual(result, expected, places=6)

    def test_against_sympy(self):
        # Compare to sympy.physics.wigner.wigner_d for random beta
        import math
        angles = [math.pi/6, math.pi/4, math.pi/3]
        for l in [1, 2]:
            for m1 in range(-l, l+1):
                for m2 in range(-l, l+1):
                    for beta in angles:
                        result = wigner_d(l=l, m1=m1, m2=m2, beta=beta)
                        reference = float(sympy_wigner_d(l, beta)[l + m1, l + m2].evalf())
                        self.assertTrue(
                            math.isclose(result, reference, rel_tol=1e-6, abs_tol=1e-8),
                            f"Mismatch for l={l}, m1={m1}, m2={m2}, beta={beta}: {result} vs {reference}"
                        )


class TestSpinOps(unittest.TestCase):
    def test_create_single_spin_Ix(self):
        # Example input
        i_times_2 = [2, 2]  # Spin-1 system
        result = create_single_spin_Ix(0, i_times_2)

        # Expected output (manually computed or verified)
        expected_shape = (9, 9)  # For a spin-1 system
        self.assertEqual(result.shape, expected_shape)
        self.assertTrue(np.allclose(result, result.T.conj()))  # Hermitian check

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            create_single_spin_Ix(0, [])  # Empty i_times_2 list should raise an error



if __name__ == "__main__":
    unittest.main()
