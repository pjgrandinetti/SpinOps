import unittest
import numpy as np
from spinOps import createIx
import unittest
from spinOps import clebsch

class TestClebsch(unittest.TestCase):
    def test_valid_clebsch(self):
        # Test valid Clebsch-Gordan coefficients
        result = clebsch(j1=1, m1=1, j2=1, m2=-1, j=1, m=0)
        expected = 1/np.sqrt(2)  # Known result for these quantum numbers
        self.assertAlmostEqual(result, expected, places=6)

        result = clebsch(j1=1, m1=0, j2=1, m2=0, j=1, m=0)
        expected = 0.0  # Known result for these quantum numbers
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
        expected = 1/np.sqrt(2)  # Known result for these quantum numbers
        self.assertAlmostEqual(result, expected, places=6)

        result = clebsch(j1=0.5, m1=-0.5, j2=0.5, m2=0.5, j=1, m=0)
        expected = 1/np.sqrt(2)  # Known result for these quantum numbers
        self.assertAlmostEqual(result, expected, places=6)



class TestSpinOps(unittest.TestCase):
    def test_createIx(self):
        # Example input
        i_times_2 = [2, 2]  # Spin-1 system
        result = createIx(0, i_times_2)

        # Expected output (manually computed or verified)
        expected_shape = (9, 9)  # For a spin-1 system
        self.assertEqual(result.shape, expected_shape)
        self.assertTrue(np.allclose(result, result.T.conj()))  # Hermitian check

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            createIx(0, [])  # Empty i_times_2 list should raise an error

if __name__ == "__main__":
    unittest.main()
    
