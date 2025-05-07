import unittest

import numpy as np

from spinOps import (create_Imf, create_Ipf, create_Ixf, create_Iyf,
                     create_Izf, createEf, number_of_states)


class TestFictitiousSpinHalfOperators(unittest.TestCase):
    def setUp(self):
        # A few test systems: single spin-½, two spin-½’s, and a spin-1
        self.systems = [
            [1],  # 2 states
            [1, 1],  # 4 states
            [2],  # 3 states
        ]

    def test_invalid_inputs(self):
        for fn in (
            createEf,
            create_Ixf,
            create_Iyf,
            create_Izf,
            create_Ipf,
            create_Imf,
        ):
            with self.subTest(fn=fn.__name__):
                # empty list → ValueError
                with self.assertRaises(ValueError):
                    fn(0, 0, [])
                # negative r or s → IndexError
                with self.assertRaises(IndexError):
                    fn(-1, 0, [1])
                with self.assertRaises(IndexError):
                    fn(0, -1, [1])

    def test_createEf(self):
        """Ef sets diag[r,r]=1 and diag[s,s]=1."""
        for spins in self.systems:
            n = number_of_states(spins)
            for r in range(n):
                for s in range(n):
                    M = createEf(r, s, spins)
                    E = np.zeros((n, n), dtype=complex)
                    E[r, r] = 1
                    E[s, s] = 1
                    np.testing.assert_array_equal(
                        M, E, err_msg=f"createEf(r={r}, s={s}, spins={spins})"
                    )

    def test_create_Ixf(self):
        """Ixf sets (r,s) and (s,r) to 0.5."""
        for spins in self.systems:
            n = number_of_states(spins)
            for r in range(n):
                for s in range(n):
                    M = create_Ixf(r, s, spins)
                    E = np.zeros((n, n), dtype=complex)
                    E[r, s] = 0.5
                    E[s, r] = 0.5
                    np.testing.assert_array_equal(
                        M, E, err_msg=f"create_Ixf(r={r}, s={s}, spins={spins})"
                    )

    def test_create_Iyf(self):
        """Iyf sets (r,s)=+0.5j and (s,r)=-0.5j (but if r==s only +0.5j)."""
        for spins in self.systems:
            n = number_of_states(spins)
            for r in range(n):
                for s in range(n):
                    M = create_Iyf(r, s, spins)
                    E = np.zeros((n, n), dtype=complex)
                    if r == s:
                        E[r, s] = 0.5j
                    else:
                        E[r, s] = 0.5j
                        E[s, r] = -0.5j
                    np.testing.assert_allclose(
                        M, E, atol=0, err_msg=f"create_Iyf(r={r}, s={s}, spins={spins})"
                    )

    def test_create_Izf(self):
        """Izf sets diag[s,s]=+0.5; diag[r,r]=-0.5 (but if r==s only +0.5)."""
        for spins in self.systems:
            n = number_of_states(spins)
            for r in range(n):
                for s in range(n):
                    M = create_Izf(r, s, spins)
                    E = np.zeros((n, n), dtype=complex)
                    # always set the s,s element to +0.5
                    E[s, s] = 0.5
                    # if r!=s, also set r,r to -0.5
                    if r != s:
                        E[r, r] = -0.5
                    np.testing.assert_allclose(
                        M, E, atol=0, err_msg=f"create_Izf(r={r}, s={s}, spins={spins})"
                    )

    def test_create_Ipf_Imf(self):
        """
        Ipf sets (s,r)=1; Imf sets (r,s)=1.
        All other entries zero.
        """
        for spins in self.systems:
            n = number_of_states(spins)
            for r in range(n):
                for s in range(n):
                    Mpf = create_Ipf(r, s, spins)
                    Epf = np.zeros((n, n), dtype=complex)
                    Epf[s, r] = 1
                    np.testing.assert_array_equal(
                        Mpf, Epf, err_msg=f"create_Ipf(r={r}, s={s}, spins={spins})"
                    )

                    Mmf = create_Imf(r, s, spins)
                    Emf = np.zeros((n, n), dtype=complex)
                    Emf[r, s] = 1
                    np.testing.assert_array_equal(
                        Mmf, Emf, err_msg=f"create_Imf(r={r}, s={s}, spins={spins})"
                    )
