# flake8: noqa: E741
import unittest

import numpy as np

from spinOps import (create_single_spin_Im, create_single_spin_Ip,
                     create_single_spin_Ix, create_single_spin_Iy,
                     create_single_spin_Iz, create_single_spin_Tlm,
                     create_single_spin_Tlm_unit, number_of_states)


def spin_basis(i_times_2):
    """
    Return for each spin the tuple (Ix, Iy, Iz, Ip, Im, E),
    constructed via the usual ladder formulas.
    """
    ops = []
    for it2 in i_times_2:
        dim = it2 + 1
        I = it2 / 2.0
        m = np.array([-I + k for k in range(dim)], dtype=float)

        # ladder operators
        I_plus = np.zeros((dim, dim), dtype=complex)
        I_minus = np.zeros((dim, dim), dtype=complex)
        for k in range(dim - 1):
            coeff = np.sqrt(I * (I + 1) - m[k] * (m[k] + 1))
            I_plus[k + 1, k] = coeff
        for k in range(1, dim):
            coeff = np.sqrt(I * (I + 1) - m[k] * (m[k] - 1))
            I_minus[k - 1, k] = coeff

        Ix = 0.5 * (I_plus + I_minus)
        Iy = (I_plus - I_minus) / (2j)
        Iz = np.diag(m)
        Ip = I_plus
        Im = I_minus
        E = np.eye(dim, dtype=complex)

        ops.append((Ix, Iy, Iz, Ip, Im, E))
    return ops


def expected_single_spin(op_name, spin_index, i_times_2):
    """
    Build the full-system operator named op_name via Kronecker products.
    op_name in {"Ix","Iy","Iz","Ip","Im"}
    """
    basis = spin_basis(i_times_2)
    mats = []
    for idx, (Ix, Iy, Iz, Ip, Im, E) in enumerate(basis):
        mats.append(
            {"Ix": Ix, "Iy": Iy, "Iz": Iz, "Ip": Ip, "Im": Im}[op_name]
            if idx == spin_index
            else E
        )
    result = mats[0]
    for M in mats[1:]:
        result = np.kron(result, M)
    return result


class TestSpinSystem(unittest.TestCase):

    def setUp(self):
        # sample spin systems
        self.systems = [
            [1],  # spin-1/2
            [2],  # spin-1
            [1, 1],  # two spin-1/2
            [3, 1],  # spin-3/2 + spin-1/2
        ]

    def test_number_of_states_invalid(self):
        with self.assertRaises(ValueError):
            number_of_states([])

    def test_invalid_spin_index(self):
        # Ix, Iy, Iz, Ip, Im
        simple_fns = [
            create_single_spin_Ix,
            create_single_spin_Iy,
            create_single_spin_Iz,
            create_single_spin_Ip,
            create_single_spin_Im,
        ]
        for fn in simple_fns:
            with self.subTest(fn=fn.__name__):
                with self.assertRaises(ValueError):
                    fn(0, [])
                with self.assertRaises(IndexError):
                    fn(-1, [1])
                with self.assertRaises(IndexError):
                    fn(1, [1])

        # Tlm constructors: signature (L, M, spin_index, i_times_2)
        tensor_fns = [create_single_spin_Tlm, create_single_spin_Tlm_unit]
        for fn in tensor_fns:
            with self.subTest(fn=fn.__name__):
                with self.assertRaises(ValueError):
                    fn(0, 0, 0, [])
                with self.assertRaises(IndexError):
                    fn(0, 0, -1, [1])
                with self.assertRaises(IndexError):
                    fn(0, 0, 1, [1])

    def test_Ix_Iy_Iz_Hermitian(self):
        """I_x, I_y, I_z must be Hermitian."""
        for spins in self.systems:
            n = number_of_states(spins)
            for op_name, ctor in [
                ("Ix", create_single_spin_Ix),
                ("Iy", create_single_spin_Iy),
                ("Iz", create_single_spin_Iz),
            ]:
                for idx in range(len(spins)):
                    M = ctor(idx, spins)
                    self.assertEqual(M.shape, (n, n))
                    np.testing.assert_allclose(
                        M,
                        M.conj().T,
                        atol=1e-12,
                        err_msg=f"{op_name} not Hermitian for spins={spins}, idx={idx}",
                    )

    def test_Ip_Im_adjoint(self):
        """I_+ = I_x + i I_y, I_- = I_x - i I_y, and I_+† = I_-"""
        for spins in self.systems:
            for idx in range(len(spins)):
                Ix = create_single_spin_Ix(idx, spins)
                Iy = create_single_spin_Iy(idx, spins)
                Ip = create_single_spin_Ip(idx, spins)
                Im = create_single_spin_Im(idx, spins)
                np.testing.assert_allclose(Ip, Ix + 1j * Iy, atol=1e-12)
                np.testing.assert_allclose(Im, Ix - 1j * Iy, atol=1e-12)
                np.testing.assert_allclose(Ip.conj().T, Im, atol=1e-12)

    def test_against_explicit_construction(self):
        """Compare ladder operators against explicit mixed-base helper."""
        for spins in self.systems:
            n = number_of_states(spins)
            for op_name, ctor in [
                ("Ix", create_single_spin_Ix),
                ("Iy", create_single_spin_Iy),
                ("Iz", create_single_spin_Iz),
                ("Ip", create_single_spin_Ip),
                ("Im", create_single_spin_Im),
            ]:
                for idx in range(len(spins)):
                    M = ctor(idx, spins)
                    exp = expected_single_spin(op_name, idx, spins)
                    self.assertEqual(M.shape, (n, n))
                    np.testing.assert_allclose(
                        M,
                        exp,
                        atol=1e-12,
                        err_msg=f"{op_name} mismatch for spins={spins}, idx={idx}",
                    )

    def test_Tlm_00_identity(self):
        """T_{0,0} = I, unit T_{0,0} = 1/√(2I+1)·I."""
        for spins in self.systems:
            for idx in range(len(spins)):
                twoI = spins[idx]
                n = number_of_states(spins)
                # raw tensor should be identity
                M_raw = create_single_spin_Tlm(0, 0, idx, spins)
                np.testing.assert_allclose(M_raw, np.eye(n, dtype=complex), atol=1e-12)
                # unit tensor scaled by 1/√(2I+1)
                M_unit = create_single_spin_Tlm_unit(0, 0, idx, spins)
                factor = 1 / np.sqrt(twoI + 1)
                np.testing.assert_allclose(
                    M_unit, factor * np.eye(n, dtype=complex), atol=1e-12
                )

    def test_Tlm_out_of_range_zero(self):
        """If L>2I or |M|>L, raw and unit tensors are zero."""
        for spins in self.systems:
            for idx in range(len(spins)):
                twoI = spins[idx]
                n = number_of_states(spins)
                # L > 2I
                Z_raw = create_single_spin_Tlm(twoI + 1, 0, idx, spins)
                Z_unit = create_single_spin_Tlm_unit(twoI + 1, 0, idx, spins)
                np.testing.assert_allclose(
                    Z_raw, np.zeros((n, n), dtype=complex), atol=1e-12
                )
                np.testing.assert_allclose(
                    Z_unit, np.zeros((n, n), dtype=complex), atol=1e-12
                )
                # |M| > L
                Z2_raw = create_single_spin_Tlm(0, twoI + 1, idx, spins)
                Z2_unit = create_single_spin_Tlm_unit(0, twoI + 1, idx, spins)
                np.testing.assert_allclose(
                    Z2_raw, np.zeros((n, n), dtype=complex), atol=1e-12
                )
                np.testing.assert_allclose(
                    Z2_unit, np.zeros((n, n), dtype=complex), atol=1e-12
                )

    def test_Tlm_Hermitian_conjugation(self):
        """T_{L,M}^† = (-)^M T_{L,-M} for raw and unit tensor."""
        for spins in self.systems:
            for idx in range(len(spins)):
                twoI = spins[idx]
                for L in range(1, twoI + 1):
                    for M_ in range(-L, L + 1):
                        T = create_single_spin_Tlm(L, M_, idx, spins)
                        Tm = create_single_spin_Tlm(L, -M_, idx, spins)
                        sign = (-1) ** M_
                        np.testing.assert_allclose(T.conj().T, sign * Tm, atol=1e-12)
                        U = create_single_spin_Tlm_unit(L, M_, idx, spins)
                        Um = create_single_spin_Tlm_unit(L, -M_, idx, spins)
                        np.testing.assert_allclose(U.conj().T, sign * Um, atol=1e-12)
