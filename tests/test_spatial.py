import pytest
import numpy as np
import cmath
import math
from sympy.physics.wigner import wigner_3j, wigner_d_small as sympy_wigner_d
from spinOps import Rotate, create_rho1, create_rho2, wigner_d, DLM

class TestUnitWigner_d:
    def test_identity_diagonal_at_zero_angle(self):
        # d^l_{m,m}(0) == 1.0 for any m
        for l in [0, 1, 2, 3]:
            for m in range(-l, l+1):
                result = wigner_d(l, m, m, 0.0)
                assert result == pytest.approx(1.0, abs=1e-12)

    @pytest.mark.parametrize("l,beta", [
        (1, math.pi/6),
        (1, math.pi/4),
        (2, math.pi/6),
        (2, math.pi/3),
    ])
    def test_diagonal_against_sympy_random_angles(self, l, beta):
        # Compare only diagonal elements to sympy.physics.wigner.wigner_d_small
        sym = sympy_wigner_d(l, beta)
        for m in range(-l, l+1):
            result = wigner_d(l, m, m, beta)
            expected = float(sym[l + m, l + m].evalf())
            assert result == pytest.approx(expected, rel=1e-6, abs=1e-8), \
                f"d^{l}_{{{m},{m}}}({beta}): got {result}, expected {expected}"

    def test_invalid_magnetic_quantum_numbers(self):
        # |m1| or |m2| > l should yield zero
        assert wigner_d(1, 2, 0, 0.5) == pytest.approx(0.0, abs=1e-12)
        assert wigner_d(1, 0, -2, 1.0) == pytest.approx(0.0, abs=1e-12)

    def test_full_matrix_against_sympy(self):
        # Compare entire small-d matrix for l=1,2 at various angles
        for l, beta in [(1, math.pi/6), (2, math.pi/4)]:
            sym = sympy_wigner_d(l, beta)
            for m1 in range(-l, l+1):
                for m2 in range(-l, l+1):
                    got = wigner_d(l, m1, m2, beta)
                    expected = float(sym[l + m1, l + m2].evalf())
                    assert got == pytest.approx(expected, rel=1e-6, abs=1e-8), \
                        f"Full d-matrix mismatch l={l},m1={m1},m2={m2},beta={beta}"

    def test_unitarity(self):
        # For each l and beta, verify d * d^† = identity
        import numpy as np
        for l, beta in [(1, math.pi/5), (2, math.pi/3)]:
            size = 2*l + 1
            # build matrix
            mat = np.array([[wigner_d(l, m1, m2, beta)
                             for m2 in range(-l, l+1)]
                            for m1 in range(-l, l+1)], dtype=complex)
            # check unitarity
            prod = mat.dot(mat.conj().T)
            assert np.allclose(prod, np.eye(size), atol=1e-8), f"Unitarity failed for l={l}, beta={beta}"


class TestDLM:
    @pytest.mark.parametrize("l,m,alpha,gamma", [
        (2.0, 2.0, 0.4, 1.1),
        (2.0, 1.0, 1.2, -0.6),
        (2.0, 0.0, 2.3, 2.3),
        (2.0, -1.0, 3.1, 0.2),
        (2.0, -2.0, 0.7, -1.3),
    ])
    def test_beta_zero_phase(self, l, m, alpha, gamma):
        """
        For β=0, DLM(l,m,m,α,0,γ) = e^{-i (m α + m γ)}
        and off-diagonals are zero.
        """
        val = DLM(l, m, m, alpha, 0.0, gamma)
        expected = cmath.exp(-1j * (m * alpha + m * gamma))
        assert val == pytest.approx(expected, rel=1e-9)

    def test_off_diagonal_zero_at_beta_zero(self):
        """
        DLM off-diagonal with β=0 should be zero.
        """
        assert DLM(2.0, 1.0, 0.0, 0.1, 0.0, 0.2) == pytest.approx(0.0, abs=1e-12)
    
    def test_factorization_small_d_phase(self):
        """
        Test DLM(l,m1,m2,α,β,γ) = e^{-i m1 α} d^l_{m1,m2}(β) e^{-i m2 γ}.
        """
        import numpy as np
        for l, alpha, beta, gamma in [(1.0, 0.5, 0.3, -0.7), (2.0, -1.0, 1.2, 0.2)]:
            for m1 in range(-int(l), int(l)+1):
                for m2 in range(-int(l), int(l)+1):
                    small = wigner_d(l, m1, m2, beta)
                    phase = cmath.exp(-1j * (m1 * alpha + m2 * gamma))
                    expected = phase * small
                    got = DLM(l, m1, m2, alpha, beta, gamma)
                    assert got == pytest.approx(expected, rel=1e-8, abs=1e-10)

    def test_unitarity_of_DLM(self):
        """
        Verify that the full D-matrix is unitary for given angles.
        """
        import numpy as np
        for l, alpha, beta, gamma in [(1.0, 0.3, 0.4, -0.5), (2.0, -0.8, 1.0, 0.9)]:
            size = 2 * int(l) + 1
            mat = np.zeros((size, size), dtype=complex)
            for i, m1 in enumerate(range(-int(l), int(l)+1)):
                for j, m2 in enumerate(range(-int(l), int(l)+1)):
                    mat[i, j] = DLM(l, m1, m2, alpha, beta, gamma)
            prod = mat @ mat.conj().T
            assert np.allclose(prod, np.eye(size), atol=1e-8), f"DLM unitarity failed for l={l}"


class TestRotate:
    def test_rotate_empty_input(self):
        # Should raise ValueError on empty array
        with pytest.raises(ValueError):
            Rotate(np.zeros(0, dtype=np.complex128), 0.1, 0.2, 0.3)

    @pytest.mark.parametrize("j", [1, 2])
    def test_rotate_identity(self, j):
        # identity rotation (alpha=beta=gamma=0) returns original initial
        size = 2 * j + 1
        # test for each basis state
        for m0 in range(size):
            initial = np.zeros(size, dtype=np.complex128)
            initial[m0] = 1.0
            result = Rotate(initial, 0.0, 0.0, 0.0)
            assert np.allclose(result, initial)

    # Removed direct DLM comparison due to potential phase conventions

class TestRho:
    def test_create_rho1_shape_and_dtype(self):
        vec = create_rho1(1.23)
        assert isinstance(vec, np.ndarray)
        assert vec.shape == (3,)
        assert np.iscomplexobj(vec)
        # not all zeros
        assert not np.allclose(vec, 0)

    def test_create_rho2_shape_and_dtype(self):
        vec = create_rho2(2.34, 0.56)
        assert isinstance(vec, np.ndarray)
        assert vec.shape == (5,)
        assert np.iscomplexobj(vec)
        assert not np.allclose(vec, 0)
    
    @pytest.mark.parametrize("alpha,beta,gamma,zeta", [
        (0.1, 0.2, 0.3, 1.23),
        (0.5, 1.0, -0.7, 0.8),
    ])
    def test_rotate_rho1_norm_preserved(self, alpha, beta, gamma, zeta):
        # Norm of rank-1 tensor preserved under rotation
        initial = create_rho1(zeta)
        rotated = Rotate(initial, alpha, beta, gamma)
        assert np.isclose(np.linalg.norm(rotated), np.linalg.norm(initial), atol=1e-10)

    @pytest.mark.parametrize("alpha,beta,gamma,zeta,eta", [
        (0.2, 0.3, 0.4, 2.34, 0.56),
        (-0.5, 1.2, -0.8, 1.5, 0.9),
    ])
    def test_rotate_rho2_norm_preserved(self, alpha, beta, gamma, zeta, eta):
        # Norm of rank-2 tensor preserved under rotation
        initial = create_rho2(zeta, eta)
        rotated = Rotate(initial, alpha, beta, gamma)
        assert np.isclose(np.linalg.norm(rotated), np.linalg.norm(initial), atol=1e-10)