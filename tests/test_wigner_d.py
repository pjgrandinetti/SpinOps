import pytest
import math
from sympy.physics.wigner import wigner_d_small as sympy_wigner_d
from spinOps import wigner_d

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
