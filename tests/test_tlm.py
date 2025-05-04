import unittest
import numpy as np
from sympy.physics.wigner import clebsch_gordan as sympy_cg
from sympy.physics.wigner import wigner_d_small as sympy_wigner_d
from sympy.physics.wigner import wigner_3j
from sympy import S, factorial, sqrt  # make sure sqrt is included
import math
from spinOps import create_single_spin_Ix, clebsch, wigner_d
from spinOps import tlm, unit_tlm

def reduced_matrix_element(I, l):
    """Compute reduced matrix element ⟨I||T^l||I⟩ from Eq. S.63 of SI."""
    I = S(I)
    l = S(l)
    numerator = factorial(l) * factorial(l) * factorial(2 * I + l + 1)
    denominator = 2**l * factorial(2 * l) * factorial(2 * I - l)
    return float(sqrt(numerator / denominator))

class TestUnitTlm:

    @pytest.mark.parametrize("l, m, j, m1, m2", [
        (1,  1,   1,    1,   0),
        (2,  0,   1.5, -0.5, 0.5),
        (1,  0,   1,    1,  -1),
    ])
    def test_unit_tensor_scaling(self, l, m, j, m1, m2):
        # skip if selection rule not satisfied
        if abs(m1 + m2 - m) > 1e-12:
            pytest.skip("m₁ + m₂ != m → unit_tlm raises ValueError")
        nonunit = tlm(l, m, j, m1, j, m2)
        red     = reduced_matrix_element(j, l)            # Eq. S.62
        scale   = math.sqrt(2*l + 1) / (math.sqrt(2*j + 1) * red)
        expected = nonunit * scale
        actual   = unit_tlm(l, m, j, m1, j, m2)
        assert actual == pytest.approx(expected, rel=1e-12)

    def test_trace_orthonormality(self):
        """Tr[𝒯ₗₘ† 𝒯ₗₘ] = ∑ₘ₁,ₘ₂ |𝒯ₗₘ(m₁,m₂)|² = 1  for each (l,m)."""
        j   = 1
        tol = 1e-12

        for l in range(3):              # ℓ = 0,1,2
            for m in range(-l, l+1):
                acc = 0.0
                for m1 in [-1, 0, 1]:
                    for m2 in [-1, 0, 1]:
                        try:
                            t = unit_tlm(l, m, j, m1, j, m2)
                        except ValueError:
                            continue
                        acc += t*t
                assert acc == pytest.approx(1.0, abs=tol)

    def test_off_diagonal_orthogonality(self):
        """Different (l,m) channels are orthogonal:
           ∑ₘ₁,ₘ₂ 𝒯_{l₁,m₁}⋅𝒯_{l₂,m₂} = 0 for (l₁,m₁)≠(l₂,m₂)."""
        j   = 1
        tol = 1e-12

        for l1 in range(3):
            for m1 in range(-l1, l1+1):
                for l2 in range(3):
                    for m2 in range(-l2, l2+1):
                        if (l1, m1) == (l2, m2):
                            continue
                        acc = 0.0
                        for a in [-1, 0, 1]:
                            for b in [-1, 0, 1]:
                                try:
                                    t1 = unit_tlm(l1, m1, j, a, j, b)
                                    t2 = unit_tlm(l2, m2, j, a, j, b)
                                except ValueError:
                                    continue
                                acc += t1 * t2
                        assert acc == pytest.approx(0.0, abs=tol)
                        
