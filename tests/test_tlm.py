import math
import pytest
from sympy import Rational
from sympy.physics.wigner import wigner_3j
from spinOps import tlm, unit_tlm


class TestUnitTlm:
    @staticmethod
    def tlm_ref(l: int, m: int, I: float, m1: float, m2: float) -> float:
        """
        Reference implementation of <I,m1| T_{l,m} |I,m2> using the Wigner-Eckart theorem (SI § S2),
        via SymPy's wigner_3j.
        """
        # Reduced matrix element ⟨I‖T^(l)‖I⟩ from Eq. S.62
        num = math.factorial(l) * math.factorial(l) * math.factorial(int(2*I + l + 1))
        den = (2**l) * math.factorial(2*l) * math.factorial(int(2*I - l))
        red = math.sqrt(num / den)

        # Prepare rational spins for wigner_3j
        j   = Rational(int(2*I), 2)
        mj1 = Rational(-int(2*m1), 2)
        mj  = Rational(int(2*m), 2)
        mj2 = Rational(int(2*m2), 2)

        threej = wigner_3j(j, l, j, mj1, mj, mj2)
        phase  = (-1) ** int(I - m1)

        return float(phase * threej * red)

    @pytest.mark.parametrize(
        "l, m, I, m1, m2", [
            # integer-spin cases
            (1,  0, 1.0,  1.0,  0.0),
            (2, -1, 3.0, -2.0, -1.0),
            (3,  2, 2.5,  1.0,  3.0),
            # half-integer spin example
            (1,  1, 0.5, -0.5,  0.5),
        ]
    )
    def test_tlm_against_reference(self, l, m, I, m1, m2):
        """
        Compare the Cython tlm() against the SymPy-based reference implementation.
        """
        got = tlm(l, m, I, m1, m2)
        exp = self.tlm_ref(l, m, I, m1, m2)
        assert got == pytest.approx(exp, rel=1e-9), (
            f"tlm({l},{m},{I},{m1},{m2}) = {got} ≠ {exp}"
        )

    @pytest.mark.parametrize(
        "l, m, I, m1, m2", [
            # reuse the same “happy‐path” examples
            (1,  0, 1.0,  1.0,  0.0),
            (2, -1, 3.0, -2.0, -1.0),
            (3,  2, 2.5,  1.0,  3.0),
            (1,  1, 0.5, -0.5,  0.5),
        ]
    )
    def test_unit_tlm_against_reference(self, l, m, I, m1, m2):
        """
        Compare unit_tlm() against the analytic unit‐tensor scaling of tlm_ref.
        """
        # 1) get the raw T_{l,m} reference
        raw = self.tlm_ref(l, m, I, m1, m2)

        # 2) compute the unit‐tensor scaling:
        #    scale = 1/l! * sqrt[(2l+1)(2I-l)! 2^l (2l)! / (2I+l+1)!]
        inv_l_fact = 1.0 / math.factorial(l)
        num = (2*l + 1) \
            * math.factorial(int(2*I - l)) \
            * (2**l) \
            * math.factorial(2*l)
        den = math.factorial(int(2*I + l + 1))
        scale = inv_l_fact * math.sqrt(num / den)

        expected = raw * scale
        got = unit_tlm(l, m, I, m1, m2)

        assert got == pytest.approx(expected, rel=1e-9), (
            f"unit_tlm({l},{m},{I},{m1},{m2}) = {got} ≠ {expected}"
        )