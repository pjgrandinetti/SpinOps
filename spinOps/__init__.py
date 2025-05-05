"""
SpinOps top-level package
...
"""

from importlib import import_module

try:
    _spinops = import_module("spinOps._spinOps")
except ModuleNotFoundError as err:
    raise ImportError(
        "The compiled extension module 'spinOps._spinOps' was not found. "
        "Build the package in place (e.g. `pip install -e .`) or install it "
        "from a wheel before importing SpinOps."
    ) from err

# ---------------------------------------------------------------------------
# Public symbols re-exported from the extension
# ---------------------------------------------------------------------------
(
    clebsch,
    tlm,
    unit_tlm,
    number_of_states,
    create_single_spin_Ix,
    create_single_spin_Iy,
    create_single_spin_Iz,
    create_single_spin_Ip,
    create_single_spin_Im,
    create_single_spin_Tlm,
    create_single_spin_Tlm_unit,
    create_rho1,
    create_rho2,
    wigner_d,
    DLM,
    Rotate,
    createEf,
    create_Ixf,
    create_Iyf,
    create_Izf,
    create_Ipf,
    create_Imf,
) = (
    _spinops.clebsch,
    _spinops.tlm,
    _spinops.unit_tlm,
    _spinops.number_of_states,
    _spinops.create_single_spin_Ix,
    _spinops.create_single_spin_Iy,
    _spinops.create_single_spin_Iz,
    _spinops.create_single_spin_Ip,
    _spinops.create_single_spin_Im,
    _spinops.create_single_spin_Tlm,
    _spinops.create_single_spin_Tlm_unit,
    _spinops.create_rho1,
    _spinops.create_rho2,
    _spinops.wigner_d,
    _spinops.DLM,
    _spinops.Rotate,
    _spinops.createEf,
    _spinops.create_Ixf,
    _spinops.create_Iyf,
    _spinops.create_Izf,
    _spinops.create_Ipf,
    _spinops.create_Imf,
)

# --------------------- Reset __module__ attributes -------------------------
for func in (
    clebsch,
    tlm,
    unit_tlm,
    number_of_states,
    create_single_spin_Ix,
    create_single_spin_Iy,
    create_single_spin_Iz,
    create_single_spin_Ip,
    create_single_spin_Im,
    create_single_spin_Tlm,
    create_single_spin_Tlm_unit,
    create_rho1,
    create_rho2,
    wigner_d,
    DLM,
    Rotate,
    createEf,
    create_Ixf,
    create_Iyf,
    create_Izf,
    create_Ipf,
    create_Imf,
):
    func.__module__ = __name__

__all__ = [
    "number_of_states",
    "create_single_spin_Ix",
    "create_single_spin_Iy",
    "create_single_spin_Iz",
    "create_single_spin_Ip",
    "create_single_spin_Im",
    "create_single_spin_Tlm",
    "tlm",
    "create_single_spin_Tlm_unit",
    "unit_tlm",
    "createEf",
    "create_Ixf",
    "create_Iyf",
    "create_Izf",
    "create_Ipf",
    "create_Imf",
    "clebsch",
    "create_rho1",
    "create_rho2",
    "wigner_d",
    "DLM",
    "Rotate",
]

__version__ = "0.1.1"
