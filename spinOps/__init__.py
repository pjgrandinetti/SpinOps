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
    numberOfStates,
    createIx,
    createIy,
    createIz,
    createIp,
    createIm,
    createTLM,
    createTLM_unit,
    createRho1,
    createRho2,
    wigner_d,
    DLM,
    Rotate,
    createEf,
    createIxf,
    createIyf,
    createIzf,
    createIpf,
    createImf,
) = (
    _spinops.clebsch,
    _spinops.tlm,
    _spinops.unit_tlm,
    _spinops.numberOfStates,
    _spinops.createIx,
    _spinops.createIy,
    _spinops.createIz,
    _spinops.createIp,
    _spinops.createIm,
    _spinops.createTLM,
    _spinops.createTLM_unit,
    _spinops.createRho1,
    _spinops.createRho2,
    _spinops.wigner_d,
    _spinops.DLM,
    _spinops.Rotate,
    _spinops.createEf,
    _spinops.createIxf,
    _spinops.createIyf,
    _spinops.createIzf,
    _spinops.createIpf,
    _spinops.createImf,
)

# --------------------- Reset __module__ attributes -------------------------
for func in (
    clebsch,
    tlm,
    unit_tlm,
    numberOfStates,
    createIx,
    createIy,
    createIz,
    createIp,
    createIm,
    createTLM,
    createTLM_unit,
    createRho1,
    createRho2,
    wigner_d,
    DLM,
    Rotate,
    createEf,
    createIxf,
    createIyf,
    createIzf,
    createIpf,
    createImf,
):
    func.__module__ = __name__

__all__ = [
    "numberOfStates",
    "createIx",
    "createIy",
    "createIz",
    "createIp",
    "createIm",
    "createTLM",
    "tlm",
    "createTLM_unit",
    "unit_tlm",
    "createEf",
    "createIxf",
    "createIyf",
    "createIzf",
    "createIpf",
    "createImf",
    "clebsch",
    "createRho1",
    "createRho2",
    "wigner_d",
    "DLM",
    "Rotate",
]

__version__ = "0.1.0"

