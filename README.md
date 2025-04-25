# spinOps

|              |                                                                                                                                                                                                                                                                                                                                                                            |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://img.shields.io/pypi/v/spinops.svg?style=flat&logo=pypi&logoColor=white)](https://pypi.org/project/spinops/) [![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spinops.svg)](https://pypi.org/project/spinops/) |
| Build Status | [![CI](https://github.com/pjgrandinetti/spinOps/actions/workflows/wheels.yml/badge.svg?branch=master)](https://github.com/pjgrandinetti/spinOps/actions/workflows/wheels.yml) [![Read the Docs](https://img.shields.io/readthedocs/spinOps)](https://spinops.readthedocs.io/en/latest/) |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                  |
| Metrics      | [![codecov](https://codecov.io/gh/pjgrandinetti/spinOps/branch/master/graph/badge.svg)](https://codecov.io/gh/pjgrandinetti/spinOps) [![CodeFactor](https://www.codefactor.io/repository/github/pjgrandinetti/spinops/badge)](https://www.codefactor.io/repository/github/pjgrandinetti/spinOps)                                                                           |


**spinOps** is a Python package for performing operations on quantum spin systems. It provides tools for creating and manipulating spin operators, irreducible spherical tensors, and performing quantum state rotations using Wigner D-matrix elements.

[![Documentation Status](https://readthedocs.org/projects/spinops/badge/?version=latest)](https://spinops.readthedocs.io/en/latest/)

## Features

- Generate Clebsch-Gordan coefficients and tensor operators.
- Create spin operator matrices (e.g., \(I_x\), \(I_y\), \(I_z\), \(I_+\), \(I_-\)).
- Handle fictitious spin-1/2 systems.
- Compute irreducible spherical tensors in PAS (\(\rho_1\), \(\rho_2\)).
- Perform rotations using Wigner D-matrix elements.

## Installation

Install the package from PyPI:

```bash
pip install spinOps
```

## Example

```python
import spinOps

# Create a spin operator matrix for a system of two coupled spin-1/2 nuclei.
# In this example:
#   - The first argument (0) specifies the operator index (e.g., I_x for the first nucleus).
#   - The second argument [1, 1] indicates the spin value for each nucleus, where
#     each entry is twice the actual spin (1 corresponds to 1/2 in this context).
operator = spinOps.createIx(0, [1, 1])

# Print the generated operator matrix to verify its structure and values.
print(operator)
```
