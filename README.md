# spinOps

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
