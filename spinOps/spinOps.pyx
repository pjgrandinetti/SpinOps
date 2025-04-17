# cython: language_level=3
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

from numpy cimport ndarray
cimport numpy as cnp
import numpy as np

# Importing external C functions
cdef extern from "complex.h":
    pass

cdef extern from "spin.h":
    # Clebsch-Gordan coefficients and tensor operators
    double clebsch_(double j1, double m1, double j2, double m2, double j, double m)
    double tlm_(double l, double m, double j1, double m1, double j2, double m2)
    double unit_tlm_(double l, double m, double j1, double m1, double j2, double m2)

    # State and operator-related functions
    int numberOfStates_(int spinCount, int *spinsTimesTwo)
    void getIx_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIy_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIz_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIp_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)

    # Tensor operators
    void getTlm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)
    void getTlm_unit_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)

    # Fictitious spin-1/2 operators
    void getEf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
    void getIxf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
    void getIyf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
    void getIzf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
    void getIpf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
    void getImf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)

cpdef double clebsch(double j1, double m1, double j2, double m2, double j, double m):
    """
    Computes the Clebsch-Gordan coefficient for the given quantum numbers.

    Parameters:
        j1 (double): Total angular momentum of the first particle.
        m1 (double): Magnetic quantum number of the first particle.
        j2 (double): Total angular momentum of the second particle.
        m2 (double): Magnetic quantum number of the second particle.
        j (double): Total angular momentum of the combined system.
        m (double): Magnetic quantum number of the combined system.

    Returns:
        double: The Clebsch-Gordan coefficient for the specified quantum numbers.

    Raises:
        ValueError: If the input quantum numbers do not satisfy the selection rules.
    """
    # Validate input quantum numbers (optional, based on selection rules)
    if m1 + m2 != m:
        raise ValueError("The magnetic quantum numbers m1 and m2 must sum to m.")
    if abs(j1 - j2) > j or j > j1 + j2:
        raise ValueError("The total angular momentum j must satisfy |j1 - j2| <= j <= j1 + j2.")

    # Call the external C function
    return clebsch_(j1, m1, j2, m2, j, m)


cpdef double tlm(double l, double m, double j1, double m1, double j2, double m2):
    """
    Computes the tensor operator matrix element for the given quantum numbers.

    Parameters:
        l (double): Rank of the tensor operator.
        m (double): Magnetic quantum number of the tensor operator.
        j1 (double): Total angular momentum of the first particle.
        m1 (double): Magnetic quantum number of the first particle.
        j2 (double): Total angular momentum of the second particle.
        m2 (double): Magnetic quantum number of the second particle.

    Returns:
        double: The tensor operator matrix element for the specified quantum numbers.

    Raises:
        ValueError: If the input quantum numbers do not satisfy the selection rules.
    """
    # Optional validation for quantum number selection rules
    if m1 + m2 != m:
        raise ValueError("The magnetic quantum numbers m1 and m2 must sum to m.")
    if abs(j1 - j2) > l or l > j1 + j2:
        raise ValueError("The rank l must satisfy |j1 - j2| <= l <= j1 + j2.")

    # Call the external C function
    return tlm_(l, m, j1, m1, j2, m2)

cpdef double unit_tlm(double l, double m, double j1, double m1, double j2, double m2):
    """
    Computes the unit tensor operator matrix element for the given quantum numbers.

    Parameters:
        l (double): Rank of the tensor operator.
        m (double): Magnetic quantum number of the tensor operator.
        j1 (double): Total angular momentum of the first particle.
        m1 (double): Magnetic quantum number of the first particle.
        j2 (double): Total angular momentum of the second particle.
        m2 (double): Magnetic quantum number of the second particle.

    Returns:
        double: The unit tensor operator matrix element for the specified quantum numbers.

    Raises:
        ValueError: If the input quantum numbers do not satisfy the selection rules.
    """
    # Optional validation for quantum number selection rules
    if m1 + m2 != m:
        raise ValueError("The magnetic quantum numbers m1 and m2 must sum to m.")
    if abs(j1 - j2) > l or l > j1 + j2:
        raise ValueError("The rank l must satisfy |j1 - j2| <= l <= j1 + j2.")

    # Call the external C function
    return unit_tlm_(l, m, j1, m1, j2, m2)


cpdef int numberOfStates(list spinsTimesTwo):
    """
    Computes the total number of quantum states in a spin system.

    Parameters:
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        int: The total number of quantum states in the spin system.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")

    # Convert the input list to a NumPy array
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)

    # Call the external C function
    return numberOfStates_(spinCount, &spins[0])


cpdef ndarray[double complex, ndim=2] createIx(int spinIndex, list spinsTimesTwo):
    """
    Creates the Ix operator matrix for a given spin in a spin system.

    Parameters:
        spinIndex (int): The index of the spin for which the Ix operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Ix operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIx_(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIy(int spinIndex, list spinsTimesTwo):
    """
    Creates the Iy operator matrix for a given spin in a spin system.

    Parameters:
        spinIndex (int): The index of the spin for which the Iy operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Iy operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIy_(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIz(int spinIndex, list spinsTimesTwo):
    """
    Creates the Iz operator matrix for a given spin in a spin system.

    Parameters:
        spinIndex (int): The index of the spin for which the Iz operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Iz operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIz_(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIp(int spinIndex, list spinsTimesTwo):
    """
    Creates the I+ (Iplus) operator matrix for a given spin in a spin system.

    Parameters:
        spinIndex (int): The index of the spin for which the I+ operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The I+ operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIp_(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIm(int spinIndex, list spinsTimesTwo):
    """
    Creates the I- (Iminus) operator matrix for a given spin in a spin system.

    Parameters:
        spinIndex (int): The index of the spin for which the I- operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The I- operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIm_(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createTLM(int L, int M, int spinIndex, list spinsTimesTwo):
    """
    Creates the TLM operator matrix for a given spin in a spin system.

    Parameters:
        L (int): The rank of the tensor operator.
        M (int): The magnetic quantum number of the tensor operator.
        spinIndex (int): The index of the spin for which the TLM operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The TLM operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getTlm_(&myOp[0, 0], spinIndex, &spins[0], spinCount, L, M)

    return myOp


cpdef ndarray[double complex, ndim=2] createTLM_unit(int L, int M, int spinIndex, list spinsTimesTwo):
    """
    Creates the unit TLM operator matrix for a given spin in a spin system.

    Parameters:
        L (int): The rank of the tensor operator.
        M (int): The magnetic quantum number of the tensor operator.
        spinIndex (int): The index of the spin for which the unit TLM operator is being created.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The unit TLM operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `spinIndex` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getTlm_unit_(&myOp[0, 0], spinIndex, &spins[0], spinCount, L, M)

    return myOp


cpdef ndarray[double complex, ndim=2] createEf(int r, int s, list spinsTimesTwo):
    """
    Creates the Ef (identity) operator matrix for a fictitious spin-1/2 system.

    Parameters:
        r (int): The index of the first state.
        s (int): The index of the second state.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Ef operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `r` or `s` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getEf_(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIxf(int r, int s, list spinsTimesTwo):
    """
    Creates the Ixf operator matrix for a fictitious spin-1/2 system.

    Parameters:
        r (int): The index of the first state.
        s (int): The index of the second state.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Ixf operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `r` or `s` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIxf_(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIyf(int r, int s, list spinsTimesTwo):
    """
    Creates the Iyf operator matrix for a fictitious spin-1/2 system.

    Parameters:
        r (int): The index of the first state.
        s (int): The index of the second state.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Iyf operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `r` or `s` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIyf_(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIzf(int r, int s, list spinsTimesTwo):
    """
    Creates the Izf operator matrix for a fictitious spin-1/2 system.

    Parameters:
        r (int): The index of the first state.
        s (int): The index of the second state.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The Izf operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `r` or `s` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIzf_(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIpf(int r, int s, list spinsTimesTwo):
    """
    Creates the I+f (Iplus fictitious) operator matrix for a fictitious spin-1/2 system.

    Parameters:
        r (int): The index of the first state.
        s (int): The index of the second state.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The I+f operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `r` or `s` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getIpf_(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createImf(int r, int s, list spinsTimesTwo):
    """
    Creates the I-f (Iminus fictitious) operator matrix for a fictitious spin-1/2 system.

    Parameters:
        r (int): The index of the first state.
        s (int): The index of the second state.
        spinsTimesTwo (list): A list of integers representing `2 * I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The I-f operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `spinsTimesTwo` is empty.
        IndexError: If `r` or `s` is out of bounds.
    """
    # Validate input
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    getImf_(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cdef extern from "spatial.h":
    # Functions for creating irreducible spherical tensors
    void getrho1_pas_(double complex *tensor, double zeta)  # Creates the rank-1 irreducible tensor (rho1) in PAS
    void getrho2_pas_(double complex *tensor, double zeta, double eta)  # Creates the rank-2 irreducible tensor (rho2) in PAS

    # Wigner rotation functions
    double wigner_d_(double l, double m1, double m2, double beta)  # Computes the Wigner d-matrix element
    double complex DLM_(double l, double m1, double m2, double alpha, double beta, double gamma)  # Computes the Wigner D-matrix element

    # Rotation of quantum states
    void Rot_(double j, double complex *initial, double alpha, double beta, double gamma, double complex *final)  # Rotates a quantum state


cpdef cnp.ndarray[double complex, ndim=1] createRho1(double zeta):
    """
    Creates the rank-1 irreducible tensor (rho1) in the principal axis system (PAS).

    Parameters:
        zeta (double): The asymmetry parameter for the tensor.

    Returns:
        cnp.ndarray[double complex, ndim=1]: A 1D NumPy array representing the rank-1 irreducible tensor.

    Raises:
        ValueError: If the input parameter `zeta` is invalid (optional validation can be added if needed).
    """
    # Allocate memory for the tensor
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(3, dtype=np.complex128)

    # Call the external C function to populate the tensor
    getrho1_pas_(<double complex *> cnp.PyArray_DATA(myOp), zeta)

    return myOp


cpdef cnp.ndarray[double complex, ndim=1] createRho2(double zeta, double eta):
    """
    Creates the rank-2 irreducible tensor (rho2) in the principal axis system (PAS).

    Parameters:
        zeta (double): The asymmetry parameter for the tensor.
        eta (double): The eta parameter for the tensor.

    Returns:
        cnp.ndarray[double complex, ndim=1]: A 1D NumPy array representing the rank-2 irreducible tensor.

    Raises:
        ValueError: If the input parameters `zeta` or `eta` are invalid (optional validation can be added if needed).
    """
    # Allocate memory for the tensor
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(5, dtype=np.complex128)

    # Call the external C function to populate the tensor
    getrho2_pas_(<double complex *> cnp.PyArray_DATA(myOp), zeta, eta)

    return myOp


cpdef double wigner_d(double l, double m1, double m2, double beta):
    """
    Computes the Wigner d-matrix element for the given quantum numbers and rotation angle.

    Parameters:
        l (double): The rank of the rotation operator.
        m1 (double): The initial magnetic quantum number.
        m2 (double): The final magnetic quantum number.
        beta (double): The rotation angle in radians.

    Returns:
        double: The Wigner d-matrix element for the specified parameters.

    Raises:
        ValueError: If the input quantum numbers do not satisfy the selection rules (optional validation can be added).
    """
    # Call the external C function to compute the Wigner d-matrix element
    return wigner_d_(l, m1, m2, beta)


cpdef double complex DLM(double l, double m1, double m2, double alpha, double beta, double gamma):
    """
    Computes the Wigner D-matrix element for the given quantum numbers and rotation angles.

    Parameters:
        l (double): The rank of the rotation operator.
        m1 (double): The initial magnetic quantum number.
        m2 (double): The final magnetic quantum number.
        alpha (double): The first Euler angle (rotation about the z-axis) in radians.
        beta (double): The second Euler angle (rotation about the y-axis) in radians.
        gamma (double): The third Euler angle (rotation about the z-axis) in radians.

    Returns:
        double complex: The Wigner D-matrix element for the specified parameters.

    Raises:
        ValueError: If the input quantum numbers do not satisfy the selection rules (optional validation can be added).
    """
    # Call the external C function to compute the Wigner D-matrix element
    return DLM_(l, m1, m2, alpha, beta, gamma)


cpdef cnp.ndarray[double complex, ndim=1] Rotate(cnp.ndarray[double complex, ndim=1] initial, double alpha, double beta, double gamma):
    """
    Rotates a quantum state using the Wigner D-matrix.

    Parameters:
        initial (cnp.ndarray[double complex, ndim=1]): A 1D NumPy array representing the initial quantum state.
        alpha (double): The first Euler angle (rotation about the z-axis) in radians.
        beta (double): The second Euler angle (rotation about the y-axis) in radians.
        gamma (double): The third Euler angle (rotation about the z-axis) in radians.

    Returns:
        cnp.ndarray[double complex, ndim=1]: A 1D NumPy array representing the rotated quantum state.

    Raises:
        ValueError: If the input array `initial` is empty.
    """
    # Validate input
    if len(initial) == 0:
        raise ValueError("The input array 'initial' cannot be empty.")

    # Allocate memory for the rotated state
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(len(initial), dtype=np.complex128)

    # Call the external C function to perform the rotation
    Rot_((len(initial) - 1) / 2, 
         <double complex *> cnp.PyArray_DATA(initial), 
         alpha, beta, gamma, 
         <double complex *> cnp.PyArray_DATA(myOp))

    return myOp
