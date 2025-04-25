# File: spinOps/spinOps.pyx
# cython: language_level=3
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

from numpy cimport ndarray
cimport numpy as cnp
import numpy as np

from .spinOps cimport clebsch_ as _clebsch
from .spinOps cimport tlm_ as _tlm
from .spinOps cimport unit_tlm_ as _unit_tlm
from .spinOps cimport numberOfStates_ as _numberOfStates

from .spinOps cimport getIx_ as _getIx
from .spinOps cimport getIy_ as _getIy
from .spinOps cimport getIz_ as _getIz
from .spinOps cimport getIp_ as _getIp
from .spinOps cimport getIm_ as _getIm

from .spinOps cimport getTlm_ as _getTlm
from .spinOps cimport getTlm_unit_ as _getTlm_unit

from .spinOps cimport getEf_ as _getEf
from .spinOps cimport getIxf_ as _getIxf
from .spinOps cimport getIyf_ as _getIyf
from .spinOps cimport getIzf_ as _getIzf
from .spinOps cimport getIpf_ as _getIpf
from .spinOps cimport getImf_ as _getImf

from .spinOps cimport mypow, fac

from .spinOps cimport getrho1_pas_ as _getrho1_pas
from .spinOps cimport getrho2_pas_ as _getrho2_pas
from .spinOps cimport wigner_d_ as _wigner_d
from .spinOps cimport DLM_ as _DLM
from .spinOps cimport Rot_ as _Rot

cpdef double clebsch(double j1, double m1, double j2, double m2, double j, double m):
    """
    Computes the Clebsch-Gordan coefficient, :math:`\langle j,m|j_1,m_1,j_2,m_2\\rangle` for the specified quantum numbers.

    Parameters
    ----------
    j1 : double
        Total angular momentum of the first particle.
    m1 : double
        Magnetic quantum number of the first particle.
    j2 : double
        Total angular momentum of the second particle.
    m2 : double
        Magnetic quantum number of the second particle.
    j : double
        Total angular momentum of the combined system.
    m : double
        Magnetic quantum number of the combined system.

    Returns
    -------
    double
        The Clebsch-Gordan coefficient for the specified quantum numbers.

    Raises
    ------
    ValueError
        If the input quantum numbers do not satisfy the required selection rules.
    """
    if m1 + m2 != m:
        raise ValueError("The magnetic quantum numbers m1 and m2 must sum to m.")
    if abs(j1 - j2) > j or j > j1 + j2:
        raise ValueError("The total angular momentum j must satisfy |j1 - j2| <= j <= j1 + j2.")

    return _clebsch(j1, m1, j2, m2, j, m)


cpdef double tlm(double l, double m, double j1, double m1, double j2, double m2):
    """
    Computes the matrix element

    .. math::

        \langle j_1,m_1|\:\hat{T}_{l,m}\:|j_2,m_2\\rangle

    of the irreducible spherical tensor operator :math:`\hat{T}_{l,m}`.

    Parameters
    ----------
    l : double
        Rank of the irreducible spherical tensor operator.
    m : double
        Order of the irreducible spherical tensor operator.
    j1 : double
        Total angular momentum quantum number of the first particle.
    m1 : double
        Angular momentum component quantum number of the first particle.
    j2 : double
        Total angular momentum quantum number of the second particle.
    m2 : double
        Angular momentum component quantum number of the second particle.

    Returns
    -------
    double
        The matrix element :math:`\langle j_1,m_1|\:\hat{T}_{l,m}\:|j_2,m_2\\rangle`.

    Raises
    ------
    ValueError
        If the quantum numbers do not satisfy selection rules.
    """
    if m1 + m2 != m:
        raise ValueError("Quantum numbers m1 and m2 must sum to m.")
    if abs(j1 - j2) > l or l > j1 + j2:
        raise ValueError("Rank l must satisfy |j1 - j2| <= l <= j1 + j2.")

    return _tlm(l, m, j1, m1, j2, m2)

cpdef double unit_tlm(double l, double m, double j1, double m1, double j2, double m2):
    """
    Computes the matrix element

    .. math::

        \langle j_1,m_1|\:\hat{\mathcal{T}}_{l,m}\:|j_2,m_2\\rangle

    of the unit irreducible spherical tensor operator 
    :math:`\hat{\mathcal{T}}_{l,m}` between the specified quantum states.

    Parameters:
        l (double): The rank of the irreducible spherical tensor operator operator.
        m (double): The order of the irreducible spherical tensor operator operator.
        j1 (double): Total angular momentum quantum number of the first particle.
        m1 (double): Angular momentum component quantum number of the first particle.
        j2 (double): Total angular momentum quantum number of the second particle.
        m2 (double): Angular momentum component of the second particle.

    Returns:
        double: The :math:`\langle j_1,m_1|\:\hat{\mathcal{T}}_{l,m}\:|j_2,m_2\\rangle` matrix element.

    Raises:
        ValueError: If the input quantum numbers do not satisfy the selection rules.
    """
    # Optional validation for quantum number selection rules
    if m1 + m2 != m:
        raise ValueError("The magnetic quantum numbers m1 and m2 must sum to m.")
    if abs(j1 - j2) > l or l > j1 + j2:
        raise ValueError("The rank l must satisfy |j1 - j2| <= l <= j1 + j2.")

    # Call the external C function
    return _unit_tlm(l, m, j1, m1, j2, m2)


cpdef int numberOfStates(list spinsTimesTwo):
    """
    Computes the total number of quantum states for a given spin system.

    Parameters
    ----------
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    int
        Total number of quantum states in the spin system.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")

    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)

    return _numberOfStates(spinCount, &spins[0])
cpdef ndarray[double complex, ndim=2] createIx(int spinIndex, list spinsTimesTwo):
    """
    Generates the single-spin :math:`\hat{I}_x` operator matrix for the specified spin within a spin system.

    Parameters
    ----------
    spinIndex : int
        Index of the spin for which the :math:`\hat{I}_x` operator is constructed.
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{I}_x` operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `spinIndex` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getIx(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createIy(int spinIndex, list spinsTimesTwo):
    """
    Generates the single-spin :math:`\hat{I}_y` operator matrix for a specified spin within a spin system.

    Parameters
    ----------
    spinIndex : int
        Index of the spin for which the :math:`\hat{I}_y` operator is constructed.
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{I}_y` operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `spinIndex` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getIy(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createIz(int spinIndex, list spinsTimesTwo):
    """
    Creates the single-spin :math:`\hat{I}_z` operator matrix for a single spin in a spin system.

    Parameters:
        spinIndex (int): The index of the spin for which the :math:`\hat{I}_z` operator is being created.
        spinsTimesTwo (list): A list of integers representing :math:`2 I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The :math:`\hat{I}_z` operator matrix as a 2D NumPy array.

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
    _getIz(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIp(int spinIndex, list spinsTimesTwo):
    """
    Generates the single-spin raising operator (:math:`\hat{I}_+`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    spinIndex : int
        Index of the spin for which the :math:`\hat{I}_+` operator is constructed.
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the raising (:math:`\hat{I}_+`) operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `spinIndex` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getIp(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createIm(int spinIndex, list spinsTimesTwo):
    """
    Generates the single-spin lowering operator (:math:`\hat{I}_-`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    spinIndex : int
        Index of the spin for which the :math:`\hat{I}_-` operator is constructed.
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the lowering (:math:`\hat{I}_-`) operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `spinIndex` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getIm(&myOp[0, 0], spinIndex, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createTLM(int L, int M, int spinIndex, list spinsTimesTwo):
    """
    Generates the single-spin irreducible spherical tensor operator (:math:`\hat{T}_{L,M}`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    L : int
        Rank of the tensor operator.
    M : int
        Order of the tensor operator.
    spinIndex : int
        Index of the spin for which the :math:`\hat{T}_{L,M}` operator is constructed.
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{T}_{L,M}` operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `spinIndex` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getTlm(<double complex *> cnp.PyArray_DATA(myOp), spinIndex, &spins[0], spinCount, L, M)

    return myOp


cpdef ndarray[double complex, ndim=2] createTLM_unit(int L, int M, int spinIndex, list spinsTimesTwo):
    """
    Generates the single-spin unit-normalized irreducible spherical tensor operator (:math:`\hat{\mathcal{T}}_{L,M}`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    L : int
        Rank of the tensor operator.
    M : int
        Order of the tensor operator.
    spinIndex : int
        Index of the spin for which the unit-normalized :math:`\hat{\mathcal{T}}_{L,M}` operator is constructed.
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the unit-normalized :math:`\hat{\mathcal{T}}_{L,M}` operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `spinIndex` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if spinIndex < 0 or spinIndex >= len(spinsTimesTwo):
        raise IndexError("The spinIndex is out of bounds.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getTlm_unit(&myOp[0, 0], spinIndex, &spins[0], spinCount, L, M)

    return myOp


cpdef ndarray[double complex, ndim=2] createEf(int r, int s, list spinsTimesTwo):
    """
    Generates the operator matrix :math:`\hat{E}^{r-s}` corresponding to the transition from state :math:`s` to :math:`r`
    in a fictitious spin-1/2 system.

    Parameters
    ----------
    r : int
        Index of the first quantum state (row index).
    s : int
        Index of the second quantum state (column index).
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{E}^{r-s}` operator matrix.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `r` or `s` is out of the valid range.
    """
    if not spinsTimesTwo:
        raise ValueError("The input list 'spinsTimesTwo' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    cdef int nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins = np.array(spinsTimesTwo, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getEf(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp


cpdef ndarray[double complex, ndim=2] createIxf(int r, int s, list spinsTimesTwo):
    """
    Generates the fictitious spin-1/2 operator matrix :math:`\hat{I}_x^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_x^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `r` or `s` is negative.
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
    _getIxf(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createIyf(int r, int s, list spinsTimesTwo):
    """
    Generates the fictitious spin-1/2 operator matrix :math:`\hat{I}_y^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_y^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `r` or `s` is negative.
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
    _getIyf(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createIzf(int r, int s, list spinsTimesTwo):
    """
    Generates the fictitious spin-1/2 operator matrix :math:`\hat{I}_z^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_z^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `r` or `s` is negative.
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
    _getIzf(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createIpf(int r, int s, list spinsTimesTwo):
    """
    Generates the fictitious spin-1/2 raising operator matrix :math:`\hat{I}_+^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_+^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `r` or `s` is negative.
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
    _getIpf(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp

cpdef ndarray[double complex, ndim=2] createImf(int r, int s, list spinsTimesTwo):
    """
    Generates the fictitious spin-1/2 lowering operator matrix :math:`\hat{I}_-^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    spinsTimesTwo : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_-^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `spinsTimesTwo` is empty.
    IndexError
        If `r` or `s` is negative.
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
    _getImf(&myOp[0, 0], r, s, &spins[0], spinCount)

    return myOp

cpdef cnp.ndarray[double complex, ndim=1] createRho1(double zeta):
    """
    Constructs the rank-1 irreducible spherical tensor :math:`\\rho_{1,m}` in the principal axis system (PAS)
    according to the Haeberlen convention.

    Parameters
    ----------
    zeta : double
        The anisotropy parameter :math:`\zeta` for the tensor.

    Returns
    -------
    cnp.ndarray[double complex, ndim=1]
        A 1D NumPy array containing the components of the rank-1 irreducible tensor.

    Raises
    ------
    ValueError
        If the input parameter `zeta` is invalid (validation not currently enforced).
    """
    # Allocate memory for the tensor
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(3, dtype=np.complex128)

    # Call the external C function to populate the tensor
    _getrho1_pas(<double complex *> cnp.PyArray_DATA(myOp), zeta)

    return myOp

cpdef cnp.ndarray[double complex, ndim=1] createRho2(double zeta, double eta):
    """
    Constructs the rank-2 irreducible spherical tensor :math:`\\rho_{2,m}` in the principal axis system (PAS)
    according to the Haeberlen convention.

    Parameters
    ----------
    zeta : double
        The anisotropy parameter :math:`\zeta` for the tensor.
    eta : double
        The asymmetry parameter :math:`\eta` for the tensor.

    Returns
    -------
    cnp.ndarray[double complex, ndim=1]
        A 1D NumPy array containing the components of the rank-2 irreducible tensor.

    Raises
    ------
    ValueError
        If the input parameters `zeta` or `eta` are invalid (validation not currently enforced).
    """
    # Allocate memory for the tensor
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(5, dtype=np.complex128)

    # Call the external C function to populate the tensor
    _getrho2_pas(<double complex *> cnp.PyArray_DATA(myOp), zeta, eta)

    return myOp

cpdef double wigner_d(double l, double m1, double m2, double beta):
    """
    Computes the reduced Wigner d-matrix element :math:`d^{(l)}_{m_1,m_2}[\\beta]` for the given quantum numbers
    and rotation angle.

    Parameters
    ----------
    l : double
        Rank of the rotation operator.
    m1 : double
        Initial magnetic quantum number.
    m2 : double
        Final magnetic quantum number.
    beta : double
        Rotation angle in radians.

    Returns
    -------
    double
        The reduced Wigner d-matrix element :math:`d^{(l)}_{m_1,m_2}[\\beta]`.

    Raises
    ------
    ValueError
        If the input quantum numbers do not satisfy the required selection rules (validation not currently enforced).
    """
    # Call the external C function to compute the Wigner d-matrix element
    return _wigner_d(l, m1, m2, beta)

cpdef double complex DLM(double l, double m1, double m2, double alpha, double beta, double gamma):
    """
    Computes the Wigner D-matrix element :math:`\mathcal{D}^{(l)}_{m_1,m_2}[\\alpha,\\beta,\gamma]`
    for the given quantum numbers and Euler angles.

    Parameters
    ----------
    l : double
        Rank of the rotation operator.
    m1 : double
        Initial magnetic quantum number.
    m2 : double
        Final magnetic quantum number.
    alpha : double
        First Euler angle (rotation about the z-axis) in radians.
    beta : double
        Second Euler angle (rotation about the y-axis) in radians.
    gamma : double
        Third Euler angle (rotation about the z-axis) in radians.

    Returns
    -------
    double complex
        The Wigner D-matrix element :math:`\mathcal{D}^{(l)}_{m_1,m_2}[\\alpha,\\beta,\gamma]`.

    Raises
    ------
    ValueError
        If the input quantum numbers do not satisfy the required selection rules (validation not enforced).
    """
    return _DLM(l, m1, m2, alpha, beta, gamma)


cpdef cnp.ndarray[double complex, ndim=1] Rotate(cnp.ndarray[double complex, ndim=1] initial, double alpha, double beta, double gamma):
    """
    Rotates a spherical tensor :math:`\\rho_{l,m}` using the Wigner D-matrix and specified Euler angles.

    Parameters
    ----------
    initial : cnp.ndarray[double complex, ndim=1]
        A 1D NumPy array representing the spherical tensor components.
    alpha : double
        First Euler angle (rotation about the z-axis) in radians.
    beta : double
        Second Euler angle (rotation about the y-axis) in radians.
    gamma : double
        Third Euler angle (rotation about the z-axis) in radians.

    Returns
    -------
    cnp.ndarray[double complex, ndim=1]
        A 1D NumPy array representing the rotated spherical tensor.

    Raises
    ------
    ValueError
        If the input array `initial` is empty.
    """
    # Validate input
    if len(initial) == 0:
        raise ValueError("The input array 'initial' cannot be empty.")

    # Allocate memory for the rotated state
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(len(initial), dtype=np.complex128)

    # Call the external C function to perform the rotation
    _Rot((len(initial) - 1) / 2, 
         <double complex *> cnp.PyArray_DATA(initial), 
         alpha, beta, gamma, 
         <double complex *> cnp.PyArray_DATA(myOp))

    return myOp
