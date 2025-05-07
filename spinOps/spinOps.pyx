# File: spinOps/spinOps.pyx
# cython: language_level=3
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

cimport numpy as cnp
from numpy cimport ndarray

import numpy as np

from .spinOps cimport DLM_ as _DLM
from .spinOps cimport Rot_ as _Rot
from .spinOps cimport clebsch_ as _clebsch
from .spinOps cimport get_single_spin_Im_ as _get_single_spin_Im
from .spinOps cimport get_single_spin_Ip_ as _get_single_spin_Ip
from .spinOps cimport get_single_spin_Ix_ as _get_single_spin_Ix
from .spinOps cimport get_single_spin_Iy_ as _get_single_spin_Iy
from .spinOps cimport get_single_spin_Iz_ as _get_single_spin_Iz
from .spinOps cimport get_single_spin_Tlm_ as _get_single_spin_Tlm
from .spinOps cimport get_single_spin_Tlm_unit_ as _get_single_spin_Tlm_unit
from .spinOps cimport getEf_ as _getEf
from .spinOps cimport getImf_ as _get_Imf
from .spinOps cimport getIpf_ as _get_Ipf
from .spinOps cimport getIxf_ as _get_Ixf
from .spinOps cimport getIyf_ as _get_Iyf
from .spinOps cimport getIzf_ as _get_Izf
from .spinOps cimport getrho1_pas_ as _getrho1_pas
from .spinOps cimport getrho2_pas_ as _getrho2_pas
from .spinOps cimport number_of_states_ as _number_of_states
from .spinOps cimport tlm_ as _tlm
from .spinOps cimport unit_tlm_ as _unit_tlm
from .spinOps cimport wigner_d_ as _wigner_d


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

    """

    return _clebsch(int(2*j1), int(2*m1), int(2*j2), int(2*m2), int(2*j), int(2*m))

cpdef double tlm(int l, int m, double I, double m1, double m2):
    """
    Computes the matrix element

    .. math::

        \langle I,m_1|\:\hat{T}_{l,m}\:|I,m_2\\rangle

    of the irreducible spherical tensor operator :math:`\hat{T}_{l,m}`.

    Parameters
    ----------
    l : double
        Rank of the irreducible spherical tensor operator.
    m : double
        Order of the irreducible spherical tensor operator.
    I : double
        Total angular momentum quantum number of the spin.
    m1 : double
        Angular momentum component quantum number of the spin.
    m2 : double
        Angular momentum component quantum number of the spin.

    Returns
    -------
    double
        The matrix element :math:`\langle I,m_1|\:\hat{T}_{l,m}\:|I,m_2\\rangle`.
    """

    return _tlm(l, m, int(2*I), int(2*m1), int(2*m2))


cpdef double unit_tlm(int l, int m, double I, double m1, double m2):
    """
    Computes the matrix element

    .. math::

        \langle j_1,m_1|\:\hat{\mathcal{T}}_{l,m}\:|j_2,m_2\\rangle

    of the unit irreducible spherical tensor operator 
    :math:`\hat{\mathcal{T}}_{l,m}` between the specified quantum states.

    Parameters:
        l (double): Rank of the irreducible spherical tensor operator.
        m (double): Order of the irreducible spherical tensor operator.
        j1 (double): Total angular momentum quantum number of the first particle.
        m1 (double): Angular momentum component of the first particle.
        j2 (double): Total angular momentum quantum number of the second particle.
        m2 (double): Angular momentum component of the second particle.

    Returns:
        double: The :math:`\langle j_1,m_1|\:\hat{\mathcal{T}}_{l,m}\:|j_2,m_2\\rangle` matrix element.
    """

    return _unit_tlm(l, m, int(2*I), int(2*m1), int(2*m2))


cpdef int number_of_states(list i_times_2):
    """
    Computes the total number of quantum states for a given spin system.

    Parameters
    ----------
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    int
        Total number of quantum states in the spin system.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")

    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)

    return _number_of_states(total_spin_count, &spins[0])


cpdef ndarray[double complex, ndim=2] create_single_spin_Ix(int spin_index, list i_times_2):
    """
    Generates the single-spin :math:`\hat{I}_x` operator matrix for the specified spin within a spin system.

    Parameters
    ----------
    spin_index : int
        Index of the spin for which the :math:`\hat{I}_x` operator is constructed.
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{I}_x` operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `spin_index` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _get_single_spin_Ix(&myOp[0, 0], spin_index, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_single_spin_Iy(int spin_index, list i_times_2):
    """
    Generates the single-spin :math:`\hat{I}_y` operator matrix for a specified spin within a spin system.

    Parameters
    ----------
    spin_index : int
        Index of the spin for which the :math:`\hat{I}_y` operator is constructed.
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{I}_y` operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `spin_index` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _get_single_spin_Iy(&myOp[0, 0], spin_index, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_single_spin_Iz(int spin_index, list i_times_2):
    """
    Creates the single-spin :math:`\hat{I}_z` operator matrix for a single spin in a spin system.

    Parameters:
        spin_index (int): The index of the spin for which the :math:`\hat{I}_z` operator is being created.
        i_times_2 (list): A list of integers representing :math:`2 I` values for each spin in the system,
                              where `I` is the spin quantum number.

    Returns:
        ndarray[double complex, ndim=2]: The :math:`\hat{I}_z` operator matrix as a 2D NumPy array.

    Raises:
        ValueError: If the input list `i_times_2` is empty.
        IndexError: If `spin_index` is out of bounds.
    """
    # Validate input
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    _get_single_spin_Iz(&myOp[0, 0], spin_index, &spins[0], total_spin_count)

    return myOp


cpdef ndarray[double complex, ndim=2] create_single_spin_Ip(int spin_index, list i_times_2):
    """
    Generates the single-spin raising operator (:math:`\hat{I}_+`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    spin_index : int
        Index of the spin for which the :math:`\hat{I}_+` operator is constructed.
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the raising (:math:`\hat{I}_+`) operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `spin_index` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _get_single_spin_Ip(&myOp[0, 0], spin_index, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_single_spin_Im(int spin_index, list i_times_2):
    """
    Generates the single-spin lowering operator (:math:`\hat{I}_-`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    spin_index : int
        Index of the spin for which the :math:`\hat{I}_-` operator is constructed.
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the lowering (:math:`\hat{I}_-`) operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `spin_index` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _get_single_spin_Im(&myOp[0, 0], spin_index, &spins[0], total_spin_count)

    return myOp


cpdef ndarray[double complex, ndim=2] create_single_spin_Tlm(int L, int M, int spin_index, list i_times_2):
    """
    Generates the single-spin irreducible spherical tensor operator (:math:`\hat{T}_{L,M}`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    L : int
        Rank of the tensor operator.
    M : int
        Order of the tensor operator.
    spin_index : int
        Index of the spin for which the :math:`\hat{T}_{L,M}` operator is constructed.
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{T}_{L,M}` operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `spin_index` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _get_single_spin_Tlm(<double complex *> cnp.PyArray_DATA(myOp), spin_index, &spins[0], total_spin_count, L, M)

    return myOp


cpdef ndarray[double complex, ndim=2] create_single_spin_Tlm_unit(int L, int M, int spin_index, list i_times_2):
    """
    Generates the single-spin unit-normalized irreducible spherical tensor operator (:math:`\hat{\mathcal{T}}_{L,M}`) matrix for a specified spin within a spin system.

    Parameters
    ----------
    L : int
        Rank of the tensor operator.
    M : int
        Order of the tensor operator.
    spin_index : int
        Index of the spin for which the unit-normalized :math:`\hat{\mathcal{T}}_{L,M}` operator is constructed.
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the unit-normalized :math:`\hat{\mathcal{T}}_{L,M}` operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `spin_index` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if spin_index < 0 or spin_index >= len(i_times_2):
        raise IndexError("The spin_index is out of bounds.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _get_single_spin_Tlm_unit(&myOp[0, 0], spin_index, &spins[0], total_spin_count, L, M)

    return myOp


cpdef ndarray[double complex, ndim=2] createEf(int r, int s, list i_times_2):
    """
    Generates the operator matrix :math:`\hat{E}^{r-s}` corresponding to the transition from state :math:`s` to :math:`r`
    in a fictitious spin-1/2 system.

    Parameters
    ----------
    r : int
        Index of the first quantum state (row index).
    s : int
        Index of the second quantum state (column index).
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the :math:`\hat{E}^{r-s}` operator matrix.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `r` or `s` is out of the valid range.
    """
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    _getEf(&myOp[0, 0], r, s, &spins[0], total_spin_count)

    return myOp


cpdef ndarray[double complex, ndim=2] create_Ixf(int r, int s, list i_times_2):
    """
    Generates the fictitious spin-1/2 operator matrix :math:`\hat{I}_x^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_x^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `r` or `s` is negative.
    """
    # Validate input
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    _get_Ixf(&myOp[0, 0], r, s, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_Iyf(int r, int s, list i_times_2):
    """
    Generates the fictitious spin-1/2 operator matrix :math:`\hat{I}_y^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_y^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `r` or `s` is negative.
    """
    # Validate input
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    _get_Iyf(&myOp[0, 0], r, s, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_Izf(int r, int s, list i_times_2):
    """
    Generates the fictitious spin-1/2 operator matrix :math:`\hat{I}_z^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_z^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `r` or `s` is negative.
    """
    # Validate input
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    _get_Izf(&myOp[0, 0], r, s, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_Ipf(int r, int s, list i_times_2):
    """
    Generates the fictitious spin-1/2 raising operator matrix :math:`\hat{I}_+^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_+^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `r` or `s` is negative.
    """
    # Validate input
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    _get_Ipf(&myOp[0, 0], r, s, &spins[0], total_spin_count)

    return myOp

cpdef ndarray[double complex, ndim=2] create_Imf(int r, int s, list i_times_2):
    """
    Generates the fictitious spin-1/2 lowering operator matrix :math:`\hat{I}_-^{r-s}` for a transition
    from state :math:`s` to state :math:`r`.

    Parameters
    ----------
    r : int
        Index of the target quantum state (row index).
    s : int
        Index of the source quantum state (column index).
    i_times_2 : list of int
        List of integers representing :math:`2I` values for each spin in the system,
        where :math:`I` is the spin quantum number.

    Returns
    -------
    ndarray[double complex, ndim=2]
        A 2D NumPy array representing the operator :math:`\hat{I}_-^{r-s}`.

    Raises
    ------
    ValueError
        If the input list `i_times_2` is empty.
    IndexError
        If `r` or `s` is negative.
    """
    # Validate input
    if not i_times_2:
        raise ValueError("The input list 'i_times_2' cannot be empty.")
    if r < 0 or s < 0:
        raise IndexError("State indices 'r' and 's' must be non-negative.")

    # Compute the number of states and prepare the operator matrix
    cdef int nstates = number_of_states(i_times_2)
    cdef int total_spin_count = len(i_times_2)
    cdef ndarray[int] spins = np.array(i_times_2, dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates, nstates), dtype=np.complex128)

    # Call the external C function to populate the operator matrix
    _get_Imf(&myOp[0, 0], r, s, &spins[0], total_spin_count)

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
    """
    return _wigner_d(int(2*l), int(2*m1), int(2*m2), beta)

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

    return _DLM(int(2*l), int(2*m1), int(2*m2), alpha, beta, gamma)


cpdef cnp.ndarray[double complex, ndim=1] Rotate(cnp.ndarray[double complex, ndim=1] initial,
                                                 double alpha, double beta, double gamma):
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

    # Determine rank L from array length: len(initial) = 2*l + 1
    cdef int two_l = len(initial) - 1

    # Allocate memory for the rotated state
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(len(initial), dtype=np.complex128)

    # Call the external C function with spin j
    _Rot(two_l,
             <double complex *> cnp.PyArray_DATA(initial),
             alpha, beta, gamma,
             <double complex *> cnp.PyArray_DATA(myOp))

    return myOp

cpdef cnp.ndarray[double complex, ndim=1] create_rho1(double zeta):
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

cpdef cnp.ndarray[double complex, ndim=1] create_rho2(double zeta, double eta):
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

