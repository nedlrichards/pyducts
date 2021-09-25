from numpy cimport ndarray
from numpy import empty

cdef extern:
    void c_pe_pade(int *mp, int *np, int *ns, int *ip, double *k0, double *dr,
                 complex *pd1, complex *pd2)

def pe_pade(int np, int ns, int ip, double k0, double dr, int mp=20):
    cdef ndarray[complex, mode="c"] pd1 = empty(mp, dtype=complex)
    cdef ndarray[complex, mode="c"] pd2 = empty(mp, dtype=complex)
    c_pe_pade(&mp, &np, &ns, &ip, &k0, &dr, &pd1[0], &pd2[0])
    return pd1, pd2
