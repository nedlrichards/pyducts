from numpy cimport ndarray, int64_t, float64_t, complex128_t
from numpy import empty, int64, float64, complex128

cdef extern:
    void c_pe_pade(int64_t *mp, int64_t *np, int64_t *ns, int64_t *ip, float64_t *k0, float64_t *dr,
                   complex128_t *pd1, complex128_t *pd2)

def pe_pade(int64_t np, int64_t ns, int64_t ip, float64_t k0, float64_t dr, int64_t mp=20):
    cdef ndarray[complex128_t, mode="c"] pd1 = empty(mp, dtype=complex128)
    cdef ndarray[complex128_t, mode="c"] pd2 = empty(mp, dtype=complex128)
    c_pe_pade(&mp, &np, &ns, &ip, &k0, &dr, &pd1[0], &pd2[0])
    return pd1, pd2
