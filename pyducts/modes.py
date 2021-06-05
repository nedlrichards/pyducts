import numpy as np
from math import pi, sqrt
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, newton
from scipy.integrate import solve_ivp
import numba as nb


# einsum implimentation of diagonal multiply
n = 1000
d0 = np.random.randn(n)
d1 = np.random.randn(n-1)
x = np.random.randn(n)

A_full = np.diag(d0) + np.diag(d1, k=1) + np.diag(d1, k=-1)
A_sparse = np.zeros((3, n))
test = np.zeros((3, n+2))

A_sparse[0, 1:] = d1
A_sparse[1, ] = d0
A_sparse[2, :-1] = d1

ref = A_full @ x
test[:, 1:-1] = A_sparse * x
ref_sparse = test[2, :-2] + test[1, 1:-1] + test[0, 2:]
assert(np.allclose(ref, ref_sparse))
1/0

def kr_modes(omega, k_min, k_max, z_max, c_ier, max_pow2=2):
    """Compute eigenvalues between k_min and k_max"""
    # bisection search for crossing regions
    # determine max number of crossings
    det, num_cross = det_estimate(omega, k_min, k_max, z_max, c_ier, max_pow2)

    # make array to save bisection results
    regions = np.full((3, 2, num_cross + 1), np.nan)
    regions[0, :, -1] = k_min
    regions[1, :, -1] = det
    regions[2, :, -1] = num_cross

    # evaluate determinant at start of region
    det, num_cross = det_estimate(omega, k_max, k_max, z_max, c_ier, max_pow2)
    regions[0, :, 0] = k_max
    regions[1, :, 0] = det
    regions[2, :, 0] = num_cross

    # bisect and split loop
    _add_bisection(regions, omega, k_max, z_max, c_ier, max_pow2)
    sr = _split_regions(regions)
    while sr is not None:
        _add_bisection(sr, omega, k_max, z_max, c_ier, max_pow2)
        sr = _split_regions(regions)

    # compute eigen-values with bisection algorithm
    kr_eig = []
    rooter = lambda kr: det_estimate(omega, kr, k_max, z_max, c_ier, max_pow2)[0]
    for i in range(regions.shape[-1] - 1):
        kr_eig.append(brentq(rooter, regions[0, 1, i], regions[0, 0, i + 1]))

    kr_eig = np.array(kr_eig)
    return kr_eig

def det_estimate(omega, kr, k_max, z_max, c_ier, max_pow2):
    """combine sturm sequence and richardson extrapolation"""
    # divide integration axis in integer factors of fundimental length

    kz_ref = sqrt(k_max ** 2 - kr ** 2)
    if kz_ref > 0:
        lz = 2 * pi / kz_ref
        dz_fundimental = z_max / max(np.ceil(z_max / lz), 1)
    else:
        dz_fundimental = z_max

    # number of iterations is set by finest sampling
    pow2s = np.arange(max_pow2 + 1)
    dzs = dz_fundimental / (8 * 2 ** pow2s)

    num_steps = np.array(np.round(z_max / dzs), dtype=np.int, ndmin=1)
    k_interval = np.array(num_steps[-1] / num_steps, dtype=np.int)
    s_interval = np.array(num_steps / num_steps[0], dtype=np.int)

    zaxis = np.arange(num_steps[-1] + 1) * dzs[-1]
    c_up = c_ier(zaxis)
    kz2 = (omega / c_up) ** 2 - kr ** 2

    det = []
    zs = []
    for k_step, s_step, dz in zip(k_interval, s_interval, dzs):
        # pressure release bottom formulation
        sturm = _sturm_sequence(kz2[::k_step], dz)
        det.append(sturm[s_step:][::s_step])
    det = np.array(det)

    for i in range(det.shape[0] - 1):
        det = det[1:, :] + np.diff(det, axis=0) / (4 ** (i + 1) - 1)

    # count number of zero crossings
    sgn = np.sign(det)
    # assume floating point will take care of 0 edge case
    num_cross = np.sum(np.abs(np.diff(sgn)) > 1)

    return det[-1, -1], num_cross

@nb.jit(nopython=True)
def _sturm_sequence(kz2, dz):
    """Scaled sturm sequence designed for mode finding"""
    y_out = np.empty_like(kz2)
    lam = dz ** 2 * kz2
    y_out[0] = 0
    y_out[1] = 1

    for i in range(2, y_out.size):
        y_out[i] = (2 - lam[i]) * y_out[i - 1] - y_out[i - 2]
        # scale sequence
        w = max(abs(y_out[i]), abs(y_out[i-1]))
        if w > 1e10:
            scale = 1e10 / w
            y_out[i] *= scale
            y_out[i - 1] *= scale
        elif w < 1e-10:
            scale = 1e-10 / w
            y_out[i] *= scale
            y_out[i - 1] *= scale
    return y_out * dz


def _add_bisection(sr, omega, k_max, z_max, c_ier, max_pow2):
    # first bisection
    k_test = (sr[0, 1, 0] + sr[0, 0, -1]) / 2
    det, num_cross = det_estimate(omega, k_test, k_max, z_max, c_ier, max_pow2)
    cross_index = int(num_cross - sr[2, 0, 0])
    if np.isnan(sr[0, 0, cross_index]):
        sr[0, :, cross_index] = k_test
        sr[1, :, cross_index] = det
        sr[2, :, cross_index] = num_cross
    else:
        if sr[0, 0, cross_index] > k_test:
            sr[0, 0, cross_index] = k_test
            sr[1, 0, cross_index] = det
        elif sr[0, 1, cross_index] < k_test:
            sr[0, 1, cross_index] = k_test
            sr[1, 1, cross_index] = det

def _split_regions(sr):
    # are there any regions yet undefined?
    isunset = np.isnan(sr[0, 0, :])
    if isunset.sum() == 0:
        # break condition
        return

    # focus on first region
    roi_start = np.argmax(isunset)
    roi_end = np.argmax(np.bitwise_not(isunset[roi_start: ]))
    sub_region = sr[:, :, roi_start - 1: roi_start + roi_end + 1]
    return sub_region
