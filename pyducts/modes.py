import numpy as np
from math import pi, sqrt
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, newton
from scipy.integrate import solve_ivp
import numba as nb

plt.ion()

D = 100.
fc = 400.

ct = 0.
rhot = 0.

c = 1500.
rho = 1000.

cb = 1800.
rhob = 0.


# compute finite difference estimate of the modal equations
k_water = 2 * pi * fc / c
k_bottom = 2 * pi * fc / cb

@nb.jit(nopython=True)
def sturm_sequence(kr, dz, d0, y0):
    """Scaled sturm sequence designed for mode finding"""
    lam = dz ** 2 * (k_water ** 2 - kr ** 2)
    y0[0] = 0
    y0[1] = 1
    for i in range(2, d0.size):
        y0[i] = (2 - lam) * y0[i - 1] - y0[i - 2]

        # scale sequence
        w = max(abs(y0[i]), abs(y0[i-1]))
        if w > 1e10:
            scale = 1e10 / w
            y0[i] *= scale
            y0[i - 1] *= scale
        elif w < 1e-10:
            scale = 1e-10 / w
            y0[i] *= scale
            y0[i - 1] *= scale
    return y0 * dz

def det_estimate(kr):
    """estimate the determinant with finite resolution sturm sequence"""
    # divide integration axis in integer factors of fundimental length
    kz = sqrt(k_water ** 2 - kr ** 2)
    if kz > 0:
        lz = 2 * pi / kz
        dz_fundimental = D / max(np.ceil(D / lz), 1)
    else:
        dz_fundimental = D

    det = []
    zs = []
    for i, exp in enumerate(range(3, 5)):
        dz = dz_fundimental / 2 ** exp
        # pressure release bottom formulation
        k_ = np.full(int(np.round(D / dz)) + 1, k_water)
        zaxis = np.arange(np.round(D / dz) + 1) * dz
        d0 = -2 + (dz * k_) ** 2
        y0 = np.empty(d0.size)
        sturm = sturm_sequence(kr, dz, d0, y0)
        det.append(sturm)
    return det
        det.append(sturm[2 ** i:][::2 ** i])
    det = np.array(det)

    for i in range(det.shape[0] - 1):
        det = det[1:, :] + np.diff(det, axis=0) / (4 ** (i + 1) - 1)

    # count number of zero crossings
    sgn = np.sign(det)
    # assume floating point will take care of 0 edge case
    num_cross = np.sum(np.abs(np.diff(sgn)) > 1)

    return det[0, -1], num_cross

def kr_modes(k_min, k_max):
    """Compute eigenvalues between k_min and k_max"""
    # bisection search for crossing regions
    # determine max number of crossings
    det, num_cross = det_estimate(k_bottom)

    # make array to save bisection results
    regions = np.full((3, 2, num_cross + 1), np.nan)
    regions[0, :, -1] = k_bottom
    regions[1, :, -1] = det
    regions[2, :, -1] = num_cross

    # evaluate determinant at start of region
    det, num_cross = det_estimate(k_water)
    regions[0, :, 0] = k_water
    regions[1, :, 0] = det
    regions[2, :, 0] = num_cross

    # bisect and split loop
    _add_bisection(regions)
    sr = _split_regions(regions)
    while sr is not None:
        _add_bisection(sr)
        sr = _split_regions(regions)

    # compute eigen-values with bisection algorithm
    kr_eig = []
    rooter = lambda kr: det_estimate(kr)[0]
    for i in range(regions.shape[-1] - 1):
        kr_eig.append(brentq(rooter, regions[0, 1, i], regions[0, 0, i + 1]))

    kr_eig = np.array(kr_eig)
    return kr_eig

def _add_bisection(sr):
    # first bisection
    k_test = (sr[0, 1, 0] + sr[0, 0, -1]) / 2
    det, num_cross = det_estimate(k_test)
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
    isunset = np.isnan(regions[0, 0, :])
    if isunset.sum() == 0:
        # break condition
        return

    # focus on first region
    roi_start = np.argmax(isunset)
    roi_end = np.argmax(np.bitwise_not(isunset[roi_start: ]))
    sub_region = regions[:, :, roi_start - 1: roi_start + roi_end + 1]
    return sub_region


