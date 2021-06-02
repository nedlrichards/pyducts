import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, newton
from scipy.integrate import solve_ivp
import numba as nb

plt.ion()

D = 100.
fc = 35.

ct = 0.
rhot = 0.

c = 1500.
rho = 1000.

cb = 1800.
rhob = 0.

# divide integration axis in integer factors of fundimental length
dz_fundimental = D / max(np.floor(D *  fc / c), 1)

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
    det = []
    for i in range(1, 4):
        dz = dz_fundimental / 2 ** i
        # pressure release bottom formulation
        k_ = np.full(int(np.round(D / dz)) + 1, k_water)
        d0 = -2 + (dz * k_) ** 2
        y0 = np.empty(d0.size)
        sturm = sturm_sequence(kr, dz, d0, y0)
        det.append(sturm[-1])
    det = np.array(det)

    for i in range(det.size - 1):
        det = det[1:] + np.diff(det) / (4 ** (i + 1) - 1)

    # count number of zero crossings
    sgn = np.sign(sturm)
    # assume floating point will take care of 0 edge case
    num_cross = np.sum(np.abs(np.diff(sgn)) > 1)

    return det, num_cross

dz = dz_fundimental / 2 ** 3
# pressure release bottom formulation
k_ = np.full(int(np.round(D / dz)) + 1, k_water)
d0 = -2 + (dz * k_) ** 2
y0 = np.empty(d0.size)


det, num_cross = det_estimate(k_water)
#f_bot = sturm_sequence(k_bottom, dz, d0, y0)


