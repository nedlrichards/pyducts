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
#rhob = 1800.
rhob = 0.

# divide integration axis in integer factors of fundimental length
dz_fundimental = D / max(np.floor(D *  fc / c), 1)

# compute finite difference estimate of the modal equations
k = 2 * pi * fc / c

#@nb.jit(nopython=True)
def sturm_sequence(kr, dz, d0):
    lam = dz ** 2 * (k ** 2 - kr ** 2)
    y = np.empty(d0.size)
    y[0] = 0
    y[1] = 1
    for i in range(2, d0.size):
        y[i] = (2 - lam) * y[i - 1] - y[i - 2]

        # scale sequence
        w = max(abs(y[i]), abs(y[i-1]))
        if w > 1e10:
            scale = 1e10 / w
            y[i] *= scale
            y[i - 1] *= scale
        elif w < 1e-10:
            scale = 1e-10 / w
            y[i] *= scale
            y[i - 1] *= scale

    return y * dz

def det_estimate(kmin, kmax, n):
    """estimate the determinant with finite resolution sturm sequence"""

    dz = dz_fundimental / 2 ** n

    # pressure release bottom formulation
    k_ = np.full(int(np.round(D / dz)) + 1, k)
    d0 = -2 + (dz * k_) ** 2
    rooter = lambda kr: sturm_sequence(kr, dz, d0)[-1]

    return brentq(rooter, kmin, kmax)

def det_estimate_test(kmin, kmax, n):
    """estimate the determinant with finite resolution sturm sequence"""

    dz = dz_fundimental / 2 ** (n - 1)
    # pressure release bottom formulation
    k_ = np.full(int(np.round(D / dz) + 1), k)
    d0 = -2 + (dz * k_) ** 2
    rooter = lambda kr: sturm_sequence(kr, dz, d0)[-1]
    r0 = brentq(rooter, kmin, kmax)

    dz = dz_fundimental / 2 ** n
    # pressure release bottom formulation
    k_ = np.full(int(np.round(D / dz) + 1), k)
    d0 = -2 + (dz * k_) ** 2
    rooter = lambda kr: sturm_sequence(kr, dz, d0)[-1]
    r1 = brentq(rooter, kmin, kmax)

    return r1 + (r1 - r0) / 3

def det_estimate_richardson(kr):
    """estimate the determinant with finite resolution sturm sequence"""
    det = []
    for i in range(1, 4):
        dz = dz_fundimental / 2 ** i
        # pressure release bottom formulation
        k_ = np.full(int(np.round(D / dz)) + 1, k)
        d0 = -2 + (dz * k_) ** 2
        det.append(sturm_sequence(kr, dz, d0)[-1])
    det = np.array(det)

    for i in range(det.size - 1):
        det = det[1:] + np.diff(det) / (4 ** (i + 1) - 1)

    return det

k_ana = np.sqrt(k ** 2 - (pi / D) ** 2)
y_est = []
y_est_test = []
for i in range(2, 10):
    y_est.append(det_estimate(0.14, 0.145, i))
    y_est_test.append(det_estimate_test(0.14, 0.145, i))

brentq(det_estimate_richardson, 0.14, 0.145)
fig, ax = plt.subplots()
ax.semilogy(np.abs(np.array(y_est) - k_ana))
ax.semilogy(np.abs(np.array(y_est_test) - k_ana))
N1 = y_est[1:] + np.diff(y_est) / 3
N2 = N1[1:] + np.diff(N1) / (4 ** 2 - 1)
N3 = N2[1:] + np.diff(N2) / (4 ** 3 - 1)
N4 = N3[1:] + np.diff(N3) / (4 ** 4 - 1)
