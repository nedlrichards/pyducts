import numpy as np
from math import pi
import matplotlib.pyplot as plt
import modes

plt.ion()

fc = 100.
D = 100.
c = 1500.

decimation = 8
dz = c / (decimation * fc)
numz = np.ceil(D / dz)
zaxis = np.arange(numz) * D / (numz - 1)
dz = (zaxis[-1] - zaxis[0]) / (zaxis.size - 1)

kz = 2 * np.pi / D
k = 2 * pi * fc / c
kr_exact = np.sqrt(k ** 2 - kz ** 2)

def sturm_sequence(kr, kz2, dz):
    """Scaled sturm sequence designed for mode finding"""
    y_out = np.empty_like(kz2)
    d_out = np.empty_like(kz2)
    d2_out = np.empty_like(kz2)
    lam = 2 - dz ** 2 * kz2
    y_out[0] = 0
    y_out[1] = 1
    d_out[0] = 0
    d_out[1] = 0
    d2_out[0] = 0
    d2_out[1] = 0

    for i in range(2, y_out.size):
        y_out[i] = lam[i] * y_out[i - 1] - y_out[i - 2]
        d_out[i] = 2 * kr * dz ** 2 * y_out[i - 1] + lam[i] * d_out[i - 1] - d_out[i - 2]
        d2_out[i] = dz ** 2 * (2 * y_out[i - 1] + 4 * kr * d_out[i - 1]) \
                  + lam[i] * d2_out[i - 1] - d2_out[i - 2]
        # scale sequence
        w = max(abs(y_out[i]), abs(y_out[i-1]))
        if w > 1e10:
            scale = 1e10 / w
            y_out[i] *= scale
            y_out[i - 1] *= scale
            d_out[i] *= scale
            d_out[i - 1] *= scale
        elif w < 1e-10:
            scale = 1e-10 / w
            y_out[i] *= scale
            y_out[i - 1] *= scale
            d_out[i] *= scale
            d_out[i - 1] *= scale
    return y_out, d_out, d2_out

# compute sturm sequence for eigenvalue
det_test, d_det, d2_det = sturm_sequence(kr_exact, np.full_like(zaxis, kz ** 2), dz)

fd_kr = 1e-6
kr = (kr_exact - fd_kr)
kz2 = k ** 2 - kr ** 2

dt_test = sturm_sequence(kr, np.full_like(zaxis, kz2), dz)[0]

kr = (kr_exact + fd_kr)
kz2 = k ** 2 - kr ** 2
dt_test2 = sturm_sequence(kr, np.full_like(zaxis, kz2), dz)[0]

fd_est = (det_test - dt_test) / fd_kr
fd_est2 = (dt_test2 - 2 * det_test + dt_test) / fd_kr ** 2

#fig, ax = plt.subplots()
#ax.plot(zaxis, np.sin(zaxis * kz))
