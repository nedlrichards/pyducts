import numpy as np
from math import pi
import matplotlib.pyplot as plt
from modes import kr_modes


plt.ion()

# Recreate Munk profile results from COA
omega = 2 * pi * 50
# munk ssp
zref = 1300.
eps = 0.00737
cmin = 1500.
zhat = lambda z: 2 * (z - zref) / zref
c_munk = lambda z: cmin * (1 + eps * (zhat(z) - 1 + np.exp(-zhat(z))))
zmax = 5e3

cb = 1600.
rhob = 0.

kr_eig = kr_modes(omega, omega / cb, omega / cmin, zmax, c_munk, 2)


1 / 0

D = 100.
fc = 400.

ct = 0.
rhot = 0.

c = 1500.
rho = 1000.

cb = 1800.
rhob = 0.

# compute finite difference estimate of the modal equations
omega = 2 * pi * fc
k_water = omega / c
k_bottom = omega / cb


