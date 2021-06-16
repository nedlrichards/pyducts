import numpy as np
import matplotlib.pyplot as plt
from read_mod import read_mod
from math import pi
from scipy.linalg import get_lapack_funcs, get_blas_funcs
from modes import kr_vectors, kr_modes

plt.ion()

phi_krak, kr_krak, z_krak = read_mod('./envs/MunkK')

# Recreate Munk profile results from COA
fc = 50
omega = 2 * pi * fc
# munk ssp
zref = 1300.
eps = 0.00737
cmin = 1500.
zhat = lambda z: 2 * (z - zref) / zref
c_munk = lambda z: cmin * (1 + eps * (zhat(z) - 1 + np.exp(-zhat(z))))
zmax = 5e3
c_bottom = 1600.

dz = cmin / (10 * fc)
numz = int(np.ceil(zmax / dz))
dz = zmax / numz

zaxis = np.arange(numz) * dz
kr = kr_modes(omega, omega / c_bottom, omega / cmin, zmax, c_munk, max_pow2=2)
phi = kr_vectors(omega, kr, zaxis, c_munk)

fig, ax = plt.subplots()
ax.plot(np.real(phi_krak[0]), z_krak)
ax.plot(-np.real(phi[0]), zaxis, '--')

ax.plot(np.real(phi_krak[1]), z_krak)
ax.plot(-np.real(phi[1]), zaxis, '--')
