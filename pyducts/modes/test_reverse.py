import numpy as np
import matplotlib.pyplot as plt
from read_mod import read_mod
from math import pi
from scipy.linalg import get_lapack_funcs, get_blas_funcs
from modes import kr_vectors, kr_modes

plt.ion()

phi, kr, z = read_mod('./envs/MunkK')

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
vectors = kr_vectors(omega, kr, zaxis, c_munk)

1/0



c_profile = c_munk(zaxis)

test_num = 55

k_profile = omega / c_profile

# allocate C matrix as complex valued
d = -2 + dz ** 2 * (k_profile ** 2 - kr[test_num] ** 2) + 0j
e = np.ones(numz-1, dtype=np.complex128)
#C_mult = np.array([d, np.ones_like(d)])
C_mult = np.array([np.ones_like(d), d + 1j, np.ones_like(d)])

# tridiagonal matrix solver
gtsv = get_lapack_funcs('gtsv', (d,))
gbmv = get_blas_funcs('gbmv', (d,))
hbmv = get_blas_funcs('hbmv', (d,))

phi_0 = np.ones(k_profile.size)
phi_1 = gtsv(e, d, e, phi_0)[3]

phi_0 /= np.sqrt(np.trapz(np.abs(phi_1) ** 2) * dz)
error = gbmv(numz, numz, 1, 1, 1., C_mult, phi_0)
error_test = hbmv(1, 1., C_mult[:2, :], phi_0)

1/0

normals=[]
phis = []
for i in range(3):
    phi_test = solve(C, phi_test)
    normals.append(np.sqrt(np.trapz(np.abs(phi_test) ** 2)))
    phi_test /= np.sqrt(np.trapz(np.abs(phi_test) ** 2) * dz)
    phis.append(phi_test)

fig, ax = plt.subplots()
ax.plot(-np.real(phis[-1]), zaxis)
ax.plot(np.real(phi[test_num]), z)
