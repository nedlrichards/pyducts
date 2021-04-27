import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve, solve_banded
from math import pi

plt.ion()

# ******************************************************
# PE
# ******************************************************

# *** set up Munk SSP ***

zs = 1000.0 # source depth
f = 50
omega = 2 * pi * f  # frequency in Hertz
eps = 0.00737
c0 = 1500.

d = 5000  # bottom depth
nz = 250  # number of finite-difference points
h = d / nz
h2 = h * h  # mesh spacing
z = np.linspace( 0, d, nz )  # grid coordinates

x = 2 * ( z - 1300 ) / 1300
c = c0 * (1 + eps * (x - 1 + np.exp(-x)))

fig, ax = plt.subplots()
ax.plot(c, -z / 1e3)
ax.set_ylabel( 'Depth (km)' )
ax.set_xlabel( 'Sound Speed (m/s)' )

# ******************************************************
# Gaussian starter
# ******************************************************

# \psi(0,z) = \sqrt{ k_0 } \,
# e^{ -{ k_0^2 \over 2 } ( z - z_s )^2 } (6.100)

k0 = omega / c0
zs = 1000.0
isd = zs / d * nz

nr     = 2000
rmax   = 100000.0
deltar = rmax / nr
r = np.linspace( 0.0, rmax, nr )

psi = np.zeros((nz, nr), dtype=np.complex128)

# we have broadened the Gaussian to make it more
# narrow angle
# the usual formula has fac = 1
fac = 10
psi[:, 0] = np.sqrt(k0 / fac) * np.exp(-( k0 / fac) ** 2 \
          * ((z - zs * np.ones(nz))) ** 2 )

# ******************************************************
# Form marching matrix
# ******************************************************
# SPE: 2ik_0 { \pa \psi \over {\pa r} }
#          + { \pa^2 \psi \over {\pa z^2} }
#          + k_0^2 ( n^2 - 1 ) \psi = 0  (6.8)

n = c0 / c

E = np.diag(np.ones(nz - 1), 1) / h2
D1 = np.diag(np.full(nz, -2), 0) / h2
D2 = np.diag(k0 ** 2 * (n ** 2 - 1), 0)

A = D1 + D2 + E + E.T
B = 2j * k0 / deltar * np.identity(nz) - A / 2
C = 2j * k0 / deltar * np.identity(nz) + A / 2

# make a banded version of c
C_diag = 2j * k0 / deltar * np.full(nz, 1. + 0j) \
       + (np.full(nz, -2. + 0j) / h2 + k0 ** 2 * (n ** 2 - 1)) / 2

C_banded = np.array([np.full(nz, 1. + 0j) / (2 * h2),
                     C_diag,
                     np.full(nz, 1. + 0j) / (2 * h2)])

# ******************************************************
# March out in range
# ******************************************************
for ir in range(nr - 1):
   # equivalent to C psi( :, ir+1 ) = B * psi( :, ir );
    #psi[:, ir+1] = solve(C, B @ psi[:, ir])
    psi[:, ir+1] = solve_banded((1, 1), C_banded, B @ psi[:, ir])

# ******************************************************
# PE
# ******************************************************

# --- plot field

# put back Hankel function
hank = np.sqrt(2 / (pi * k0)) * np.exp(1j * ( k0 * r - pi / 4)) \
     * np.diag(1.0 / np.sqrt(r))

tl = 20 * np.log10(np.abs(psi * np.diag(hank)))

fig, ax = plt.subplots()
ax.pcolormesh(r / 1e3, -z / 1e3, tl, vmax=-60, vmin=-100, cmap=plt.cm.gray_r)

ax.set_xlabel( 'Range (km)' )
ax.set_ylabel( 'Depth (km)' )
ax.set_title('PE intensity')
