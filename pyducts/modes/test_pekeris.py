import numpy as np
from math import pi
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from modes import Modes
from kraken_wrapper import write_env, run_kraken, read_mod

plt.ion()

D = 100.
fc = 100.
zs = 25.
zr = 50.

c = 1500.
rho = 1000.

cb = 2000.
rho_b = 1000.


# Analytical solution
k_water = 2 * pi * fc / c
k_bottom = 2 * pi * fc / cb

kz_water = lambda kr: np.sqrt(k_water ** 2 - kr ** 2)
# gamma is -1j * kz_bottom
gamma = lambda kr: np.sqrt(kr ** 2 - k_bottom ** 2)

rooter = lambda kr: np.tan(kz_water(kr) * D) \
                    + rho_b * kz_water(kr) / (rho * gamma(kr))

num_test = 500
test_kr = np.linspace(k_bottom, k_water, num_test + 2)[1:]
kr_bounds = test_kr[: -1][np.abs(np.diff(np.sign(rooter(test_kr)))) > 1]

kr_eig = []
for kr0, kr1 in zip(kr_bounds[:-1], kr_bounds[1:]):
    x0, r = brentq(rooter, kr0, kr1, full_output=True)
    # make sure we drop tan branch cuts
    if np.abs(rooter(x0)) < 1e-3:
        kr_eig.append(x0)
kr_eig = np.array(kr_eig)

numz = 300
z_extension = 1.5

dz = (D * z_extension) / numz
zaxis = np.arange(numz) * dz
water_i = zaxis <= D
bottom_i = zaxis > D

phi = np.sin(kz_water(kr_eig)[:, None] * zaxis) + 0j
phi[:, bottom_i] = np.sin(kz_water(kr_eig)[:, None] * D) \
                 * np.exp(-gamma(kr_eig)[:, None] * (zaxis[bottom_i] - D))

# Kraken solution
write_env('./envs/pekeris.env', fc, [0, D], [c, c], bottom_HS=[cb, rho_b / 1000.])
run_kraken('./envs/pekeris.env')
phi_krak, k_krak, z_krak = read_mod('./envs/pekeris.env')

c_ier = lambda z: np.full_like(np.array(z, ndmin=1, dtype=np.float64), c)
fd_modes = Modes(2 * pi * fc, c_ier, [c, cb], D, bottom_HS=[cb, rho_b])
kr_modes = fd_modes.kr_modes()
1/0


kr_test = kr_eig[-1]
#kr_test = np.real(k_krak[0])
#kr_test = np.real(kr_modes[0])
f_p = lambda t, y: [y[1], -y[0] * (k_water ** 2 - kr_test ** 2)]
phi_out = solve_ivp(f_p, [0, D], [0, 1], method='BDF')

import modes
dets = []
ncs = []
k_test = np.linspace(k_bottom, k_water, num_test + 2)
for k in k_test:
    det, nc = modes._det_estimate(k, 2  * pi * fc, k_water, D, c_ier, 2, [cb, rho_b / 1000.])
    dets.append(det)
    ncs.append(nc)
dets = np.array(dets)

acc = phi_out.y[0, -1] + phi_out.y[1, -1] * rho_b / (rho * np.sqrt(kr_test ** 2 - k_bottom ** 2))

fig, ax = plt.subplots()
ax.plot(phi_out.y[0], phi_out.t)

ax.set_ylim(105, -2)
1/0

1/0

fig, axes = plt.subplots(1, 2, sharey=True)
axes[0].plot(phi[-1, :], zaxis)
axes[0].plot(phi_krak[-1, :], z_krak)
axes[1].plot(phi[-4, :], zaxis)
axes[1].plot(phi_krak[-4, :], z_krak)
axes[0].plot([-10, 10],[D, D], 'k--')
axes[0].set_ylim(120, -5)
axes[0].set_xlim(-1.2, 1.2)
