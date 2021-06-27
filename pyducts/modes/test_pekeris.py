import numpy as np
from math import pi
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import modes
from kraken_wrapper import write_env, run_kraken, read_mod

plt.ion()

D = 100.
fc = 100.
zs = 25.
zr = 50.

c = 1500.
rho = 1000.

cb = 1600.
rho_b = 1800.
#rho_b = 1000.

# Analytical solution
k_water = 2 * pi * fc / c
k_bottom = 2 * pi * fc / cb
# numerical implimentation
bottom_HS = [cb, rho_b]

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
kr_eig = np.array(kr_eig)[::-1]

numz = 300
z_extension = 1.5

dz = (D * z_extension) / numz
zaxis = np.arange(numz) * dz
water_i = zaxis <= D
bottom_i = zaxis > D

phi_ana = np.sin(kz_water(kr_eig)[:, None] * zaxis) + 0j
phi_ana[:, bottom_i] = np.sin(kz_water(kr_eig)[:, None] * D) \
                 * np.exp(-gamma(kr_eig)[:, None] * (zaxis[bottom_i] - D))

# normalize the modes
rho_vec = np.full_like(zaxis, rho_b)
rho_vec[water_i] = rho
norm = np.sqrt(np.trapz(np.abs(phi_ana) ** 2 / rho_vec, axis=-1) * dz)

norm_ana = np.sqrt(np.trapz(np.abs(phi_ana[:, water_i]) ** 2 / rho, axis=-1) * dz) \
        + np.abs(phi_ana[:, water_i][:, -1]) ** 2 / (2 * gamma(kr_eig) * rho_b)

phi_ana /= norm[:, None]

# Kraken solution
write_env('./envs/pekeris.env', fc, [0, D], [c, c], bottom_HS=[cb, rho_b / rho])
run_kraken('./envs/pekeris.env')
phi_krak, k_krak, z_krak = read_mod('./envs/pekeris.env')

# focus on water-born modes
wb_i = np.abs(np.imag(k_krak)) < 1e-6
k_krak = np.real(k_krak[wb_i]).astype(dtype=np.float64)
phi_krak = np.real(phi_krak[wb_i, :])
# match units of denisty
phi_krak *= np.sqrt(1000.)

# test reverse iteration generation of modes
decimation = 10
dz_modes = c / (decimation * fc)
numz = int(np.ceil(D / dz_modes))
dz_modes = D / numz
z_modes = np.arange(numz + 1) * dz_modes
c_modes = np.full_like(z_modes, c)

phi_init = np.ones_like(z_modes)
# normalize phi0
norm = np.sqrt(np.trapz(np.abs(phi_init) ** 2 / rho, axis=-1) * dz) \
     + np.abs(phi_init[-1]) ** 2 / (2 * gamma(k_water) * rho_b)

phi_init /= norm

phi_modes = []
for i, k in enumerate(k_krak):
    phi = modes.reverse_iteration(2 * pi * fc, k, z_modes,
                                            c_modes, phi_init,
                                            bottom_HS=bottom_HS)
    # Gram-Schmidt orthogonalization
    phi_init -= phi * np.dot(phi, phi_init) / np.dot(phi, phi)

    norm = np.sqrt(np.trapz(np.abs(phi_init) ** 2 / rho, axis=-1) * dz) \
        + np.abs(phi_init[-1]) ** 2 / (2 * gamma(k) * rho_b)

    phi_init /= norm

    # make sign of first peak positive
    peaki = find_peaks(np.abs(phi_init))[0]
    phi_modes.append(np.sign(phi[peaki[0]]) * phi)
    #phi_modes.append(phi)

phi_modes = np.array(phi_modes)

p_i = 0
fig, ax = plt.subplots()
ax.plot(phi_ana[p_i], zaxis)
ax.plot(phi_krak[p_i], z_krak, ':')
ax.plot(phi_modes[p_i], z_modes, '--')

ax.set_ylim(110, -1)

1/0
# mode finding is not worked out yet

c_ier = lambda z: np.full_like(np.array(z, ndmin=1, dtype=np.float64), c)
fd_modes = Modes(2 * pi * fc, c_ier, [c, cb], D, bottom_HS=[cb, rho_b / 1000.])
kr_modes = fd_modes.kr_modes()


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

fig, ax = plt.subplots()
ax.plot(k_test, dets)
ax.plot(kr_modes[0, 0, :], kr_modes[1, 0, :], '.')

1/0

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
