import numpy as np
from math import pi
from scipy.linalg import sqrtm, expm

def solve_pe(psi_start, freq, c, c0, rho, dz, dr, numr):
    """Solve PE by stepping intial pressure field forward in range"""
    k0 = 2 * pi * freq / c0

    n = c0 / c
    n1 = np.hstack([n[0], n[:-1]])
    n2 = np.hstack([n[1:], n[-1]])

    rho1 = np.hstack([rho[0], rho[:-1]])
    rho2 = np.hstack([rho[1:], rho[-1]])

    # rho ratios
    r2_rsum = rho2 / (rho1 + rho2)
    r1_r2 = rho1 / rho2

    # off diagonal term
    nz = psi_start.size

    # direct implimentation of PE operator, Eq. (6.169)
    diff = np.diag(np.full(nz - 1, 1 + 0j), -1) \
         + np.diag(r2_rsum) \
         + np.diag(r1_r2[1: ], 1)
    diff *= 2 * r2_rsum / (k0 * dz) ** 2

    # density modified index of refraction, eta. Eq. (6.170)
    eta = r2_rsum * (n1 ** 2 + r1_r2 * n2 ** 2)

    di = np.diag_indices_from(diff)

    oper = sqrtm(diff + np.diag(eta, 0))
    oper[di] -= 1
    oper = expm(1j * k0 * dr * oper)

    # March out in range
    psi = [psi_start]
    for ir in range(numr - 1):
        psi.append(oper @ psi[ir])

    return np.array(psi)

def solve_pe_energy(psi_start, freq, c, c0, rho, dz, dr, numr):
    """Solve PE by stepping intial pressure field forward in range"""
    k0 = 2 * pi * freq / c0
    k = 2 * pi * freq / c

    # energy conservation quanity
    #alpha = np.sqrt(c * rho)
    alpha = np.ones_like(c)

    # index of refraction term
    nt = k ** 2 / k0 ** 2

    kappa_n1 = np.hstack([nt[0], nt[:-1]])
    kappa_0 = nt
    kappa_p1 = np.hstack([nt[1:], nt[-1]])

    gal_1 = np.diag(kappa_n1[1:] + kappa_0[1:], k=-1) \
          + np.diag(kappa_n1 + 6 * kappa_0 + kappa_p1, k=0) \
          + np.diag(kappa_0[:-1] + kappa_p1[:-1], k=1)
    gal_1 /= 12

    #gal_1 = np.diag(kappa_0, k=0)

    # differential term
    t1 = rho / alpha
    t2 = 1 / rho
    t3 = alpha

    gal_2 = np.diag(t3[:-1] * t1[1:], k=-1) \
          - np.diag(t3 * t1, k=0) \
          + np.diag(t3[1:] * t1[:-1], k=1)

    # finite difference estimate of derivative
    t2_n1 = np.hstack([t2[0], t2[:-1]])
    t2_0 = t2
    t2_p1 = np.hstack([t2[1:], t2[-1]])

    gal_2 *= np.diag(t2_n1[:-1] + t2_0[1:], k=-1) \
           + np.diag(t2_n1 + 2 * t2_0 + t2_p1, k=0) \
           + np.diag(t2_0[:-1] + t2_p1[1:], k=1)
    gal_2 /= (2 * (k0 * dz) ** 2)

    di = np.diag_indices_from(gal_1)

    oper = sqrtm(gal_1 + gal_2)
    oper[di] -= 1
    oper = expm(1j * k0 * dr * oper)

    # March out in range
    field = [psi_start / alpha]
    for ir in range(numr - 1):
        field.append(oper @ field[ir])

    field = np.array(field) * alpha

    return field
