import numpy as np
from math import pi, sqrt, copysign
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, newton
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import numba as nb
from scipy.linalg import get_lapack_funcs, get_blas_funcs

class Modes:
    """Compute eigenvalues and vectors for acoustic waveguides"""
    def __init__(self, omega, c_ier, c_bounds, z_bottom, bottom_HS=None):
        """basic setup"""
        self.omega = omega
        self.c_ier = c_ier
        self.c_bounds = c_bounds
        self.z_bottom = z_bottom
        self.bottom_HS = bottom_HS

    def kr_modes(self, max_pow2=2):
        """Compute eigenvalues between k_min and k_max"""

        k_min = self.omega / max(self.c_bounds)
        k_max = self.omega / min(self.c_bounds)

        # bisection search for crossing regions
        # determine max number of crossings
        args = self.omega, k_max, self.z_bottom, self.c_ier, max_pow2, self.bottom_HS
        det, num_cross = _det_estimate(k_min, *args)

        # make array to save bisection results
        regions = np.full((3, 2, num_cross + 1), np.nan)
        regions[0, :, -1] = k_min
        regions[1, :, -1] = det
        regions[2, :, -1] = num_cross

        # evaluate determinant at start of region
        det, num_cross = _det_estimate(k_max, *args)
        regions[0, :, 0] = k_max
        regions[1, :, 0] = det
        regions[2, :, 0] = num_cross

        # bisect and split loop
        _add_bisection(regions, *args)
        sr = _split_regions(regions)
        while sr is not None:
            _add_bisection(sr, *args)
            sr = _split_regions(regions)

        if self.bottom_HS is None:
            return regions
        else:
            # acoustic HS mean changes in number of zero crossings do not
            # occur at eigen_values
            _bound_zeros(regions)
            return regions

        # compute eigen-values with bisection algorithm
        kr_eig = []
        rooter = lambda kr: _det_estimate(kr, *args)[0]

        for i in range(regions.shape[-1] - 1):
            kr_eig.append(brentq(rooter, regions[0, 1, i], regions[0, 0, i + 1]))

        kr_eig = np.array(kr_eig)
        return kr_eig

    def kr_vectors(self, krs, zaxis):
        """compute eigenvectors from eigenvalues"""
        phi_init = np.ones_like(zaxis)
        dz = (zaxis[-1]- zaxis[0]) / (zaxis.size - 1)
        c = self.c_ier(zaxis)

        vecs = []

        for k_eig in krs:
            eig_vec = reverse_iteration(self.omega, k_eig, zaxis, c, phi_init,
                                         bottom_HS=self.bottom_HS)
            vecs.append(eig_vec)
            # get next initial vector by spliting last one at largest peak
            phi_init = eig_vec.copy()
            peaki = find_peaks(np.abs(phi_init))[0]
            max_i = np.argmax(np.abs(phi_init[peaki]))
            phi_init[:peaki[max_i]] *= -1

        vecs = np.array(vecs)
        return vecs

def reverse_iteration(omega, kr, zaxis, c, phi_0, bottom_HS=None):
    """Use reverse iteration to calculate eigenvector from value"""
    dz = (zaxis[-1] - zaxis[0]) / (zaxis.size - 1)
    if bottom_HS is not None:
        numz = zaxis.size - 1
        # assume pressure release top
        k_profile = omega / c[1:]
        phi_0 = np.real(phi_0[1:])
        f = 1
        gamma = np.sqrt(kr ** 2 - (omega / bottom_HS[0]) ** 2)

        g = bottom_HS[1] / (1000. * gamma)
    else:
        # assume pressure release top and bottom
        numz = zaxis.size - 2
        k_profile = omega / c[1: -1]
        phi_0 = np.real(phi_0[1:-1])

    # allocate C matrix as real valued
    d = -2 + dz ** 2 * (k_profile ** 2 - np.real(kr) ** 2) + 0j
    C_band = np.ones((3, numz), dtype=np.complex128)
    C_band[1, :] = d

    if bottom_HS is not None:
        C_band[0, -1] = -2
        # assume a water column density of 1000
        C_band[1, -1] = 2 * dz * 1000. * gamma / bottom_HS[1] + 2 \
                      - dz ** 2 * (k_profile[-1] ** 2 - np.real(kr) ** 2)

    # tridiagonal matrix solver
    gtsv = get_lapack_funcs('gtsv', (d,))
    gbmv = get_blas_funcs('gbmv', (d,))

    phi_out = np.zeros_like(zaxis)
    flag = True
    norm_previous = None

    for i in range(50):
        phi_1 = gtsv(C_band[0, 1:], C_band[1, :], C_band[2, 1:], phi_0)[3]
        norm = np.sqrt(np.trapz(np.abs(phi_1) ** 2) * dz)
        # assume a water column density of 1000
        norm /= np.sqrt(1000.)

        phi_0 = phi_1 / norm
        #print(norm)

        if norm_previous is not None and abs(norm - norm_previous) / norm < 1e-3:
            flag = False
            break
        else:
            norm_previous = norm

    if i > 10:
        print(i)

    if flag:
        import ipdb; ipdb.set_trace()
        pass

    if bottom_HS is None:
        phi_out[1: -1] = np.real(phi_0)
        norm += np.abs(phi_1[-1]) ** 2 / (2 * np.real(gamma) * bottom_HS[1])
        phi_out /= norm
    else:
        phi_out[1:] = np.real(phi_0)
    return phi_out


def _det_estimate(kr, omega, k_max, z_max, c_ier, max_pow2, bottom_HS):
    """combine sturm sequence and richardson extrapolation"""
    # divide integration axis in integer factors of fundimental length
    # add a guard against kr=k_min case
    if bottom_HS is not None and np.abs(kr - omega / bottom_HS[0]) < 1e-10:
        kr += 1e-10

    kz_ref = sqrt(k_max ** 2 - kr ** 2)
    if kz_ref > 0:
        lz = 2 * pi / kz_ref
        dz_fundimental = z_max / max(np.ceil(z_max / lz), 1)
    else:
        dz_fundimental = z_max

    # number of iterations is set by finest sampling
    pow2s = np.arange(max_pow2 + 1)
    dzs = dz_fundimental / (8 * 2 ** pow2s)

    num_steps = np.array(np.round(z_max / dzs), dtype=np.int, ndmin=1)
    k_interval = np.array(num_steps[-1] / num_steps, dtype=np.int)
    s_interval = np.array(num_steps / num_steps[0], dtype=np.int)

    zaxis = np.arange(num_steps[-1] + 1) * dzs[-1]
    c_up = c_ier(zaxis)
    kz2 = (omega / c_up) ** 2 - kr ** 2

    shot = []
    zs = []
    d_phi_D = []
    for k_step, s_step, dz in zip(k_interval, s_interval, dzs):
        # pressure release bottom formulation
        sturm = _sturm_sequence(kz2[::k_step], dz)
        shot.append(sturm[s_step:][::s_step])

        # use second order estimate of first derivative
        dphi_dz = (sturm[-1] - sturm[-2]) / dz - kz2[-1] * sturm[-1] * dz / 2
        d_phi_D.append(dphi_dz)

    shot = np.array(shot)
    d_phi_D = np.array(d_phi_D)

    for i in range(shot.shape[0] - 1):
        shot = shot[1:, :] + np.diff(shot, axis=0) / (4 ** (i + 1) - 1)
        d_phi_D = d_phi_D[1:] + np.diff(d_phi_D) / (4 ** (i + 1) - 1)

    # count number of zero crossings
    sgn = np.sign(shot)
    # assume floating point will take care of 0 edge case
    num_cross = np.sum(np.abs(np.diff(sgn)) > 1)


    # contribution from half space
    phi_D = shot[-1, -1]
    if bottom_HS is None:
        det = phi_D
    else:
        k_bottom = omega / bottom_HS[0]
        gamma = np.sqrt(kr ** 2 - k_bottom ** 2)
        det = phi_D + d_phi_D * bottom_HS[1] / gamma

    return det, num_cross

@nb.jit(nopython=True)
def _sturm_sequence(kz2, dz):
    """Scaled sturm sequence designed for mode finding"""
    y_out = np.empty_like(kz2)
    lam = dz ** 2 * kz2
    y_out[0] = 0
    y_out[1] = 1

    for i in range(2, y_out.size):
        y_out[i] = (2 - lam[i]) * y_out[i - 1] - y_out[i - 2]
        # scale sequence
        w = max(abs(y_out[i]), abs(y_out[i-1]))
        if w > 1e10:
            scale = 1e10 / w
            y_out[i] *= scale
            y_out[i - 1] *= scale
        elif w < 1e-10:
            scale = 1e-10 / w
            y_out[i] *= scale
            y_out[i - 1] *= scale
    return y_out * dz

def _bound_zeros(regions):
    """Work up from zero-ith mode to bound eigenvalues"""
    # trival case of a (0 | 1) eigenvalues
    return regions
    if regions.shape[-1] == 2:
        return regions

    hs_regions = regions.copy()
    hs_regions[:, 0, 0] = regions[:, 0, 0]

    # figure the crossing direction from first mode value
    cross_sign = copysign(1, hs_regions[1, 0, 0])

    for i in np.arange(regions.shape[-1] - 1):
        # check if region has been split by value of determinant
        if copysign(1, regions[1, 1, i]) == -cross_sign:
            hs_regions[:, 1, i] = regions[:, 1, i]
        # check if next region has needed sign
        elif copysign(1, regions[1, 0, i + 1]) == -cross_sign:
            hs_regions[:, 1, i] = regions[1, 0, i + 1]
        # need to find needed sign change between sampled points
        else:
            # finding sign change isn't trival, acutally
            first_guess = regions[:, 0, i + 1]


def _add_bisection(sr, *args):
    # first bisection
    k_test = (sr[0, 1, 0] + sr[0, 0, -1]) / 2
    det, num_cross = _det_estimate(k_test, *args)
    current_index = sr[2, 0, 0]
    # adjust crossing index bases on sign
    cross_index = int(num_cross - current_index)
    if np.isnan(sr[0, 0, cross_index]):
        sr[0, :, cross_index] = k_test
        sr[1, :, cross_index] = det
        sr[2, :, cross_index] = num_cross
    else:
        if sr[0, 0, cross_index] > k_test:
            sr[0, 0, cross_index] = k_test
            sr[1, 0, cross_index] = det
        elif sr[0, 1, cross_index] < k_test:
            sr[0, 1, cross_index] = k_test
            sr[1, 1, cross_index] = det

def _split_regions(sr):
    # are there any regions yet undefined?
    isunset = np.isnan(sr[0, 0, :])
    if isunset.sum() == 0:
        # break condition
        return

    # focus on first region
    roi_start = np.argmax(isunset)
    roi_end = np.argmax(np.bitwise_not(isunset[roi_start: ]))
    sub_region = sr[:, :, roi_start - 1: roi_start + roi_end + 1]
    return sub_region
