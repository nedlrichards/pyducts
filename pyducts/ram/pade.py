import numpy as np
import rampy
import refpy
from roots import cmplx_roots_gen
from scipy.special import binom, factorial
from scipy.linalg import solve

class Enviornment:
    """Basic enviornment for RAM run"""
    def __init__(self, ns=1, np=8):
        """setup common parameters for all ranges"""
        self.facous = facous  # acoustic frequency, Hz
        self.z_src = z_src  # source depth, m
        self.z_rcr = z_rcr  # receiver depth used for tl.line, m
        self.rmax = rmax  # max range of computed pressure, m
        self.dr = dr  # range step, m
        self.zmax = zmax  # max depth of computed pressure, m
        self.dz = dz  # depth step, m
        self.c0 = c0  # reference sound speed, m/s
        self.ns = ns  # number of stablity constraints, {1 | 2}
        self.np = np  # number of pade terms
        self.rs = rs  # max range of stabilty constraints

        # bottom depth in linearly interpolated between bathymetry points
        self.rb rb  # range of bathymetry points
        self.zb = zb  # depth of bathymetry points

        # profiles are treated as discontinuities at each range
        self.r_profiles = r_profiles
        self.c_w = c_w  # water sound speed, m/s
        self.c_b = c_b  # bottom sound speed, m/s
        self.rhob = rhob  # bottom density, g/cc
        self.attn = attn  # bottom attenuation, dB/lambda

#def self_starter(


"""
    def matrix_c(k0, dz):
    """Setup diagonal matrix"""
      a1 = k0 ** 2 / 6.0
      a2 = 2.0 * k0 ** 2 / 3.0
      a3 = k0 ** 2 / 6.0
      cfact = 0.5 / dz ** 2
      dfact = 1.0 / 12.0

      f1 = np.ones(mz)
      f2 = np.ones(mz)
      f3 = np.ones(mz)

      if iz == jz:
        i1 = 2
        i2 = nz + 1
        f1[: iz] = 1.0 / alpw
        f3 = alpw
        ksq=ksqw

c
c     New matrices when iz.eq.jz.
c
      if(iz.eq.jz)then
      i1=2
      i2=nz+1
      do 1 i=1,iz
      f1(i)=1.0/alpw(i)
      f2(i)=1.0
      f3(i)=alpw(i)
      ksq(i)=ksqw(i)
    1 continue
      do 2 i=iz+1,nz+2
      f1(i)=rhob(i)/alpb(i)
      f2(i)=1.0/rhob(i)
      f3(i)=alpb(i)
      ksq(i)=ksqb(i)
    2 continue
      end if
"""

def pade_coeffs(k0, dr, mpr, npr, nsr, ipr):
    """Setup pade coefficents"""
    sig = k0 * dr
    n = 2 * npr

    if ipr == 1:
        nu = 0.0
        alp = 0.0
    else:
        nu = 1.0
        alp = -0.25

    binomial = binom(np.arange(n + 1)[:, None],
                     np.arange(n + 1)[None, :]) * 1.0

    dg = rampy.deriv(n+1, sig, alp, binomial, nu)
    bp = dg[1:]

    ap = np.zeros((n, 2 * n), dtype=np.complex)
    fact = 1
    for i in np.arange(n):
        tout = np.zeros(2 * (i + 1), dtype=np.complex)
        tvec = -binom(i + 1, np.arange(i + 1) + 1.) \
               * factorial(np.arange(i + 1) + 1, exact=True) \
               + 0j
        tvec = tvec[:n] if tvec.size > n else tvec
        tvec *= dg[i::-1]
        tout[1 :: 2] = tvec

        ap[i, :2 * (i + 1)] = tout

    # factorial term
    ap = ap[:, :n]
    flati = np.arange(n - 1) * n + 2 * np.arange(n - 1)

    fi = np.unravel_index(flati[:n // 2], ap.shape)
    ap[fi] = factorial(np.arange(1, n // 2 + 1), exact=True)

    # stability constraints
    if nsr >= 1:
        bp[-1] = -1.0 + 0j
        ap[-1, ::2] = (-3) ** (np.arange(n // 2) + 1) + 0j
        ap[-1, 1::2] = 0 + 0j

    if nsr >= 2:
        bp[-2] = -1.0 + 0j
        ap[-2, ::2] = (-1.5) ** (np.arange(n // 2) + 1) + 0j
        ap[-2, 1::2] = 0 + 0j

    cp = solve(ap, bp)
    # find roots of pade polynomials
    dh1 = np.ones(npr + 1, dtype=np.complex)
    dh1[1:] = cp[0:: 2]
    roots = np.zeros(npr, dtype=np.complex)
    cmplx_roots_gen(roots, dh1, 1, 0)
    pd1 = -1 / roots

    dh1 = np.ones(npr + 1, dtype=np.complex)
    dh1[1:] = cp[1:: 2]
    roots = np.zeros(npr, dtype=np.complex)
    cmplx_roots_gen(roots, dh1, 1, 0)
    pd2 = -1 / roots

    return pd1, pd2

# test outputs
kt = 2 * np.pi
drt = 50.

mp_r = 30  # fortran size estimate, must be larger than np
np_r = 8
ns_r = 1
ip_r = 2

dgr = np.zeros(40, dtype=np.complex)
dh1r = np.zeros(40, dtype=np.complex)

pd1r, pd2r = refpy.epade(mp_r, np_r, ns_r, ip_r, kt, 1500., drt)
pd1, pd2 = pade_coeffs(kt, drt, mp_r, np_r, ns_r, ip_r)

rmax,dr,ndr,dz,ndz,c0,np,ns,rs,omega,ir,dirk0,nz,nzplt,iz,rb,zb,r1,r2,u,v,lz = refpy.setup()
