import numpy as np
import rampy
import refpy
from roots import cmplx_roots_gen
from scipy.special import binom, factorial
from scipy.linalg import solve

# fortran preallocation
mr = 1000
mz = 80000
mp = 30
r1 = np.zeros((mz,mp))
r2 = np.zeros((mz,mp))
r3 = np.zeros((mz,mp))
u = np.zeros(mz)
v = np.zeros(mz)
pd1 = np.zeros(mp, dtype=np.complex128)
pd2 = np.zeros(mp, dtype=np.complex128)

with open("ram.in", 'r') as f:
    _ = f.readline()
    [facous, z_src, z_rcr] = np.array(f.readline().split()).astype(float)
    [rmax, dr, ndr] = np.array(f.readline().split()).astype(float)
    [zmax, dz, ndz, zmplt] = np.array(f.readline().split()).astype(float)
    [c0, npr, ns, rs] = np.array(f.readline().split()).astype(float)

    # read bottom position
    bathy = [np.array(f.readline().split()).astype(float)]
    while bathy[-1][0] >= 0:
        bathy.append(np.array(f.readline().split()).astype(float))
    bottom_position = np.array(bathy[:-1])

    # read profiles one at a time
    profile = []
    rp = 0.0
    while isinstance(rp, float):
        cw = [np.array(f.readline().split()).astype(float)]
        while cw[-1][0] >= 0:
            cw.append(np.array(f.readline().split()).astype(float))
        cw = np.array(cw[:-1])

        cb = [np.array(f.readline().split()).astype(float)]
        while cb[-1][0] >= 0:
            cb.append(np.array(f.readline().split()).astype(float))
        cb = np.array(cb[:-1])

        rhob = [np.array(f.readline().split()).astype(float)]
        while rhob[-1][0] >= 0:
            rhob.append(np.array(f.readline().split()).astype(float))
        rhob = np.array(rhob[:-1])

        attn = [np.array(f.readline().split()).astype(float)]
        while attn[-1][0] >= 0:
            attn.append(np.array(f.readline().split()).astype(float))
        attn = np.array(attn[:-1])

        profile.append([rp, cw, cb, rhob, attn])

        rp = f.readline().split()
        if rp:
            rp = float(rp[0])
        else:
            break

    # constants defined in setup
    eta = 1.0 / (40.0 * np.pi * np.log10(np.e))
    eps = 1.0e-20
    omega = 2.0 * np.pi * facous
    ib = 1
    mdr = 0
    r = dr
    ri = 1.0 + z_rcr / dz
    ir = np.floor(ri)
    dir = ri - ir
    k0 = omega/c0
    nz = zmax / dz - 0.5
    nzplt = zmplt / dz - 0.5
    iz = 1.0 + bottom_position[0, 0] / dz
    iz = max(2, iz)
    iz = min(nz, iz)
    lz = int(np.floor(nzplt / ndz))

    # interpolate profile to create variables used by ram
    zaxis = np.arange(int(nz)) * dz

    if cw.shape[0] == 1:
        cw_grid = np.full_like(zaxis, cw[0, 1])
    if cb.shape[0] == 1:
        cb_grid = np.full_like(zaxis, cb[0, 1])

    ksqw = (omega / cw) ** 2 - k0 ** 2
    ksqb = ((omega / cb) * (1.0 + 1j * eta * attn)) ** 2 - k0 ** 2
    alpw = np.sqrt(cw / c0)
    alpb = np.sqrt(rhob * cb / c0)

    # insert source into pressure field
    dis_src = z_src / dz
    src_i = int(dis_src)
    rem_src = dis_src - src_i
    u[src_i] = (1.0 - rem_src) * np.sqrt(2.0 * np.pi / k0) / (dz * alpw[src_i])
    u[src_i + 1] = rem_src * np.sqrt(2.0 * np.pi / k0) / (dz * alpw[src_i])

    # Divide the delta function by (1-X)**2 to get a smooth rhs.

    pd1[0] = 0.0
    pd2[0] = -1.0

    f1,f2,f3,r1,r2,r3,s1,s2,s3  = rampy.matrc(mz,nz,mp,1,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,pd1,pd2)


"""
    def matrix_c(k0, dz):
    #Setup diagonal matrix
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

#pd1r, pd2r = refpy.epade(mp_r, np_r, ns_r, ip_r, kt, 1500., drt)
#pd1, pd2 = pade_coeffs(kt, drt, mp_r, np_r, ns_r, ip_r)

#rmax,dr,ndr,dz,ndz,c0,np,ns,rs,omega,ir,dirk0,nz,nzplt,iz,rb,zb,r1,r2,u,v,lz = refpy.setup()
