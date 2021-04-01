import rampy
import numpy as np
from scipy.interpolate import interp1d
from enviornment import Enviornment


def run_ram():
    """replaces the main loop in ram"""
    env = Enviornment()
    env.read_in("ram.in")

    # FORTRAN preallocation parameters
    mz = len(env.zaxis)
    mr = len(env.raxis)
    mp = 2 * env.ram_info["num_pade"]

    pline = []
    pgrid = []

    (nz,nmp,ns,ndr,ndz,iz,nzplt,lz,ib,ir,deir,dr,dz,omega,rmax,
     c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,
     ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2) = rampy.setup(mr,mz,mp)

    pline.append(rampy.outln(ir,deir,r,f3,u))
    pgrid.append(rampy.outgr(ndz,nzplt,lz,r,f3,u))

    if len(env.profiles) > 2:
        rp = env.profiles[1]['rp']
    else:
        # range independent enviornment
        rp = rmax * 2

    # interpolate enviornment to RAM grid
    prof_i = 0
    profs = env.populate_quantities(prof_i)
    (cw_p, cb_p, rhob_p, attn_p, ksqw_p, ksqb_p, alpw_p, alpb_p) = profs

    mdr=0
    for r in env.raxis:
        (iz,ib) = rampy.updat(nz,nmp,dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,
                              rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,
                              s1,s2,s3,pd1,pd2)

        """
        # update changing profiles here
        if r > rp:
            # interpolate enviornment to RAM grid
            prof_i += 1
            profs = env.populate_quantities(prof_i)
            (cw_p, cb_p, rhob_p, attn_p, ksqw_p, ksqb_p, alpw_p, alpb_p) = profs

            if len(env.profiles) > prof_i + 1:
                rp = env.profiles[prof_i + 1]['rp']
            else:
                # range independent enviornment
                rp = rmax * 2

            (f1,f2,f3,r1,r2,r3,s1,s2,s3) = rampy.matrc(nz,nmp,iz,iz,dz,
                                                       k0,rhob,alpw,alpb,ksq,
                                                       ksqw,ksqb,pd1,pd2)

        """
        u, v = rampy.solve(nz,nmp,iz,r1,r2,r3,s1,s2,s3)

        pline.append(rampy.outln(ir,deir,r,f3,u))

        mdr=mdr+1
        if mdr == ndr:
            mdr=0
            pgrid.append(rampy.outgr(ndz,nzplt,lz,r,f3,u))


    pline = np.array(pline)
    pgrid = np.array(pgrid)

    return pline, pgrid


def read_in(file_name):
    """read RAM imput from file_name"""
    with open(file_name, 'r') as f:
        title = f.readline()  # title
        [facous, z_src, z_rcr] = np.array(f.readline().split()).astype(float)
        omega = 2 * np.pi * facous
        [rmax, dr, ndr] = np.array(f.readline().split()).astype(float)
        ndr = int(ndr)
        [zmax, dz, ndz, zmplt] = np.array(f.readline().split()).astype(float)
        ndr = int(ndz)
        [c0, num_pade, num_stab, rs] = np.array(f.readline().split()).astype(float)
        num_pade = int(num_pade)
        num_stab = int(num_stab)

        # read bottom position
        bathy = [np.array(f.readline().split()).astype(float)]
        while bathy[-1][0] >= 0:
            bathy.append(np.array(f.readline().split()).astype(float))
        bathy[-1] = np.array([2 * rmax, bathy[-2][-1]])
        bathy = np.array(bathy)

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

            profile.append({"rp":rp,
                            "cw":cw,
                            "cb":cb,
                            "rhob":rhob,
                            "attn":attn})

            rp = f.readline().split()
            if rp:
                rp = float(rp[0])
            else:
                break

        return profile



def read_line(file_name, num_bytes=8, is_cmpx=True):
    """read RAM line outpult from file_name"""
    tl = np.loadtxt(file_name)
    r = tl[:, 0]
    if not is_cmpx:
        tl = tl[:, 1]
    else:
        tl = -20 * np.log10(np.abs(tl[:, 1] + 1j * tl[:, 2]) + np.spacing(1))
    return r, tl

def read_grid(file_name, num_bytes=8, is_cmpx=True):
    """read RAM grid outpult from file_name"""

    if num_bytes == 8:
        cdt=np.complex128
        fdt=np.float64
        idt=np.int64
    elif num_bytes == 4:
        cdt=np.complex64
        fdt=np.float32
        idt=np.int32
    else:
        raise(ValueError("num_bytes must be 4 or 8"))

    with open(file_name, 'rb') as f:
        _ = np.frombuffer(f.read(4), dtype='int32' )
        # fortran adds a field at every write
        freq = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        z_src = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        z_rcr = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        r_max = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        dr = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        ndr = np.frombuffer(f.read(num_bytes), dtype=idt)[0]
        z_max = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        dz = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        ndz = np.frombuffer(f.read(num_bytes), dtype=idt)[0]
        zmplt = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        c0 = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        num_pade = np.frombuffer(f.read(num_bytes), dtype=idt)[0]
        ns = np.frombuffer(f.read(num_bytes), dtype=idt)[0]
        rs = np.frombuffer(f.read(num_bytes), dtype=fdt)[0]
        num_zplot = np.frombuffer(f.read(num_bytes), dtype=idt)[0]
        _ = np.frombuffer(f.read(4), dtype='int32' )

        zplot = np.arange(num_zplot) * zmplt / (num_zplot-1);
        num_rplot = int(np.floor(r_max / (ndr * dr)))
        rplot = (np.arange(num_rplot) + 1) * dr

        tl_grid = []

        for _ in rplot:
            _ = np.frombuffer(f.read(4), dtype='int32' )
            if is_cmpx:
                TL = np.frombuffer(f.read(2 * num_bytes * num_zplot), dtype=cdt)
                TL = -20 * np.log10(np.abs(TL) + np.spacing(1))
            else:
                TL = np.frombuffer(f.read(num_bytes * num_zplot), dtype=fdt)
            _ = np.frombuffer(f.read(4), dtype='int32' )
            tl_grid.append(TL)

        tl_grid = np.array(tl_grid)

    return rplot, zplot, tl_grid

"""
# preallocation parameters required by FORTRAN
mr = 1000
mz = 80000
mp = 30

# dir_r = dir
# np_r = np

(nz,np_r,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,dir_r,dr,
 dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,
 attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
 pd1,pd2) = rampy.setup(mr,mz,mp)

1/0
pd1 = np.zeros_like(pd1)
pd2 = np.zeros_like(pd2)
zs = 50.
rampy.selfs(nz, np_r, ns, iz, zs, dr, dz, k0, rhob,
            alpw, alpb, ksq, ksqw, ksqb, f1, f2, f3, u, v, r1, r2, r3,
            s1, s2, s3, pd1, pd2)

omega = 2 * np.pi * facous
k0 = omega / c0
eta = 1.0 / (40.0 * np.pi * np.log10(np.e))
eps = 1.0e-20

assert(np.sum(ksqw != 0) == int(np.floor(zmax / dz)) + 1)
assert(lz == int(np.floor(zmplt / (ndz * dz))) - 1)

# interpolate profile to create variables used by ram
zaxis = np.arange(int(np.floor(zmax / dz)) + 1) * dz

cp_water = np.interp(zaxis, profile[0][1][:, 0], profile[0][1][:, 1])
cp_bottom = np.interp(zaxis, profile[0][2][:, 0], profile[0][2][:, 1])
rho_bottom = np.interp(zaxis, profile[0][3][:, 0], profile[0][3][:, 1])
attn_bottom = np.interp(zaxis, profile[0][4][:, 0], profile[0][4][:, 1])

ksqw_test = (omega / cp_water) ** 2 - k0 ** 2
ksqb_test = ((omega / cp_bottom) * (1.0 + 1j * eta * attn_bottom)) ** 2 - k0 ** 2
alpw_test = np.sqrt(cp_water / c0)
alpb_test = np.sqrt(rho_bottom * cp_bottom / c0)

ib = 1
mdr = 0
r = dr  # start range step at dr

# receiver interpolation specifications
rcr_i = int(np.floor(z_rcr / dz))
rcr_rem = (z_rcr / dz) % 1
nzplt = int(zmplt / dz - 0.5)

# best guess at start possition for decompostion
iz = int(1.0 + bathy[0, 1] / dz)
iz = max(2, iz)
iz = min(nz, iz)

# keep stabalization on for entire run
if rs < dr:
    rs = 2.0 * rmax


alloc_z = zaxis.size
alloc_np = num_pade
#alloc_z = mz
#alloc_np = mp

ksq = np.zeros(alloc_z, dtype=np.complex128)

f1 = np.zeros(alloc_z, dtype=np.float64)
f2 = np.zeros(alloc_z, dtype=np.float64)
f3 = np.zeros(alloc_z, dtype=np.float64)

u = np.zeros(alloc_z, dtype=np.complex128)
v = np.zeros(alloc_z, dtype=np.complex128)

r1 = np.zeros((alloc_z, alloc_np), dtype=np.complex128)
r2 = np.zeros((alloc_z, alloc_np), dtype=np.complex128)
r3 = np.zeros((alloc_z, alloc_np), dtype=np.complex128)

s1 = np.zeros((alloc_z, alloc_np), dtype=np.complex128)
s2 = np.zeros((alloc_z, alloc_np), dtype=np.complex128)
s3 = np.zeros((alloc_z, alloc_np), dtype=np.complex128)

pd1 = np.zeros(alloc_np, dtype=np.complex128)
pd2 = np.zeros(alloc_np, dtype=np.complex128)

rampy.selfs(zaxis.size, num_pade, num_stab, iz, z_src, dr, dz, k0, rho_bottom,
            alpw_test, alpb_test, ksq, ksqw_test, ksqb_test, f1, f2, f3, u, v, r1, r2, r3,
            s1, s2, s3, pd1, pd2)

#rampy.selfs(alloc_z, alloc_np, num_stab, iz, z_src, dr, dz, k0, rho_bottom,
            #alpw_test, alpb_test, ksq, ksqw_test, ksqb_test, f1, f2, f3, u, v, r1, r2, r3,
            #s1, s2, s3, pd1, pd2)

#ksqw = (omega / cw) ** 2 - k0 ** 2
#ksqb = ((omega / cb) * (1.0 + 1j * eta * attn)) ** 2 - k0 ** 2
#alpw = np.sqrt(cw / c0)
#alpb = np.sqrt(rhob * cb / c0)

#     open(unit=1,status='old',file='ram.in')
#     open(unit=2,status='unknown',file='tl.line')
#     open(unit=3,status='unknown',file='tl.grid',form='unformatted')
#
#     call setup(mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,
#    >   dz,pi,eta,eps,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,
#    >   attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
#    >   pd1,pd2)
#     March the acoustic field out in range.
#
#   1 call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,eta,omega,rmax,c0,k0,r,
#    >   rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,
#    >   r1,r2,r3,s1,s2,s3,pd1,pd2)
#     call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
#     r=r+dr
#     call outpt(mz,mdr,ndr,ndz,nzplt,lz,ir,dir,eps,r,f3,u,tlg)
#     if(r.lt.rmax)go to 1
#
#     close(1)
#     close(2)
#     close(3)
#

"""

