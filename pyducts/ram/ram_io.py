import numpy as np
from math import pi

class RamIn:
    """Read in from ram.in"""

    def __init__(self):
        """Read inputs from ram.on"""
        with open("ram.in", "r") as fid:
            fid.readline()

            freq, z_src, z_rcr = fid.readline().split()
            rmax, dr, range_decimation = fid.readline().split()
            zmax, dz, z_decimation, zmax_plot = fid.readline().split()
            ref_c, num_pade, num_stability, range_stability = fid.readline().split()

            def read_bathy(fid):
                """Read bathymetric data line by line"""
                rb, zb = fid.readline().split()
                last_zb = zb
                while rb > 0.:
                    last_zb = zb
                    yield rb, zb
                yield 2 * rmax, last_zb

            bathy = np.array([read_bathy(fid)])
            # TODO: check for no bathy case

            ib = 1
            r_current = dr
            omega = 2.0 * pi * freq
            receiver_bins = 1.0 + z_rcr / dz
            receiver_i = int(ri)
            receiver_mod = ri - ir
            k0 = omega / ref_c
            num_z = int(zmax / dz - 0.5)
            num_z_plot = int(zmax_plot / dz - 0.5)
            num_z_save = int(num_z_plot / z_decimation)

            # find closest index of the bottom
            z_bottom = zb(1)
            bottom_index = int(1.0 + z_bottom / dz)
            bottom_index = max(2, bottom_index)
            bottom_index = min(num_z, bottom_index)

            if range_stability < dr:
                range_stability = 2.0 * rmax

            #r3 = np.zeros((num_z + 2, num_pade), dtype=np.complex128)
            #r1 = np.zeros((num_z + 2, num_pade), dtype=np.complex128)
            #u = np.zeros(num_z + 2, dtype=np.complex128)
            #v = np.zeros(num_z + 2, dtype=np.complex128)


"""
    subroutine setup(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz, &
                    omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,   &
                    alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
        ! Initialize the parameters, acoustic field, and matrices.
        integer*8  ,intent(in)  :: mr,mz,mp
        integer*8  ,intent(out) :: nz,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir
        real*8     ,intent(out) :: dir,dr,dz,omega,rmax,c0,k0,r,rp,rs

        rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),alpb(mz),ksqw(mz),f1(mz),f2(mz),f3(mz)

        complex*16 ,intent(out) :: u(mz),v(mz),ksq(mz),ksqb(mz),
                                   r1(mz,mp),r2(mz,mp),r3(mz,mp),
                                   s1(mz,mp),s2(mz,mp),s3(mz,mp),
                                   pd1(mp),pd2(mp)

      integer*8               :: i,j
      real*8                  :: z,zs,zmax,ri,freq,zmplt,zr

    end subroutine

    subroutine profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,alpw,alpb,ksqw,ksqb)
        ! Set up the profiles.
        integer*8  ,intent(in)  :: mz,nz
        real*8     ,intent(in)  :: dz,omega,rmax,c0,k0
        real*8     ,intent(out) :: rp,cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),alpb(mz),ksqw(mz)
        complex*16 ,intent(out) :: ksqb(mz)

        integer*8               :: i
        real*8                  :: eta

        eta=0.01832338997198569352181968569348d0

        call zread(mz,nz,dz,cw)
        call zread(mz,nz,dz,cb)
        call zread(mz,nz,dz,rhob)
        call zread(mz,nz,dz,attn)
        rp=2.0*rmax
        read(1,*,end=1)rp

        do i=1,nz+2
            ksqw(i)=(omega/cw(i))**2-k0**2
            ksqb(i)=((omega/cb(i))*(1.0+i_*eta*attn(i)))**2-k0**2
            alpw(i)=sqrt(cw(i)/c0)
            alpb(i)=sqrt(rhob(i)*cb(i)/c0)
        continue
    end subroutine

    subroutine zread(mz,nz,dz,prof)
        ! Profile reader and interpolator.

        integer*8 ,intent(in)  :: mz,nz
        real*8    ,intent(in)  :: dz
        real*8    ,intent(out) :: prof(mz)

        integer*8              :: i,j,k,iold
        real*8                 :: zi,profi

        do i=1,nz+2
            prof(i)=-1.0
        continue

        read(1,*) zi,profi
        prof(1) = profi
        i = int(1.5 + zi / dz, 8)
        prof(i) = profi
        iold=i
        read(1,*) zi,profi

        if(zi.lt.0.0)go to 3
        i=int(1.5+zi/dz, 8)
        if(i.eq.iold)i=i+1
        prof(i)=profi
        iold=i
        go to 2
    3 prof(nz+2)=prof(i)
      i=1
      j=1
    4 i=i+1
      if(prof(i).lt.0.0)go to 4
      if(i-j.eq.1)go to 6
      do 5 k=j+1,i-1
      prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)
    5 continue
        j=i
        if(j.lt.nz+2)go to 4

    end subroutine
"""
