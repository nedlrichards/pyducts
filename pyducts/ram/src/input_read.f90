module input_read
    implicit none
    private
    use constants     ,only : i_, pi

contains
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

        read(1,*)
        read(1,*) freq, zs, zr
        read(1,*) rmax, dr, ndr
        read(1,*) zmax, dz, ndz, zmplt
        read(1,*) c0, np, ns, rs

        i=1
        read(1,*) rb(i), zb(i)

        do while rb(i) >= 0.0
            i = i + 1
            read(1,*) rb(i), zb(i)
        end do

        rb(i) = 2.0 * rmax
        zb(i) = zb(i - 1)

        ib = 1
        r = dr
        omega = 2.0 * pi * freq
        ri = 1.0 + zr / dz
        ir = int(ri, 8)
        dir = ri - float(ir)
        k0 = omega / c0
        nz = int(zmax / dz - 0.5, 8)
        nzplt = int(zmplt / dz - 0.5, 8)
        z = zb(1)
        iz = int(1.0 + z / dz, 8)
        iz = max(2, iz)
        iz = min(nz, iz)

        if(rs < dr) rs = 2.0 * rmax

        if(nz + 2 > mz) then
        write(*,*)'   Need to increase parameter mz to ',nz+2
        stop
        end if

        if(np > mp) then
        write(*,*)'   Need to increase parameter mp to ',np
        stop
        end if

        if(i > mr) then
        write(*,*)'   Need to increase parameter mr to ',i
        stop
        end if

        do j = 1, mp
            r3(1, j) = 0.0
            r1(nz + 2, j) = 0.0
        continue

        do i = 1, nz+2
            u(i) = 0.0
            v(i) = 0.0
        continue

        lz = 0
        do i = ndz, nzplt, ndz
            lz = lz + 1
        continue

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
