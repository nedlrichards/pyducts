module ram_v3
    use pade_coeffs   ,only : pe_pade
    use constants     ,only : pi, i_, e
    implicit none

contains
    subroutine inram(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz,    &
                     omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksq,ksqw,   &
                     ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,zs,zr,zmax,zmplt)

        ! Initialize the parameters, acoustic field, and matrices.
        integer*8  ,intent(in)  :: mr,mz,mp
        integer*8  ,intent(out) :: nz,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir
        real*8     ,intent(out) :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt

        real*8   ,allocatable    ,dimension(:)    :: rp_tmp
        real*8   ,allocatable    ,dimension(:,:)  :: cw,cb,attn,rho_tmp

        real*8     ,intent(out) ,allocatable ,dimension(:)  :: rp
        real*8     ,intent(out) ,allocatable ,dimension(:,:):: rhob,alpw,alpb,ksqw
        complex*16 ,intent(out) ,allocatable ,dimension(:,:):: ksqb

        complex*16 ,intent(out) ,dimension(mz)    :: u,v,ksq
        complex*16 ,intent(out) ,dimension(mz,mp) :: r1,r2,r3,s1,s2,s3

        real*8     ,intent(out) ,dimension(mr)    :: rb,zb
        real*8     ,intent(out) ,dimension(mz)    :: f1,f2,f3

        integer*8               :: i,j,max_nprof,nprof,iostat
        real*8                  :: z,ri,freq,eta,r_read

        eta=1 / (40. * pi * log10(e))

        read(1,*)
        read(1,*)freq,zs,zr
        read(1,*)rmax,dr,ndr
        read(1,*)zmax,dz,ndz,zmplt
        read(1,*)c0,np,ns,rs

        i=1
        read(1,*)rb(i),zb(i)
        do while (rb(i) >= 0.0)
            i=i+1
            read(1,*)rb(i),zb(i)
        end do
        rb(i)=2.0*rmax
        zb(i)=zb(i-1)

        max_nprof=int(rmax / dr, 8)
        omega=2.0*pi*freq
        ri=1.0+zr/dz
        ir=int(ri,8)
        dir=ri-float(ir)
        k0=omega/c0
        nz=int(zmax/dz-0.5,8)
        nzplt=int(zmplt/dz-0.5,8)
        z=zb(1)
        iz=int(1.0+z/dz,8)
        iz=max(2,iz)
        iz=min(nz,iz)
        if(rs.lt.dr)rs=2.0*rmax

        if(nz+2.gt.mz)then
        write(*,*)'   Need to increase parameter mz to ',nz+2
        stop
        end if
        if(np.gt.mp)then
        write(*,*)'   Need to increase parameter mp to ',np
        stop
        end if
        if(i.gt.mr)then
        write(*,*)'   Need to increase parameter mr to ',i
        stop
        end if

        lz=0
        do i=ndz,nzplt,ndz
            lz=lz+1
        end do

        allocate(rp_tmp(max_nprof))
        allocate(cw(mz, max_nprof))
        allocate(cb(mz, max_nprof))
        allocate(attn(mz, max_nprof))
        allocate(rho_tmp(mz, max_nprof))

        iostat = 0
        ! set range of 0th profile before loop
        r_read = 0
        nprof = 1

        do while (iostat == 0)
            rp_tmp(nprof) = r_read
            ! interpolated profiles read from file

            call zread(nz,dz,cw(:, nprof))
            call zread(nz,dz,cb(:, nprof))
            call zread(nz,dz,rho_tmp(:, nprof))
            call zread(nz,dz,attn(:, nprof))

            read(1,*,iostat=iostat)r_read
            nprof = nprof + 1
        end do

        rp_tmp(nprof) = 2.0 * rmax

        allocate(rp(nprof))
        allocate(ksqw(nz+2, nprof-1))
        allocate(ksqb(nz+2, nprof-1))
        allocate(alpw(nz+2, nprof-1))
        allocate(alpb(nz+2, nprof-1))
        allocate(rhob(nz+2, nprof-1))

        rp = rp_tmp(1:nprof)
        ! computed quantities from interpolated profiles
        ksqw=(omega / cw(1:nz+2, 1:nprof-1)) ** 2 - k0 ** 2
        ksqb=((omega / cb(1:nz+2, 1:nprof-1))                                  &
            * (1.0 + i_ * eta * attn(1:nz+2, 1:nprof-1))) ** 2 - k0 ** 2
        alpw=sqrt(cw(1:nz+2, 1:nprof-1) / c0)
        alpb=sqrt(rho_tmp(1:nz+2, 1:nprof-1) * cb(1:nz+2, 1:nprof-1) / c0)
        rhob=rho_tmp(1:nz+2, 1:nprof-1)

    end subroutine

    subroutine zread(nz,dz,prof)
        ! Profile reader and interpolator.
        integer*8 ,intent(in)  :: nz
        real*8    ,intent(in)  :: dz
        real*8    ,intent(out) :: prof(nz+2)

        integer*8              :: i,j,k,iold
        real*8                 :: zi,profi


        do i=1, nz+2
            prof(i) = -1.0
        end do

        read(1,*) zi,profi
        prof(1) = profi
        i=int(1.5 + zi / dz,8)
        prof(i) = profi
        iold=i

        read(1,*)zi, profi
        do while (zi >= 0.0)
            i=int(1.5 + zi / dz, 8)
            if (i == iold) i = i + 1
            prof(i) = profi
            iold = i
            read(1,*)zi, profi
        end do

        prof(nz+2)=prof(i)
        i=1
        j=1

        do while (j < nz + 2)
            i=i+1
            if (prof(i) < 0.0) cycle

            if(i-j > 1) then
            ! interpolat between defined values
                do k = j + 1, i - 1
                    prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)
                end do
            end if
            j=i
        end do
    end subroutine

    subroutine outln(mz,ir,dir,r,f3,u)
        ! Output transmission loss line
        integer*8  ,intent(in) :: mz,ir
        real*8     ,intent(in) :: dir,r,f3(mz)
        complex*16 ,intent(in) :: u(mz)

        real*8                 :: eps
        complex*16             :: ur,pout

        eps=1e-20

        ur=(1.0-dir)*f3(ir)*u(ir)+dir*f3(ir+1)*u(ir+1)
        pout=ur/sqrt(r+eps)
        write(2,*)r,real(pout),imag(pout)
    end subroutine

    subroutine outgr(mz,ndz,nzplt,lz,r,f3,u)
        ! Output transmission loss.
        integer*8  ,intent(in) :: mz,ndz,nzplt,lz
        real*8     ,intent(in) :: r,f3(mz)
        complex*16 ,intent(in) :: u(mz)

        integer*8              :: i,j
        real*8                 :: eps
        complex*16             :: ur,poutg(mz)

        eps=1e-20

        j=0
        do i=ndz,nzplt,ndz
            ur=u(i)*f3(i)
            j=j+1
            poutg(j)=ur/sqrt(r+eps)
        end do

        write(3)(poutg(j),j=1,lz)

    end subroutine
end module
