module ram_v4
    use pade_coeffs   ,only : pe_pade
    use constants     ,only : pi, i_
    use ram_io ,only      : outln, outgr, storeln, storegr
    implicit none

contains
    subroutine main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,  &
                    dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,  &
                    ksqb,zs,zr,zmax,zmplt,pline,pgrid,tofile)

        integer*8  :: mr, mz, mp, nprof
        integer*8  :: nr,nz,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,mdr,prof_i,i,j
        real*8     :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt
        logical    :: tofile

        complex*16 ,dimension (mp)    :: pd1,pd2
        real*8     ,dimension (mr)    :: rb, zb
        real*8     ,dimension (mz)    :: f1,f2,f3
        complex*16 ,dimension (mz)    :: ksq,u,v
        complex*16 ,dimension (mz,mp) :: r1,r2,r3,s1,s2,s3

        real*8     ,dimension (nprof)          :: rp
        real*8     ,dimension (nz+2,nprof-1)   :: rhob,alpw,alpb,ksqw
        complex*16 ,dimension (nz+2,nprof-1)   :: ksqb

        complex*16 ,dimension (nr)    ,intent(out) :: pline
        complex*16 ,dimension (lz,nr) ,intent(out) :: pgrid

        do j=1,mp
            r3(1,j)=0.0
            r1(nz+2,j)=0.0
        end do

        do i=1,nz+2
            u(i)=0.0
            v(i)=0.0
        end do

        call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,alpb,ksq,           &
                ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)

        i = 1
        j=1
        r = i * dr

        if (tofile) then
            open(unit=2,status='unknown',file='tl.line')
            open(unit=3,status='unknown',file='tl.grid',form='unformatted')
            write(3) omega/(2*pi),zs,zr,rmax,dr,ndr,zmax,dz,ndz,zmplt,c0,np,ns,rs,lz
        end if

        if (tofile) then
            call outln(mz,ir,dir,r,f3,u)
            call outgr(mz,ndz,nzplt,lz,r,f3,u)
        end if

        call storeln(mz,ir,dir,r,f3,u,pline(i))
        call storegr(mz,ndz,nzplt,lz,r,f3,u,pgrid(:,j))
        j = j + 1

        !! The propagation matrices.
        call pe_pade(mp,np,ns,1_8,k0,dr,pd1,pd2)
        call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob(:,1),alpw(:,1),alpb(:,1),      &
                ksq,ksqw(:,1),ksqb(:,1),f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)

        ! March the acoustic field out in range.
        mdr=0
        prof_i = 1
        ib=1
        do i = 2, nr
            call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,rmax,k0,              &
                    r,rp,rs,rb,zb,rhob,alpw,alpb,ksq,                         &
                    ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,prof_i)

            call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
            r = i * dr

            if (tofile) then
                call outln(mz,ir,dir,r,f3,u)
            end if

            call storeln(mz,ir,dir,r,f3,u,pline(i))

            mdr=mdr+1
            if(mdr == ndr)then
                mdr=0

                if (tofile) then
                    call outgr(mz,ndz,nzplt,lz,r,f3,u)
                end if

                call storegr(mz,ndz,nzplt,lz,r,f3,u,pgrid(:,j))
                j = j + 1

            end if
        end do

        if (tofile) then
            close(2)
            close(3)
        end if

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

    subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,ksqw,          &
                    ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
        ! The tridiagonal matrices.

        integer*8     ,intent(in)    :: mz,nz,mp,np,iz,jz
        real*8        ,intent(in)    :: dz,k0,rhob(:),alpw(:),alpb(:),ksqw(:)
        real*8        ,intent(out)   :: f1(mz),f2(mz),f3(mz)
        complex*16    ,intent(in)    :: ksqb(:),pd1(mp),pd2(mp)
        complex*16    ,intent(out)   :: ksq(:), r1(mz,mp),r2(mz,mp),           &
                                        r3(mz,mp),s1(mz,mp),                   &
                                        s2(mz,mp),s3(mz,mp)

        integer*8                    :: i,j,i1,i2
        real*8                       :: a1,a2,a3,c1,c2,c3,                     &
                                        cfact,dfact
        complex*16                   :: d1,d2,d3,rfact

        a1=k0**2/6.0
        a2=2.0*k0**2/3.0
        a3=k0**2/6.0
        cfact=0.5/dz**2
        dfact=1.0/12.0

        ! New matrices when iz.eq.jz.

        if(iz == jz)then
            i1=2
            i2=nz+1
            do i=1,iz
                f1(i)=1.0/alpw(i)
                f2(i)=1.0
                f3(i)=alpw(i)
                ksq(i)=ksqw(i)
            end do
            do i=iz+1,nz+2
                f1(i)=rhob(i)/alpb(i)
                f2(i)=1.0/rhob(i)
                f3(i)=alpb(i)
                ksq(i)=ksqb(i)
            end do
        end if

        ! Updated matrices when iz.ne.jz.

        if(iz > jz)then
            i1=jz
            i2=iz+1
            do i=jz+1,iz
                f1(i)=1.0/alpw(i)
                f2(i)=1.0
                f3(i)=alpw(i)
                ksq(i)=ksqw(i)
            end do
        end if

        if(iz < jz)then
            i1=iz
            i2=jz+1
            do i=iz+1,jz
                f1(i)=rhob(i)/alpb(i)
                f2(i)=1.0/rhob(i)
                f3(i)=alpb(i)
                ksq(i)=ksqb(i)
            end do
        end if

        do i=i1,i2

        ! Discretization by Galerkins method.

            c1=cfact*f1(i)*(f2(i-1)+f2(i))*f3(i-1)
            c2=-cfact*f1(i)*(f2(i-1)+2.0*f2(i)+f2(i+1))*f3(i)
            c3=cfact*f1(i)*(f2(i)+f2(i+1))*f3(i+1)
            d1=c1+dfact*(ksq(i-1)+ksq(i))
            d2=c2+dfact*(ksq(i-1)+6.0*ksq(i)+ksq(i+1))
            d3=c3+dfact*(ksq(i)+ksq(i+1))

            do j=1,np
                r1(i,j)=a1+pd2(j)*d1
                r2(i,j)=a2+pd2(j)*d2
                r3(i,j)=a3+pd2(j)*d3
                s1(i,j)=a1+pd1(j)*d1
                s2(i,j)=a2+pd1(j)*d2
                s3(i,j)=a3+pd1(j)*d3
            end do
        end do

        ! The matrix decomposition.
        do j=1,np
            do i=i1,iz
                rfact=1.0/(r2(i,j)-r1(i,j)*r3(i-1,j))
                r1(i,j)=r1(i,j)*rfact
                r3(i,j)=r3(i,j)*rfact
                s1(i,j)=s1(i,j)*rfact
                s2(i,j)=s2(i,j)*rfact
                s3(i,j)=s3(i,j)*rfact
            end do

            do i=i2,iz+2,-1
                rfact=1.0/(r2(i,j)-r3(i,j)*r1(i+1,j))
                r1(i,j)=r1(i,j)*rfact
                r3(i,j)=r3(i,j)*rfact
                s1(i,j)=s1(i,j)*rfact
                s2(i,j)=s2(i,j)*rfact
                s3(i,j)=s3(i,j)*rfact
            end do

            r2(iz+1,j)=r2(iz+1,j)-r1(iz+1,j)*r3(iz,j)
            r2(iz+1,j)=r2(iz+1,j)-r3(iz+1,j)*r1(iz+2,j)
            r2(iz+1,j)=1.0/r2(iz+1,j)

        end do
    end subroutine

    subroutine solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
        ! The tridiagonal solver.
        integer*8   ,intent(in)    :: mz,nz,mp,np,iz
        complex*16  ,intent(in)    :: r1(mz,mp),r2(mz,mp),r3(mz,mp),           &
                                      s1(mz,mp),s2(mz,mp),s3(mz,mp)
        complex*16  ,intent(inout) :: u(mz),v(mz)

        real*8                     :: eps
        integer*8                  :: i,j

        eps=1.0e-30

        do j=1,np
            ! The right side.

            do i=2,nz+1
                v(i)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)+eps
            end do

            ! The elimination steps.

            do i=3,iz
                v(i)=v(i)-r1(i,j)*v(i-1)+eps
            end do

            do i=nz,iz+2,-1
                v(i)=v(i)-r3(i,j)*v(i+1)+eps
            end do

            u(iz+1)=(v(iz+1)-r1(iz+1,j)*v(iz)-r3(iz+1,j)*v(iz+2))* &
                    r2(iz+1,j)+eps

            ! The back substitution steps.

            do i=iz,2,-1
                u(i)=v(i)-r3(i,j)*u(i+1)+eps
            end do

            do i=iz+2,nz+1
                u(i)=v(i)-r1(i,j)*u(i-1)+eps
            end do
        end do
    end subroutine

    subroutine updat(mr,mz,nz,mp,np,iz,ib,dr,dz,rmax,k0,              &
                     r,rp,rs,rb,zb,rhob,alpw,alpb,ksq,                         &
                     ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,prof_i)
        ! Matrix updates.
        integer*8 ,intent(in)    :: mr,mz,nz,mp,np
        integer*8 ,intent(inout) :: iz,ib,prof_i
        real*8                   :: dr,dz,rmax,k0,r,rb(mr),           &
                                    zb(mr),f1(mz),f2(mz),f3(mz)
        real*8                   :: rp(:),rs
        real*8                   :: rhob(:,:),alpw(:,:),alpb(:,:),ksqw(:,:)
        complex*16               :: ksq(:),ksqb(:,:),                          &
                                    r1(mz,mp),r2(mz,mp),r3(mz,mp),             &
                                    s1(mz,mp),s2(mz,mp),s3(mz,mp)
        complex*16               :: pd1(mp),pd2(mp)

        integer*8                :: jz,ns
        real*8                   :: z
        logical                  :: isup

        ! Varying bathymetry.
        isup = .false.

        if(r >= rb(ib+1)) ib=ib+1
        jz=iz
        z=zb(ib)+(r+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
        iz=int(1.0+z/dz,8)
        iz=max(2,iz)
        iz=min(nz,iz)

        if(iz /= jz) isup = .true.


        if(r >= rp(prof_i+1)) then
            prof_i=prof_i+1
            isup = .true.
        end if

        ! Turn off the stability constraints.
        if(r >= rs)then
            ns=0
            rs=2.0*rmax
            call pe_pade(mp,np,ns,1_8,k0,dr,pd1,pd2)
            isup = .true.
        end if

        ! Varying profiles.
        if(isup)then
            call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob(:,prof_i),alpw(:,prof_i),  &
                       alpb(:,prof_i),ksq,ksqw(:,prof_i),ksqb(:,prof_i),       &
                       f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
        end if

    end subroutine

    subroutine selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,alpb,             &
                    ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
        ! The self-starter.
        integer*8 ,intent(in)     :: mz,nz,mp,np,ns,iz
        real*8    ,intent(in)     :: zs,dr,dz,k0,rhob(:,:),alpw(:,:),          &
                                     alpb(:,:),ksqw(:,:)
        real*8    ,intent(out)    :: f1(mz),f2(mz),f3(mz)
        complex*16 ,intent(in)    :: ksqb(:,:)
        complex*16 ,intent(inout) :: ksq(:),v(mz),u(mz)
        complex*16 ,intent(out)   :: pd1(mp),pd2(mp),                          &
                                     r1(mz,mp),r2(mz,mp),r3(mz,mp),            &
                                     s1(mz,mp),s2(mz,mp),s3(mz,mp)

        integer*8                 :: is
        real*8                    :: dis,si

        ! Conditions for the delta function.

        si=1.0+zs/dz
        is=int(si,8)
        dis=si-float(is)
        u(is)=(1.0-dis)*sqrt(2.0*pi/k0)/(dz*alpw(is,1))
        u(is+1)=dis*sqrt(2.0*pi/k0)/(dz*alpw(is,1))

        ! Divide the delta function by (1-X)**2 to get a smooth rhs.
        pd1(1)=0.0
        pd2(1)=-1.0
        call matrc(mz,nz,mp,1_8,iz,iz,dz,k0,rhob(:,1),alpw(:,1),alpb(:,1),ksq, &
                   ksqw(:,1),ksqb(:,1),f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
        call solve(mz,nz,mp,1_8,iz,u,v,r1,r2,r3,s1,s2,s3)
        call solve(mz,nz,mp,1_8,iz,u,v,r1,r2,r3,s1,s2,s3)

        ! Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
        call pe_pade(mp,np,ns,2_8,k0,dr,pd1,pd2)
        call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob(:,1),alpw(:,1),alpb(:,1),ksq,  &
                   ksqw(:,1),ksqb(:,1),f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
        call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)

    end subroutine
end module
