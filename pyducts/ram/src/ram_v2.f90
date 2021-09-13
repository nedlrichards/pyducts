program ram
    implicit none

    integer*8 mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir
    complex*16 ksq,ksqb,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2
    real*8 dir,dr,dz,omega,rmax,c0,k0,r,rp,rs
    real*8 rb,zb,cw,cb,rhob,attn,alpw,alpb,ksqw,f1,f2,f3
    ! mr=bathymetry points, mz=depth grid, mp=pade terms.

    parameter (mr=1000,mz=80000,mp=20)
    dimension rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),   &
         alpb(mz),f1(mz),f2(mz),f3(mz),ksq(mz),ksqw(mz),ksqb(mz),u(mz), &
         v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),                 &
         s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)

    open(unit=1,status='old',file='ram.in')
    open(unit=2,status='unknown',file='tl.line')
    open(unit=3,status='unknown',file='tl.grid',form='unformatted')

    call setup(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,    &
         dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb, &
         ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)

    ! March the acoustic field out in range.
    mdr=0
    do while (r < rmax)
        call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,omega,rmax,c0,k0,r, rp,rs,rb, &
                   zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,     &
                   r1,r2,r3,s1,s2,s3,pd1,pd2)
        call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
        r=r+dr

        call outln(mz,ir,dir,r,f3,u)

        mdr=mdr+1
        if(mdr == ndr)then
            mdr=0
            call outgr(mz,ndz,nzplt,lz,r,f3,u)
        end if
    end do

    close(1)
    close(2)
    close(3)
    stop
end program

subroutine setup(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz,    &
                 omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb, &
                 ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
    ! Initialize the parameters, acoustic field, and matrices.

    use pade_coeffs ,only : pe_pade
    use constants     ,only : pi
    implicit none

    integer*8  ,intent(in)  :: mr,mz,mp
    integer*8  ,intent(out) :: nz,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir
    real*8     ,intent(out) :: dir,dr,dz,omega,rmax,c0,k0,r,rp,rs,rb(mr), &
                               zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),    &
                               alpw(mz),alpb(mz),ksqw(mz),f1(mz),f2(mz),f3(mz)

    complex*16 ,intent(out) :: u(mz),v(mz),ksq(mz),ksqb(mz),   &
                               r1(mz,mp),r2(mz,mp),r3(mz,mp),  &
                               s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)

    integer*8               :: i,j
    real*8                  :: eps,z,zs,zmax,ri,freq,zmplt,zr


    eps=1.0e-20

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

    ib=1
    r=dr
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

    do j=1,mp
        r3(1,j)=0.0
        r1(nz+2,j)=0.0
    end do

    do i=1,nz+2
        u(i)=0.0
        v(i)=0.0
    end do

    lz=0
    do i=ndz,nzplt,ndz
        lz=lz+1
    end do
    write(3) freq,zs,zr,rmax,dr,ndr,zmax,dz,ndz,zmplt,c0,np,ns,rs,lz

    !  The initial profiles and starting field.

    call profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,     &
               alpw,alpb,ksqw,ksqb)
    call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,alpb,ksq, &
               ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
    call outln(mz,ir,dir,r,f3,u)
    call outgr(mz,ndz,nzplt,lz,r,f3,u)

    !! The propagation matrices.
    call pe_pade(mp,np,ns,1_8,k0,dr,pd1,pd2)
    call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,           &
               ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
end subroutine

subroutine profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob, &
                 attn,alpw,alpb,ksqw,ksqb)
    ! Set up the profiles.
    use constants ,only : i_
    implicit none

    integer*8  ,intent(in)  :: mz,nz
    real*8     ,intent(in)  :: dz,omega,rmax,c0,k0
    real*8     ,intent(out) :: rp,cw(mz),cb(mz),rhob(mz),attn(mz), &
                               alpw(mz),alpb(mz),ksqw(mz)
    complex*16 ,intent(out) :: ksqb(mz)

    integer*8               :: i, iostat
    real*8                  :: eta

    eta=0.01832338997198569352181968569348d0

    call zread(mz,nz,dz,cw)
    call zread(mz,nz,dz,cb)
    call zread(mz,nz,dz,rhob)
    call zread(mz,nz,dz,attn)
    rp=2.0*rmax
    read(1,*,iostat=iostat)rp
    if (iostat /= 0) then
        print *,'end of file'
    end if

    do i=1,nz+2
        ksqw(i)=(omega/cw(i))**2-k0**2
        ksqb(i)=((omega/cb(i))*(1.0+i_*eta*attn(i)))**2-k0**2
        alpw(i)=sqrt(cw(i)/c0)
        alpb(i)=sqrt(rhob(i)*cb(i)/c0)
    end do
end subroutine

subroutine zread(mz,nz,dz,prof)
    ! Profile reader and interpolator.
    implicit none
    integer*8 ,intent(in)  :: mz,nz
    real*8    ,intent(in)  :: dz
    real*8    ,intent(out) :: prof(mz)

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

subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,ksqw, &
                 ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
    ! The tridiagonal matrices.
    implicit none

    integer*8     ,intent(in)    :: mz,nz,mp,np,iz,jz
    real*8        ,intent(in)    :: dz,k0,rhob(mz),alpw(mz),      &
                                    alpb(mz),ksqw(mz)
    real*8        ,intent(out)   :: f1(mz),f2(mz),f3(mz)
    complex*16    ,intent(in)    :: ksqb(mz),pd1(mp),pd2(mp)
    complex*16    ,intent(out)   :: ksq(mz), r1(mz,mp),r2(mz,mp), &
                                    r3(mz,mp),s1(mz,mp),         &
                                    s2(mz,mp),s3(mz,mp)

    integer*8                    :: i,j,i1,i2
    real*8                       :: a1,a2,a3,c1,c2,c3,            &
                                    cfact,dfact,rfact
    complex*16                   :: d1,d2,d3

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
    implicit none

    integer*8   ,intent(in)    :: mz,nz,mp,np,iz
    complex*16  ,intent(in)    :: r1(mz,mp),r2(mz,mp),r3(mz,mp), &
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

subroutine updat(mr,mz,nz,mp,np,iz,ib,dr,dz,omega,rmax,c0,k0, &
                 r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq, &
                 ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)

    ! Matrix updates.

    use pade_coeffs ,only : pe_pade
    implicit none

    integer*8 ,intent(in)    :: mr,mz,nz,mp,np
    integer*8 ,intent(out)   :: iz,ib
    real*8                   :: dr,dz,omega,rmax,c0,k0,r,rb(mr), &
                                zb(mr),f1(mz),f2(mz),f3(mz)
    real*8                   :: rp,rs
    real*8                   :: cw(mz),cb(mz),rhob(mz),attn(mz), &
                                alpw(mz),alpb(mz),ksqw(mz)
    complex*16               :: ksq(mz),ksqb(mz),                &
                                r1(mz,mp),r2(mz,mp),r3(mz,mp),   &
                                s1(mz,mp),s2(mz,mp),s3(mz,mp)
    complex*16               :: pd1(mp),pd2(mp)

    integer*8                :: jz,ns
    real*8                   :: z

    ! Varying bathymetry.

    if(r >= rb(ib+1)) ib=ib+1
    jz=iz
    z=zb(ib)+(r+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
    iz=int(1.0+z/dz,8)
    iz=max(2,iz)
    iz=min(nz,iz)

    if(iz /= jz) then
        call matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,   &
                   ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
    end if

    ! Varying profiles.

    if(r >= rp)then
        call profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob,attn, &
                alpw,alpb,ksqw,ksqb)
        call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,   &
                   ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
    end if

    ! Turn off the stability constraints.
    if(r >= rs)then
        ns=0
        rs=2.0*rmax
        call pe_pade(mp,np,ns,1_8,k0,dr,pd1,pd2)
        call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,  &
                   ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
    end if

end subroutine

subroutine selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,alpb,  &
                 ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
    ! The self-starter.
    use pade_coeffs ,only : pe_pade
    use constants ,only : i_, pi
    implicit none
    integer*8 ,intent(in)     :: mz,nz,mp,np,ns,iz
    real*8    ,intent(in)     :: zs,dr,dz,k0,rhob(mz),alpw(mz), &
                                  alpb(mz),ksqw(mz)
    real*8    ,intent(out)    :: f1(mz),f2(mz),f3(mz)
    complex*16 ,intent(in)    :: ksqb(mz)
    complex*16 ,intent(inout) :: ksq(mz),v(mz),u(mz)
    complex*16 ,intent(out)   :: pd1(mp),pd2(mp),               &
                                 r1(mz,mp),r2(mz,mp),r3(mz,mp), &
                                 s1(mz,mp),s2(mz,mp),s3(mz,mp)

    integer*8                 :: is
    real*8                    :: dis,si

    ! Conditions for the delta function.

    si=1.0+zs/dz
    is=int(si,8)
    dis=si-float(is)
    u(is)=(1.0-dis)*sqrt(2.0*pi/k0)/(dz*alpw(is))
    u(is+1)=dis*sqrt(2.0*pi/k0)/(dz*alpw(is))

    ! Divide the delta function by (1-X)**2 to get a smooth rhs.

    pd1(1)=0.0
    pd2(1)=-1.0
    call matrc(mz,nz,mp,1_8,iz,iz,dz,k0,rhob,alpw,alpb,ksq,     &
               ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
    call solve(mz,nz,mp,1_8,iz,u,v,r1,r2,r3,s1,s2,s3)
    call solve(mz,nz,mp,1_8,iz,u,v,r1,r2,r3,s1,s2,s3)
    ! Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).

    call pe_pade(mp,np,ns,2_8,k0,dr,pd1,pd2)
    call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw, &
               ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
    call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)

end subroutine

subroutine outln(mz,ir,dir,r,f3,u)
    ! Output transmission loss line
    implicit none
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
    implicit none
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
