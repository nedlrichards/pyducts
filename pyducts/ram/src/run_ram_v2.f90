program ram
    use ram_io ,only      : inram
    implicit none

    ! mr=bathymetry points, mz=depth grid, mp=pade terms.
    integer*8 ,parameter  :: mr=1000, mz=80000, mp=20

    integer*8  :: np,nz,ns,ndr,ndz,iz,nzplt,lz,ib,ir,mdr,nprof
    real*8     :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt

    complex*16 ,dimension (mp)    :: pd1,pd2
    real*8     ,dimension (mr)    :: rb, zb

    real*8     ,allocatable ,dimension (:)    :: rp
    real*8     ,allocatable ,dimension (:,:)  :: rhob,alpw,alpb,ksqw
    complex*16 ,allocatable ,dimension (:,:)  :: ksqb

    call inram(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz,    &
               omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,       &
               ksqb,zs,zr,zmax,zmplt,nprof)

    call main(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz,           &
               omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,             &
               ksqb,zs,zr,zmax,zmplt,nprof,.true.)


end program

subroutine main(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz,         &
                omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,            &
                ksqb,zs,zr,zmax,zmplt,nprof,tofile)
    use ram_v3 ,only      : setup, updat, solve, selfs, pe_pade, matrc
    use ram_io ,only      : outln, outgr, storeln, storegr
    use constants ,only : pi

    implicit none
    integer*8  :: mr, mz, mp, nprof
    integer*8  :: np,nz,ns,ndr,ndz,iz,nzplt,lz,ib,ir,mdr,prof_i,i,j,numr
    real*8     :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt
    logical    :: tofile

    complex*16 ,dimension (mp)    :: pd1,pd2
    real*8     ,dimension (mr)    :: rb, zb
    real*8     ,dimension (mz)    :: f1,f2,f3
    complex*16 ,dimension (mz)    :: ksq,u,v
    complex*16 ,dimension (mz,mp) :: r1,r2,r3,s1,s2,s3

    real*8     ,dimension (nprof)     :: rp
    real*8     ,dimension (nz+2,nprof-1)   :: rhob,alpw,alpb,ksqw
    complex*16 ,dimension (nz+2,nprof-1)   :: ksqb

    complex*16 ,allocatable ,dimension (:)    :: pline
    complex*16 ,allocatable ,dimension (:,:)  :: pgrid

    do j=1,mp
        r3(1,j)=0.0
        r1(nz+2,j)=0.0
    end do

    do i=1,nz+2
        u(i)=0.0
        v(i)=0.0
    end do

    call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,alpb,ksq,               &
               ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)

    i = 1
    j=1
    r = i * dr

    numr = int(rmax / dr, 8)
    allocate(pline(numr))
    allocate(pgrid(lz,numr))

    if (tofile) then
        open(unit=2,status='unknown',file='tl.line')
        open(unit=3,status='unknown',file='tl.grid',form='unformatted')
        write(3) omega/(2*pi),zs,zr,rmax,dr,ndr,zmax,dz,ndz,zmplt,c0,np,ns,rs,lz
    end if

    !! The propagation matrices.
    call pe_pade(mp,np,ns,1_8,k0,dr,pd1,pd2)
    call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob(:,1),alpw(:,1),alpb(:,1),          &
               ksq,ksqw(:,1),ksqb(:,1),f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)

    if (tofile) then
        call outln(mz,ir,dir,r,f3,u)
        call outgr(mz,ndz,nzplt,lz,r,f3,u)
    else
        call storeln(mz,ir,dir,r,f3,u,i,pline)
        call storegr(mz,ndz,nzplt,lz,r,f3,u,j,pgrid)
    end if
    j = j + 1

    ! March the acoustic field out in range.
    mdr=0
    prof_i = 1
    ib=1
    do i = 2, numr
        call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,omega,rmax,c0,k0,              &
                   r,rp,rs,rb,zb,rhob,alpw,alpb,ksq,                         &
                   ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2,prof_i)

        call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
        r = i * dr

        if (tofile) then
            call outln(mz,ir,dir,r,f3,u)
        else
            call storeln(mz,ir,dir,r,f3,u,i,pline)
        end if

        mdr=mdr+1
        if(mdr == ndr)then
            mdr=0

            if (tofile) then
                call outgr(mz,ndz,nzplt,lz,r,f3,u)
            else
                call storegr(mz,ndz,nzplt,lz,r,f3,u,j,pgrid)
            end if
            j = j + 1

        end if
    end do

    if (tofile) then
        close(2)
        close(3)
    end if

end subroutine
