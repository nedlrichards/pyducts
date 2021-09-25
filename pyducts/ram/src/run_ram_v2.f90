program ram
    use ram_io ,only      : inram
    use ram_v4 ,only      : main
    implicit none

    ! mr=bathymetry points, mz=depth grid, mp=pade terms.
    integer*8 ,parameter  :: mr=1000, mz=80000, mp=20

    integer*8  :: nr,np,nz,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof
    real*8     :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt

    complex*16 ,dimension (mp)    :: pd1,pd2
    real*8     ,dimension (mr)    :: rb, zb

    real*8     ,allocatable ,dimension (:)    :: rp
    real*8     ,allocatable ,dimension (:,:)  :: rhob,alpw,alpb,ksqw
    complex*16 ,allocatable ,dimension (:)    :: pline
    complex*16 ,allocatable ,dimension (:,:)  :: ksqb,pgrid

    call inram(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,dr,dz, &
               omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,             &
               ksqb,zs,zr,zmax,zmplt)

    allocate(pline(nr))
    allocate(pgrid(lz,nr))

    call main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,dr,dz,  &
               omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,             &
               ksqb,zs,zr,zmax,zmplt,pline,pgrid,.true.)
end program
