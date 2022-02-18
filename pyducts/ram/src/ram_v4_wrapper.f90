module ram_v4_wrapper

use iso_c_binding ,only : c_double, c_int64_t, c_double_complex, c_bool
use ram_v4        ,only : main
use ram_io        ,only : inram

implicit none

contains

subroutine c_main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,    &
                  dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,    &
                  ksqb,zs,zr,zmax,zmplt,pline,pgrid) bind(c)

integer(c_int64_t)         ,intent(in)  :: mr, mz, mp, nr, np, nz, ns, ndr,    &
                                           ndz, iz, nzplt, lz, ib, ir, nprof
real(c_double)             ,intent(in)  :: dir, dr, dz, omega, rmax, c0, k0, r,&
                                           rs, zs, zr, zmax, zmplt

real(c_double)            ,intent(in) ,dimension (nz+2,nprof-1)   :: rhob,alpw,alpb,ksqw
complex(c_double_complex) ,intent(in) ,dimension (nz+2,nprof-1)   :: ksqb
logical                                                           :: f_tf

real(c_double)             ,intent(out) ,dimension (nprof):: rp
real(c_double)             ,intent(out) ,dimension (mr)   :: rb, zb
complex(c_double_complex)  ,intent(out)                   :: pline(nr), pgrid(lz,nr)

f_tf = .false.

call main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,            &
          dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,            &
          ksqb,zs,zr,zmax,zmplt,pline,pgrid,f_tf)


end subroutine

subroutine c_inram(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,   &
                    dr,dz,omega,rmax,c0,k0,r,rp_np,rs,rb,zb,rhob_np,alpw_np,   &
                    alpb_np,ksqw_np,ksqb_np,zs,zr,zmax,zmplt) bind(c)

integer(c_int64_t) :: mr,mz,mp
integer(c_int64_t) :: nr,nz,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof
real(c_double)     :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt

real(c_double)     ,allocatable ,dimension(:)  :: rp
real(c_double)     ,allocatable ,dimension(:,:):: rhob,alpw,alpb,ksqw
complex(c_double_complex) ,allocatable ,dimension(:,:):: ksqb

real(c_double)     intent(out),dimension(mr)    :: rb,zb

real(c_double)     ,intent(out) ,dimension(nprof)  :: rp_np
real(c_double)     ,intent(out) ,dimension(nz+2,nprof-1) :: rhob_np,alpw_np,alpb_np,ksqw_np
complex(c_double_complex) ,intent(out) ,dimension(nz+2,nprof-1):: ksqb_np


call inram(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,           &
           dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,           &
           ksqb,zs,zr,zmax,zmplt)

rp_np=rp
rhob=rhob_np
alpw=alpw_np
alpb=alpb_np
ksqw=ksqw_np
ksqb=ksqb_np

end subroutine

subroutine c_insize(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,  &
                    dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt) bind(c)

integer(c_int64_t) ,intent(in)  :: mr,mz,mp
integer(c_int64_t) ,intent(out) :: nr,nz,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof
real(c_double)     ,intent(out) :: dir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt

real(c_double)     ,allocatable ,dimension(:)  :: rp
real(c_double)     ,allocatable ,dimension(:,:):: rhob,alpw,alpb,ksqw
complex(c_double_complex) ,allocatable ,dimension(:,:):: ksqb

real(c_double)     ,dimension(mr)    :: rb,zb

call inram(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,           &
           dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,           &
           ksqb,zs,zr,zmax,zmplt)

end subroutine


end module
