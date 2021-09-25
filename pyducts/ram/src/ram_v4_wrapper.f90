module ram_v4_wrapper

use iso_c_binding ,only : c_double, c_int64_t, c_double_complex, c_bool
use ram_v4   ,only      : main

implicit none

contains

subroutine c_main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,    &
                  dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,    &
                  ksqb,zs,zr,zmax,zmplt,pline,pgrid,tofile) bind(c)

integer(c_int64_t)         ,intent(in)  :: mr, mz, mp, nr, np, nz, ns, ndr,    &
                                           ndz, iz, nzplt, lz, ib, ir, nprof
real(c_double)             ,intent(in)  :: dir, dr, dz, omega, rmax, c0, k0, r,&
                                           rs, zs, zr, zmax, zmplt

real(c_double)            ,intent(in) ,dimension (nz+2,nprof-1)   :: rhob,alpw,alpb,ksqw
complex(c_double_complex) ,intent(in) ,dimension (nz+2,nprof-1)   :: ksqb
logical(c_bool)           ,intent(in)                             :: tofile
logical                                                           :: f_tf

real(c_double)             ,intent(out) ,dimension (nprof):: rp
real(c_double)             ,intent(out) ,dimension (mr)   :: rb, zb
complex(c_double_complex)  ,intent(out)                   :: pline(nr), pgrid(lz,nr)

f_tf = tofile

call main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,            &
          dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,            &
          ksqb,zs,zr,zmax,zmplt,pline,pgrid,f_tf)


end subroutine

end module
