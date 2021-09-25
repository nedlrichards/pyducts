module ram_v4_wrapper

use iso_c_binding ,only : c_double, c_int64_t, c_double_complex, c_bool
use ram_v4   ,only      : main

implicit none

contains

subroutine c_main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,dir,dr, &
                  dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,       &
                  ksqb,zs,zr,zmax,zmplt,pline,pgrid,tofile) bind(c)

integer(c_int64_t)         ,intent(in)  :: mr, mz, mp, rn, np, nz, ns, ndr,    &
                                           ndz, iz, nzplt, lz, ib, ir, nprof
real(c_double)             ,intent(in)  :: dir, dr, dz, omega, rmax, c0, k0, r,&
                                           rs, zs, zr, zmax, zmplt

real(c_double)            ,intent(in) ,dimension (nz+2,nprof-1)   :: rhob,alpw,alpb,ksqw
complex(c_double_complex) ,intent(in) ,dimension (nz+2,nprof-1)   :: ksqb

complex(c_double_complex)  ,intent(out) :: pline(nr), pgrid(lz,nr)
logical(c_bool)            ,intent(in)  :: tofile


call main(mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,dz,omega,rmax, &
          c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,ksqb,zs,zr,zmax,zmplt,nprof, &
          pline,pgrid,tofile)

end subroutine

end module
