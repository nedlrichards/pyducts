from numpy import empty, float64, complex128, asarray
from numpy cimport ndarray, float64_t, complex128_t

cdef extern:
    void c_main(long *mr,long *nr,long *mz,long *nz,long *mp,long *np,long *ns,long *ndr,
                long *ndz,long *iz,long *nzplt,long *lz,long *ib,long *ir,long *nprof,
                double *rdir,double *dr,double *dz,double *omega,double *rmax,
                double *c0,double *k0,double *r,double *rp,double *rs,
                double *rb,double *zb,double *rhob,double *alpw,double *alpb,
                double *ksqw,complex *ksqb,double *zs,double *zr,double *zmax,
                double *zmplt,complex *pline,complex *pgrid)

    void c_inram(long *mr,long *nr,long *mz,long *nz,long *mp,long *np,
                 long *ns,long *ndr,long *ndz,long *iz,long *nzplt,long *lz,
                 long *ib,long *ir,long *nprof,double *rdir,double *dr,
                 double *dz,double *omega,double *rmax,double *c0,
                 double *k0,double *r,double *rp,double *rs,double *rb,
                 double *zb,double *rhob,double *alpw,double *alpb,double *ksqw,
                 complex *ksqb,double *zs,double *zr,double *zmax,double *zmplt)


def main(long nr,long nz,long np,long ns,long ndr,
         long ndz,long iz,long nzplt,long lz,long ib,long ir,long nprof,
         double rdir,double dr,double dz,double omega,double rmax,
         double c0,double k0,double r,
         ndarray[float64_t,ndim=1] rp,
         double rs,
         ndarray[float64_t,ndim=1] rb,
         ndarray[float64_t,ndim=1] zb,
         ndarray[float64_t,ndim=2] rhob,
         ndarray[float64_t,ndim=2] alpw,
         ndarray[float64_t,ndim=2] alpb,
         ndarray[float64_t,ndim=2] ksqw,
         ndarray[complex128_t,ndim=2] ksqb,
         double zs,double zr,double zmax,double zmplt,
         long mr=1000,long mz=80000, long mp=20):

    cdef ndarray[complex128_t,ndim=1] pline = empty(nr, dtype=complex128, order="F")
    cdef ndarray[complex128_t,ndim=2] pgrid = empty((lz,nr), dtype=complex128, order="F")
    cdef ndarray[float64_t,ndim=2] rhob_f = asarray(rhob, order="F")
    cdef ndarray[float64_t,ndim=2] alpw_f = asarray(alpw, order="F")
    cdef ndarray[float64_t,ndim=2] alpb_f = asarray(alpb, order="F")
    cdef ndarray[float64_t,ndim=2] ksqw_f = asarray(ksqw, order="F")
    cdef ndarray[complex128_t,ndim=2] ksqb_f = asarray(ksqb, order="F")

    c_main(&mr,&nr,&mz,&nz,&mp,&np,&ns,&ndr,&ndz,&iz,&nzplt,&lz,&ib,&ir,&nprof,
           &rdir,&dr,&dz,&omega,&rmax,&c0,&k0,&r,&rp[0],&rs,&rb[0],&zb[0],
           &rhob_f[0,0],
           &alpw_f[0,0],
           &alpb_f[0,0],
           &ksqw_f[0,0],
           &ksqb_f[0,0],
           &zs,&zr,&zmax,&zmplt,
           &pline[0],&pgrid[0,0])

    return pline, pgrid

def inram(long mr=1000,long mz=80000, long mp=20, long *nz, long *num_prof):

    cdef long   nr,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,
    cdef double rdir,dr,dz,omega,rmax,c0,k0,r,rs,zs,zr,zmax,zmplt

    cdef ndarray[float64_t,ndim=1] rp = np.empty(num_prof)
    cdef ndarray[float64_t,ndim=1] rb = np.empty(mr)
    cdef ndarray[float64_t,ndim=1] zb = np.empty(mr)

    cdef ndarray[float64_t,ndim=2] rhob,alpw,alpb,ksqw
    cdef ndarray[complex128_t,ndim=2] ksqb

    c_inram(&mr,&nr,&mz,&nz,&mp,&np,&ns,&ndr,&ndz,&iz,&nzplt,&lz,&ib,&ir,&nprof,
           &rdir,&dr,&dz,&omega,&rmax,&c0,&k0,&r,&rp[0],&rs,&rb[0],&zb[0],
           &rhob[0,0],&alpw[0,0],&alpb[0,0],&ksqw[0,0],&ksqb[0,0],
           &zs,&zr,&zmax,&zmplt)

    return mr,nr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,nprof,rdir,dr,dz, \
           omega,rmax,c0,k0,r,rp,rs,rb,zb,rhob,alpw,alpb,ksqw,ksqb,zs,zr,zmax,zmplt
