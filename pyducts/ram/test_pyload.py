import numpy as np
import matplotlib.pyplot as plt

from ram_io import RamIn
from ram_main_v4 import main

ram_in = RamIn()

pline, pgrid = main(ram_in.ram_args())

long nr,long nz,long np,long ns,long ndr,
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

