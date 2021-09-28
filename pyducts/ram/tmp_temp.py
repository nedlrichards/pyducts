import numpy as np
import tmp

rhob=np.empty((2,12), dtype=np.float64)
alpw=np.empty((2,12), dtype=np.float64)
alpb=np.empty((2,12), dtype=np.float64)
ksqw=np.empty((2,12), dtype=np.float64)
ksqb=np.empty((2,12), dtype=np.complex128)

nprof=3
mr=4
rp = np.arange(nprof, dtype=np.float64)
rb = np.arange(mr, dtype=np.float64)
zb = np.arange(mr, dtype=np.float64)

nr=2
lz=3
pline = np.empty(1, dtype=np.complex128, order="F")
pgrid = np.empty((lz,1), dtype=np.complex128, order="F")

jane = tmp.test(nr,2,3,4,5,6,7,8,lz,10,11,nprof,13.,14.,15.,16.,17.,18,19.,20.,
                rp,21.,rb,zb,rhob,alpw,alpb,ksqw,ksqb,22.,23.,24.,25.,mr=mr)

print(jane)
