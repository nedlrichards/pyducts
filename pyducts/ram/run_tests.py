import numpy as np

import refpy

mr = 1000
mz = 80000
mp = 30

(nz,np_r,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,dir_r,dr,dz,pi,eta,eps,
 omega,rmax,c0,k0,ci,r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,
 ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2,tlg) = refpy.setup(mr,mz,mp)

tl = -20.0 * np.log10(np.abs(f3 * u) + eps) + 10.0 * np.log10(r + eps)

