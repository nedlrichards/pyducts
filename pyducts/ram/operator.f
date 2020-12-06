      subroutine setup(rmax,dr,ndr,dz,ndz,c0,np,ns,rs,omega,ir,dir,k0,nz,
     >   nzplt,iz,rb,zb,r1,r2,u,v,lz)
c    subroutine setup(nz,np,ns,ndr,ndz,iz,nzplt,lz,ir,
c    >   dir,dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,
c    >   rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,
c    >   s3,pd1,pd2,tlg)
c
      integer mr,mz,mp
      integer nz,np,ns,ndr,ndz,iz,nzplt,lz,ir

      complex*16 u(80000),v(80000),ksq(80000),ksqb(80000),r1(80000,30),
     >   r2(80000,30),r3(80000,30),s1(80000,30),s2(80000,30),
     >   s3(80000,30),pd1(30),pd2(30)
      real*8 k0,rb(1000),zb(1000),cw(80000),cb(80000),rhob(80000),
     >   attn(80000),alpw(80000),alpb(80000),f1(80000),f2(80000),
     >   f3(80000),ksqw(80000),tlg(80000),freq,zs,zr
      real*8 eta, eps
      complex*16 i_
      real*8 pi
c
c     Initialize the parameters, acoustic field, and matrices.
c
cf2py intent(out) rmax,dr,ndr,dz,ndz,c0,np,ns,rs
cf2py intent(out) omega,ir,dir,k0,nz,nzplt,iz,rb,zb
cf2py intent(out) r1,r2,u,v,lz
c
      mr=1000
      mz=80000
      mp=30

      i_=(0.0d0, 1.0d0)
      pi=3.1415926535897932384626433832795d0
      eta=1.0/(40.0*pi*alog10(exp(1.0)))
      eps=1.0e-20

      open(unit=1,status='old',file='ram.in')
      read(1,*)
      read(1,*)freq,zs,zr
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,ns,rs

      r=dr
      omega=2.0*pi*freq
      ri=1.0+zr/dz
      ir=ifix(ri)
      dir=ri-float(ir)
      k0=omega/c0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5
      z=zb(1)
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      if(rs.lt.dr)rs=2.0*rmax


      i=1
    1 read(1,*)rb(i),zb(i)
      if(rb(i).lt.0.0)go to 2
      i=i+1
      go to 1
    2 rb(i)=2.0*rmax
      zb(i)=zb(i-1)



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
c
      do 3 j=1,mp
      r3(1,j)=0.0
      r1(nz+2,j)=0.0
    3 continue
      do 4 i=1,nz+2
      u(i)=0.0
      v(i)=0.0
    4 continue
      lz=0
      do 5 i=ndz,nzplt,ndz
      lz=lz+1
    5 continue

c
c     The initial profiles and starting field.
c
      call profl(mz,nz,i_,dz,eta,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,
     >   alpw,alpb,ksqw,ksqb)
c     call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,c0,k0,rhob,alpw,alpb,ksq,
c    >   ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
c     call outpt(mz,0,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,f3,u,tlg)
c
c     The propagation matrices.
c
c     call epade(mp,np,ns,1,k0,c0,dr,pd1,pd2)
c     call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
c    >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c

c     close(1)
      return
      end

      subroutine profl(mz,nz,ci,dz,eta,omega,rmax,c0,k0,rp,cw,cb,rhob,
     >   attn,alpw,alpb,ksqw,ksqb)
      complex ci,ksqb(mz)
      real k0,cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),alpb(mz),ksqw(mz)
c
      call zread(mz,nz,dz,cw)
      call zread(mz,nz,dz,cb)
      call zread(mz,nz,dz,rhob)
      call zread(mz,nz,dz,attn)
      rp=2.0*rmax
      read(1,*,end=1)rp
c
    1 do 2 i=1,nz+2
      ksqw(i)=(omega/cw(i))**2-k0**2
      ksqb(i)=((omega/cb(i))*(1.0+ci*eta*attn(i)))**2-k0**2
      alpw(i)=sqrt(cw(i)/c0)
      alpb(i)=sqrt(rhob(i)*cb(i)/c0)
    2 continue
c
      return
      end
c
c     Profile reader and interpolator.
c
      subroutine zread(mz,nz,dz,prof)
      real prof(mz)
c
      do 1 i=1,nz+2
      prof(i)=-1.0
    1 continue
      read(1,*)zi,profi
      prof(1)=profi
      i=1.5+zi/dz
      prof(i)=profi
      iold=i
    2 read(1,*)zi,profi
      if(zi.lt.0.0)go to 3
      i=1.5+zi/dz
      if(i.eq.iold)i=i+1
      prof(i)=profi
      iold=i
      go to 2
    3 prof(nz+2)=prof(i)
      i=1
      j=1
    4 i=i+1
      if(prof(i).lt.0.0)go to 4
      if(i-j.eq.1)go to 6
      do 5 k=j+1,i-1
      prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)
    5 continue
    6 j=i
      if(j.lt.nz+2)go to 4
c
      return
      end

      subroutine selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,c0,k0,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
      complex*16 u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
      real*8 k0,rhob(mz),alpw(mz),alpb(mz),f1(mz),f2(mz),f3(mz),
     >   ksqw(mz)
      complex*16 i_
      real*8 pi

c
c     The self-starter.
c     Conditions for the delta function.
c
      i_=(0.0d0, 1.0d0)
      pi=3.1415926535897932384626433832795d0

      si=1.0+zs/dz
      is=ifix(si)
      dis=si-float(is)
      u(is)=(1.0-dis)*sqrt(2.0*pi/k0)/(dz*alpw(is))
      u(is+1)=dis*sqrt(2.0*pi/k0)/(dz*alpw(is))
c
c     Divide the delta function by (1-X)**2 to get a smooth rhs.
c
      pd1(1)=0.0
      pd2(1)=-1.0
c     call matrc(mz,nz,mp,1,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
c    >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c     call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,s1,s2,s3)
c     call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,s1,s2,s3)
c
c     Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(i_*k0*r*sqrt(1+X)).
c
c     call epade(mp,np,ns,2,k0,c0,dr,pd1,pd2)
c     call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
c    >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c     call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
c
      return
      end

      subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,ksqw,
     >   ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      integer mz,nz,mp,np,iz,jz
      complex*16 d1,d2,d3,rfact,ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
      real*8 dz,k0,rhob(mz),alpw(mz),alpb(mz),ksqw(mz),f1(mz),f2(mz),
     >   f3(mz)
c
      a1=k0**2/6.0
      a2=2.0*k0**2/3.0
      a3=k0**2/6.0
      cfact=0.5/dz**2
      dfact=1.0/12.0
c
c     New matrices when iz.eq.jz.
c
      if(iz.eq.jz)then
      i1=2
      i2=nz+1
      do 1 i=1,iz
      f1(i)=1.0/alpw(i)
      f2(i)=1.0
      f3(i)=alpw(i)
      ksq(i)=ksqw(i)
    1 continue
      do 2 i=iz+1,nz+2
      f1(i)=rhob(i)/alpb(i)
      f2(i)=1.0/rhob(i)
      f3(i)=alpb(i)
      ksq(i)=ksqb(i)
    2 continue
      end if
c
c     Updated matrices when iz.ne.jz.
c
      if(iz.gt.jz)then
      i1=jz
      i2=iz+1
      do 3 i=jz+1,iz
      f1(i)=1.0/alpw(i)
      f2(i)=1.0
      f3(i)=alpw(i)
      ksq(i)=ksqw(i)
    3 continue
      end if
c
      if(iz.lt.jz)then
      i1=iz
      i2=jz+1
      do 4 i=iz+1,jz
      f1(i)=rhob(i)/alpb(i)
      f2(i)=1.0/rhob(i)
      f3(i)=alpb(i)
      ksq(i)=ksqb(i)
    4 continue
      end if
c
      do 6 i=i1,i2
c
c     Discretization by Galerkins method.
c
      c1=cfact*f1(i)*(f2(i-1)+f2(i))*f3(i-1)
      c2=-cfact*f1(i)*(f2(i-1)+2.0*f2(i)+f2(i+1))*f3(i)
      c3=cfact*f1(i)*(f2(i)+f2(i+1))*f3(i+1)
      d1=c1+dfact*(ksq(i-1)+ksq(i))
      d2=c2+dfact*(ksq(i-1)+6.0*ksq(i)+ksq(i+1))
      d3=c3+dfact*(ksq(i)+ksq(i+1))
c
      do 5 j=1,np
      r1(i,j)=a1+pd2(j)*d1
      r2(i,j)=a2+pd2(j)*d2
      r3(i,j)=a3+pd2(j)*d3
      s1(i,j)=a1+pd1(j)*d1
      s2(i,j)=a2+pd1(j)*d2
      s3(i,j)=a3+pd1(j)*d3
    5 continue
    6 continue
c
c     The matrix decomposition.
c
      do 9 j=1,np
      do 7 i=i1,iz
      rfact=1.0/(r2(i,j)-r1(i,j)*r3(i-1,j))
      r1(i,j)=r1(i,j)*rfact
      r3(i,j)=r3(i,j)*rfact
      s1(i,j)=s1(i,j)*rfact
      s2(i,j)=s2(i,j)*rfact
      s3(i,j)=s3(i,j)*rfact
    7 continue
c
      do 8 i=i2,iz+2,-1
      rfact=1.0/(r2(i,j)-r3(i,j)*r1(i+1,j))
      r1(i,j)=r1(i,j)*rfact
      r3(i,j)=r3(i,j)*rfact
      s1(i,j)=s1(i,j)*rfact
      s2(i,j)=s2(i,j)*rfact
      s3(i,j)=s3(i,j)*rfact
    8 continue
c
      r2(iz+1,j)=r2(iz+1,j)-r1(iz+1,j)*r3(iz,j)
      r2(iz+1,j)=r2(iz+1,j)-r3(iz+1,j)*r1(iz+2,j)
      r2(iz+1,j)=1.0/r2(iz+1,j)
c
    9 continue
c
      return
      end

      subroutine g_func(g,sig,x,alp,nu)
c
c     The operator function.
c
      complex*16 g
      real*8 sig,x,alp,nu
      complex*16 i_
      parameter (i_=(0.0d0, 1.0d0))

cf2py intent(in) sig,x,alp,nu
cf2py intent(out) g
      g=(1.0d0-nu*x)**2*cdexp(alp*dlog(1.0d0+x)+
     >   i_*sig*(-1.0d0+dsqrt(1.0d0+x)))
      return
      end

      subroutine deriv(n,sig,alp,dg,bin,nu)
      integer n
      complex*16 dg(n),dh1(n),dh2(n)
      real*8 bin(n,n), nu

      complex*16 dh3(n)
      real*8 exp1, exp2, exp3
      complex*16 i_/(0.0d0, 1.0d0)/

cf2py intent(in) n,sig,alp,bin,nu
cf2py intent(out) dg
cf2py depend(n) bin, dg, dh1, dh2
c
c     The derivatives of the operator function at x=0.
c
      dh1(1)=0.5d0*i_*sig
      exp1=-0.5d0
      dh2(1)=alp
      exp2=-1.0d0
      dh3(1)=-2.0d0*nu
      exp3=-1.0d0
      do 1 i=2,n
      dh1(i)=dh1(i-1)*exp1
      exp1=exp1-1.0d0
      dh2(i)=dh2(i-1)*exp2
      exp2=exp2-1.0d0
      dh3(i)=-nu*dh3(i-1)*exp3
      exp3=exp3-1.0d0
    1 continue
c
      dg(1)=1.0d0
      dg(2)=dh1(1)+dh2(1)+dh3(1)
      do 3 i=2,n
      dg(i+1)=dh1(i)+dh2(i)+dh3(i)
      do 2 j=1,i-1
      dg(i+1)=dg(i+1)+bin(i,j)*(dh1(j)+dh2(j)+dh3(j))*dg(i-j+1)
    2 continue
    3 continue
c
      return
      end
