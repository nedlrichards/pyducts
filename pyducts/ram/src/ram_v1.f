      program ram
c
c     ******************************************************************
c     ***** Range-dependent Acoustic Model, Version 1.5, 13-Sep-00 *****
c     ******************************************************************
c
c     This code was developed by Michael D. Collins at the Naval
c     Research Laboratory in Washington, DC. It solves range-dependent
c     ocean acoustics problems with the split-step Pade algorithm
c     [M. D. Collins, J. Acoust. Soc. Am. 93, 1736-1742 (1993)]. A
c     users guide and updates of the code are available via anonymous
c     ftp from ram.nrl.navy.mil.
c
c     Version 1.5 contains a correction to a bug in the dimension of
c     quantities passed to subroutines fndrt and guerre that Laurie
c     Fialkowski noticed.
c
c     Version 1.4 contains a correction to a minor bug in subroutine
c     guerre that Dave King noticed (amp1 and amp2 were declared
c     twice) and a few other minor improvements.
c
c     Version 1.3 contains a new root-finding subroutine.
c
c     Version 1.2 contains a minor modification. The output to tl.grid
c     is no longer zeroed out along the ocean bottom. This was done in
c     previous versions so that the ocean bottom would be highlighted
c     in graphical displays. The graphics codes ramclr, ramctr, and
c     ramcc read in the bathymetry from ram.in and plot the ocean
c     bottom directly.
c
c     Version 1.1 contains two improvements:
c
c     (1) An improved self starter. Stability is improved by using the
c     factor (1-X)**2 instead of (1+X)**2 to smooth the delta function.
c     The factor (1+X)**2 is nearly singular for some problems involving
c     deep water and/or weak attenuation. Numerical problems associated
c     with this singularity were detected by Eddie Scheer of Woods Hole
c     Oceanographic Institute.
c
c     (2) Elimination of underflow problems. A very small number is
c     added to the solution in subroutine solve to prevent underflow,
c     which can adversely affect run time on some computers. This
c     improvement was suggested by Ed McDonald of the SACLANT Undersea
c     Research Centre.
c
      implicit none
      integer*8 mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir
      complex*16 ksq,ksqb,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2
      real*8 dir,dr,dz,omega,rmax,c0,k0,r,rp,rs
      real*8 rb,zb,cw,cb,rhob,attn,alpw,alpb,ksqw,f1,f2,f3
c
c     mr=bathymetry points, mz=depth grid, mp=pade terms.
c
      parameter (mr=1000,mz=80000,mp=20)
      dimension rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),
     >   alpb(mz),f1(mz),f2(mz),f3(mz),ksq(mz),ksqw(mz),ksqb(mz),u(mz),
     >   v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
c
      open(unit=1,status='old',file='ram.in')
      open(unit=2,status='unknown',file='tl.line')
      open(unit=3,status='unknown',file='tl.grid',form='unformatted')
c
      call setup(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,
     >   dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,rhob,
     >   attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
     >   pd1,pd2)
c
c     March the acoustic field out in range.
c

      mdr=0
    1 call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,omega,rmax,c0,k0,r,
     >   rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,
     >   r1,r2,r3,s1,s2,s3,pd1,pd2)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
      r=r+dr

      call outln(mz,ir,dir,r,f3,u)

      mdr=mdr+1
      if(mdr.eq.ndr)then
      mdr=0
      call outgr(mz,ndz,nzplt,lz,r,f3,u)
      end if

      if(r.lt.rmax)go to 1
c
      close(1)
      close(2)
      close(3)
c
      stop
      end

      subroutine setup(mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,
     >   dir,dr,dz,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,cb,
     >   rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,
     >   s3,pd1,pd2)
c
c     Initialize the parameters, acoustic field, and matrices.
c
      implicit none
      integer*8 mr,mz,nz,mp,np,ns,ndr,ndz,iz,nzplt,lz,ib,ir,i,j
      real*8 dir,dr,dz,pi,eps,omega,rmax,c0,k0,r,rp,rs,z,zs,
     >   zmax,ri,freq,zmplt,zr,rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),
     >   attn(mz),alpw(mz),alpb(mz),ksqw(mz),f1(mz),f2(mz),f3(mz)
      complex*16 ci,u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
c
      ci=cmplx(0.0,1.0,8)
      eps=1.0e-20
      pi=3.1415926535897932384626433832795d0

      read(1,*)
      read(1,*)freq,zs,zr
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,ns,rs

C     mbp: adding grid size info to the TL output file
C      write( 3 ) dz, ndz, zmplt, dr, ndr, rmax, freq, zs

c
      i=1
    1 read(1,*)rb(i),zb(i)
      if(rb(i).lt.0.0)go to 2
      i=i+1
      go to 1
    2 rb(i)=2.0*rmax
      zb(i)=zb(i-1)
c
      ib=1
      r=dr
      omega=2.0*pi*freq
      ri=1.0+zr/dz
      ir=int(ri,8)
      dir=ri-float(ir)
      k0=omega/c0
      nz=int(zmax/dz-0.5,8)
      nzplt=int(zmplt/dz-0.5,8)
      z=zb(1)
      iz=int(1.0+z/dz,8)
      iz=max(2,iz)
      iz=min(nz,iz)
      if(rs.lt.dr)rs=2.0*rmax
c
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
c     write(3)lz
c     header containing a more detailed description of the problem (jcp)
      write(3) freq,zs,zr,rmax,dr,ndr,zmax,dz,ndz,zmplt,c0,np,ns,rs,lz

c
c     The initial profiles and starting field.
c
      call profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,
     >   alpw,alpb,ksqw,ksqb)
      call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,alpb,ksq,
     >   ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
      call outln(mz,ir,dir,r,f3,u)
      call outgr(mz,ndz,nzplt,lz,r,f3,u)

c
c     The propagation matrices.
c
      call epade(mp,np,ns,1_8,k0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
      return
      end

      subroutine profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob,
     >   attn,alpw,alpb,ksqw,ksqb)
c
c     Set up the profiles.
c
      implicit none
      integer*8 mz,nz,i
      real*8 dz,eta,omega,rmax,c0,k0,rp,cw(mz),cb(mz),rhob(mz),
     >   attn(mz),alpw(mz),alpb(mz),ksqw(mz)
      complex*16 ci,ksqb(mz)

      ci=cmplx(0.0,1.0,8)
      eta=0.01832338997198569352181968569348d0
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

      subroutine zread(mz,nz,dz,prof)
c
c     Profile reader and interpolator.
c
      implicit none
      integer*8 mz,nz,i,j,k,iold
      real*8 dz,zi,profi,prof(mz)
c
      do 1 i=1,nz+2
      prof(i)=-1.0
    1 continue
      read(1,*)zi,profi
      prof(1)=profi
      i=int(1.5+zi/dz,8)
      prof(i)=profi
      iold=i
    2 read(1,*)zi,profi
      if(zi.lt.0.0)go to 3
      i=int(1.5+zi/dz, 8)
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

      subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,ksqw,
     >   ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
c     The tridiagonal matrices.
c
      implicit none
      integer*8 mz,nz,mp,np,iz,jz,i,j,i1,i2
      real*8 dz,k0,rhob(mz),f1(mz),f2(mz),f3(mz),alpw(mz),alpb(mz),
     >   ksqw(mz),a1,a2,a3,c1,c2,c3,cfact,dfact
      complex*16 d1,d2,d3,rfact,ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
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

      subroutine solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
c
c     The tridiagonal solver.
c
      implicit none
      integer*8 mz,nz,mp,np,iz,i,j
      real*8 eps
      complex*16 u(mz),v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),
     >   s1(mz,mp),s2(mz,mp),s3(mz,mp)
      eps=1.0e-30
c
      do 6 j=1,np
c
c     The right side.
c
      do 1 i=2,nz+1
      v(i)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)+eps
    1 continue
c
c     The elimination steps.
c
      do 2 i=3,iz
      v(i)=v(i)-r1(i,j)*v(i-1)+eps
    2 continue
      do 3 i=nz,iz+2,-1
      v(i)=v(i)-r3(i,j)*v(i+1)+eps
    3 continue
c
      u(iz+1)=(v(iz+1)-r1(iz+1,j)*v(iz)-r3(iz+1,j)*v(iz+2))*
     >   r2(iz+1,j)+eps
c
c     The back substitution steps.
c
      do 4 i=iz,2,-1
      u(i)=v(i)-r3(i,j)*u(i+1)+eps
    4 continue
      do 5 i=iz+2,nz+1
      u(i)=v(i)-r1(i,j)*u(i-1)+eps
    5 continue
    6 continue
c
      return
      end

      subroutine updat(mr,mz,nz,mp,np,iz,ib,dr,dz,omega,rmax,c0,k0,
     >   r,rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,
     >   f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
c     Matrix updates.
c
      implicit none
      integer*8 mr,mz,nz,mp,np,iz,ib,jz,ns
      real*8 dr,dz,omega,rmax,c0,k0,r,rp,rs,rb(mr),z,zb(mr),
     >   attn(mz),cb(mz),rhob(mz),cw(mz),ksqw(mz),f1(mz),f2(mz),f3(mz),
     >   alpw(mz),alpb(mz)
      complex*16 ci,ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),
     >   s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
c
c     Varying bathymetry.
c

      ci=cmplx(0.0d0,1.0d0,8)
      if(r.ge.rb(ib+1))ib=ib+1
      jz=iz
      z=zb(ib)+(r+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
      iz=int(1.0+z/dz,8)
      iz=max(2,iz)
      iz=min(nz,iz)
      if(iz.ne.jz)call matrc(mz,nz,mp,np,iz,jz,dz,k0,rhob,alpw,alpb,ksq,
     >   ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
c     Varying profiles.
c
      if(r.ge.rp)then
      call profl(mz,nz,dz,omega,rmax,c0,k0,rp,cw,cb,rhob,attn,
     >   alpw,alpb,ksqw,ksqb)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      end if
c
c     Turn off the stability constraints.
c
      if(r.ge.rs)then
      ns=0
      rs=2.0*rmax
      call epade(mp,np,ns,1_8,k0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      end if
c
      return
      end
      subroutine selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,k0,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
c     The self-starter.
c
      implicit none
      integer*8 mz,nz,mp,np,ns,iz,is
      real*8 zs,dr,dz,pi,k0,rhob(mz),alpw(mz),alpb(mz),dis,si,
     >   f1(mz),f2(mz),f3(mz),ksqw(mz)
      complex*16 u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
c
c     Conditions for the delta function.
c
      pi=3.1415926535897932384626433832795d0
      si=1.0+zs/dz
      is=int(si,8)
      dis=si-float(is)
      u(is)=(1.0-dis)*sqrt(2.0*pi/k0)/(dz*alpw(is))
      u(is+1)=dis*sqrt(2.0*pi/k0)/(dz*alpw(is))
c
c     Divide the delta function by (1-X)**2 to get a smooth rhs.
c
      pd1(1)=0.0
      pd2(1)=-1.0
      call matrc(mz,nz,mp,1_8,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      call solve(mz,nz,mp,1_8,iz,u,v,r1,r2,r3,s1,s2,s3)
      call solve(mz,nz,mp,1_8,iz,u,v,r1,r2,r3,s1,s2,s3)

c
c     Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
c
      call epade(mp,np,ns,2_8,k0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhob,alpw,alpb,ksq,ksqw,ksqb,
     >   f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
c
      return
      end

      subroutine outln(mz,ir,dir,r,f3,u)
c
c     Output transmission loss line
c
      implicit none
      integer*8 mz,ir
      real*8 dir,eps,r,f3(mz)
      complex*16 ur,u(mz),pout
      eps=1e-20
c
      ur=(1.0-dir)*f3(ir)*u(ir)+dir*f3(ir+1)*u(ir+1)
c     tl=-20.0*log10(abs(ur)+eps)+10.0*log10(r+eps)
      pout=ur/sqrt(r+eps)
      write(2,*)r,real(pout),imag(pout)
      return
      end

      subroutine outgr(mz,ndz,nzplt,lz,r,f3,u)
c
c     Output transmission loss.
c
      implicit none
      integer*8 mz,ndz,nzplt,lz,i,j
      real*8 eps,r,f3(mz)
      complex*16 ur,u(mz),poutg(mz)
      eps=1e-20

      j=0
      do 1 i=ndz,nzplt,ndz
      ur=u(i)*f3(i)
      j=j+1
c     tlg(j)=-20.0*log10(abs(ur)+eps)+10.0*log10(r+eps)
      poutg(j)=ur/sqrt(r+eps)
    1 continue
      write(3)(poutg(j),j=1,lz)

      return
      end

      subroutine epade(mp,np,ns,ip,k0,dr,pd1,pd2)
c
c     The coefficients of the rational approximation.
c
      implicit none
      integer*8 mp,np,ns,ip,i,j,n
      real*8 nu,k0,dr,alp,bin(2*mp,2*mp),fact(2*mp),pi,sig
      complex*16 ci,z1,dg(2*mp),dh1(2*mp),dh2(2*mp),dh3(2*mp),
     >   a(2*mp,2*mp),b(2*mp),pd1(mp),pd2(mp)

      pi=3.1415926535897932384626433832795d0
      ci=cmplx(0.0d0,1.0d0,8)
      sig=k0*dr
      n=2*np
c
      if(ip.eq.1)then
      nu=0.0d0
      alp=0.0d0
      else
      nu=1.0d0
      alp=-0.25d0
      end if
c
c     The factorials.
c
      fact(1)=1.0d0
      do 1 i=2,n
      fact(i)=dfloat(i)*fact(i-1)
    1 continue
c
c     The binomial coefficients.
c
      do 2 i=1,n+1
      bin(i,1)=1.0d0
      bin(i,i)=1.0d0
    2 continue
      do 4 i=3,n+1
      do 3 j=2,i-1
      bin(i,j)=bin(i-1,j-1)+bin(i-1,j)
    3 continue
    4 continue
c
      do 6 i=1,n
      do 5 j=1,n
      a(i,j)=0.0d0
    5 continue
    6 continue
c
c     The accuracy constraints.
c
      call deriv(2*mp,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
c
      do 7 i=1,n
      b(i)=dg(i+1)
    7 continue
      do 9 i=1,n
      if(2*i-1.le.n)a(i,2*i-1)=fact(i)
      do 8 j=1,i
      if(2*j.le.n)a(i,2*j)=-bin(i+1,j+1)*fact(j)*dg(i-j+1)
    8 continue
    9 continue
c
c     The stability constraints.
c
      if(ns.ge.1)then
      z1=-3.0d0
      b(n)=-1.0d0
      do 10 j=1,np
      a(n,2*j-1)=z1**j
      a(n,2*j)=0.0d0
   10 continue
      end if
c
      if(ns.ge.2)then
      z1=-1.5d0
      b(n-1)=-1.0d0
      do 11 j=1,np
      a(n-1,2*j-1)=z1**j
      a(n-1,2*j)=0.0d0
   11 continue
      end if
c
      call gauss(2*mp,n,a,b)
c
      dh1(1)=1.0d0
      do 12 j=1,np
      dh1(j+1)=b(2*j-1)
   12 continue
      call fndrt(dh1,np,dh2,2*mp)
      do 13 j=1,np
      pd1(j)=-1.0d0/dh2(j)
   13 continue
c
      dh1(1)=1.0d0
      do 14 j=1,np
      dh1(j+1)=b(2*j)
   14 continue
      call fndrt(dh1,np,dh2,2*mp)
      do 15 j=1,np
      pd2(j)=-1.0d0/dh2(j)
   15 continue
c
      return
      end

      function g(sig,x,alp,nu)
c
c     The operator function.
c
      implicit none
      real*8 alp,sig,x,nu
      complex*16 ci,g

      ci=cmplx(0.0d0,1.0d0,8)
      g=(1.0d0-nu*x)**2*cdexp(alp*log(1.0d0+x)+
     >   ci*sig*(-1.0d0+sqrt(1.0d0+x)))
      return
      end

      subroutine deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
c
c     The derivatives of the operator function at x=0.
c
      implicit none
      integer*8 m,n,i,j
      real*8 sig,alp,bin(m,m),nu,exp1,exp2,exp3
      complex*16 ci,dg(m),dh1(m),dh2(m),dh3(m)
c
      ci=cmplx(0.0d0,1.0d0,8)
c
      dh1(1)=0.5d0*ci*sig
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

      subroutine gauss(m,n,a,b)
c
c     Gaussian elimination.
c
      implicit none
      integer*8 m,n,i,j,k
      complex*16 a(m,m),b(m)
c
c     Downward elimination.
c
      do 4 i=1,n
      if(i.lt.n)call pivot(m,n,i,a,b)
      a(i,i)=1.0d0/a(i,i)
      b(i)=b(i)*a(i,i)
      if(i.lt.n)then
      do 1 j=i+1,n
      a(i,j)=a(i,j)*a(i,i)
    1 continue
      do 3 k=i+1,n
      b(k)=b(k)-a(k,i)*b(i)
      do 2 j=i+1,n
      a(k,j)=a(k,j)-a(k,i)*a(i,j)
    2 continue
    3 continue
      end if
    4 continue
c
c     Back substitution.
c
      do 6 i=n-1,1,-1
      do 5 j=i+1,n
      b(i)=b(i)-a(i,j)*b(j)
    5 continue
    6 continue
c
      return
      end

      subroutine pivot(m,n,i,a,b)
c
c     Rows are interchanged for stability.
c
      implicit none
      integer*8 m,n,i,j,i0
      real*8 amp0,amp
      complex*16 temp,a(m,m),b(m)
c
      i0=i
      amp0=cdabs(a(i,i))
      do 1 j=i+1,n
      amp=cdabs(a(j,i))
      if(amp.gt.amp0)then
      i0=j
      amp0=amp
      end if
    1 continue
      if(i0.eq.i)return
c
      temp=b(i)
      b(i)=b(i0)
      b(i0)=temp
      do 2 j=i,n
      temp=a(i,j)
      a(i,j)=a(i0,j)
      a(i0,j)=temp
    2 continue
c
      return
      end

      subroutine fndrt(a,n,z,m)
c
c     The root-finding subroutine.
c
      implicit none
      integer*8 n,m,i,k
      real*8 err
      complex*16 a(m),z(m),root
c
      if(n.eq.1)then
      z(1)=-a(1)/a(2)
      return
      end if
      if(n.eq.2)go to 4
c
      do 3 k=n,3,-1
c
c     Obtain an approximate root.
c
      root=0.0d0
      err=1.0d-12
      call guerre(a,k,m,root,err,1000_8)
c
c     Refine the root by iterating five more times.
c
      err=0.0d0
      call guerre(a,k,m,root,err,5_8)
      z(k)=root
c
c     Divide out the factor (z-root).
c
      do 1 i=k,1,-1
      a(i)=a(i)+root*a(i+1)
    1 continue
      do 2 i=1,k
      a(i)=a(i+1)
    2 continue
c
    3 continue
c
c     Solve the quadratic equation.
c
    4 z(2)=0.5*(-a(2)+sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
      z(1)=0.5*(-a(2)-sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
c
      return
      end

      subroutine guerre(a,n,m,z,err,nter)
c
c     This subroutine finds a root of a polynomial of degree n > 2
c     by Laguerres method.
c
      implicit none
      integer*8 m,n,i,nter,iter,jter
      real*8 amp1,amp2,rn,eps,err
      complex*16 a(m),az(50),azz(50),z,dz,p,pz,pzz,f,g,h,ci
c
      ci=cmplx(0.0d0,1.0d0,8)
      eps=1.0d-20
      rn=real(n)
c
c     The coefficients of p_der(z) and p_der_der(z).
c
      do 1 i=1,n
      az(i)=float(i)*a(i+1)
    1 continue
      do 2 i=1,n-1
      azz(i)=float(i)*az(i+1)
    2 continue
c
      iter=0
    3 p=a(n)+a(n+1)*z
      do 4 i=n-1,1,-1
      p=a(i)+z*p
    4 continue
      if(abs(p).lt.eps)return
c
      pz=az(n-1)+az(n)*z
      do 5 i=n-2,1,-1
      pz=az(i)+z*pz
    5 continue
c
      pzz=azz(n-2)+azz(n-1)*z
      do 6 i=n-3,1,-1
      pzz=azz(i)+z*pzz
    6 continue
c
c     The Laguerre perturbation.
c
      f=pz/p
      g=f**2-pzz/p
      h=sqrt((rn-1.0d0)*(rn*g-f**2))
      amp1=abs(f+h)
      amp2=abs(f-h)
      if(amp1.gt.amp2)then
      dz=-rn/(f+h)
      else
      dz=-rn/(f-h)
      end if
c
      iter=iter+1
c
c     Rotate by 90 degrees to avoid limit cycles.
c
      jter=jter+1
      if(jter.eq.10)then
      jter=1
      dz=dz*ci
      end if
      z=z+dz
c
      if(iter.eq.100)then
      write(*,*)' '
      write(*,*)'   Laguerre method not converging.'
      write(*,*)'   Try a different combination of DR and NP.'
      write(*,*)' '
      stop
      end if
c
      if((abs(dz).gt.err).and.(iter.lt.nter))go to 3
c
      return
      end
