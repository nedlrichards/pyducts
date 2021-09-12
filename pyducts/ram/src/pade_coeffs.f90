module pade_coeffs
    use constants ,only : i_, pi
    implicit none

contains
    subroutine pe_pade(mp,np,ns,ip,k0,dr,pd1,pd2)
        ! The coefficients of the rational approximation.
        integer*8  ,intent(in)  :: mp,np,ns,ip
        real*8     ,intent(in)  :: k0,dr
        complex*16 ,intent(out) :: pd1(mp),pd2(mp)

        real*8                  :: nu,alp,bin(2*mp,2*mp),fact(2*mp),sig
        integer*8               :: i,j,n
        complex*16              :: z1,dg(2*mp),dh1(2*mp),dh2(2*mp),dh3(2*mp),&
                                   a(2*mp,2*mp),b(2*mp)

        sig=k0*dr
        n=2*np

        if(ip.eq.1)then
            nu=0.0d0
            alp=0.0d0
        else
            nu=1.0d0
            alp=-0.25d0
        end if

        ! The factorials.
        fact(1)=1.0d0
        do i=2,n
            fact(i)=dfloat(i)*fact(i-1)
        end do

        ! The binomial coefficients.
        do i=1,n+1
            bin(i,1)=1.0d0
            bin(i,i)=1.0d0
        end do
        do i=3,n+1
            do j=2,i-1
                bin(i,j)=bin(i-1,j-1)+bin(i-1,j)
            end do
        end do

        do i=1,n
            do j=1,n
                a(i,j)=0.0d0
            end do
        end do

        ! The accuracy constraints.
        call deriv(2*mp,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)

        do i=1,n
            b(i)=dg(i+1)
        end do
        do i=1,n
            if(2*i-1.le.n)a(i,2*i-1)=fact(i)
            do j=1,i
                if(2*j.le.n)a(i,2*j)=-bin(i+1,j+1)*fact(j)*dg(i-j+1)
            end do
        end do

        ! The stability constraints.

        if(ns.ge.1)then
            z1=-3.0d0
            b(n)=-1.0d0
            do j=1,np
                a(n,2*j-1)=z1**j
                a(n,2*j)=0.0d0
            end do
        end if

        if(ns.ge.2)then
        z1=-1.5d0
        b(n-1)=-1.0d0
            do j=1,np
                a(n-1,2*j-1)=z1**j
                a(n-1,2*j)=0.0d0
            end do
        end if

        call gauss(2*mp,n,a,b)

        dh1(1)=1.0d0
        do j=1,np
            dh1(j+1)=b(2*j-1)
        end do
        !call fndrt(dh1,np,dh2,2*mp)
        call cmplx_roots_gen(dh2,dh1,np,.True.,.False.)
        do j=1,np
            pd1(j)=-1.0d0/dh2(j)
        end do

        dh1(1)=1.0d0
        do j=1,np
            dh1(j+1)=b(2*j)
        end do
        !call fndrt(dh1,np,dh2,2*mp)
        call cmplx_roots_gen(dh2,dh1,np,.True.,.False.)
        do j=1,np
            pd2(j)=-1.0d0/dh2(j)
        end do

    end subroutine

    complex function g(sig,x,alp,nu)
        ! The operator function.
        real*8 ,intent(in) :: sig,x,alp,nu

        g=(1.0d0-nu*x)**2*cdexp(alp*log(1.0d0+x)+i_*sig*(-1.0d0+sqrt(1.0d0+x)))
    end function

    subroutine deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
        ! The derivatives of the operator function at x=0.
        integer*8  ,intent(in)  :: m,n
        real*8     ,intent(in)  :: sig,alp,bin(m,m),nu
        complex*16 ,intent(out) :: dg(m),dh1(m),dh2(m),dh3(m)

        integer*8               :: i,j
        real*8                  :: exp1,exp2,exp3

        dh1(1)=0.5d0*i_*sig
        exp1=-0.5d0
        dh2(1)=alp
        exp2=-1.0d0
        dh3(1)=-2.0d0*nu
        exp3=-1.0d0
        do i=2,n
            dh1(i)=dh1(i-1)*exp1
            exp1=exp1-1.0d0
            dh2(i)=dh2(i-1)*exp2
            exp2=exp2-1.0d0
            dh3(i)=-nu*dh3(i-1)*exp3
            exp3=exp3-1.0d0
        end do

        dg(1)=1.0d0
        dg(2)=dh1(1)+dh2(1)+dh3(1)
        do i=2,n
            dg(i+1)=dh1(i)+dh2(i)+dh3(i)
            do j=1,i-1
                dg(i+1)=dg(i+1)+bin(i,j)*(dh1(j)+dh2(j)+dh3(j))*dg(i-j+1)
            end do
        end do

    end subroutine

    subroutine gauss(m,n,a,b)
        ! Gaussian elimination.
        integer*8  ,intent(in)    :: m,n
        complex*16 ,intent(inout) :: a(m,m),b(m)

        integer*8                 :: i,j,k
        ! Downward elimination.
        do i=1,n
            if(i.lt.n)call pivot(m,n,i,a,b)
            a(i,i)=1.0d0/a(i,i)
            b(i)=b(i)*a(i,i)
            if(i.lt.n)then
                do j=i+1,n
                    a(i,j)=a(i,j)*a(i,i)
                end do
                do k=i+1,n
                    b(k)=b(k)-a(k,i)*b(i)
                    do j=i+1,n
                        a(k,j)=a(k,j)-a(k,i)*a(i,j)
                    end do
                end do
            end if
        end do
        ! Back substitution.
        do i=n-1,1,-1
            do j=i+1,n
                b(i)=b(i)-a(i,j)*b(j)
            end do
        end do
    end subroutine

    subroutine pivot(m,n,i,a,b)
        ! Rows are interchanged for stability.
        integer*8  ,intent(in)    :: m,n,i
        complex*16 ,intent(inout) :: a(m,m),b(m)

        integer*8                 :: j,i0
        real*8                    :: amp0,amp
        complex*16                :: temp

        i0=i
        amp0=cdabs(a(i,i))
        do j=i+1,n
            amp=cdabs(a(j,i))
            if(amp.gt.amp0)then
                i0=j
                amp0=amp
            end if
        end do

        if(i0.eq.i)return

        temp=b(i)
        b(i)=b(i0)
        b(i0)=temp
        do j=i,n
            temp=a(i,j)
            a(i,j)=a(i0,j)
            a(i0,j)=temp
        end do
    end subroutine

end module
