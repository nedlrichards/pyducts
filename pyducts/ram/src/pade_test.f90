program main
    use pade_coeffs ,only : pe_pade
    implicit none

    integer*8 :: np,ns,ip,mp
    real*8 :: k0,dr
    complex*16 ,dimension(:) ,allocatable :: pd1,pd2

    mp=20
    np=8
    ns=2
    ip=2
    k0=1/10.
    dr=200.
    allocate(pd1(mp), pd2(mp))

    call pe_pade(mp,np,ns,ip,k0,dr,pd1,pd2)
    print '(F8.3,SP,F8.3,"i")', pd1(1:np)
    print *,''
    print '(F8.3,SP,F8.3,"i")', pd2(1:np)
end program
