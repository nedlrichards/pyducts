module pade_coeffs_wrapper

use iso_c_binding ,only : c_double, c_int64_t, c_double_complex
use pade_coeffs   ,only : pe_pade

implicit none

contains

subroutine c_pe_pade(mp,np,ns,ip,k0,dr,pd1,pd2) bind(c)
    integer(c_int64_t)         ,intent(in)  :: mp,np,ns,ip
    real(c_double)             ,intent(in)  :: k0,dr
    complex(c_double_complex)  ,intent(out) :: pd1(mp), pd2(mp)
    call pe_pade(mp,np,ns,ip,k0,dr,pd1,pd2)
end subroutine

end module
