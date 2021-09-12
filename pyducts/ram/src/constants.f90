module constants
    ! Constants contain more digits than double precision, so that
    ! they are rounded correctly. Single letter constants contain underscore so
    ! that they do not clash with user variables ("e" and "i" are frequently used as
    ! loop variables)

    real*8 pi, e_
    complex*16 i_

    parameter (pi=3.1415926535897932384626433832795d0, &
               e_=2.7182818284590452353602874713527d0, &
               i_=(0.0d0, 1.0d0))
end module
