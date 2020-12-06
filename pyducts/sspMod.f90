MODULE sspmod

  ! holds SSP input by user and associated variables

  USE FatalError

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER, PRIVATE :: ENVFile = 5
  INTEGER, PARAMETER          :: MaxSSP = 20001, MaxMedia = 501
  INTEGER, PRIVATE            :: N, iz, ILoc, Lay
  INTEGER                     :: ISSP
  REAL (KIND=8)               :: alphaR = 1500, betaR = 0, alphaI = 0, betaI = 0, rhoR = 1
  REAL (KIND=8), PRIVATE      :: h, z, R

  ! SSP
  INTEGER           :: SLoc( MaxMedia ), SNPts( MaxMedia ), SNMedia
  REAL     (KIND=8) :: Sz( MaxSSP ), SalphaR( MaxSSP ), SalphaI( MaxSSP ), Srho( MaxSSP ), SbetaR( MaxSSP ), SbetaI( MaxSSP )
  REAL     (KIND=8) :: SDepth( MaxMedia ), Ssigma( MaxMedia ), Sbeta( MaxMedia ), SfT( MaxMedia )
  COMPLEX  (KIND=8) :: Scp( MaxSSP ), Scs( MaxSSP ), Sn2( MaxSSP ),  &
                       ScpSpline( 4, MaxSSP ), ScsSpline( 4, MaxSSP ), SrhoSpline( 4, MaxSSP )
  COMPLEX  (KIND=8) :: ScpCoef( 4, MaxSSP ), ScsCoef( 4, MaxSSP ), SrhoCoef( 4, MaxSSP ), ScsWork( 4, MaxSSP )   ! for the PCHIP coefficients
  COMPLEX  (KIND=8) :: SrhoT( MaxSSP )   ! temporary values for calling PCHIP (ordinate argument must be complex valued!)
  CHARACTER (LEN=1) :: SType
  CHARACTER (LEN=2) :: SAttenUnit
  CHARACTER (LEN=8) :: SMaterial( MaxMedia )

  CHARACTER (LEN=1) :: TBC                            ! Boundary condition type
  REAL     (KIND=8) :: TalphaR, TalphaI, TbetaR, TbetaI  ! P-wave, S-wave speeds (user units)
  REAL     (KIND=8) :: Tbeta, TfT                      ! power law and transition frequency
  COMPLEX  (KIND=8) :: TcP, TcS                        ! P-wave, S-wave speeds (neper/m loss)
  REAL     (KIND=8) :: Trho, TBumpDensity, Teta, Txi     ! density, boss parameters

  CHARACTER (LEN=1) :: BBC                            ! Boundary condition type
  REAL     (KIND=8) :: BalphaR, BalphaI, BbetaR, BbetaI  ! P-wave, S-wave speeds (user units)
  REAL     (KIND=8) :: Bbeta, BfT                      ! power law and transition frequency
  COMPLEX  (KIND=8) :: BcP, BcS                        ! P-wave, S-wave speeds (neper/m loss)
  REAL     (KIND=8) :: Brho, BBumpDensity, Beta, Bxi     ! density, boss parameters


  SUBROUTINE ReadSSP( Medium, N1 )

    ! reads the SSP data from the environmental file for a given medium
    INTEGER, INTENT( IN    ) :: Medium
    INTEGER, INTENT( INOUT ) :: N1
    INTEGER                  :: iSSP

    SNPts( Medium ) = N1

    ! The variable SSP%Loc( Medium ) points to the starting point for the
    ! data in the arrays z, alpha, beta and rho
    IF ( Medium == 1 ) THEN
       SLoc( Medium ) = 0
    ELSE
       SLoc( Medium ) = SLoc( Medium - 1 ) + SNPts( Medium - 1 )
    ENDIF
    ILoc = SLoc( Medium )

    !  Read in data and convert attenuation to Nepers/m
    N1 = 1
    DO iSSP = 1, MaxSSP
       iz = SLoc( Medium ) + iSSP

       READ(  ENVFile, *    ) Sz( iz ), alphaR, betaR, rhoR, alphaI, betaI

       ! Verify that the depths are monotonically increasing
       IF ( iSSP > 1 ) THEN
          IF ( Sz( iz ) .LE. Sz( iz - 1 ) ) THEN
              CALL ERROUT( 'ReadSSP', 'The depths in the SSP must be monotonically increasing' )
          END IF
       END IF

       SalphaR( iz ) = alphaR
       SalphaI( iz ) = alphaI
       Srho(    iz ) = rhoR
       SbetaR(  iz ) = betaR
       SbetaI(  iz ) = betaI

       ! Did we read the last point?
       IF ( ABS( Sz( iz ) - SDepth( Medium + 1 ) ) < 100. * EPSILON( 1.0e0 ) ) THEN
          SNPts( Medium ) = N1
          IF ( Medium == 1 ) SDepth( 1 ) = Sz( 1 )
          IF ( SNPts( Medium ) == 1 ) THEN
              CALL ERROUT( 'ReadSSP', 'The SSP must have at least 2 points in each layer' )
          END IF

          RETURN
       ENDIF

       N1 = N1 + 1
    END DO

    ! Fall through means too many points in the profile
    CALL ERROUT( 'ReadSSP', 'Number of SSP points exceeds limit' )

  END SUBROUTINE ReadSSP
