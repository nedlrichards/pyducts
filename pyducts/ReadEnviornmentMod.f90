MODULE ReadEnvironmentMod

  ! mbp 12/2018, based on much older subroutine

  USE sspMod
  USE FatalError

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: ENVFile = 5

CONTAINS
  SUBROUTINE ReadEnvironment(FileRoot, Title, freq0, MaxMedium, TopOpt, NG, BotOpt, cLow, cHigh, RMax)

    ! Reads in the info in ENVFile

    !USE SubTabulate

    INTEGER,             INTENT( IN  ) :: MaxMedium
    INTEGER,             INTENT( OUT ) :: NG( * )
    REAL       (KIND=8), INTENT( OUT ) :: freq0, cLow, cHigh, RMax
    CHARACTER ( LEN=* ), INTENT( IN  ) :: FileRoot
    CHARACTER,           INTENT( OUT ) :: TopOpt*( * ), BotOpt*( * ), Title*( * )
    INTEGER             :: NElts, Medium, Nneeded, iostat
    REAL       (KIND=8) :: rho( MaxSSP ), c, deltaz
    COMPLEX    (KIND=8) :: cP( MaxSSP ), cS( MaxSSP )
    CHARACTER ( LEN=2 ) :: AttenUnit
    CHARACTER ( LEN=8 ) :: Task

    ! Open the environmental file
    OPEN( UNIT = ENVFile, FILE = TRIM( FileRoot ) // '.env', STATUS = 'OLD',     IOSTAT = iostat, ACTION = 'READ' )

    IF ( IOSTAT /= 0 ) THEN   ! successful open?
       CALL ERROUT( 'READIN', 'Unable to open the environmental file' )
    END IF

    ! set default values for the ssp (in case this is the nth profile after one has already been read)
    alphaR = 1500
    betaR  = 0
    alphaI = 0
    betaI  = 0
    rhoR   = 1

    NElts  = 0         ! this is a dummy variable, passed to profil during read of SSP

    READ(  ENVFile, *, END = 9999 ) Title( 9 : 80 )
    READ(  ENVFile, *, END = 9999  ) freq0
    READ(  ENVFile, *, END = 9999  ) SNMedia

    IF ( SNMedia > MaxMedium ) THEN
       CALL ERROUT( 'READIN', 'Too many Media' )
    ENDIF

    CALL ReadTopOpt( TopOpt, TBC, AttenUnit )
    ! read top BC
    CALL TopBot(TBC, TalphaR, TalphaI, TbetaR, TbetaI, Tbeta, Tft, TcP, TcS, Trho, TBumpDensity, Teta, Txi)

    !  *** Internal media ***

    MediumLoop: DO Medium = 1, SNMedia
       IF ( AttenUnit( 1 : 1 ) == 'm' ) THEN   ! this particular attenuation unit needs to get a power law and transition frequency
          READ(  ENVFile, *, END = 9999 ) NG( Medium ), Ssigma( Medium ), &
               SDepth( Medium + 1 ), Sbeta( Medium ), SfT( Medium )
       ELSE
          READ(  ENVFile, *, END = 9999 ) NG( Medium ), Ssigma( Medium ), SDepth( Medium + 1 )
       END IF

       !  Call EvaluateSSP to read in SSP
       Task = 'INIT'
       CALL EvaluateSSP( cP, cS, rho, Medium, NElts, freq0, Task )

       ! estimate number of points needed (can be a bad estimate if the layer has a big change in sound speed)
       c = alphar   ! this is the last sound speed value that was read
       IF ( betar > 0.0 ) c = betar     ! shear?
       deltaz  = c / freq0 / 20         ! default sampling: 20 points per wavelength
       Nneeded = INT( ( SDepth( Medium + 1 ) - SDepth( Medium ) ) / deltaz )
       Nneeded = MAX( Nneeded, 10 )     ! require a minimum of 10 points

       IF ( NG( Medium ) == 0 ) THEN    ! automatic calculation of f.d. mesh
          NG( Medium ) = Nneeded
       ELSE IF ( NG( Medium ) < Nneeded / 2 ) THEN
          CALL ERROUT( 'ReadEnvironment', 'Mesh is too coarse' )
       END IF

    END DO MediumLoop

    ! *** Bottom properties ***

    IF ( AttenUnit( 1 : 1 ) == 'm' ) THEN   ! this particular attenuation unit needs to get a power law and transition frequency
       READ( ENVFile, *, END = 9999 ) BotOpt( 1 : 8 ), Ssigma( SNMedia + 1 ), Bbeta, BfT
    ELSE
       READ( ENVFile, *, END = 9999 ) BotOpt( 1 : 8 ), Ssigma( SNMedia + 1 )
    END IF

    BBC = BotOpt( 1 : 1 )
    ! Read bottom BC
    CALL TopBot(BBC, BalphaR, BalphaI, BbetaR, BbetaI, Bbeta, Bft, BcP, BcS, Brho, BBumpDensity, Beta, Bxi)


    ! Read phase speed limits
    READ(  ENVFile, *    ) cLow, cHigh                 ! Spectral limits (m/s)
    IF ( cLow >= cHigh ) CALL ERROUT( 'GetPar', 'Need phase speeds cLow < cHigh'  )

    ! Read maximum range
    READ(  ENVFile, * ) RMax          ! Maximum range for calculations (km)
    IF ( RMax < 0.0 ) CALL ERROUT( ' ', 'RMax must be non-negative'  )

    RETURN

    CLOSE( ENVFile )
    STOP

  END SUBROUTINE ReadEnvironment

  !**********************************************************************!

  SUBROUTINE ReadTopOpt( TopOpt, BC, AttenUnit )

    USE AttenMod
    INTEGER, PARAMETER :: SSPFile = 40
    CHARACTER (LEN= 8), INTENT( OUT ) :: TopOpt
    CHARACTER (LEN= 1), INTENT( OUT ) :: BC                     ! Boundary condition type
    CHARACTER (LEN= 2), INTENT( OUT ) :: AttenUnit

    TopOpt = '      '   ! initialize to blanks
    READ( ENVFile, * ) TopOpt

    SType      = TopOpt( 1 : 1 )
    BC            = TopOpt( 2 : 2 )
    AttenUnit     = TopOpt( 3 : 4 )
    SAttenUnit = AttenUnit

    !  Added volume attenuation

    SELECT CASE ( AttenUnit( 2 : 2 ) )
    CASE ( 'F' )
       READ(  ENVFile, * ) T, Salinity, pH, z_bar
    CASE ( 'B' )
       CALL ERROUT( 'Bio Layers are not yet implimented' )
       READ( ENVFile, *  ) NBioLayers
       IF ( NBioLayers > MaxBioLayers ) THEN
          CALL ERROUT( 'READIN', 'Too many biolayers' )
       END IF

       !DO iBio = 1, NBioLayers
       !   READ( ENVFile, *  ) bio( iBio )%Z1, bio( iBio )%Z2, bio( iBio )%f0, bio( iBio )%Q, bio( iBio )%a0
       !END DO
    CASE ( ' ' )
    CASE DEFAULT
       CALL ERROUT( 'ReadTopOpt', 'Unknown top option letter in fourth position' )
    END SELECT

  END SUBROUTINE ReadTopOpt

  !**********************************************************************!
  SUBROUTINE TopBot(BC, alphaR, alphaI, betaR, betaI, beta, ft, cP, cS, rho, BumpDensity, eta, xi)

    ! Handles top and bottom boundary conditions

    ! Input:
    !     HS%BC:   Boundary condition type
    !
    ! Output:
    !    HS%cP:    P-wave speed in halfspace
    !    HS%cS:    S-wave speed in halfspace
    !    HS%rho:   density in halfspace

    TYPE( HSInfo ), INTENT( INOUT ) :: HS
    REAL     (KIND=8) :: zTemp

    ! Echo to PRTFile user's choice of boundary condition
    SELECT CASE ( HS%BC )
    CASE ( 'V' )
       WRITE( PRTFile, * ) '    VACUUM'
    CASE ( 'R' )
       WRITE( PRTFile, * ) '    Perfectly RIGID'
    CASE ( 'A' )
       WRITE( PRTFile, * ) '    ACOUSTO-ELASTIC half-space'
    CASE ( 'F' )
       WRITE( PRTFile, * ) '    FILE used for reflection loss'
    CASE ( 'W' )
       WRITE( PRTFile, * ) '    Writing an IRC file'
    CASE ( 'P' )
       WRITE( PRTFile, * ) '    reading PRECALCULATED IRC'
    CASE DEFAULT
       CALL ERROUT( 'TopBot', 'Unknown boundary condition type' )
    END SELECT

    ! Read in BC parameters depending on particular choice
    cP  = 0.0
    cS  = 0.0
    rho = 0.0

    SELECT CASE ( HS%BC )
    CASE ( 'A' )                   !  Half-space properties
       zTemp = 0.0
       READ(  ENVFile, *    ) zTemp, alphaR, betaR, rhoR, alphaI, betaI
       IF ( alphaR == 0.0 .OR. rhoR == 0.0 ) &
            CALL ERROUT( 'TopBot', 'Sound speed or density vanishes in halfspace' )
    END SELECT

    RETURN
  END SUBROUTINE TopBot
END MODULE ReadEnvironmentMod
