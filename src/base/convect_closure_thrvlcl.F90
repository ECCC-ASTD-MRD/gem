!######################################################################
 SUBROUTINE CONVECT_CLOSURE_THRVLCL( KLON, KLEV,                      &
                              &  PPRES, PTH, PRV, PZ, OWORK1,         &
                              &  PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
                              &  KLCL, KDPL, KPBL )
!######################################################################

!!**** Determine thermodynamic properties at new LCL
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the thermodynamic
!!      properties at the new lifting condensation level LCL
!!
!!
!!
!!**  METHOD
!!    ------
!!    see CONVECT_TRIGGER_FUNCT
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! Reference pressure
!!          RD, RV           ! Gaz  constants for dry air and water vapor
!!          RCPD               ! Cpd (dry air)
!!          RTT                ! triple point temperature
!!          RBETW, RGAMW      ! constants for vapor saturation pressure
!!
!!      Module YOE_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! lowest allowed pressure difference between
!!                             ! surface and LCL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!
!!      Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAR
USE YOE_CONVPAREXT

IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _ONE_    1.0

!*       0.1   Declarations of dummy arguments :

integer,                    INTENT(IN) :: KLON  ! horizontal dimension
integer,                    INTENT(IN) :: KLEV  ! vertical dimension
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! theta
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! vapor mixing ratio
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of grid point (m)
integer, DIMENSION(KLON),   INTENT(IN) :: KDPL  ! contains vert. index of DPL
integer, DIMENSION(KLON),   INTENT(IN) :: KPBL  ! " vert. index of source layer top
LOGICAL,   DIMENSION(KLON),   INTENT(IN) :: OWORK1! logical mask

real, DIMENSION(KLON),     INTENT(OUT):: PTHLCL ! theta at LCL
real, DIMENSION(KLON),     INTENT(OUT):: PRVLCL ! vapor mixing ratio at  LCL
real, DIMENSION(KLON),     INTENT(OUT):: PZLCL  ! height at LCL (m)
real, DIMENSION(KLON),     INTENT(OUT):: PTLCL  ! temperature at LCL (m)
real, DIMENSION(KLON),     INTENT(OUT):: PTELCL ! environm. temp. at LCL (K)
integer, DIMENSION(KLON),  INTENT(OUT):: KLCL   ! contains vert. index of LCL

!*       0.2   Declarations of local variables :

integer :: JK, JKM                    ! vertical loop index
integer :: JI                         ! horizontal loop index
integer :: IIE, IKB, IKE              ! horizontal + vertical loop bounds
real    :: ZEPS, ZEPSA    ! R_d / R_v, R_v / R_d
real    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd

real, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
real, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
real, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure
real, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
real, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
real, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
real, DIMENSION(KLON) :: ZWORK1, ZWORK2     ! work arrays


!-------------------------------------------------------------------------------

!*       0.3    Compute array bounds
!               --------------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Initialize local variables
!               --------------------------

ZEPS      = RD  / RV
ZEPSA     = RV  / RD
ZCPORD    = RCPD / RD
ZRDOCP    = RD  / RCPD

ZDPTHMIX(:) = _ZERO_
ZPRESMIX(:) = _ZERO_
PTHLCL(:)   = 300.
PTLCL(:)    = 300.
PTELCL(:)   = 300.
PRVLCL(:)   = _ZERO_
PZLCL(:)    = PZ(:,IKB)
ZTMIX(:)    = 230.
ZPLCL(:)    = 1.E4
KLCL(:)     = IKB + 1


!*       2.     Construct a mixed layer as in TRIGGER_FUNCT
!               -------------------------------------------

     DO JI = 1, IIE
        DO JK = KDPL(JI), KPBL(JI)
           JKM = JK + 1
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            PTHLCL(JI)   = PTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            PRVLCL(JI)   = PRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
        ENDDO
     ENDDO


WHERE ( OWORK1(:) )

        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        PTHLCL(:)   = PTHLCL(:)   / ZDPTHMIX(:)
        PRVLCL(:)   = PRVLCL(:)   / ZDPTHMIX(:)

!*       3.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               --------------------------------------------------


        ZTMIX(:)  = PTHLCL(:) * ( ZPRESMIX(:) / RATM ) ** ZRDOCP
        ZEVMIX(:) = PRVLCL(:) * ZPRESMIX(:) / ( PRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3 )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / &
                  & ( 17.502 - ZWORK1(:) )
              ! adiabatic saturation temperature
        PTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - RTT )   &
                  & - 4.36E-4 * ( ZTMIX(:) - RTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        PTLCL(:)  = MIN( PTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = RATM * ( PTLCL(:) / PTHLCL(:) ) ** ZCPORD

END WHERE

     ZPLCL(:) = MIN( 2.E5, MAX( 10., ZPLCL(:) ) ) ! bound to avoid overflow


!*       3.2    Correct PTLCL in order to be completely consistent
!               with MNH saturation formula
!               --------------------------------------------------

     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, PTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / PTLCL(:) * ( RBETW / PTLCL(:) - RGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) / &
                  &  ( _ONE_ + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        PTLCL(:)  = PTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)

     END WHERE


!*       3.3    If PRVLCL is oversaturated set humidity and temperature
!               to saturation values.
!               -------------------------------------------------------

    CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) .AND. PRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( RBETW / ZTMIX(:) - RGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( PRVLCL(:)- ZWORK1(:) ) / &
                  &  ( _ONE_ + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        PTLCL(:)  = ZTMIX(:) + ZLV(:) / ZCPH(:) * ZWORK2(:)
        PRVLCL(:) = PRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        PTHLCL(:) = PTLCL(:) * ( RATM / ZPLCL(:) ) ** ZRDOCP
     END WHERE


!*        4.1   Determine  vertical loop index at the LCL
!               -----------------------------------------

     DO JI = 1, IIE
        DO JK = KDPL(JI), IKE - 1
        IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. OWORK1(JI) ) THEN
            KLCL(JI)  = JK + 1
            PZLCL(JI) = PZ(JI,JK+1)
        ENDIF
        ENDDO
     ENDDO


!*        4.2   Estimate height and environmental temperature at LCL
!               ----------------------------------------------------

    DO JI = 1, IIE
        JK   = KLCL(JI)
        JKM  = JK - 1
        ZDP(JI)     = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /   &
                    & LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI)  = PTH(JI,JK)  * ( PPRES(JI,JK)  / RATM ) ** ZRDOCP
        ZWORK2(JI)  = PTH(JI,JKM) * ( PPRES(JI,JKM) / RATM ) ** ZRDOCP
        ZWORK1(JI)  = ZWORK2(JI) + ( ZWORK1(JI) - ZWORK2(JI) ) * ZDP(JI)
           ! we compute the precise value of the LCL
           ! The precise height is between the levels KLCL and KLCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    ENDDO
    WHERE( OWORK1(:) )
       PTELCL(:) = ZWORK1(:)
       PZLCL(:)  = ZWORK2(:)
    END WHERE



END SUBROUTINE CONVECT_CLOSURE_THRVLCL

