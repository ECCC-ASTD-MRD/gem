!     ##########################################################################
      SUBROUTINE CONVECT_CONDENS( KLON,                                        &
                            &  KICE, PPRES, PTHL, PRW, PRCO, PRIO, PZ, OWORK1, &
                            &  PT, PEW, PRC, PRI, PLV, PLS, PCPH   )
!     ##########################################################################

!!**** Compute temperature cloud and ice water content from enthalpy and r_w
!!
!!
!!    PURPOSE
!!    -------
!!     The purpose of this routine is to determine cloud condensate
!!     and to return values for L_v, L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!     Condensate is extracted iteratively
!!
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module YOMCST
!!          RG                   ! gravity constant
!!          RALPW, RBETW, RGAMW ! constants for water saturation pressure
!!          RALPS, RBETS, RGAMS ! constants for ice saturation pressure
!!          RATM                 ! reference pressure
!!          RD, RV             ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV           ! specific heat for dry air and water vapor
!!          RCW, RCS             ! specific heat for liquid water and ice
!!          RTT                  ! triple point temperature
!!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR
!!          XTFRZ1               ! begin of freezing interval
!!          XTFRZ2               ! end of freezing interval
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CONDENS)
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

IMPLICIT NONE
#include <arch_specific.hf>
#define _ZERO_   0.0
#define _ONE_    1.0

!*       0.1   Declarations of dummy arguments :

integer, INTENT(IN)                :: KLON    ! horizontal loop index
integer, INTENT(IN)                :: KICE    ! flag for ice ( 1 = yes,
                                              !                0 = no ice )
real, DIMENSION(KLON),   INTENT(IN) :: PPRES  ! pressure
real, DIMENSION(KLON),   INTENT(IN) :: PTHL   ! enthalpy (J/kg)
real, DIMENSION(KLON),   INTENT(IN) :: PRW    ! total water mixing ratio
real, DIMENSION(KLON),   INTENT(IN) :: PRCO   ! cloud water estimate (kg/kg)
real, DIMENSION(KLON),   INTENT(IN) :: PRIO   ! cloud ice   estimate (kg/kg)
real, DIMENSION(KLON),   INTENT(IN) :: PZ     ! level height (m)
LOGICAL, DIMENSION(KLON),INTENT(IN) :: OWORK1 ! logical mask


real, DIMENSION(KLON),   INTENT(OUT):: PT     ! temperature
real, DIMENSION(KLON),   INTENT(OUT):: PRC    ! cloud water mixing ratio(kg/kg)
real, DIMENSION(KLON),   INTENT(OUT):: PRI    ! cloud ice mixing ratio  (kg/kg)
real, DIMENSION(KLON),   INTENT(OUT):: PLV    ! latent heat L_v
real, DIMENSION(KLON),   INTENT(OUT):: PLS    ! latent heat L_s
real, DIMENSION(KLON),   INTENT(OUT):: PCPH   ! specific heat C_ph
real, DIMENSION(KLON),   INTENT(OUT):: PEW    ! water saturation mixing ratio

!*       0.2   Declarations of local variables KLON

integer :: JITER          ! iteration index
real    :: ZEPS, ZEPSA    ! R_d / R_v, 1 / ZEPS
real    :: ZCVOCD         ! RCPV / RCPD
real    :: ZRDOCP         ! R_d / C_pd

real, DIMENSION(KLON)    :: ZEI           ! ice saturation mixing ratio
real, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZT ! work arrays


!-------------------------------------------------------------------------------

!*       1.     Initialize temperature and Exner function
!               -----------------------------------------

ZRDOCP      = RD   / RCPD
ZEPS        = RD   / RV
ZEPSA       = _ONE_ / ZEPS
ZCVOCD      = RCPV  / RCPD


    ! Make a first temperature estimate, based e.g. on values of
    !  r_c and r_i at lower level

      !! Note that the definition of ZCPH is not the same as used in
      !! routine CONVECT_SATMIXRATIO
     PCPH(:)   = RCPD + RCPV * PRW(:)
     ZWORK1(:) = ( _ONE_ + PRW(:) ) * RG * PZ(:)
     PT(:)     = ( PTHL(:) + PRCO(:) * RLVTT + PRIO(:) * RLSTT - ZWORK1(:) )   &
               & / PCPH(:)
     PT(:)     = MAX(180., MIN( 330., PT(:) ) ) ! set overflow bounds in
                                                          ! case that PTHL=0


!*       2.     Enter the iteration loop
!               ------------------------

DO JITER = 1,6
     PEW(:) = EXP( RALPW - RBETW / PT(:) - RGAMW * LOG( PT(:) ) )
     ZEI(:) = EXP( RALPS - RBETS / PT(:) - RGAMS * LOG( PT(:) ) )
     PEW(:) = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
     ZEI(:) = ZEPS * ZEI(:) / ( PPRES(:) - ZEI(:) )

     PLV(:)    = RLVTT + ( RCPV - RCW ) * ( PT(:) - RTT ) ! compute L_v
     PLS(:)    = RLSTT + ( RCPV - RCS ) * ( PT(:) - RTT ) ! compute L_i

     ZWORK2(:) = ( PT(:) - XTFRZ2 ) / ( XTFRZ1 - XTFRZ2 ) ! freezing interval
     ZWORK2(:) = MAX( _ZERO_, MIN(_ONE_, ZWORK2(:) ) ) 
     IF ( KICE==0 ) ZWORK2(:) = _ONE_
     ZWORK2(:) = ZWORK2(:) * ZWORK2(:)
     ZWORK3(:) = ( _ONE_ - ZWORK2(:) ) * ZEI(:) + ZWORK2(:) * PEW(:)
     PRC(:)    = MAX( _ZERO_, ZWORK2(:) * ( PRW(:) - ZWORK3(:) ) )
     PRI(:)    = MAX( _ZERO_, ( _ONE_ - ZWORK2(:) ) * ( PRW(:) - ZWORK3(:) ) )
     ZT(:)     = ( PTHL(:) + PRC(:) * PLV(:) + PRI(:) * PLS(:) - ZWORK1(:) )   &
               & / PCPH(:)
     PT(:) = PT(:) + ( ZT(:) - PT(:) ) * 0.4  ! force convergence
     PT(:) = MAX( 175., MIN( 330., PT(:) ) )
ENDDO


END SUBROUTINE CONVECT_CONDENS

