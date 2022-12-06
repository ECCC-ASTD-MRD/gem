!     #######################################################
      SUBROUTINE CONVECT_MIXING_FUNCT( KLON,                &
                                     & PMIXC, KMF, PER, PDR )
!     #######################################################

!!**** Determine the area under the distribution function
!!     KMF = 1 : gaussian  KMF = 2 : triangular distribution function
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the entrainment and
!!      detrainment rate by evaluating the are under the distribution
!!      function. The integration interval is limited by the critical
!!      mixed fraction PMIXC
!!
!!
!!
!!**  METHOD
!!    ------
!!      Use handbook of mathemat. functions by Abramowitz and Stegun, 1968
!!
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine MIXING_FUNCT)
!!      Abramovitz and Stegun (1968), handbook of math. functions
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

IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0

!*       0.1   Declarations of dummy arguments :

integer,               INTENT(IN) :: KLON   ! horizontal dimension
integer,               INTENT(IN) :: KMF    ! switch for dist. function
real, DIMENSION(KLON), INTENT(IN) :: PMIXC  ! critical mixed fraction

real, DIMENSION(KLON), INTENT(OUT):: PER    ! normalized entrainment rate
real, DIMENSION(KLON), INTENT(OUT):: PDR    ! normalized detrainment rate

!*       0.2   Declarations of local variables :

real    :: ZSIGMA = 0.166666667         ! standard deviation
real    :: ZFE    = 4.931813949         ! integral normalization

real    :: ZSQRTP = 2.506628 ,  ZP  = 0.33267   ! constants
real    :: ZA1    = 0.4361836 , ZA2 =-0.1201676 ! constants
real    :: ZA3    = 0.9372980 , ZT1 = 0.500498  ! constants
real    :: ZE45   = 0.01111                          ! constant

real, DIMENSION(KLON) :: ZX, ZY, ZW1, ZW2         ! work variables
real    :: ZW11


!-------------------------------------------------------------------------------

!       1.     Use gaussian function for KMF=1
!              -------------------------------

IF( KMF == 1 ) THEN
    ! ZX(:)  = ( PMIXC(:) - _HALF_ ) / ZSIGMA
      ZX(:)  = 6. * PMIXC(:) - 3.
      ZW1(:) = _ONE_ / ( _ONE_+ ZP * ABS ( ZX(:) ) )
      ZY(:)  = EXP( -_HALF_ * ZX(:) * ZX(:) )
      ZW2(:) = ZA1 * ZW1(:) + ZA2 * ZW1(:) * ZW1(:) +                   &
             &   ZA3 * ZW1(:) * ZW1(:) * ZW1(:)
      ZW11   = ZA1 * ZT1 + ZA2 * ZT1 * ZT1 + ZA3 * ZT1 * ZT1 * ZT1
ENDIF

WHERE ( KMF == 1 .AND. ZX(:) >= _ZERO_ )
        PER(:) = ZSIGMA * ( _HALF_ * ( ZSQRTP - ZE45 * ZW11              &
               & - ZY(:) * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )        &
               & - _HALF_ * ZE45 * PMIXC(:) * PMIXC(:)
        PDR(:) = ZSIGMA*( _HALF_ * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )    &
               & + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
               & - ZE45 * ( _HALF_ + _HALF_ * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE
WHERE ( KMF == 1 .AND. ZX(:) < _ZERO_ )
        PER(:) = ZSIGMA*( _HALF_ * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )    &
               & + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
               & - _HALF_ * ZE45 * PMIXC(:) * PMIXC(:)
        PDR(:) = ZSIGMA * ( _HALF_ * ( ZSQRTP - ZE45 * ZW11 - ZY(:)      &
               & * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )                &
               & - ZE45 * ( _HALF_ + _HALF_ * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE

      PER(:) = PER(:) * ZFE
      PDR(:) = PDR(:) * ZFE


!       2.     Use triangular function KMF=2
!              -------------------------------

!     not yet released


END SUBROUTINE CONVECT_MIXING_FUNCT

