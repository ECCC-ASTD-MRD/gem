!     ######################################################################
      SUBROUTINE CONVECT_TSTEP_PREF( KLON, KLEV,                           &
                                   & PU, PV, PPRES, PZ, PDXDY, KLCL, KCTL, &
                                   & PTIMEA, PPREF )
!     ######################################################################

!!**** Routine to compute convective advection time step and precipitation
!!     efficiency
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      advection time step PTIMEC as a function of the mean ambient
!!      wind as well as the precipitation efficiency as a function
!!      of wind shear and cloud base height.
!!
!!
!!**  METHOD
!!    ------
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
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
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

USE YOE_CONVPAREXT
IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :

integer, INTENT(IN)                    :: KLON   ! horizontal dimension
integer, INTENT(IN)                    :: KLEV   ! vertical dimension
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PU     ! grid scale horiz. wind u
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PV     ! grid scale horiz. wind v
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m)
real, DIMENSION(KLON),      INTENT(IN) :: PDXDY  ! grid area (m^2)
integer, DIMENSION(KLON),   INTENT(IN) :: KLCL   ! lifting condensation level index
integer, DIMENSION(KLON),   INTENT(IN) :: KCTL   ! cloud top level index

real, DIMENSION(KLON),      INTENT(OUT):: PTIMEA ! advective time period
real, DIMENSION(KLON),      INTENT(OUT):: PPREF  ! precipitation efficiency


!*       0.2   Declarations of local variables KLON

integer :: IIE, IKB, IKE                      ! horizontal + vertical loop bounds
integer :: JI                                 ! horizontal loop index
integer :: JK, JKLC, JKP5, JKCT               ! vertical loop index

integer, DIMENSION(KLON)  :: IP500       ! index of 500 hPa levels
real, DIMENSION(KLON)     :: ZCBH        ! cloud base height
real, DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3  ! work arrays


!-------------------------------------------------------------------------------

!        0.3   Set loop bounds
!              ---------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Determine vertical index for 500 hPa levels
!               ------------------------------------------


IP500(:) = IKB
DO JK = IKB, IKE
    WHERE ( PPRES(:,JK) >= 500.E2 ) IP500(:) = JK
ENDDO


!*       2.     Compute convective time step
!               ----------------------------

            ! compute wind speed at LCL, 500 hPa, CTL

DO JI = 1, IIE
   JKLC = KLCL(JI)
   JKP5 = IP500(JI)
   JKCT = KCTL(JI)
   ZWORK1(JI) = SQRT( PU(JI,JKLC) * PU(JI,JKLC) +           &
              &       PV(JI,JKLC) * PV(JI,JKLC)  )
   ZWORK2(JI) = SQRT( PU(JI,JKP5) * PU(JI,JKP5) +           &
              &       PV(JI,JKP5) * PV(JI,JKP5)  )
   ZWORK3(JI) = SQRT( PU(JI,JKCT) * PU(JI,JKCT) +           &
              &       PV(JI,JKCT) * PV(JI,JKCT)  )
ENDDO

ZWORK2(:) = MAX( 0.1, 0.5 * ( ZWORK1(:) + ZWORK2(:) ) )

PTIMEA(:) = SQRT( PDXDY(:) ) / ZWORK2(:)


!*       3.     Compute precipitation efficiency
!               -----------------------------------

!*       3.1    Precipitation efficiency as a function of wind shear
!               ----------------------------------------------------

ZWORK2(:) = SIGN( 1., ZWORK3(:) - ZWORK1(:) )
DO JI = 1, IIE
    JKLC = KLCL(JI)
    JKCT = KCTL(JI)
    ZWORK1(JI) = ( PU(JI,JKCT) - PU(JI,JKLC) )  *          &
               & ( PU(JI,JKCT) - PU(JI,JKLC) )  +          &
               & ( PV(JI,JKCT) - PV(JI,JKLC) )  *          &
               & ( PV(JI,JKCT) - PV(JI,JKLC) )
    ZWORK1(JI) = 1.E3 * ZWORK2(JI) * SQRT( ZWORK1(JI) ) /  &
               & MAX( 1.E-2, PZ(JI,JKCT) - PZ(JI,JKLC) )
ENDDO

PPREF(:)  = 1.591 + ZWORK1(:) * ( -.639 + ZWORK1(:)       &
          &              * (  9.53E-2 - ZWORK1(:) * 4.96E-3 ) )
PPREF(:)  = MAX( .4, MIN( PPREF(:), .70 ) )

!*       3.2    Precipitation efficiency as a function of cloud base height
!               ----------------------------------------------------------

DO JI = 1, IIE
   JKLC = KLCL(JI)
   ZCBH(JI)   = MAX( 3., ( PZ(JI,JKLC) - PZ(JI,IKB) ) * 3.281E-3 )
ENDDO
ZWORK1(:) = .9673 + ZCBH(:) * ( -.7003 + ZCBH(:) * ( .1622 + &
          &   ZCBH(:) *  ( -1.2570E-2 + ZCBH(:) * ( 4.2772E-4 -   &
          &   ZCBH(:) * 5.44E-6 ) ) ) )
ZWORK1(:) = MAX( .4, MIN( .70, _ONE_/ ( _ONE_ + ZWORK1(:) ) ) )

!*       3.3    Mean precipitation efficiency is used to compute rainfall
!               ----------------------------------------------------------

PPREF(:) = 0.5 * ( PPREF(:) + ZWORK1(:) )


END SUBROUTINE CONVECT_TSTEP_PREF

