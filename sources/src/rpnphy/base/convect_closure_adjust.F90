!    #########################################################################
     SUBROUTINE CONVECT_CLOSURE_ADJUST( KLON, KLEV, PADJ,                    &
                                  &   PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR, &
                                  &   PDMF, PZDMF, PDER, PZDER, PDDR, PZDDR, &
                                  &   PPRMELT, PZPRMELT, PDTEVR, PZDTEVR,    &
                                  &   PTPR, PZTPR,                           &
                                  &   PPRLFLX, PZPRLFL, PPRSFLX, PZPRSFL     )
!    #########################################################################

!!**** Uses closure adjustment factor to adjust mass flux and to modify
!!     precipitation efficiency  when necessary. The computations are
!!     similar to routine CONVECT_PRECIP_ADJUST.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust the mass flux using the
!!      factor PADJ computed in CONVECT_CLOSURE
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!    EXTERNAL
!!    --------
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
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


integer,                    INTENT(IN) :: KLON     ! horizontal dimension
integer,                    INTENT(IN) :: KLEV     ! vertical dimension
real, DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor


real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF  ! downdraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDMF ! initial value of  "
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER  ! downdraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDER ! initial value of  "
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR  ! downdraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDDR ! initial value of  "
real, DIMENSION(KLON),   INTENT(INOUT):: PTPR     ! total precipitation (kg/s)
real, DIMENSION(KLON),   INTENT(INOUT):: PZTPR    ! initial value of "
real, DIMENSION(KLON),   INTENT(INOUT):: PDTEVR   ! donwndraft evapor. (kg/s)
real, DIMENSION(KLON),   INTENT(INOUT):: PZDTEVR  ! initial value of "
real, DIMENSION(KLON),   INTENT(INOUT):: PPRMELT  ! melting of precipitation
real, DIMENSION(KLON),   INTENT(INOUT):: PZPRMELT ! initial value of "
real, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRLFLX! liquid precip flux
real, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRLFL! initial value "
real, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRSFLX! solid  precip flux
real, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRSFL! initial value "


!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE                 ! horiz. + vert. loop bounds
integer :: JK                            ! vertical loop index


!-------------------------------------------------------------------------------

!*       0.3   Compute loop bounds
!              -------------------

IIE  = KLON
IKB  = 1 + JCVEXB
IKE  = KLEV - JCVEXT


!*       1.     Adjust mass flux by the factor PADJ to converge to
!               specified degree of stabilization
!               ----------------------------------------------------

          PPRMELT(:)  = PZPRMELT(:)   * PADJ(:)
          PDTEVR(:)   = PZDTEVR(:)    * PADJ(:)
          PTPR(:)     = PZTPR(:)      * PADJ(:)

     DO JK = IKB + 1, IKE
          PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
          PDMF(:,JK)  = PZDMF(:,JK)   * PADJ(:)
          PDER(:,JK)  = PZDER(:,JK)   * PADJ(:)
          PDDR(:,JK)  = PZDDR(:,JK)   * PADJ(:)
          PPRLFLX(:,JK) = PZPRLFL(:,JK) * PADJ(:)
          PPRSFLX(:,JK) = PZPRSFL(:,JK) * PADJ(:)
     ENDDO

END SUBROUTINE CONVECT_CLOSURE_ADJUST

