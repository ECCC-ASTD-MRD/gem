!     ######################################################################
      SUBROUTINE CONVECT_PRECIP_ADJUST( KLON, KLEV,                        &
                                      & PPRES, PUMF, PUER, PUDR,           &
                                      & PUPR, PUTPR, PURW,                 &
                                      & PDMF, PDER, PDDR, PDTHL, PDRW,     &
                                      & PPREF, PTPR, PMIXF, PDTEVR,        &
                                      & KLFS, KDBL, KLCL, KCTL, KETL,      &
                                      & PDTEVRF )
!     ######################################################################

!!**** Adjust up- and downdraft mass fluxes to be consistent with the
!!     mass transport at the LFS given by the precipitation efficiency
!!     relation.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust up- and downdraft mass
!!      fluxes below the LFS to be consistent with the precipitation
!!      efficiency relation
!!
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
!!     Module YOE_CONVPAR
!!        XUSRDPTH             ! pressure depth to compute updraft humidity
!!                             ! supply rate for downdraft
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_PRECIP_ADJUST)
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
USE YOE_CONVPAR

IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :


integer,                    INTENT(IN) :: KLON  ! horizontal dimension
integer,                    INTENT(IN) :: KLEV  ! vertical dimension
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
real, DIMENSION(KLON),      INTENT(IN) :: PUTPR ! updraft  total precipit. (kg/s
real, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
real, DIMENSION(KLON),      INTENT(IN) :: PMIXF ! critical mixed fraction at LCL
integer, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
integer, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
integer, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of equilibrium
                                                  ! (zero buoyancy) level
integer, DIMENSION(KLON),  INTENT(INOUT) :: KLFS ! contains vert. index of LFS
integer, DIMENSION(KLON),  INTENT(INOUT) :: KDBL ! contains vert. index of DBL

real, DIMENSION(KLON),      INTENT(INOUT) :: PDTEVR ! total downdraft evaporation
                                                      ! rate at LFS
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTEVRF! downdraft evaporation rate
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER   ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   ! updraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUPR   ! updraft  precipit. (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER   ! downdraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   ! downdraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTHL  ! downdraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDRW   ! downdraft total water (kg/kg)

real, DIMENSION(KLON),     INTENT(OUT)   :: PTPR    ! total precipitation (kg/s)
                                                      ! = downdraft precipitation

!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE        ! horizontal + vertical loop bounds
integer :: JK, JKT1, JKT2, JKT3 ! vertical loop index
integer :: JI                   ! horizontal loop index

integer, DIMENSION(KLON) :: IPRL
real, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, &
                           &  ZWORK4, ZWORK5, ZWORK6    ! work arrays


!-------------------------------------------------------------------------------

!        0.3   Set loop bounds
!              ---------------

IKB  = 1 + JCVEXB
IKE  = KLEV - JCVEXT
IIE  = KLON
JKT1 = MAXVAL( KLFS(:) )
JKT2 = MAXVAL( KCTL(:) )
JKT3 = MINVAL( KLCL(:) )


!        1.    Set some output variables for columns where no downdraft
!              exists. Exit if there is no downdraft at all.
!              --------------------------------------------------------

IPRL(:) = IKB
PTPR(:) = _ZERO_

WHERE ( PDTEVR(:) == _ZERO_ )
     PTPR(:)    = PUTPR(:)  ! no downdraft evaporation => no downdraft, all
                            ! precipitation occurs in updraft
END WHERE
IF ( COUNT( PDTEVR(:) > _ZERO_ ) == 0 ) RETURN ! exit routine if no downdraft exists

!*       2.     The total mass transported from the updraft to the down-
!               draft at the LFS must be consistent with the three water
!               budget terms :
!               ---------------------------------------------------------

!*       2.1    Downdraft evaporation rate at the DBL. The evaporation
!               rate in downdraft must be consistent with precipitation
!               efficiency relation.
!               --------------------------------------------------------


DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI) = PDTEVR(JI) / MIN( -1.E-1, PDMF(JI,JK) )
     ZWORK6(JI) = PDMF(JI,JK)
ENDDO

!*       2.2    Some preliminar computations for downdraft = total
!               precipitation rate. The precipitation is evaluated in
!               a layer thickness DP=XUSRDPTH=165 hPa above the LCL.
!               The difference between updraft precipitation and downdraft
!               precipitation (updraft supply rate) is used to drive the
!               downdraft through evaporational cooling.
!               --------------------------------------------------------

DO JI = 1, IIE
     JK = KLCL(JI)
     ZWORK5(JI) = PPRES(JI,JK)
ENDDO

PTPR(:) = _ZERO_
DO JK = JKT3, JKT2
    WHERE ( JK >= KLCL(:) .AND. PPRES(:,JK) >= ZWORK5(:) - XUSRDPTH )
        PTPR(:) = PTPR(:) + PUPR(:,JK)
        IPRL(:) = JK
    END WHERE
ENDDO
IPRL(:) = MIN( KETL(:), IPRL(:) )

DO JI = 1, IIE
     JK = IPRL(JI)
     PTPR(JI) = PUMF(JI,JK+1) * PURW(JI,JK+1) + PTPR(JI)
ENDDO

PTPR(:)   = PPREF(:) * MIN( PUTPR(:), PTPR(:) )
ZWORK4(:) = PUTPR(:) - PTPR(:)


!*       2.3    Total amount of precipitation that falls out of the up-
!               draft between the LCL and the LFS.
!               Condensate transfer from up to downdraft at LFS
!               ---------------------------------------------------------

ZWORK5(:) = _ZERO_
DO JK = JKT3, JKT1
     WHERE ( JK >= KLCL(:) .AND. JK <= KLFS(:) )
           ZWORK5(:) = ZWORK5(:) +  PUPR(:,JK)
     END WHERE
ENDDO

DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK2(JI) = ( _ONE_ - PPREF(JI) ) * ZWORK5(JI) *                     &
                & ( _ONE_ - PMIXF(JI) ) / MAX( 1.E-1, PUMF(JI,JK) )
ENDDO


!*       2.4    Increase the first guess downdraft mass flux to satisfy
!               precipitation efficiency relation.
!               If downdraft does not evaporate any water at the DBL for
!               the specified relative humidity, or if the corrected mass
!               flux at the LFS is positive no downdraft is allowed
!               ---------------------------------------------------------


ZWORK1(:) = ZWORK4(:) / ( ZWORK1(:) + ZWORK2(:) + 1.E-8 )
ZWORK2(:) = ZWORK1(:) / MIN( -1.E-1, ZWORK6(:) ) ! ratio of budget consistent to actual DMF

ZWORK3(:) = _ONE_
ZWORK6(:) = _ONE_
WHERE ( ZWORK1(:) > _ZERO_ .OR. PDTEVR(:) < _ONE_ )
   KDBL(:)   = IKB
   KLFS(:)   = IKB
   PDTEVR(:) = _ZERO_
   ZWORK2(:) = _ZERO_
   ZWORK3(:) = _ZERO_
   ZWORK6(:) = _ZERO_
END WHERE

DO JK = IKB, JKT1
     PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
     PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:)
     PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:)
   PDTEVRF(:,JK) = PDTEVRF(:,JK)* ZWORK2(:)
     PDRW(:,JK)  = PDRW(:,JK)  * ZWORK3(:)
     PDTHL(:,JK) = PDTHL(:,JK) * ZWORK3(:)
ENDDO
ZWORK4(:) = ZWORK2(:)


!*       3.     Increase updraft mass flux, mass detrainment rate, and water
!               substance detrainment rates to be consistent with the transfer
!               of the estimated mass from the up- to the downdraft at the LFS
!               --------------------------------------------------------------

DO JI = 1, IIE
    JK = KLFS(JI)
    ZWORK2(JI) = ( _ONE_ - ZWORK6(JI) ) + ZWORK6(JI) *                   &
               &  ( PUMF(JI,JK) - ( _ONE_ - PMIXF(JI) ) * ZWORK1(JI) ) / &
               &  MAX( 1.E-1, PUMF(JI,JK) )
ENDDO


JKT1  = MAXVAL( KLFS(:) )  ! value of KLFS might have been reset to IKB above
DO JK = IKB, JKT1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) ) THEN
        PUMF(JI,JK)  = PUMF(JI,JK)  * ZWORK2(JI)
        PUER(JI,JK)  = PUER(JI,JK)  * ZWORK2(JI)
        PUDR(JI,JK)  = PUDR(JI,JK)  * ZWORK2(JI)
        PUPR(JI,JK)  = PUPR(JI,JK)  * ZWORK2(JI)
      ENDIF
    ENDDO
ENDDO


!*       4.     Increase total = downdraft precipitation and evaporation rate
!               -------------------------------------------------------------

WHERE ( PDTEVR(:) > _ZERO_ )
    PDTEVR(:)  = PDTEVR(:) * ZWORK4(:)
    PTPR(:)    = PTPR(:) + PPREF(:) * ZWORK5(:) * ( ZWORK2(:) - _ONE_ )
ELSEWHERE
    PTPR(:)    = PUTPR(:)
END WHERE


END SUBROUTINE CONVECT_PRECIP_ADJUST

