!#######################################################################
  SUBROUTINE CONVECT_CHEM_TRANSPORT( KLON, KLEV, KCH, PCH1, PCH1C,     &
                                 & KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                 & PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                 & PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                 & KFTSTEPS )
!#######################################################################

!!**** Compute  modified chemical tracer values due to convective event
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!      environmental values of the chemical tracers
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      main deep convection code
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    11/12/97
!!
!-------------------------------------------------------------------------------

!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAREXT

IMPLICIT NONE

#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0

!*       0.1   Declarations of dummy arguments :

integer,                INTENT(IN) :: KLON     ! horizontal dimension
integer,                INTENT(IN) :: KLEV     ! vertical dimension
integer,                INTENT(IN) :: KCH      ! number of passive tracers

real,DIMENSION(KLON,KLEV,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
real,DIMENSION(KLON,KLEV,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.

integer, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
integer, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
integer, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
integer, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
integer, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
integer, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level

real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)

real, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
real, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
real, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
integer,                   INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps

!*       0.2   Declarations of local variables :

integer :: INCH1          ! number of chemical tracers
integer :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
integer :: IKS            ! vertical dimension
integer :: JI             ! horizontal loop index
integer :: JK, JKP        ! vertical loop index
integer :: JN             ! chemical tracer loop index
integer :: JSTEP          ! fractional time loop index
integer :: JKLD, JKLP, JKMAX ! loop index for levels
integer, DIMENSION(KLON) :: ITSTEP !fractional convective time step

real, DIMENSION(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
real, DIMENSION(KLON,KLEV,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
real, DIMENSION(KLON)          :: ZTIMEC! adjust fractional convective time step
real, DIMENSION(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
real, DIMENSION(KLON,KLEV,KCH) :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
real, DIMENSION(KLON,KCH)      :: ZWORK1, ZWORK2, ZWORK3

!-------------------------------------------------------------------------------

!*       0.3   Compute loop bounds
!              -------------------

INCH1  = KCH
IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
!JKMAX  = MAXVAL( KCTL(:) ) -> KCTL(JI)


!*      2.      Updraft computations
!               --------------------

ZUCH1(:,:,:) = _ZERO_

!*      2.1     Initialization  at LCL
!               ----------------------------------

DO JI = 1, IIE
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,:) = _HALF_ * ( PCH1(JI,JKLD,:) + PCH1(JI,JKLP,:) )
ENDDO

!*      2.2     Final updraft loop
!               ------------------

DO JI = 1, IIE
   !DO JK = MINVAL( KDPL(:) ), KCTL(JI)
   DO JK = KDPL(JI), KCTL(JI)
      JKP = JK + 1

      DO JN = 1, INCH1
       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK ) ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)

       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
                     !if you have reactive i.e. non-passive tracers
                     ! add the corresponding sink term in the following equation
           ZUCH1(JI,JKP,JN) = ( PUMF(JI,JK) * ZUCH1(JI,JK,JN) +              &
                            &   PUER(JI,JKP) * PCH1(JI,JK,JN) )  /           &
                            & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 )
       ENDIF
      ENDDO
   ENDDO
ENDDO

!*      3.      Downdraft computations
!               ----------------------

ZDCH1(:,:,:) = _ZERO_

!*      3.1     Initialization at the LFS
!               -------------------------

ZWORK1(:,:) = SPREAD( PMIXF(:), DIM=2, NCOPIES=INCH1 )
DO JI = 1, IIE
     JK = KLFS(JI)
     ZDCH1(JI,JK,:) = ZWORK1(JI,:) * PCH1(JI,JK,:) +                          &
                    &                  ( _ONE_ - ZWORK1(JI,:) ) * ZUCH1(JI,JK,:)
ENDDO

!*      3.2     Final downdraft loop
!               --------------------

!DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
DO JI = 1, IIE
    DO JK = KLFS(JI), IKB + 1, -1
       JKP = JK - 1
    DO JN = 1, INCH1
      IF ( JK <= KLFS(JI) .AND. JKP >= KDBL(JI) ) THEN
       ZDCH1(JI,JKP,JN) = ( ZDCH1(JI,JK,JN) * PDMF(JI,JK) -              &
                        &   PCH1(JI,JK,JN) *  PDER(JI,JKP) ) /           &
                        & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 )
      ENDIF
    ENDDO
    ENDDO
ENDDO


!*      4.      Final closure (environmental) computations
!               ------------------------------------------

PCH1C(:,IKB:IKE,:) = PCH1(:,IKB:IKE,:) ! initialize adjusted envir. values

DO JK = IKB, IKE
   ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG ! environmental subsidence
ENDDO

!PLEASE CHECK HERE!!!
!Recalculate the fractional time step like in convect_closure_shal -vlee
!ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                           ! to be an integer multiple of PTIMEC
ZTIMEC(:) = PTIMEC(:)
ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1
ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) )
WHERE ( PTIMEC(:) < _ONE_ ) ZTIMEC(:) = _ZERO_
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )

!These variables are inside the loop in the convect_uv_transport_shal
!perhaps due to the fractional time step loop that I have implemented
!in here with itstep(:) instead of kftsteps -vlee
!PLEASE CHECK HERE
!ZCH1MFIN(:,:,:)   = _ZERO_
!ZCH1MFOUT(:,:,:)  = _ZERO_


DO JI = 1, IIE
   FRACTIONAL_STEPS:DO JSTEP = 1, ITSTEP(JI)! Enter the fractional time step loop
!DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
      ZCH1MFIN(:,:,:)   = _ZERO_
      ZCH1MFOUT(:,:,:)  = _ZERO_

      DO JK = IKB + 1, KCTL(JI)
         JKP = MAX( IKB + 1, JK - 1 )
          ZWORK3(JI,1) = ZOMG(JI,JK)
          ZWORK1(JI,1) = SIGN( _ONE_, ZWORK3(JI,1) )
          ZWORK2(JI,1) = _HALF_ * ( 1. + ZWORK1(JI,1) )
          ZWORK1(JI,1) = _HALF_ * ( 1. - ZWORK1(JI,1) )
          ZCH1MFIN(JI,JK,:)  = - ZWORK3(JI,1) * PCH1C(JI,JKP,:) * ZWORK1(JI,1)
          ZCH1MFOUT(JI,JK,:) =   ZWORK3(JI,1) * PCH1C(JI,JK,:)  * ZWORK2(JI,1)
          ZCH1MFIN(JI,JKP,:) = ZCH1MFIN(JI,JKP,:) + ZCH1MFOUT(JI,JK,:) * ZWORK2(JI,1)
          ZCH1MFOUT(JI,JKP,:)= ZCH1MFOUT(JI,JKP,:) + ZCH1MFIN(JI,JK,:) * ZWORK1(JI,1)
      END DO
!
      DO JK = IKB + 1, KCTL(JI)
      DO JN = 1, INCH1
         PCH1C(JI,JK,JN) = PCH1C(JI,JK,JN) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (    &
                      ZCH1MFIN(JI,JK,JN) + PUDR(JI,JK) * ZUCH1(JI,JK,JN) +        &
                      PDDR(JI,JK) * ZDCH1(JI,JK,JN) - ZCH1MFOUT(JI,JK,JN) -       &
                      ( PUER(JI,JK) + PDER(JI,JK) ) * PCH1(JI,JK,JN)    )
      !  PCH1C(JI,JK,JN) = MAX( 0., PCH1C(JI,JK,JN) )
      END DO
      END DO

   ENDDO FRACTIONAL_STEPS! final values
ENDDO

END SUBROUTINE CONVECT_CHEM_TRANSPORT

