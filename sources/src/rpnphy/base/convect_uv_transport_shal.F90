!#######################################################################
  SUBROUTINE CONVECT_UV_TRANSPORT_SHAL( KLON, KLEV, PU, PV, PUC, PVC,       & 
                                 & KDPL, KPBL, KLCL, KCTL,      &  
                                 & PUMF, PUER, PUDR,            &
                                 & PTIMEC, PDXDY, PLMASS, PWSUB,     &
                                 & KFTSTEPS, PTIMC, ITSTEP,GWORK1 ) 
!#######################################################################

!!**** Compute  modified horizontal wind components due to convective event
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine convective adjusted
!!      horizontal wind components u and v
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PUC-PU)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      tracer routine but includes pressure term
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
!!      Original    11/02/02
!!
!!Jing, May 2017, delete all parts realted to downdraft to clean the code  
!!
!-------------------------------------------------------------------------------


!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAR
USE YOE_CONVPAREXT

IMPLICIT NONE
#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0

!*       0.1   Declarations of dummy arguments :

integer,                INTENT(IN) :: KLON     ! horizontal dimension
integer,                INTENT(IN) :: KLEV     ! vertical dimension

real,DIMENSION(KLON,KLEV),INTENT(IN) :: PU     ! horizontal wind in x (m/s)
real,DIMENSION(KLON,KLEV),INTENT(IN) :: PV     ! horizontal wind in x (m/s)
real,DIMENSION(KLON,KLEV),INTENT(OUT):: PUC    ! convective adjusted value of u (m/s)
real,DIMENSION(KLON,KLEV),INTENT(OUT):: PVC    ! convective adjusted value of v (m/s)

integer, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
integer, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
integer, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
integer, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level

real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)

real, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step ,modify to fractional timestep
real, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
integer,                   INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
real, DIMENSION(KLON),     INTENT(IN) :: PTIMC   ! fractional convective time step
integer,DIMENSION(KLON),   INTENT(IN) :: ITSTEP  ! # of fractional convective timesteps
LOGICAL, DIMENSION(KLON),  INTENT(IN)   :: GWORK1  ! flag for newly activated columns including MASS CONSERVATION criteria from closure

!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
integer :: IKS            ! vertical dimension
integer :: JI             ! horizontal loop index
integer :: JK, JKP        ! vertical loop index
integer :: JSTEP          ! fractional time loop index
integer :: JKLD, JKLP, JKMAX ! loop index for levels

integer, PARAMETER             :: IUV = 2    ! for u and v
real, DIMENSION(KLON,KLEV)     :: ZOMG       ! compensat. subsidence (Pa/s)
real, DIMENSION(KLON,KLEV,IUV) :: ZUUV        ! updraft values
real, DIMENSION(KLON)          :: ZTIMEC     ! fractional convective time step
real, DIMENSION(KLON,KLEV,IUV) :: ZUVMFIN, ZUVMFOUT
                                   ! work arrays for environm. compensat. mass
real, DIMENSION(KLON,IUV)      :: ZWORK1, ZWORK2, ZWORK3

integer, DIMENSION(KLON)  :: ICOUNT    ! timestep counter
LOGICAL, DIMENSION(KLON)  :: GWORK3      

!-------------------------------------------------------------------------------

!*       0.3   Compute loop bounds
!              -------------------

IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
!JKMAX  = MAXVAL( KCTL(:) )


!*      2.      Updraft computations
!               --------------------

ZUUV(:,:,:) = _ZERO_

!*      2.1     Initialization  at LCL
!               ----------------------------------

DO JI = 1, IIE
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,1) = _HALF_ * ( PU(JI,JKLD) + PU(JI,JKLP) )
    ZWORK1(JI,2) = _HALF_ * ( PV(JI,JKLD) + PV(JI,JKLP) )
ENDDO

!*      2.2     Final updraft loop
!               ------------------

DO JI = 1, IIE
   !DO JK = MINVAL(KDPL(:)), JKMAX
   DO JK = KDPL(JI),KCTL(JI)
      JKP = JK + 1

       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK ) THEN
            ZUUV(JI,JK,1) = ZWORK1(JI,1)
            ZUUV(JI,JK,2) = ZWORK1(JI,2)
       END IF

       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
                            ! instead of passive tracers equations
                            ! wind equations also include pressure term
           ZUUV(JI,JKP,1) = ( PUMF(JI,JK) * ZUUV(JI,JK,1) +                   &
                            &   PUER(JI,JKP) * PU(JI,JK) )  /                 &
                            & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 ) +  &
                            &   XUVDP * ( PU(JI,JKP) - PU(JI,JK) ) 
           ZUUV(JI,JKP,2) = ( PUMF(JI,JK) * ZUUV(JI,JK,2) +                   &
                            &   PUER(JI,JKP) * PV(JI,JK) )  /                 &
                            & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 ) +  &
                            &   XUVDP * ( PV(JI,JKP) - PV(JI,JK) ) 
       ENDIF
     ENDDO
ENDDO


!*      4.      Final closure (environmental) computations
!               ------------------------------------------

PUC(:,IKB:IKE) = PU(:,IKB:IKE) ! initialize adjusted envir. values
PVC(:,IKB:IKE) = PV(:,IKB:IKE) ! initialize adjusted envir. values

 DO JK = IKB, IKE
    ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG ! environmental subsidence
 ENDDO

ZTIMEC(:) = PTIMC(:)  
                                         ! to be an integer multiple of PTIMEC
WHERE ( PTIMC(:) < _ONE_ ) ZTIMEC(:) = _ZERO_
!ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )

 ILOOP: DO JI = 1, IIE
     IF ( GWORK1(JI) )  THEN
       FRACTIONAL_SUBSTEPS: DO JSTEP = 1, ITSTEP(JI)
       ! Enter the fractional time step loop for each column

        ZUVMFIN(:,:,:)   = _ZERO_         
        ZUVMFOUT(:,:,:)  = _ZERO_        

        !Do this K loop first to sum up
        !DO JK = IKB + 1, JKMAX
         DO JK = IKB + 1, KCTL(JI)
          JKP = MAX( IKB + 1, JK - 1 )
          ZWORK3(JI,1) = ZOMG(JI,JK)
          ZWORK1(JI,1) = SIGN( _ONE_, ZWORK3(JI,1) )
          ZWORK2(JI,1) = _HALF_ * ( _ONE_ + ZWORK1(JI,1) )
          ZWORK1(JI,1) = _HALF_ * ( _ONE_ - ZWORK1(JI,1) )
          ZUVMFIN(JI,JK,1)  = - ZWORK3(JI,1) * PUC(JI,JKP) * ZWORK1(JI,1)
          ZUVMFOUT(JI,JK,1) =   ZWORK3(JI,1) * PUC(JI,JK)  * ZWORK2(JI,1)
          ZUVMFIN(JI,JK,2)  = - ZWORK3(JI,1) * PVC(JI,JKP) * ZWORK1(JI,1)
          ZUVMFOUT(JI,JK,2) =   ZWORK3(JI,1) * PVC(JI,JK)  * ZWORK2(JI,1)
          ZUVMFIN(JI,JKP,1) = ZUVMFIN(JI,JKP,1) + ZUVMFOUT(JI,JK,1) * ZWORK2(JI,1)
          ZUVMFIN(JI,JKP,2) = ZUVMFIN(JI,JKP,2) + ZUVMFOUT(JI,JK,2) * ZWORK2(JI,1)
          ZUVMFOUT(JI,JKP,1)= ZUVMFOUT(JI,JKP,1)+ ZUVMFIN(JI,JK,1)  * ZWORK1(JI,1)
          ZUVMFOUT(JI,JKP,2)= ZUVMFOUT(JI,JKP,2)+ ZUVMFIN(JI,JK,2)  * ZWORK1(JI,1)
         END DO

        !Then add on for the tendencies
        !DO JK = IKB + 1, JKMAX
         DO JK = IKB + 1, KCTL(JI)
         PUC(JI,JK) = PUC(JI,JK) + ZTIMEC(JI) / PLMASS(JI,JK) *  (       &
                   &   ZUVMFIN(JI,JK,1) + PUDR(JI,JK) * ZUUV(JI,JK,1)     &
                   &   - ZUVMFOUT(JI,JK,1) -   &
                   &   PUER(JI,JK) * PU(JI,JK)    )
         PVC(JI,JK) = PVC(JI,JK) + ZTIMEC(JI) / PLMASS(JI,JK) *  (       &
                   &   ZUVMFIN(JI,JK,2) + PUDR(JI,JK) * ZUUV(JI,JK,2)     &
                   &    - ZUVMFOUT(JI,JK,2) -   &
                   &    PUER(JI,JK) * PV(JI,JK)    )
         ENDDO
       ENDDO FRACTIONAL_SUBSTEPS !exit the fractional time step loop
     END IF
  ENDDO ILOOP

END SUBROUTINE CONVECT_UV_TRANSPORT_SHAL

