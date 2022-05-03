!#######################################################################
 SUBROUTINE CONVECT_DOWNDRAFT2( KLON, KLEV,                             &
                          & KICE, PPRES, PDPRES, PZ, PTH, PTHES,       &
                          & PRW, PRC, PRI,                             &
                          & PPREF, KLCL, KCTL, KETL,                   &
                          & PUTHL, PURW, PURC, PURI,                   &
                          & PDMF, PDER, PDDR, PDTHL, PDRW,             &
                  !Jingorg& PMIXF, PDTEVR, KLFS, KDBL, KML,            &
                          & PMIXF, PDTEVR, KLFS, KDBL, KML,IDDT,       &
                          & PDTEVRF )
!#######################################################################

!!**** Compute downdraft properties from LFS to DBL.
!!
!!
!!    PDRPOSE
!!    -------
!!      The purpose of this routine is to determine downdraft properties
!!      ( mass flux, thermodynamics )
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from top.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
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
!!
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RPI                ! Pi
!!          RATM               ! reference pressure
!!          RD, RV           ! gaz  constants for dry air and water vapor
!!          RCPD               ! Cpd (dry air)
!!          RCPV, RCW, RCS     ! Cp of water vapor, liquid water and ice
!!          RTT                ! triple point temperature
!!          RLVTT, RLSTT       ! vaporisation/sublimation heat at RTT
!!
!!      Module YOE_CONVPAR
!!          XCRAD              ! cloud radius
!!          XZPBL              ! thickness of downdraft detrainment layer
!!          XENTR              ! entrainment constant in pressure coordinates
!!          XRHDBC             ! relative humidity in downdraft below cloud
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_DOWNDRAFT)
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

USE YOMCST
USE YOE_CONVPAR
USE YOE_CONVPAREXT

IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :


integer,                    INTENT(IN) :: KLON  ! horizontal dimension
integer,                    INTENT(IN) :: KLEV  ! vertical dimension
integer,                    INTENT(IN) :: KICE  ! flag for ice ( 1 = yes,
                                                  !                0 = no ice )
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! grid scale theta
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water
                                                  ! mixing ratio
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRC   ! grid scale r_c (cloud water)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRI   ! grid scale r_i (cloud ice)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between
                                                  ! bottom and top of layer (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! level height (m)
integer, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
integer, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
integer, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of
                                                  ! equilibrium (zero buoyancy) level
integer, DIMENSION(KLON),   INTENT(IN) :: KML   ! " vert. index of melting level
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUTHL ! updraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PURC  ! updraft r_c (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PURI  ! updraft r_i (kg/kg)
real, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency


real, DIMENSION(KLON,KLEV), INTENT(OUT):: PDMF   ! downdraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PDER   ! downdraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PDDR   ! downdraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTHL  ! downdraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PDRW   ! downdraft total water (kg/kg)
real, DIMENSION(KLON),      INTENT(OUT):: PMIXF  ! mixed fraction at LFS
real, DIMENSION(KLON),      INTENT(OUT):: PDTEVR ! total downdraft evaporation
                                                   ! rate at LFS (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTEVRF! downdraft evaporation rate
integer, DIMENSION(KLON),  INTENT(OUT):: KLFS    ! contains vert. index of LFS
integer, DIMENSION(KLON),  INTENT(OUT):: KDBL    ! contains vert. index of DBL
integer, DIMENSION(KLON),  INTENT(OUT):: IDDT    ! contains vert. index of DDT, !Jing

!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE     ! horizontal + vertical loop bounds
integer :: JK, JKP, JKM, JKT ! vertical loop index
integer :: JI, JL            ! horizontal loop index
integer :: JITER             ! iteration loop index
real    :: ZCPORD, ZRDOCP    ! C_pd / R_d, R_d / C_pd
real    :: ZEPS              ! R_d / R_v
real    :: ZEPSA, ZCVOCD     ! R_v / R_d, C_pv / C_pd

!Jingorg integer, DIMENSION(KLON) :: IDDT      ! top level of detrainm. layer
real, DIMENSION(KLON)    :: ZTHE      ! environm. theta_e (K)
real, DIMENSION(KLON)    :: ZDT, ZDTP ! downdraft temperature (K)
real, DIMENSION(KLON)    :: ZCPH      ! specific heat C_ph
real, DIMENSION(KLON)    :: ZLV, ZLS  ! latent heat of vaporis., sublim.
real, DIMENSION(KLON)    :: ZDDT      ! thickness (hPa) of detrainm. layer
real, DIMENSION(KLON)    :: ZPI       ! Pi=(P0/P)**(Rd/Cpd)
real, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
real, DIMENSION(KLON)    :: ZWORK4, ZWORK5          ! work arrays
LOGICAL, DIMENSION(KLON)   :: GWORK1                  ! work array


!-------------------------------------------------------------------------------

!        0.3    Set loop bounds
!               ---------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Initialize downdraft properties
!               -------------------------------

ZCPORD     = RCPD / RD
ZRDOCP     = RD  / RCPD
ZEPS       = RD  / RV
ZEPSA      = RV  / RD
ZCVOCD     = RCPV / RCPD

PDMF(:,:)  = _ZERO_
PDER(:,:)  = _ZERO_
PDDR(:,:)  = _ZERO_
PDRW(:,:)  = _ZERO_
PDTHL(:,:) = _ZERO_
PDTEVR(:)  = _ZERO_
PMIXF(:)   = _ZERO_
ZTHE(:)    = _ZERO_
ZDDT(:)    = PDPRES(:,IKB+2)
KDBL(:)    = IKB + 1
KLFS(:)    = IKB + 1
IDDT(:)    = KDBL(:) + 1


!*       2.     Determine the LFS by looking for minimum of environmental
!               saturated theta_e
!               ----------------------------------------------------------

ZWORK1(:) = 900.   ! starting value for search of minimum envir. theta_e
DO JK = MINVAL( KLCL(:) ) + 2, MAXVAL( KETL(:) )
   DO JI = 1, IIE
      GWORK1(JI) = JK >= KLCL(JI) + 2 .AND. JK < KETL(JI)
      IF ( GWORK1(JI) .AND. ZWORK1(JI) > PTHES(JI,JK) ) THEN
         KLFS(JI)   = JK
         ZWORK1(JI) = MIN( ZWORK1(JI), PTHES(JI,JK) )
      ENDIF
   ENDDO
ENDDO


!*       3.     Determine the mixed fraction using environmental and updraft
!               values of theta_e at LFS
!               ---------------------------------------------------------

DO JI = 1, IIE
    JK = KLFS(JI)
    ZPI(JI)    = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
      ! compute updraft theta_e
    ZWORK3(JI) = PURW(JI,JK) - PURC(JI,JK) - PURI(JI,JK)
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI)
    ZLV(JI)    = RLVTT + ( RCPV - RCW ) * ( ZDT(JI) - RTT )
    ZLS(JI)    = RLSTT + ( RCPV - RCS ) * ( ZDT(JI) - RTT )
    ZCPH(JI)   = RCPD + RCPV * PURW(JI,JK)
    ZDT(JI)    = ( PUTHL(JI,JK) - ( _ONE_ + PURW(JI,JK) ) * RG * PZ(JI,JK) &
               & + ZLV(JI) * PURC(JI,JK) + ZLS(JI) * PURI(JI,JK) ) / ZCPH(JI)
    ZWORK1(JI) = ZDT(JI) * ZPI(JI) ** ( _ONE_ - 0.28 * ZWORK3(JI) )   &
               & * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )         &
               & * ZWORK3(JI) * ( _ONE_ + 0.81 * ZWORK3(JI) ) )
      ! compute environmental theta_e
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI)
    ZLV(JI)    = RLVTT + ( RCPV - RCW ) * ( ZDT(JI) - RTT )
    ZLS(JI)    = RLSTT + ( RCPV - RCS ) * ( ZDT(JI) - RTT )
    ZWORK3(JI) = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
    ZCPH(JI)   = RCPD + RCPV * PRW(JI,JK)
    ZWORK2(JI) = ZDT(JI) * ZPI(JI) ** ( _ONE_ - 0.28 * ZWORK3(JI) ) &
               & * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )       &
               & * ZWORK3(JI) * ( _ONE_ + 0.81 * ZWORK3(JI) ) )
      ! compute mixed fraction
    PMIXF(JI)  = MAX( _ZERO_, ( ZWORK1(JI) - PTHES(JI,JK) ) )            &
               &  / ( ZWORK1(JI) - ZWORK2(JI) + 1.E-10 )
    PMIXF(JI)  = MAX(_ZERO_, MIN( _ONE_, PMIXF(JI) ) )
    ZWORK4(JI) = PPRES(JI,JK)
ENDDO


!*       4.     Estimate the effect of melting on the downdraft
!               ---------------------------------------------

ZWORK1(:) = _ZERO_
      ! use total solid precipitation
!DO JK = IKB + 1, IKE
!    ZWORK1(:) = ZWORK1(:) + PURS(:,JK) ! total snow/hail content
!END DO

DO JI = 1, IIE
     JK  = KLCL(JI)
     JKP = KCTL(JI)
     ZWORK1(JI) = _HALF_ * ( PURW(JI,JK) - PURW(JI,JKP) )
ENDDO

      ! temperature perturbation due to melting at LFS
ZWORK3(:) = _ZERO_
WHERE( KML(:) > IKB + 2 )
          ZWORK3(:) = ZWORK1(:) * ( ZLS(:) - ZLV(:) ) / ZCPH(:)
          ZDT(:)    = ZDT(:) - ZWORK3(:) * REAL(KICE)
END WHERE


!*       5.     Initialize humidity at LFS as a saturated mixture of
!               updraft and environmental air
!               -----------------------------------------------------

DO JI = 1, IIE
     JK = KLFS(JI)
     PDRW(JI,JK)  = PMIXF(JI) * PRW(JI,JK) + ( _ONE_ - PMIXF(JI) ) * PURW(JI,JK)
     ZWORK2(JI)   = PDRW(JI,JK) - ( _ONE_ - PMIXF(JI) )                          &
                  &                  * ( PURC(JI,JK) + PURI(JI,JK) )
ENDDO


!*       6.1    Determine the DBL by looking for level where the envir.
!               theta_es at the LFS corrected by melting effects  becomes
!               larger than envir. value
!               ---------------------------------------------------------

! compute satur. mixing ratio for melting corrected temperature
    CALL CONVECT_SATMIXRATIO( KLON, ZWORK4, ZDT, ZWORK3, ZLV, ZLS, ZCPH )

      ! compute envir. saturated theta_e for melting corrected temperature
    ZWORK1(:) = MIN( ZWORK2(:), ZWORK3(:) )
    ZWORK3(:) = ZWORK3(:) * ZWORK4(:) / ( ZWORK3(:) + ZEPS ) ! sat. pressure
    ZWORK3(:) = LOG( ZWORK3(:) / 613.3 )

      ! dewp point temperature
    ZWORK3(:) = ( 4780.8 - 32.19 * ZWORK3(:) ) / ( 17.502 - ZWORK3(:) )

      ! adiabatic saturation temperature
    ZWORK3(:) = ZWORK3(:) - ( .212 + 1.571E-3 * ( ZWORK3(:) - RTT )          &
              &   - 4.36E-4 * ( ZDT(:) - RTT ) ) * ( ZDT(:) - ZWORK3(:) )
    ZWORK4(:) = SIGN(_HALF_, ZWORK2(:) - ZWORK3(:) )
    ZDT(:)    = ZDT(:) * ( _HALF_ + ZWORK4(:) ) + ( _HALF_ - ZWORK4(:) ) * ZWORK3(:)
    ZWORK2(:) = ZDT(:) * ZPI(:) ** ( _ONE_ - 0.28 * ZWORK2(:) )                   &
              &                   * EXP( ( 3374.6525 / ZDT(:) - 2.5403 )     &
              &                   * ZWORK1(:) * ( _ONE_ + 0.81 * ZWORK1(:) ) )

GWORK1(:) = .TRUE.
JKM = MAXVAL( KLFS(:) )
DO JK = JKM - 1, IKB + 1, -1
  DO JI = 1, IIE
     IF ( JK < KLFS(JI) .AND. ZWORK2(JI) > PTHES(JI,JK) .AND. GWORK1(JI) ) THEN
          KDBL(JI) = JK
          GWORK1(JI) = .FALSE.
     ENDIF
  ENDDO
ENDDO


!*       7.     Define mass flux and entr/detr. rates at LFS
!               -------------------------------------------

DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI)  = PPRES(JI,JK) /                                            &
                 & ( RD * ZDT(JI) * ( _ONE_ + ZEPS * ZWORK1(JI) ) ) ! density
     PDMF(JI,JK) = - ( _ONE_ - PPREF(JI) ) * ZWORK1(JI) * RPI * XCRAD * XCRAD
     PDTHL(JI,JK)= ZWORK2(JI)   ! theta_l is here actually theta_e
     ZWORK2(JI)  = PDMF(JI,JK)
     PDDR(JI,JK) = _ZERO_
     PDER(JI,JK) = - PMIXF(JI) * PDMF(JI,JK)
ENDDO


!         7.1   Downdraft detrainment is assumed to occur in a layer
!               of 60 hPa, determine top level IDDT of this layer
!               ---------------------------------------------------------

ZWORK1(:) = _ZERO_
DO JK = IKB + 2, JKM
      ZWORK1(:) = ZWORK1(:) + PDPRES(:,JK)
      WHERE ( JK > KDBL(:) .AND. ZWORK1(:) <= XZPBL )
           ZDDT(:) = ZWORK1(:)
           IDDT(:) = JK
      END WHERE
ENDDO


!*       8.     Enter loop for downdraft computations. Make a first guess
!               of initial downdraft mass flux.
!               In the downdraft computations we use theta_es instead of
!               enthalpy as it allows to better take into account evaporation
!               effects. As the downdraft detrainment rate is zero apart
!               from the detrainment layer, we just compute enthalpy
!               downdraft from theta_es in this layer.
!               ----------------------------------------------------------


ZWORK5(:) = _ZERO_

DO JK =  JKM - 1, IKB + 1, -1
  JKP = JK + 1
  DO JI = 1, IIE
    IF ( JK < KLFS(JI) .AND. JK >= IDDT(JI) )  THEN
      PDER(JI,JK)  = - ZWORK2(JI) * XENTR * PDPRES(JI,JKP) / XCRAD
                                               ! DER and DPRES are positive
      PDMF(JI,JK)  = PDMF(JI,JKP) - PDER(JI,JK)
      ZPI(JI)      = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
      ZDT(JI)      = PTH(JI,JK) / ZPI(JI)
      ZWORK1(JI)   = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
      ZTHE(JI)     = ZDT(JI) * ZPI(JI) ** ( _ONE_ - 0.28 * ZWORK1(JI) )       &
                   &           * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )   &
                   &           * ZWORK1(JI) * ( _ONE_ + 0.81 * ZWORK1(JI) ) )
         ! PDTHL is here theta_es, later on in this routine this table is
         ! reskipped to enthalpy
      PDTHL(JI,JK) = ( PDTHL(JI,JKP) * PDMF(JI,JKP) - ZTHE(JI) * PDER(JI,JK) ) /   &
                   & ( PDMF(JI,JK) - 1.E-7 )
      PDRW(JI,JK)  = ( PDRW(JI,JKP) * PDMF(JI,JKP) - PRW(JI,JK) * PDER(JI,JK) ) /  &
                   & ( PDMF(JI,JK) - 1.E-7 )
    ENDIF
    IF ( JK < IDDT(JI) .AND. JK >= KDBL(JI) )   THEN
      JL = IDDT(JI)
      PDDR(JI,JK)  = - PDMF(JI,JL) * PDPRES(JI,JKP) / ZDDT(JI)
      PDMF(JI,JK)  = PDMF(JI,JKP) + PDDR(JI,JK)
      PDTHL(JI,JK) = PDTHL(JI,JKP)
      PDRW(JI,JK)  = PDRW(JI,JKP)
    ENDIF
  ENDDO
ENDDO


!*       9.     Calculate total downdraft evaporation
!               rate for given mass flux (between DBL and IDDT)
!               -----------------------------------------------

PDTEVRF(:,:) = _ZERO_

JKT = MAXVAL( IDDT(:) )
DO JK = IKB + 1, JKT

       ZPI(:) = ( RATM / PPRES(:,JK) ) ** ZRDOCP
       ZDT(:) = PTH(:,JK) / ZPI(:)

!*       9.1    Determine wet bulb temperature at DBL from theta_e.
!               The iteration algoritm is similar to that used in
!               routine CONVECT_CONDENS
!               --------------------------------------------------

   DO JITER = 1, 4
       CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK1, ZLV, ZLS, ZCPH )
       ZDTP(:) = PDTHL(:,JK) / ( ZPI(:) ** ( _ONE_ - 0.28 * ZWORK1(:) )      &
               &      * EXP( ( 3374.6525 / ZDT(:) - 2.5403 )            &
               &             * ZWORK1(:) * ( _ONE_ + 0.81 * ZWORK1(:) ) ) )
       ZDT(:)  = 0.4 * ZDTP(:) + 0.6 * ZDT(:) ! force convergence
   ENDDO


!*       9.2    Sum total downdraft evaporation rate. No evaporation
!               if actual humidity is larger than specified one.
!               -----------------------------------------------------

   ZWORK2(:) = ZWORK1(:) / ZDT(:) * ( RBETW / ZDT(:) - RGAMW ) ! dr_sat/dT
   ZWORK2(:) = ZLV(:) / ZCPH(:) * ZWORK1(:) * ( _ONE_ - XRHDBC ) /              &
             &      ( _ONE_ + ZLV(:) / ZCPH(:) * ZWORK2(:) ) ! temperature perturb                                                           ! due to evaporation
   ZDT(:)    = ZDT(:) + ZWORK2(:)

   CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK3, ZLV, ZLS, ZCPH )

   ZWORK3(:)    = ZWORK3(:) * XRHDBC
   ZWORK1(:)    = MAX( _ZERO_, ZWORK3(:) - PDRW(:,JK) )
   PDTEVR(:)    = PDTEVR(:) + ZWORK1(:) * PDDR(:,JK)
   PDTEVRF(:,JK)= PDTEVRF(:,JK) + ZWORK1(:) * PDDR(:,JK)
      ! compute enthalpie and humidity in the detrainment layer
   PDRW(:,JK)   = MAX( PDRW(:,JK), ZWORK3(:) )
   PDTHL(:,JK)  = ( ( RCPD + PDRW(:,JK) * RCPV ) * ZDT(:)                    &
                &   + ( _ONE_ + PDRW(:,JK) ) * RG * PZ(:,JK) )

ENDDO


!*      12.     If downdraft does not evaporate any water for specified
!               relative humidity, no downdraft is allowed
!               ---------------------------------------------------------

ZWORK2(:) = _ONE_
WHERE ( PDTEVR(:) < _ONE_ .OR. KLFS(:) == IKB + 1 ) ZWORK2(:) = _ZERO_
DO JK = IKB, JKM
      KDBL(:)     = KDBL(:) * INT( ZWORK2(:) ) + ( 1 - INT( ZWORK2(:) ) ) * IKB
      KLFS(:)     = KLFS(:) * INT( ZWORK2(:) ) + ( 1 - INT( ZWORK2(:) ) ) * IKB
      PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
      PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:)
      PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:)
      ZWORK1(:)   = REAL( KLFS(:) - JK )              ! use this to reset thl_d
      ZWORK1(:)   = MAX( _ZERO_,MIN(_ONE_,ZWORK1(:) ) ) ! and rv_d to zero above LFS
      PDTHL(:,JK) = PDTHL(:,JK) * ZWORK2(:) * ZWORK1(:)
      PDRW(:,JK)  = PDRW(:,JK)  * ZWORK2(:) * ZWORK1(:)
      PDTEVR(:)   = PDTEVR(:)   * ZWORK2(:)
      PDTEVRF(:,JK)= PDTEVRF(:,JK) * ZWORK2(:)
ENDDO

END SUBROUTINE CONVECT_DOWNDRAFT2

