!-------------------------------------------------------------------------------
!   ############################################################################
    SUBROUTINE BKF_DEEP3 ( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,             &
                         & PDTCONV,  ODOWN, KICE,                              &
                         & KENSM,                                              &
                         & PPABS, PZZ, PDXDY, PHSFLX,                          &
                         & PT, PRV, PRC, PRI, PU, PV, PW,                      &
                         & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,              &
                         & ZPRLTEN, PPRSTEN,                                   &
                         & PUMF, PDMF, PPRLFLX, PPRSFLX, PCAPE, ZTOP, ZBAS,    &
                         & PURV, PURC, PURI, PCFD,                             &
                         & PAREAUP,  PWUMAXOUT,                                &
                         & OUVTRANS, PUTEN, PVTEN,                             &
                         & OCHTRANS, KCH1, PCH1, PCH1TEN,                      &
                         & PUDR, PDDR, PKKFC, PSP, XLAT, MG, MLAC              ) ! for ERA40
!   ############################################################################
       use phy_status, only: phy_error_L
!
!!**** Interface routine to the fast Meso-NH convection code developed for ECMWF/ARPEGE
!!     having a structure typical for operational routines
!!
!!     Transformations necessary to call deep+ shallow code
!!     - skip input vertical arrays/levels : bottom=1, top=KLEV
!!     - transform specific humidities in mixing ratio
!!
!!
!!    PURPOSE
!!    -------
!!      The routine interfaces the MNH convection code as developed for operational
!!      forecast models like ECMWF/ARPEGE or HIRLAM with the typical Meso-NH array structure
!!      Calls the deep and/or shallow convection routine
!!
!!
!!**  METHOD
!!    ------
!!     Returns one tendency for shallow+deep convection but each part can
!!     be activated/desactivated separately
!!     For deep convection one can enable up to 3 additional ensemble members
!!     - this substantially improves the smoothness of the scheme and reduces
!!       allows for runs with different cloud radii (entrainment rates) and
!!       reduces the arbitrariness inherent to convective trigger condition
!!
!!
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_DEEP
!!    SU_CONVPAR, SU_CONVPAR1
!!    SUCST:   ECMWF/ARPEGE routine
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/12/98
!!      modified    20/03/2002 by P. Marquet : transformed for ARPEGE/Climat
!!                             (tsmbkind.h, REAL_B, INTEGER_M, _JPRB, "& &",
!!                              _ZERO_, _ONE_, _HALF_)
!!      modified    11/04/O2 allow for ensemble of deep updrafts/downdrafts
!!
!!    REFERENCE
!!    ---------
!!    Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc., Vol 127, pp 869-886:
!!           A mass flux convection scheme for regional and global models.
!!
!-------------------------------------------------------------------------------
!Revisions in RPN
! 001   -PV-jan2015 - do not aggregate deep and shallow tendencies; add output of shallow tendencies
! 002   -PV-apr2015 - output cloud top and base heights for deep only
! 003   -JY-Jul2015 - separate from convection.ftn90 (which has both deep and shallow convection scheme)
!                   - add output of UA(area of updraft), K6(maximum upward velocity in Bechtold scheme)
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :
!
!
integer,                    INTENT(IN)   :: KLON   ! horizontal dimension
integer,                    INTENT(IN)   :: KLEV   ! vertical dimension
integer,                    INTENT(IN)   :: KIDIA  ! value of the first point in x
integer,                    INTENT(IN)   :: KFDIA  ! value of the last point in x
integer,                    INTENT(IN)   :: KBDIA  ! vertical  computations start at
!                                                    ! KBDIA that is at least 1
integer,                    INTENT(IN)   :: KTDIA  ! vertical computations can be
                                                     ! limited to KLEV + 1 - KTDIA
                                                     ! default=1
real,                       INTENT(IN)   :: PDTCONV! Interval of time between two
                                                     ! calls of the deep convection
                                                     ! scheme
LOGICAL,                      INTENT(IN)   :: ODOWN  ! take or not convective
                                                     ! downdrafts into account
integer,                    INTENT(IN)   :: KICE   ! flag for ice ( 1 = yes,
                                                     !                0 = no ice )
integer,                    INTENT(IN)   :: KENSM  ! number of additional deep convection calls
                                                     ! for ensemble (presently limited to 3)
                                                     ! KENSM=0 corresponds to base run with
                                                     ! 1 deep and 1 shallow call
!
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PT     ! grid scale T at time t  (K)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PRV    ! grid scale water vapor  (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PRC    ! grid scale r_c (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PRI    ! grid scale r_i (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PU     ! grid scale horiz. wind u (m/s)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PV     ! grid scale horiz. wind v (m/s)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PW     ! grid scale vertical velocity (m/s)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PPABS  ! grid scale pressure (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PZZ    ! geopotential (m2/s2)
real, DIMENSION(KLON),      INTENT(IN)   :: PDXDY  ! grid area (m2)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PHSFLX ! turbulent sensible heat flux (W/m^2)
integer, DIMENSION(KLON), INTENT(INOUT)  :: KCOUNT ! convective counter(recompute
                                                     ! tendency or keep it
!
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperat. tendency (K/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
!real, DIMENSION(KLON),      INTENT(INOUT):: PPRTEN ! total surf precipitation tendency (m/s)
real, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf precipitation tendency (kg/s m2)
real, DIMENSION(KLON),      INTENT(INOUT):: ZPRLTEN! liquid surf precipitation tendency (kg/s m2)
!
! convective U,V transport (Cu friction)
LOGICAL,                      INTENT(IN)        :: OUVTRANS ! flag to compute convective
                                                            ! transport for chemical tracer
real, DIMENSION(KLON,KLEV), INTENT(INOUT)     :: PUTEN    ! convecctive u tendency (m/s^2)
real, DIMENSION(KLON,KLEV), INTENT(INOUT)     :: PVTEN    ! convecctive v tendency (m/s^2)

! convective Chemical Tracers:
LOGICAL,                      INTENT(IN)        :: OCHTRANS ! flag to compute convective
                                                            ! transport for chemical tracer
integer,                    INTENT(IN)        :: KCH1     ! number of species
real, DIMENSION(KLON,KLEV,KCH1), INTENT(IN)   :: PCH1     ! grid scale chemical species
real, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN  ! chemical convective tendency
                                                            ! (1/s)
!					
! Diagnostic variables:
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux   (kg/s m2)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s m2)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PPRLFLX! liquid precip flux  (m/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PPRSFLX! solid precip flux   (m/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURV   ! water vapor in updraft (kg/kg)
!real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCI  ! total condensate in updraft (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURC   !  grid-scale liq condensate in updraft (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURI   !  grid-scale ice condensate in updraft (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PCFD   ! Cloud Fraction - deep convection
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PAREAUP   ! area of updraft - deep convection (m^2)

!integer, DIMENSION(KLON),   INTENT(INOUT) :: KCLTOP ! cloud top level (number of model level)
!integer, DIMENSION(KLON),   INTENT(INOUT) :: KCLBAS ! cloud base level(number of model level)
                                                      ! they are given a value of
                                                      ! 0 if no convection

real, DIMENSION(KLON),       INTENT(INOUT) :: ZTOP ! cloud top height
real, DIMENSION(KLON),       INTENT(INOUT) :: ZBAS ! cloud base height
real, DIMENSION(KLON),       INTENT(OUT) :: PCAPE  ! CAPE (J/kg)
real, DIMENSION(KLON),       INTENT(OUT) :: PWUMAXOUT  ! maximum upward velocity in BKF (m/s)
!
! special for ERA40
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   ! updraft detrainment rate   (kg/s m3)
real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   ! downdraft detrainment rate (kg/s m3)
real, DIMENSION(KLON),      INTENT(OUT)   :: PKKFC  ! deep convective counter
!
real, DIMENSION(KLON),      INTENT(IN)   :: PSP   ! surface pressure
real, DIMENSION(KLON),      INTENT(IN)   :: XLAT  ! latitude
real, DIMENSION(KLON),      INTENT(IN)   :: MG    ! land_sea mask
real, DIMENSION(KLON),      INTENT(IN)   :: MLAC  ! lake mask



!*       0.2   Declarations of local variables :
!
integer  :: JI, JK, JKP, JN  ! loop index
!
real, DIMENSION(KLON)               :: PPRTEN
!
! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
! increase one extra level for Bechtold scheme
real, DIMENSION(KLON,KLEV+1) :: ZKKFC
real, DIMENSION(KLON,KLEV+1) :: ZT     ! grid scale T at time t  (K)
real, DIMENSION(KLON,KLEV+1) :: ZRV    ! grid scale water vapor  (kg/kg)
real, DIMENSION(KLON,KLEV+1) :: ZRC    ! grid scale r_c mixing ratio (kg/kg)
real, DIMENSION(KLON,KLEV+1) :: ZRI    ! grid scale r_i mixing ratio (kg/kg)
real, DIMENSION(KLON,KLEV+1) :: ZU     ! grid scale horiz. wind u (m/s)
real, DIMENSION(KLON,KLEV+1) :: ZV     ! grid scale horiz. wind v (m/s)
real, DIMENSION(KLON,KLEV+1) :: ZW     ! grid scale vertical velocity (m/s)
real, DIMENSION(KLON,KLEV+1) :: ZW1    ! perturbed vertical velocity for ensemble (m/s)
real, DIMENSION(KLON,KLEV+1) :: ZPABS  ! grid scale pressure (Pa)
real, DIMENSION(KLON,KLEV+1) :: ZZZ    ! height of model layer (m)
real, DIMENSION(KLON,KLEV+1) :: ZHSFLX ! turbulent sensible heat flux (W/m^2)
!
real, DIMENSION(KLON,KLEV+1) :: ZTTEN  ! convective temperat. tendency (K/s)
real, DIMENSION(KLON,KLEV+1) :: ZRVTEN ! convective r_v tendency (1/s)
real, DIMENSION(KLON,KLEV+1) :: ZRCTEN ! convective r_c tendency (1/s)
real, DIMENSION(KLON,KLEV+1) :: ZRITEN ! convective r_i tendency (1/s)
real, DIMENSION(KLON,KLEV+1) :: ZUTEN  ! convective u tendency (m/s^2)
real, DIMENSION(KLON,KLEV+1) :: ZVTEN  ! convective m tendency (m/s^2)
real, DIMENSION(KLON,KLEV+1) :: ZCFD   ! cloud fraction - deep
real, DIMENSION(KLON,KLEV+1) :: ZAREAUP   ! cloud area -deep  (m^2)
real, DIMENSION(KLON,KLEV+1) :: ZUMF   ! updraft mass flux   (kg/s m2)
real, DIMENSION(KLON,KLEV+1) :: ZDMF   ! downdraft mass flux (kg/s m2)
real, DIMENSION(KLON,KLEV+1) :: ZURV   ! water vapor in updrafts (kg/kg)
!real, DIMENSION(KLON,KLEV+1) :: ZURCI  ! total condensate in updrafts (kg/kg)
real, DIMENSION(KLON,KLEV+1) :: ZURC   ! liq condensate in updrafts (kg/kg)
real, DIMENSION(KLON,KLEV+1) :: ZURI   ! ice condensate in updrafts (kg/kg)
real, DIMENSION(KLON,KLEV+1) :: ZPRLFLX! liquid precip flux  (m/s)
real, DIMENSION(KLON,KLEV+1) :: ZPRSFLX! solid precip flux   (m/s)
!integer, DIMENSION(KLON)   :: ICLTOP ! cloud top level (number of model level)
!integer, DIMENSION(KLON)   :: ICLBAS ! cloud base level(number of model level)
real, DIMENSION(KLON,KLEV+1,KCH1):: ZCH1     ! grid scale chemical species
real, DIMENSION(KLON,KLEV+1,KCH1):: ZCH1TEN  ! chemical convective tendency

common/cnvmain/ZTTENE,  ZRVTENE, ZRCTENE,       &
ZRITENE, ZPRLTENE, ZPRSTENE, ZUMFE, ZDMFE,  ZURVE, ZURCE, ZURIE, ZPRLFLXE, &
ZPRSFLXE, ICLTOPE, ICLBASE,  ZUTENE, ZVTENE, ZCH1TENE, ZEDUMMY, IEDUMMY,     &
ZWEIGHT, ZUDRE, ZDDRE
!$OMP THREADPRIVATE(/cnvmain/)
#define ALLOCATABLE POINTER
!
!
!*       0.5   Declarations of additional Ensemble fields:
!
integer                             :: KENS     ! number of allowed
                                                  ! additional deep convection calls
real, DIMENSION(:,:,:), ALLOCATABLE :: ZTTENE   ! convective temperat. tendency (K/s)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZRVTENE  ! convective r_v tendency (1/s)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZRCTENE  ! convective r_c tendency (1/s)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZRITENE  ! convective r_i tendency (1/s)
real, DIMENSION(:,:),   ALLOCATABLE :: ZPRLTENE ! liquid surf precipitation tendency (kg/s m2)
real, DIMENSION(:,:),   ALLOCATABLE :: ZPRSTENE ! solid surf precipitation tendency (kg/s m2)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZUMFE    ! updraft mass flux   (kg/s m2)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZDMFE    ! downdraft mass flux (kg/s m2)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZURVE    ! updraft water vapor  (kg/kg)
!real, DIMENSION(:,:,:), ALLOCATABLE :: ZURCIE   ! updraft condensate   (kg/kg)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZURCE    ! updraft liq condensate   (kg/kg)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZURIE    ! updraft ice condensate   (kg/kg)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZPRLFLXE ! liquid precip flux  (kg/s m2)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZPRSFLXE ! solid precip flux   (kg/s m2)
real, DIMENSION(:,:),ALLOCATABLE :: ICLTOPE  ! cloud top level (m)
real, DIMENSION(:,:),ALLOCATABLE :: ICLBASE  ! cloud base level(m)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZUTENE   ! convective u tendency (m/s^2)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZVTENE   ! convective u tendency (m/s^2)
real, DIMENSION(:,:,:,:),ALLOCATABLE:: ZCH1TENE ! chemical convective tendency
real, DIMENSION(:),     ALLOCATABLE :: ZEDUMMY  ! field not to be recomputed by ensemble
integer, DIMENSION(:),  ALLOCATABLE :: IEDUMMY  ! field not to be recomputed by ensemble
real, DIMENSION(:),     ALLOCATABLE :: ZWEIGHT  ! weighting factor for ensemble members
real                                :: ZSUM     ! sum of weighting factors
!
integer                      :: KLEV1   ! vertical dimension (add sfc as first level)
! special for ERA40
real, DIMENSION(KLON,KLEV+1)          :: ZUDR, ZDDR  ! updraft/downdraft detrainment rates (kg/s m3)
real, DIMENSION(:,:,:), ALLOCATABLE :: ZUDRE, ZDDRE
!-------------------------------------------------------------------------------
!
!
!*       .9   Setup fundamental thermodunamical/physical constants using ECMWF/ARPEGE routine
!             ------------------------------------------------------------------------------
!
     CALL SUCST(54,20020211,0,0)
!
!
!*       1.   Allocate 2D (horizontal, vertical) arrays and additional ensemble arrays
!             ------------------------------------------------------------------------
!

!    KCLTOP(:)  = 1 ! set default value when no convection
!    KCLBAS(:)  = 1 ! can be changed  depending on user
!    ZTOP(:)=0.
!    ZBAS(:)=0.
    PUDR = 0.
    PDDR = 0.

!
!
KENS = MIN( KENSM, 3 )
IF ( KENS > 0 ) THEN
    ALLOCATE( ZTTENE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZRVTENE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZRCTENE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZRITENE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZUMFE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZDMFE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZURVE(KLON,KLEV+1,KENS) )
!    ALLOCATE( ZURCIE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZURCE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZURIE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZUTENE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZVTENE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZCH1TENE(KLON,KLEV+1,KCH1,KENS) )
    ALLOCATE( ZPRLFLXE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZPRSFLXE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZPRLTENE(KLON,KENS) )
    ALLOCATE( ZPRSTENE(KLON,KENS) )
    ALLOCATE( ICLTOPE(KLON,KENS) )
    ALLOCATE( ICLBASE(KLON,KENS) )
    ALLOCATE( ZEDUMMY(KLON) )
    ALLOCATE( IEDUMMY(KLON) )
    ALLOCATE( ZWEIGHT(KENS) )
    IEDUMMY(:)   = 0
    ICLTOPE(:,:) = 1 ! set default value when no convection
    ICLBASE(:,:) = 1

    ALLOCATE( ZUDRE(KLON,KLEV+1,KENS) )
    ALLOCATE( ZDDRE(KLON,KLEV+1,KENS) )
END IF
!
!
!*       2.   Flip arrays upside-down as  first vertical level in convection is 1
!             --------------------------------------------------------------------
!
! modifications to move whole column one level up, and the first level value is surface value
! level 1 --> level klev+1 ;  level klev --> level 2 ;
DO JK = 1, KLEV
   JKP = KLEV - JK + 2
   DO JI = KIDIA, KFDIA
      ZPABS(JI,JKP) = PPABS(JI,JK)
      ZZZ(JI,JKP)   = PZZ(JI,JK)
      ZT(JI,JKP)    = PT(JI,JK)
      ZRV(JI,JKP)   = PRV(JI,JK) / ( _ONE_ - PRV(JI,JK) ) ! transform specific humidity
      ZRC(JI,JKP)   = PRC(JI,JK) / ( _ONE_ - PRC(JI,JK) ) ! in mixing ratio
      ZRI(JI,JKP)   = PRI(JI,JK) / ( _ONE_ - PRI(JI,JK) )
      ZU(JI,JKP)    = PU(JI,JK)
      ZV(JI,JKP)    = PV(JI,JK)
      ZW(JI,JKP)    = PW(JI,JK)
      ZHSFLX(JI,JKP)= PHSFLX(JI,JK)
   END DO
END DO
DO JK = 1, KLEV
   JKP = KLEV - JK + 2
   DO JI = KIDIA, KFDIA
      ZTTEN(JI,JKP) = PTTEN(JI,JK)
      ZRVTEN(JI,JKP) = PRVTEN(JI,JK) * ( _ONE_ + ZRV(JI,JKP) ) ** 2  ! specific humidity tendency to mixing ratio tendency
      ZRCTEN(JI,JKP) = PRCTEN(JI,JK) * ( _ONE_ + ZRC(JI,JKP) ) ** 2
      ZRITEN(JI,JKP) = PRITEN(JI,JK) * ( _ONE_ + ZRI(JI,JKP) ) ** 2
      ZUTEN(JI,JKP)  = PUTEN(JI,JK)
      ZVTEN(JI,JKP)  = PVTEN(JI,JK)
      ZCFD(JI,JKP)   = PCFD(JI,JK)
      ZAREAUP(JI,JKP)= PAREAUP(JI,JK)
      ZUMF(JI,JKP)   = PUMF(JI,JK)
      ZDMF(JI,JKP)   = PDMF(JI,JK)
      !         ZURV(JI,JKP)   = PURV(JI,JK)  / ( _ONE_ - PURV(JI,JK) )
      !         ZURCI(JI,JKP)  = PURCI(JI,JK) / ( _ONE_ - PURCI(JI,JK) )
      ZURC(JI,JKP)   = PURC(JI,JK) / ( _ONE_ - PURC(JI,JK) )
      ZURI(JI,JKP)   = PURI(JI,JK) / ( _ONE_ - PURI(JI,JK) )
      ZPRLFLX(JI,JKP)= PPRLFLX(JI,JK)
      ZPRSFLX(JI,JKP)= PPRSFLX(JI,JK)
      ZUDR(JI,JKP)   = PUDR(JI,JK)
      ZDDR(JI,JKP)   = PDDR(JI,JK)
   END DO
END DO
DO JI = KIDIA, KFDIA
   ZTTEN(JI,1)  = 0.
   ZRVTEN(JI,1) = 0.
   ZRCTEN(JI,1) = 0.
   ZRITEN(JI,1) = 0.
   ZUTEN(JI,1)  = 0.
   ZVTEN(JI,1)  = 0.
   ZCFD(JI,1)   = 0.
   ZAREAUP(JI,1)= 0.
   ZUMF(JI,1)   = 0.
   ZDMF(JI,1)   = 0.
   ZURC(JI,1)   = 0.
   ZURI(JI,1)   = 0.
   ZPRLFLX(JI,1)= 0.
   ZPRSFLX(JI,1)= 0.
   ZUDR(JI,1)   = 0.
   ZDDR(JI,1)   = 0.
END DO
!PV:  sert a quoi????
!   ZPRLTEN=PPRTEN

! set first level value as surface values
 DO JI = KIDIA, KFDIA
   ZPABS(JI,1) = PSP(JI)
   ZZZ(JI,1)=0.
   ZT(JI,1)=PT(JI,KLEV)
   ZRV(JI,1)= 0.
   ZRC(JI,1)= 0.
   ZRI(JI,1)= 0.
!   ZRV(JI,1)= PRV(JI,KLEV)
!   ZRC(JI,1)= PRC(JI,KLEV)
!   ZRI(JI,1)= PRI(JI,KLEV)
   ZU(JI,1)=0.
   ZV(JI,1)=0.
   ZW(JI,1)=0.
   ZHSFLX(JI,1)=0.
 END DO


IF ( OCHTRANS ) THEN
   DO JK = 1, KLEV
      JKP = KLEV - JK + 2
      DO JN = 1, KCH1
         DO JI = KIDIA, KFDIA
            ZCH1(JI,JKP,JN) = PCH1(JI,JK,JN)
         END DO
      END DO
   END DO
END IF
!
!*       4.a  Call deep convection routine
!             ----------------------------
!
!
! 1. Base version
!
    CALL SU_CONVPAR
!
    KLEV1=KLEV+1
    CALL CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                 &
                        & PDTCONV, KICE, ODOWN,                                 &
                        & ZPABS, ZZZ, PDXDY, ZHSFLX,                            &
                        & ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW,                        &
                        & KCOUNT, ZTTEN, ZRVTEN, ZRCTEN, ZRITEN,                &
                        & ZPRLTEN, PPRSTEN,                                     &
                        & ZTOP, ZBAS, ZPRLFLX, ZPRSFLX,                         &
                        & ZUMF, ZDMF, ZURV, ZURC, ZURI, PCAPE, ZCFD,            &
                        & ZAREAUP, PWUMAXOUT,                                   &
                        & OUVTRANS, ZUTEN, ZVTEN,                               &
                        & OCHTRANS, KCH1, ZCH1, ZCH1TEN                         &
                        &,ZUDR, ZDDR, PKKFC, XLAT, MG, MLAC                     ) ! for ERA40
    if (phy_error_L) return

!  2. Additional Ensemble members
!
    IF ( KENS > 0 ) THEN
!
    CALL SU_CONVPAR1
!
!* first member - changes in YOE_CONVPAR (cloud radius of 500 m)
!                                          specified in SU_CONVPAR1
!
    CALL CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                                  &
                        & PDTCONV, KICE, ODOWN,                                                  &
                        & ZPABS, ZZZ, PDXDY, ZHSFLX,                                             &
                        & ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW,                                         &
                        & IEDUMMY, ZTTENE(:,:,1), ZRVTENE(:,:,1), ZRCTENE(:,:,1), ZRITENE(:,:,1),&
                        & ZPRLTENE(:,1), ZPRSTENE(:,1),                                          &
                        & ICLTOPE(:,1), ICLBASE(:,1), ZPRLFLXE(:,:,1), ZPRSFLXE(:,:,1),          &
                        & ZUMFE(:,:,1), ZDMFE(:,:,1), ZURVE(:,:,1), ZURCE(:,:,1), ZURIE(:,:,1), ZEDUMMY,ZCFD,      &
                        & ZAREAUP, PWUMAXOUT,                                                    &
                        & OUVTRANS, ZUTENE(:,:,1), ZVTENE(:,:,1),                                &
                        & OCHTRANS, KCH1, ZCH1, ZCH1TENE(:,:,:,1)                                &
                        &,ZUDRE(:,:,1), ZDDRE(:,:,1), ZKKFC, XLAT, MG, MLAC                      ) ! for ERA40
    if (phy_error_L) return
    END IF
!
    IF ( KENS > 1 ) THEN
!
    CALL SU_CONVPAR
!
!* second member (positive vertical velocity perturb for Trigger)
!
    ZW1=ZW*1.5+1.E-4
    CALL CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                                  &
                        & PDTCONV, KICE, ODOWN,                                                  &
                        & ZPABS, ZZZ, PDXDY, ZHSFLX,                                             &
                        & ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW1,                                        &
                        & IEDUMMY, ZTTENE(:,:,2), ZRVTENE(:,:,2), ZRCTENE(:,:,2), ZRITENE(:,:,2),&
                        & ZPRLTENE(:,2), ZPRSTENE(:,2),                                          &
                        & ICLTOPE(:,2), ICLBASE(:,2), ZPRLFLXE(:,:,2), ZPRSFLXE(:,:,2),          &
                        & ZUMFE(:,:,2), ZDMFE(:,:,2), ZURVE(:,:,2), ZURCE(:,:,2), ZURIE(:,:,2), ZEDUMMY,ZCFD,      &
                        & ZAREAUP, PWUMAXOUT,                                                    &
                        & OUVTRANS, ZUTENE(:,:,2), ZVTENE(:,:,2),                                &
                        & OCHTRANS, KCH1, ZCH1, ZCH1TENE(:,:,:,2)                                &
                        &,ZUDRE(:,:,2), ZDDRE(:,:,2), ZKKFC, XLAT, MG, MLAC                      )  ! for ERA40
    if (phy_error_L) return
    END IF
!
    IF ( KENS > 2 ) THEN
!
!* third member (negative vertical velocity perturb for Trigger)
!
    ZW1=ZW*0.5-1.E-4
    CALL CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                                  &
                        & PDTCONV, KICE, ODOWN,                                                  &
                        & ZPABS, ZZZ, PDXDY, ZHSFLX,                                             &
                        & ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW1,                                        &
                        & IEDUMMY, ZTTENE(:,:,3), ZRVTENE(:,:,3), ZRCTENE(:,:,3), ZRITENE(:,:,3),&
                        & ZPRLTENE(:,3), ZPRSTENE(:,3),                                          &
                        & ICLTOPE(:,3), ICLBASE(:,3), ZPRLFLXE(:,:,3), ZPRSFLXE(:,:,3),          &
                        & ZUMFE(:,:,3), ZDMFE(:,:,3), ZURVE(:,:,3), ZURCE(:,:,3), ZURIE(:,:,3), ZEDUMMY,ZCFD,      &
                        & ZAREAUP, PWUMAXOUT,                                                    &
                        & OUVTRANS, ZUTENE(:,:,3), ZVTENE(:,:,3),                                &
                        & OCHTRANS, KCH1, ZCH1, ZCH1TENE(:,:,:,3)                                &
                        &,ZUDRE(:,:,3), ZDDRE(:,:,3), ZKKFC, XLAT, MG, MLAC                      )  ! for ERA40
    if (phy_error_L) return
    END IF
!
!
!             -------------------------------
!
!
!*       5.  Add  - if activated - ensemble average values for deep
!            convective tendencies
!            ---------------------------------------------------------
!
ZSUM = _ONE_
IF ( KENS > 0 ) THEN
    IF ( KENS == 1 ) ZWEIGHT(:) = .5
    IF ( KENS >  1 ) ZWEIGHT(:) = _ONE_
    DO JN = 1, KENS
    DO JK = 1, KLEV1
    DO JI = KIDIA, KFDIA
       ZTTEN(JI,JK)  = ZTTEN(JI,JK)  + ZWEIGHT(JN) * ZTTENE(JI,JK,JN)
       ZRVTEN(JI,JK) = ZRVTEN(JI,JK) + ZWEIGHT(JN) * ZRVTENE(JI,JK,JN)
       ZRCTEN(JI,JK) = ZRCTEN(JI,JK) + ZWEIGHT(JN) * ZRCTENE(JI,JK,JN)
       ZRITEN(JI,JK) = ZRITEN(JI,JK) + ZWEIGHT(JN) * ZRITENE(JI,JK,JN)
       ZPRLFLX(JI,JK)= ZPRLFLX(JI,JK)+ ZWEIGHT(JN) * ZPRLFLXE(JI,JK,JN)
       ZPRSFLX(JI,JK)= ZPRSFLX(JI,JK)+ ZWEIGHT(JN) * ZPRSFLXE(JI,JK,JN)
       ZUMF(JI,JK)   = ZUMF(JI,JK)   + ZWEIGHT(JN) * ZUMFE(JI,JK,JN)
       ZDMF(JI,JK)   = ZDMF(JI,JK)   + ZWEIGHT(JN) * ZDMFE(JI,JK,JN)
       ZURV(JI,JK)   = ZURV(JI,JK)   + ZWEIGHT(JN) * ZURVE(JI,JK,JN)
!       ZURCI(JI,JK)  = ZURCI(JI,JK)  + ZWEIGHT(JN) * ZURCIE(JI,JK,JN)
       ZURC(JI,JK)   = ZURC(JI,JK)   + ZWEIGHT(JN) * ZURCE(JI,JK,JN)
       ZURI(JI,JK)   = ZURI(JI,JK)   + ZWEIGHT(JN) * ZURIE(JI,JK,JN)

       ZUDR(JI,JK)   = ZUDR(JI,JK)   + ZWEIGHT(JN) * ZUDRE(JI,JK,JN) ! for ERA40
       ZDDR(JI,JK)   = ZDDR(JI,JK)   + ZWEIGHT(JN) * ZDDRE(JI,JK,JN)
    END DO
    END DO
    DO JI = KIDIA, KFDIA
       ZPRLTEN(JI)  = ZPRLTEN(JI)  + ZWEIGHT(JN) * ZPRLTENE(JI,JN)
       PPRSTEN(JI)  = PPRSTEN(JI)  + ZWEIGHT(JN) * ZPRSTENE(JI,JN)
       ZTOP(JI)   = MAX(ZTOP(JI), ICLTOPE(JI,JN))
       ZBAS(JI)   = MAX(ZBAS(JI), ICLBASE(JI,JN))
    END DO
       IF ( OUVTRANS ) THEN
         DO JK = 1, KLEV1
         DO JI = KIDIA, KFDIA
            ZUTEN(JI,JK) = ZUTEN(JI,JK) + ZWEIGHT(JN) * ZUTENE(JI,JK,JN)
            ZVTEN(JI,JK) = ZVTEN(JI,JK) + ZWEIGHT(JN) * ZVTENE(JI,JK,JN)
         END DO
         END DO
       END IF
       IF ( OCHTRANS )  THEN
         DO JK = 1, KLEV1
         DO JI = KIDIA, KFDIA
          ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) + ZWEIGHT(JN) * ZCH1TENE(JI,JK,:,JN)
         END DO
         END DO
       END IF
    END DO
!
    ZSUM = _ONE_ / ( _ONE_ + SUM( ZWEIGHT(:) ) )
END IF
!

    DO JI = KIDIA, KFDIA
       PPRTEN(JI)   = ( ZPRLTEN(JI) + PPRSTEN(JI) ) * ZSUM
       PPRSTEN(JI)  = PPRSTEN(JI)  * ZSUM
    END DO
       IF ( OUVTRANS ) THEN
         DO JK = 1, KLEV1
         DO JI = KIDIA, KFDIA
            ZUTEN(JI,JK) = ZUTEN(JI,JK) * ZSUM
            ZVTEN(JI,JK) = ZVTEN(JI,JK) * ZSUM
         END DO
         END DO
       END IF
       IF ( OCHTRANS ) THEN
         DO JK = 1, KLEV1
         DO JI = KIDIA, KFDIA
!            ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) * ZSUM + ZCH1TENS(JI,JK,:)
            ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) * ZSUM
         END DO
         END DO
       END IF
!
!
!*       6.  Reflip arrays to ECMWF/ARPEGE vertical structure
!            change mixing ratios to specific humidity
!
! shift the whole column one level down
! level 1 <-- level klev+1 ;  level klev <-- level 2 ;
DO JK = 1, KLEV
   JKP = KLEV - JK + 2
   DO JI = KIDIA, KFDIA
      PTTEN(JI,JK)  = ZTTEN(JI,JKP)
   ! don't transform back to specific hum, does not conserve integrals
      PRVTEN(JI,JK) = ZRVTEN(JI,JKP) / ( _ONE_ + ZRV(JI,JKP) ) ** 2
      PRCTEN(JI,JK) = ZRCTEN(JI,JKP) / ( _ONE_ + ZRC(JI,JKP) ) ** 2
      PRITEN(JI,JK) = ZRITEN(JI,JKP) / ( _ONE_ + ZRI(JI,JKP) ) ** 2
      PUTEN(JI,JK)  = ZUTEN(JI,JKP)
      PVTEN(JI,JK)  = ZVTEN(JI,JKP)
      PCFD(JI,JK)   = ZCFD(JI,JKP)
      PAREAUP(JI,JK)= ZAREAUP(JI,JKP)
      PUMF(JI,JK)   = ZUMF(JI,JKP)
      PDMF(JI,JK)   = ZDMF(JI,JKP)
!      PURV(JI,JK)   = ZURV(JI,JKP) / ( _ONE_ + ZURV(JI,JKP) )
      PURC(JI,JK)   = ZURC(JI,JKP)/ ( _ONE_ + ZURC(JI,JKP) )
      PURI(JI,JK)   = ZURI(JI,JKP)/ ( _ONE_ + ZURI(JI,JKP) )
      PPRLFLX(JI,JK)= ZPRLFLX(JI,JKP)
      PPRSFLX(JI,JK)= ZPRSFLX(JI,JKP)
      PUDR(JI,JK)   = ZUDR(JI,JKP)
      PDDR(JI,JK)   = ZDDR(JI,JKP)
   END DO
END DO


!DO JI = KIDIA, KFDIA
!   JK = ICLTOP(JI)
!   KCLTOP(JI) = KLEV - JK + 1
!   JK = ICLBAS(JI)
!   KCLBAS(JI) = KLEV - JK + 1
!   IF ( ICLTOP(JI) == 1 ) KCLTOP(JI) = 1
!   IF ( ICLBAS(JI) == 1 ) KCLBAS(JI) = 1
!PV use deep precip as flag for deep conv activity
!PV output heights of cloud top and base
!    IF( ZPRLTEN(JI) > _ZERO_ ) THEN
!       ZTOP(JI)=ZZZ(JI,ICLTOP(JI))
!       ZBAS(JI)=ZZZ(JI,ICLBAS(JI))
!    ELSE
!      ZTOP(JI)=0.
!      ZBAS(JI)=0.
!    ENDIF
!END DO

IF ( OCHTRANS ) THEN
   DO JK = 1, KLEV
      JKP = KLEV - JK + 2
      DO JN = 1, KCH1
         DO JI = KIDIA, KFDIA
            PCH1TEN(JI,JK,JN) = ZCH1TEN(JI,JKP,JN)
         END DO
      END DO
   END DO
END IF

!*       7.  Deallocate local arrays
!


IF ( KENS > 0 ) THEN
       DEALLOCATE( ZTTENE )
       DEALLOCATE( ZRVTENE )
       DEALLOCATE( ZRCTENE )
       DEALLOCATE( ZRITENE )
       DEALLOCATE( ZUTENE )
       DEALLOCATE( ZVTENE )
       DEALLOCATE( ZUMFE )
       DEALLOCATE( ZDMFE )
       DEALLOCATE( ZURVE )
!       DEALLOCATE( ZURCIE )
       DEALLOCATE( ZURCE )
       DEALLOCATE( ZURIE )
       DEALLOCATE( ZCH1TENE )
       DEALLOCATE( ZPRLFLXE )
       DEALLOCATE( ZPRSFLXE )
       DEALLOCATE( ZPRLTENE )
       DEALLOCATE( ZPRSTENE )
       DEALLOCATE( ZEDUMMY )
       DEALLOCATE( IEDUMMY )
       DEALLOCATE( ICLTOPE )
       DEALLOCATE( ICLBASE )
       DEALLOCATE( ZWEIGHT )

       DEALLOCATE( ZUDRE )
       DEALLOCATE( ZDDRE )
END IF
!
!
END SUBROUTINE BKF_DEEP3
!
!----------------------------------------------------------------------------
