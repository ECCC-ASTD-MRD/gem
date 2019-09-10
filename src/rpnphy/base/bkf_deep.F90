subroutine bkf_deep4(KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA, &
     & PDTCONV,  ODOWN, KICE,                              &
     & KENSM,                                              &
     & PPABS, PZZ, PDXDY, PHSFLX,                          &
     & PT, PRV, PRC, PRI, PU, PV, PW,                      &
     & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,              &
     & ZPRLTEN, PPRSTEN,                                   &
     & PUMF, PDMF, PPRLFLX, PPRSFLX, PCAPE, ZTOP, ZBAS,    &
     & PURC, PURI, PCFD,                                   &
     & PAREAUP,  PWUMAXOUT,                                &
     & OUVTRANS, PUTEN, PVTEN,                             &
     & OCHTRANS, KCH1, PCH1, PCH1TEN,                      &
     & PUDR, PDDR, PKKFC, PSP, XLAT, MG, MLAC)
   use phy_status, only: phy_error_L
   implicit none
!!!#include <arch_specific.hf>

   !!**** Interface routine to the fast Meso-NH convection code developed for ECMWF/ARPEGE
   !!     having a structure typical for operational routines
   !
   !!     Transformations necessary to call deep+ shallow code
   !!     - skip input vertical arrays/levels : bottom=1, top=KLEV
   !!     - transform specific humidities in mixing ratio
   !
   !
   !!    PURPOSE
   !!    -------
   !!      The routine interfaces the MNH convection code as developed for operational
   !!      forecast models like ECMWF/ARPEGE or HIRLAM with the typical Meso-NH array structure
   !!      Calls the deep and/or shallow convection routine
   !
   !
   !!**  METHOD
   !!    ------
   !!     Returns one tendency for shallow+deep convection but each part can
   !!     be activated/desactivated separately
   !!     For deep convection one can enable up to 3 additional ensemble members
   !!     - this substantially improves the smoothness of the scheme and reduces
   !!       allows for runs with different cloud radii (entrainment rates) and
   !!       reduces the arbitrariness inherent to convective trigger condition
   !
   !
   !
   !!    EXTERNAL
   !!    --------
   !!    CONVECT_DEEP
   !!    SU_CONVPAR, SU_CONVPAR1
   !!    SUCST:   ECMWF/ARPEGE routine
   !
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !
   !
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    11/12/98
   !!      modified    20/03/2002 by P. Marquet : transformed for ARPEGE/Climat
   !!                             (tsmbkind.h, REAL_B, INTEGER_M, _JPRB, "& &",
   !!                              _ZERO_, _ONE_, _HALF_)
   !!      modified    11/04/O2 allow for ensemble of deep updrafts/downdrafts
   !
   !!    REFERENCE
   !!    ---------
   !!    Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc., Vol 127, pp 869-886:
   !!           A mass flux convection scheme for regional and global models.
   !
   !-------------------------------------------------------------------------------
   !Revisions in RPN
   ! 001   -PV-jan2015 - do not aggregate deep and shallow tendencies; add output of shallow tendencies
   ! 002   -PV-apr2015 - output cloud top and base heights for deep only
   ! 003   -JY-Jul2015 - separate from convection.ftn90 (which has both deep and shallow convection scheme)
   !                   - add output of UA(area of updraft), K6(maximum upward velocity in Bechtold scheme)

   !*       0.    DECLARATIONS
   !              ------------

#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

   !*       0.1   Declarations of dummy arguments :

   integer,                    intent(IN)   :: KLON   ! horizontal dimension
   integer,                    intent(IN)   :: KLEV   ! vertical dimension
   integer,                    intent(IN)   :: KIDIA  ! value of the first point in x
   integer,                    intent(IN)   :: KFDIA  ! value of the last point in x
   integer,                    intent(IN)   :: KBDIA  ! vertical  computations start at
   !                                                    ! KBDIA that is at least 1
   integer,                    intent(IN)   :: KTDIA  ! vertical computations can be
   ! limited to KLEV + 1 - KTDIA
   ! default=1
   real,                       intent(IN)   :: PDTCONV! Interval of time between two
   ! calls of the deep convection
   ! scheme
   logical,                      intent(IN)   :: ODOWN  ! take or not convective
   ! downdrafts into account
   integer,                    intent(IN)   :: KICE   ! flag for ice ( 1 = yes,
   !                0 = no ice )
   integer,                    intent(IN)   :: KENSM  ! number of additional deep convection calls
   ! for ensemble (presently limited to 3)
   ! KENSM=0 corresponds to base run with
   ! 1 deep and 1 shallow call

   real, dimension(KLON,KLEV), intent(IN)   :: PT     ! grid scale T at time t  (K)
   real, dimension(KLON,KLEV), intent(IN)   :: PRV    ! grid scale water vapor  (kg/kg)
   real, dimension(KLON,KLEV), intent(IN)   :: PRC    ! grid scale r_c (kg/kg)
   real, dimension(KLON,KLEV), intent(IN)   :: PRI    ! grid scale r_i (kg/kg)
   real, dimension(KLON,KLEV), intent(IN)   :: PU     ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV), intent(IN)   :: PV     ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV), intent(IN)   :: PW     ! grid scale vertical velocity (m/s)
   real, dimension(KLON,KLEV), intent(IN)   :: PPABS  ! grid scale pressure (Pa)
   real, dimension(KLON,KLEV), intent(IN)   :: PZZ    ! geopotential (m2/s2)
   real, dimension(KLON),      intent(IN)   :: PDXDY  ! grid area (m2)
   real, dimension(KLON,KLEV), intent(IN)   :: PHSFLX ! turbulent sensible heat flux (W/m^2)
   integer, dimension(KLON), intent(INOUT)  :: KCOUNT ! convective counter(recompute
   ! tendency or keep it

   real, dimension(KLON,KLEV), intent(INOUT):: PTTEN  ! convective temperat. tendency (K/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRVTEN ! convective r_v tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRCTEN ! convective r_c tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRITEN ! convective r_i tendency (1/s)
   !real, DIMENSION(KLON),      INTENT(INOUT):: PPRTEN ! total surf precipitation tendency (m/s)
   real, dimension(KLON),      intent(INOUT):: PPRSTEN! solid surf precipitation tendency (kg/s m2)
   real, dimension(KLON),      intent(INOUT):: ZPRLTEN! liquid surf precipitation tendency (kg/s m2)

   ! convective U,V transport (Cu friction)
   logical,                      intent(IN)        :: OUVTRANS ! flag to compute convective
   ! transport for chemical tracer
   real, dimension(KLON,KLEV), intent(INOUT)     :: PUTEN    ! convecctive u tendency (m/s^2)
   real, dimension(KLON,KLEV), intent(INOUT)     :: PVTEN    ! convecctive v tendency (m/s^2)

   ! convective Chemical Tracers:
   logical,                      intent(IN)        :: OCHTRANS ! flag to compute convective
   ! transport for chemical tracer
   integer,                    intent(IN)        :: KCH1     ! number of species
   real, dimension(KLON,KLEV,KCH1), intent(IN)   :: PCH1     ! grid scale chemical species
   real, dimension(KLON,KLEV,KCH1), intent(INOUT):: PCH1TEN  ! chemical convective tendency
   ! (1/s)

   ! Diagnostic variables:
   real, dimension(KLON,KLEV), intent(INOUT) :: PUMF   ! updraft mass flux   (kg/s m2)
   real, dimension(KLON,KLEV), intent(INOUT) :: PDMF   ! downdraft mass flux (kg/s m2)
   real, dimension(KLON,KLEV), intent(INOUT) :: PPRLFLX! liquid precip flux  (m/s)
   real, dimension(KLON,KLEV), intent(INOUT) :: PPRSFLX! solid precip flux   (m/s)
!!$   real, dimension(KLON,KLEV), intent(INOUT) :: PURV   ! water vapor in updraft (kg/kg)   !#TODO: never used
   !real, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCI  ! total condensate in updraft (kg/kg)
   real, dimension(KLON,KLEV), intent(INOUT) :: PURC   !  grid-scale liq condensate in updraft (kg/kg)
   real, dimension(KLON,KLEV), intent(INOUT) :: PURI   !  grid-scale ice condensate in updraft (kg/kg)
   real, dimension(KLON,KLEV), intent(INOUT) :: PCFD   ! Cloud Fraction - deep convection
   real, dimension(KLON,KLEV), intent(INOUT) :: PAREAUP   ! area of updraft - deep convection (m^2)

   !integer, DIMENSION(KLON),   INTENT(INOUT) :: KCLTOP ! cloud top level (number of model level)
   !integer, DIMENSION(KLON),   INTENT(INOUT) :: KCLBAS ! cloud base level(number of model level)
   ! they are given a value of
   ! 0 if no convection

   real, dimension(KLON),       intent(INOUT) :: ZTOP ! cloud top height
   real, dimension(KLON),       intent(INOUT) :: ZBAS ! cloud base height
   real, dimension(KLON),       intent(OUT) :: PCAPE  ! CAPE (J/kg)
   real, dimension(KLON),       intent(OUT) :: PWUMAXOUT  ! maximum upward velocity in BKF (m/s)

   ! special for ERA40
   real, dimension(KLON,KLEV), intent(INOUT) :: PUDR   ! updraft detrainment rate   (kg/s m3)
   real, dimension(KLON,KLEV), intent(INOUT) :: PDDR   ! downdraft detrainment rate (kg/s m3)
   real, dimension(KLON),      intent(OUT)   :: PKKFC  ! deep convective counter

   real, dimension(KLON),      intent(IN)   :: PSP   ! surface pressure
   real, dimension(KLON),      intent(IN)   :: XLAT  ! latitude
   real, dimension(KLON),      intent(IN)   :: MG    ! land_sea mask
   real, dimension(KLON),      intent(IN)   :: MLAC  ! lake mask


   !*       0.2   Declarations of local variables :

   integer  :: JI, JK, JKP, JN  ! loop index

   real, dimension(KLON)               :: PPRTEN

   ! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
   ! increase one extra level for Bechtold scheme
   real, dimension(KLON,KLEV+1) :: ZKKFC
   real, dimension(KLON,KLEV+1) :: ZT     ! grid scale T at time t  (K)
   real, dimension(KLON,KLEV+1) :: ZRV    ! grid scale water vapor  (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZRC    ! grid scale r_c mixing ratio (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZRI    ! grid scale r_i mixing ratio (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZU     ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV+1) :: ZV     ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV+1) :: ZW     ! grid scale vertical velocity (m/s)
   real, dimension(KLON,KLEV+1) :: ZW1    ! perturbed vertical velocity for ensemble (m/s)
   real, dimension(KLON,KLEV+1) :: ZPABS  ! grid scale pressure (Pa)
   real, dimension(KLON,KLEV+1) :: ZZZ    ! height of model layer (m)
   real, dimension(KLON,KLEV+1) :: ZHSFLX ! turbulent sensible heat flux (W/m^2)

   real, dimension(KLON,KLEV+1) :: ZTTEN  ! convective temperat. tendency (K/s)
   real, dimension(KLON,KLEV+1) :: ZRVTEN ! convective r_v tendency (1/s)
   real, dimension(KLON,KLEV+1) :: ZRCTEN ! convective r_c tendency (1/s)
   real, dimension(KLON,KLEV+1) :: ZRITEN ! convective r_i tendency (1/s)
   real, dimension(KLON,KLEV+1) :: ZUTEN  ! convective u tendency (m/s^2)
   real, dimension(KLON,KLEV+1) :: ZVTEN  ! convective m tendency (m/s^2)
   real, dimension(KLON,KLEV+1) :: ZCFD   ! cloud fraction - deep
   real, dimension(KLON,KLEV+1) :: ZAREAUP   ! cloud area -deep  (m^2)
   real, dimension(KLON,KLEV+1) :: ZUMF   ! updraft mass flux   (kg/s m2)
   real, dimension(KLON,KLEV+1) :: ZDMF   ! downdraft mass flux (kg/s m2)
   real, dimension(KLON,KLEV+1) :: ZURV   ! water vapor in updrafts (kg/kg)
   !real, DIMENSION(KLON,KLEV+1) :: ZURCI  ! total condensate in updrafts (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZURC   ! liq condensate in updrafts (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZURI   ! ice condensate in updrafts (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZPRLFLX! liquid precip flux  (m/s)
   real, dimension(KLON,KLEV+1) :: ZPRSFLX! solid precip flux   (m/s)
   !integer, DIMENSION(KLON)   :: ICLTOP ! cloud top level (number of model level)
   !integer, DIMENSION(KLON)   :: ICLBAS ! cloud base level(number of model level)
   real, dimension(KLON,KLEV+1,KCH1):: ZCH1     ! grid scale chemical species
   real, dimension(KLON,KLEV+1,KCH1):: ZCH1TEN  ! chemical convective tendency

   common/cnvmain/ZTTENE,  ZRVTENE, ZRCTENE,       &
        ZRITENE, ZPRLTENE, ZPRSTENE, ZUMFE, ZDMFE,  ZURVE, ZURCE, ZURIE, ZPRLFLXE, &
        ZPRSFLXE, ICLTOPE, ICLBASE,  ZUTENE, ZVTENE, ZCH1TENE, ZEDUMMY, IEDUMMY,     &
        ZWEIGHT, ZUDRE, ZDDRE
!$OMP THREADPRIVATE(/cnvmain/)


   !*       0.5   Declarations of additional Ensemble fields:

   integer                             :: KENS     ! number of allowed
   ! additional deep convection calls
   real, dimension(:,:,:), pointer :: ZTTENE   ! convective temperat. tendency (K/s)
   real, dimension(:,:,:), pointer :: ZRVTENE  ! convective r_v tendency (1/s)
   real, dimension(:,:,:), pointer :: ZRCTENE  ! convective r_c tendency (1/s)
   real, dimension(:,:,:), pointer :: ZRITENE  ! convective r_i tendency (1/s)
   real, dimension(:,:),   pointer :: ZPRLTENE ! liquid surf precipitation tendency (kg/s m2)
   real, dimension(:,:),   pointer :: ZPRSTENE ! solid surf precipitation tendency (kg/s m2)
   real, dimension(:,:,:), pointer :: ZUMFE    ! updraft mass flux   (kg/s m2)
   real, dimension(:,:,:), pointer :: ZDMFE    ! downdraft mass flux (kg/s m2)
   real, dimension(:,:,:), pointer :: ZURVE    ! updraft water vapor  (kg/kg)
   !real, DIMENSION(:,:,:), POINTER :: ZURCIE   ! updraft condensate   (kg/kg)
   real, dimension(:,:,:), pointer :: ZURCE    ! updraft liq condensate   (kg/kg)
   real, dimension(:,:,:), pointer :: ZURIE    ! updraft ice condensate   (kg/kg)
   real, dimension(:,:,:), pointer :: ZPRLFLXE ! liquid precip flux  (kg/s m2)
   real, dimension(:,:,:), pointer :: ZPRSFLXE ! solid precip flux   (kg/s m2)
   real, dimension(:,:),pointer :: ICLTOPE  ! cloud top level (m)
   real, dimension(:,:),pointer :: ICLBASE  ! cloud base level(m)
   real, dimension(:,:,:), pointer :: ZUTENE   ! convective u tendency (m/s^2)
   real, dimension(:,:,:), pointer :: ZVTENE   ! convective u tendency (m/s^2)
   real, dimension(:,:,:,:),pointer:: ZCH1TENE ! chemical convective tendency
   real, dimension(:),     pointer :: ZEDUMMY  ! field not to be recomputed by ensemble
   integer, dimension(:),  pointer :: IEDUMMY  ! field not to be recomputed by ensemble
   real, dimension(:),     pointer :: ZWEIGHT  ! weighting factor for ensemble members
   real                                :: ZSUM     ! sum of weighting factors

   integer                      :: KLEV1   ! vertical dimension (add sfc as first level)
   ! special for ERA40
   real, dimension(KLON,KLEV+1)          :: ZUDR, ZDDR  ! updraft/downdraft detrainment rates (kg/s m3)
   real, dimension(:,:,:), pointer :: ZUDRE, ZDDRE
   !-------------------------------------------------------------------------------


   !*       .9   Setup fundamental thermodunamical/physical constants using ECMWF/ARPEGE routine
   !             ------------------------------------------------------------------------------

   call SUCST(54,20020211,0,0)


   !*       1.   Allocate 2D (horizontal, vertical) arrays and additional ensemble arrays
   !             ------------------------------------------------------------------------


   !    KCLTOP(:)  = 1 ! set default value when no convection
   !    KCLBAS(:)  = 1 ! can be changed  depending on user
   !    ZTOP(:)=0.
   !    ZBAS(:)=0.
   PUDR = 0.
   PDDR = 0.



   KENS = min( KENSM, 3 )
   if ( KENS > 0 ) then
      allocate( ZTTENE(KLON,KLEV+1,KENS) )
      allocate( ZRVTENE(KLON,KLEV+1,KENS) )
      allocate( ZRCTENE(KLON,KLEV+1,KENS) )
      allocate( ZRITENE(KLON,KLEV+1,KENS) )
      allocate( ZUMFE(KLON,KLEV+1,KENS) )
      allocate( ZDMFE(KLON,KLEV+1,KENS) )
      allocate( ZURVE(KLON,KLEV+1,KENS) )
      !    ALLOCATE( ZURCIE(KLON,KLEV+1,KENS) )
      allocate( ZURCE(KLON,KLEV+1,KENS) )
      allocate( ZURIE(KLON,KLEV+1,KENS) )
      allocate( ZUTENE(KLON,KLEV+1,KENS) )
      allocate( ZVTENE(KLON,KLEV+1,KENS) )
      allocate( ZCH1TENE(KLON,KLEV+1,KCH1,KENS) )
      allocate( ZPRLFLXE(KLON,KLEV+1,KENS) )
      allocate( ZPRSFLXE(KLON,KLEV+1,KENS) )
      allocate( ZPRLTENE(KLON,KENS) )
      allocate( ZPRSTENE(KLON,KENS) )
      allocate( ICLTOPE(KLON,KENS) )
      allocate( ICLBASE(KLON,KENS) )
      allocate( ZEDUMMY(KLON) )
      allocate( IEDUMMY(KLON) )
      allocate( ZWEIGHT(KENS) )
      IEDUMMY(:)   = 0
      ICLTOPE(:,:) = 1 ! set default value when no convection
      ICLBASE(:,:) = 1

      allocate( ZUDRE(KLON,KLEV+1,KENS) )
      allocate( ZDDRE(KLON,KLEV+1,KENS) )
   end if


   !*       2.   Flip arrays upside-down as  first vertical level in convection is 1
   !             --------------------------------------------------------------------

   ! modifications to move whole column one level up, and the first level value is surface value
   ! level 1 --> level klev+1 ;  level klev --> level 2 ;
   do JK = 1, KLEV
      JKP = KLEV - JK + 2
      do JI = KIDIA, KFDIA
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
      end do
   end do
   do JK = 1, KLEV
      JKP = KLEV - JK + 2
      do JI = KIDIA, KFDIA
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
      end do
   end do
   do JI = KIDIA, KFDIA
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
   end do
   !PV:  sert a quoi????
   !   ZPRLTEN=PPRTEN

   ! set first level value as surface values
   do JI = KIDIA, KFDIA
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
   end do


   if ( OCHTRANS ) then
      do JK = 1, KLEV
         JKP = KLEV - JK + 2
         do JN = 1, KCH1
            do JI = KIDIA, KFDIA
               ZCH1(JI,JKP,JN) = PCH1(JI,JK,JN)
            end do
         end do
      end do
   end if

   !*       4.a  Call deep convection routine
   !             ----------------------------


   ! 1. Base version

   call SU_CONVPAR

   KLEV1=KLEV+1
   call CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                 &
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

   if ( KENS > 0 ) then

      call SU_CONVPAR1

      !* first member - changes in YOE_CONVPAR (cloud radius of 500 m)
      !                                          specified in SU_CONVPAR1

      call CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                                  &
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
   end if

   if ( KENS > 1 ) then

      call SU_CONVPAR

      !* second member (positive vertical velocity perturb for Trigger)

      ZW1=ZW*1.5+1.E-4
      call CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                                  &
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
   end if

   if ( KENS > 2 ) then

      !* third member (negative vertical velocity perturb for Trigger)

      ZW1=ZW*0.5-1.E-4
      call CONVECT_DEEP2(KLON, KLEV1, KIDIA, KFDIA, KBDIA, KTDIA,                                  &
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
   end if


   !             -------------------------------


   !*       5.  Add  - if activated - ensemble average values for deep
   !            convective tendencies
   !            ---------------------------------------------------------

   ZSUM = _ONE_
   if ( KENS > 0 ) then
      if ( KENS == 1 ) ZWEIGHT(:) = .5
      if ( KENS >  1 ) ZWEIGHT(:) = _ONE_
      do JN = 1, KENS
         do JK = 1, KLEV1
            do JI = KIDIA, KFDIA
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
            end do
         end do
         do JI = KIDIA, KFDIA
            ZPRLTEN(JI)  = ZPRLTEN(JI)  + ZWEIGHT(JN) * ZPRLTENE(JI,JN)
            PPRSTEN(JI)  = PPRSTEN(JI)  + ZWEIGHT(JN) * ZPRSTENE(JI,JN)
            ZTOP(JI)   = max(ZTOP(JI), ICLTOPE(JI,JN))
            ZBAS(JI)   = max(ZBAS(JI), ICLBASE(JI,JN))
         end do
         if ( OUVTRANS ) then
            do JK = 1, KLEV1
               do JI = KIDIA, KFDIA
                  ZUTEN(JI,JK) = ZUTEN(JI,JK) + ZWEIGHT(JN) * ZUTENE(JI,JK,JN)
                  ZVTEN(JI,JK) = ZVTEN(JI,JK) + ZWEIGHT(JN) * ZVTENE(JI,JK,JN)
               end do
            end do
         end if
         if ( OCHTRANS )  then
            do JK = 1, KLEV1
               do JI = KIDIA, KFDIA
                  ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) + ZWEIGHT(JN) * ZCH1TENE(JI,JK,:,JN)
               end do
            end do
         end if
      end do

      ZSUM = _ONE_ / ( _ONE_ + sum( ZWEIGHT(:) ) )
   end if


   do JI = KIDIA, KFDIA
      PPRTEN(JI)   = ( ZPRLTEN(JI) + PPRSTEN(JI) ) * ZSUM
      PPRSTEN(JI)  = PPRSTEN(JI)  * ZSUM
   end do
   if ( OUVTRANS ) then
      do JK = 1, KLEV1
         do JI = KIDIA, KFDIA
            ZUTEN(JI,JK) = ZUTEN(JI,JK) * ZSUM
            ZVTEN(JI,JK) = ZVTEN(JI,JK) * ZSUM
         end do
      end do
   end if
   if ( OCHTRANS ) then
      do JK = 1, KLEV1
         do JI = KIDIA, KFDIA
            !            ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) * ZSUM + ZCH1TENS(JI,JK,:)
            ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) * ZSUM
         end do
      end do
   end if


   !*       6.  Reflip arrays to ECMWF/ARPEGE vertical structure
   !            change mixing ratios to specific humidity

   ! shift the whole column one level down
   ! level 1 <-- level klev+1 ;  level klev <-- level 2 ;
   do JK = 1, KLEV
      JKP = KLEV - JK + 2
      do JI = KIDIA, KFDIA
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
      end do
   end do


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

   if ( OCHTRANS ) then
      do JK = 1, KLEV
         JKP = KLEV - JK + 2
         do JN = 1, KCH1
            do JI = KIDIA, KFDIA
               PCH1TEN(JI,JK,JN) = ZCH1TEN(JI,JKP,JN)
            end do
         end do
      end do
   end if

   !*       7.  Deallocate local arrays



   if ( KENS > 0 ) then
      deallocate( ZTTENE )
      deallocate( ZRVTENE )
      deallocate( ZRCTENE )
      deallocate( ZRITENE )
      deallocate( ZUTENE )
      deallocate( ZVTENE )
      deallocate( ZUMFE )
      deallocate( ZDMFE )
      deallocate( ZURVE )
      !       DEALLOCATE( ZURCIE )
      deallocate( ZURCE )
      deallocate( ZURIE )
      deallocate( ZCH1TENE )
      deallocate( ZPRLFLXE )
      deallocate( ZPRSFLXE )
      deallocate( ZPRLTENE )
      deallocate( ZPRSTENE )
      deallocate( ZEDUMMY )
      deallocate( IEDUMMY )
      deallocate( ICLTOPE )
      deallocate( ICLBASE )
      deallocate( ZWEIGHT )

      deallocate( ZUDRE )
      deallocate( ZDDRE )
   end if

end subroutine bkf_deep4
