!############################################################################
  SUBROUTINE CONVECT_DEEP2(KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,          & !#TODO: rename deep3 after changing args list
                           & PDTCONV, KICE, ODOWN,                          &
                           & PPABST, PZZ, PDXDY, PHSFLX,                    &
                           & PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                           & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                           & PPRLTEN, PPRSTEN,                              &
                           & KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              &
                           & PUMF, PDMF, PURV, PURC, PURI, PCAPE, PCF,      &
                           & PAREAUP, PWUMAXOUT,                            &
                           & OUVCONV, PUTEN, PVTEN,                         &
                           & OCH1CONV, KCH1, PCH1, PCH1TEN,                 &
                           & PUDR, PDDR, PKKFC, XLAT, MG, MLAC              ) ! for ERA40
!############################################################################

!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays.
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!
!!
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_FUNCT
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_TSTEP_PREF
!!    CONVECT_DOWNDRAFT
!!    CONVECT_PRECIP_ADJUST
!!    CONVECT_CLOSURE
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                   ! gravity constant
!!          RPI                  ! number Pi
!!          RATM                 ! reference pressure
!!          RD, RV               ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV           ! specific heat for dry air and water vapor
!!          RALPW, RBETW, RGAMW  ! constants for water saturation pressure
!!          RTT                  ! triple point temperature
!!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
!!          RCW, RCS             ! specific heat for liquid water and ice
!!          RHOH2O               ! density of liquid water (normally defined in Module YOETHF )
!!                               ! but redefine dhere as RATM/100
!!
!!      Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module YOE_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc. :
!!           A mass flux convection scheme for regional and global models.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 04/10/97 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!-------------------------------------------------------------------------------


!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAREXT
USE YOE_CONVPAR
use cnv_options, deep_name=>deep

IMPLICIT NONE
#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :


integer,                    INTENT(IN) :: KLON     ! horizontal dimension
integer,                    INTENT(IN) :: KLEV     ! vertical dimension
integer,                    INTENT(IN) :: KIDIA    ! value of the first point in x
integer,                    INTENT(IN) :: KFDIA    ! value of the last point in x
integer,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                    ! KBDIA that is at least 1
integer,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                     ! limited to KLEV + 1 - KTDIA
                                                     ! default=1
real,                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                     ! calls of the deep convection
                                                     ! scheme
integer,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes,
                                                     !                0 = no ice )
LOGICAL,                      INTENT(IN) :: ODOWN    ! take or not convective
                                                     ! downdrafts into account
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature (K)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor(kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u (m/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v (m/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical
                                                     ! velocity (m/s)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
real, DIMENSION(KLON),      INTENT(IN) :: PDXDY    ! horizontal grid area (m2)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PHSFLX   ! turbulent heat flux (W/m2)

integer, DIMENSION(KLON),   INTENT(INOUT):: KCOUNT ! convective counter (recompute
                                                     ! tendency or keep it)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature
                                                     ! tendency (K/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
real, DIMENSION(KLON),      INTENT(INOUT):: PPRLTEN! liquid surf. precipitation
                                                     ! tendency (kg/s m2)
real, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf. precipitation
                                                     ! tendency (kg/s m2)
real, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level (m)
real, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level (m)
                                                     ! they are given a value of
                                                     ! 0 if no convection
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRLFLX! liquid precip flux (m/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRSFLX! solid  precip flux (m/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF   ! downdraft mass flux (kg/s m2)
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURV   ! updraft water vapor (kg/kg)
!real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURCI  ! updraft liquid+ice condensate (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURC   ! updraft liquid condensate (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURI   ! updraft ice condensate (kg/kg)
real, DIMENSION(KLON),      INTENT(INOUT):: PCAPE  ! maximum CAPE (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(INOUT)  :: PCF    ! Cloud Fraction
real, DIMENSION(KLON,KLEV), INTENT(INOUT)  :: PAREAUP! area of updraft, cloud coverage area (m^2)= PCF x dxdy
real, DIMENSION(KLON),      INTENT(INOUT)  :: PWUMAXOUT ! maximum upward velocity in bechtold deep convection

LOGICAL,                      INTENT(IN)   :: OUVCONV! include wind transport (Cu friction)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUTEN  ! convective u tendency (m/s^2)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PVTEN  ! convective v tendency (m/s^2)

LOGICAL,                      INTENT(IN) :: OCH1CONV ! include tracer transport
integer,                    INTENT(IN) :: KCH1     ! number of species
real, DIMENSION(KLON,KLEV,KCH1), INTENT(IN)   :: PCH1    ! grid scale chemical species
real, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN ! species conv. tendency (1/s)

! for ERA40
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PUDR   ! updraft detrainment rate (kg/s m2)
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PDDR   ! downdraft detrainment rate (kg/s m2)
real, DIMENSION(KLON),      INTENT(OUT)  :: PKKFC  ! shallow convective counter


real, DIMENSION(KLON),      INTENT(IN)  :: XLAT    ! Latitude
real, DIMENSION(KLON),      INTENT(IN)  :: MG      ! land-sea mask
real, DIMENSION(KLON),      INTENT(IN)  :: MLAC    ! lake mask


!*       0.2   Declarations of local fixed memory variables :

integer  :: ITEST, ICONV, ICONV1    ! number of convective columns
integer  :: IIB, IIE                ! horizontal loop bounds
integer  :: IKB, IKE                ! vertical loop bounds
integer  :: IKS                     ! vertical dimension
integer  :: JI, JL                  ! horizontal loop index
integer  :: JN                      ! number of tracers
integer  :: JK, JKP, JKM            ! vertical loop index
integer  :: IFTSTEPS                ! only used for chemical tracers
real     :: RHOH2O                  ! density of liquid water (normally defined in Module YOETHF )
real     :: ZEPS, ZEPSA, ZEPSB      ! R_d / R_v, R_v / R_d, RCPV / RCPD - ZEPSA
real     :: ZCPORD, ZRDOCP          ! C_p/R_d,  R_d/C_p
LOGICAL, PARAMETER :: CONSERV_CORRECT=.true. !apply conservation correction

LOGICAL, DIMENSION(KLON, KLEV)        :: GTRIG3 ! 3D logical mask for convection
LOGICAL, DIMENSION(KLON)              :: GTRIG  ! 2D logical mask for trigger test
real,  DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta,
                                                              ! theta_v, theta_es
real,  DIMENSION(KLON)              :: ZTIME  ! convective time period
real,  DIMENSION(KLON)              :: ZTIMC   ! fractional convective time step
integer,  DIMENSION(KLON)           :: ITSTEP  ! # of fractional convective timesteps
LOGICAL,  DIMENSION(KLON)             :: GWORK1  ! flag for newly activated columns including MASS CONSERVATION criteria from closure
real,  DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array
real                                :: ZW1 ! work

common/deep/IDPL    ,IPBL    ,ILCL    ,IETL    ,ICTL    ,ILFS    ,IDBL, IDDT , &
IML     ,ISDPL   ,ISPBL   ,ISLCL   ,ZSTHLCL ,ZSTLCL  ,ZSRVLCL ,ZSWLCL  ,&
ZSZLCL  ,ZSTHVELCL,ZSDXDY  ,ZZ      ,ZPRES   ,ZDPRES  ,ZU      ,ZV      ,&
ZW      ,ZHSFLX  ,ZRHO   ,ZTT     ,ZTH     ,ZTHV    ,ZTHL    ,ZTHES, ZTHEST ,ZRW     ,&
ZRV     ,ZRC     ,ZRI     ,ZDXDY   ,ZUMF    ,ZUER    ,ZUDR    ,ZUPR    ,&
ZUTHL   ,ZUTHV   ,ZWU     ,ZURW    ,ZURC    ,ZURI    ,ZURR    ,ZURS    ,ZUTPR   ,&
ZMFLCL  ,ZCAPE   ,ZTHLCL  ,ZTLCL   ,ZRVLCL  ,ZWLCL   ,ZZLCL   ,ZTHVELCL,&
ZDMF    ,ZDER    ,ZDDR    ,ZDTHL   ,ZDRW    ,ZMIXF   ,ZTPR    ,ZSPR    ,&
ZDTEVR  ,ZPREF   ,ZDTEVRF ,ZPRLFLX ,ZPRSFLX ,ZLMASS  ,ZTIMEA  ,ZTIMEC  ,&
ZTHC    ,ZRVC    ,ZRCC    ,ZRIC    ,ZWSUB   ,GTRIG1  ,GWORK   ,         &
IINDEX, IJINDEX, IJSINDEX, IJPINDEX,ZCPH    ,ZLV, ZLS,ZUC     ,ZVC     ,&
ZCH1    , ZCH1C   , ZWORK3  , GTRIG4
!$OMP THREADPRIVATE(/deep/)
#define ALLOCATABLE POINTER
!*       0.2   Declarations of local allocatable  variables :

integer, DIMENSION(:),ALLOCATABLE  :: IDPL    ! index for parcel departure level
integer, DIMENSION(:),ALLOCATABLE  :: IPBL    ! index for source layer top
integer, DIMENSION(:),ALLOCATABLE  :: ILCL    ! index for lifting condensation level
integer, DIMENSION(:),ALLOCATABLE  :: IETL    ! index for zero buoyancy level
integer, DIMENSION(:),ALLOCATABLE  :: ICTL    ! index for cloud top level
integer, DIMENSION(:),ALLOCATABLE  :: ILFS    ! index for level of free sink
integer, DIMENSION(:),ALLOCATABLE  :: IDBL    ! index for downdraft base level
integer, DIMENSION(:),ALLOCATABLE  :: IDDT    ! index for downdraft detrainment top level
integer, DIMENSION(:),ALLOCATABLE  :: IML     ! melting level

integer, DIMENSION(:), ALLOCATABLE :: ISDPL   ! index for parcel departure level
integer, DIMENSION(:), ALLOCATABLE :: ISPBL   ! index for source layer top
integer, DIMENSION(:), ALLOCATABLE :: ISLCL   ! index for lifting condensation level

real, DIMENSION(:), ALLOCATABLE    :: ZSTHLCL ! updraft theta at LCL
real, DIMENSION(:), ALLOCATABLE    :: ZSTLCL  ! updraft temp. at LCL
real, DIMENSION(:), ALLOCATABLE    :: ZSRVLCL ! updraft rv at LCL
real, DIMENSION(:), ALLOCATABLE    :: ZSWLCL  ! updraft w at LCL
real, DIMENSION(:), ALLOCATABLE    :: ZSZLCL  ! LCL height
real, DIMENSION(:), ALLOCATABLE    :: ZSTHVELCL! envir. theta_v at LCL
real, DIMENSION(:), ALLOCATABLE    :: ZSDXDY  ! grid area (m^2)

! grid scale variables
real, DIMENSION(:,:), ALLOCATABLE  :: ZZ      ! height of model layer (m)
real, DIMENSION(:,:), ALLOCATABLE  :: ZPRES   ! grid scale pressure
real, DIMENSION(:,:), ALLOCATABLE  :: ZDPRES  ! pressure difference between
                                              ! bottom and top of layer (Pa)
real, DIMENSION(:,:), ALLOCATABLE  :: ZU      ! grid scale horiz. u component on theta grid
real, DIMENSION(:,:), ALLOCATABLE  :: ZV      ! grid scale horiz. v component on theta grid
real, DIMENSION(:,:), ALLOCATABLE  :: ZW      ! grid scale vertical velocity on theta grid
real, DIMENSION(:,:), ALLOCATABLE  :: ZHSFLX  ! turbulent sensible heat flux (W/m2)
real, DIMENSION(:,:), ALLOCATABLE  :: ZRHO    ! air density
real, DIMENSION(:,:), ALLOCATABLE  :: ZTT     ! temperature
real, DIMENSION(:,:), ALLOCATABLE  :: ZTH     ! grid scale theta
real, DIMENSION(:,:), ALLOCATABLE  :: ZTHV    ! grid scale theta_v
real, DIMENSION(:,:), ALLOCATABLE  :: ZTHL    ! grid scale enthalpy (J/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZTHES, ZTHEST ! grid scale saturated theta_e
real, DIMENSION(:,:), ALLOCATABLE  :: ZRW     ! grid scale total water (kg/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZRV     ! grid scale water vapor (kg/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZRC     ! grid scale cloud water (kg/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZRI     ! grid scale cloud ice (kg/kg)
real, DIMENSION(:),   ALLOCATABLE  :: ZDXDY   ! grid area (m^2)

! updraft variables
real, DIMENSION(:,:), ALLOCATABLE  :: ZUMF    ! updraft mass flux (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZUER    ! updraft entrainment (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZUDR    ! updraft detrainment (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZUPR    ! updraft precipitation in
                                              ! flux units (kg water / s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
real, DIMENSION(:,:), ALLOCATABLE  :: ZWU     ! updraft vertical vel (m/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZURR    ! liquid precipit. (kg/kg)
                                              ! produced in  model layer
real, DIMENSION(:,:), ALLOCATABLE  :: ZURS    ! solid precipit. (kg/kg)
                                              ! produced in  model layer
real, DIMENSION(:),   ALLOCATABLE  :: ZUTPR   ! total updraft precipitation (kg/s)
real, DIMENSION(:),   ALLOCATABLE  :: ZMFLCL  ! cloud base unit mass flux(kg/s)
real, DIMENSION(:),   ALLOCATABLE  :: ZCAPE   ! available potent. energy
real, DIMENSION(:),   ALLOCATABLE  :: ZTHLCL  ! updraft theta at LCL
real, DIMENSION(:),   ALLOCATABLE  :: ZTLCL   ! updraft temp. at LCL
real, DIMENSION(:),   ALLOCATABLE  :: ZRVLCL  ! updraft rv at LCL
real, DIMENSION(:),   ALLOCATABLE  :: ZWLCL   ! updraft w at LCL
real, DIMENSION(:),   ALLOCATABLE  :: ZZLCL   ! LCL height
real, DIMENSION(:),   ALLOCATABLE  :: ZTHVELCL! envir. theta_v at LCL

! downdraft variables
real, DIMENSION(:,:), ALLOCATABLE  :: ZDMF    ! downdraft mass flux (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZDER    ! downdraft entrainment (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZDDR    ! downdraft detrainment (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZDTHL   ! downdraft enthalpy (J/kg)
real, DIMENSION(:,:), ALLOCATABLE  :: ZDRW    ! downdraft total water (kg/kg)
real, DIMENSION(:),   ALLOCATABLE  :: ZMIXF   ! mixed fraction at LFS
real, DIMENSION(:),   ALLOCATABLE  :: ZTPR    ! total surf precipitation (kg/s)
real, DIMENSION(:),   ALLOCATABLE  :: ZSPR    ! solid surf precipitation (kg/s)
real, DIMENSION(:),   ALLOCATABLE  :: ZDTEVR  ! donwndraft evapor. (kg/s)
real, DIMENSION(:),   ALLOCATABLE  :: ZPREF   ! precipitation efficiency
real, DIMENSION(:,:), ALLOCATABLE  :: ZDTEVRF ! donwndraft evapor. (kg/s)
real, DIMENSION(:,:), ALLOCATABLE  :: ZPRLFLX ! liquid precip flux
real, DIMENSION(:,:), ALLOCATABLE  :: ZPRSFLX ! solid precip flux

! closure variables
real, DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
real, DIMENSION(:),   ALLOCATABLE  :: ZTIMEA  ! advective time period
real, DIMENSION(:),   ALLOCATABLE  :: ZTIMEC  ! time during which convection is
                                                ! active at grid point (as ZTIME)

real, DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
real, DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w
real, DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c
real, DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i
real, DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)

LOGICAL,   DIMENSION(:),ALLOCATABLE  :: GTRIG1  ! logical mask for convection
LOGICAL,   DIMENSION(:),ALLOCATABLE  :: GWORK   ! logical work array
integer, DIMENSION(:),ALLOCATABLE  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index

real, DIMENSION(:),   ALLOCATABLE  :: ZCPH    ! specific heat C_ph
real, DIMENSION(:),   ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
real                               :: ZES     ! saturation vapor mixng ratio

! for U, V transport:
real,  DIMENSION(:,:), ALLOCATABLE :: ZUC     ! horizontal wind u (m/s)
real,  DIMENSION(:,:), ALLOCATABLE :: ZVC     ! horizontal wind v (m/s)

! for Chemical Tracer transport:
real,  DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
real,  DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
real,  DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! work arrays
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE:: GTRIG4  ! logical mask

!-------------------------------------------------------------------------------
RHOH2O = RATM *1.E-2

!*       0.3    Compute loop bounds
!               -------------------

IIB    = KIDIA
IIE    = KFDIA
JCVEXB = MAX( 0, KBDIA - 1 )
IKB    = 1 + JCVEXB
IKS    = KLEV
JCVEXT = MAX( 0, KTDIA - 1 )
IKE    = IKS - JCVEXT


!*       0.5    Update convective counter ( where KCOUNT > 0
!               convection is still active ).
!               ---------------------------------------------

KCOUNT(IIB:IIE) = max(KCOUNT(IIB:IIE) - 1,0)

GTRIG(:)  = KCOUNT(:) <= 0
ITEST     = COUNT( GTRIG(:) )

WHERE ( KCOUNT(:) > 0 )
  PKKFC(:) = 1.
END WHERE


!PV - must undo change in units done in ccdiagnostics.ftn90 to all precip rates
!PV - or else precip rates will divided by 1000 at each timestep during convective timesacle
PPRLTEN(:) = PPRLTEN(:)*1000.0
PPRSTEN(:) = PPRSTEN(:)*1000.0

IF ( ITEST == 0 ) RETURN  ! if convection is already active at every grid point
                          ! exit CONVECT_DEEP



!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------

GTRIG3(:,:) = SPREAD( GTRIG(:), DIM=2, NCOPIES=IKS )
WHERE ( GTRIG3(:,:) )
    PTTEN(:,:)  = _ZERO_
    PRVTEN(:,:) = _ZERO_
    PRCTEN(:,:) = _ZERO_
    PRITEN(:,:) = _ZERO_
    PPRLFLX(:,:)= _ZERO_
    PPRSFLX(:,:)= _ZERO_
    PUTEN(:,:)  = _ZERO_
    PVTEN(:,:)  = _ZERO_
    PCF(:,:)   = _ZERO_
    PAREAUP(:,:)   = _ZERO_
    PUMF(:,:)   = _ZERO_
    PDMF(:,:)   = _ZERO_
    PURV(:,:)   = _ZERO_
!    PURCI(:,:)  = _ZERO_
    PURC(:,:)   = _ZERO_
    PURI(:,:)   = _ZERO_
    PUDR(:,:)   = _ZERO_
    PDDR(:,:)   = _ZERO_
END WHERE

WHERE ( GTRIG(:) )
   PPRLTEN(:) = _ZERO_
   PPRSTEN(:) = _ZERO_
   KCLTOP(:)  = _ZERO_
   KCLBAS(:)  = _ZERO_
   PCAPE(:)   = _ZERO_
   PWUMAXOUT(:)=_ZERO_
END WHERE

IF ( OCH1CONV ) THEN
   ALLOCATE( GTRIG4(KLON,KLEV,KCH1) )
   GTRIG4(:,:,:) = SPREAD( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
   WHERE( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = _ZERO_
   DEALLOCATE( GTRIG4 )
ENDIF


!*       1.     Initialize  local variables
!               ----------------------------

ZEPS   = RD / RV
ZEPSA  = RV / RD
ZEPSB  = RCPV / RCPD - ZEPSA
ZCPORD = RCPD / RD
ZRDOCP = RD / RCPD


!*       1.1    Set up grid scale theta, theta_v, theta_es
!               ------------------------------------------

ZTHT(:,:)   = 300.
ZSTHV(:,:)  = 300.
ZSTHES(:,:) = 400.

DO JK = IKB, IKE
DO JI = IIB, IIE
   IF ( PPABST(JI,JK) > 40.E2 ) THEN
      ZTHT(JI,JK)  = PTT(JI,JK) * ( RATM / PPABST(JI,JK) ) ** ZRDOCP
      ZSTHV(JI,JK) = ZTHT(JI,JK) * ( _ONE_ + ZEPSA * PRVT(JI,JK) ) /           &
                   & ( _ONE_ + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )

          ! use conservative Bolton (1980) formula for theta_e
          ! it is used to compute CAPE for undilute parcel ascent
          ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here

      ZES = EXP( RALPW - RBETW / PTT(JI,JK) - RGAMW * LOG( PTT(JI,JK) ) )
      ZES = MIN( _ONE_, ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
      ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **         &
                    & ( _ONE_ - 0.28 * ZES ) *                        &
                    &   EXP( ( 3374.6525 / PTT(JI,JK) - 2.5403 ) &
                    &   * ZES * ( _ONE_ + 0.81 * ZES ) )
   ENDIF
ENDDO
ENDDO



!*       2.     Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------

!*       2.1    Allocate arrays depending on number of model columns that need
!               to be tested for convection (i.e. where no convection is present
!               at the moment.
!               --------------------------------------------------------------

     ALLOCATE( ZPRES(ITEST,IKS) )
     ALLOCATE( ZZ(ITEST,IKS) )
     ALLOCATE( ZW(ITEST,IKS) )
     ALLOCATE( ZTH(ITEST,IKS) )
     ALLOCATE( ZTHV(ITEST,IKS) )
     ALLOCATE( ZTHEST(ITEST,IKS) )
     ALLOCATE( ZRV(ITEST,IKS) )
     ALLOCATE( ZSTHLCL(ITEST) )
     ALLOCATE( ZSTLCL(ITEST) )
     ALLOCATE( ZSRVLCL(ITEST) )
     ALLOCATE( ZSWLCL(ITEST) )
     ALLOCATE( ZSZLCL(ITEST) )
     ALLOCATE( ZSTHVELCL(ITEST) )
     ALLOCATE( ISDPL(ITEST) )
     ALLOCATE( ISPBL(ITEST) )
     ALLOCATE( ISLCL(ITEST) )
     ALLOCATE( ZSDXDY(ITEST) )
     ALLOCATE( GTRIG1(ITEST) )
     ALLOCATE( ZCAPE(ITEST) )
     ALLOCATE( IINDEX(KLON) )
     ALLOCATE( IJSINDEX(ITEST) )
     ALLOCATE( ZHSFLX(ITEST,IKS) )
     DO JI = 1, KLON
        IINDEX(JI) = JI
     ENDDO
     IJSINDEX(:) = PACK( IINDEX(:), MASK=GTRIG(:) )

  DO JK = IKB, IKE
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZPRES(JI,JK)  = PPABST(JL,JK)
     ZZ(JI,JK)     = PZZ(JL,JK)
     ZTH(JI,JK)    = ZTHT(JL,JK)
     ZTHV(JI,JK)   = ZSTHV(JL,JK)
     ZTHEST(JI,JK) = ZSTHES(JL,JK)
     ZRV(JI,JK)    = MAX( _ZERO_, PRVT(JL,JK) )
     ZW(JI,JK)     = PWT(JL,JK)
     ZHSFLX(JI,JK) = PHSFLX(JL,JK)
  ENDDO
  ENDDO
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZSDXDY(JI)    = PDXDY(JL)
  ENDDO

!*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c
!               and envir. saturation theta_e
!               ------------------------------------------------------------


!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------

     ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL
     ISDPL(:) = IKB
     ISPBL(:) = IKB


     CALL CONVECT_TRIGGER_FUNCT2(ITEST, KLEV,                              &
                               & ZPRES, ZTH, ZTHV, ZTHEST,                 &
                               & ZRV, ZW, ZZ, ZSDXDY, ZHSFLX,              &
                               & ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                               & ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1,   &
                               & ZCAPE, XLAT, MG, MLAC )

!PV cape in columns rejected by criteria in trigger must be put to zero
     WHERE ( .NOT. GTRIG1(:) )
         ZCAPE(:)   = _ZERO_
     END WHERE

     DO JI = 1, ITEST
        JL = IJSINDEX(JI)
        PCAPE(JL) = ZCAPE(JI)
     ENDDO

     DEALLOCATE( ZPRES )
     DEALLOCATE( ZZ )
     DEALLOCATE( ZTH )
     DEALLOCATE( ZTHV )
     DEALLOCATE( ZTHEST )
     DEALLOCATE( ZRV )
     DEALLOCATE( ZW )
     DEALLOCATE( ZCAPE )
     DEALLOCATE( ZHSFLX )


!*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
!               arrays used in the convection scheme using the mask GTRIG, i.e.
!               we do calculus only in convective columns. This corresponds to
!               a GATHER operation.
!               --------------------------------------------------------------

     ICONV = COUNT( GTRIG1(:) )
     IF ( ICONV == 0 )  THEN
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( GTRIG1 )
         DEALLOCATE( IINDEX )
         DEALLOCATE( IJSINDEX )
         RETURN   ! no convective column has been found, exit CONVECT_DEEP
     ENDIF

     ! vertical index variables

         ALLOCATE( IDPL(ICONV) )
         ALLOCATE( IPBL(ICONV) )
         ALLOCATE( ILCL(ICONV) )
         ALLOCATE( ICTL(ICONV) )
         ALLOCATE( IETL(ICONV) )

         ! grid scale variables

         ALLOCATE( ZZ(ICONV,IKS) )
         ALLOCATE( ZPRES(ICONV,IKS) )
         ALLOCATE( ZDPRES(ICONV,IKS) )
         ALLOCATE( ZU(ICONV,IKS) )
         ALLOCATE( ZV(ICONV,IKS) )
         ALLOCATE( ZRHO(ICONV, IKS) )
         ALLOCATE( ZTT(ICONV, IKS) )
         ALLOCATE( ZTH(ICONV,IKS) )
         ALLOCATE( ZTHV(ICONV,IKS) )
         ALLOCATE( ZTHL(ICONV,IKS) )
         ALLOCATE( ZTHES(ICONV,IKS) )
         ALLOCATE( ZRV(ICONV,IKS) )
         ALLOCATE( ZRC(ICONV,IKS) )
         ALLOCATE( ZRI(ICONV,IKS) )
         ALLOCATE( ZRW(ICONV,IKS) )
         ALLOCATE( ZDXDY(ICONV) )

         ! updraft variables

         ALLOCATE( ZUMF(ICONV,IKS) )
         ALLOCATE( ZUER(ICONV,IKS) )
         ALLOCATE( ZUDR(ICONV,IKS) )
         ALLOCATE( ZUPR(ICONV,IKS) )
         ALLOCATE( ZUTHL(ICONV,IKS) )
         ALLOCATE( ZUTHV(ICONV,IKS) )
         ALLOCATE( ZWU(ICONV,IKS) )
         ALLOCATE( ZURW(ICONV,IKS) )
         ALLOCATE( ZURC(ICONV,IKS) )
         ALLOCATE( ZURI(ICONV,IKS) )
         ALLOCATE( ZURR(ICONV,IKS) )
         ALLOCATE( ZURS(ICONV,IKS) )
         ALLOCATE( ZUTPR(ICONV) )
         ALLOCATE( ZTHLCL(ICONV) )
         ALLOCATE( ZTLCL(ICONV) )
         ALLOCATE( ZRVLCL(ICONV) )
         ALLOCATE( ZWLCL(ICONV) )
         ALLOCATE( ZMFLCL(ICONV) )
         ALLOCATE( ZZLCL(ICONV) )
         ALLOCATE( ZTHVELCL(ICONV) )
         ALLOCATE( ZCAPE(ICONV) )

         ! work variables

         ALLOCATE( IJINDEX(ICONV) )
         ALLOCATE( IJPINDEX(ICONV) )
         ALLOCATE( ZCPH(ICONV) )
         ALLOCATE( ZLV(ICONV) )
         ALLOCATE( ZLS(ICONV) )


!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------

         GTRIG(:)      = UNPACK( GTRIG1(:), MASK=GTRIG(:), FIELD=.FALSE. )
         IJINDEX(:)    = PACK( IINDEX(:), MASK=GTRIG(:) )

    DO JK = IKB, IKE
    DO JI = 1, ICONV
         JL = IJINDEX(JI)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZTT(JI,JK)    = PTT(JL,JK)
         ZRHO(JI,JK)   = ZPRES(JI,JK)/(ZTT(JI,JK)*RD)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHES(JI,JK)  = ZSTHES(JL,JK)
         ZRV(JI,JK)    = MAX( _ZERO_, PRVT(JL,JK) )
         ZRC(JI,JK)    = MAX( _ZERO_, PRCT(JL,JK) )
         ZRI(JI,JK)    = MAX( _ZERO_, PRIT(JL,JK) )
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
         ZU(JI,JK)     = PUT(JL,JK)
         ZV(JI,JK)     = PVT(JL,JK)
    ENDDO
    ENDDO

    DO JI = 1, ITEST
       IJSINDEX(JI) = JI
    ENDDO
    IJPINDEX(:) = PACK( IJSINDEX(:), MASK=GTRIG1(:) )
    DO JI = 1, ICONV
         JL = IJPINDEX(JI)
         IDPL(JI)      = ISDPL(JL)
         IPBL(JI)      = ISPBL(JL)
         ILCL(JI)      = ISLCL(JL)
         ZTHLCL(JI)    = ZSTHLCL(JL)
         ZTLCL(JI)     = ZSTLCL(JL)
         ZRVLCL(JI)    = ZSRVLCL(JL)
         ZWLCL(JI)     = ZSWLCL(JL)
         ZZLCL(JI)     = ZSZLCL(JL)
         ZTHVELCL(JI)  = ZSTHVELCL(JL)
         ZDXDY(JI)     = ZSDXDY(JL)
    ENDDO
         ALLOCATE( GWORK(ICONV) )
         GWORK(:)      = PACK( GTRIG1(:),  MASK=GTRIG1(:) )
         DEALLOCATE( GTRIG1 )
         ALLOCATE( GTRIG1(ICONV) )
         GTRIG1(:)     = GWORK(:)

         DEALLOCATE( GWORK )
!         DEALLOCATE( IJPINDEX )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )


!*           3.2    Compute pressure difference
!                   ---------------------------------------------------

        ZDPRES(:,IKB) = _ZERO_
        DO JK = IKB + 1, IKE
            ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK)
        ENDDO

!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c
!                  ----------------------------------------------------------

        DO JK = IKB, IKE, 1
            ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
            ZCPH(:)    = RCPD + RCPV * ZRW(:,JK)
            ZLV(:)     = RLVTT + ( RCPV - RCW ) * ( ZTT(:,JK) - RTT ) ! compute L_v
            ZLS(:)     = RLSTT + ( RCPV - RCS ) * ( ZTT(:,JK) - RTT ) ! compute L_i
            ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( _ONE_ + ZRW(:,JK) ) * RG * ZZ(:,JK) &
                       & - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
        ENDDO


!*           4.     Compute updraft properties
!                   ----------------------------

!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s )
!                   -------------------------------------------------------------

         DO JI = 1, ICONV
               JK = ILCL(JI) - 1
               ZMFLCL(JI) = ZPRES(JI,JK) / ( RD * ZTT(JI,JK) *                    &
                          & ( _ONE_ + ZEPS * ZRVLCL(JI) ) ) * RPI * XCRAD * XCRAD &
                          & * MAX( _ONE_, ZDXDY(JI)/XA25 )
         ENDDO

         DEALLOCATE( ZCPH )
         DEALLOCATE( ZLV )
         DEALLOCATE( ZLS )


     CALL CONVECT_UPDRAFT( ICONV, KLEV,                                     &
                         & KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                         & ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   &
                         & ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                         & ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZWU,             &
                         & ZURW, ZURC, ZURI, ZURR, ZURS, ZUPR,              &
                         & ZUTPR, ZCAPE, ICTL, IETL                         )

! print '(a50,5i3,x,l8)','UPDRAFT LEVELS:lc,kpbl,klcl,let,ltop:',idpl,ipbl,ilcl,ietl,ictl,gtrig1






!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------


     ICONV1 = COUNT(GTRIG1)

!PV CAPE in columns rejected by criteria in updraft must be put to zero
     DO JI = 1, ICONV
          JL = IJPINDEX(JI)
          IF(.NOT. GTRIG1(JI) ) PCAPE(JL)   = _ZERO_
     ENDDO
     DEALLOCATE( IJPINDEX )


     IF ( ICONV1 > 0 )  THEN


!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------

! downdraft variables

        ALLOCATE( ILFS(ICONV) )
        ALLOCATE( IDBL(ICONV) )
        ALLOCATE( IDDT(ICONV) )
        ALLOCATE( IML(ICONV) )
        ALLOCATE( ZDMF(ICONV,IKS) )
        ALLOCATE( ZDER(ICONV,IKS) )
        ALLOCATE( ZDDR(ICONV,IKS) )
        ALLOCATE( ZDTHL(ICONV,IKS) )
        ALLOCATE( ZDRW(ICONV,IKS) )
        ALLOCATE( ZLMASS(ICONV,IKS) )
        DO JK = IKB, IKE
           ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / RG  ! mass of model layer
        ENDDO
        ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
        ALLOCATE( ZMIXF(ICONV) )
        ALLOCATE( ZTPR(ICONV) )
        ALLOCATE( ZSPR(ICONV) )
        ALLOCATE( ZDTEVR(ICONV) )
        ALLOCATE( ZPREF(ICONV) )
        ALLOCATE( ZDTEVRF(ICONV,IKS) )
        ALLOCATE( ZPRLFLX(ICONV,IKS) )
        ALLOCATE( ZPRSFLX(ICONV,IKS) )

! closure variables

        ALLOCATE( ZTIMEA(ICONV) )
        ALLOCATE( ZTIMEC(ICONV) )
        ALLOCATE( ZTHC(ICONV,IKS) )
        ALLOCATE( ZRVC(ICONV,IKS) )
        ALLOCATE( ZRCC(ICONV,IKS) )
        ALLOCATE( ZRIC(ICONV,IKS) )
        ALLOCATE( ZWSUB(ICONV,IKS) )


!*           5.     Compute downdraft properties
!                   ----------------------------

!*           5.1    Compute advective time period and precipitation
!                   efficiency as a function of mean ambient wind (shear)
!                   --------------------------------------------------------

        CALL CONVECT_TSTEP_PREF( ICONV, KLEV,                          &
                               & ZU, ZV, ZPRES, ZZ, ZDXDY, ILCL, ICTL, &
                               & ZTIMEA, ZPREF )

          ! exclude convective downdrafts if desired
        IF ( .NOT. ODOWN ) ZPREF(:) = _ONE_

        ! Compute the period during which convection is active
        if (deep_timeconv_sec > 0.) then
           ztimec(:) = max(pdtconv,deep_timeconv_sec)
        elseif (deep_timeconv == 'ADVECTIVE') then
           ZTIMEC(:) = MAX( 1800., MIN( 3600., ZTIMEA(:) ) )
           ZTIMEC(:) = REAL( NINT( ZTIMEC(:) / PDTCONV ) ) * PDTCONV
           ZTIMEC(:) = MAX( PDTCONV, ZTIMEC(:) ) ! necessary if PDTCONV > 1800
        else
           call physeterror('convect_deep', 'adjustment time '//trim(deep_timeconv)//' unsupported by '//trim(deep_name))
           return
        endif

!*           5.2    Compute melting level
!                   ----------------------

        IML(:) = IKB
        DO JK = IKE, IKB, -1
          WHERE( ZTT(:,JK) <= RTT )  IML(:) = JK
        ENDDO

        CALL CONVECT_DOWNDRAFT2( ICONV, KLEV,                               &
                              & KICE, ZPRES, ZDPRES, ZZ, ZTH, ZTHES,       &
                              & ZRW, ZRC, ZRI,                             &
                              & ZPREF, ILCL, ICTL, IETL,                   &
                              & ZUTHL, ZURW, ZURC, ZURI,                   &
                              & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,             &
                              & ZMIXF, ZDTEVR, ILFS, IDBL, IML,IDDT,       &
                              & ZDTEVRF                                    )


!*           6.     Adjust up and downdraft mass flux to be consistent
!                   with precipitation efficiency relation.
!                   ---------------------------------------------------

       CALL CONVECT_PRECIP_ADJUST( ICONV, KLEV,                              &
                                 & ZPRES,ZUMF, ZUER, ZUDR, ZUPR, ZUTPR, ZURW,&
                                 & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,            &
                                 & ZPREF, ZTPR, ZMIXF, ZDTEVR,               &
                                 & ILFS, IDBL, ILCL, ICTL, IETL,             &
                                 & ZDTEVRF                                   )


!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------

       CALL CONVECT_CLOSURE2( ICONV, KLEV,                                &
                           & ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,           &
                           & ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,           &
                           & ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,              &
                           & ILCL, IDPL, IPBL, ILFS, ICTL, IML,          &
                           & ZUMF, ZUER, ZUDR, ZUTHL, ZURW,              &
                           & ZURC, ZURI, ZUPR,                           &
                           & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,              &
                           & ZTPR, ZSPR, ZDTEVR,                         &
                           & ZCAPE, ZTIMEC,                              &
                           & IFTSTEPS,                                   &
                           & ZDTEVRF, ZPRLFLX, ZPRSFLX, ZTIMC,ITSTEP,GWORK1 )



!*           8.     Determine the final grid-scale (environmental) convective
!                   tendencies and set convective counter
!                   --------------------------------------------------------


!*           8.1    Grid scale tendencies
!                   ---------------------

          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values

      DO JK = IKB, IKE
         ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
                    & * ( ZPRES(:,JK) / RATM ) ** ZRDOCP  ! change theta in temperature
         ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
                    &                            / ZTIMEC(:)

         ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
         ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:)

         ZPRLFLX(:,JK) = ZPRLFLX(:,JK) / ( RHOH2O * ZDXDY(:) )
         ZPRSFLX(:,JK) = ZPRSFLX(:,JK) / ( RHOH2O * ZDXDY(:) )
      ENDDO


      ZPRLFLX(:,IKB) = ZPRLFLX(:,IKB+1)
      ZPRSFLX(:,IKB) = ZPRSFLX(:,IKB+1)


!*           8.2    Apply conservation correction
!                   -----------------------------
      FORCE_CONSERVATION: if (CONSERV_CORRECT) then

         ! Compute moisture integral
       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = _ZERO_
       ZWORK2B(:)= _ZERO_
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
               ZW1 = _HALF_ *  (ZPRES(JI,MAX(JK-1,IKB+1)) - ZPRES(JI,JKP)) / RG
               ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) *   &
              &                   ZW1
               ZWORK2B(JI) = ZWORK2B(JI) + ZW1
         ENDDO
       ENDDO
         ! Compute moisture correction factor (integrates to surface precip.)
       DO JI = 1, ICONV
         IF ( ZTPR(JI) > _ZERO_) THEN
               ZWORK2(JI) = ( ZTPR(JI) / ZDXDY(JI) + ZWORK2(JI) ) / ZWORK2B(JI)
         ENDIF
       ENDDO
         ! Apply moisture correction
       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > _ZERO_ .AND. JK <= ICTL(JI) ) THEN
                  ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)
           ENDIF
         ENDDO
       ENDDO
         ! Compute enthalpy integral
         ZWORK2(:) = _ZERO_
         ZWORK2B(:)= _ZERO_
         DO JK = IKB+1, JKM
            JKP = JK + 1
            DO JI = 1, ICONV
               ZW1 = _HALF_ *  (ZPRES(JI,MAX(JK-1,IKB+1)) - ZPRES(JI,JKP)) / RG
               ZWORK2(JI) = ZWORK2(JI) + (                                               &
                    &  ( RCPD + RCPV * ZRW(JI,JK) )* ZTHC(JI,JK)                           &
                    &  + RCPV * ZTT(JI,JK) * (ZRVC(JI,JK)+ZRCC(JI,JK)+ZRIC(JI,JK))         &
                    &  - ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK)   &
                    &  - ( RLSTT + ( RCPV - RCS ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK)   &
                    &                        ) * ZW1
               ZWORK2B(JI) = ZWORK2B(JI) + ZW1
            ENDDO
         ENDDO
         ! Compute enthalpy correction factor (integrates to surface precip.)
         DO JI = 1, ICONV
            IF ( ZTPR(JI) > _ZERO_) THEN
               ZWORK2(JI) = ( ZTPR(JI) / ZDXDY(JI) *                                    &
                    &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,IKB) - RTT ) ) - ZWORK2(JI) )     &
                    &                           / ZWORK2B(JI)
            ENDIF
         ENDDO
         ! Apply enthalpy correction
         DO JK = JKM, IKB+1, -1
            DO JI = 1, ICONV
               IF ( ZTPR(JI) > _ZERO_ .AND. JK <= ICTL(JI) ) THEN
                  ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2(JI) / ( RCPD + RCPV * ZRW(JI,JK) )
               ENDIF
            ENDDO
         ENDDO

      endif FORCE_CONSERVATION

          ! extend tendencies to first model level

    ! DO JI = 1, ICONV
    !    ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
    !    ZTHC(JI,IKB)  = ZTHC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZTHC(JI,IKB+1)= ZTHC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    !    ZRVC(JI,IKB)  = ZRVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZRVC(JI,IKB+1)= ZRVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    ! ENDDO

              ! execute a "scatter"= pack command to store the tendencies in
              ! the final 2D tables

      DO JK = IKB, IKE
      DO JI = 1, ICONV
         JL = IJINDEX(JI)
         PTTEN(JL,JK)   = ZTHC(JI,JK)
         PRVTEN(JL,JK)  = ZRVC(JI,JK)
         PRCTEN(JL,JK)  = ZRCC(JI,JK)
         PRITEN(JL,JK)  = ZRIC(JI,JK)

         PURV(JL,JK)    = ZURW(JI,JK) - ZURC(JI,JK) - ZURI(JI,JK)
!         PURCI(JL,JK)   = ZURC(JI,JK) + ZURI(JI,JK)
         PURC(JL,JK)     = ZURC(JI,JK)
         PURI(JL,JK)     = ZURI(JI,JK)
!PV modified to output in same units as KFC i.e. kgm-2s-1
         PPRLFLX(JL,JK) = ZPRLFLX(JI,JK)*RHOH2O
         PPRSFLX(JL,JK) = ZPRSFLX(JI,JK)*RHOH2O
      ENDDO
      ENDDO




!*           8.3    Convective rainfall tendency
!                   ----------------------------

                 ! liquid and solid surface rainfall tendency in m/s
       ZTPR(:)   = ZTPR(:) / ( RHOH2O * ZDXDY(:) ) ! total surf precip
       ZSPR(:)   = ZSPR(:) / ( RHOH2O * ZDXDY(:) ) ! solid surf precip
       ZTPR(:)   = ZTPR(:) - ZSPR(:) ! compute liquid part

!PV remet en kgm-2s-1 etant donne que c'est fait dans ccdiagnostics
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        PPRLTEN(JL) = ZTPR(JI)*RHOH2O
        PPRSTEN(JL) = ZSPR(JI)*RHOH2O
     ENDDO

!                   Cloud base and top levels
!                   -------------------------
!PV change kcltop and kclbas from index to height
     ILCL(:) = MIN( ILCL(:), ICTL(:) )
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        KCLTOP(JL) = ZZ(JI,ICTL(JI))
        KCLBAS(JL) = ZZ(JI,ILCL(JI))
        IF(GWORK1(JI)) PKKFC(JL) = 1.
     ENDDO


!*           8.4    Set convective counter
!                   ----------------------

         ! compute convective counter for just activated convective
         ! grid points
         ! If the advective time period is less than specified
         ! minimum for convective period, allow feedback to occur only
         ! during advective time
     if (deep_timerefresh_sec > 0.) then
        zwork2 = min(deep_timerefresh_sec,ztimec)
     elseif (deep_timerefresh == 'TIMECONV') then
        zwork2 = ztimec
     elseif (deep_timerefresh == 'ADVECTIVE') then
        zwork2 = min(ztimea,ztimec)
     else
           call physeterror('convect_deep', 'refresh time '//trim(deep_timeconv)//' unsupported by '//trim(deep_name))
        return
     endif
     ZTIME(:) = _ONE_
     DO JI = 1, ICONV
       JL = IJINDEX(JI)
       ZTIME(JL)  =  ZWORK2(JI)
       ZWORK2(JL) =  MAX( ZTIME(JL), PDTCONV )
       IF ( GTRIG(JL) )  KCOUNT(JL) = NINT( ZWORK2(JL) / PDTCONV )
       IF ( GTRIG(JL) .AND. PPRLTEN(JL)<1.E-14 ) KCOUNT(JL) = 0
     ENDDO


!*           8.6    Compute convective tendencies for Tracers
!                   ------------------------------------------

  IF ( OCH1CONV ) THEN

       ALLOCATE( ZCH1(ICONV,IKS,KCH1) )
       ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )
       ALLOCATE( ZWORK3(ICONV,KCH1) )

       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          ZCH1(JI,JK,:) = PCH1(JL,JK,:)
       ENDDO
       ENDDO

      CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,      &
                                 & IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
                                 & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
                                 & ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB, &
                                 & IFTSTEPS )

       DO JK = IKB, IKE
       DO JN = 1, KCH1
          ZCH1C(:,JK,JN) = ( ZCH1C(:,JK,JN)- ZCH1(:,JK,JN) ) / ZTIMEC(:)
       ENDDO
       ENDDO


!*           8.7    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals

       JKM = MAXVAL( ICTL(:) )
       ZWORK3(:,:) = _ZERO_
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
           ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) *                    &
                        &   _HALF_ * (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
         ENDDO
       ENDDO

          ! Mass error (integral must be zero)

       DO JI = 1, ICONV
         IF ( ZTPR(JI) > _ZERO_) THEN
           JKP = ICTL(JI)
           ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
                        & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
           ZWORK3(JI,:) = ZWORK3(JI,:) * ZW1
         ENDIF
       ENDDO

          ! Apply uniform correction but assure positive mass at each level

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > _ZERO_ .AND. JK <= ICTL(JI) ) THEN
                ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
              ! ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
           ENDIF
         ENDDO
       ENDDO
!
          ! extend tendencies to first model level

    !  DO JI = 1, ICONV
    !     ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
    !  ENDDO
    !  DO JN = 1, KCH1
    !  DO JI = 1, ICONV
    !    ZCH1(JI,IKB,JN)  = ZCH1(JI,IKB+1,JN) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZCH1(JI,IKB+1,JN)= ZCH1(JI,IKB+1,JN) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    !  ENDDO
    !  ENDDO

       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PCH1TEN(JL,JK,:) = ZCH1C(JI,JK,:)
       ENDDO
       ENDDO
  ENDIF
!
!*           8.8    Compute convective tendencies for wind
!                   --------------------------------------

  IF ( OUVCONV ) THEN

       ALLOCATE( ZUC(ICONV,IKS) )
       ALLOCATE( ZVC(ICONV,IKS) )

       CALL CONVECT_UV_TRANSPORT_DEEP( ICONV, KLEV, ZU, ZV, ZUC, ZVC,   &
                                & IDPL, IPBL, ILCL, ICTL, ILFS, IDBL, IDDT, &
                                & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,   &
                                & ZTIMEC,  ZDXDY, ZMIXF, ZLMASS, ZWSUB,  &
                                & IFTSTEPS, ZTIMC, ITSTEP,GWORK1 )

       DO JK = IKB, IKE
          ZUC(:,JK) = ( ZUC(:,JK)- ZU(:,JK) ) / ZTIMEC(:)
          ZVC(:,JK) = ( ZVC(:,JK)- ZV(:,JK) ) / ZTIMEC(:)
       ENDDO

!*           8.9    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals

       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = _ZERO_
       ZWORK2B(:)= _ZERO_
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
           ZW1 = _HALF_ *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
           ZWORK2(JI) = ZWORK2(JI) + ZUC(JI,JK) * ZW1
           ZWORK2B(JI)= ZWORK2B(JI)+ ZVC(JI,JK) * ZW1
         ENDDO
       ENDDO

          !  error (integral must be zero)

       DO JI = 1, ICONV
         IF ( ZTPR(JI) > _ZERO_) THEN
           JKP = ICTL(JI)
           ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
                       & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
           ZWORK2(JI) = ZWORK2(JI) * ZW1
           ZWORK2B(JI)= ZWORK2B(JI)* ZW1
         ENDIF
       ENDDO

          ! Apply uniform correction

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > _ZERO_ .AND. JK <= ICTL(JI) ) THEN
                ZUC(JI,JK) = ZUC(JI,JK) - ZWORK2(JI)
                ZVC(JI,JK) = ZVC(JI,JK) - ZWORK2B(JI)
           ENDIF
         ENDDO
       ENDDO
!
          ! extend tendencies to first model level

    ! DO JI = 1, ICONV
    !    ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
    !    ZUC(JI,IKB)  = ZUC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZUC(JI,IKB+1)= ZUC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    !    ZVC(JI,IKB)  = ZVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZVC(JI,IKB+1)= ZVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    ! ENDDO


       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PUTEN(JL,JK)   = ZUC(JI,JK)
          PVTEN(JL,JK)   = ZVC(JI,JK)
       ENDDO
       ENDDO
  ENDIF

!
!*           9.     Write up- and downdraft mass fluxes
!                   ------------------------------------
!
    DO JK = IKB, IKE
       ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
       ZDMF(:,JK)  = ZDMF(:,JK) / ZDXDY(:)
    ENDDO
    DO JK = IKB+1, IKE
       ZUDR(:,JK)  = ZUDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )! detrainment for ERA40
       ZDDR(:,JK)  = ZDDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )
    ENDDO
    ZWORK2(:) = _ONE_
    DO JK = IKB, IKE
    DO JI = 1, ICONV
       JL = IJINDEX(JI)
       IF ( PPRLTEN(JL) < 1.E-14 ) ZWORK2(JL) = _ZERO_
       PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)
       PDMF(JL,JK) = ZDMF(JI,JK) * ZWORK2(JL)
       PUDR(JL,JK) = ZUDR(JI,JK) * ZWORK2(JL)
       PDDR(JL,JK) = ZDDR(JI,JK) * ZWORK2(JL)
       PURV(JL,JK) = PURV(JL,JK) * ZWORK2(JL)
!       PURCI(JL,JK)= PURCI(JL,JK)* ZWORK2(JL)
!PV calculate cloud fraction -same recipe as in  KFC
       IF ( ZWU(JI,JK) > 1.E-1 ) THEN
         PCF(JL,JK) = ZUMF(JI,JK)/(ZRHO(JI,JK)*ZWU(JI,JK)) * ZWORK2(JL)
       ELSE
         PCF(JL,JK) = 0.0
       ENDIF
! calculate updraft areas from cloud fractions
       PAREAUP(JL,JK)=PCF(JL,JK)*ZDXDY(JI)
       PWUMAXOUT(JL) = MAX( ZWU(JI,JK), PWUMAXOUT(JL) )

!PV output "grid-scale" condensate rather than "in-cloud"  condensate
       PURC(JL,JK)= PURC(JL,JK)* PCF(JL,JK) *ZWORK2(JL)
       PURI(JL,JK)= PURI(JL,JK)* PCF(JL,JK) *ZWORK2(JL)
    ENDDO
    ENDDO
!
!
!*           10.    Deallocate all local arrays
!                   ---------------------------
!
! downdraft variables
!
      DEALLOCATE( ZDMF )
      DEALLOCATE( ZDER )
      DEALLOCATE( ZDDR )
      DEALLOCATE( ZDTHL )
      DEALLOCATE( ZDRW )
      DEALLOCATE( ZLMASS )
      DEALLOCATE( ZMIXF )
      DEALLOCATE( ZTPR )
      DEALLOCATE( ZSPR )
      DEALLOCATE( ZDTEVR )
      DEALLOCATE( ZPREF )
      DEALLOCATE( IML )
      DEALLOCATE( ILFS )
      DEALLOCATE( IDBL )
      DEALLOCATE( IDDT )
      DEALLOCATE( ZDTEVRF )
      DEALLOCATE( ZPRLFLX )
      DEALLOCATE( ZPRSFLX )
!
!   closure variables
!
      DEALLOCATE( ZTIMEA )
      DEALLOCATE( ZTIMEC )
      DEALLOCATE( ZTHC )
      DEALLOCATE( ZRVC )
      DEALLOCATE( ZRCC )
      DEALLOCATE( ZRIC )
      DEALLOCATE( ZWSUB )
!
       IF ( OCH1CONV ) THEN
           DEALLOCATE( ZCH1 )
           DEALLOCATE( ZCH1C )
           DEALLOCATE( ZWORK3 )
       ENDIF
       IF ( OUVCONV ) THEN
           DEALLOCATE( ZUC )
           DEALLOCATE( ZVC )
       ENDIF
!
    ENDIF
!
!    vertical index
!
    DEALLOCATE( IDPL )
    DEALLOCATE( IPBL )
    DEALLOCATE( ILCL )
    DEALLOCATE( ICTL )
    DEALLOCATE( IETL )
!
! grid scale variables
!
    DEALLOCATE( ZZ )
    DEALLOCATE( ZPRES )
    DEALLOCATE( ZDPRES )
    DEALLOCATE( ZU )
    DEALLOCATE( ZV )
    DEALLOCATE( ZRHO )
    DEALLOCATE( ZTT )
    DEALLOCATE( ZTH )
    DEALLOCATE( ZTHV )
    DEALLOCATE( ZTHL )
    DEALLOCATE( ZTHES )
    DEALLOCATE( ZRW )
    DEALLOCATE( ZRV )
    DEALLOCATE( ZRC )
    DEALLOCATE( ZRI )
    DEALLOCATE( ZDXDY )
!
! updraft variables
!
    DEALLOCATE( ZUMF )
    DEALLOCATE( ZUER )
    DEALLOCATE( ZUDR )
    DEALLOCATE( ZUTHL )
    DEALLOCATE( ZUTHV )
    DEALLOCATE( ZWU )
    DEALLOCATE( ZURW )
    DEALLOCATE( ZURC )
    DEALLOCATE( ZURI )
    DEALLOCATE( ZURR )
    DEALLOCATE( ZURS )
    DEALLOCATE( ZUPR )
    DEALLOCATE( ZUTPR )
    DEALLOCATE( ZTHLCL )
    DEALLOCATE( ZTLCL )
    DEALLOCATE( ZRVLCL )
    DEALLOCATE( ZWLCL )
    DEALLOCATE( ZZLCL )
    DEALLOCATE( ZTHVELCL )
    DEALLOCATE( ZMFLCL )
    DEALLOCATE( ZCAPE )
!
! work arrays
!
    DEALLOCATE( IINDEX )
    DEALLOCATE( IJINDEX )
    DEALLOCATE( IJSINDEX )
    DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE CONVECT_DEEP2
