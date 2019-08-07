!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine CONVECT_DEEP2(KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,          &
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

   !!**** Monitor routine to compute all convective tendencies by calls
   !!     of several subroutines.
   !
   !
   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine the convective
   !!      tendencies. The routine first prepares all necessary grid-scale
   !!      variables. The final convective tendencies are then computed by
   !!      calls of different subroutines.
   !
   !
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
   !
   !
   !
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
   !
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
   !
   !!      Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
   !
   !!      Module YOE_CONVPAR
   !!          XA25                 ! reference grid area
   !!          XCRAD                ! cloud radius
   !
   !
   !
   !!    REFERENCE
   !!    ---------
   !
   !!      Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc. :
   !!           A mass flux convection scheme for regional and global models.
   !!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
   !!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
   !
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    26/03/96
   !!   Peter Bechtold 04/10/97 replace theta_il by enthalpy
   !!         "        10/12/98 changes for ARPEGE
   !----------------------------------------------------------------------------


   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAREXT
   use YOE_CONVPAR
   use cnv_options, deep_name=>deep

   implicit none
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

   !*       0.1   Declarations of dummy arguments :


   integer,                    intent(IN) :: KLON     ! horizontal dimension
   integer,                    intent(IN) :: KLEV     ! vertical dimension
   integer,                    intent(IN) :: KIDIA    ! value of the first point in x
   integer,                    intent(IN) :: KFDIA    ! value of the last point in x
   integer,                    intent(IN) :: KBDIA    ! vertical  computations start at
   !                                                    ! KBDIA that is at least 1
   integer,                    intent(IN) :: KTDIA    ! vertical computations can be
   ! limited to KLEV + 1 - KTDIA
   ! default=1
   real,                       intent(IN) :: PDTCONV  ! Interval of time between two
   ! calls of the deep convection
   ! scheme
   integer,                    intent(IN) :: KICE     ! flag for ice ( 1 = yes,
   !                0 = no ice )
   logical,                      intent(IN) :: ODOWN    ! take or not convective
   ! downdrafts into account
   real, dimension(KLON,KLEV), intent(IN) :: PTT      ! grid scale temperature (K)
   real, dimension(KLON,KLEV), intent(IN) :: PRVT     ! grid scale water vapor(kg/kg)
   real, dimension(KLON,KLEV), intent(IN) :: PRCT     ! grid scale r_c (kg/kg)
   real, dimension(KLON,KLEV), intent(IN) :: PRIT     ! grid scale r_i (kg/kg)
   real, dimension(KLON,KLEV), intent(IN) :: PUT      ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV), intent(IN) :: PVT      ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV), intent(IN) :: PWT      ! grid scale vertical
   ! velocity (m/s)
   real, dimension(KLON,KLEV), intent(IN) :: PPABST   ! grid scale pressure at (Pa)
   real, dimension(KLON,KLEV), intent(IN) :: PZZ      ! height of model layer (m)
   real, dimension(KLON),      intent(IN) :: PDXDY    ! horizontal grid area (m2)
   real, dimension(KLON,KLEV), intent(IN) :: PHSFLX   ! turbulent heat flux (W/m2)

   integer, dimension(KLON),   intent(INOUT):: KCOUNT ! convective counter (recompute
   ! tendency or keep it)
   real, dimension(KLON,KLEV), intent(INOUT):: PTTEN  ! convective temperature
   ! tendency (K/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRVTEN ! convective r_v tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRCTEN ! convective r_c tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRITEN ! convective r_i tendency (1/s)
   real, dimension(KLON),      intent(INOUT):: PPRLTEN! liquid surf. precipitation
   ! tendency (kg/s m2)
   real, dimension(KLON),      intent(INOUT):: PPRSTEN! solid surf. precipitation
   ! tendency (kg/s m2)
   real, dimension(KLON),   intent(INOUT):: KCLTOP ! cloud top level (m)
   real, dimension(KLON),   intent(INOUT):: KCLBAS ! cloud base level (m)
   ! they are given a value of
   ! 0 if no convection
   real, dimension(KLON,KLEV), intent(INOUT):: PPRLFLX! liquid precip flux (m/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PPRSFLX! solid  precip flux (m/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
   real, dimension(KLON,KLEV), intent(INOUT):: PDMF   ! downdraft mass flux (kg/s m2)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURV   ! updraft water vapor (kg/kg)
   !real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURCI  ! updraft liquid+ice condensate (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURC   ! updraft liquid condensate (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURI   ! updraft ice condensate (kg/kg)
   real, dimension(KLON),      intent(INOUT):: PCAPE  ! maximum CAPE (J/kg)
   real, dimension(KLON,KLEV), intent(INOUT)  :: PCF    ! Cloud Fraction
   real, dimension(KLON,KLEV), intent(INOUT)  :: PAREAUP! area of updraft, cloud coverage area (m^2)= PCF x dxdy
   real, dimension(KLON),      intent(INOUT)  :: PWUMAXOUT ! maximum upward velocity in bechtold deep convection

   logical,                      intent(IN)   :: OUVCONV! include wind transport (Cu friction)
   real, dimension(KLON,KLEV), intent(INOUT):: PUTEN  ! convective u tendency (m/s^2)
   real, dimension(KLON,KLEV), intent(INOUT):: PVTEN  ! convective v tendency (m/s^2)

   logical,                      intent(IN) :: OCH1CONV ! include tracer transport
   integer,                    intent(IN) :: KCH1     ! number of species
   real, dimension(KLON,KLEV,KCH1), intent(IN)   :: PCH1    ! grid scale chemical species
   real, dimension(KLON,KLEV,KCH1), intent(INOUT):: PCH1TEN ! species conv. tendency (1/s)

   ! for ERA40
   real, dimension(KLON,KLEV), intent(OUT)  :: PUDR   ! updraft detrainment rate (kg/s m2)
   real, dimension(KLON,KLEV), intent(OUT)  :: PDDR   ! downdraft detrainment rate (kg/s m2)
   real, dimension(KLON),      intent(OUT)  :: PKKFC  ! shallow convective counter


   real, dimension(KLON),      intent(IN)  :: XLAT    ! Latitude
   real, dimension(KLON),      intent(IN)  :: MG      ! land-sea mask
   real, dimension(KLON),      intent(IN)  :: MLAC    ! lake mask


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
   logical, parameter :: CONSERV_CORRECT=.true. !apply conservation correction

   logical, dimension(KLON, KLEV)        :: GTRIG3 ! 3D logical mask for convection
   logical, dimension(KLON)              :: GTRIG  ! 2D logical mask for trigger test
   real,  dimension(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta,
   ! theta_v, theta_es
   real,  dimension(KLON)              :: ZTIME  ! convective time period
   real,  dimension(KLON)              :: ZTIMC   ! fractional convective time step
   integer,  dimension(KLON)           :: ITSTEP  ! # of fractional convective timesteps
   logical,  dimension(KLON)             :: GWORK1  ! flag for newly activated columns including MASS CONSERVATION criteria from closure
   real,  dimension(KLON)              :: ZWORK2, ZWORK2B ! work array
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
   !*       0.2   Declarations of local allocatable  variables :

   integer, dimension(:),pointer  :: IDPL    ! index for parcel departure level
   integer, dimension(:),pointer  :: IPBL    ! index for source layer top
   integer, dimension(:),pointer  :: ILCL    ! index for lifting condensation level
   integer, dimension(:),pointer  :: IETL    ! index for zero buoyancy level
   integer, dimension(:),pointer  :: ICTL    ! index for cloud top level
   integer, dimension(:),pointer  :: ILFS    ! index for level of free sink
   integer, dimension(:),pointer  :: IDBL    ! index for downdraft base level
   integer, dimension(:),pointer  :: IDDT    ! index for downdraft detrainment top level
   integer, dimension(:),pointer  :: IML     ! melting level

   integer, dimension(:), pointer :: ISDPL   ! index for parcel departure level
   integer, dimension(:), pointer :: ISPBL   ! index for source layer top
   integer, dimension(:), pointer :: ISLCL   ! index for lifting condensation level

   real, dimension(:), pointer    :: ZSTHLCL ! updraft theta at LCL
   real, dimension(:), pointer    :: ZSTLCL  ! updraft temp. at LCL
   real, dimension(:), pointer    :: ZSRVLCL ! updraft rv at LCL
   real, dimension(:), pointer    :: ZSWLCL  ! updraft w at LCL
   real, dimension(:), pointer    :: ZSZLCL  ! LCL height
   real, dimension(:), pointer    :: ZSTHVELCL! envir. theta_v at LCL
   real, dimension(:), pointer    :: ZSDXDY  ! grid area (m^2)

   ! grid scale variables
   real, dimension(:,:), pointer  :: ZZ      ! height of model layer (m)
   real, dimension(:,:), pointer  :: ZPRES   ! grid scale pressure
   real, dimension(:,:), pointer  :: ZDPRES  ! pressure difference between
   ! bottom and top of layer (Pa)
   real, dimension(:,:), pointer  :: ZU      ! grid scale horiz. u component on theta grid
   real, dimension(:,:), pointer  :: ZV      ! grid scale horiz. v component on theta grid
   real, dimension(:,:), pointer  :: ZW      ! grid scale vertical velocity on theta grid
   real, dimension(:,:), pointer  :: ZHSFLX  ! turbulent sensible heat flux (W/m2)
   real, dimension(:,:), pointer  :: ZRHO    ! air density
   real, dimension(:,:), pointer  :: ZTT     ! temperature
   real, dimension(:,:), pointer  :: ZTH     ! grid scale theta
   real, dimension(:,:), pointer  :: ZTHV    ! grid scale theta_v
   real, dimension(:,:), pointer  :: ZTHL    ! grid scale enthalpy (J/kg)
   real, dimension(:,:), pointer  :: ZTHES, ZTHEST ! grid scale saturated theta_e
   real, dimension(:,:), pointer  :: ZRW     ! grid scale total water (kg/kg)
   real, dimension(:,:), pointer  :: ZRV     ! grid scale water vapor (kg/kg)
   real, dimension(:,:), pointer  :: ZRC     ! grid scale cloud water (kg/kg)
   real, dimension(:,:), pointer  :: ZRI     ! grid scale cloud ice (kg/kg)
   real, dimension(:),   pointer  :: ZDXDY   ! grid area (m^2)

   ! updraft variables
   real, dimension(:,:), pointer  :: ZUMF    ! updraft mass flux (kg/s)
   real, dimension(:,:), pointer  :: ZUER    ! updraft entrainment (kg/s)
   real, dimension(:,:), pointer  :: ZUDR    ! updraft detrainment (kg/s)
   real, dimension(:,:), pointer  :: ZUPR    ! updraft precipitation in
   ! flux units (kg water / s)
   real, dimension(:,:), pointer  :: ZUTHL   ! updraft enthalpy (J/kg)
   real, dimension(:,:), pointer  :: ZUTHV   ! updraft theta_v (K)
   real, dimension(:,:), pointer  :: ZWU     ! updraft vertical vel (m/s)
   real, dimension(:,:), pointer  :: ZURW    ! updraft total water (kg/kg)
   real, dimension(:,:), pointer  :: ZURC    ! updraft cloud water (kg/kg)
   real, dimension(:,:), pointer  :: ZURI    ! updraft cloud ice   (kg/kg)
   real, dimension(:,:), pointer  :: ZURR    ! liquid precipit. (kg/kg)
   ! produced in  model layer
   real, dimension(:,:), pointer  :: ZURS    ! solid precipit. (kg/kg)
   ! produced in  model layer
   real, dimension(:),   pointer  :: ZUTPR   ! total updraft precipitation (kg/s)
   real, dimension(:),   pointer  :: ZMFLCL  ! cloud base unit mass flux(kg/s)
   real, dimension(:),   pointer  :: ZCAPE   ! available potent. energy
   real, dimension(:),   pointer  :: ZTHLCL  ! updraft theta at LCL
   real, dimension(:),   pointer  :: ZTLCL   ! updraft temp. at LCL
   real, dimension(:),   pointer  :: ZRVLCL  ! updraft rv at LCL
   real, dimension(:),   pointer  :: ZWLCL   ! updraft w at LCL
   real, dimension(:),   pointer  :: ZZLCL   ! LCL height
   real, dimension(:),   pointer  :: ZTHVELCL! envir. theta_v at LCL

   ! downdraft variables
   real, dimension(:,:), pointer  :: ZDMF    ! downdraft mass flux (kg/s)
   real, dimension(:,:), pointer  :: ZDER    ! downdraft entrainment (kg/s)
   real, dimension(:,:), pointer  :: ZDDR    ! downdraft detrainment (kg/s)
   real, dimension(:,:), pointer  :: ZDTHL   ! downdraft enthalpy (J/kg)
   real, dimension(:,:), pointer  :: ZDRW    ! downdraft total water (kg/kg)
   real, dimension(:),   pointer  :: ZMIXF   ! mixed fraction at LFS
   real, dimension(:),   pointer  :: ZTPR    ! total surf precipitation (kg/s)
   real, dimension(:),   pointer  :: ZSPR    ! solid surf precipitation (kg/s)
   real, dimension(:),   pointer  :: ZDTEVR  ! donwndraft evapor. (kg/s)
   real, dimension(:),   pointer  :: ZPREF   ! precipitation efficiency
   real, dimension(:,:), pointer  :: ZDTEVRF ! donwndraft evapor. (kg/s)
   real, dimension(:,:), pointer  :: ZPRLFLX ! liquid precip flux
   real, dimension(:,:), pointer  :: ZPRSFLX ! solid precip flux

   ! closure variables
   real, dimension(:,:), pointer  :: ZLMASS  ! mass of model layer (kg)
   real, dimension(:),   pointer  :: ZTIMEA  ! advective time period
   real, dimension(:),   pointer  :: ZTIMEC  ! time during which convection is
   ! active at grid point (as ZTIME)

   real, dimension(:,:), pointer  :: ZTHC    ! conv. adj. grid scale theta
   real, dimension(:,:), pointer  :: ZRVC    ! conv. adj. grid scale r_w
   real, dimension(:,:), pointer  :: ZRCC    ! conv. adj. grid scale r_c
   real, dimension(:,:), pointer  :: ZRIC    ! conv. adj. grid scale r_i
   real, dimension(:,:), pointer  :: ZWSUB   ! envir. compensating subsidence (Pa/s)

   logical,   dimension(:),pointer  :: GTRIG1  ! logical mask for convection
   logical,   dimension(:),pointer  :: GWORK   ! logical work array
   integer, dimension(:),pointer  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index

   real, dimension(:),   pointer  :: ZCPH    ! specific heat C_ph
   real, dimension(:),   pointer  :: ZLV, ZLS! latent heat of vaporis., sublim.
   real                               :: ZES     ! saturation vapor mixng ratio

   ! for U, V transport:
   real,  dimension(:,:), pointer :: ZUC     ! horizontal wind u (m/s)
   real,  dimension(:,:), pointer :: ZVC     ! horizontal wind v (m/s)

   ! for Chemical Tracer transport:
   real,  dimension(:,:,:), pointer:: ZCH1    ! grid scale chemical specy (kg/kg)
   real,  dimension(:,:,:), pointer:: ZCH1C   ! conv. adjust. chemical specy 1
   real,  dimension(:,:),   pointer:: ZWORK3  ! work arrays
   logical, dimension(:,:,:), pointer:: GTRIG4  ! logical mask

   !-------------------------------------------------------------------------------
   RHOH2O = RATM *1.E-2

   !*       0.3    Compute loop bounds
   !               -------------------

   IIB    = KIDIA
   IIE    = KFDIA
   JCVEXB = max( 0, KBDIA - 1 )
   IKB    = 1 + JCVEXB
   IKS    = KLEV
   JCVEXT = max( 0, KTDIA - 1 )
   IKE    = IKS - JCVEXT


   !*       0.5    Update convective counter ( where KCOUNT > 0
   !               convection is still active ).
   !               ---------------------------------------------

   KCOUNT(IIB:IIE) = max(KCOUNT(IIB:IIE) - 1,0)

   GTRIG(:)  = KCOUNT(:) <= 0
   ITEST     = count( GTRIG(:) )

   where ( KCOUNT(:) > 0 )
      PKKFC(:) = 1.
   end where


   !PV - must undo change in units done in ccdiagnostics.ftn90 to all precip rates
   !PV - or else precip rates will divided by 1000 at each timestep during convective timesacle
   PPRLTEN(:) = PPRLTEN(:)*1000.0
   PPRSTEN(:) = PPRSTEN(:)*1000.0

   if ( ITEST == 0 ) return  ! if convection is already active at every grid point
   ! exit CONVECT_DEEP



   !*       0.7    Reset convective tendencies to zero if convective
   !               counter becomes negative
   !               -------------------------------------------------

   GTRIG3(:,:) = spread( GTRIG(:), DIM=2, NCOPIES=IKS )
   where ( GTRIG3(:,:) )
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
   end where

   where ( GTRIG(:) )
      PPRLTEN(:) = _ZERO_
      PPRSTEN(:) = _ZERO_
      KCLTOP(:)  = _ZERO_
      KCLBAS(:)  = _ZERO_
      PCAPE(:)   = _ZERO_
      PWUMAXOUT(:)=_ZERO_
   end where

   if ( OCH1CONV ) then
      allocate( GTRIG4(KLON,KLEV,KCH1) )
      GTRIG4(:,:,:) = spread( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
      where( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = _ZERO_
      deallocate( GTRIG4 )
   endif


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

   do JK = IKB, IKE
      do JI = IIB, IIE
         if ( PPABST(JI,JK) > 40.E2 ) then
            ZTHT(JI,JK)  = PTT(JI,JK) * ( RATM / PPABST(JI,JK) ) ** ZRDOCP
            ZSTHV(JI,JK) = ZTHT(JI,JK) * ( _ONE_ + ZEPSA * PRVT(JI,JK) ) /    &
                 & ( _ONE_ + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )

            ! use conservative Bolton (1980) formula for theta_e
            ! it is used to compute CAPE for undilute parcel ascent
            ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here

            ZES = exp( RALPW - RBETW / PTT(JI,JK) - RGAMW * log( PTT(JI,JK) ) )
            ZES = min( _ONE_, ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
            ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **      &
                 & ( _ONE_ - 0.28 * ZES ) *                        &
                 &   exp( ( 3374.6525 / PTT(JI,JK) - 2.5403 ) &
                 &   * ZES * ( _ONE_ + 0.81 * ZES ) )
         endif
      enddo
   enddo



   !*       2.     Test for convective columns and determine properties at the LCL
   !               --------------------------------------------------------------

   !*       2.1    Allocate arrays depending on number of model columns that need
   !               to be tested for convection (i.e. where no convection is present
   !               at the moment.
   !               --------------------------------------------------------------

   allocate( ZPRES(ITEST,IKS) )
   allocate( ZZ(ITEST,IKS) )
   allocate( ZW(ITEST,IKS) )
   allocate( ZTH(ITEST,IKS) )
   allocate( ZTHV(ITEST,IKS) )
   allocate( ZTHEST(ITEST,IKS) )
   allocate( ZRV(ITEST,IKS) )
   allocate( ZSTHLCL(ITEST) )
   allocate( ZSTLCL(ITEST) )
   allocate( ZSRVLCL(ITEST) )
   allocate( ZSWLCL(ITEST) )
   allocate( ZSZLCL(ITEST) )
   allocate( ZSTHVELCL(ITEST) )
   allocate( ISDPL(ITEST) )
   allocate( ISPBL(ITEST) )
   allocate( ISLCL(ITEST) )
   allocate( ZSDXDY(ITEST) )
   allocate( GTRIG1(ITEST) )
   allocate( ZCAPE(ITEST) )
   allocate( IINDEX(KLON) )
   allocate( IJSINDEX(ITEST) )
   allocate( ZHSFLX(ITEST,IKS) )
   do JI = 1, KLON
      IINDEX(JI) = JI
   enddo
   IJSINDEX(:) = pack( IINDEX(:), MASK=GTRIG(:) )

   do JK = IKB, IKE
      do JI = 1, ITEST
         JL = IJSINDEX(JI)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
         ZTHEST(JI,JK) = ZSTHES(JL,JK)
         ZRV(JI,JK)    = max( _ZERO_, PRVT(JL,JK) )
         ZW(JI,JK)     = PWT(JL,JK)
         ZHSFLX(JI,JK) = PHSFLX(JL,JK)
      enddo
   enddo
   do JI = 1, ITEST
      JL = IJSINDEX(JI)
      ZSDXDY(JI)    = PDXDY(JL)
   enddo

   !*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c
   !               and envir. saturation theta_e
   !               ------------------------------------------------------------


   !*       2.3    Test for convective columns and determine properties at the LCL
   !               ------------------------------------------------------------

   ISLCL(:) = max( IKB, 2 )   ! initialize DPL PBL and LCL
   ISDPL(:) = IKB
   ISPBL(:) = IKB


   call CONVECT_TRIGGER_FUNCT3(ITEST, KLEV,         &
        & ZPRES, ZTH, ZTHV, ZTHEST,                 &
        & ZRV, ZW, ZZ,                              &
        & ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
        & ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1,   &
        & ZCAPE, XLAT, MG, MLAC )

   !PV cape in columns rejected by criteria in trigger must be put to zero
   where ( .not. GTRIG1(:) )
      ZCAPE(:)   = _ZERO_
   end where

   do JI = 1, ITEST
      JL = IJSINDEX(JI)
      PCAPE(JL) = ZCAPE(JI)
   enddo

   deallocate( ZPRES )
   deallocate( ZZ )
   deallocate( ZTH )
   deallocate( ZTHV )
   deallocate( ZTHEST )
   deallocate( ZRV )
   deallocate( ZW )
   deallocate( ZCAPE )
   deallocate( ZHSFLX )


   !*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
   !               arrays used in the convection scheme using the mask GTRIG, i.e.
   !               we do calculus only in convective columns. This corresponds to
   !               a GATHER operation.
   !               ------------------------------------------------------------

   ICONV = count( GTRIG1(:) )
   if ( ICONV == 0 )  then
      deallocate( ZSTHLCL )
      deallocate( ZSTLCL )
      deallocate( ZSRVLCL )
      deallocate( ZSWLCL )
      deallocate( ZSZLCL )
      deallocate( ZSTHVELCL )
      deallocate( ZSDXDY )
      deallocate( ISLCL )
      deallocate( ISDPL )
      deallocate( ISPBL )
      deallocate( GTRIG1 )
      deallocate( IINDEX )
      deallocate( IJSINDEX )
      return   ! no convective column has been found, exit CONVECT_DEEP
   endif

   ! vertical index variables

   allocate( IDPL(ICONV) )
   allocate( IPBL(ICONV) )
   allocate( ILCL(ICONV) )
   allocate( ICTL(ICONV) )
   allocate( IETL(ICONV) )

   ! grid scale variables

   allocate( ZZ(ICONV,IKS) )
   allocate( ZPRES(ICONV,IKS) )
   allocate( ZDPRES(ICONV,IKS) )
   allocate( ZU(ICONV,IKS) )
   allocate( ZV(ICONV,IKS) )
   allocate( ZRHO(ICONV, IKS) )
   allocate( ZTT(ICONV, IKS) )
   allocate( ZTH(ICONV,IKS) )
   allocate( ZTHV(ICONV,IKS) )
   allocate( ZTHL(ICONV,IKS) )
   allocate( ZTHES(ICONV,IKS) )
   allocate( ZRV(ICONV,IKS) )
   allocate( ZRC(ICONV,IKS) )
   allocate( ZRI(ICONV,IKS) )
   allocate( ZRW(ICONV,IKS) )
   allocate( ZDXDY(ICONV) )

   ! updraft variables

   allocate( ZUMF(ICONV,IKS) )
   allocate( ZUER(ICONV,IKS) )
   allocate( ZUDR(ICONV,IKS) )
   allocate( ZUPR(ICONV,IKS) )
   allocate( ZUTHL(ICONV,IKS) )
   allocate( ZUTHV(ICONV,IKS) )
   allocate( ZWU(ICONV,IKS) )
   allocate( ZURW(ICONV,IKS) )
   allocate( ZURC(ICONV,IKS) )
   allocate( ZURI(ICONV,IKS) )
   allocate( ZURR(ICONV,IKS) )
   allocate( ZURS(ICONV,IKS) )
   allocate( ZUTPR(ICONV) )
   allocate( ZTHLCL(ICONV) )
   allocate( ZTLCL(ICONV) )
   allocate( ZRVLCL(ICONV) )
   allocate( ZWLCL(ICONV) )
   allocate( ZMFLCL(ICONV) )
   allocate( ZZLCL(ICONV) )
   allocate( ZTHVELCL(ICONV) )
   allocate( ZCAPE(ICONV) )

   ! work variables

   allocate( IJINDEX(ICONV) )
   allocate( IJPINDEX(ICONV) )
   allocate( ZCPH(ICONV) )
   allocate( ZLV(ICONV) )
   allocate( ZLS(ICONV) )


   !*           3.1    Gather grid scale and updraft base variables in
   !                   arrays using mask GTRIG
   !                   ---------------------------------------------------

   GTRIG(:)      = unpack( GTRIG1(:), MASK=GTRIG(:), FIELD=.false. )
   IJINDEX(:)    = pack( IINDEX(:), MASK=GTRIG(:) )

   do JK = IKB, IKE
      do JI = 1, ICONV
         JL = IJINDEX(JI)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZTT(JI,JK)    = PTT(JL,JK)
         ZRHO(JI,JK)   = ZPRES(JI,JK)/(ZTT(JI,JK)*RD)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHES(JI,JK)  = ZSTHES(JL,JK)
         ZRV(JI,JK)    = max( _ZERO_, PRVT(JL,JK) )
         ZRC(JI,JK)    = max( _ZERO_, PRCT(JL,JK) )
         ZRI(JI,JK)    = max( _ZERO_, PRIT(JL,JK) )
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
         ZU(JI,JK)     = PUT(JL,JK)
         ZV(JI,JK)     = PVT(JL,JK)
      enddo
   enddo

   do JI = 1, ITEST
      IJSINDEX(JI) = JI
   enddo
   IJPINDEX(:) = pack( IJSINDEX(:), MASK=GTRIG1(:) )
   do JI = 1, ICONV
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
   enddo
   allocate( GWORK(ICONV) )
   GWORK(:)      = pack( GTRIG1(:),  MASK=GTRIG1(:) )
   deallocate( GTRIG1 )
   allocate( GTRIG1(ICONV) )
   GTRIG1(:)     = GWORK(:)

   deallocate( GWORK )
   !         DEALLOCATE( IJPINDEX )
   deallocate( ISDPL )
   deallocate( ISPBL )
   deallocate( ISLCL )
   deallocate( ZSTHLCL )
   deallocate( ZSTLCL )
   deallocate( ZSRVLCL )
   deallocate( ZSWLCL )
   deallocate( ZSZLCL )
   deallocate( ZSTHVELCL )
   deallocate( ZSDXDY )


   !*           3.2    Compute pressure difference
   !                   ---------------------------------------------------

   ZDPRES(:,IKB) = _ZERO_
   do JK = IKB + 1, IKE
      ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK)
   enddo

   !*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c
   !                  ----------------------------------------------------------

   do JK = IKB, IKE, 1
      ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
      ZCPH(:)    = RCPD + RCPV * ZRW(:,JK)
      ZLV(:)     = RLVTT + ( RCPV - RCW ) * ( ZTT(:,JK) - RTT ) ! compute L_v
      ZLS(:)     = RLSTT + ( RCPV - RCS ) * ( ZTT(:,JK) - RTT ) ! compute L_i
      ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( _ONE_ + ZRW(:,JK) ) * RG * ZZ(:,JK) &
           & - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
   enddo


   !*           4.     Compute updraft properties
   !                   ----------------------------

   !*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s )
   !                   ---------------------------------------------------------

   do JI = 1, ICONV
      JK = ILCL(JI) - 1
      ZMFLCL(JI) = ZPRES(JI,JK) / ( RD * ZTT(JI,JK) *                    &
           & ( _ONE_ + ZEPS * ZRVLCL(JI) ) ) * RPI * XCRAD * XCRAD &
           & * max( _ONE_, ZDXDY(JI)/XA25 )
   enddo

   deallocate( ZCPH )
   deallocate( ZLV )
   deallocate( ZLS )


   call CONVECT_UPDRAFT( ICONV, KLEV,                                     &
        & KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
        & ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   &
        & ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
        & ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZWU,             &
        & ZURW, ZURC, ZURI, ZURR, ZURS, ZUPR,              &
        & ZUTPR, ZCAPE, ICTL, IETL                         )

   ! print '(a50,5i3,x,l8)','UPDRAFT LEVELS:lc,kpbl,klcl,let,ltop:',idpl,ipbl,ilcl,ietl,ictl,gtrig1






   !*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud
   !                   thickness is smaller than 3 km
   !                   --------------------------------------------------------


   ICONV1 = count(GTRIG1)

   !PV CAPE in columns rejected by criteria in updraft must be put to zero
   do JI = 1, ICONV
      JL = IJPINDEX(JI)
      if(.not. GTRIG1(JI) ) PCAPE(JL)   = _ZERO_
   enddo
   deallocate( IJPINDEX )


   if ( ICONV1 > 0 )  then


      !*       4.3    Allocate memory for downdraft variables
      !               ---------------------------------------

      ! downdraft variables

      allocate( ILFS(ICONV) )
      allocate( IDBL(ICONV) )
      allocate( IDDT(ICONV) )
      allocate( IML(ICONV) )
      allocate( ZDMF(ICONV,IKS) )
      allocate( ZDER(ICONV,IKS) )
      allocate( ZDDR(ICONV,IKS) )
      allocate( ZDTHL(ICONV,IKS) )
      allocate( ZDRW(ICONV,IKS) )
      allocate( ZLMASS(ICONV,IKS) )
      do JK = IKB, IKE
         ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / RG  ! mass of model layer
      enddo
      ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
      allocate( ZMIXF(ICONV) )
      allocate( ZTPR(ICONV) )
      allocate( ZSPR(ICONV) )
      allocate( ZDTEVR(ICONV) )
      allocate( ZPREF(ICONV) )
      allocate( ZDTEVRF(ICONV,IKS) )
      allocate( ZPRLFLX(ICONV,IKS) )
      allocate( ZPRSFLX(ICONV,IKS) )

      ! closure variables

      allocate( ZTIMEA(ICONV) )
      allocate( ZTIMEC(ICONV) )
      allocate( ZTHC(ICONV,IKS) )
      allocate( ZRVC(ICONV,IKS) )
      allocate( ZRCC(ICONV,IKS) )
      allocate( ZRIC(ICONV,IKS) )
      allocate( ZWSUB(ICONV,IKS) )


      !*           5.     Compute downdraft properties
      !                   ----------------------------

      !*           5.1    Compute advective time period and precipitation
      !                   efficiency as a function of mean ambient wind (shear)
      !                   -----------------------------------------------------

      call CONVECT_TSTEP_PREF( ICONV, KLEV,                          &
           & ZU, ZV, ZPRES, ZZ, ZDXDY, ILCL, ICTL, &
           & ZTIMEA, ZPREF )

      ! exclude convective downdrafts if desired
      if ( .not. ODOWN ) ZPREF(:) = _ONE_

      ! Compute the period during which convection is active
      if (deep_timeconv_sec > 0.) then
         ztimec(:) = max(pdtconv,deep_timeconv_sec)
      elseif (deep_timeconv == 'ADVECTIVE') then
         ZTIMEC(:) = max( 1800., min( 3600., ZTIMEA(:) ) )
         ZTIMEC(:) = real( nint( ZTIMEC(:) / PDTCONV ) ) * PDTCONV
         ZTIMEC(:) = max( PDTCONV, ZTIMEC(:) ) ! necessary if PDTCONV > 1800
      else
         call physeterror('convect_deep', 'adjustment time '//trim(deep_timeconv)//' unsupported by '//trim(deep_name))
         return
      endif

      !*           5.2    Compute melting level
      !                   ----------------------

      IML(:) = IKB
      do JK = IKE, IKB, -1
         where( ZTT(:,JK) <= RTT )  IML(:) = JK
      enddo

      call CONVECT_DOWNDRAFT2( ICONV, KLEV,                               &
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

      call CONVECT_PRECIP_ADJUST( ICONV, KLEV,                              &
           & ZPRES,ZUMF, ZUER, ZUDR, ZUPR, ZUTPR, ZURW,&
           & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,            &
           & ZPREF, ZTPR, ZMIXF, ZDTEVR,               &
           & ILFS, IDBL, ILCL, ICTL, IETL,             &
           & ZDTEVRF                                   )


      !*           7.     Determine adjusted environmental values assuming
      !                   that all available buoyant energy must be removed
      !                   within an advective time step ZTIMEC.
      !                   ---------------------------------------------------

      call CONVECT_CLOSURE2( ICONV, KLEV,                                &
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
      !                   -----------------------------------------------------


      !*           8.1    Grid scale tendencies
      !                   ---------------------

      ! in order to save memory, the tendencies are temporarily stored
      ! in the tables for the adjusted grid-scale values

      do JK = IKB, IKE
         ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
              & * ( ZPRES(:,JK) / RATM ) ** ZRDOCP  ! change theta in temperature
         ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
              &                            / ZTIMEC(:)

         ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
         ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:)

         ZPRLFLX(:,JK) = ZPRLFLX(:,JK) / ( RHOH2O * ZDXDY(:) )
         ZPRSFLX(:,JK) = ZPRSFLX(:,JK) / ( RHOH2O * ZDXDY(:) )
      enddo


      ZPRLFLX(:,IKB) = ZPRLFLX(:,IKB+1)
      ZPRSFLX(:,IKB) = ZPRSFLX(:,IKB+1)


      !*           8.2    Apply conservation correction
      !                   -----------------------------
      FORCE_CONSERVATION: if (CONSERV_CORRECT) then

         ! Compute moisture integral
         JKM = maxval( ICTL(:) )
         ZWORK2(:) = _ZERO_
         ZWORK2B(:)= _ZERO_
         do JK = IKB+1, JKM
            JKP = JK + 1
            do JI = 1, ICONV
               ZW1 = _HALF_ *  (ZPRES(JI,max(JK-1,IKB+1)) - ZPRES(JI,JKP)) / RG
               ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) *   &
                    &                   ZW1
               ZWORK2B(JI) = ZWORK2B(JI) + ZW1
            enddo
         enddo
         ! Compute moisture correction factor (integrates to surface precip.)
         do JI = 1, ICONV
            if ( ZTPR(JI) > _ZERO_) then
               ZWORK2(JI) = ( ZTPR(JI) / ZDXDY(JI) + ZWORK2(JI) ) / ZWORK2B(JI)
            endif
         enddo
         ! Apply moisture correction
         do JK = JKM, IKB+1, -1
            do JI = 1, ICONV
               if ( ZTPR(JI) > _ZERO_ .and. JK <= ICTL(JI) ) then
                  ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)
               endif
            enddo
         enddo
         ! Compute enthalpy integral
         ZWORK2(:) = _ZERO_
         ZWORK2B(:)= _ZERO_
         do JK = IKB+1, JKM
            JKP = JK + 1
            do JI = 1, ICONV
               ZW1 = _HALF_ *  (ZPRES(JI,max(JK-1,IKB+1)) - ZPRES(JI,JKP)) / RG
               ZWORK2(JI) = ZWORK2(JI) + (                                    &
                    &  ( RCPD + RCPV * ZRW(JI,JK) )* ZTHC(JI,JK)              &
                    &  + RCPV * ZTT(JI,JK) * (ZRVC(JI,JK)+ZRCC(JI,JK)+ZRIC(JI,JK))         &
                    &  - ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK)   &
                    &  - ( RLSTT + ( RCPV - RCS ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK)   &
                    &                        ) * ZW1
               ZWORK2B(JI) = ZWORK2B(JI) + ZW1
            enddo
         enddo
         ! Compute enthalpy correction factor (integrates to surface precip.)
         do JI = 1, ICONV
            if ( ZTPR(JI) > _ZERO_) then
               ZWORK2(JI) = ( ZTPR(JI) / ZDXDY(JI) *                          &
                    &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,IKB) - RTT ) ) - ZWORK2(JI) )     &
                    &                           / ZWORK2B(JI)
            endif
         enddo
         ! Apply enthalpy correction
         do JK = JKM, IKB+1, -1
            do JI = 1, ICONV
               if ( ZTPR(JI) > _ZERO_ .and. JK <= ICTL(JI) ) then
                  ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2(JI) / ( RCPD + RCPV * ZRW(JI,JK) )
               endif
            enddo
         enddo

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

      do JK = IKB, IKE
         do JI = 1, ICONV
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
         enddo
      enddo




      !*           8.3    Convective rainfall tendency
      !                   ----------------------------

      ! liquid and solid surface rainfall tendency in m/s
      ZTPR(:)   = ZTPR(:) / ( RHOH2O * ZDXDY(:) ) ! total surf precip
      ZSPR(:)   = ZSPR(:) / ( RHOH2O * ZDXDY(:) ) ! solid surf precip
      ZTPR(:)   = ZTPR(:) - ZSPR(:) ! compute liquid part

      !PV remet en kgm-2s-1 etant donne que c'est fait dans ccdiagnostics
      do JI = 1, ICONV
         JL = IJINDEX(JI)
         PPRLTEN(JL) = ZTPR(JI)*RHOH2O
         PPRSTEN(JL) = ZSPR(JI)*RHOH2O
      enddo

      !                   Cloud base and top levels
      !                   -------------------------
      !PV change kcltop and kclbas from index to height
      ILCL(:) = min( ILCL(:), ICTL(:) )
      do JI = 1, ICONV
         JL = IJINDEX(JI)
         KCLTOP(JL) = ZZ(JI,ICTL(JI))
         KCLBAS(JL) = ZZ(JI,ILCL(JI))
         if(GWORK1(JI)) PKKFC(JL) = 1.
      enddo


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
      do JI = 1, ICONV
         JL = IJINDEX(JI)
         ZTIME(JL)  =  ZWORK2(JI)
         ZWORK2(JL) =  max( ZTIME(JL), PDTCONV )
         if ( GTRIG(JL) )  KCOUNT(JL) = nint( ZWORK2(JL) / PDTCONV )
         if ( GTRIG(JL) .and. PPRLTEN(JL)<1.E-14 ) KCOUNT(JL) = 0
      enddo


      !*           8.6    Compute convective tendencies for Tracers
      !                   ------------------------------------------

      if ( OCH1CONV ) then

         allocate( ZCH1(ICONV,IKS,KCH1) )
         allocate( ZCH1C(ICONV,IKS,KCH1) )
         allocate( ZWORK3(ICONV,KCH1) )

         do JK = IKB, IKE
            do JI = 1, ICONV
               JL = IJINDEX(JI)
               ZCH1(JI,JK,:) = PCH1(JL,JK,:)
            enddo
         enddo

         call CONVECT_CHEM_TRANSPORT1(ICONV, KLEV, KCH1, ZCH1, ZCH1C,      &
              & IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
              & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
              & ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB)

         do JK = IKB, IKE
            do JN = 1, KCH1
               ZCH1C(:,JK,JN) = ( ZCH1C(:,JK,JN)- ZCH1(:,JK,JN) ) / ZTIMEC(:)
            enddo
         enddo


         !*           8.7    Apply conservation correction
         !                   -----------------------------

         ! Compute vertical integrals

         JKM = maxval( ICTL(:) )
         ZWORK3(:,:) = _ZERO_
         do JK = IKB+1, JKM
            JKP = JK + 1
            do JI = 1, ICONV
               ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) *                 &
                    &   _HALF_ * (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
            enddo
         enddo

         ! Mass error (integral must be zero)

         do JI = 1, ICONV
            if ( ZTPR(JI) > _ZERO_) then
               JKP = ICTL(JI)
               ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
                    & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
               ZWORK3(JI,:) = ZWORK3(JI,:) * ZW1
            endif
         enddo

         ! Apply uniform correction but assure positive mass at each level

         do JK = JKM, IKB+1, -1
            do JI = 1, ICONV
               if ( ZTPR(JI) > _ZERO_ .and. JK <= ICTL(JI) ) then
                  ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
                  ! ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
               endif
            enddo
         enddo

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

         do JK = IKB, IKE
            do JI = 1, ICONV
               JL = IJINDEX(JI)
               PCH1TEN(JL,JK,:) = ZCH1C(JI,JK,:)
            enddo
         enddo
      endif

      !*           8.8    Compute convective tendencies for wind
      !                   --------------------------------------

      if ( OUVCONV ) then

         allocate( ZUC(ICONV,IKS) )
         allocate( ZVC(ICONV,IKS) )

         call CONVECT_UV_TRANSPORT_DEEP1(ICONV, KLEV, ZU, ZV, ZUC, ZVC, &
              & IDPL, IPBL, ILCL, ICTL, ILFS, IDBL, IDDT, &
              & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,   &
              & ZDXDY, ZMIXF, ZLMASS, ZWSUB,  &
              & IFTSTEPS, ZTIMC, ITSTEP, GWORK1)

         do JK = IKB, IKE
            ZUC(:,JK) = ( ZUC(:,JK)- ZU(:,JK) ) / ZTIMEC(:)
            ZVC(:,JK) = ( ZVC(:,JK)- ZV(:,JK) ) / ZTIMEC(:)
         enddo

         !*           8.9    Apply conservation correction
         !                   -----------------------------

         ! Compute vertical integrals

         JKM = maxval( ICTL(:) )
         ZWORK2(:) = _ZERO_
         ZWORK2B(:)= _ZERO_
         do JK = IKB+1, JKM
            JKP = JK + 1
            do JI = 1, ICONV
               ZW1 = _HALF_ *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
               ZWORK2(JI) = ZWORK2(JI) + ZUC(JI,JK) * ZW1
               ZWORK2B(JI)= ZWORK2B(JI)+ ZVC(JI,JK) * ZW1
            enddo
         enddo

         !  error (integral must be zero)

         do JI = 1, ICONV
            if ( ZTPR(JI) > _ZERO_) then
               JKP = ICTL(JI)
               ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
                    & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
               ZWORK2(JI) = ZWORK2(JI) * ZW1
               ZWORK2B(JI)= ZWORK2B(JI)* ZW1
            endif
         enddo

         ! Apply uniform correction

         do JK = JKM, IKB+1, -1
            do JI = 1, ICONV
               if ( ZTPR(JI) > _ZERO_ .and. JK <= ICTL(JI) ) then
                  ZUC(JI,JK) = ZUC(JI,JK) - ZWORK2(JI)
                  ZVC(JI,JK) = ZVC(JI,JK) - ZWORK2B(JI)
               endif
            enddo
         enddo

         ! extend tendencies to first model level

         ! DO JI = 1, ICONV
         !    ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
         !    ZUC(JI,IKB)  = ZUC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
         !    ZUC(JI,IKB+1)= ZUC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
         !    ZVC(JI,IKB)  = ZVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
         !    ZVC(JI,IKB+1)= ZVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
         ! ENDDO


         do JK = IKB, IKE
            do JI = 1, ICONV
               JL = IJINDEX(JI)
               PUTEN(JL,JK)   = ZUC(JI,JK)
               PVTEN(JL,JK)   = ZVC(JI,JK)
            enddo
         enddo
      endif


      !*           9.     Write up- and downdraft mass fluxes
      !                   ------------------------------------

      do JK = IKB, IKE
         ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
         ZDMF(:,JK)  = ZDMF(:,JK) / ZDXDY(:)
      enddo
      do JK = IKB+1, IKE
         ZUDR(:,JK)  = ZUDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )! detrainment for ERA40
         ZDDR(:,JK)  = ZDDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )
      enddo
      ZWORK2(:) = _ONE_
      do JK = IKB, IKE
         do JI = 1, ICONV
            JL = IJINDEX(JI)
            if ( PPRLTEN(JL) < 1.E-14 ) ZWORK2(JL) = _ZERO_
            PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)
            PDMF(JL,JK) = ZDMF(JI,JK) * ZWORK2(JL)
            PUDR(JL,JK) = ZUDR(JI,JK) * ZWORK2(JL)
            PDDR(JL,JK) = ZDDR(JI,JK) * ZWORK2(JL)
            PURV(JL,JK) = PURV(JL,JK) * ZWORK2(JL)
            !       PURCI(JL,JK)= PURCI(JL,JK)* ZWORK2(JL)
            !PV calculate cloud fraction -same recipe as in  KFC
            if ( ZWU(JI,JK) > 1.E-1 ) then
               PCF(JL,JK) = ZUMF(JI,JK)/(ZRHO(JI,JK)*ZWU(JI,JK)) * ZWORK2(JL)
            else
               PCF(JL,JK) = 0.0
            endif
            ! calculate updraft areas from cloud fractions
            PAREAUP(JL,JK)=PCF(JL,JK)*ZDXDY(JI)
            PWUMAXOUT(JL) = max( ZWU(JI,JK), PWUMAXOUT(JL) )

            !PV output "grid-scale" condensate rather than "in-cloud"  condensate
            PURC(JL,JK)= PURC(JL,JK)* PCF(JL,JK) *ZWORK2(JL)
            PURI(JL,JK)= PURI(JL,JK)* PCF(JL,JK) *ZWORK2(JL)
         enddo
      enddo


      !*           10.    Deallocate all local arrays
      !                   ---------------------------

      ! downdraft variables

      deallocate( ZDMF )
      deallocate( ZDER )
      deallocate( ZDDR )
      deallocate( ZDTHL )
      deallocate( ZDRW )
      deallocate( ZLMASS )
      deallocate( ZMIXF )
      deallocate( ZTPR )
      deallocate( ZSPR )
      deallocate( ZDTEVR )
      deallocate( ZPREF )
      deallocate( IML )
      deallocate( ILFS )
      deallocate( IDBL )
      deallocate( IDDT )
      deallocate( ZDTEVRF )
      deallocate( ZPRLFLX )
      deallocate( ZPRSFLX )

      !   closure variables

      deallocate( ZTIMEA )
      deallocate( ZTIMEC )
      deallocate( ZTHC )
      deallocate( ZRVC )
      deallocate( ZRCC )
      deallocate( ZRIC )
      deallocate( ZWSUB )

      if ( OCH1CONV ) then
         deallocate( ZCH1 )
         deallocate( ZCH1C )
         deallocate( ZWORK3 )
      endif
      if ( OUVCONV ) then
         deallocate( ZUC )
         deallocate( ZVC )
      endif

   endif

   !    vertical index

   deallocate( IDPL )
   deallocate( IPBL )
   deallocate( ILCL )
   deallocate( ICTL )
   deallocate( IETL )

   ! grid scale variables

   deallocate( ZZ )
   deallocate( ZPRES )
   deallocate( ZDPRES )
   deallocate( ZU )
   deallocate( ZV )
   deallocate( ZRHO )
   deallocate( ZTT )
   deallocate( ZTH )
   deallocate( ZTHV )
   deallocate( ZTHL )
   deallocate( ZTHES )
   deallocate( ZRW )
   deallocate( ZRV )
   deallocate( ZRC )
   deallocate( ZRI )
   deallocate( ZDXDY )

   ! updraft variables

   deallocate( ZUMF )
   deallocate( ZUER )
   deallocate( ZUDR )
   deallocate( ZUTHL )
   deallocate( ZUTHV )
   deallocate( ZWU )
   deallocate( ZURW )
   deallocate( ZURC )
   deallocate( ZURI )
   deallocate( ZURR )
   deallocate( ZURS )
   deallocate( ZUPR )
   deallocate( ZUTPR )
   deallocate( ZTHLCL )
   deallocate( ZTLCL )
   deallocate( ZRVLCL )
   deallocate( ZWLCL )
   deallocate( ZZLCL )
   deallocate( ZTHVELCL )
   deallocate( ZMFLCL )
   deallocate( ZCAPE )

   ! work arrays

   deallocate( IINDEX )
   deallocate( IJINDEX )
   deallocate( IJSINDEX )
   deallocate( GTRIG1 )


end subroutine CONVECT_DEEP2
