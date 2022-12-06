!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module cnv_options
   implicit none 
   public
   save
   
   !# Select closures for shallow convection
   !# * 'CAPE'
   !# * 'EQUILIBRIUM'
   character(len=16) :: bkf_closures    = 'CAPE'
   namelist /convection_cfgs/ bkf_closures
   character(len=*), parameter :: BKF_CLOSURES_OPT(2) = (/ &
        'CAPE       ', &
        'EQUILIBRIUM'  &
        /)

   !# Select formulation of fractional detrainment rate for shallow convection
   !# * 'BECHTOLD01'
   !# * 'CUIJPERS95'
   !# * 'DEROOY10'
   character(len=16) :: bkf_detrains    = 'BECHTOLD01'
   namelist /convection_cfgs/ bkf_detrains
   character(len=*), parameter :: BKF_DETRAINS_OPT(3) = (/ &
        'BECHTOLD01', &
        'CUIJPERS95', &
        'DEROOY10  '  &
        /)

   !# Select formulation of fractional entrainment rate for shallow convection
   !# * 'BECHTOLD01'
   !# * 'BECHTOLD08'
   !# * 'DEROOY11'
   !# * 'SIEBESMA03'
   character(len=16) :: bkf_entrains    = 'BECHTOLD01'
   namelist /convection_cfgs/ bkf_entrains
   character(len=*), parameter :: BKF_ENTRAINS_OPT(4) = (/ &
        'BECHTOLD01', &
        'BECHTOLD08', &
        'DEROOY11  ', &
        'SIEBESMA03'  &
        /)

   !# Evaporate detrained condensate in shallow convection
   logical           :: bkf_evaps       = .false.
   namelist /convection_cfgs/ bkf_evaps

   !# Number of species for convective transport (never tested)
   integer           :: bkf_kch         = 0
   namelist /convection_cfgs/ bkf_kch

   !# Number of additional ensemble members (max 3) for deep bkf convection
   integer           :: bkf_kens        = 0
   namelist /convection_cfgs/ bkf_kens

   !# Temperature perturbations at LCL for triggering bkf_shallow
   !# An ensemble of shall. cumuli will be generated using perturbations
   !# starting from  bkf_tperts(1) to bkf_tperts(2) 
   !# with increment bkf_tperts(3) 
   real           :: bkf_tperts(3) = (/0.2, 0.2, 0./)
   namelist /convection_cfgs/ bkf_tperts

   !# Cloud radii at LCL for bkf_shallow
   !# from bkf_rads(1) to bkf_rads(2) 
   !# with increment bkf_rads(3)
   real           :: bkf_rads(3)  = (/50., 50., 0./)
   namelist /convection_cfgs/ bkf_rads

   !# Take ice phase into account in deep bkf (yes=1)
   integer           :: bkf_kice         = 1
   namelist /convection_cfgs/ bkf_kice

   !# Limit vertical computation by ktdia-1 levels
   integer           :: bkf_ktdia        = 1
   namelist /convection_cfgs/ bkf_ktdia

   !# Activate convective transport of species for deep and shallow bkf
   logical           :: bkf_lch1conv     = .false.
   namelist /convection_cfgs/ bkf_lch1conv

   !# Allow downdrafts in deep bkf
   logical           :: bkf_ldown        = .true.
   namelist /convection_cfgs/ bkf_ldown

   !# Activate shallow convective momentum transport
   logical           :: bkf_lshalm       = .false.
   namelist /convection_cfgs/ bkf_lshalm

   !# Deep convection scheme name
   !# * 'NIL     ' :
   !# * 'SEC     ' :
   !# * 'KFC     ' :
   !# * 'KFC2    ' :
   !# * 'BECHTOLD' :
   character(len=16) :: deep            = 'nil'
   namelist /convection_cfgs/ deep
   character(len=*), parameter :: DEEP_OPT(5) = (/ &
        'NIL     ', &
        'SEC     ', &
        'KFC     ', &
        'KFC2    ', &
        'BECHTOLD'  &
        /)

   !# Treat convective clouds as cloud objects
   logical :: deep_cloudobj = .false.
   namelist /convection_cfgs/ deep_cloudobj

   !# Decay timescale for convective cloud objects (seconds)
   real :: deep_codecay = 600.
   namelist /convection_cfgs/ deep_codecay

   !# Conservation corrections for deep convective scheme
   !# * 'NIL   ' : No conservation correction applied
   !# * 'TEND  ' : Temperature and moisture tendencies corrected
   !# * 'PRECIP' : Surface precipitation rate corrected
   character(len=16) :: deep_conserve = 'PRECIP'
   namelist /convection_cfgs/ deep_conserve
   namelist /convection_cfgs_p/ deep_conserve
   character(len=*), parameter :: DEEP_CONSERVE_OPT(3) = (/ &
        'NIL   ', &
        'TEND  ', &
        'PRECIP'  &
        /)

   !# Time over which all LCL air can be drawn into up/downdrafts
   !# *  USER DEFINED (e.g. "1h") (units D,H,M,S,P)
   !# * 'NIL        ' : No time limit on entraniment into the up/downdrafts
   !# * 'TIMECONV   ' : Use convective time scale for entrainment limit
   !# * 'TIMEREFRESH' : Use refresh time scale for entrainment limit
   character(len=16) :: deep_timeent     = 'timeconv'
   real              :: deep_timeent_sec = -1.
   namelist /convection_cfgs/ deep_timeent
   character(len=*), parameter :: DEEP_TIMEENT_OPT(3) = (/ &
        'NIL        ', &
        'TIMECONV   ', &
        'TIMEREFRESH' &
        /)

   !# Convective adjustment time-scale for deep convection
   !# * USER DEFINED (e.g. "1h") (units D,H,M,S,P)
   !# * 'ADVECTIVE'
   !# * 'KFCTAUCAPE'
   !# * 'BECHTOLD09' : Adjustment described in Bechtold (2009; ECMWF Training)
   character(len=16) :: deep_timeconv     = '1h'
   real              :: deep_timeconv_sec = -1.
   namelist /convection_cfgs/ deep_timeconv
   character(len=*), parameter :: DEEP_TIMECONV_OPT(4) = (/ &
        'USER DEFINED (e.g. "1h")', &
        'ADVECTIVE               ', &
        'KFCTAUCAPE              ', &
        'BECHTOLD09              ' &
        /) !#TODO: remove 'USER DEFINED'?

   !# Refresh time-scale for deep convection
   !# max = deep_timeconv imposed in code
   !# min = model timestep imposed in code
   !# * USER DEFINED (e.g. "1h") (units D,H,M,S,P)
   !# * 'ADVECTIVE'
   !# * 'TIMECONV'
   character(len=16) :: deep_timerefresh     = 'timeconv'
   real              :: deep_timerefresh_sec = -1.
   namelist /convection_cfgs/ deep_timerefresh
   character(len=*), parameter :: DEEP_TIMEREFRESH_OPT(3) = (/ &
        'USER DEFINED (e.g. "1h")', &
        'ADVECTIVE               ', &
        'TIMECONV                '  &
        /) !#TODO: remove 'USER DEFINED'?

   !# Minimum depth of conv. updraft for KFC  trigger (m)
   real              :: kfcdepth        = 4000.
   namelist /convection_cfgs/ kfcdepth

   !# Maximum depth of the downdraft detrainment layer (Pa) for 'kfc2'
   real              :: kfcdpdd         = 10000.
   namelist /convection_cfgs/ kfcdpdd

   !# 
   !# Generate wind tendencies (CMT) in deep=KFC,KFC2
   !# * 'NIL      ': No wind tendencies applied
   !# * 'ECMWF_PH2': ECMWF approach over anvil only, as implemented in phase 2 with bug 
   !# * 'ECMWF    ': ECMWF approach, debugged and applied over whole cloud (KFC2 only)
   !# * 'GKI      ': GKI approach (KFC2 only)
   character(len=16) :: cmt_type        = 'NIL'
   namelist /convection_cfgs/ cmt_type
   character(len=*), parameter :: CMT_TYPE_OPT(4) = (/ &
        'NIL      ', &
        'ECMWF_PH2', &
        'ECMWF    ', &
        'GKI      ' &
        /)
   integer, parameter :: CMT_NONE = 0
   integer, parameter :: CMT_ECMWF_PH2 = 1
   integer, parameter :: CMT_ECMWF = 2
   integer, parameter :: CMT_GKI = 3
   integer :: cmt_type_i = CMT_NONE

   !# CMT lambda factor for cmt_type=ECMWF* formulations
   real              :: cmt_ecmwf_lambda = 2.
   namelist /convection_cfgs/ cmt_ecmwf_lambda
   
   !# CMT GKI for updraft coef for cmt_type=GKI
   real              :: cmt_gki_cup      = 0.7
   namelist /convection_cfgs/ cmt_gki_cup
   
   !#  CMT GKI for downdraft coef for cmt_type=GKI
   real              :: cmt_gki_cdown    = 0.
   namelist /convection_cfgs/ cmt_gki_cdown

   !# Compute production terms for Kain-Fritsch scheme
   logical           :: kfcprod         = .false.
   namelist /convection_cfgs/ kfcprod

   !# Initial convective updraft radius in KFC scheme(m)
   real              :: kfcrad          = 1500.
   namelist /convection_cfgs/ kfcrad

   !# Convective updraft radius over water in KFC scheme(m)
   real              :: kfcradw         = -1.
   namelist /convection_cfgs/ kfcradw

   !# Varies convective timescale as a function of CAPE for Kain-Fritsch scheme
   !# KFCTAUCAPE = time1, time2, cmean, dcape
   !# * time1 (s): max kfctimec 
   !# * time2 (s): min kfctimec 
   !# * cmean (J/Kg): cape value at which kfctimec will be mean of time1 and time2 
   !# * dcape (J/Kg): decrease in kfctimec from time1 to time2 will occur over range cmean-dcape to cmean+dcape 
   real              :: kfctaucape(4)   = (/-1., -1., -1., -1./)
   namelist /convection_cfgs/ kfctaucape

   !# Trigger parameter of Kain-Fritsch convection scheme (WKLCL).
   !# Trigger parameter will increase from kfctrigw(3) to kfctrigw(4) [m/s]
   !# between wstar values kfctrigw(1) and kfctrigw(2)
   real              :: kfctrigw(4)     = (/0., 0., 0., 0./)
   namelist /convection_cfgs/ kfctrigw

   !# Trigger parameter of Kain-Fritsch convection scheme (WKLCL).
   !# Trigger parameter will increase from kfctrig4(3) to kfctrig4(4) [m/s]
   !# between timestep kfctrig4(1) and timestep kfctrig4(2)
   real              :: kfctrig4(4)     = (/0., 0., 0.05, 0.05/)
   namelist /convection_cfgs/ kfctrig4

   !# Nominal resolution for which KFCTRIG4 is set.
   !# This is inactive if value <= 0.
   real              :: kfctriga        = -1.0
   namelist /convection_cfgs/ kfctriga

   !# Over land and lakes we keep the value set by the "ramp" above over sea water:
   !# * for |lat| >= TRIGLAT(2) we keep value set by the "ramp" KFCTRIG4
   !# * for |lat| <= TRIGLAT(1) we use the new value KFCTRIGL [m/s]
   !# * and linear interpolation in between TRIGLAT(1) and TRIGLAT(2)
   real              :: kfctrigl        = 0.05
   namelist /convection_cfgs/ kfctrigl

   !# Logical key for variation of the trigger function depending on latitude and land-sea-lake mask
   logical            :: kfctriglat      = .false.
   namelist /convection_cfgs/ kfctriglat

   !# Relaxation timescale for trigger velocity
   real               :: kfctrigtau      = -1.
   namelist /convection_cfgs/ kfctrigtau

   !# Switch for mid-level convection
   !# * 'NIL' : No mid-level convective scheme
   !# * 'KF ' : Kain-Fritsch-based mid-level convective scheme
   character(len=16) :: mid            = 'nil'
   namelist /convection_cfgs/ mid
   character(len=*), parameter :: MID_OPT(2) = (/ &
        'NIL', &
        'KF ' &
        /)

   !# Conservation corrections for mid-level convective scheme
   !# * 'NIL   ' : No conservation correction applied
   !# * 'TEND  ' : Temperature and moisture tendencies corrected
   !# * 'PRECIP' : Surface precipitation rate corrected
   character(len=16) :: mid_conserve = 'PRECIP'
   namelist /convection_cfgs/ mid_conserve
   namelist /convection_cfgs_p/ mid_conserve
   character(len=*), parameter :: MID_CONSERVE_OPT(3) = (/ &
        'NIL   ', &
        'TEND  ', &
        'PRECIP'  &
        /)

   !# Minimum cloud depth for mid-level convection (m)
   real              :: mid_depth = 2000.
   namelist /convection_cfgs/ mid_depth
   namelist /convection_cfgs_p/ mid_depth

   !# Downdraft detrainment depth for mid-level convection (Pa)
   real              :: mid_dpdd    = 6000.
   namelist /convection_cfgs/ mid_dpdd
   namelist /convection_cfgs_p/ mid_dpdd

   !# Fraction of environmental mass flux that enters updrafts
   character(len=16) :: mid_emffrac = 'all'
   namelist /convection_cfgs/ mid_emffrac
   namelist /convection_cfgs_p/ mid_emffrac
   character(len=*), parameter :: MID_EMFFRAC_OPT(4) = (/ &
        'ALL ', &
        'TLCL', &
        'TDD ', &
        'QLCL' &
        /)

   !# Modulation of the minimum environmental mass flux for mid-level convection
   character(len=16) :: mid_emfmod = 'nil'
   namelist /convection_cfgs/ mid_emfmod
   namelist /convection_cfgs_p/ mid_emfmod
   character(len=*), parameter :: MID_EMFMOD_OPT(4) = (/ &
        'NIL     ', &
        'LATITUDE', &
        'LATMOD1 ', &
        'LATMOD2 ' &
        /)

   !# Minimum environmental mass flux for mid-level convection (kg/s)
   real              :: mid_minemf = 1e7
   namelist /convection_cfgs/ mid_minemf
   namelist /convection_cfgs_p/ mid_minemf

   !# Maximum deep CAPE (J/kg/m2) for mid-level convective triggering
   real              :: mid_maxcape = -1
   namelist /convection_cfgs/ mid_maxcape
   namelist /convection_cfgs_p/ mid_maxcape

   !# Minimum parcel departure level for mid-level convection (m)
   real              :: mid_minbase  = 500.
   namelist /convection_cfgs/ mid_minbase
   namelist /convection_cfgs_p/ mid_minbase

   !# Precipitation efficiency for the mid-level convective scheme (fractional)
   !# * 'USER DEFINED (e.g. 0.4)' : Use a fixed value prescribed in the namelist
   !# * 'BLUESTEIN83            ' : Use subcloud precipitable water deficit 
   !# * 'CONDSE                 ' : Use efficiency estimate from gridscale scheme
   character(len=16)  :: mid_peff       = '0.4'
   real               :: mid_peff_const = -1.
   namelist /convection_cfgs/ mid_peff
   namelist /convection_cfgs_p/ mid_peff
   character(len=*), parameter :: MID_PEFF_OPT(5) = (/ &
        'BLUESTEIN83', &
        'CONDSE     ', &
        'TLCL       ', &
        'TDD        ', &
        'QLCL       ' &
        /)

   !# Switch for shallow convection
   !# * 'NIL'
   !# * 'KTRSNT'
   !# * 'BECHTOLD'
   character(len=16) :: shal            = 'nil'
   namelist /convection_cfgs/ shal
   character(len=*), parameter :: SHAL_OPT(3) = (/ &
        'KTRSNT   ', &
        'BECHTOLD ', &
        'NIL      '  &
        /)

   !# Conservation corrections for shallow convective scheme
   !# * 'NIL ' : No conservation correction applied
   !# * 'TEND' : Temperature and moisture tendencies corrected
   character(len=16) :: shal_conserve = 'NIL'
   namelist /convection_cfgs/ shal_conserve
   namelist /convection_cfgs_p/ shal_conserve
   character(len=*), parameter :: SHAL_CONSERVE_OPT(2) = (/ &
        'NIL ', &
        'TEND'  &
        /)

   !# Convective adjustment time-scale for shallow convection
   !# 'USER DEFINED (e.g. "3h") (units D,H,M,S,P)
   character(len=16) :: shal_timeconv     = '3h'
   real              :: shal_timeconv_sec = -1.
   namelist /convection_cfgs/ shal_timeconv
   character(len=*), parameter :: SHAL_TIMECONV_OPT(1) = (/ &
        'USER DEFINED (e.g. "3h")' &
        /) !#TODO: remove 'USER DEFINED'?

   !# Over land and lakes we keep the value set by the "ramp" above over sea water:
   !# * for |lat| >= TRIGLAT(2) we keep value set by the "ramp" KFCTRIG4
   !# * for |lat| <= TRIGLAT(1) we use the new value KFCTRIGL
   !# * and linear interpolation in between TRIGLAT(1) and TRIGLAT(2)
   real              :: triglat(2)      = 0.0
   namelist /convection_cfgs/ triglat

contains

   function cnv_options_init() result(F_istat)
      implicit none 
      integer :: F_istat
#include <rmnlib_basics.hf>
      F_istat = RMN_OK
      return
   end function cnv_options_init

end module cnv_options
