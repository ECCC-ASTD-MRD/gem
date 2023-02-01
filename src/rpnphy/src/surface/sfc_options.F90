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

module sfc_options
   use, intrinsic :: iso_fortran_env, only: INT64
   implicit none
   public
   save
 
   !#
   integer, parameter :: CLASS_URB = 21 !# Class "urban"
   integer, parameter :: NCLASS    = 26 !# NUMBER OF CLASSES FOR NATURAL COVERS
   integer, parameter :: NCLASSURB = 12 !# NUMBER OF URBAN CLASSES
   integer, parameter :: NL = 3 !# nombre de niveaux dans la glace marine


   !# ----  FOR SVS LAND SCHEME ------
   integer, parameter :: MAX_NL_SVS = 50 !# maximum number of soil layers specified by user
   integer, parameter :: NL_SVS_DEFAULT = 7 ! default number of soil layers
   ! default depth of default soil layers in [METERS]
   real, dimension(NL_SVS_DEFAULT):: DP_SVS_DEFAULT =  (/ 0.05, 0.1, 0.2, 0.4, 1.0, 2.0, 3.0 /)
   !#----------------------------------

   real, parameter :: CRITEXTURE = 0.1
   real, parameter :: CRITLAC    = 0.01
   real, parameter :: CRITMASK   = 0.001
   real, parameter :: CRITSNOW   = 0.0001
   real, parameter :: CRITWATER  = 0.001
   real, parameter :: HIMIN      = 0.001
   real, parameter :: MINICEDP   = 0.05
   real, parameter :: N0RIB = 1.0E-5
   real, parameter :: SNOH0 = 0.1
   real, parameter :: VAMIN = 1.0e-4
   real, parameter :: Z0GLA = 0.0003

   !#
   logical           :: atm_external = .false.
   logical           :: atm_tplus  = .false.
   logical           :: climat     = .false.
   logical           :: cplocn     = .false.
   real              :: delt       = 0.
   integer(INT64)  :: jdateo     = 0
   logical           :: rad_off    = .false.
   logical           :: radslope   = .false.
   logical           :: update_alwater = .false.
   logical           :: z0veg_only = .false.
   logical           :: thermal_stress = .false.
   integer           :: kntveg     = -1
   logical           :: timings_L  = .false.
   integer           :: nphyoutlist = -1
   character(len=32), pointer :: phyoutlist_S(:) => NULL()

   !# Surface layer coefficients for exchange with PBL
   real :: bh91_a, bh91_b, bh91_c, bh91_d, d97_as, dg92_ci, l07_ah, l07_am

   !# Adjust surface temperature over snow after reading (coherency check)
   logical           :: adj_i0_snow = .true.
   namelist /surface_cfgs/ adj_i0_snow

   !# Prandtl number for neutral stability (initialized by SL module)
   real              :: beta        = 0.
   namelist /surface_cfgs/ beta

   !# Diurnal SST scheme
   !# * 'NIL    ' : No Diurnal SST scheme
   !# * 'FAIRALL' : #TODO: define
   character(len=16) :: diusst      = 'NIL'
   namelist /surface_cfgs/ diusst
   character(len=*), parameter :: DIUSST_OPT(2) = (/ &
        'NIL    ', &
        'FAIRALL'  &
        /)

   !# Diurnal SST scheme active over freshwater lakes if .true.
   logical           :: diusst_lakes = .true.
   namelist /surface_cfgs/ diusst_lakes

   !# Diurnal SST scheme active over ocean if .true.
   logical           :: diusst_ocean = .true.
   namelist /surface_cfgs/ diusst_ocean

   !# Diurnal SST scheme active warmlayer over freshwater lakes if .true.
   logical           :: diusst_warmlayer_lakes = .true.
   namelist /surface_cfgs/ diusst_warmlayer_lakes

   !# Depth of soil layers in [METERS] in SVS land surface scheme (schmsol=SVS)
   real :: dp_svs(MAX_NL_SVS) = -1.0
   namelist /surface_cfgs/ dp_svs

   !# Emissivity for ice (glacier and sea ice)
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: ice_emiss = '0.99'
   real              :: ice_emiss_const = -1.
   namelist /surface_cfgs/ ice_emiss

   !# Set water temperature of ice-covered lakes to 0C for points north of
   !# ice line if .true.
   !# needs an initialization file otherwise the model stops
   logical           :: icelac      = .false.
   namelist /surface_cfgs/ icelac

   !# Sea ice melting
   logical           :: icemelt     = .false.
   namelist /surface_cfgs/ icemelt

   !# Implicit surface fluxes if .true.; explicit fluxes if .false.
   logical           :: impflx      = .false.
   namelist /surface_cfgs/ impflx

   !# If .true. apply temporary fix to ISBA
   !# * timestep dependent KCOEF
   !# * No PSN factor for meting and freezing
   logical           :: isba_melting_fix = .false.
   namelist /surface_cfgs/ isba_melting_fix

   !# Use the vegetation-only roughness length to compute vegetation snow fraction
   logical           :: isba_snow_z0veg = .false.
   namelist /surface_cfgs/ isba_snow_z0veg

   !# Emissivity for bare soil (ISBA scheme only)
   !# * '_constant_' : A fixed floating point value used as a constant
   !# * 'CLIMATO'    : Value read from an input climatology file (EMIB)
   character(len=16) :: isba_soil_emiss = '0.95'
   real              :: isba_soil_emiss_const = -1.
   namelist /surface_cfgs/ isba_soil_emiss
   character(len=*), parameter :: ISBA_SOIL_EMISS_OPT(1) = (/ &
        'CLIMATO' &
        /)

   !# If .true., freeze precipitation reaching the ground in sub-zero conditions
   logical           :: isba_zr_freeze = .false.
   namelist /surface_cfgs/ isba_zr_freeze

   !# OBSOLETE, REPLACED by KHYD !!! WILL BE EVENTUALLY REMOVED
   !# Deepest active (permeable) soil layer in SVS land surface scheme (schmsol=SVS)
   integer           :: kdp    = -1
   namelist /surface_cfgs/ kdp

   !# Last/Deepest soil layer considered during the accumulation of
   !# lateral flow and drainage. Drainage is taken as the vertical flux
   !# leaving layer KHYD, and lateral flow as the sum of lateral flows from
   !# layers 1 to KHYD
   integer           :: khyd    = -1
   namelist /surface_cfgs/ khyd

   !# Vegetation field update frequency (units D,H,M,S,P)
   character(len=16) :: kntveg_S     = ''
   namelist /surface_cfgs/ kntveg_S

   !# Lead fraction for ice-covered lakes
   real              :: lake_leadfrac = 0.
   namelist /surface_cfgs/ lake_leadfrac

   !# Minimum fraction of leads in sea ice.&nbsp; Multiply ice fraction by (1.-leadfrac)
   real              :: leadfrac    = 0.03
   namelist /surface_cfgs/ leadfrac

   !# Limit snow depth to 10 cm for calculation of heat conductivity of snow
   !# over sea-ice and glacier if .true.
   logical           :: limsnodp    = .false.
   namelist /surface_cfgs/ limsnodp

   !# (coupling) fluxes over ocean are taken from ocean model if .true.
   logical           :: owflux      = .false.
   namelist /surface_cfgs/ owflux

   !# read-in land surface emissivity if .true.
   logical           :: read_emis     = .false.
   namelist /surface_cfgs/ read_emis

   !# read-in high vegetation roughness for SVS if .true.
   logical           :: read_z0vh     = .false.
   namelist /surface_cfgs/ read_z0vh

   !# Takes into account effect of ocean salinity on saturation specific
   !# humidity at ocean surface (boundary condition for LH flux calculation)
   logical           :: salty_qsat  = .false.
   namelist /surface_cfgs/ salty_qsat

   !# Land surface processes
   !# * 'NIL ' : No Land surface processes
   !# * 'ISBA' : Interaction Soil Biosphere Atmosphere (ISBA) land sfc scheme
   !# * 'SVS ' : Soil, Vegetation, and Snow (SVS) (Multibudget) land sfc scheme
   character(len=16) :: schmsol     = 'ISBA'
   namelist /surface_cfgs/ schmsol
   character(len=*), parameter :: SCHMSOL_OPT(3) = (/ &
        'NIL ', &
        'ISBA', &
        'SVS '  &
        /)

   !# Urban surface processes
   !# * 'NIL' : No Urban surface processes
   !# * 'TEB' : Town Energy Balance (TEB) urban scheme
   character(len=16) :: schmurb     = 'NIL'
   namelist /surface_cfgs/ schmurb
   character(len=*), parameter :: SCHMURB_OPT(2) = (/ &
        'NIL', &
        'TEB'  &
        /)

   !# Lake surface processes
   !# * 'NIL' :
   !# * 'FLAKE' :
   !# * 'CSLM' :
   character(len=16) :: schmlake    = 'NIL'
   namelist /surface_cfgs/ schmlake
   character(len=*), parameter :: SCHMLAKE_OPT(3) = (/ &
        'NIL  ', &
        'FLAKE' , &
        'CSLM ' &
        /)

   !# River surface processes
   !# * 'NIL' :
   character(len=16) :: schmriver    = 'NIL'
   namelist /surface_cfgs/ schmriver
   character(len=*), parameter :: SCHMRIVER_OPT(2) = (/ &
        'NIL   ', &
        'RIVERS' &
        /)

   !# Minimum Obukhov length (L) for glaciers
   real           :: sl_Lmin_glacier = -1.
   namelist /surface_cfgs/ sl_Lmin_glacier
   
   !# Minimum Obukhov length (L) for sea ice
   real           :: sl_Lmin_seaice = -1.
   namelist /surface_cfgs/ sl_Lmin_seaice
   
   !# Mimimum Obukhov length (L) for soil surfaces
   real           :: sl_Lmin_soil = -1.
   namelist /surface_cfgs/ sl_Lmin_soil

   !# Minimum Obukhov length (L) for water
   real           :: sl_Lmin_water = -1.
   namelist /surface_cfgs/ sl_Lmin_water

   !# Minimum Obukhov length (L) for town
   real           :: sl_Lmin_town = -1.
   namelist /surface_cfgs/ sl_Lmin_town
   
   !# Define bulk Ri values for near-neutral regime in the surface layer
   real           :: sl_rineutral = 0.
   namelist /surface_cfgs/ sl_rineutral

   !# Class of stability functions (stable case) to use in the surface layer
   !# * 'DELAGE97  ' : Use functions described by Delage (1997; BLM)
   !# * 'BELJAARS91' : Use functions described by Beljaars and Holtslag (1991; JAM)
   !# * 'LOCK07    ' : Use functions described by Lock (2007; Tech Report) employed at UKMO
   character(len=16) :: sl_func_stab = 'DELAGE97'
   namelist /surface_cfgs/ sl_func_stab
   character(len=*), parameter :: SL_FUNC_STAB_OPT(3) = (/ &
        'DELAGE97  ', &
        'BELJAARS91', &
        'LOCK07    '  &
        /)

   !# Class of stability functions (unstable case) to use in the surface layer
   !# * 'DELAGE92' : Use functions described by Delage and Girard (1992; BLM)
   !# * 'DYER74  ' : Use functions described by Dyer (1974; BLM)
   character(len=16) :: sl_func_unstab = 'DELAGE92'
   namelist /surface_cfgs/ sl_func_unstab
   character(len=*), parameter :: SL_FUNC_UNSTAB_OPT(2) = (/ &
        'DELAGE92', &
        'DYER74  '  &
        /)
 
   !# Use a reference roughness for surface layer calculations
   logical :: sl_z0ref = .false.
   namelist /surface_cfgs/ sl_z0ref

   !# Emissivity for snow
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: snow_emiss = '1.'
   real              :: snow_emiss_const = -1.
   namelist /surface_cfgs/ snow_emiss

   !#  Soil texture database/calculations for SVS land surface scheme
   !# * 'GSDE   '   : 8 layers of sand & clay info from Global Soil Dataset for ESMs (GSDE)
   !# * 'SLC    '   : 5 layers of sand & clay info from Soil Landscape of Canada (SLC)
   !# * 'SOILGRIDS' : 7 layers of sand & clay info from ISRIC ? World Soil Information
   character(len=16) :: soiltext    = 'GSDE'
   namelist /surface_cfgs/ soiltext
   character(len=*), parameter :: SOILTEXT_OPT(3) = (/ &
        'GSDE     ',  &
        'SLC      ',  &
        'SOILGRIDS' &
        /)

   !# If .true., SVS1 simulates soil freezing and thawing and its impact on hydrology
   logical           :: lsoil_freezing_svs1 = .false.
   namelist /surface_cfgs/ lsoil_freezing_svs1 

   !# If .true., SVS1 simulates water ponding at the surface
   logical           :: lwater_ponding_svs1 = .false.
   namelist /surface_cfgs/ lwater_ponding_svs1 


   !# Use snow albedo "I6" directly if .true.;
   !# Use snow age "XA" to calculate snow albedo if .false.
   logical           :: snoalb_anl  = .true.
   namelist /surface_cfgs/ snoalb_anl

   !# use dynamic calculation of z0h for bare ground + vegetation  for SVS if .true.
   logical           :: svs_dynamic_z0h     = .false.
   namelist /surface_cfgs/ svs_dynamic_z0h 

   
   !# use hrsurf based on soil texture for SVS if .true.
   logical           :: svs_hrsurf_sltext     = .false.
   namelist /surface_cfgs/ svs_hrsurf_sltext

   !# use local momentum (no snow) roughness for SVS if .true.
   logical           :: svs_local_z0m     = .false.
   namelist /surface_cfgs/ svs_local_z0m  
   
   !# Option to deal with melt due to rain on snow in SVS snow scheme
   !# BELAIR03: Original formulation from Belair et al. (2003)
   !# BELAIR03_DTGEM: Revised formulation of Belair et al. (2003) using the time step of GEM at the time when the parameterisation has been
   !developped (12 min). Recommended for coupled applications with short time steps (< 5 min). 
   !# NONE: Snow melt due to rain deactivated
   !# LEONARDINI21:  Assume all the energy brought by rain is used to melt the snowpack (see Leonardini et al. (2021)
   ! Options REV and L21 remvoe the large sensitivity to the time step found with DEF.  
   character(len=16) :: svs_snow_rain   = 'BELAIR03'
   namelist /surface_cfgs/ svs_snow_rain
   character(len=*), parameter :: SVS_SNOW_RAIN_OPT(4) = (/ &
        'BELAIR03         ',  &
        'BELAIR03_DTGEM   ',  &
        'NONE             ',  &
        'LEONARDINI21     ' &
        /)

   
   !# Limit temperature inversions to 8K/40m in surface layer if .true.
   logical           :: tdiaglim    = .false.
   namelist /surface_cfgs/ tdiaglim

   !# OPTION FOR CALCULATION of AVERAGE LAND SURFACE TEMPERATURE AND HUMIDITY IN SVS
   !# .FALSE. :  Area-average only calculation for sfc T and Hum.
   !# .TRUE.  :  Option that uses effective surface temperature  and specific humidity instead
   !             of composite (area-averaged only) counterparts in surface flux calculations (D. Deacu)
   logical           :: use_eff_surf_tq    = .false.
   namelist /surface_cfgs/ use_eff_surf_tq

   !# OPTION TO USE PHOTOSYNTHESIS CODE FOR STOMATAL RESISTANCE in SVS
   logical           :: use_photo = .true.
   namelist /surface_cfgs/ use_photo
 
   !# Factor multiplying stomatal resistance in ISBA
   real              :: veg_rs_mult = 1.
   namelist /surface_cfgs/ veg_rs_mult

   !#  VF definitions and mapping in SVS 
   !# * 'CLASSIC' : Same VF definitions as ISBA
   !# * 'CCILCECO' : New VF definitions used in SVS only
   !# with geo. fields generated using CCILC 2015 + Ecobiomes
   character(len=16) :: vf_type    = 'CLASSIC'
   namelist /surface_cfgs/ vf_type
   character(len=*), parameter :: VFTYPE_OPT(2) = (/ &
        'CLASSIC  ',  &
        'CCILCECO '   &
        /)

   !# Emissivity for water
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: water_emiss = '1.'
   real              :: water_emiss_const = -1.
   namelist /surface_cfgs/ water_emiss

   !# Use directional roughness length if .true.
   logical           :: z0dir       = .false.
   namelist /surface_cfgs/ z0dir

   !# Constant value of thermal roughness length (m) applied over water within
   !# latitudinal band defined by z0tlat
   real              :: z0hcon      = 4.0e-5
   namelist /surface_cfgs/ z0hcon

   !# Minimum value of momentum roughness length (m)
   real              :: z0min       = 1.5e-5
   namelist /surface_cfgs/ z0min

   !# Momentum roughness length formulation over water
   !# * 'CHARNOCK' : Standard Charnock clipped at high wind speed
   !# * 'BELJAARS' : #TODO: define
   !# * 'WRF1'     : ISFTCFLX=1 from WRF (Green and Zhang 2013)
   !# * 'WRF2'     : ISFTCFLX=2 from WRF (Green and Zhang 2013)
   character(len=16) :: z0mtype     = 'CHARNOCK'
   namelist /surface_cfgs/ z0mtype
   character(len=*), parameter :: Z0MTYPE_OPT(4) = (/ &
        'CHARNOCK', &
        'BELJAARS', &
        'WRF1    ', &
        'WRF2    ' &
        /)

   !# Roughness length for sea ice
   real              :: z0seaice    = 1.6e-4
   namelist /surface_cfgs/ z0seaice
   
   !# Thermal roughness length formulation over water
   !# * 'MOMENTUM' : Uses z0h = z0m (replaces key z0trdps300=.false.)
   !# * 'DEACU12'  : #TODO: define  (replaces key z0trdps300=.true.)
   !# * 'ECMWF'    : #TODO: define  (New formulation used by ECMWF)
   !# * 'WRF1'     : ISFTCFLX=1 from WRF (Green and Zhang 2013)
   !# * 'WRF2'     : ISFTCFLX=2 from WRF (Green and Zhang 2013)
   character(len=16) :: z0ttype     = 'MOMENTUM'
   namelist /surface_cfgs/ z0ttype
   character(len=*), parameter :: Z0TTYPE_OPT(5) = (/ &
        'MOMENTUM', &
        'DEACU12 ', &
        'ECMWF   ', &
        'WRF1    ', &
        'WRF2    ' &
        /)

   !# Thermal roughness length formulation over vegetation
   !# * 'FIXED' : Uses z0h = z0m 
   !# * 'ZILI95': evolves with u*
   character(len=16) :: z0tevol     = 'FIXED'
   namelist /surface_cfgs/ z0tevol
   character(len=*), parameter :: Z0TEVOL_OPT(2) = (/ &
        'FIXED ', &
        'ZILI95'  &
        /)

   !# Latitude (2 elements, in degrees) used to specify Z0T over water
   !# * If |lat| <= Z0TLAT(1) constant Z0T.
   !# * If |lat| >= Z0TLAT(2) Charnock's relation.
   !# * In between, linear interpolation is used.
   real              :: z0tlat(2)   = 0.
   namelist /surface_cfgs/ z0tlat

   !# Height at which to compute screen-level temperature (m)
   real              :: zt = 1.5
   namelist /surface_cfgs/ zt

   !# Height at which to compute anemomenter-level winds (m)
   real              :: zu = 10.
   namelist /surface_cfgs/ zu

   !# New urban surface parameters within SVS only (not used in TEB)
   logical           :: svs_urban_params = .false.
   namelist /surface_cfgs/ svs_urban_params
 
   !# Adjust wind diagnostic in TEB in the street  if .true.
   logical           :: urb_diagwind = .false.
   namelist /surface_cfgs/ urb_diagwind

   !# Adjust temperature diagnostic in TEB in the street  if .true.
   logical           :: urb_diagtemp = .false.
   namelist /surface_cfgs/ urb_diagtemp

contains

   function sfc_options_init() result(F_istat)
      use sfclayer, only: sl_get, SL_OK
      implicit none
      integer :: F_istat
#include <rmn/msg.h>
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      F_istat = sl_get('beta',beta)
      F_istat = min(sl_get('bh91_a',bh91_a),F_istat)
      F_istat = min(sl_get('bh91_b',bh91_b),F_istat)
      F_istat = min(sl_get('bh91_c',bh91_c),F_istat)
      F_istat = min(sl_get('bh91_d',bh91_d),F_istat)
      F_istat = min(sl_get('d97_as',d97_as),F_istat)
      F_istat = min(sl_get('dg92_ci',dg92_ci),F_istat)
      F_istat = min(sl_get('l07_ah',l07_ah),F_istat)
      F_istat = min(sl_get('l07_am',l07_am),F_istat)
      F_istat = min(sl_get('rineutral',sl_rineutral),F_istat)
      if (.not.RMN_IS_OK(F_istat)) &
           call msg(MSG_ERROR,'(sfc_options_init) cannot retrieve AS, BETA, CI or NEUTRAL')
      return
   end function sfc_options_init

end module sfc_options
