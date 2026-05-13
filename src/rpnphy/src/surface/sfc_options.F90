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

   !# ----  FOR TEB BUILT-UP LAND SCHEME ------
   integer, parameter :: NB_TEB_SNOWPAR = 13    !  number of parameters for snow in TEB
   ! values of (1) WCRN,(2)tau_a (cold snow alb),(3)Almin, (4)almax, (5)rhomin, (6)rhomax, (7)tau_r (compaction),
   !           (8)emissivity, (9)tau_f (warm snow alb), (10) tau_anthropo, (11) zom, (12) z0h  (13) snfrac_max
   !                                                           (1)  (2)    (3)   (4)  (5)   (6)   (7)   (8)    (9)  (10)  (11)   (12)  (13)
    real, dimension(NB_TEB_SNOWPAR):: TEB_SNROOF_DEFAULT =  (/ 1.0, 0.04, 0.50, 0.85, 100., 750., 0.01, 1.00, 0.01, 0.0, 0.001, 0.0001, 1.0 /)
    real, dimension(NB_TEB_SNOWPAR):: TEB_SNROAD_DEFAULT =  (/ 1.0, 0.08, 0.20, 0.85, 100., 750., 0.10, 1.00, 0.01, 1.0, 0.001, 0.0001, 0.7 /)
   !
   integer, parameter :: NB_TEB_HYDROPAR = 2    !  number of parameters for hydrology in TEB
   ! values of (1) WS_ROOF_MAX,(2) WS_ROAD_MAX
    real, dimension(NB_TEB_HYDROPAR):: TEB_HYDROPAR_DEFAULT =  (/ 1.0, 1.0 /)
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
   logical           :: thermal_stress_utci = .false.
   logical           :: thermal_stress_shade = .false.
   logical           :: thermal_stress_roof = .false.
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

   !# (SVS2) Options for transfer of the forcing in the canopy module of SVS2
   !# * 'FOR'   :  forcing is below the canopy
   !# * 'O2F'   :  met forcing is below canopy height in the open and is transfered below canopy
   !# * 'ABV'   :  met forcing is above canopy
   character(len=3) :: cano_ref_forcing    = 'ABV'
   namelist /surface_cfgs/ cano_ref_forcing 
   character(len=*), parameter :: CANO_REF_FORCING_OPT(3) = (/ &
        'FOR',  &
        'O2F',  &
        'ABV'   &
        /)

   !# chi_min parameter: the minimal efficiency factor for phase change in the soil_freezing scheme (the same concept as in Surfex and described in Pitman et al., 1991, Slater et al., 1998.
   real              :: chi_min    = 0.6
   namelist /surface_cfgs/ chi_min
   
   !# CSLM: maximum number of lake layers 
   integer           :: cslm_nlakmax   = 20
   namelist /surface_cfgs/ cslm_nlakmax
       
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

   !# (SVS2) Option for the compaction scheme for Crocus
   !# * 'B92' :  default Crocus from Brun et al. 1992 or Vionnet et al. 2012 (Default in SVS2)
   !# * 'S14' : use the settling param of Schleef et al. (2014) for fresh snow (less than 2 days) 
   !# * 'T11' : param of snow viscosity from Teufelsbauer (2011)
   !# * 'R21' :  Reduction in snow compaction due to basal vegetation and no snowdrift below veg height (Royer 2021 & Lackner 2022)   
   !# * 'R2V' : Reduction in snow compaction due to basal vegetation only  (Royer 2021 & Lackner 2022) 
   !# * 'R2D' : No snowdrift below veg height (Royer 2021)
   character(len=16) :: hsnowcomp = 'B92'
   namelist /surface_cfgs/ hsnowcomp
   character(len=*), parameter :: HSNOWCOMP_OPT(6) = (/ &
        'B92',  &
        'S14',  &  
        'T11',  &
        'R21',   &
        'R2V',  &
        'R2D'   &
         /)


   !# (SVS2) Option for the thermal conductivity scheme for Crocus
   !# * 'Y81': default Crocus from Yen et al. 1981 (Default in SVS2)
   !# * 'I02': ISBA_ES snow conductivity parametrization (Boone et al. 2002)
   !# * 'C11': Calonne et al. 2011 snow conductivity parametrization
   !# * 'S97': Sturm et al. 1997
   !# * 'J91': Jordan et al. 1991
   !# * 'F21': Fourteau et al. 2021
   character(len=16) :: hsnowcond = 'Y81'
   namelist /surface_cfgs/ hsnowcond
   character(len=*), parameter :: HSNOWCOND_OPT(6) = (/ &
        'Y81',  &
        'I02',  &  
        'C11',  &
        'S97',  &  
        'J91',  &
        'F21'   &
        /)

   
   !# (SVS2)  Option for the snowdrift scheme for Crocus: Mechanical transformation of snow grain and compaction 
   !#  and effect of wind  on falling snow properties
   !# * 'NONE': No snowdrift scheme
   !# * 'DFLT': falling snow falls as purely dendritic
   !# * 'GA01': Gallee et al 2001
   !# * 'VI13': Vionnet et al 2013 (Default in SVS2)
   !# * 'R21R': Royer et al 2021 (Increase in Maximum Density)
   !# * 'R21W': Royer et al 2021 (Increase in Wind_Effect)
   !# * 'R21F': Royer et al 2021 (Increase in Max Density and Wind Effect)
   character(len=16) :: hsnowdrift_cro = 'VI13'
   namelist /surface_cfgs/ hsnowdrift_cro
   character(len=*), parameter :: HSNOWDRIFT_CRO_OPT(7) = (/ &
        'NONE',  &
        'DFLT',  &
        'GA01',  &  
        'VI13',  &
        'R21F', &
        'R21W', &
        'R21R' &  
        /)
   
   !# (SVS2) Option for the snowdrift scheme for ES: Mechanical transformation of snow grain and compaction 
   !# * 'DFLT': Snowdrift scheme activated  (Default in SVS2)
   !# * 'NONE': No snowdrift scheme
   character(len=16) :: hsnowdrift_es = 'DFLT'
   namelist /surface_cfgs/ hsnowdrift_es
   character(len=*), parameter :: HSNOWDRIFT_ES_OPT(2) = (/ &
        'NONE',  &
        'DFLT'  &
        /)

   !# (SVS2) Option for the falling snow density for Crocus
   !# * 'V12': Vionnet et al. 2012 from Brun et al. 1989 (Default in SVS2)
   !# * 'A76': Anderson et al. 1976
   !# * 'S02': Lehning el al. 2002
   !# * 'P75': Pahaut 1975
   !# * 'NZE': Constant density 200 kg/m3 (defined snowcro.F90 )
   !# * 'R21': Royer et al. 2021
   !# * 'L22': Lackner et al. 2022
   !# * 'GW1': 
   !# * 'GW2':         
   character(len=16) :: hsnowfall = 'V12'
   namelist /surface_cfgs/ hsnowfall
   character(len=*), parameter :: HSNOWFALL_OPT(9) = (/ &
        'V12',  &
        'A76',  &  
        'S02',  &  
        'P75',  &  
        'NZE',  &
        'R21',  &
        'L22',  &  
        'GW1',  &
        'GW2'   &   
         /)
   
   !# (SVS2) Option for the liquid water content scheme for Crocus
   !# * 'B92': default Crocus from Brun et al. 1992 or Vionnet et al. 2012 
   !# * 'B02': ISBA_ES  parametrization (Boone et al. 2002)
   !# * 'CLM': parametrization (Oleson et al 2004)
   !# * 'SPK': SNOWPACK parametrization (Lehning et al 2002)    
   character(len=16) :: hsnowhold = 'B92'
   namelist /surface_cfgs/ hsnowhold
   character(len=*), parameter :: HSNOWHOLD_OPT(4) = (/ &
        'B92',  &
        'B02',  &  
        'O04',  &  
        'SPK'   &   
         /)
 
   !# (SVS2) Option for the metamorphism scheme for Crocus
   !#  * 'B21': Correction of C13 to correctly handle the conversion from 
   !            dendricity/sphericity/grain size to optical diameter/sphericity
   !# * 'T07': Taillandier et al 2007
   !# * 'F06': Flanner et al 2006
   !# * 'S-F': Schlef et al 2014
   !# * 'S-B': Schlef et al 2014
   character(len=16) :: hsnowmetamo = 'B21'
   namelist /surface_cfgs/ hsnowmetamo
   character(len=*), parameter :: HSNOWMETAMO_OPT(5) = (/ &
        'T07',  &  
        'F06',  &  
        'B21',  &
        'S-F',  &
        'S-B'   &  
         /)

   !# (SVS2) Option for the radiative transfer scheme for Crocus
   !# * 'B92':  Brun et al 1992  (Default in SVS-Cro)
   !# * 'T17': (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme
   character(len=16) :: hsnowrad = 'B92'
   namelist /surface_cfgs/ hsnowrad
   character(len=*), parameter :: HSNOWRAD_OPT(2) = (/ &
        'B92',  &
        'T17'   &  
         /)

   !# (SVS2) Option for the turbulent fluxes in Crocus
   !# * 'DEF' : Default: Louis (ISBA: Noilhan and Mahfouf 1996)
   !# * 'RIL' : Limit Richarson number under very stable conditions to 0.2 (currently testing)
   !# * 'RI1' : Limit Richarson number under very stable conditions to 0.1
   !# * 'RI2' : Limit Richarson number under very stable conditions to 0.026
   !# * 'M98' : Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus   
   character(len=16) :: hsnowres = 'RIL'
   namelist /surface_cfgs/ hsnowres
   character(len=*), parameter :: HSNOWRES_OPT(5) = (/ &
        'DEF',  &
        'RIL',  &  
        'RI1',  &  
        'RI2',  &  
        'M98'   &  
         /)


   !# (SVS2) Option for the snowpack representation in SVS-2
   !# * 'ES' : Explicit snow scheme 
   !# * 'CRO': Crocus
   character(len=16) :: hsnowscheme = 'CRO '
   namelist /surface_cfgs/ hsnowscheme
   character(len=*), parameter :: HSNOWSCHEME_OPT(2) = (/ &
        'ES ',  &
        'CRO'  &  
         /)
   
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

   !# Specify which method is used to compute hrsurf
   !# * ALPHA_JN90  : (default) [pending description]
   !# * BETA_ECMWF12: [pending description]
   character(len=16) :: isba_hrsurf_method    = 'ALPHA_JN90'
   namelist /surface_cfgs/ isba_hrsurf_method
   character(len=*), parameter :: ISBA_HRSURF_METHOD_OPT(2) = (/ &
        'ALPHA_JN90   ', &
        'BETA_ECMWF12 ' &
        /)  

   !# Specify the exponent applied to the ratio WD/WFC when computing Beta (only used by BETA_ECMWF12 method) 
   real              :: isba_hrsurf_power = 1.
   namelist /surface_cfgs/ isba_hrsurf_power   

   !# Specify the typical soil resistance for hrsurf computation (only used by BETA_ECMWF12 method) [s/m]
   real              :: isba_hrsurf_rs = 50.
   namelist /surface_cfgs/ isba_hrsurf_rs
   
   !# If .true., field capacity for bare ground evaporation is calculated with respect to 
   !             liquid water using modified porosity in presence of frozen soil
   logical           :: isba_modwsat_ice = .false.
   namelist /surface_cfgs/ isba_modwsat_ice

   !# If .true. apply temporary fix to ISBA
   !# * timestep dependent KCOEF
   !# * No PSN factor for meting and freezing
   logical           :: isba_melting_fix = .false.
   namelist /surface_cfgs/ isba_melting_fix
   
   !# Snow melt/freeze under vegetation impact T2 instead of TST
   logical           :: isba_snow_melt_t2veg = .false.
   namelist /surface_cfgs/ isba_snow_melt_t2veg
   
   !# Use the vegetation-only roughness length to compute vegetation snow fraction
   logical           :: isba_snow_z0veg = .false.
   namelist /surface_cfgs/ isba_snow_z0veg

   !# Surface weights / Poids; LEGACY as of sps_6.3.0.
   character(len=6)           :: sfc_poids = 'LEGACY'
   namelist /surface_cfgs/ sfc_poids
   character(len=*), parameter :: SFC_POIDS_OPT(3) = (/ &
        'LEGACY', &
        'CONSER', &
        'MODMG ' &
        /)

   !# Computation of bare ground snow fraction
   !# * 'NIL'   : Legacy approach implemented in ISBA (Belair et al. 003)
   !# * 'PHY98' : Use the same definition as in Physics 1998 documentation, Douville 1995 and Pitman 1991
   !# * 'SVS1'  : Use the same definition as in SVS1
   character(len=16)           :: isba_snowfrac_bare = 'NIL'
   namelist /surface_cfgs/ isba_snowfrac_bare
   character(len=*), parameter :: ISBA_SNOWFRAC_BARE_OPT(3) = (/ &
        'NIL   ', &
        'PHY98 ', &
        'SVS1  '  &
        /)
   
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

   !# Last/Deepest soil layer where ploughing has an effet
   !# when activating the option "svs_tdrains_plough".
   !# Ploughing generally affects soils down to 20/40 cm
   integer           :: kplough    = -1
   namelist /surface_cfgs/ kplough

   !# Number of the soil layer containing the Tile Drains,
   !# when activating the option "svs_tdrains_plough".
   !# Tile Drains generally located between 90 and 200 cm deep. 
   !# (generally 90/100cm in Canada, sometimes deeper in the US).
   integer           :: ktdrains    = -1
   namelist /surface_cfgs/ ktdrains

   !# Vegetation field update frequency (units D,H,M,S,P)
   character(len=16) :: kntveg_S     = ''
   namelist /surface_cfgs/ kntveg_S

   !# Lead fraction for ice-covered lakes
   real              :: lake_leadfrac = 0.
   namelist /surface_cfgs/ lake_leadfrac

   !# (SVS2) Lower boundary conditions for heat equation in SVS2
   !# * 'TPERM'   :  prescribed temperature at bottom of the soil column 
   !# * '0FLUX'   :  assumes heat flux at the bottom of the soil column is negligible
   character(len=5) :: lbcheat_svs2    = 'TPERM'
   namelist /surface_cfgs/ lbcheat_svs2 
   character(len=*), parameter :: LBCHEAT_SVS2_OPT(2) = (/ &
        'TPERM',  &
        '0FLUX'  &
        /)

   !# Minimum fraction of leads in sea ice.&nbsp; Multiply ice fraction by (1.-leadfrac)
   real              :: leadfrac    = 0.03
   namelist /surface_cfgs/ leadfrac

   !# Limit snow depth to 10 cm for calculation of heat conductivity of snow
   !# over sea-ice and glacier if .true.
   logical           :: limsnodp    = .false.
   namelist /surface_cfgs/ limsnodp

   !# (SVS2) activate forest litter in SVS2 if .true.
   logical           :: lforlit     = .false.
   namelist /surface_cfgs/ lforlit

   !# (SVS1) If .true., macropores are activated and enhance infiltration in a frozen soil
   logical           :: lmacropores_svs = .false.
   namelist /surface_cfgs/ lmacropores_svs

   !#  (SVS1) If .true., we use an efficiency factor phase change (chi) in the soil freezing scheme
   logical           :: lphase_change_eff_svs1 = .true.
   namelist /surface_cfgs/ lphase_change_eff_svs1 

   !#  (SVS1) If .true., field capacity for bare ground evaporation and wilting point for bare ground albedo are 
   !              calculated with respect to liquid water using a modified porosity in presence of frozen soil
   logical           :: lmodwsat_ice_svs1 = .false.
   namelist /surface_cfgs/ lmodwsat_ice_svs1
   
   !# (SVS2) If .true., SVS2 simulates interception of snow by canopy, sublimation and inloading of intercepted snow
   logical           :: lsnow_interception_svs2 = .false.
   namelist /surface_cfgs/ lsnow_interception_svs2

   !# (SVS2) If .true., SVS2 uses the aerodynamic resistance from sfc_layer
   logical           :: lsfclayer_crocus_svs2 = .false.
   namelist /surface_cfgs/ lsfclayer_crocus_svs2


   !# (SVS2) Activate spatially variable snow aging coefficient in Crocus snow albedo parameterization
   logical           :: lsnowaging_var = .false. 
   namelist /surface_cfgs/ lsnowaging_var

   !# (SVS2) Activate mass loss due to blowing snow sublimation in ES and Crocus (default TRUE in SVS2)
   logical           :: lsnowdrift_sublim = .true.
   namelist /surface_cfgs/ lsnowdrift_sublim
   
   !# (SVS2) If .true., SVS2 uses only a unique column for the soil
   logical           :: lunique_profile_svs2 = .true.
   namelist /surface_cfgs/ lunique_profile_svs2 

   !# If .true., SVS1 or SVS2 simulates water ponding at the surface
   logical           :: lwater_ponding_svs = .false.
   namelist /surface_cfgs/ lwater_ponding_svs

   ! mp_alpha is a parameter related to the macropore activation function in SVS 
   !# mp_alpha parameter in the macropore activation option (based on a sensitivity analysis)
   real              :: mp_alpha    = 0.55
   namelist /surface_cfgs/ mp_alpha

   !# (SVS2) Number of snow layers in multi-layer snowpack scheme:
   !# In Crocus it refers to the maximal number of snow layers
   integer           :: nsl    = 12
   namelist /surface_cfgs/ nsl

   !# (coupling) fluxes over ocean are taken from ocean model if .true.
   logical           :: owflux      = .false.
   namelist /surface_cfgs/ owflux
   
   !# (CSLM) read-in light extinction coef. for CSLM  if .true.
   logical           :: read_blak_cslm     = .false.
   namelist /surface_cfgs/ read_blak_cslm
   
   !# read-in land surface emissivity if .true.
   logical           :: read_emis     = .false.
   namelist /surface_cfgs/ read_emis
   
   !# (SVS2) read-in height of polar low vegetation for SVS2 if .true.
   logical           :: read_hveglpol     = .false.
   namelist /surface_cfgs/ read_hveglpol

   !# (SVS2) read-in organic content for SVS2 if .true.
   logical           :: read_oc     = .false.
   namelist /surface_cfgs/ read_oc
   
   !# read-in high vegetation density for SVS2 if .true.
   !logical           :: read_vgh_dens     = .false.
   !namelist /surface_cfgs/ read_vgh_dens
   
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
   !# * 'SVS2' : Advanced version of the SVS land sfc scheme
   character(len=16) :: schmsol     = 'ISBA'
   namelist /surface_cfgs/ schmsol
   character(len=*), parameter :: SCHMSOL_OPT(4) = (/ &
        'NIL ', &
        'ISBA', &
        'SVS ', &
        'SVS2'  &
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


   !# Coefficient for calculation of displacement height from roughness (recommended value=8)
   real              :: sl_afd = -1.
   namelist /surface_cfgs/ sl_afd
   
   !# Type of diagnostic-level calculations to use
   !# * 'FLUX      ' : Use turbulent surface fluxes to estimate diagnostic-level fields
   !# * 'INTERP    ' : Use interpolation to estimate diagnostic-level fields
   character(len=16) :: sl_diag_type = 'FLUX'
   namelist /surface_cfgs/ sl_diag_type
   character(len=*), parameter :: SL_DIAG_TYPE_OPT(2) = (/ &
        'FLUX  ', &
        'INTERP' &
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

   !# Type of Lmin calculation to use
   !# * 'SHEAR     ' : Impose a minimum lowest-level wind speed to guarantee L>L_min
   !# * 'SFO       ' : Impose a minimum L_min directly in stability function calculations
   character(len=16) :: sl_Lmin_type = 'SHEAR'
   namelist /surface_cfgs/ sl_Lmin_type
   character(len=*), parameter :: SL_LMIN_TYPE_OPT(2) = (/ &
        'SHEAR', &
        'SFO  ' &
        /)
   
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

   !# Coefficient for scalar interpolation to the diagnostic level for Rib-based adjustment
   real              :: sl_re2 = -1.
   namelist /surface_cfgs/ sl_re2

   !# Set maximum scale height for surface layer stability functions
   real              :: sl_ximax = -1.
   namelist /surface_cfgs/ sl_ximax
   
   !# Use a reference roughness for surface layer calculations
   logical :: sl_z0ref = .false.
   namelist /surface_cfgs/ sl_z0ref

   !# Emissivity for snow
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: snow_emiss = '1.'
   real              :: snow_emiss_const = -1.
   namelist /surface_cfgs/ snow_emiss

   ! Two options to compute the thermal conductivity of soil (should be used only with soil freezing) 
   ! * PL1998       : Use the model from Peters-Lidard et al. (1998) for frozen soil that involve the Kersten number (Johanssen, 1975) 
   ! * TIAN2016     : Use the physical model from Tian et al. (2016) [https://doi.org/10.1111/ejss.12366]
   character(len=16) :: soil_cond    = 'PL1998'
   namelist /surface_cfgs/ soil_cond
   character(len=*), parameter :: SOIL_COND_OPT(2) = (/ &
        'PL1998  ', &
        'TIAN2016' &
        /)
   
   !#  Soil texture database/calculations for SVS land surface scheme
   !# * 'GSDE   '   : 8 layers of sand & clay info from Global Soil Dataset for ESMs (GSDE)
   !# * 'SLC    '   : 5 layers of sand & clay info from Soil Landscape of Canada (SLC)
   !# * 'SOILGRIDS' : 7 layers of sand & clay info from ISRIC World Soil Information
   !# * 'SOILGRIDSV2' : 6 layers of sand, clay, bulk density, & oc info from ISRIC World Soil Information
   character(len=16) :: soiltext    = 'GSDE'
   namelist /surface_cfgs/ soiltext
   character(len=*), parameter :: SOILTEXT_OPT(4) = (/ &
        'GSDE       ',  &
        'SLC        ',  &
        'SOILGRIDS  ',  &
        'SOILGRIDSV2'   &
        /)

   !# If .true., SVS1 simulates soil freezing and thawing and its impact on hydrology
   logical           :: lsoil_freezing_svs1 = .false.
   namelist /surface_cfgs/ lsoil_freezing_svs1 

   !#  Option to compute the heat flux between the soil and the snowpack in the module soil_freezing.F90
   character(len=16) :: soilsnowhf_svs1     = 'ST_D_DD'
   namelist /surface_cfgs/ soilsnowhf_svs1
   character(len=*), parameter :: SOILSNOWHF_SVS1_OPT(4) = (/ &
        'DST_HD  ', &   ! DST_HD: use the Deep Snow Temperature and Half the snow Depth
        'DST_FD  ', &   ! DST_FD: use the Deep Snow Temperature and the Full snow Depth
        'DST_MAXD', &   ! DST_MAXD (DEFAULT): use the Deep Snow Temp. and the MAX of half the snow Depth and the full snow depth minus the damping depth
        'ST_D_DD '  &   ! ST-D-DD: if snow depth is lower than damping depth, use the full snow depth and the skin Snow Temperature if snow Depth is greater than Damping Depth, use half snow depth and deep snow temp.
        /)

   !#  Option to compute HFLUX for bareground in the module soil_freezing.F90 (Use the soil thermal restiance (default) or the skin conductivity from Table 2.1 of the sup. mat. of Boussetta et al., 2021)
   character(len=16) :: soilgrndhf_svs1     = 'LAM_BOU2021'
   namelist /surface_cfgs/ soilgrndhf_svs1
   character(len=*), parameter :: SOILGRNDHF_SVS1_OPT(2) = (/ &
        'RTH_GRND   ', &   ! RTH_GRND: use a soil thermal resistance in the computation of surface heat flux of bare ground
        'LAM_BOU2021' &   ! LAM_BOU2021: use the skin thermal conductivity from Table 2.1 of the sup. mat. of Boussetta et al., 2021 in the computation of surface heat flux of bare ground
        /)

   !#  Option to set the depth of the lower boundrary condition below the soil column (module soil_freezing.F90)
   character(len=16) :: soildbtm_svs1     = 'DEEP'
   namelist /surface_cfgs/ soildbtm_svs1
   character(len=*), parameter :: SOILDBTM_SVS1_OPT(3) = (/ &
        'MID ', &   ! default condition. The temperature is set at a depth of half the thickness of the lower layer below the bottom of the soil column
        'DEEP', &   ! the temperature is set at a depth of 2.5 times the total thickness of the soil colum with a minimum (minimum 7.5 m and maximum 12.5 m below the surface)
        'NOFL'  &   ! Zero flux condition at the bottom of the soil column
        /)

   ! Three options to compute hydraulic conductivity in the presence of ice 
   ! * NONE       : No modification of hydraulic conductivity in presence of ice 
   ! * ZHANGGRAY97: Correction factor taken from Zhang and Gray (1997). Same as CLASS 3.6
   ! * BOONE2000  : Impedance factor taken from SURFEX (Boone et al., 2000)
   character(len=16) :: soil_ksat_ice    = 'ZHANGGRAY97'
   namelist /surface_cfgs/ soil_ksat_ice
   character(len=*), parameter :: SOIL_KSAT_ICE_OPT(3) = (/ &
        'NONE       ', &
        'ZHANGGRAY97', &
        'BOONE2000  ' &
        /)

   !# Use snow albedo "I6" directly if .true.;
   !# Use snow age "XA" to calculate snow albedo if .false.
   logical           :: snoalb_anl  = .true.
   namelist /surface_cfgs/ snoalb_anl

   !# determine depth over which soil thermal inertia is computed [m]
   !# when set to zero, only 1st layer is used (default option for backward compatibility)
   !# default should eventually be changed to 0.12 based on the fact that cg = sqrt(tau*k/pi)
   !# where k is the thermal diffusivity of the soil [m2/s] and tau = 86400 [s]
   !# if no information is known about the soil, IEC 60853-2 recommends k=5·10-7 m2/s
   !# hence svs_cg_depth = 0.12 m
   real              :: svs_cg_depth = 0.
   namelist /surface_cfgs/ svs_cg_depth

   !# specify which equation to use to take ice into account when computing the heat capacity of bare ground
   !# * BELAIR2003: (default) Adjustment to heat capacity based on total water content (WD+WF)
   !# * BOONE2000:  Adjustment to heat capacity that takes into account the reduced heat capacity of ice vs liquid water
   character(len=16) :: svs_cg_ice    = 'BELAIR2003'
   namelist /surface_cfgs/ svs_cg_ice
   character(len=*), parameter :: SVS_CG_ICE_OPT(2) = (/ &
        'BELAIR2003', &
        'BOONE2000 '  &
        /)   

   !# use dynamic calculation of z0h for bare ground + vegetation  for SVS if .true.
   logical           :: svs_dynamic_z0h     = .false.
   namelist /surface_cfgs/ svs_dynamic_z0h

   !# Specify which method is used to compute hrsurf
   !# * ALPHA_JN90    : (default) "alpha" method proposed by Jacquemin and Noilhan (1990), BLM
   !# * BETA_ECMWF12  : resistance method of Albergel et al. (2012) modified by Fortin et al. (2026)
   !#                   in order to behave like a "beta" method in dry conditions
   character(len=16) :: svs_hrsurf_method    = 'ALPHA_JN90'
   namelist /surface_cfgs/ svs_hrsurf_method
   character(len=*), parameter :: SVS_HRSURF_METHOD_OPT(2) = (/ &
        'ALPHA_JN90   ', &
        'BETA_ECMWF12 ' &
        /)

   !# Specify the exponent applied to the ratio WD/WFC when computing Beta (only used by BETA_ECMWF12 method)
   !# With default value of 1, BETA_ECMWF12 method becomes equivalent to resistance method of Albergel et al. (2012)
   real              :: svs_hrsurf_power = 1.
   namelist /surface_cfgs/ svs_hrsurf_power

   !# Specify the typical soil resistance for hrsurf computation (only used by BETA_ECMWF12 method) [s/m]
   !# Default value of 50 s/m proposed by Albergel et al. (2012)
   real              :: svs_hrsurf_rs = 50.
   namelist /surface_cfgs/ svs_hrsurf_rs

   !------------------------------------------------------------------
   ! NUMERICAL METHOD HYDRO_SVS
   ! Numerical method used to solve Richards' equations in hydro_svs
   ! 'EULER':   Euler, forward, 1st order
   ! 'RK4_IC4': Runge-Kutta 4th order (RK4) (incorrect implementation, but used operationnally in NSRPS IC4)
   ! 'RK4':     Runge-Kutta 4th order (RK4) (improved implementation as part of IC5)
   character(len=16) :: svs_hydro_nummethod = 'RK4_IC4'
   namelist /surface_cfgs/ svs_hydro_nummethod
   character(len=*), parameter :: SVS_HYDRO_NUMMETHOD_OPT(3) = (/ &
        'EULER  ', &
        'RK4_IC4', &
        'RK4    ' &
        /)
        
   !# Option to use average normalized water content rather than average stress to compute transpiration
   logical           :: svs_etr_avg_beta = .false.
   namelist /surface_cfgs/ svs_etr_avg_beta

   !# Proportion of static vs dynamic root distribution as a function of soil moisture
   !# Number betweeen 0 and 1 that reflects the confidence that we have in the prescribed root density profile
   !# 0 = static root distribution given by lookup tables, 1 = fully dynamic root distribution depending on soil moisture
   real           :: svs_etr_fcd_dyn = 0.
   namelist /surface_cfgs/ svs_etr_fcd_dyn

   !# Maximum fraction of wilted roots that can be ignored when computing average stress of vegetation
   real              :: svs_etr_max_roots_ignored = 0.0
   namelist /surface_cfgs/ svs_etr_max_roots_ignored

   !# Use local momentum (no snow) roughness for SVS if .true.
   logical           :: svs_local_z0m     = .false.
   namelist /surface_cfgs/ svs_local_z0m  
   
   !# Option to deal with melt due to rain on snow in SVS snow scheme
   !# * BELAIR03: Original formulation from Belair et al. (2003)
   !# * BELAIR03_DTGEM: Revised formulation of Belair et al. (2003) using the time step of GEM at the time when the parameterisation has been
   !developped (12 min). Recommended for coupled applications with short time steps (< 5 min). 
   !# * NONE: Snow melt due to rain deactivated
   !# * LEONARDINI21:  Assume all the energy brought by rain is used to melt the snowpack (see Leonardini et al. (2021)
   ! Options REV and L21 remvoe the large sensitivity to the time step found with DEF.  
   character(len=16) :: svs_snow_rain   = 'BELAIR03'
   namelist /surface_cfgs/ svs_snow_rain
   character(len=*), parameter :: SVS_SNOW_RAIN_OPT(4) = (/ &
        'BELAIR03         ',  &
        'BELAIR03_DTGEM   ',  &
        'NONE             ',  &
        'LEONARDINI21     ' &
        /)

   !# Computation of snow cover fraction over ground in SVS1
   !# * 'NIL'   : Legacy approach implemented in SVS (Belair et al. 003)
   !# * 'LA23' : Adjusted method from Lalande et al. (2023) accoungint for hysterisis of snow cover fraction and subgrid topo
   !# * 'AR25' : Method from Abolafia-Rosenzweig et al. (2025) applied to resolution of 2.5 km
   character(len=16)           :: svs_snowfrac_ground = 'NIL'
   namelist /surface_cfgs/ svs_snowfrac_ground
   character(len=*), parameter :: SVS_SNOWFRAC_GROUND_OPT(3) = (/ &
        'NIL   ', &
        'LA23  ', &
        'AR25  '  &
        /)

   ! Parameters used in the param. of snow cover fraction LA23 when VF_TYPE is set to CCILC_WE
   !# z0min_la23 is the minimal roughness parameter used in LA23 in presence of low vegetation of roughness (Z0VL) below 0.01 m
   !# z0max_la23 is the maximal roughness parameter used in LA23 in presence of low vegetation of roughness (Z0VL) below 0.15 m
   !# z0exp_la23 is the exponent used in the power law giving the roughness parameter as a function of Z0VL
   real              :: z0min_la23    = 0.01
   namelist /surface_cfgs/ z0min_la23   
   real              :: z0max_la23    = 0.02
   namelist /surface_cfgs/ z0max_la23   
   real              :: z0exp_la23    = 6.
   namelist /surface_cfgs/ z0exp_la23   

   !# Computation of snow albedo in SVS1
   !# * 'NIL'   : Legacy approach implemented in SVS (Leonardini et al., 2021)
   !# * 'CLASS' : Method from the land surface scheme CLASS
   character(len=16)           :: svs_snowalb = 'NIL'
   namelist /surface_cfgs/ svs_snowalb
   character(len=*), parameter :: SVS_SNOWALB_OPT(2) = (/ &
        'NIL     ' , &
        'CLASS   ' &
        /)

   !# specify functions to use to compute wsat, wfc and wwilt
   !# SURFEX: equations used in SURFEX v8.1, obtained from Giordani (1993) and Noilhan & Lacarrère (1995)
   !# USDA2006: equations proposed by USDA, see Saxton and Rawls (2006)
   character(len=16) :: svs_soiltext2prop = 'SURFEXV8'
   namelist /surface_cfgs/ svs_soiltext2prop
   character(len=*), parameter :: SVS_SOILTEXT2PROP_OPT(2) = (/ &
        'SURFEXV8     ', &
        'USDA2006     ' &
        /)

   !# Option to activate the effect of ploughing and Tile Drains in SVS
   logical           :: svs_tdrains_plough = .false.
   namelist /surface_cfgs/ svs_tdrains_plough

   !# Option to change the lookup table ALDAT in inicover_svs.F90
   logical :: svs_read_aldat = .false.
   namelist /surface_cfgs/ svs_read_aldat

   !# Values (1:NCLASS) used for the lookup table ALDAT if svs_read_aldat=.T.
   real    :: svs_aldat(NCLASS) = -999.
   namelist /surface_cfgs/ svs_aldat

   !# Option to change the lookup table D2DAT in inicover_svs.F90
   logical :: svs_read_d2dat = .false.
   namelist /surface_cfgs/ svs_read_d2dat

   !# Values (1:NCLASS) used for the lookup table D2DAT if svs_read_d2dat=.T.
   real    :: svs_d2dat(NCLASS) = -999.
   namelist /surface_cfgs/ svs_d2dat

   !# Option to change the lookup table D50DAT in inicover_svs.F90
   logical :: svs_read_d50dat = .false.
   namelist /surface_cfgs/ svs_read_d50dat

   !# Values (1:NCLASS) used for the lookup table D50DAT if svs_read_d50dat=.T.
   real    :: svs_d50dat(NCLASS) = -999.
   namelist /surface_cfgs/ svs_d50dat

   !# Option to change the lookup table D95DAT in inicover_svs.F90
   logical :: svs_read_d95dat = .false.
   namelist /surface_cfgs/ svs_read_d95dat

   !# Values (1:NCLASS) used for the lookup table D95DAT if svs_read_d95dat=.T.
   real    :: svs_d95dat(NCLASS) = -999.
   namelist /surface_cfgs/ svs_d95dat

   !# Option to change the lookup table VEGDAT in inicover_svs.F90
   logical :: svs_read_vegdat = .false.
   namelist /surface_cfgs/ svs_read_vegdat

   !# Values (1:NCLASS) used for the lookup table VEGDAT if svs_read_vegdat=.T.
   real    :: svs_vegdat(NCLASS) = -999.
   namelist /surface_cfgs/ svs_vegdat

   !# Option to change the lookuptable Z0MDAT in inicover_svs.F90
   logical :: svs_read_z0mdat = .false.
   namelist /surface_cfgs/ svs_read_z0mdat

   !# Values (1:NCLASS) used for the lookup table Z0MDAT if svs_read_z0mdat=.T.
   real    :: svs_z0mdat(NCLASS) = -999.
   namelist /surface_cfgs/ svs_z0mdat

   !# Option to change the lookup table VF2CTEMDAT in phtsyn_svs.F90
   logical :: svs_read_vf2ctemdat = .false.
   namelist /surface_cfgs/ svs_read_vf2ctemdat

   !# Values (1:NCLASS) used for the lookup table VF2CTEMDAT if svs_read_vf2ctemdat=.T.
   integer :: svs_vf2ctemdat(NCLASS) = -999
   namelist /surface_cfgs/ svs_vf2ctemdat
   
   
   
   !# Option to change the lookup table VEGCROPS in inicover_svs.F90
   logical :: svs_read_vegcrops = .false.
   namelist /surface_cfgs/ svs_read_vegcrops

   !# Values (1:13) used for the lookup table VEGCROPS if svs_read_vegcrops=.T.
   real    :: svs_vegcrops(13) = -999.
   namelist /surface_cfgs/ svs_vegcrops

   !# Option to read a lookup table of monthly values for D50 for class 15
   logical :: svs_read_d50veg15 = .false.
   namelist /surface_cfgs/ svs_read_d50veg15

   !# Option to read a lookup table of monthly values for D50 for class 16
   logical :: svs_read_d50veg16 = .false.
   namelist /surface_cfgs/ svs_read_d50veg16

   !# Values (1:13) used for the lookup table D50VEG15 if svs_read_d50veg15=.T.
   real    :: svs_d50veg15(13) = -999.
   namelist /surface_cfgs/ svs_d50veg15

   !# Values (1:13) used for the lookup table D50VEG16 if svs_read_d50veg16=.T.
   real    :: svs_d50veg16(13) = -999.
   namelist /surface_cfgs/ svs_d50veg16

   !# Option to read a lookup table of monthly values for D95 for class 15 
   logical :: svs_read_d95veg15 = .false.
   namelist /surface_cfgs/ svs_read_d95veg15

   !# Option to read a lookup table of monthly values for D95 for class 16 
   logical :: svs_read_d95veg16 = .false.
   namelist /surface_cfgs/ svs_read_d95veg16

   !# Values (1:13) used for the lookup table D95VEG15 if svs_read_d95veg15=.T.
   real    :: svs_d95veg15(13) = -999.
   namelist /surface_cfgs/ svs_d95veg15

   !# Values (1:13) used for the lookup table D95VEG16 if svs_read_d95veg16=.T.
   real    :: svs_d95veg16(13) = -999.
   namelist /surface_cfgs/ svs_d95veg16


   !# Option to read a lookup table of monthly values for D50 for class 16
   logical :: svs_read_gexpveg = .false.
   namelist /surface_cfgs/ svs_read_gexpveg

   !# Exponent in function defining vegetation stress when estimating transpiration
   !# Transpiration decreases more slowly with soil moisture when svs_gexp is high
   !# it does not start decreasing until about soil moisture is half way between wilting
   !# point and field capacity when svs_gexp=10.
   !# Prior to introducing this key, a value of two was used in phtsyn_svs.F90 but a value
   !# of one was assumed in vegi_svs.F90, leading to an inconsistency in the code when
   !# the CTEM parameterization is used.
   !# This behaviour is preserved if svs_read_gexpveg = .false. for
   !# backward compatibility purposes.
   !# Since this is a bugfix, eventually the default value of the key should be changed
   !# to a positive value.
   !# Values (26) used for the lookup table D95VEG16 if svs_read_gexpveg=.T.
   real    :: svs_gexpveg(26) = -999.
   namelist /surface_cfgs/ svs_gexpveg

   !# Option to read a lookup table of monthly values for LAI for class 11
   logical :: svs_read_lai11 = .false.
   namelist /surface_cfgs/ svs_read_lai11

   !# Values (1:13) used for the lookup table LAI11 if svs_read_lai11=.T.
   real    :: svs_lai11(13) = -999.
   namelist /surface_cfgs/ svs_lai11
   
   !# Option to read a lookup table of monthly values for LAI for class 14
   logical :: svs_read_lai14 = .false.
   namelist /surface_cfgs/ svs_read_lai14

   !# Values (1:13) used for the lookup table LAI14 if svs_read_lai14=.T.
   real    :: svs_lai14(13) = -999.
   namelist /surface_cfgs/ svs_lai14
   
   !# Option to read a lookup table of monthly values for LAI for class 15
   logical :: svs_read_lai15 = .false.
   namelist /surface_cfgs/ svs_read_lai15

   !# Values (1:13) used for the lookup table LAI15 if svs_read_lai15=.T.
   real    :: svs_lai15(13) = -999.
   namelist /surface_cfgs/ svs_lai15
   
   !# Option to read a lookup table of monthly values for LAI for class 16
   logical :: svs_read_lai16 = .false.
   namelist /surface_cfgs/ svs_read_lai16

   !# Values (1:13) used for the lookup table LAI16 if svs_read_lai16=.T.
   real    :: svs_lai16(13) = -999.
   namelist /surface_cfgs/ svs_lai16
   
   !# Option to read a lookup table of monthly values for LAI for class 16
   logical :: svs_read_lai17 = .false.
   namelist /surface_cfgs/ svs_read_lai17

   !# Values (1:13) used for the lookup table LAI17 if svs_read_lai17=.T.
   real    :: svs_lai17(13) = -999.
   namelist /surface_cfgs/ svs_lai17

   !# Option to read a lookup table of monthly values for VEGDAT for class 14
   logical :: svs_read_vegdat14 = .false.
   namelist /surface_cfgs/ svs_read_vegdat14

   !# Values (1:13) used for the lookup table VEGDAT(13) if svs_read_vegdat13=.T.
   real    :: svs_vegdat14(13) = -999.
   namelist /surface_cfgs/ svs_vegdat14

   !# Option to read a lookup table of monthly values for VEGDAT for class 16
   logical :: svs_read_vegdat16 = .false.
   namelist /surface_cfgs/ svs_read_vegdat16

   !# Values (1:13) used for the lookup table VEGDAT(16) if svs_read_vegdat16=.T.
   real    :: svs_vegdat16(13) = -999.
   namelist /surface_cfgs/ svs_vegdat16

   !# Option to read a lookup table of monthly values for VEGDAT for class 17
   logical :: svs_read_vegdat17 = .false.
   namelist /surface_cfgs/ svs_read_vegdat17

   !# Values (1:13) used for the lookup table VEGDAT(17) if svs_read_vegdat17=.T.
   real    :: svs_vegdat17(13) = -999.
   namelist /surface_cfgs/ svs_vegdat17

   !# Logical - Activate Crocus debugging options  if .true. (see README_crodebug.md)
   logical           :: svs2_crodebug   = .false.
   namelist /surface_cfgs/ svs2_crodebug
   
   !# Limit temperature inversions to dlrate in surface layer if .true.
   logical           :: tdiaglim    = .false.
   namelist /surface_cfgs/ tdiaglim

   !# Set tdiaglim lapse rate (K/m) in surface layer if .true.
   real           :: tdlrate    = 0.2
   namelist /surface_cfgs/ tdlrate   

   !# OPTION FOR CALCULATION of AVERAGE LAND SURFACE TEMPERATURE AND HUMIDITY IN SVS
   !# * .FALSE. :  Area-average only calculation for sfc T and Hum.
   !# * .TRUE.  :  Option that uses effective surface temperature  and specific humidity instead
   !               of composite (area-averaged only) counterparts in surface flux calculations (D. Deacu)
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
   character(len=*), parameter :: VFTYPE_OPT(3) = (/ &
        'CLASSIC   ',  &
        'CCILCECO  ',  &
        'CCILC_WE  '   &
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

   !# (SVS2) Methodology to compute roughness length over snow in SVS2
   !# *  'SURFEX'  : Same method as in SURFEX/Crocus using : value for snow only 
   !# *  'SVS1  '  : Same method as in SVS1: grid-averaged value for momentum, local value for snow for heat
   !# *  'HYB   '  : Hybrid method:  grid-averaged value for momentum (as in SVS1), local value for snow from SURFEX for heat 
   character(len=6) :: z0snow_svs2 = 'SURFEX'
   character(len=*), parameter :: Z0SNOW_SVS2_OPT(3) = (/ &
        'SURFEX',  &
        'SVS1  ',  &
        'HYB   '   &
        /)
   namelist /surface_cfgs/ z0snow_svs2
   
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

   !# (TEB) Options for parameters governing the snow evolution on roofs nd roads
   !# * 'MA00'   :  Default Masson 2000
   !# * 'JA14'   :    opt Järvi et al 2014
   !# * 'LE10'   :    opt Lemonsu et al 2010
   character(len=6) :: teb_snow    = 'MA00  '
   namelist /surface_cfgs/ teb_snow
   character(len=*), parameter :: TEB_SNOW_OPT(4) = (/ &
        'MA00  ',  &
        'JA14  ',  &
        'LE10  ',  &
        'CUSTOM' &
        /)

  !# read the TEB snow parameters (if schmsol=TEB and teb_snow=custom)
   real :: teb_snroof(NB_TEB_SNOWPAR) = -1.0
   real :: teb_snroad(NB_TEB_SNOWPAR) = -1.0
   namelist /surface_cfgs/ teb_snroof
   namelist /surface_cfgs/ teb_snroad

  !# (TEB) Options for parameters governing the hydrology evolution on roofs nd roads
   !# * 'MA00'   :  Default Masson 2000
   !# * 'CUSTOM'   :    opt for customization 
   character(len=6) :: teb_hydro    = 'MA00  '
   namelist /surface_cfgs/ teb_hydro
   character(len=*), parameter :: TEB_HYDRO_OPT(2) = (/ &
        'MA00  ',  &
        'CUSTOM' &
        /)

  !# read the TEB hydro parameters (if schmsol=TEB and teb_hydro=custom)
   real :: teb_hydropar(NB_TEB_HYDROPAR) = -1.0
   namelist /surface_cfgs/ teb_hydropar


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
