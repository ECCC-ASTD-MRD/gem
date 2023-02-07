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

module phy_options
   use, intrinsic :: iso_fortran_env, only: INT64
   implicit none
   public
   save

   integer, parameter :: OPT_OPTIX_OLD     = 1
   integer, parameter :: OPT_OPTIX_NEW     = 2
   integer, parameter :: RAD_NUVBRANDS     = 6 !#TODO: move to a radiation specific module/cdk
   integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: move to mu_jdate
   logical           :: chemistry    = .false.
   logical           :: climat       = .false.
   logical           :: cmt_comp_diag = .false.
   character(len=16) :: conv_mid     = 'NIL'
   character(len=16) :: conv_shal    = 'NIL'
   character(len=16) :: convec       = 'NIL'
   integer(INT64)    :: jdateo       = 0
   real              :: delt         = 0.
   character(len=4),  pointer :: dyninread_list_s(:) => NULL()
   logical           :: dynout       = .false.
   logical           :: ebdiag       = .false.
   logical           :: ecdiag       = .false.
   logical           :: etccdiag     = .false.
   logical           :: impflx       = .false.
   logical           :: inincr       = .false.
   integer           :: ioptix       = OPT_OPTIX_OLD
   integer           :: kntrad       = 1
   integer           :: kntraduv     = -1
   logical           :: llight       = .false.
   logical           :: llinoz       = .false.
   logical           :: llingh       = .false.
   integer           :: lhn_ramp     = 10
   integer           :: lhn_start    = 10
   integer           :: lhn_stop     = 360
   logical           :: lrefract     = .false.
   logical           :: out_linoz    = .false.
   logical           :: age_linoz    = .false.
   character(len=1024) :: ozone_file_s = 'NIL'
   character(len=16) :: schmsol      = 'ISBA'
   logical           :: slt_winds    = .false.
   logical           :: tdiaglim     = .false.
   integer           :: tlift        = 0
   integer           :: nphyoutlist  = 0
   integer           :: ilongmel     = -1
   character(len=32), pointer :: phyoutlist_S(:) => NULL()
   character(len=32) :: vgrid_M_S = 'ref-m'
   character(len=32) :: vgrid_T_S = 'ref-t'

   !# Time length (hours) for special time accumulated physics variables
   integer           :: acchr        = 0
   namelist /physics_cfgs/ acchr
   namelist /physics_cfgs_p/ acchr

   !# Turbulent kinetic energy advect. is active if .true.
   logical           :: advectke     = .false.
   namelist /physics_cfgs/ advectke
   namelist /physics_cfgs_p/ advectke

   !# Clip tracers negative values
   logical           :: clip_tr_L      = .true.
   namelist /physics_cfgs/ clip_tr_L
   namelist /physics_cfgs_p/ clip_tr_L

   !# Conservation corrections for gridscale condensation
   !# * 'NIL ' : No conservation correction applied
   !# * 'TEND' : Temperature and moisture tendencies corrected
   character(len=16) :: cond_conserve = 'NIL'
   namelist /physics_cfgs/ cond_conserve
   namelist /physics_cfgs_p/ cond_conserve
   character(len=*), parameter :: COND_CONSERVE_OPT(2) = (/ &
        'NIL ', &
        'TEND'  &
        /)

   !# Evaporation parameter for Sunqvist gridscale condensation
   real           :: cond_evap       = 2.e-4
   namelist /physics_cfgs/ cond_evap
   namelist /physics_cfgs_p/ cond_evap

   !# Dry accretion of ice crystals by snow in Sundqvist (factor)
   real           :: cond_iceacc    = 5.
   namelist /physics_cfgs/ cond_iceacc
   namelist /physics_cfgs_p/ cond_iceacc
   
   !# Minimum cloud mixing ratio (kg/kg) for autoconversion in
   !# Sunqvist gridscale condensation
   real           :: cond_hmrst     = 3.e-4
   namelist /physics_cfgs/ cond_hmrst
   namelist /physics_cfgs_p/ cond_hmrst

   !# Min allowed values of modified hu00 (threshold relative humidity
   !# for stratiform condensation, Sunqvist gridscale condensation)
   real           :: cond_hu0min       = 0.85
   namelist /physics_cfgs/ cond_hu0min
   namelist /physics_cfgs_p/ cond_hu0min

   !# Max allowed values of modified hu00 (threshold relative humidity
   !# for stratiform condensation, Sunqvist gridscale condensation)
   real           :: cond_hu0max       = 0.975
   namelist /physics_cfgs/ cond_hu0max
   namelist /physics_cfgs_p/ cond_hu0max
   
   !# Activate computing of all diags, requested for output or not.
   logical           :: debug_alldiag_L     = .false.
   namelist /physics_cfgs/ debug_alldiag_L
   namelist /physics_cfgs_p/ debug_alldiag_L
   
   !# Run only the physics nml+init (skip input and step)
   logical           :: debug_initonly_L     = .false.
   namelist /physics_cfgs/ debug_initonly_L
   namelist /physics_cfgs_p/ debug_initonly_L

   !# Activate Debug memory mode
   logical           :: debug_mem_L      = .false.
   namelist /physics_cfgs/ debug_mem_L
   namelist /physics_cfgs_p/ debug_mem_L

   !# Print a trace of the phy functions (MSG verbosity = debug)
   logical           :: debug_trace_L      = .false.
   namelist /physics_cfgs/ debug_trace_L
   namelist /physics_cfgs_p/ debug_trace_L

   !# Diffuse vertical motion if .true.
   logical           :: diffuw       = .false.
   namelist /physics_cfgs/ diffuw
   namelist /physics_cfgs_p/ diffuw

   !# Minimal value for TKE in stable case (for 'CLEF')
   real              :: etrmin2      = 1.E-4
   namelist /physics_cfgs/ etrmin2
   namelist /physics_cfgs_p/ etrmin2

   !# Boundary layer processes
   !# * 'NIL    ': no vertical diffusion
   !# * 'CLEF   ': non-cloudy boundary layer formulation
   !# * 'MOISTKE': cloudy boundary layer formulation
   !# * 'SURFACE': TODO
   !# * 'SIMPLE ': a very simple mixing scheme for neutral PBLs
   !# * 'YSU    ': Yonsei University PBL scheme (from WRF 4.2.1)
   character(len=16) :: fluvert      = 'NIL'
   namelist /physics_cfgs/ fluvert
   namelist /physics_cfgs_p/ fluvert
   character(len=*), parameter :: FLUVERT_OPT(6) = (/ &
        'NIL    ', &
        'CLEF   ', &
        'MOISTKE', &
        'SURFACE', &
        'SIMPLE ', &
        'YSU    ' &
        /)

   !# (MOISTKE only) Apply factor fnn_reduc
   !# * .false.: everywhere
   !# * .true.: over water only
   logical           :: fnn_mask     = .false.
   namelist /physics_cfgs/ fnn_mask
   namelist /physics_cfgs_p/ fnn_mask

   !# (MOISTKE only) Reduction factor (between 0. and 1.) to be applied to the
   !# parameter FNN (turbulent flux enhancement due to boundary layer clouds)
   real              :: fnn_reduc    = 1.
   namelist /physics_cfgs/ fnn_reduc
   namelist /physics_cfgs_p/ fnn_reduc

   !# (CLEF+CONRES only) Non-dimensional parameter (must be >= 1.) that controls
   !# the value of the flux enhancement factor in CONRES
   real              :: fnnmod       = 2.
   namelist /physics_cfgs/ fnnmod
   namelist /physics_cfgs_p/ fnnmod

   !# Gravity wave drag formulation
   !# * 'NIL  ': no Gravity wave drag
   !# * 'SGO16': gravity wave drag + low-level blocking (new formulation 2016)
   character(len=16) :: gwdrag       = 'NIL'
   namelist /physics_cfgs/ gwdrag
   namelist /physics_cfgs_p/ gwdrag
   character(len=*), parameter :: GWDRAG_OPT(2) = (/ &
        'NIL  ', &
        'SGO16'  &
        /)

   !# Number of times the 3-point filter will be applied to smooth the GW flux profiles
   integer           :: hines_flux_filter = 0
   namelist /physics_cfgs/ hines_flux_filter
   namelist /physics_cfgs_p/ hines_flux_filter

   !# Consider heating from non-orog. drag if = 1
   integer           :: iheatcal     = 0
   namelist /physics_cfgs/ iheatcal
   namelist /physics_cfgs_p/ iheatcal

   !# Comma-separated list of diagnostic level inputs to read.
   !# Default: indiag_list_s(1) = 'DEFAULT LIST',
   !# expanded to: UU, VV, TT, HU + all dynamic Tracers
   character(len=32) :: indiag_list_s(128) = ' '
   namelist /physics_cfgs/ indiag_list_s

   !# Initialize water content and cloud fraction seen by radiation for time 0 if .true.
   logical           :: inilwc       = .false.
   namelist /physics_cfgs/ inilwc
   namelist /physics_cfgs_p/ inilwc

   !# Type of input system used
   !# * 'DIST   ' : GEM 5.0 input system, RPN_COMM_IO/RPN_COMM_ezshuf_dist based
   !# * 'BLOC   ' : GEM 5.0 input system, RPN_COMM_bloc based
   !# * 'GEM_4.8' : GEM 4.8 input system, RPN_COMM_bloc based
   character(len=16) :: input_type = 'DIST'
   namelist /physics_cfgs/ input_type
   namelist /physics_cfgs_p/ input_type
   character(len=*), parameter :: INPUT_TYPE_OPT(3) = (/ &
        'DIST   ', &
        'BLOC   ', &
        'GEM_4.8'  &
        /)

   !# Update ozone climatology during the run
   logical           :: intozot      = .false.
   namelist /physics_cfgs/ intozot
   namelist /physics_cfgs_p/ intozot

   !# Time between full radiation calculation (units D,H,M,S,P)
   character(len=16) :: kntrad_S     = ''
   namelist /physics_cfgs/ kntrad_S
   namelist /physics_cfgs_p/ kntrad_S

   !# Time between UV radiation calculation (units D,H,M,S,P)
   character(len=16) :: kntraduv_S     = ''
   namelist /physics_cfgs/ kntraduv_S
   namelist /physics_cfgs_p/ kntraduv_S
   
   !# Add methane oxydation as source of humidity in the stratosphere if .true.
   logical           :: lmetox       = .false.
   namelist /physics_cfgs/ lmetox
   namelist /physics_cfgs_p/ lmetox

   !# Mixing length calc. scheme
   !# * 'BLAC62  ': mixing length calc. using Blackadar
   !# * 'BOUJO   ': mixing length calc. using Bougeault
   !# * 'TURBOUJO': mixing length calc. using Bougeault in turbulent regimes (otherwise Blackadar)
   !# * 'LH      ': mixing length calc. using Lenderink and Holtslag
   character(len=16) :: longmel      = 'BLAC62'
   namelist /physics_cfgs/ longmel
   namelist /physics_cfgs_p/ longmel

   !# Time length (hours) for special time averaged physics variables
   integer           :: moyhr = 0
   namelist /physics_cfgs/ moyhr
   namelist /physics_cfgs_p/ moyhr

   !# Number of ice-phase hydrometeor categories to use in the P3 microphysics
   !# scheme (currently limited to <5)
   integer           :: p3_ncat = 1
   namelist /physics_cfgs/ p3_ncat
   namelist /physics_cfgs_p/ p3_ncat

   !# Maximum time step (s) to be taken by the microphysics (P3) scheme, with time-splitting
   !# used to reduce step to below this value if necessary
   real           :: p3_dtmax = 60.
   namelist /physics_cfgs/ p3_dtmax
   namelist /physics_cfgs_p/ p3_dtmax

   !# calibration factor for ice deposition in microphysics (P3)
   real           :: p3_depfact = 1.0
   namelist /physics_cfgs/ p3_depfact
   namelist /physics_cfgs_p/ p3_depfact

   !# calibration factor for ice sublimation in microphysics (P3)
   real           :: p3_subfact = 1.0
   namelist /physics_cfgs/ p3_subfact
   namelist /physics_cfgs_p/ p3_subfact

   !# switch for real-time debugging in microphysics (P3)
   logical         :: p3_debug = .false.
   namelist /physics_cfgs/ p3_debug
   namelist /physics_cfgs_p/ p3_debug

   !# switch for subgrid cloud/precipitation fraction scheme (SCPF) in microphysics (P3)
   logical         :: p3_scpf_on = .false.
   namelist /physics_cfgs/ p3_scpf_on
   namelist /physics_cfgs_p/ p3_scpf_on

   !# precipitation fraction factor used by SCPF in microphysics (P3)
   real           :: p3_pfrac = 1.0
   namelist /physics_cfgs/ p3_pfrac
   namelist /physics_cfgs_p/ p3_pfrac

   !# model resolution factor used by SCPF in microphysics (P3)
   real           :: p3_resfact = 1.0
   namelist /physics_cfgs/ p3_resfact
   namelist /physics_cfgs_p/ p3_resfact

   !# Switch for aerosol activation scheme (1 = default, 2 = ARG + Aerosol climatology)
   integer           :: mp_aeroact = 1
   namelist /physics_cfgs/ mp_aeroact
   namelist /physics_cfgs_p/ mp_aeroact

   !# Switch for airmass type (1 = maritime, 2 = continental)
   integer           :: my_ccntype   = 1
   namelist /physics_cfgs/ my_ccntype
   namelist /physics_cfgs_p/ my_ccntype

   !# Compute MY Diagnostic fields if .true.
   logical           :: my_diagon    = .true.
   namelist /physics_cfgs/ my_diagon
   namelist /physics_cfgs_p/ my_diagon

   !# Ice-phase switched on if .true.
   logical           :: my_iceon     = .true.
   namelist /physics_cfgs/ my_iceon
   namelist /physics_cfgs_p/ my_iceon

   !# Autoconversion (cloud to rain) switched on
   logical           :: my_rainon    = .true.
   namelist /physics_cfgs/ my_rainon
   namelist /physics_cfgs_p/ my_rainon

   !# Sedimentation switched on
   logical           :: my_sedion    = .true.
   namelist /physics_cfgs/ my_sedion
   namelist /physics_cfgs_p/ my_sedion

   !# Snow initiation switched on
   logical           :: my_snowon    = .true.
   namelist /physics_cfgs/ my_snowon
   namelist /physics_cfgs_p/ my_snowon

   !# Warm-phase switched on
   logical           :: my_warmon    = .true.
   namelist /physics_cfgs/ my_warmon
   namelist /physics_cfgs_p/ my_warmon

   !# Physic input blocking along X
   integer           :: ninblocx     = 1
   namelist /physics_cfgs/ ninblocx
   namelist /physics_cfgs_p/ ninblocx

   !# Physic input blocking along Y
   integer           :: ninblocy     = 1
   namelist /physics_cfgs/ ninblocy
   namelist /physics_cfgs_p/ ninblocy

   !# Hines non-orographic GWD scheme is active if .true.
   logical           :: non_oro      = .false.
   namelist /physics_cfgs/ non_oro
   namelist /physics_cfgs_p/ non_oro

   !# Pressure (in Pa) that defines the bottom emission level for gravity waves
   real              :: non_oro_pbot = 61000.0
   namelist /physics_cfgs/ non_oro_pbot
   namelist /physics_cfgs_p/ non_oro_pbot

   !# Number of timesteps for which surface fluxes "FC" and "FV" are
   !# gradually set from 0 to their full value in a "slow start fashion"
   !# at the beginning of a time integration
   integer           :: nsloflux     = 0
   namelist /physics_cfgs/ nsloflux
   namelist /physics_cfgs_p/ nsloflux
   
   !# Vectoc lenght physics memory space folding for openMP
   integer           :: p_runlgt     = -1
   namelist /physics_cfgs/ p_runlgt
   namelist /physics_cfgs_p/ p_runlgt

   !# Time-averaging of transfer coefficient for momentum to reduce 2-dt 
   !# oscillations in fluxes
   logical           :: pbl_cmu_timeavg = .false.
   namelist /physics_cfgs/ pbl_cmu_timeavg
   namelist /physics_cfgs_p/ pbl_cmu_timeavg
   
   !# Diffuse condensate fields
   logical           :: pbl_diff_condens = .false.
   namelist /physics_cfgs/ pbl_diff_condens
   namelist /physics_cfgs_p/ pbl_diff_condens

   !# Run with a modified closure for the dissipation length scale
   !# * 'NIL  ' : No modified closure for the dissipation length scale
   !# * 'LIM50' : A maximum value of 50m is imposed on dissipation length
   character(len=16) :: pbl_diss     = 'NIL'
   namelist /physics_cfgs/ pbl_diss
   namelist /physics_cfgs_p/ pbl_diss
   character(len=*), parameter :: PBL_DISS_OPT(2) = (/ &
        'LIM50', &
        'NIL  '  &
        /)

   !# Dissipative heating tendencies are computed for the PBL scheme such
   !# that total energy (kinetic + internal) is conserved
   !# * 'NIL       ' : No dissipative heating is computed
   !# * 'LOCAL_K   ' : Local total energy conservation based on diffusion coefficients
   !# * 'LOCAL_TEND' : Local total energy conservation based on wind tendencies
   character(len=16) :: pbl_dissheat = 'NIL'
   namelist /physics_cfgs/ pbl_dissheat
   namelist /physics_cfgs_p/ pbl_dissheat
   character(len=*), parameter :: PBL_DISSHEAT_OPT(3) = (/ &
        'LOCAL_TEND', &
        'LOCAL_K   ', &
        'NIL       '  &
        /)

   !# Conservation corrections for PBL scheme
   !# * 'NIL ' : No conservation correction applied
   !# * 'TEND' : Temperature and moisture tendencies corrected
   character(len=16) :: pbl_conserve = 'NIL'
   namelist /physics_cfgs/ pbl_conserve
   namelist /physics_cfgs_p/ pbl_conserve
   character(len=*), parameter :: PBL_CONSERVE_OPT(2) = (/ &
        'NIL ', &
        'TEND'  &
        /)

   !# Include the turbulent effects of trade wind cumulus clouds
   logical           :: pbl_cucloud  = .true.
   namelist /physics_cfgs/ pbl_cucloud
   namelist /physics_cfgs_p/ pbl_cucloud
 
   !# Class of stability functions (stable case) to use in the PBL
   !# * 'DELAGE97  ' : Use functions described by Delage (1997; BLM)
   !# * 'BELJAARS91' : Use functions described by Beljaars and Holtslag (1991; JAM)
   !# * 'LOCK07    ' : Use functions described by Lock (2007; Tech Report) employed at UKMO
   character(len=16) :: pbl_func_stab = 'DELAGE97'
   namelist /physics_cfgs/ pbl_func_stab
   namelist /physics_cfgs_p/ pbl_func_stab
   character(len=*), parameter :: PBL_FUNC_STAB_OPT(3) = (/ &
        'DELAGE97  ', &
        'BELJAARS91', &
        'LOCK07    '  &
        /)

   !# Class of stability functions (unstable case) to use in the PBL
   !# * 'DELAGE92' : Use functions described by Delage and Girard (1992; BLM)
   !# * 'DYER74  ' : Use functions described by Dyer (1974; BLM)
   character(len=16) :: pbl_func_unstab = 'DELAGE92'
   namelist /physics_cfgs/ pbl_func_unstab
   namelist /physics_cfgs_p/ pbl_func_unstab
   character(len=*), parameter :: PBL_FUNC_UNSTAB_OPT(2) = (/ &
        'DELAGE92', &
        'DYER74  '  &
        /)

   !# Choose form of asymptotic mixing length for Blacadar-type estimates
   !# * 'BLAC62' : Asymptotic 200 m proposed by Blackadar (1962; JGR) with clipping
   !# * 'LOCK07' : Diagnosed asymptotic scale of Lock (2007; Tech Report) used at UKMO
   character(len=16) :: pbl_mlblac_max = 'BLAC62'
   namelist /physics_cfgs/ pbl_mlblac_max
   namelist /physics_cfgs_p/ pbl_mlblac_max
   character(len=*), parameter :: PBL_MLBLAC_MAX_OPT(2) = (/ &
        'BLAC62', &
        'LOCK07'  &
        /)

   !# Apply "turboujo" turbulence conditions to dissipation length scale
   logical           :: pbl_mlturb_diss = .false.
   namelist /physics_cfgs/ pbl_mlturb_diss
   namelist /physics_cfgs_p/ pbl_mlturb_diss

   !# Run with legacy moistke clouds (no limits on cloud effects)
   logical           :: pbl_moistke_legacy_cloud = .false.
   namelist /physics_cfgs/ pbl_moistke_legacy_cloud
   namelist /physics_cfgs_p/ pbl_moistke_legacy_cloud

   !# Use the non-local PBL cloud formulation
   !# * 'NIL   ' : no non-local PBL cloud formulation
   !# * 'LOCK06' : Non-local cloud scheme of Lock and Mailhot (2006)
   character(len=16) :: pbl_nonloc   = 'NIL'
   namelist /physics_cfgs/ pbl_nonloc
   namelist /physics_cfgs_p/ pbl_nonloc
   character(len=*), parameter :: PBL_NONLOC_OPT(2) = (/ &
        'NIL   ',  &
        'LOCK06'  &
        /)

   !# Use the mixing length to average the Richardson number profile of (potentially)
   !# many layers to derive a "background" Ri estimate
   logical           :: pbl_ribkg    = .false.
   namelist /physics_cfgs/ pbl_ribkg
   namelist /physics_cfgs_p/ pbl_ribkg

   !# Richardson num. critical values for hysteresis
   real              :: pbl_ricrit(2)= 1.
   namelist /physics_cfgs/ pbl_ricrit
   namelist /physics_cfgs_p/ pbl_ricrit

   !# PBL representation of boundary layer clouds
   !# * 'NIL     ': No Shallow convection
   !# * 'CONRES  ': Bulk Richardson number-based turbulent enhancement
   !# * 'SHALOW  ': Deprecated (see 1998 RPN physics doc)
   !# * 'SHALODQC': Deprecated (see 1998 RPN physics doc)
   !# * 'GELEYN  ': Deprecated (see 1998 RPN physics doc)
   character(len=16) :: pbl_shal     = 'NIL'
   namelist /physics_cfgs/ pbl_shal
   namelist /physics_cfgs_p/ pbl_shal
   character(len=*), parameter :: PBL_SHAL_OPT(5) = (/ &
        'NIL     ', &
        'CONRES  ', &
        'SHALOW  ', &
        'SHALODQC', &
        'GELEYN  '  &
        /)

   !# Layer over which to adjust from SL to PBL stability functions [(bot,top) in m]
   real, dimension(2) :: pbl_slblend_layer = (/-1.,-1./)
   namelist /physics_cfgs/ pbl_slblend_layer
   namelist /physics_cfgs_p/ pbl_slblend_layer

   !# Adjustment to coefficient for TKE diffusion
   real              :: pbl_tkediff  = 1.
   namelist /physics_cfgs/ pbl_tkediff
   namelist /physics_cfgs_p/ pbl_tkediff

   !# Control of time scale for TKE diffusion
   logical           :: pbl_tkediff2dt  = .false.
   namelist /physics_cfgs/ pbl_tkediff2dt
   namelist /physics_cfgs_p/ pbl_tkediff2dt

   !# Depth (Pa) of the always-turbulent near-surface layer in the PBL
   real              :: pbl_turbsl_depth  = 3000.
   namelist /physics_cfgs/ pbl_turbsl_depth
   namelist /physics_cfgs_p/ pbl_turbsl_depth

   !# Use RPNphy solver for diffusion equations from YSU coefficients
   logical           :: pbl_ysu_rpnsolve = .false.
   namelist /physics_cfgs/ pbl_ysu_rpnsolve
   namelist /physics_cfgs_p/ pbl_ysu_rpnsolve

   !# Relaxation timescale (s) for mixing length smoothing
   real              :: pbl_zntau    = 7200.
   namelist /physics_cfgs/ pbl_zntau
   namelist /physics_cfgs_p/ pbl_zntau

   !# Use true (motionless) surface boundary conditions for TKE diffusion
   logical           :: pbl_zerobc   = .false.
   namelist /physics_cfgs/ pbl_zerobc
   namelist /physics_cfgs_p/ pbl_zerobc

   !# Scheme to determine precipitation type
   !# * 'NIL     ': no call to bourge
   !# * 'BOURGE  ': use Bourgouin algorithm (bourge1) to determine precip. types.
   !# * 'BOURGE3D':
   !# * 'SPS_W19 ': phase separation based on near-surface wet-bulb temperature (from Wang et al., 2019). Only for SPS
   !# * 'SPS_FRC ': fraction of each precipitation type is read directly in the atmospheric forcing (for SPS only)
   !# * 'SPS_H13 ': phase separation based on near-surface hydrometeor temperature (from Harder and Pomeroy, 2013). Only for SPS
   character(len=16) :: pcptype      = 'NIL'
   namelist /physics_cfgs/ pcptype
   namelist /physics_cfgs_p/ pcptype
   character(len=*), parameter :: PCPTYPE_OPT(6) = (/ &
        'NIL     ', &
        'BOURGE  ', &
        'BOURGE3D', &
        'SPS_W19 ', &
        'SPS_FRC ', &
        'SPS_H13 ' &
        /)

   !# Use diagnostic pcp types for surface for MP scheme when .true.
   !# otherwise use MP scheme pcp types
   logical           :: mpdiag_for_sfc = .false.
   namelist /physics_cfgs/ mpdiag_for_sfc
   namelist /physics_cfgs_p/ mpdiag_for_sfc

   !# Print stats for phy_input read var
   logical           :: phystat_input_l = .false.
   namelist /physics_cfgs/ phystat_input_l
   namelist /physics_cfgs_p/ phystat_input_l

   !# Use double presision for physic statistics output
   logical           :: phystat_dble_l = .false.
   namelist /physics_cfgs/ phystat_dble_l
   namelist /physics_cfgs_p/ phystat_dble_l

   !# Physic statistics output for 3d varables:
   !# * .false. : mean, var, min and max for the whole 3d fiels
   !# * .true.  : mean, var, min and max are done for each levels independently
   logical           :: phystat_2d_l = .false.
   namelist /physics_cfgs/ phystat_2d_l
   namelist /physics_cfgs_p/ phystat_2d_l

   !# Physic statistics output Frequency
   character(len=16) :: phystat_freq_S = '0h'
   namelist /physics_cfgs/ phystat_freq_S
   namelist /physics_cfgs_p/ phystat_freq_S

   !# Physic statistics output: bus variable list that should be included in physics
   !# "block" stats. Possible values:
   !# * Long varnames
   !# * Short varnames
   !# * 'ALLVARS=EDPV': all variables from E, D, P, V buses (any combination of the 4 letters);
   character(len=32) :: phystat_list_s(2048) = ' '
   namelist /physics_cfgs/ phystat_list_s
!!$   namelist /physics_cfgs_p/ phystat_list_s


   !# CFC11 bckgrnd atmospheric concentration (PPMV)
   real              :: qcfc11       = -1.
   namelist /physics_cfgs/ qcfc11
   namelist /physics_cfgs_p/ qcfc11
   real, parameter :: QCFC11_DEFAULT = 0.280

   !# CFC12 bckgrnd atmospheric concentration (PPMV)
   real              :: qcfc12       = -1
   namelist /physics_cfgs/ qcfc12
   namelist /physics_cfgs_p/ qcfc12
   real, parameter :: QCFC12_DEFAULT = 0.530


   !# CH4 bckgrnd atmospheric concentration (PPMV)
   real              :: qch4         = -1.
   namelist /physics_cfgs/ qch4
   namelist /physics_cfgs_p/ qch4
   real, parameter :: QCH4_DEFAULT   = 1.783

   !# CO2 bckgrnd atmospheric concentration (PPMV)
   real              :: qco2         = -1.
   namelist /physics_cfgs/ qco2
   namelist /physics_cfgs_p/ qco2
   real, parameter :: QCO2_DEFAULT   = 380.

   !# N2O bckgrnd atmospheric concentration (PPMV)
   real              :: qn2o         = -1.
   namelist /physics_cfgs/ qn2o
   namelist /physics_cfgs_p/ qn2o
   real, parameter :: QN2O_DEFAULT   = 0.3186

   !# Atmospheric path length for solar radiation
   !# * 'RODGERS67' : Formulation used by Li and Barker (2005)
   !# * 'LI06' : Estimate of  Li and Shibata (2006)
   character(len=16) :: rad_atmpath = 'RODGERS67'
   namelist /physics_cfgs/ rad_atmpath
   namelist /physics_cfgs_p/ rad_atmpath
   character(len=*), parameter :: RAD_ATMPATH_OPT(2) = (/ &
        'RODGERS67', &
        'LI06     '  &
        /)

   !# Phase partition of total water content for radiation when CONSUN is used
   !# * 'BOUOPS' : Boudala et al. (2004), QJRMS, 130, pp. 2919-2931 - bugged
   !# * 'BOUDALA' :Boudala et al. (2004), QJRMS, 130, pp. 2919-2931
   !# * 'ECMWF' : IFS docu CY25R1
   !# * 'Rockel' : Rockel et al. Beitr. Atmos. Phy. 1991
   character(len=16) :: rad_part_nomp  = 'BOUOPS'
   namelist /physics_cfgs/ rad_part_nomp
   namelist /physics_cfgs_p/ rad_part_nomp
   character(len=*), parameter :: RAD_PART_NOMP_OPT(4) = (/ &
        'BOUOPS  ', &
        'BOUDALA ',  &
        'ECMWF   ',  &
        'ROCKEL  '  &
        /)

   !# Optical properties of ice cloud from condensation scheme for radiation
   !# * 'CCCMA'      : Radius estimate from CCCMA
   !# * 'SIGMA'      : Altitude-dependent ice radius
   !# * 'ECMWF'      : Radius estimate from ecmwf
   character(len=16) :: rad_cond_rei = '15.'
   real              :: rei_const    = -1.
   namelist /physics_cfgs/ rad_cond_rei
   namelist /physics_cfgs_p/ rad_cond_rei
   character(len=*), parameter :: RAD_COND_REI_OPT(3) = (/ &
        'CCCMA', &
        'SIGMA',  &
        'ECMWF'  &
        /)

   !# Optical properties of liquid cloud from condensation scheme for radiation
   !# * 'BARKER'     : Radii based on aircraft data (Slingo), range 4-17 microns
   !# * 'NEWRAD'     : The "new" optical properties used in the "newrad" scheme
   !# * 'ROTSTAYN03' : Radii based on Rotstayn and Liu (2003)
   character(len=16) :: rad_cond_rew = 'BARKER'
   real              :: rew_const    = -1.
   namelist /physics_cfgs/ rad_cond_rew
   namelist /physics_cfgs_p/ rad_cond_rew
   character(len=*), parameter :: RAD_COND_REW_OPT(3) = (/ &
        'BARKER    ', &
        'NEWRAD    ', &
        'ROTSTAYN03'  &
        /)

   !# Conservation corrections for radiation scheme
   !# * 'NIL ' : No conservation correction applied
   !# * 'TEND' : Temperature and moisture tendencies corrected
   character(len=16) :: rad_conserve = 'NIL'
   namelist /physics_cfgs/ rad_conserve
   namelist /physics_cfgs_p/ rad_conserve
   character(len=*), parameter :: RAD_CONSERVE_OPT(2) = (/ &
        'NIL ', &
        'TEND'  &
        /)

   !# Use emissivity computed by the surface schemes
   logical           :: rad_esfc     = .false.
   namelist /physics_cfgs/ rad_esfc
   namelist /physics_cfgs_p/ rad_esfc

   !# Compute and apply tendencies from longwave radiation
   logical :: rad_lw = .true.
   namelist /physics_cfgs/ rad_lw
   namelist /physics_cfgs_p/ rad_lw
   
   !# For calculation of DIAGNOSTIC low, mid and high TRUE and EFFECTIVE cloud covers in cldoppro and cldoppro_mp
   !# TRUE:      rad_siglim(1)=limit between low and mid clouds in sigma; rad_siglim(2)=limit between mid and high clouds in sigma; 
   !# EFFECTIVE: rad_siglim(3)=limit between low and mid clouds in sigma; rad_siglim(4)=limit between mid and high clouds in sigma; 
   real, dimension(4) :: rad_siglim = (/0.7,0.4,0.7,0.4/)
   namelist /physics_cfgs/ rad_siglim
   namelist /physics_cfgs_p/ rad_siglim

   !# Fix use of effective solar zenith angle
   logical           :: rad_sun_angle_fix_l = .false.
   namelist /physics_cfgs/ rad_sun_angle_fix_l
   namelist /physics_cfgs_p/ rad_sun_angle_fix_l

   !# Compute and apply tendencies from shortwave radiation
   logical :: rad_sw = .true.
   namelist /physics_cfgs/ rad_sw
   namelist /physics_cfgs_p/ rad_sw
   
   !# For calculation of DIAGNOSTIC low, mid and high TRUE cloud covers in cldoppro and cldoppro_mp with height criteria for Calipso-GOCCP
   !# TRUE:      rad_zlim(1)=limit between low and mid clouds in height; rad_zlim(2)=limit between mid and high clouds in height; 
   real, dimension(2) :: rad_zlim = (/3200.,6500./)
   namelist /physics_cfgs/ rad_zlim
   namelist /physics_cfgs_p/ rad_zlim

   !# Use climatological values of GHG in radiation (CCCMARAD2 only)
   logical           :: radghg_L     = .false.
   namelist /physics_cfgs/ radghg_L
   namelist /physics_cfgs_p/ radghg_L

   !# Use LINOZ prognostic Ozone in radiation (CCCMARAD2 .and. LINOZ only)
   logical           :: rad_linoz_L  = .false.
   namelist /physics_cfgs/ rad_linoz_L
   namelist /physics_cfgs_p/ rad_linoz_L

   !# LINOZ prognostic stratospheric ozone
   !# * 'NIL     ' :
   !# * 'OZONE   ' :
   !# * 'GHG     ' :
   !# * 'OZONEGHG' :
   character(len=10) :: linoz_chm    = 'NIL'
   namelist /physics_cfgs/ linoz_chm
   namelist /physics_cfgs_p/ linoz_chm
   character(len=*), parameter :: LINOZ_CHM_OPT(4) = (/ &
        'NIL     ', &
        'OZONE   ', &
        'GHG     ', &
        'OZONEGHG'  &
        /)

   !# UV Method Optical properties of liquid cloud from condensation scheme for radiation
   !# * 'NIL'        : No UV calculation
   !# * 'IntegFit'   : Integration over the four UV bands with integrand weighting
   !# * 'LinearFit'  : Fitted sum of the four UV broadband irradiances
   !# * 'BandRatio'  : Scaling of input clear-sky UV index
   character(len=10) :: iuv_method   = 'IntegFit'
   namelist /physics_cfgs/ iuv_method
   namelist /physics_cfgs_p/ iuv_method
   character(len=*), parameter :: IUV_METHOD_OPT(4) = (/ &
        'NIL      ', &
        'INTEGFIT ', &
        'LINEARFIT', &
        'BANDRATIO'  &
        /)

   !# Radiation scheme
   !# * 'NIL      ': no radiation scheme
   !# * 'CCCMARAD ': most advanced radiation scheme
   !# * 'CCCMARAD2': most advanced radiation scheme v2
   character(len=16) :: radia        = 'NIL'
   namelist /physics_cfgs/ radia
   namelist /physics_cfgs_p/ radia
   character(len=*), parameter :: RADIA_OPT(3) = (/ &
        'NIL      ', &
        'CCCMARAD ', &
        'CCCMARAD2'  &
        /)

   !# Key for activation of the radiation along slopes
   logical           :: radslope     = .false.
   namelist /physics_cfgs/ radslope
   namelist /physics_cfgs_p/ radslope

   !# Launching level value of GW RMS wind (m/s) from non-orographic origin
   real              :: rmscon       = 1.0
   namelist /physics_cfgs/ rmscon
   namelist /physics_cfgs_p/ rmscon

   !# water/ice phase for saturation calc. if .true.;
   !# water phase only for saturation calc. if .false.
   logical           :: satuco       = .true.
   namelist /physics_cfgs/ satuco
   namelist /physics_cfgs_p/ satuco

   !# Surface fluxes digital filter order
   !# * -1 : No filter
   !# * 2 or 4: Apply filter 2nd or 4th order respectively
   integer :: sfcflx_filter_order = -1
   namelist /physics_cfgs/ sfcflx_filter_order
   namelist /physics_cfgs_p/ sfcflx_filter_order

   !# Surface fluxes figital filter, number of iterations
   integer :: sfcflx_filter_iter = 1
   namelist /physics_cfgs/ sfcflx_filter_iter
   namelist /physics_cfgs_p/ sfcflx_filter_iter
   
   !# Tuning factor for blocking height
   real              :: sgo_bhfac    = 1.5
   namelist /physics_cfgs/ sgo_bhfac
   namelist /physics_cfgs_p/ sgo_bhfac

   !# Sets the minimum value of the drag coefficient in the orographic
   !# blocking scheme.
   real              :: sgo_cdmin    = 1.0
   namelist /physics_cfgs/ sgo_cdmin
   namelist /physics_cfgs_p/ sgo_cdmin

   !# Turns on/off the non-linear amplification factor (depending on wind
   !# direction) of the drag coefficient in the orographic blocking scheme
   logical           :: sgo_nldirfac = .true.
   namelist /physics_cfgs/ sgo_nldirfac
   namelist /physics_cfgs_p/ sgo_nldirfac

   !# Critical phase for blocking height
   real              :: sgo_phic    = 0.2
   namelist /physics_cfgs/ sgo_phic
   namelist /physics_cfgs_p/ sgo_phic

   !# Turns on/off the amplification factor (due to stability) of the drag
   !# coefficient in the orographic blocking scheme
   logical           :: sgo_stabfac  = .true.
   namelist /physics_cfgs/ sgo_stabfac
   namelist /physics_cfgs_p/ sgo_stabfac

   !# Standard deviation length scale (gridpoints) of Gaussian smoother
   !# applied to wind GWD tendencies
   real           :: sgo_tdfilter    = 1.
   namelist /physics_cfgs/ sgo_tdfilter
   namelist /physics_cfgs_p/ sgo_tdfilter

   !# Description of threshold for mean wind speed for blocking
   real, dimension(2) :: sgo_windfac = (/2.,0.01/)
   namelist /physics_cfgs/ sgo_windfac
   namelist /physics_cfgs_p/ sgo_windfac

   !# Condensation scheme name
   !# * 'NIL       ' : No explicit condensation scheme used
   !# * 'CONSUN    ' : Sunqvist type condensation scheme
   !# * 'MP_MY2    ' : Milbrandtl and Yau microphysics scheme
   !# * 'MP_P3     ' : P3 microphysics scheme
   !# * 'KESSLER   ' : Kessler warm rain scheme
   character(len=16) :: stcond       = 'NIL'
   namelist /physics_cfgs/ stcond
   namelist /physics_cfgs_p/ stcond
   character(len=*), parameter :: STCOND_OPT(5) = (/ &
        'NIL       ', &
        'CONSUN    ', &
        'MP_MY2    ', &
        'MP_P3     ', &
        'KESSLER   ' &
        /)

   !# Special treatment of stratosphere;
   !# if .true. ignore convection/condensation tendencies where pressure is lower
   !# than topc as specified in nocld.cdk
   logical           :: stratos      = .false.
   namelist /physics_cfgs/ stratos
   namelist /physics_cfgs_p/ stratos

   !# Factor used in the gwd formulation = 1/(LENGTH SCALE)
   real              :: taufac       = 8.E-6
   namelist /physics_cfgs/ taufac
   namelist /physics_cfgs_p/ taufac

   !# Run the physics in test harness mode
   logical           :: test_phy     = .false.
   namelist /physics_cfgs/ test_phy
   namelist /physics_cfgs_p/ test_phy

   !# Print runtime timings
   logical           :: timings_L    = .false.
   namelist /physics_cfgs/ timings_L
   namelist /physics_cfgs_p/ timings_L

   !# Select a turbulent orographic form drag scheme
   !# * 'NIL'        : No turbulent orographic form drag scheme
   !# * 'BELJAARS04' : Form drag scheme described by Beljaars et al. (2006; QJRMS)
   !# WARNING: This option is broken thus disabled- will be fixed in dev branch
   character(len=16) :: tofd         = 'NIL'
   namelist /physics_cfgs/ tofd
   namelist /physics_cfgs_p/ tofd
   character(len=*), parameter :: TOFD_OPT(2) = (/ &
        'NIL       ', &
        'BELJAARS04'  &
        /)

   !# BELJAARS04 turbulent orographic form drag scheme alpha parameter
   real :: tofd_alpha = 12.
   namelist /physics_cfgs/ tofd_alpha
   namelist /physics_cfgs_p/ tofd_alpha

   !# use Latent Heat Nudging for the assmilation of radar-inferred precipitation rates
   !# * 'NIL'   : No Latent Heat Nudging
   !# * 'IRPCP' : Latent heat nudging from radar-inferred precipitation rates
   character(len=16) :: lhn         = 'NIL'
   namelist /physics_cfgs/ lhn
   namelist /physics_cfgs_p/ lhn
   character(len=*), parameter :: LHN_OPT(2) = (/ &
        'NIL  ', &
        'IRPCP'  &
        /)
 
   !# Standard deviation length scale (gridpoints) of Gaussian smoother
   !# applied to RDPR, PR and TA before the application of Latent Heat Nudging
   !# No smoothing applied when -ve
   real           :: lhn_filter    = -1.
   namelist /physics_cfgs/ lhn_filter
   namelist /physics_cfgs_p/ lhn_filter

   !# Modulation factor for the magnitude of Latent Heat Nudging being applied
   !# modulated_tendencies = lhn_weight*(LHN tendencies)
   real           :: lhn_weight    = 0.
   namelist /physics_cfgs/ lhn_weight
   namelist /physics_cfgs_p/ lhn_weight

   !# first time step at which LHN is applied
   character(len=16) :: lhn_start_S = '10p'
   namelist /physics_cfgs/ lhn_start_S
   namelist /physics_cfgs_p/ lhn_start_S

   !# last time step at which LHN is applied
   character(len=16) :: lhn_stop_S = '360p'
   namelist /physics_cfgs/ lhn_stop_S
   namelist /physics_cfgs_p/ lhn_stop_S

   !# To avoid shocking the model, LHN is turned on gradually 
   !# this parameter controls how many time steps it takes (after lhn_timestep_start)
   !# for LHN to be fully active
   character(len=16) :: lhn_ramp_S = '10p'
   namelist /physics_cfgs/ lhn_ramp_S
   namelist /physics_cfgs_p/ lhn_ramp_S

contains

   function phy_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      indiag_list_s(1) = 'DEFAULT LIST'
      return
   end function phy_options_init

end module phy_options
