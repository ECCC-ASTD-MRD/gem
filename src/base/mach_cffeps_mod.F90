!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_cffeps_mod.ftn90
! Creation         : Kerry Anderson, A. Akingunola, J. Chen, and P. Makar
!                    - Fall 2018
! Description      :
!
!============================================================================
module mach_cffeps_mod
   implicit none

   save
   integer(kind=4), parameter :: MET_LEVELS = 40
   integer(kind=4), parameter :: LABEL_SIZE = 32
   integer(kind=4), parameter :: NB_CHAR_PATH = 240
   integer(kind=4), parameter :: FUEL_SIZE = 4

   integer(kind=4), parameter :: max_fuels = 20

   character(len=fuel_size), dimension(max_fuels), parameter :: &
                         cffeps_fueltypes = &
                    (/'C1  ', 'C2  ', 'C3  ', 'C4  ', 'C5  ', 'C6  ', &
                      'C7  ', 'D1  ', 'D2  ', 'M1  ', 'M2  ', 'M3  ', &
                      'M4  ', 'S1  ', 'S2  ', 'S3  ', 'O1a ', 'O1b ', &
                      'WA  ', 'NF  '/)
!   integer(kind=4), parameter :: ifuel_c2 = 2, ifuel_d1 = 8

   integer(kind=4), parameter :: cffeps_nspec = 9
   character(len=5), dimension(cffeps_nspec), parameter :: cffeps_species = &
                       (/'PM10 ', 'PM2.5', 'CO   ', 'CO2  ', 'CH4  ', &
                         'NOX  ', 'NH3  ', 'SO2  ', 'NMHC '/)

   real(kind=4), parameter :: msl_pressure = 1013.25   ! /* MSL pressure [mb] */

   integer(kind=4), dimension(12), parameter :: jdays = &
                    (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

   character(len=LABEL_SIZE) :: cffeps_method      ! method of employment: AB= Alberta Plume Study ICAO=Standard Atmosphere test anything else=cmc calculations
   character(len=LABEL_SIZE) :: cffeps_fire_shape  ! shape of entrainment cloud: line/wedge, crescent, ellipse
   character(len=LABEL_SIZE) :: cffeps_fire_type   ! type of upper air profile method used (average lapse rate, dry, wet)
   character(len=LABEL_SIZE) :: cffeps_diurnal     ! choice of diurnal approaches: CVW, BDL, EMC, OFF -- 2020-06-04
   real(kind=4)              :: cffeps_timestep    ! timestep for modelling black body radiation loss
   real(kind=4)              :: cffeps_alpha       ! entrainment half-angle (in radians)
   real(kind=4)              :: cffeps_reset       ! reset time of smoke plume calculations [decimal hours LST] (off if < 0)
   real(kind=4)              :: cffeps_thstart     ! time to start (exclusive) top-hat fire growth [decimal hours LST]
   real(kind=4)              :: cffeps_thend       ! time to end (inclusive) top-hat fire growth [decimal hours LST]

   integer(kind=4)           :: cffeps_LDT         ! local daylight time: 1=Yes, 0=No
   integer(kind=4)           :: cffeps_radiation   ! radiation method: < 0 Byram's 1200/8600 =0 A_wall/A_top > 0  value as % of Qfire
   integer(kind=4)           :: cffeps_fc_method   ! Fuel consumption method (1 = default FBP)

   logical(kind=4)           :: cffeps_sinks       ! include the sinks terms: 1=Yes, 0=No

   logical(kind=4)           :: cffeps_buieff      ! BUI effect for FBP calculations

   logical(kind=4)           :: cffeps_fbp_accel   ! .true. for point; .false. for line

   integer(kind=4)           :: aqm_aerobins       ! Number of aerosol size bins (offline use only)
   character(len=NB_CHAR_PATH) :: hotspotfile, met_input_dir

   real(kind=4),  dimension(24) :: weight

   integer(kind=4)            :: lc_hotspots        ! Number of hotspots in local tile
   integer(kind=4), parameter :: fire_parameters = 31
   integer(kind=4), parameter :: i_gilc  = 1,  i_gjlc  = 2,  i_fuel = 3,   &
                                 i_dj    = 4,  i_dtime = 5,  i_bui  = 6,   &
                                 i_ffmc  = 7,  i_dmc   = 8,  i_dc   = 9,   &
                                 i_fmc   = 10, i_area  = 11, i_cfl  = 12,  &
                                 i_sfl   = 13, i_ffl   = 14, i_wfl  = 15,  &
                                 i_cured = 16, i_pc    = 17, i_pdf  = 18,  &
                                 i_cbh   = 19, i_asp   = 20, i_gfl  = 21,  &
                                 i_gs    = 22, i_sd    = 23, i_sh   = 24,  &
                                 i_lfl   = 25, i_dfl   = 26, i_q0   = 27,  &
                                 i_emis  = 28, i_area0 = 29, i_rsmk = 30,  &
                                 i_zplm  = 31
   real(kind=4), pointer, dimension(:,:) :: lfire_info
   real(kind=4), pointer, dimension(:,:) :: lfire_emissions
   integer(kind=4), dimension(:), allocatable :: fire_me_species_index

! To be (re-)defined at initialization
   integer(kind=4) :: cffeps_steps
   integer(kind=4) :: MAX_TIMESTEPS

   type :: feps_type
! from the GEM forecast data
      integer(kind=4) :: fueltype   ! Fuel type (index of cffeps_fueltypes)

      integer(kind=4) :: dj         ! detection Julian Day
      real(kind=4)    :: dtime      ! detection time of fire (UTC in fractional hours)
                                    ! time all values are assumed to be collected

      real(kind=4)    :: lat        ! Latitude [decimal degrees]
      real(kind=4)    :: lon        ! Longitude [decimal degrees]

      real(kind=4)    :: bui        ! BUI
      real(kind=4)    :: ffmc       ! FFMC
      real(kind=4)    :: dmc        ! Duff Moisture Code
      real(kind=4)    :: dc         ! Drought Code

      real(kind=4)    :: fmc        ! Foliar Moisture Content

      real(kind=4)    :: cfl        ! Crown Fuel Load [kg/m^2]
      real(kind=4)    :: sfl        ! Surface Fuel Load [kg/m^2]
      real(kind=4)    :: ffl        ! Fine Fuel Load [kg/m^2]
      real(kind=4)    :: wfl        ! Woody Fuel Load [kg/m^2]

      real(kind=4)    :: curing     ! Percent Cured for O1a/O1b
      real(kind=4)    :: pc         ! Percent Confier for M1/M2
      real(kind=4)    :: pdf        ! Percent Dead Fir for M3/M4
      real(kind=4)    :: cbh        ! Crown to Base Height [m]
      real(kind=4)    :: aspect     ! Aspect [degrees]
      real(kind=4)    :: gfl        ! Grass Fuel Load [kg/m^2]
      real(kind=4)    :: gs         ! Slope [percent]
      real(kind=4)    :: sd         ! C6 Stand Density [stems/ha]
      real(kind=4)    :: sh         ! C6 Stand Height [m]
      real(kind=4)    :: lfl        ! Litter Fuel Load [kg/m^2]
      real(kind=4)    :: dfl        ! Duff (Fermentation) Fuel Load [kg/m^2]

      real(kind=4)    :: area       ! size of fire at dtime
      real(kind=4)    :: estarea    ! estimated fire size at time of detection

      real(kind=4)    :: Qo         ! amount of energy previously injected into atmosphere

! Currently fix the residence time to 24
      real(kind=4), dimension(24) :: qs

      real(kind=4)    :: totalemissions  ! total emissions from the fire (tonnes) (used in emissions.c)

   end type feps_type

   type :: fbp_type
! outputs
      real(kind=4) :: ros             ! Rate of Spread [m/min]
      real(kind=4) :: fros            ! Flank rate of Spread [m/min]
      real(kind=4) :: bros            ! Back Rate of Spread [m/min]
      real(kind=4) :: cfb             ! Crown Fraction Burned
      real(kind=4) :: hfi             ! Head Fire Intensity [kW/m]
      real(kind=4) :: tfc             ! Total Fuel Consumption [kg/m^2]
      real(kind=4) :: sfc             ! Surface Fuel Consumption [kg/m^2]
   end type fbp_type

   type :: emissions_type
! Currently fix the residence time to 24
      real(kind=4), dimension(24) :: f  ! cumulative emissions for flaming combustion
      real(kind=4), dimension(24) :: s  ! cumulative emissions for smoldering combustion
      real(kind=4), dimension(24) :: r  ! cumulative emissions for residual combustion
   end type emissions_type

   type :: met_type
      real(kind=4) :: ta             ! near-surface air temperature [K]
      real(kind=4) :: hus            ! near-surface specific humidity [kg/kg]
      real(kind=4) :: ps             ! Surface pressure [mb]
      real(kind=4) :: ws             ! near-surface wind speed [km/hr]
      real(kind=4) :: wd             ! near-surface wind direction [radians]
      real(kind=4) :: pr             ! Total precipitation rate (mm/hr)
      real(kind=4) :: P(MET_LEVELS)  ! Atmospheric pressure profile [Pa]
      real(kind=4) :: T(MET_LEVELS)  ! Atmospheric temperature profile [K]
      real(kind=4) :: Z(MET_LEVELS)  ! Atmospheric height (thermo) profile [m]
   end type met_type

! Consumption allocation table
! Table by: Fuel, Flaming, Smoldering, Residual, Remainder, Notes
   type :: allocations_type
      real(kind=4), dimension(4) :: sfc_l = (/0.9,  0.1,  0.0, 0.0/) ! Litter layer
      real(kind=4), dimension(4) :: sfc_f = (/0.1,  0.7,  0.2, 0.0/) ! Fermentation layer
      real(kind=4), dimension(4) :: sfc_h = (/0.0,  0.2,  0.8, 0.0/) ! Humus layer
      real(kind=4), dimension(4) :: ffc   = (/1.0,  0.0,  0.0, 0.0/) ! Slash - fine fuels
      real(kind=4), dimension(4) :: wfc   = (/0.4,  0.3,  0.3, 0.0/) ! Slash - woody fuels
      real(kind=4), dimension(4) :: gfc   = (/0.95, 0.05, 0.0, 0.0/) ! Grass fuel load
      real(kind=4), dimension(4) :: cfc_f = (/0.94, 0.06, 0.0, 0.0/) ! Canopy
      real(kind=4), dimension(4) :: cfc_b = (/0.6,  0.0,  0.0, 0.4/) ! Branchwood
      real(kind=4), dimension(4) :: cfc_s = (/0.25, 0.25, 0.0, 0.5/) ! Snags
   end type allocations_type

   namelist /cffeps_inputs_cfgs/ cffeps_method,  cffeps_fire_shape, cffeps_fire_type, &
                                 cffeps_diurnal, cffeps_timestep,   cffeps_alpha,     &
                                 cffeps_reset,   cffeps_thstart,    cffeps_thend,     &
                                 cffeps_LDT,     cffeps_radiation,  cffeps_fc_method, &
                                 cffeps_sinks,   cffeps_buieff,     cffeps_fbp_accel, &
                                 aqm_aerobins, hotspotfile, met_input_dir

   contains

   integer function mach_cffeps_nml(F_unit)
      use chm_utils_mod,  only: chm_lun_out
      implicit none
      integer(kind=4), intent(in) :: F_unit
!
!  Local variable declarations
!
      integer(kind=4)   :: read_status
      character(len=80) :: line

      !  Return +1 if mach_cffeps_nml works properly, else return -1
      mach_cffeps_nml = -1
!
!/* These are the default CFFEPS values */

      cffeps_method = "cmc"              !/* FIREWORK calculations */

      cffeps_fire_shape = "weighted"
      cffeps_fire_type = "wet"   !/* type of upper air profile method used (average lapse rate, dry, wet) */
      cffeps_alpha = 12.0        ! /* entrainment half-angle (0. = no entrainment) */
      cffeps_reset = -1.0        ! /* reset time of smoke plume calculations [decimal hours LST] (off if < 0) */
      cffeps_radiation = 0       ! /* method of calculating radiation term: < 0 Byram's 1200/8600; =0 A_wall/A_top; > 0  value as % of Qfire */
      cffeps_sinks = .true.      ! /* include heat sinks in calculation >0 = yes */
      cffeps_thstart = 9.0       ! /* time to start (exclusive) top-hat fire growth [decimal hours LST] */
      cffeps_thend = 21.0        ! /* time to end (inclusive) top-hat fire growth [decimal hours LST] */
      cffeps_ldt = 0             ! /* Local Daylight Time 1= one hour time offset */
! NOTE: The FBP Accel should be hotspot dependent, but it is fixed for now.
      cffeps_fbp_accel = .false. ! // false, no acceleration
      cffeps_buieff = .false.    ! /* Turn off BUI effect by default
      cffeps_fc_method = 1       ! Fuel consumption method (1 = default FBP)
      cffeps_diurnal = "OFF"     ! Equilibrium Moisture Content ('EMC'), or
                                 ! BDL tabular hourly FFMC ('BDL')

      hotspotfile = 'cffeps_hotspots.csv'  ! Input hotspotfile
      met_input_dir = 'meteo'    ! offline use only
      aqm_aerobins = 12          ! offline use only

!  Attempt a read
      read(F_unit, nml = cffeps_inputs_cfgs, iostat = read_status)
      if (read_status > 0) then
         if (chm_lun_out > 0) then
            write(chm_lun_out, *) ' Error in reading namelist cffeps_cfgs (s/r mach_cffeps_mod)'
            write(chm_lun_out, *) ' Read status = ', read_status
            backspace(F_unit)
            read(F_unit, '(A)') line
            write(chm_lun_out, *) 'Invalid line in cffeps_cfgs namelist is: ', trim(line)
         end if
         if (chm_lun_out > 0) &
            write(chm_lun_out, *) 'MACH_CFFEPS_NML - FAULTY CFFEPS NAMELIST -> PLEASE REVISE '
         mach_cffeps_nml = -1
         return
      else if (read_status < 0) then
         if (chm_lun_out > 0) &
            write(chm_lun_out, *) ' No cffeps_cfgs found in gem_settings --> ABORT ALL!!!!'
         mach_cffeps_nml = -1
         return
      end if
!
      mach_cffeps_nml = 1
      return
   end function mach_cffeps_nml


end module mach_cffeps_mod
