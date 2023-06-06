!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : chm_nml_mod.ftn90
! Creation       : H. Landry, Mai 2008
! Description    : Modules defining namelist keys
!
! Addition_1
!      By: A. Kallaur
!      Date: September 19 2011
!      Synopsis: Add logical variables  "chm_master_switch" and "chm_emis_master_switch" as
!                as part of high level control mechanism.
! Extra info     :
!
!============================================================================

module chm_nml_mod
   use chm_utils_mod, only: NMLKEY_LEN, NOMV_LEN

   save

   character(len=NMLKEY_LEN) chm_model_s           ! Name of chemistry model                (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_pkg_gas_s         ! Name of gas species package            (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_pkg_pm_s          ! Name of PM species package             (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_strato_s          ! Name of Stratospheric O3 scheme        (Default = 'NIL')
   character(len=NMLKEY_LEN) chm_biog_s            ! Name of emissions package              (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_gas_drydep_s      ! Name of dry deposition scheme of gas   (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_soa_s             ! Name of SOA scheme                     (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_vert_diff_s       ! Name of vertical diffusion algo        (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_mj_treatment_s    ! Name of major point treatment algo     (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_hetchem_s         ! Name of heterogeneous chemistry scheme (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_seaflux_s         ! Name of sea salt flux scheme           (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_winddust_s        ! Name of wind dust package              (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_pm_drydep_s       ! Name of dry deposition scheme of PM    (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_aqueous_s         ! Name of aqueous phase chemistry scheme (Default = 'undefined')
   character(len=NMLKEY_LEN) chm_met_modulation_s  ! Meteorological modulation of fugitive dust emissions (Default = 'OFF')
   character(len=NMLKEY_LEN) chm_wdep_scav_coef_s  ! Switch to change below-cloud scavenging parameterization for particles (Default = 'SLINN')
   character(len=NMLKEY_LEN) chm_mass_s            ! Name of mass conservation scheme       (Default = 'bc'  )
   character(len=NMLKEY_LEN) chm_mono_s            ! Name of monotonicity scheme            (Default = 'clip')
   character(len=NMLKEY_LEN) chm_ammonia_bidi_s    ! Name of ammonia bidirectional flux scheme (Default = 'OFF')

   logical(kind=4) chm_master                        ! If false, then GEM should run as if in chm_stub mode  (Def=T)
   logical(kind=4) chm_emis_master                   ! if false, reading all types (mobile,biogenic, major) of emissions are off (Def=T)
   logical(kind=4) chm_get_ae_emis_l                 ! Do we import area emissions?                (Default = .false.)
   logical(kind=4) chm_get_be_emis_l                 ! Do we import biogenic emissions?            (Default = .false.)
   logical(kind=4) chm_get_mj_emis_l                 ! Do we import major point sources emissions? (Default = .false.)
   logical(kind=4) chm_get_wf_emis_l                 ! Do we import CFFEPS wildfire emissions      (Default = .false.)
   logical(kind=4) chm_get_on_emis_l                 ! Do we import on-road area emissions      (Default = .false.)
   logical(kind=4) chm_get_fd_emis_l                 ! Do we import fugitive dust area emissions?  (Default = .false.)
   logical(kind=4) chm_htap_emis_l                   ! Are we running with the HTAP set of emissions?  (Default = .false.)
   integer(kind=4) chm_step_factor                   ! Chemistry to physics timestep factor        (Default = '-1')
   logical(kind=4) chm_hetchem_metstlb_l             ! Heterogenous metastable option              (Default = .false.)
   logical(kind=4) chm_blcld_l                       ! Do we apply below cloud process?            (Default = .false.)
   logical(kind=4) chm_split_snow_rain_l             ! whether to keep liquid and solid precipitation separate in below cloud scavenging (Default = .false.)
   logical(kind=4) chm_nosurfzoneflux_l              ! Do we exclude surf-zone sea-salt flux       (Default = .false.)
   integer(kind=4) chm_o3icedep                      ! Use different deposition rate over ice./glaciers for O3
   integer(kind=4) chm_pm_coag_step_intvl            ! Interval of time steps to calculate coagulation coefficient for PM (Default='-1')
   integer(kind=4) chm_debug_2d_i, chm_debug_3d_i    ! Number of debug variables respectively for 2D and 3D fields (default = -1, max = 18)
   integer(kind=4) chm_intrsec_ndiv                  ! Number of subdivisions per bin in mach_cam_intrsec*
   integer(kind=4) chm_intrsec_ver                   ! Version number of the code in mach_cam_intrsec*
   logical(kind=4) chm_aqhi_l                        ! Output of AQHI index: AQ25 and AQ10         (Default = .false.)
   logical(kind=4) chm_ss_ao_l                       ! Enable including sea salt to the aerosol-concentration output (e.g. AF and AC) (Default = .false.)
   logical(kind=4) chm_diag_wetdep_l                 ! Diagnostic for wet deposition               (Default = .false.)
   logical(kind=4) chm_diag_drydep_l                 ! Diagnostic for dry deposition               (Default = .false.)
   logical(kind=4) chm_diag_accum_l                  ! Diagnostic for accumulators                 (Default = .false.)
   logical(kind=4) chm_diag_aerosols_l               ! Diagnostic for CAM species output in ug/m3  (Default = .false.)
   logical(kind=4) chm_diag_colum_l                  ! Diagnostic for column values                (Default = .false.)
   logical(kind=4) chm_diag_aero_opt_l               ! Diagnostic for aerosol optical properties   (Default = .false.)
   logical(kind=4) chm_timings_l                     ! True to enable timings in chemistry         (Default = .false.)
   logical(kind=4) chm_debug_trace_l                 ! True to enable debug tracing in chemistry   (Default = .false.)
   logical(kind=4) chm_messy_jval_l                  ! If true invoke MESSY for photolysis rates   (Default = .false.)
   logical(kind=4) chm_active_ch4_l                  ! True to enable active CH4 tracer            (Default = .false.)
   logical(kind=4) chm_vit_l                         ! True to enable vehicle-induced turbulence   (Default = .false.)
   logical(kind=4) chm_direct_l                      ! True to enable direct radiative feedback    (Default = .false.)
   logical(kind=4) chm_indirect_l                    ! True to enable indirect (cloud) feedback    (Default = .false.)
   logical(kind=4) chm_canopy_shading_l              ! True to enable forest canopy shading        (Default = .false.)
   logical(kind=4) chm_dep_lai2d_l                   ! True to use the satellite-derived LAI field
                                                     ! in gas dry deposition                       (Default = .false.)
   logical(kind=4) chm_sat_seasons_l                 ! True to enable satellite LAI - based adjustment of seasons in
                                                     ! biogenic emissions modules (Default = .false.)
   logical(kind=4) chm_cffeps_online_l               ! True to enable online CFFEPS calculations   (Default = .false.)
   integer(kind=4) chm_moyhr                         ! Accumulation interval in hours              (Default = 24     )
   integer(kind=4) chm_acchr                         ! Running mean interval in hours              (Default = 24     )
   logical(kind=4) chm_mar_halo_l                    ! Use parameterization for marine halogen chemistry (Default = .false.)
   integer(kind=4) wf_case                           ! type of plumerise for wildfire emissions    (Default = 0      )

   integer(kind=4)  nk_start       ! chemistry activated only from level nk_start to level chm_nk
   integer(kind=4)  nk_start_pm    ! aerosol chemistry activated only from level nk_start_pm to level chm_nk
   real(kind=4)     dbg_lat        ! latitude of debug location within model domain
   real(kind=4)     dbg_lon        ! longitude of debug location within model domain
   integer(kind=4)  dbg_lev        ! model level number for diagnostic outputs of cam species
   integer(kind=4)  dbg_step       ! time step for diagnostic outputs of cam species
   integer(kind=4)  dbg_itr        ! tracer index for diagnostic outputs of cam species

   real(kind=4) aerosize(13)        ! aerosol bin boundaries; depends upon number of bins
   real(kind=4) chm_bkgd_ch4        ! Background CH4 concentration (Default = 1.8ppm)
   real(kind=4) chm_bkgd_co2        ! Background CO2 concentration (Default = 400.0ppm)
   real(kind=4) chm_kt_minmax(3)    ! 1:min 2:max 3:urban-min vertical diffusion coeff. override GEM input Kt in mach_input_check
   real(kind=4) chm_urban_abl_min   ! Nighttime ABL increment due to urban heat island generated using TEB, from Ren et al
   real(kind=4) aero_opt_wavel(15)  ! Mid-band wavelengths (in meters) for AOD calculation
   real(kind=4) chm_ammonia_gep(15) ! Ground emission potential values used for ammonia bidirectional flux for option chm_ammonia_bidi_s = 'GEP' (15 land-use categories)

   character(len=NOMV_LEN) gas_ppb_out_list_s(40)    ! Comma-separated list of (ppb output name) of
                                                     ! gas species to output in unit of ppb (Default = 'O3', 'N2', 'NO', 'S2')

   real(kind=4)  em_nesdt          ! Emissions preperation timestep (s)

   logical(kind=4) chm_pblh_min_l  ! Impose the minimum of 100 m for PBL height      (Default = .false.)
   logical(kind=4) chm_ae_spread_l ! Impose the spread of area emissions to 2 layers (Default = .false.)

   logical(kind=4) chm_do_mjpts_l  ! Do we have any major point sources for local grid (tile)?
   logical(kind=4) chm_do_cffeps_l ! Do we have any cffeps fire hotspot for local grid (tile)? (Def. = .false.)

   namelist /chemistry_cfgs/ chm_master,            chm_model_s,            chm_pkg_pm_s,       chm_pkg_gas_s,     &
                             chm_biog_s,            chm_strato_s,           chm_mass_s,         chm_mono_s,        &
                             chm_vert_diff_s,       chm_messy_jval_l,       chm_gas_drydep_s,   chm_pm_drydep_s,   &
                             chm_dep_lai2d_l,       chm_sat_seasons_l,      chm_o3icedep,       chm_active_ch4_l,  &
                             chm_emis_master,       chm_get_ae_emis_l,      chm_get_be_emis_l,  chm_get_mj_emis_l, &
                             chm_get_wf_emis_l,     chm_get_on_emis_l,      chm_htap_emis_l,    chm_mj_treatment_s,&
                             wf_case,               chm_canopy_shading_l,   chm_cffeps_online_l,chm_vit_l,         &
                             chm_met_modulation_s,  chm_seaflux_s,          chm_winddust_s,     chm_soa_s,         &
                             chm_hetchem_s,         chm_hetchem_metstlb_l,  chm_aqueous_s,      chm_blcld_l,       &
                             chm_split_snow_rain_l, chm_wdep_scav_coef_s,   chm_aqhi_l,         chm_ss_ao_l,       &
                             chm_pm_coag_step_intvl,chm_intrsec_ndiv,       chm_intrsec_ver,                       &
                             chm_direct_l,          chm_indirect_l,         chm_get_fd_emis_l,  chm_mar_halo_l,    &
                             chm_diag_wetdep_l,     chm_diag_drydep_l,      chm_diag_accum_l,   chm_diag_colum_l,  &
                             chm_diag_aerosols_l,   chm_diag_aero_opt_l,    chm_moyhr,          chm_acchr,         &
                             chm_step_factor,       chm_kt_minmax,          chm_timings_l,      chm_debug_trace_l, &
                             chm_debug_2d_i,        chm_debug_3d_i,         chm_ammonia_bidi_s, chm_ammonia_gep,   &
                             em_nesdt,              aerosize,               nk_start,           nk_start_pm,       &
                             gas_ppb_out_list_s,    chm_bkgd_ch4,           chm_bkgd_co2,       aero_opt_wavel,    &
                             dbg_lat,               dbg_lon,                dbg_lev,            dbg_step, dbg_itr, &
                             chm_pblh_min_l,        chm_ae_spread_l,        chm_nosurfzoneflux_l

end module chm_nml_mod
