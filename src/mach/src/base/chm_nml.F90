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
! Creation       : Alexander Kallaur, S. Gaudreault, Sept. 2010
! Description    : Read chemistry namelists and initialize chemistry configuration
!
! Extra info     :
!
! Arguments:
!           IN     :  F_namelist-> name of namelist file
!                     lun_out   -> Unit number of diagnostic message outputs
!
!
!           IN/OUT : Integer value of chm_nml function (>0 = succesful)
!
!=============================================================================
!
!!if_on
integer function chm_nml (F_namelist, lun_out)
!!if_off
   use chm_utils_mod,       only: chm_lun_out, global_debug, undefined, chm_stop
   use chm_nml_mod
   use mach_cam_utils_mod, only: isize, icob, nbnd
   use mach_hetv_mod,      only: itermax
   use mach_incld_mod,     only: incld_dtmin
   use mach_mie_data_mod,  only: nwl_aod
   use mach_cffeps_mod,    only: mach_cffeps_nml
   use phy_options,        only: stcond
   implicit none
!!if_on
   character(len=*), intent(in) :: F_namelist
   integer(kind=4),  intent(in) :: lun_out
!!if_off
!
!  Local variable declarations
!
   integer(kind=4)   :: read_status, istatus
   logical(kind=4)   :: local_dbg
   integer(kind=4)   :: unf, nrec, err_open, err
   character(len=80) :: line
   integer  fnom, fclos
   external fnom, fclos
!
   chm_lun_out = lun_out
!    Assign default values for the chem_cfgs namelist variables.
!    In the case of a garbled read from nml file, easier to detect
!    where read error happened.
!
   chm_model_s          = undefined
   chm_pkg_gas_s        = undefined
   chm_pkg_pm_s         = undefined
   chm_master           = .true.
   chm_emis_master      = .true.
   chm_get_ae_emis_l    = .false.
   chm_get_be_emis_l    = .false.
   chm_get_mj_emis_l    = .false.
   chm_get_wf_emis_l    = .false.
   chm_get_on_emis_l    = .false.
   chm_get_fd_emis_l    = .false.
   chm_htap_emis_l      = .false.
   chm_biog_s           = undefined
   chm_gas_drydep_s     = undefined
   chm_soa_s            = undefined
   chm_strato_s         = 'NIL'
   chm_vert_diff_s      = undefined
   chm_mass_s           = 'BC'
   chm_mono_s           = 'CLIP'
   chm_mj_treatment_s   = undefined
   chm_seaflux_s        = undefined
   chm_nosurfzoneflux_l = .false.
   chm_winddust_s       = undefined
   chm_hetchem_s        = undefined
   chm_hetchem_metstlb_l= .false.
   chm_o3icedep         = 0
   chm_step_factor      = -1
   em_nesdt             = 0.0
   chm_pm_drydep_s      = undefined
   chm_aqueous_s        = undefined
   chm_met_modulation_s = 'OFF'
   chm_blcld_l          = .false.
   chm_split_snow_rain_l= .false.
   chm_wdep_scav_coef_s = 'SLINN'
   chm_pm_coag_step_intvl = -1
   chm_debug_2d_i       = -1
   chm_debug_3d_i       = -1
   chm_intrsec_ndiv     = -1
   chm_intrsec_ver      = -1
   chm_aqhi_l           = .false.
   chm_ss_ao_l          = .false.
   chm_diag_wetdep_l    = .false.
   chm_diag_drydep_l    = .false.
   chm_diag_accum_l     = .false.
   chm_diag_aerosols_l  = .false.
   chm_diag_colum_l     = .false.
   chm_diag_aero_opt_l  = .false.
   chm_mar_halo_l       = .false.
   chm_active_ch4_l     = .false.
   chm_moyhr            = 24
   chm_acchr            = 24
   chm_kt_minmax        = (/0.1, 1500.0, 1.95/) ! for non-urban, max kt, min urban kt, used in chm_vert_diff_s
   chm_urban_abl_min    = 484.0
   chm_timings_l        = .false.
   chm_debug_trace_l    = .false.
   chm_messy_jval_l     = .false.
   chm_vit_l            = .false.
   chm_direct_l         = .false.
   chm_indirect_l       = .false.
   chm_canopy_shading_l = .false.
   chm_dep_lai2d_l      = .false.
   chm_sat_seasons_l    = .false.
   chm_cffeps_online_l  = .false.
   wf_case              = 0
   chm_pblh_min_l       = .false.
   chm_ae_spread_l      = .false.
   chm_ammonia_bidi_s   = 'OFF'
   chm_ammonia_gep      = (/ 1000.0, 1000.0, 1000.0, 500.0, 1000.0, 2000.0, 800.0, 0.0, 20.0, 100.0, 20.0, 0.0, 0.0, 0.0, 0.0 /)
!
   aerosize             = -1.0
   chm_bkgd_ch4         = 1.8   ! in ppm unit
   chm_bkgd_co2         = 400.0 ! in ppm unit
   gas_ppb_out_list_s   = '    '
   gas_ppb_out_list_s(1)= 'LIST'
   aero_opt_wavel       = -1.0
!
   dbg_lat              = -99.0
   dbg_lon              = -999.0
   dbg_lev              = -1
   dbg_step             = -1
   dbg_itr              = -1
   nk_start             =  1
   nk_start_pm          =  1
!
!  Set debug message
!
   local_dbg = (.false. .or. global_debug)

!  Return +1 if chm_nml works properly, else return -1 (start with the assumption of failure)
!
   chm_nml  = -1
   unf       = 0
!
!  Attempt to open gem_settings.nml file -> ABORT if cannot
!
   err_open = fnom (unf, F_namelist, 'SEQ+OLD' , nrec)
   if ((err_open.ne.0) .and. (chm_lun_out>0)) write (chm_lun_out, 1000) trim(F_namelist)
   call chm_stop ('CHM_NML - ERROR OPENING NML FILE',err_open)
!
!  Attempt a read
!
   read (unf, nml = chemistry_cfgs, iostat = read_status)
   if (read_status > 0) then
      if (chm_lun_out > 0) then
         write(chm_lun_out, *) ' Error in reading namelist chem_cfgs (s/r chm_nml)'
         write(chm_lun_out, *) ' Read status = ', read_status
         backspace(unf)
         read(unf, '(A)') line
         write(chm_lun_out, *) 'Invalid line in chemistry_cfgs namelist is: ', trim(line)
      endif
      err = fclos (unf)
      write (0, *) 'CHM_NML - FAULTY CHEMISTRY NAMELIST -> PLEASE REVISE '
      chm_nml     = -1
      return
   else if (read_status < 0) then
      write (0, *) ' No chemistry_cfgs found in gem_settings --> ABORT ALL!!!!'
      chm_nml     = -1
      err = fclos (unf)
      return
   endif
!
!  Detect Master switch. If false, NORMAL EXIT WITH MESSAGE
!
   if (.not. chm_master) then
    ! Ensure that setting that have effects outside chemistry are also disabled
      chm_direct_l   = .false.
      chm_indirect_l = .false.
      if (chm_lun_out > 0) write(chm_lun_out, 1010)
      err = fclos (unf)
      chm_nml = 1
      return
   end if
!
! Both online and standalone CFFEPS emissions processing are not supposed to run
! at the same time; if both are set in the input namelist, the online instance
! is preferred and the stand-alone instance is disabled
   if (chm_cffeps_online_l .and. chm_get_wf_emis_l) then
      chm_get_wf_emis_l = .false.
   end if

   if (.not.chm_emis_master) then
      if (chm_lun_out > 0) then
         if (chm_get_ae_emis_l .or. chm_get_be_emis_l .or. chm_get_mj_emis_l .or. chm_get_wf_emis_l) then
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> Emission master switch set to false '
            write(chm_lun_out, *) '> chm_get_ae_emis_l, chm_get_be_emis_l, chm_get_mj_emis_l and chm_get_wf_emis_l set to false'
         end if
      end if
      chm_get_ae_emis_l = .false.
      chm_get_be_emis_l = .false.
      chm_get_mj_emis_l = .false.
      chm_get_wf_emis_l = .false.
      chm_get_on_emis_l = .false.
   else
      if (.not.chm_get_ae_emis_l .and. .not.chm_get_be_emis_l .and. &
          .not.chm_get_mj_emis_l .and. .not.chm_get_wf_emis_l .and. &
          .not.chm_get_on_emis_l) &
         chm_emis_master = .false.
   endif
!
! Consistency check according to selected model (NIL, MACH, LINOZ, etc)
!
   istatus = 0
!
   select case (chm_mono_s)
     case ('CLIP', 'ILMC', 'NIL')
      continue
     case default
      write(0, *) '### Error in chm_nml ###'
      write(0, *) '# CHM_MONO_S unknown: ', chm_mono_s
      write(0, *) '###         ABORT         ###'
      istatus = min(istatus,-1)
   end select

   select case (chm_mass_s)
     case ('BC', 'NIL')
      continue
     case default
      write(0, *) '### Error in chm_nml ###'
      write(0, *) '# CHM_MASS_S unknown: ', chm_mass_s
      write(0, *) '###         ABORT         ###'
      istatus = min(istatus,-1)
   end select

   if ((chm_mono_s == 'NIL') .and. (chm_mass_s == 'NIL')) then
      write(0, *) '###              Warning in chm_nml              ###'
      write(0, *) '# When CHM_MASS_S is nil, monotonicity is applied unless adw_mono_L is set to false in adw_cfgs'
      write(0, *) '###                                              ###'
   endif
!
   model: select case (chm_model_s)
   case ('NIL') model
      if (chm_lun_out > 0) then
         write(chm_lun_out, *) '> Warning'
         write(chm_lun_out, *) '> No chemistry package selected'
      endif
      err = fclos (unf)
      chm_nml    = 1
      return
   case ('MACH') model
!
!     MACH MODEL
!     =====================================================================
!
!     Gas chemistry
!     =============
      if (chm_pkg_gas_s == 'NIL') then
         chm_biog_s          = 'NIL'      !Override chm_biog_s
         chm_soa_s           = 'NIL'      !Override chm_soa_s
         chm_strato_s        = 'NIL'      !Override chm_strato_s
         chm_aqhi_l          = .false.    !Override chm_aqhi_l
         chm_messy_jval_l    = .false.    !Override chm_messy_jval_l
         chm_diag_accum_l    = .false.    !Override chm_diag_accum_l
         chm_diag_colum_l    = .false.    !Override chm_diag_colum_l
         if (chm_lun_out > 0) then
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> No gas package selected'
            write(chm_lun_out, *) '> chm_biog_s = NIL'
            write(chm_lun_out, *) '> chm_soa_s = NIL'
            write(chm_lun_out, *) '> chm_strato_s = NIL'
            write(chm_lun_out, *) '> chm_aqhi_l = .false.'
            write(chm_lun_out, *) '> chm_messy_jval_l = .false.'
            write(chm_lun_out, *) '> chm_diag_accum_l = .false.'
            write(chm_lun_out, *) '> chm_diag_colum_l = .false.'
            if (chm_pkg_pm_s /=  'NIL') then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> chm_pkg_pm_s requires selected gas species '
               write(chm_lun_out, *) '> these gas fields will be created but gas chemistry will not be performed '
            end if
         end if
      end if
! Ensure chm_messy_jval_l is set for all SAPRC mechanisms
      if ((chm_pkg_gas_s(1:5) == 'SAPRC') .and. (.not. chm_messy_jval_l)) then
         write(chm_lun_out, *) '> Warning'
         write(chm_lun_out, *) "> MESSy module is required for chm_pkg_gas_s(1:5)='SAPRC'"
         write(chm_lun_out, *) "> reset chm_messy_jval_l='.TRUE.'"
         chm_messy_jval_l = .true.
      end if
! Set default PPB gas output list, if unset
      if (gas_ppb_out_list_s(1) == 'LIST') then
          gas_ppb_out_list_s(1) = 'O3  '; gas_ppb_out_list_s(2) = 'N2  '
          gas_ppb_out_list_s(3) = 'NO  '; gas_ppb_out_list_s(4) = 'S2  '
      end if

!     Aerosol chemistry
!     =================
      if ((chm_pkg_pm_s == 'NIL') .or. (chm_pkg_pm_s == 'GOCART_SO2')) then
         chm_soa_s            = 'NIL'      !Override chm_soa_s
         chm_seaflux_s        = 'NIL'      !Override chm_seaflux_s
         chm_winddust_s       = 'NIL'      !Override chm_winddust_s
         chm_hetchem_s        = 'NIL'      !Override chm_hetchem_s
         chm_pm_drydep_s      = 'NIL'      !Override chm_pm_drydep_s
         chm_aqueous_s        = 'NIL'      !Override chm_aqueous_s
         chm_met_modulation_s = 'NIL'      !Override chm_met_modulation_s
         chm_aqhi_l           = .false.    !Override chm_aqhi_l
         chm_blcld_l          = .false.    !Override chm_blcld_l
         chm_diag_wetdep_l    = .false.    !Override chm_diag_wetdp_l
         chm_diag_aerosols_l  = .false.    !Override chm_diag_aerosols_l
         chm_diag_aero_opt_l  = .false.    !Override chm_diag_aero_opt_l
         chm_direct_l         = .false.    !Override chm_direct_l
         chm_indirect_l       = .false.    !Override chm_indirect_l
         if (chm_lun_out > 0) then
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> No PM package selected'
            write(chm_lun_out, *) '> chm_soa_s = NIL, chm_seaflux_s = NIL, chm_winddust_s = NIL'
            write(chm_lun_out, *) '> chm_hetchem_s = NIL, chm_pm_drydep_s = NIL, chm_aqueous_s = NIL, chm_met_modulation_s = NIL'
            write(chm_lun_out, *) '> chm_aqhi_l = .false. , chm_blcld_l = .false. , chm_diag_wetdp_l = .false. '
            write(chm_lun_out, *) '> chm_diag_aero_opt_l = .false. , chm_direct_l = .false., chm_indirect_l = .false. '
         end if
      else
         if (chm_pkg_pm_s == 'CAM2BINS') then
            isize = 2
            icob  = 2
            nbnd  = isize * chm_intrsec_ndiv
            incld_dtmin = 1.0e-6
            itermax = 10
         else if (chm_pkg_pm_s == 'CAM12BINS') then
            isize = 12
            icob  = 12
            nbnd  = isize
            incld_dtmin = 1.0e-8
            itermax = 20
         else
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# Aerosol package unknown: ', chm_pkg_pm_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
         end if

         if (chm_intrsec_ndiv < 1) then
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# chm_intrsec_ndiv must be positive and greater than zero'
            write(0, *) '# you choose: ', chm_intrsec_ndiv
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
         end if
         if (chm_intrsec_ver < 1 .or. chm_intrsec_ver > 6) then
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# chm_intrsec_ver must be between 1 and 6'
            write(0, *) '# you choose: ', chm_intrsec_ver
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
         end if

!        Dry deposition for PM
!        =====================
         select case (chm_pm_drydep_s)
            case ('ZHANG', 'ZHANG_MAKAR')
               continue
            case ('NIL')
               if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> No PM dry deposition package selected'
               endif
            case default
               write(0, *) '### Error in chm_nml ###'
               write(0, *) '# PM dry deposition package unknown: ', chm_pm_drydep_s
               write(0, *) '###         ABORT         ###'
               istatus = min(istatus,-1)
         end select

!        Heterogeneous chemistry
!        =======================
         select case (chm_hetchem_s)
           case ('HETV')
              if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning for HETV'
                  write(chm_lun_out, *) '> chm_hetchem_metstlb_l =', chm_hetchem_metstlb_l
              endif
              continue
           case ('NIL')
              if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> No heterogeneous chemistry package selected'
              endif
           case default
               write(0, *) '### Error in chm_nml ###'
               write(0, *) '# Heterogeneous chemistry package unknown: ', chm_hetchem_s
               write(0, *) '###         ABORT         ###'
               istatus = min(istatus,-1)
         end select

!        Sea salt
!        ========
         select case (chm_seaflux_s)
           case ('GONG_MONAHAN', 'GONG_MONAHAN_F', 'SMITH')
              continue
           case ('NIL')
              if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> No sea salt flux treatment'
              endif
           case default
               write(0, *) '### Error in chm_nml ###'
               write(0, *) '# Sea salt flux package unknown: ', chm_seaflux_s
               write(0, *) '###         ABORT         ###'
               istatus = min(istatus,-1)
         end select

!        Wind dust
!        =========
         select case (chm_winddust_s)
           case ('CAM_WINDDUST')
              if (chm_lun_out > 0 .and. chm_pkg_pm_s == 'CAM2BINS') then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> CAM_WINDDUST inactive when chm_pkg_pm_s == CAM2BINS'
              endif
              continue
           case ('NIL')
              if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> No wind dust flux treatment'
              endif
           case default
              if (chm_pkg_pm_s == 'CAM2BINS' .or. chm_pkg_pm_s == 'CAM12BINS') then
                  if (chm_lun_out > 0) then
                      write(chm_lun_out, *) '> Warning'
                      write(chm_lun_out, *) '> CAM_WINDDUST inactive when chm_pkg_pm_s == CAM2BINS or CAM12BINS'
                  endif
              else
                  write(0, *) '### Error in chm_nml ###'
                  write(0, *) '# Wind dust flux package unknown: ', chm_winddust_s
                  write(0, *) '###         ABORT         ###'
                  istatus = min(istatus,-1)
              endif
         end select

!        Aqueous chemistry
!        =================
         select case (chm_aqueous_s)
           case ('GONG')
              continue
           case ('NIL')
              if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> No Aqueous phase chemistry package selected'
              endif
           case default
               write(0, *) '### Error in chm_nml ###'
               write(0, *) '# Aqueous phase chemistry package unknown: ', chm_aqueous_s
               write(0, *) '###         ABORT         ###'
               istatus = min(istatus,-1)
         end select

!        Fugitive dust
!        =============
         select case (chm_met_modulation_s)
           case ('ON', 'CM_ONLY')
              continue
           case ('NIL', 'OFF')
              if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '> Warning'
                  write(chm_lun_out, *) '> No meteorological dust modulation applied'
              end if
              chm_met_modulation_s = 'NIL'
           case default
              write(0, *) '### Error in chm_nml ###'
              write(0, *) '# Fugitive dust modulation method unknown: ', chm_met_modulation_s
              write(0, *) '###         ABORT         ###'
              istatus = min(istatus,-1)
         end select
         if ((.not. chm_get_fd_emis_l) .and. (chm_met_modulation_s(1:2) == 'ON')) then
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '### Error in chm_nml ###'
               write(chm_lun_out, *) '# Fugitive dust emissions must be enabled for the following setting'
               write(chm_lun_out, *) '# chm_met_modulation_s = "ON"'
               write(chm_lun_out, *) '###         ABORT         ###'
            end if
            istatus = min(istatus,-1)
         end if
         if ((chm_get_fd_emis_l) .and. (chm_met_modulation_s(1:7) == 'CM_ONLY')) then
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '### Error in chm_nml ###'
               write(chm_lun_out, *) '# chm_met_modulation_s == "CM_ONLY" is only allowed for '
               write(chm_lun_out, *) '# backward compatibility when chm_get_fd_emis_l = .false.'
               write(chm_lun_out, *) '###         ABORT         ###'
            end if
            istatus = min(istatus,-1)
         end if

!        Wet deposition for PM
!        =====================
         if (chm_wdep_scav_coef_s == 'WANG2014') then
            chm_split_snow_rain_l = .true.
         end if

!     Diagnostics for aerosol output in ug/m3
!     =============================================================
         if (chm_diag_aerosols_L) then
            if (.not.((chm_pkg_pm_s == 'CAM2BINS') .or. (chm_pkg_pm_s == 'CAM12BINS'))) then
               if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '### Error in chm_nml ###'
                  write(chm_lun_out, *) '# Diagnostic package chm_diag_aerosols_L= ', chm_diag_aerosols_L, 'must be set to .false. because chm_pkg_pm_s is neither CAM2BINS nor CAM12BINS'
                  write(chm_lun_out, *) '###         ABORT         ###'
                  istatus = min(istatus,-1)
               end if
            end if
         end if

!     Diagnostics for aerosol optical properties
!     =============================================================
         nwl_aod = count(aero_opt_wavel > 0.0)
         if (chm_direct_l) chm_diag_aero_opt_l = .true.
         if (chm_diag_aero_opt_l) then
            if (all(aero_opt_wavel < 0.0)) then
               if (chm_lun_out > 0) then
                  write(chm_lun_out, *) '### Error in chm_nml ###'
                  write(chm_lun_out, *) '# Diagnostic package chm_diag_aero_opt_l= ', chm_diag_aero_opt_l
                  write(chm_lun_out, *) '# but no mid-band wavelengths specified for evaluating the aerosol optical properties '
                  write(chm_lun_out, *) '###         ABORT         ###'
               end if
               istatus = min(istatus,-1)
            end if
         end if
         if (chm_direct_l .and. nwl_aod < 4) then
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '### Error in chm_nml ###'
               write(chm_lun_out, *) '# Mid-band wavelengths for aerosol optical properties= ', aero_opt_wavel
               write(chm_lun_out, *) '# At least 4 wavelength used for short-wave raditive transfer '
               write(chm_lun_out, *) '# used in radriv s/r are required for the direct feedback effect'
               write(chm_lun_out, *) '###         ABORT         ###'
            end if
            istatus = min(istatus,-1)
         end if
!
         if (chm_indirect_l .and. stcond(1:2) /= "MP") then
            chm_indirect_l = .false.    !Override chm_indirect_l
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> No microphysics condensation scheme selected'
               write(chm_lun_out, *) '> chm_indirect_l = .false. '
            end if
         end if

      end if ! End of if on value of chm_pkg_pm_s


!     Biogenic chemistry package
!     ==========================
      if (.not.chm_get_be_emis_l) then
        if (chm_lun_out > 0) then
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> No biogenic emissions read, biogenic chemistry package set to NIL'
        endif
        chm_biog_s = 'NIL'
      endif
      select case (chm_biog_s)
         case ('BEIS3.09')
            continue
         case ('NIL')
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '> No biogenic package selected'
            endif
         case default
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# Biogenic emission package unknown: ', chm_biog_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
      end select

!     Ammonia bidirectional flux
!     ==========================
      select case (chm_ammonia_bidi_s)
         case ('GEP', 'GEP2D')
            write(chm_lun_out, *) '> Ammonia bidirectional flux enabled'
            write(chm_lun_out, *) '> '
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> The bidirectional flux models ammonia emissions through'
            write(chm_lun_out, *) '> the specification of emission potentials. If the emission'
            write(chm_lun_out, *) '> potentials are used to model an emissions source that is'
            write(chm_lun_out, *) '> already specified in the input emissions files, then these'
            write(chm_lun_out, *) '> emissions will need to be removed from the input emissions'
            write(chm_lun_out, *) '> files to avoid double counting the emissions from this source.'
         case ('OFF')
            continue
         case default
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# CHM_AMMONIA_BIDI_S unknown: ', chm_ammonia_bidi_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
      end select

!     Dry deposition for gas
!     ======================
      select case (chm_gas_drydep_s)
         case ('ROBICHAUD','ROBICHAUD2','ROBICHAUD3')
            continue
         case ('NIL')
            chm_diag_drydep_l    = .false.
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> No dry deposition package selected'
               write(chm_lun_out, *) '> chm_diag_drydep_l set to  .false. '
            endif
            if (trim(chm_ammonia_bidi_s) /= 'OFF') then
               write(0, *) '### Error in chm_nml ###'
               write(0, *) '# Ammonia bidirectional flux requires a dry deposition scheme.'
               write(0, *) '###         ABORT         ###'
               istatus = min(istatus,-1)
            end if
         case default
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# Dry deposition package unknown: ', chm_gas_drydep_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
      end select

!     SOA chemistry package
!     =====================
      select case (chm_soa_s)
        case ('JIANG', 'ODUM')
           continue
        case ('NIL')
           if (chm_lun_out > 0) then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> No SOA chemistry package selected'
           endif
        case default
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# SOA chemistry package unknown: ', chm_soa_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
      end select

!
!     Stratospheric O3 chemistry
!     ==========================
      if (chm_strato_s == 'LINOZ') then
         if (chm_Lun_out > 0) write(chm_lun_out,*)'LINOZ set for the stratospheric O3 package'

      else if (chm_strato_s == 'NIL') then
         if (chm_lun_out > 0) then
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> No stratospheric O3 package selected'
         end if
      end if
!
!     Vertical diffusion
!     ==================
      select case (chm_vert_diff_s)
         case ('RPNPHY', 'RPNPHY_I', 'RPNPHY_U')
            continue
         case ('FLUX', 'BOUNDARY')
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# Vertical diffusion algorithm : ', chm_vert_diff_s
            write(0, *) '# not yet supported with GEM-MACH version 2'
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
         case ('NIL')
            if (chm_lun_out > 0) then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> No diffusion package selected'
               write(chm_lun_out, *) '> Area emissions will not be treated '
            endif
         case default
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# Vertical diffusion algorithm unknown: ', chm_vert_diff_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
      end select

!     Wildfire emissions treatment
!     ============================
      if (.not. chm_get_wf_emis_l) then
        if (wf_case > 3 .and. chm_get_mj_emis_l) then
           write(0, *) '###    Error in chm_nml    ###'
           write(0, *) '# CFFEPS fire emissions plume-rise parameters are'
           write(0, *) '# required for this plume-rise treatment wf_case:', wf_case
           write(0, *) '###          ABORT         ###'
           istatus = min(istatus,-1)
        end if
      else
        if (wf_case /= 4 .and. (chm_lun_out > 0)) then
           write(chm_lun_out, *) '> Warning'
           write(chm_lun_out, *) '> CFFEPS fire emissions processing is currently '
           write(chm_lun_out, *) '> only viable with plume-rise treatment wf_case=4'
           write(chm_lun_out, *) '> Resetting wf_case to 4 from initial wf_case=', wf_case
        end if
        wf_case = 4
      endif

!     Plumerise
!     =========
      if (.not. (chm_get_mj_emis_l .or. chm_get_wf_emis_l)) then
        if (chm_lun_out > 0) then
            write(chm_lun_out, *) '> Warning'
            write(chm_lun_out, *) '> No major point emissions and no CFFEPS wildfire emissions read, plumerise package set to NIL'
        endif
        chm_mj_treatment_s = 'NIL'
        chm_do_mjpts_l = .false.
      endif
      select case (chm_mj_treatment_s)
        case ('PLUMERISE', 'PLUMERISE2')
           continue
        case ('NIL')
           if (chm_get_mj_emis_l .and. (chm_lun_out > 0)) then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> Major point emissions read but no plumerise '
           endif
           if (chm_get_wf_emis_l .and. (chm_lun_out > 0)) then
               write(chm_lun_out, *) '> Warning'
               write(chm_lun_out, *) '> CFFEPS wildfire emissions read but no plumerise '
           endif
        case default
            write(0, *) '### Error in chm_nml ###'
            write(0, *) '# Major points treatment package unknown: ', chm_mj_treatment_s
            write(0, *) '###         ABORT         ###'
            istatus = min(istatus,-1)
      end select
!
!     Vehicle induced turbulence
!     =========
      if (chm_vit_l .and. (.not. chm_get_on_emis_l)) then
         if (chm_lun_out > 0) then
            write(chm_lun_out, *) '### Error in chm_nml ###'
            write(chm_lun_out, *) '# The VIT scheme requires on-road area emissions'
            write(chm_lun_out, *) '# chm_get_on_emis_l needs to be set to .true.'
            write(chm_lun_out, *) '###         ABORT         ###'
         end if
         istatus = min(istatus,-1)
      end if
!
!     Online CFFEPS - Read the CFFEPS_inputs namelist
!     =========
      if (chm_cffeps_online_l) then
         istatus = min(mach_cffeps_nml(unf), 1)
      end if

   case default model
!
!     INVALID MODEL NAME
!     =====================================================================
!
      write(0, *) '### Error in chm_nml ###'
      write(0, *) '# Chemistry package unknown: ', chm_model_s
      write(0, *) '###         ABORT         ###'
      istatus = min(istatus,-1)
   end select model

   call chm_stop('chm_nml',istatus)

! Set CHM debug_trace if global_debug is true.
   if (global_debug) chm_debug_trace_l = .true.

   if (chm_lun_out > 0) then
      write(chm_lun_out, *) 'Final read for chemistry_cfgs namelist:'
      write(chm_lun_out, nml = chemistry_cfgs)
   end if

!
!  A chemistry namelist was read in successfully.
!  Close namelist file
!
   err = fclos(unf)
   chm_nml    = 1

1000 format (/,' NAMELIST FILE ',a,' NOT AVAILABLE - ERROR->ABORT')
1010 format (/,' CHM_NML: MASTER CHEMICAL SWITCH OFF ---> ALL CHEMISTRY SHUT OFF!!!')

end function chm_nml
