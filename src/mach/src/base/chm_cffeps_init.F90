!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2018 - Air Quality Research Division &
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
! Projet / Project : GEM-MACH
! Fichier / File   : mach_cffeps_init.ftn90
! Creation         : A. Akingunola, J. Chen, P. Makar, and Kerry Anderson
!                     - Fall 2018
! Description      : Read the input namelist and initialize CFFEPS
!
!
!!if_on
integer function chm_cffeps_init(F_basedir_S, my_pe, numproc, rpi0_j0)
!!if_off
   use mach_cffeps_mod
   use mach_cffeps_speciations_mod, only: me_species_gas, me_species_pm,   &
                                          me_species_voc, fire_me_species, &
                                          nb_fire_me_species, ngas, nvocs, npm
   use chm_utils_mod,           only: chm_timestep, global_debug
   use chm_mjrpts_sortinfo_mod, only: me_species_index, nb_me_species
   use chm_species_info_mod,    only: species_master
   use chm_nml_mod,             only: chm_do_cffeps_l, chm_step_factor
   use chm_datime_mod,          only: secondsin1hour
   use chm_consphychm_mod,      only: pi
   use mach_cam_utils_mod,      only: isize
   use phygridmap,              only: mapmod2phy, phy_lcl_gid, phy_glb_gid, &
                                      phy_lcl_ni, phydim_ni
   use ezgrid_mod,              only: ezgrid_find_ij0
   use rpn_comm_itf_mod

   implicit none
!!if_on
   character(len=*),                       intent(in) :: F_basedir_S
   integer(kind=4),                        intent(in) :: my_pe, numproc
   integer(kind=4), dimension(5, numproc), intent(in) :: rpi0_j0
!!if_off

   integer(kind=4)             :: funit, ier, i, il, jk, int_x, int_y
   integer(kind=4)             :: ij(2), i0, j0, lni, lnj
   integer(kind=4)             :: lx_lbound, lx_ubound, ly_lbound, ly_ubound
   real(kind=4)                :: dt

   character(len=1) :: zbin
   integer(kind=4), parameter :: imax = 24
   integer(kind=4)  :: nhotspots, nsp, nbin, sn, emis_nb_bins, max_hours, lc_hotspots2

   character(len=8) :: fuel     ! Fueltype
   integer(kind=4)  :: ifuel    ! Fueltype
   real(kind=4)     :: dtime    ! detection time of fire HHMM (time all values are assumed to be collected)
   real(kind=4)     :: djulian  ! detection Julian Day
   real(kind=4)     :: ffmc     ! FFMC
   real(kind=4)     :: dmc      ! Duff Moisture Code
   real(kind=4)     :: dc       ! Drought Code
   real(kind=4)     :: bui      ! BUI
   real(kind=4)     :: fmc
   real(kind=4)     :: curing, gs, aspect, pc, pdf, cbh, gfl, lfl, dfl,   &
                       sfl, cfl, ffl, wfl, sh, sd, rsmk, zplm, area
   real(kind=4)     :: estarea  ! Fire size at time of detection
   real(kind=4)     :: q0       ! Fire energy at model start time
   real(kind=4)     :: total_emis  ! Total emissions at model start
   real(kind=4)   , dimension(imax)  :: qs       ! Plume energy over time
   real(kind=4)   , dimension(imax)  :: mws      ! Plume moisture content
   real(kind=4)   , dimension(imax)  :: emis_f, emis_s, emis_r
   real(kind=4)   , dimension(:), allocatable :: lat      ! Latitude [decimal degrees]
   real(kind=4)   , dimension(:), allocatable :: lon      ! Longitude [decimal degrees]
   real(kind=4)   , dimension(:), allocatable :: lat2x, lon2y
   integer(kind=4), dimension(:), allocatable :: pe_info

   real(kind=4)   , dimension(:,:), allocatable :: fire_info, emis_info
   type(rpncomm_context)                      :: mjr_context
!  Declaration of external functions and subroutines
   integer(kind=4), external :: gdxyfll
   external physeterror
   logical(kind=4)                            :: local_dbg

   local_dbg = (.false. .or. global_debug)
!
   chm_cffeps_init = -1

! Initialize weight; weight(17) adjusted to make the total = 1
   weight =  (/0.015648395, 0.012470576, 0.010137905, 0.008443153, &
               0.004159725, 0.003641906, 0.004455177, 0.005716244, &
               0.017972534, 0.032938579, 0.043543306, 0.055667556, &
               0.065778522, 0.078207403, 0.090168462, 0.104189127, &
               0.101294339, 0.091901405, 0.072798822, 0.057406612, &
               0.044004006, 0.033647664, 0.025825035, 0.019983547/)

! Over-write the hotspotfile (again)
   hotspotfile = trim(F_basedir_S)//'/MODEL_INPUT/cffeps_hotspots.csv'

   cffeps_alpha = cffeps_alpha * pi / 180.0    !/* convert to radians */

   max_hours = 24 ! Maximum number of hours for future distribution of fire energy

   select case (trim(cffeps_method))
      case ("Alberta", "AB")
        cffeps_method = "AB"
      case ("cmc", "firework")
        cffeps_method = "cmc"
      case ("bigfoot", "bf")
        cffeps_method = "bigfoot"
      case ("standard", "icao")
        cffeps_method = "icao"
   end select

   select case (trim(cffeps_fire_shape))
      case ("line", "wedge", "perimeter")
        cffeps_fire_shape = "line"
      case ("tophat", "persistence")
        cffeps_fire_shape = "tophat"
   end select

   select case (trim(cffeps_diurnal))
      case ("off", "false", "no", "0")
        cffeps_diurnal = "OFF"
      case ("cvw", "hr", "hourly")
        cffeps_diurnal = "CVW"
      case ("bdl", "lawson", "tables")
        cffeps_diurnal = "BDL"
      case ("on", "emc", "2")
        cffeps_diurnal = "EMC"
   end select

   if (trim(cffeps_fire_type) == "piecewise") cffeps_fire_type = "dry"

   if (cffeps_thstart > 24. .or. cffeps_thend > 24.) then
      cffeps_thstart = (int(cffeps_thstart) / 100) + (mod(int(cffeps_thstart), 100) / 60.)
      cffeps_thend = (int(cffeps_thend) / 100) + (mod(int(cffeps_thend), 100) / 60.)
   end if

   if (cffeps_thstart >= cffeps_thend) then
      write(*,*)"Warning: Top-hat start ", cffeps_thstart, " >= Top-hat end ", cffeps_thend
      cffeps_thstart = 9.0    !/* time to start (exclusive) top-hat fire growth [decimal hours LST] */
      cffeps_thend = 21.0     !/* time to end (inclusive) top-hat fire growth [decimal hours LST] */
   end if

   cffeps_timestep = 1.0 !chm_timestep / sngl(secondsin1hour) [In unit of hours]

   max_timesteps = nint(real(max_hours) / cffeps_timestep)
   ! How often (in model steps) do we run CFFEPS to generate emissions
   cffeps_steps = nint(sngl(secondsin1hour) / chm_timestep) *  chm_step_factor

   if (cffeps_reset > 24.) &
      cffeps_reset = (int(cffeps_reset) / 100) + (mod(int(cffeps_reset), 100) / 60.)

   if (cffeps_ldt == 100 .or. cffeps_ldt == 1) then
      cffeps_ldt = 1
   else
      cffeps_ldt = 0
   end if
!
!  Assign weight array values to top-hat approach (if used)
   if (trim(cffeps_fire_shape) == "top-hat") then
      do i = 0, 23
         if (i >= floor(cffeps_thstart) .and. i < ceiling(cffeps_thend)) then
            if ((cffeps_thend - cffeps_thstart) < 1.0) then
               dt = cffeps_thend - cffeps_thstart
            else if ((i < cffeps_thstart) .and. &
                    (cffeps_thstart - real(i)) < 1.0) then
                dt = real(i) + 1.0 - cffeps_thstart
            else if (cffeps_thend - real(i) < 1.0) then
                dt = cffeps_thend - real(i)
            else
                dt = 1.0
            end if
            weight(i + 1) = dt / (cffeps_thend - cffeps_thstart)
         else
            weight(i + 1) = 0.0
         end if
      end do
   end if
!
   if (cffeps_diurnal == "OFF") then ! /* no diurnal effect -- 2020-06-04 */
      weight = 1.0 / 24.0
   end if

   if (cffeps_fire_shape(1:4) == "line") then
      cffeps_fbp_accel = .false.
   end if
!
   nhotspots = 0
   if (my_pe == 0) then
      funit = 11
      open(unit=funit, file=trim(hotspotfile), status='old', action='read')

       ! First determine the number of fire hotspots
      read(funit, *)  ! The header
      do
        read(funit, *, end=45)
        nhotspots = nhotspots + 1
      end do

 45   allocate(lat(nhotspots), lon(nhotspots),     &
               lat2x(nhotspots), lon2y(nhotspots), &
               fire_info(nhotspots, fire_parameters + max_timesteps), &
               emis_info(nhotspots, 4 * max_timesteps))

    ! Now read the hotspots' data
      rewind(funit)

      read(funit, *) ! First read the header
      qs = 0.0   ! Temporarily
      do i = 1, nhotspots
         read(funit, 30) lat(i), lon(i), dtime, djulian, fuel, estarea,       &
                         area, ffmc, fmc, dmc, dc, curing, gs, aspect, pc,    &
                         pdf, cbh, gfl, lfl, dfl, sfl, cfl, ffl, wfl, sh, sd, &
                         q0, total_emis, rsmk, zplm, (qs(jk), jk = 1, imax),  &
                         (mws(jk), jk = 1, imax),                &
                         (emis_f(jk), emis_s(jk), emis_r(jk), jk = 1, imax)
!
! Reassign fuel type (if required) and setup to appropriate index.
!      If there is no matching fuel type, set to NF. */
!      forall(j = 1:fuel_size) fuel(j:j) = clib_toupper(fuel(j:j))
         fuel = adjustl(fuel)
         if (fuel(1:2) == "C1") then
            ifuel = 1
         else if (fuel(1:2) == "C2") then
            ifuel = 2
         else if (fuel(1:2) == "C3") then
            ifuel = 3
         else if (fuel(1:2) == "C4") then
            ifuel = 4
         else if (fuel(1:2) == "C5") then
            ifuel = 5
         else if (fuel(1:2) == "C6") then
            ifuel = 6
         else if (fuel(1:2) == "C7") then
            ifuel = 7
         else if (fuel(1:2) == "D2" .or. fuel(1:2) == "D1") then
            ifuel = 8
         else if (fuel(1:2) == "M1") then
            ifuel = 10
         else if (fuel(1:2) == "M2") then
            ifuel = 11
         else if (fuel(1:2) == "M3") then
            ifuel = 12
         else if (fuel(1:2) == "M4") then
            ifuel = 13
         else if (fuel(1:2) == "S1") then
            ifuel = 14
         else if (fuel(1:2) == "S2") then
            ifuel = 15
         else if (fuel(1:2) == "S3") then
            ifuel = 16
         else if (trim(fuel) == "CROPLAND" .or. trim(fuel) == "LOW_VEG" .or. &
                  trim(fuel) == "O1A" .or. trim(fuel) == "O1a") then
            ifuel = 17
         else if (fuel(1:2) == "O1") then
            ifuel = 18
         else if (fuel(1:2) == "WA") then
            ifuel = 19
         else
        ! ("URBAN", "BOG", "WATER", "NON-FUEL") and anything else is undefined
            ifuel = 20
         end if
         if (ifuel > 18) ifuel = 20

!/* Calculate the BUI value from DMC and DC values.  Determine whether the BUI
!      effect is being used. (Excerpt of BUICalc in Fbp sub-module) */
             !// updated fix 2015-11-09 KRA
         if ((dmc * dc) == 0.0) then
            bui = 0.0
         else
            if (dmc <= (0.4 * dc)) then       ! /* 27a */
               bui = 0.8 * dmc * dc / (dmc + 0.4 * dc)
            else                              ! /* 27b */
               bui = dmc - (1.0 - 0.8 * dc / (dmc + 0.4 * dc)) *    &
                        (0.92 + ((0.0114 * dmc)**1.7))
            end if
            if (bui < 0.0 .or. bui > 1000.0) &
               bui = 0.0    !/* This turns off BUI effect */
         end if

         fire_info(i, i_fuel)  = real(ifuel) ! Fuel type
         fire_info(i, i_dj)    = djulian     ! detection Julian date
         fire_info(i, i_dtime) = dtime       ! detection time of fire HHMM
         fire_info(i, i_bui)   = bui         ! BUI
         fire_info(i, i_area0) = estarea     ! Area of fire at detection
         fire_info(i, i_ffmc)  = ffmc        ! FFMC
         fire_info(i, i_dmc)   = dmc         ! Duff Moisture Code
         fire_info(i, i_dc)    = dc          ! Drought Code
         fire_info(i, i_fmc)   = fmc         ! Foliar Moisture Content
         fire_info(i, i_area)  = area        ! Fire hotspot size (ha)
         fire_info(i, i_cfl)   = cfl        ! Crown Fuel Load [kg/m^2]
         fire_info(i, i_sfl)   = sfl        ! Surface Fuel Load [kg/m^2]
         fire_info(i, i_ffl)   = ffl        ! Fine Fuel Load [kg/m^2]
         fire_info(i, i_wfl)   = wfl        ! Woody Fuel Load [kg/m^2]

         fire_info(i, i_cured) = curing     ! Percent Cured for O1a/O1b
         fire_info(i, i_pc)    = pc         ! Percent Confier for M1/M2
         fire_info(i, i_pdf)   = pdf        ! Percent Dead Fir for M3/M4
         fire_info(i, i_cbh)   = cbh        ! Crown to Base Height [m]
         fire_info(i, i_asp)   = aspect     ! Aspect [degrees]
         fire_info(i, i_gfl)   = gfl        ! Grass Fuel Load [kg/m^2]
         fire_info(i, i_gs)    = gs         ! Slope [percent]
         fire_info(i, i_sd)    = sd         ! C6 Stand Density [stems/ha]
         fire_info(i, i_sh)    = sh         ! C6 Stand Height [m]
         fire_info(i, i_lfl)   = lfl        ! Litter Fuel Load [kg/m^2]
         fire_info(i, i_dfl)   = dfl        ! Duff (Fermentation) Fuel Load [kg/m^2]

         fire_info(i, i_q0)    = q0          ! amount of energy previously injected into atmosphere (J)
         fire_info(i, i_emis)  = total_emis  ! total emissions from the fire (tonnes)
         fire_info(i, i_rsmk)  = rsmk
         fire_info(i, i_zplm)  = zplm
         do jk = 1, imax
            il = jk + fire_parameters
            fire_info(i, il)   = qs(jk)      ! Plume energy over 24hrs
            emis_info(i, jk)   = mws(jk)     ! Plume moisture content
            emis_info(i, jk + imax)   = emis_f(jk) ! Cumulative emissions for flaming combustion
            emis_info(i, jk + imax*2) = emis_s(jk) ! cumulative emissions for smoldering combustion
            emis_info(i, jk + imax*3) = emis_r(jk) ! cumulative emissions for residual combustion
         end do
      end do
      close(funit)
 30   format(4(F10.4, 1X), A8, 21(1X, F10.4), 124(1X, E15.7))
!
      ier = gdxyfll(phy_glb_gid, lat2x, lon2y, lat, lon, nhotspots)
      write(*,*) 'Total number of fire hotspots for CFFEPS is: ', nhotspots

      allocate(pe_info(nhotspots))
      pe_info = -1
      do i = 1, nhotspots
!        Round to the nearest i-j point
         int_x = nint(lat2x(i))
         int_y = nint(lon2y(i))
!
         do il = 1, numproc
            lx_lbound = rpi0_j0(1, il)                      ! i0
            lx_ubound = rpi0_j0(1, il) + rpi0_j0(3, il) - 1 ! in
            ly_lbound = rpi0_j0(2, il)                      ! j0
            ly_ubound = rpi0_j0(2, il) + rpi0_j0(4, il) - 1 ! jn

            if (int_x .ge. lx_lbound .and. int_x .le. lx_ubound  .and. &
                int_y .ge. ly_lbound .and. int_y .le. ly_ubound) then
               pe_info(i) = rpi0_j0(5, il) ! pe
               exit
            end if
         end do

         fire_info(i, i_gilc)  = real(int_x)  ! i-index on full physics grid
         fire_info(i, i_gjlc)  = real(int_y)  ! j-index on full physics grid
      end do
!
      deallocate(lat, lon, lat2x, lon2y)
   else ! my_pe /= 0
      allocate(pe_info(1))
      allocate(fire_info(1, 1))
      allocate(emis_info(1, 1))
   end if
!
   mjr_context=NULL_rpncomm_context
   ier = RPN_COMM_spread_context(mjr_context, RPN_COMM_GRID, 0, pe_info, nhotspots)
   if (ier /= 0) then
      call physeterror('chm_cffeps_init', 'problem in spread_context')
      return
   end if

!    From global array fire_info of dimension (nhotspots, fire_parameters) obtain
!    local POINTER lfire_Info for dimension (fire_parameters, lc_hotspots)
!    since some sources may not be part of the domain Sum(lc_hotspots) /= nhotspots

   lc_hotspots = RPN_COMM_spread(mjr_context, fire_info, nhotspots, &
                                 (fire_parameters + max_timesteps), lfire_info)
   if (lc_hotspots < 0) then
      call physeterror('chm_cffeps_init', 'problem in fire_info rpn_comm_spread')
      return
   end if

   lc_hotspots2 = RPN_COMM_spread(mjr_context, emis_info, nhotspots, &
                                 (4 * max_timesteps), lfire_emissions)
   if (lc_hotspots /= lc_hotspots2) then
      call physeterror('chm_cffeps_init', 'problem in emis_info rpn_comm_spread')
      return
   end if

   if (lc_hotspots > 0) then
      chm_do_cffeps_l = .true.
!
!    remap global i/j indices to local indices
      ier = ezgrid_find_ij0(phy_lcl_gid, phy_glb_gid, i0, j0, lni, lnj)
      do i = 1, lc_hotspots
         int_x = int(lfire_info(i_gilc, i)) - i0 + 1
         int_y = int(lfire_info(i_gjlc, i)) - j0 + 1
         ij = mapmod2phy(int_x, int_y, phy_lcl_ni, phydim_ni)
         lfire_info(i_gilc, i) = real(ij(1))
         lfire_info(i_gjlc, i) = real(ij(2))
      end do
!
! Build the list fire major point emissions speciated species
!# Speciate the emissions
      emis_nb_bins = max(2, (isize - 2)) ! no emissions for bins B and C
      nb_fire_me_species = ngas + nvocs + npm * emis_nb_bins
      allocate(fire_me_species(nb_fire_me_species))
      fire_me_species(1:ngas) = me_species_gas
      fire_me_species(ngas+1 : ngas+nvocs) = me_species_voc
      nsp = ngas + nvocs
      do sn = 1, npm
         do nbin = 1, emis_nb_bins
            write(zbin, '(Z1)') nbin
            nsp = nsp + 1
            fire_me_species(nsp) = trim(me_species_pm(sn)) // zbin
         end do
      end do

! (Indirect) address of the fire_me_species in species_master
      allocate(fire_me_species_index(nb_fire_me_species))
      do nsp = 1, nb_fire_me_species
         fire_me_species_index(nsp) = -1
         inner: do i = 1, nb_me_species
            il = me_species_index(i)
            if (trim(fire_me_species(nsp)) == trim(species_master(il) % me_name)) then
               fire_me_species_index(nsp) = il
               exit inner
            end if
         end do inner
      end do
!
   else
      chm_do_cffeps_l = .false.
      allocate(fire_me_species(1))
      allocate(fire_me_species_index(1))
   end if
   if (local_dbg) then
      write(*, *)' -- for processor ', my_pe, ' chm_do_cffeps_l is ', chm_do_cffeps_l
      write(*, *)' -- for processor ', my_pe, ' RPN_COMM_spread returned ', lc_hotspots , ' sources'
   end if

   deallocate(pe_info)
   deallocate(fire_info)
   deallocate(emis_info)
!
!***********************************************************************
   chm_cffeps_init = 1

   return
end function chm_cffeps_init
