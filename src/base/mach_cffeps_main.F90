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
! Fichier/File   : mach_cffeps_main.ftn90
! Creation       : A. Akingunola, J. Chen, P. Makar, and Kerry Anderson - 2018
!
! Description    : Online fire emission calculations based on CFFEPS.
! Arguments:
!            IN
!                slab_index    -> Slab number
!                metvar{2d,3d} -> local storage (chemistry lib only, see chm_exe) for met. fields
!                                 copied from Physics buses (see chm_load_metvar).
!
!             IN/OUT
!                chem_tr       -> Chemistry species' concentration.
!
! Extra info     :
!
!============================================================================
!
!!if_on
subroutine mach_cffeps_main(chem_tr, metvar2d, metvar3d, slab_index, istep)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
!!if_off
   use chm_metvar_mod,       only: MV2D_DLAT, MV2D_DLON, MV2D_DXDY, MV2D_PPLUS, &
                                   MV2D_TDIAG, MV2D_QDIAG, MV2D_RAINRATE,  &
                                   MV2D_WSDIAG, MV2D_WDDIAG, MV3D_TPLUS,   &
                                   MV3D_ZMOM, MV3D_ZPLUS, MV3D_SIGT, MV3D_RHO
   use chm_utils_mod,        only: chm_timestep
   use chm_datime_mod,       only: chm_dttim_s
   use chm_consphychm_mod,   only: pi
   use mach_cam_utils_mod,   only: isize
   use mach_cffeps_mod,      only: lc_hotspots, lfire_info, lfire_emissions,  &
                                   feps_type, met_type, emissions_type,       &
                                   met_levels, max_timesteps, fire_parameters,&
                                   fire_me_species_index, cffeps_steps, jdays,&
                                   i_gilc, i_gjlc, i_fuel, i_dtime, i_dj,     &
                                   i_bui, i_ffmc, i_dmc, i_dc, i_fmc, i_area, &
                                   i_cfl, i_sfl, i_ffl, i_wfl, i_cured, i_pc, &
                                   i_pdf, i_cbh, i_asp, i_gfl, i_gs, i_sd,    &
                                   i_sh, i_lfl, i_dfl, i_q0, i_emis, i_area0, &
                                   i_rsmk, i_zplm
   use mach_cffeps_speciations_mod, only: nb_fire_me_species,                &
                                          cffeps_emiss_speciations
   use mach_headers_mod,     only: mach_cffeps_calc
!   use chm_species_info_mod, only: species_master
   implicit none
!
!!if_on
   integer(kind=4), intent   (in) :: slab_index, istep
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
!!if_off
!
!  Local variables
!========================
   type(feps_type)      :: feps
   type(emissions_type) :: emissions
   type(met_type)       :: met

   integer(kind=4) :: i_prev, i_me, cur_source, fire_i, isp, iistep
   integer(kind=4) :: iyear, imonth, iday, ihour, iminute, cmcDj
   real(kind=4)    :: fcst_UTC, r_smoke, zplume
   real(kind=4), dimension(nb_fire_me_species) :: mj_emis
   real(kind=4), dimension(max_timesteps)      :: mws, ross


   integer(kind=4) :: kk, jk, zplm_indx, zplm_nlay
   real(kind=4), dimension(chm_nk+1) :: z_magl
   real(kind=4), dimension(chm_nk)   :: rho, weight, weight_f, rate_factor
   real(kind=4)                      :: cell_area   ! Area of a grid cell (m^2)
   real(kind=4)                      :: weight_f_sum, emiss_rate
   real(kind=4), parameter           :: GRAMMES_TO_MICROGRAMMES = 1.0e+06
   real(kind=4), parameter           :: conva = 180.0 / pi
   integer(kind=4), parameter        :: imax = 24

!
   read(chm_dttim_s(1:4),'(i4.4)') iyear      !YY
   read(chm_dttim_s(5:6),'(i2.2)') imonth     !MN
   read(chm_dttim_s(7:8),'(i2.2)') iday       !DD
   read(chm_dttim_s(10:11),'(i2.2)') ihour    !HH
   read(chm_dttim_s(12:13),'(i2.2)') iminute  !MM
   fcst_UTC = real(ihour) + real(iminute) / 60.0 !// UTC time of forecast hour
   !NOTE: The cmcDj here is slightly different from ijul_day from chm_datime_mod.
   cmcDj = jdays(imonth) + iday
   if (mod(iyear, 4) == 0 .and. imonth >2) cmcDj = cmcDj + 1
!
   iistep = istep + 24   ! include the already evaluated offline 'persistence' hours
! Loop over the total number of fire hotspots found: estimate emission rates
! using CFFEPS modules; calculates the weight factor for the entire column
! for each hotspot, and then apply it to the emission rate of all the species.
!
   i_prev = 0
   do cur_source = 1, lc_hotspots
      if (nint(lfire_info(i_gjlc, cur_source)) == slab_index) then
         fire_i = nint(lfire_info(i_gilc, cur_source))
      else
         cycle
      end if

      ! Only run CFFEPS at specified intervals
      if (mod(istep, cffeps_steps) == 0) then
         if (i_prev /= fire_i) then
            do kk = 1, met_levels
               jk = chm_nk - kk + 1 !
               met%P(kk) = metvar2d(fire_i, MV2D_PPLUS) * &
                           metvar3d(fire_i, jk, MV3D_SIGT) ! [Pa]
               met%T(kk) = metvar3d(fire_i, jk, MV3D_TPLUS)
               met%Z(kk) = metvar3d(fire_i, jk, MV3D_ZPLUS)
            end do

            met%ta = metvar2d(fire_i, MV2D_TDIAG)
            met%hus = metvar2d(fire_i, MV2D_QDIAG)
            met%ps = metvar2d(fire_i, MV2D_PPLUS) * 1.0e-2    ! surface pressure [mb]
            met%ws = metvar2d(fire_i, MV2D_WSDIAG) * 3.6      ! convert m/s to km/hr
            met%wd = metvar2d(fire_i, MV2D_WDDIAG) * pi / 180.0 !  Convert wind direction from degrees to radians
            met%pr = metvar2d(fire_i, MV2D_RAINRATE) * 3.6e6  ! Convert total precipitation from m/s to mm/hr
!
            i_prev = fire_i
         end if

         feps%lat = metvar2d(fire_i, MV2D_DLAT) * conva
         feps%lon = metvar2d(fire_i, MV2D_DLON) * conva - 360.0
         feps%fueltype = nint(lfire_info(i_fuel, cur_source))
         feps%dj = nint(lfire_info(i_dj, cur_source))
         feps%dtime = lfire_info(i_dtime, cur_source)
         feps%bui= lfire_info(i_bui, cur_source)
         feps%ffmc= lfire_info(i_ffmc, cur_source)
         feps%dmc = lfire_info(i_dmc, cur_source)
         feps%dc = lfire_info(i_dc, cur_source)
         feps%fmc = lfire_info(i_fmc, cur_source)
         feps%area = lfire_info(i_area, cur_source)
         feps%estarea = lfire_info(i_area0, cur_source)

         feps%cfl = lfire_info(i_cfl, cur_source)         ! Crown Fuel Load [kg/m^2]
         feps%sfl = lfire_info(i_sfl, cur_source)         ! Surface Fuel Load [kg/m^2]
         feps%ffl = lfire_info(i_ffl, cur_source)         ! Fine Fuel Load [kg/m^2]
         feps%wfl = lfire_info(i_wfl, cur_source)         ! Woody Fuel Load [kg/m^2]

         feps%curing = lfire_info(i_cured, cur_source)    ! Percent Cured for O1a/O1b
         feps%pc  = lfire_info(i_pc, cur_source)          ! Percent Confier for M1/M2
         feps%pdf = lfire_info(i_pdf, cur_source)         ! Percent Dead Fir for M3/M4
         feps%cbh = lfire_info(i_cbh, cur_source)         ! Crown to Base Height [m]
         feps%aspect = lfire_info(i_asp, cur_source)      ! Aspect [degrees]
         feps%gfl = lfire_info(i_gfl, cur_source)         ! Grass Fuel Load [kg/m^2]
         feps%gs  = lfire_info(i_gs, cur_source)          ! Slope [percent]
         feps%sd  = lfire_info(i_sd, cur_source)          ! C6 Stand Density [stems/ha]
         feps%sh  = lfire_info(i_sh, cur_source)          ! C6 Stand Height [m]
         feps%lfl = lfire_info(i_lfl, cur_source)         ! Litter Fuel Load [kg/m^2]
         feps%dfl = lfire_info(i_dfl, cur_source)         ! Duff (Fermentation) Fuel Load [kg/m^2]

         feps%Qo = lfire_info(i_q0, cur_source)
         feps%totalemissions = lfire_info(i_emis, cur_source)
      !TODO: max_timesteps currently fixed at imax (24); must fix
         do jk = 1, max_timesteps
            feps%Qs(jk) = lfire_info(jk + fire_parameters, cur_source)

            mws(jk) = lfire_emissions(jk, cur_source)
            emissions%f(jk) = lfire_emissions(imax + jk, cur_source)
            emissions%s(jk) = lfire_emissions(2*imax + jk, cur_source)
            emissions%r(jk) = lfire_emissions(3*imax + jk, cur_source)

            ross(jk) = 0.0 ! Not properly connected yet
         end do

         call mach_cffeps_calc(feps, emissions, r_smoke, zplume, met, &
                               mws, ross, cmcDj, fcst_UTC, ihour, iistep)
!
! Return updated values to the local hotspot data pointer
         lfire_info(i_area, cur_source) = feps%area
         ! amount of energy previously injected into atmosphere (J)
         lfire_info(i_q0, cur_source)   = feps%Qo
         ! total emissions from the fire (tonnes)
         lfire_info(i_emis, cur_source) = feps%totalemissions
         lfire_info(i_rsmk, cur_source) = r_smoke
         lfire_info(i_zplm, cur_source) = zplume
         do jk = 1, max_timesteps
            lfire_info(jk + fire_parameters, cur_source) = feps%Qs(jk)

            lfire_emissions(jk, cur_source) = mws(jk)
            lfire_emissions(imax + jk, cur_source) = emissions%f(jk)
            lfire_emissions(2*imax + jk, cur_source) = emissions%s(jk)
            lfire_emissions(3*imax + jk, cur_source) = emissions%r(jk)
         end do
!
      else
         do jk = 1, max_timesteps
            emissions%f(jk) = lfire_emissions(imax + jk, cur_source)
            emissions%s(jk) = lfire_emissions(2*imax + jk, cur_source)
            emissions%r(jk) = lfire_emissions(3*imax + jk, cur_source)
         end do
         r_smoke = lfire_info(i_rsmk, cur_source)
         zplume = lfire_info(i_zplm, cur_source)
      end if
!
      if (istep == 0) cycle ! Chemistry is not done at this time

      cell_area = metvar2d(fire_i, MV2D_DXDY)
      ! Compute height of the model level interface in meter above ground-level
      do kk = 1, chm_nk
         z_magl(kk) = metvar3d(fire_i, kk, MV3D_ZMOM)
         rho(kk)    = metvar3d(fire_i, kk, MV3D_RHO)
      end do
!
!  Weight of each cell in the column above the fire hotspot
      weight = 0.
      do kk = chm_nk, 1, -1
       ! Evaluate position zplm from the ground up
         zplm_indx = kk
         if (z_magl(kk) >= zplume) then
            exit
         end if
      end do
!
      zplm_nlay = chm_nk - zplm_indx + 1  ! find # of layers below ZPLM
    ! calculate smoke density below zplm
      weight_f = 0.0
      weight_f(zplm_indx:chm_nk) = rho(zplm_indx:chm_nk) * r_smoke
      weight_f_sum = sum(weight_f)
      do kk = zplm_indx, chm_nk
         if (weight_f_sum > 0.0) then
            weight(kk) = weight_f(kk) / weight_f_sum
         else
            weight(kk) = 1.0 / zplm_nlay
         end if
      end do
! normalize weighing factor by layer thickness
      do kk = 1, chm_nk-1
         weight(kk) = weight(kk) / (z_magl(kk) - z_magl(kk + 1))
      end do
      weight(chm_nk) = weight(chm_nk) / z_magl(chm_nk)

!     Estimate the mass of model species emiited from the fire
!*NOTE: We pass 'isize' as argument from here instead of using the mach_cam_utils
!       module value, so that the 'mach_cffeps_emiss_speciations' can be similarly
!       called (un-modified) in the CFFEPS stand-alone instance.
      call cffeps_emiss_speciations(emissions, mj_emis, isize)
!
!  Apply plumerise injection for species to the chemistry dynamical bus
      do kk = 1, chm_nk
         rate_factor(kk) = GRAMMES_TO_MICROGRAMMES * weight(kk) * chm_timestep / &
                           (cell_area * rho(kk))
      end do
      do i_me = 1, nb_fire_me_species
         emiss_rate = mj_emis(i_me)
         isp = fire_me_species_index(i_me)

         if (emiss_rate > 0.0 .and. isp > 0) then

!            if (local_dbg) then
!               write(*, *) 'FOUND SPECIES {DYN,mjpt} NAME : ', &
!                            species_master(isp) % dyn_name,    &
!                            species_master(isp) % me_name
!               write(*, *) 'current source value  : ', emiss_rate
!            end if

            do kk = 1, chm_nk
               chem_tr(fire_i, kk, isp) = chem_tr(fire_i, kk, isp) + &
                                          emiss_rate * rate_factor(kk)
            end do
         end if
      end do
!
   end do

   return
end subroutine mach_cffeps_main
