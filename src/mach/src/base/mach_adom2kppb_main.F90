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
! Fichier/File   : mach_adom2kppb_main.ftn90
! Creation       : Kenjiro Toyota, Jack Chen, Deji Akingunola, 2018
! Description    : This KKP-based ADOM2 solver is an alternative version to
!                  mach_adom2kpp, using the same reaction rates as mach_adom2yb
!
! Extra info     :
!               - Code structure underneath this subroutine based largely on so-called Mathew Russell version
!               - Daytime and nighttime solvers are different by cutting off all the photolysis
!                 reactions from the latter
!               - The solver can run with chm_messy_jval_l = .true. or .false.
!
! Arguments:  IN
!               p2d      -> pressure (pa) on thermodynamic levels
!               tplus    -> temperature (deg k)
!               hu_vmr   -> Atmospheric humidity (in ppmv)
!               sigt     -> Thermodynamic levels sigma coordinates
!               rjval_lookup -> photolysis rates from look-up table
!
!             OUT
!               voc_diff -> Changes in VOC species' concentrations for SOA calc
!
!             IN/OUT
!               tppmgs   -> Dynamic gas species' concentration in ppmv
!               rjval    -> MESSy photolysis rates
!               hstart   -> Solver internal timestep
!
!============================================================================
!
!!if_on
subroutine mach_adom2kppb_main(p2d, tplus, hu_vmr, sigt, rjval, tppmgs, hstart, &
                               rjval_lookup, voc_diff, gni, gnk)
   use mach_pkg_gas_mod,     only: nspec, nsp_soa_gases, njrxs
   use mach_adom2_rates_mod, only: ntype2, ntype4
!!if_off
   use chm_utils_mod,        only: chm_lun_out, global_debug, CHM_MSG_DEBUG,&
                                   chm_timestep, chm_error_l
   use chm_consphychm_mod,   only: consth, rgasi, avno
   use mach_pkg_gas_mod,     only: soa_gas_idx, njidx, minRx, minconc,      &
                                   o2ppmv, airppmv
   use mach_pkg_adom2_mod,   only: nreac, nprcf, nvar, ind_c2h6, ind_ch4,   &
                                   ind_h2o, ind_o2, ind_no, ind_m,          &
                                   ind_oh, layerboundary1, layerboundary2,  &
                                   layerconcentration1, layerconcentration2, &
                                   layerconcentration3
   use chm_nml_mod,          only: chm_soa_s, nk_start, chm_strato_s,       &
                                   chm_messy_jval_l, chm_timings_l,         &
                                   chm_active_ch4_l
   use mach_gas_headers_mod, only: mach_kppb_adom2_interface
   use mach_adom2_rates_mod, only: mach_adom2_uprate
   implicit none
!
! Declare subroutine arguments
!
!!if_on
   integer(kind=4), intent   (in) :: gni, gnk
   real(kind=4),    intent   (in) :: p2d     (gni, gnk)
   real(kind=4),    intent   (in) :: hu_vmr  (gni, gnk)
   real(kind=4),    intent   (in) :: sigt    (gni, gnk)
   real(kind=4),    intent   (in) :: tplus   (gni, gnk)
   real(kind=8),    intent(inout) :: rjval   (gni * gnk, njrxs)
   real(kind=4),    intent(inout) :: tppmgs  (nspec, gni, gnk)
   real(kind=8),    intent(inout) :: hstart  (gni, gnk)
   real(kind=8),    intent  (out) :: voc_diff(nsp_soa_gases, gni, gnk)
   real(kind=8),    intent   (in) :: rjval_lookup(gni * gnk, ntype2+ntype4)
!!if_off
!
!  Declare local variables
!
   integer(kind=4)             :: i, k, ii, isp
!
   real(kind=8)                :: rga(gni * gnk, nreac)
   real(kind=8)                :: bg(gni * gnk, nprcf)
!
   real(kind=8)                :: conc(nspec), conc1d(nspec)
   real(kind=8)                :: bgs(nprcf), rgs(nreac)

   real(kind=4)                :: cno(gni * gnk)
!  Parameters for ADOM2-KPPB RODAS3 solver
   real(kind=8)                :: tstep_kpp         ! Solver internal timestep in minutes
   integer(kind=4)             :: day_night_flag    ! flag for switching between day(1)/night(0) solvers

!  index remapping 2-d arrays to 1-d
   integer(kind=4)             :: nn
!
   real(kind=8)                :: atol  ! Absolute tolerance parameter for solver (default)
   real(kind=8)                :: rtol  ! Relative tolerance parameter for solver (default)
! OH average concentration
   real(kind=8)                :: conc_oh
!
   logical(kind=4)             :: local_dbg, kpp_solve
!
!  Declare external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp
!
!  CODE BEGINS HERE
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_adom2kppb_main [BEGIN]')
   if (chm_timings_L) call timing_start_omp(334, 'mach_adom2kppb_main', 330)

   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         nn = nn + 1
         cno(nn) = max(0.0, tppmgs(ind_no, i, k))
      end do
   end do
!
!     Calculate reaction rate constants for current slice:
   call mach_adom2_uprate(tplus, p2d, cno, rjval_lookup, rga, bg, &
                          gni * gnk, gni, gnk)
!
   if (chm_messy_jval_l) then
! update Rx rate array rga(:,:) with MESSY rates:
      do i = 1, njrxs
         rga(:, njidx(i)) = rjval(:, i)
      end do
   else
! if not MESSY, rjval is used to hold the light attenuation coefficients (for instance,
! under the canopy shading effect)
      do i = 1, njrxs
         rga(:, njidx(i)) = rga(:, njidx(i)) * rjval(:, i)
!
! Return the look-up table base j-values in to 'rjval' for possible output;
! Convert back to the original units of  s-1 (from min-1) for output
         rjval(:, i) = rga(:, njidx(i)) / 60.0d0
      end do
   end if

   voc_diff = 0.0d0

! Some constant values
   conc(ind_o2)  = o2ppmv     !  should relative amount of o2 change? upper atm?
   conc(ind_m)   = airppmv

   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         nn = nn + 1

         do isp = 1, nvar
            conc(isp) = dble(tppmgs(isp, i, k))
         end do
! 'Fixed' species in ppmv
         conc(ind_ch4)  = dble(tppmgs(ind_ch4, i, k))
         conc(ind_c2h6) = dble(layerconcentration1) ! constant concentration of ethane  as a function of height.
         if (sigt(i, k) > layerboundary1) conc(ind_c2h6) = dble(layerconcentration2)
         if (sigt(i, k) > layerboundary2) conc(ind_c2h6) = dble(layerconcentration3)
!
!  determine value of water vapour from the meteorology:
!    (a) converting humidity variable to ppmv, and
!    (b) taking the average of the met value to determine the value for the given time step.
!        for humidity in kg(h2o)/kg(air) to ppmv, the conversion (consth) is * 1.e+06/(m(h2o)/m(air));
!        28.96 used for molecular mass of air, 18.015 for molecular mass of h2o
         conc(ind_h2o) = dble(max(consth * 1.0e-10, hu_vmr(i, k)))

         do ii = 1, nreac
            rgs(ii)   = rga(nn, ii)
         end do

         do ii = 1, nprcf
            bgs(ii)   = bg(nn, ii)
         end do
!
!  Day(1) vs. night(0) flag for solver
!
         day_night_flag = 1
!** The night time flad (day_night_flag = 0) actually slows down the solver,
!** and sometimes can lead to the solver outrightly failing (Mantis #4571)
!         if (chm_messy_jval_l) then
!            ! nighttime for messy-jval, JNO2 = 0
!            if (rgs(njidx(1)) < minRx) day_night_flag = 0
!         else
!            ! nighttime for Table-J, cos(sza) = 0.00001
!            if (csza(i) <= 0.00001) day_night_flag = 0
!         end if
!
!
         atol = 1.d-20  ! Absolute tolerance parameter for solver (default)
         rtol = 5.d-2   ! Relative tolerance parameter for solver (default)
!         rtol = 1.d-1   ! Relative tolerance parameter for solver (somewhat aggressive)

         kpp_solve = .true.

! Iterate over the KPP solver with varying rtol until 'steady' solutions are achieved
         do while (kpp_solve .and. (rtol >= 5.1d-5))
            conc1d = conc
            tstep_kpp = hstart(i, k) / 60.0d0  ! Internal solver timestep in minutes
            call mach_kppb_adom2_interface(conc1d, rgs, bgs, atol, rtol, &
                                           tstep_kpp,  day_night_flag)
            if (chm_error_l) return
! Check the returned concentrations from the solver for negatives;
! If found, lower rtol for the grid and repeat solver until all conc are positive
            if (any(conc1d < 0.0d0)) then
               if (local_dbg) then
                  write(*, *) 'negative concentrations in the KPP solution'
                  write(*, *) ' rtol: ', rtol
                  write(*, *) 'i, k, conc1d: ', i, k, conc1d
               end if

               ! Reduce the solver's relative tolerance for this grid
               rtol = rtol * 0.1d0

               if (minval(conc1d) > -1.0d-3) then
                  conc1d = max(minconc, conc1d)
                  kpp_solve = .false.
               else if  (rtol < 5.1d-5) then
                  write(*,*) 'negative concentration from kppb'
                  do isp = 1, nvar
                     write(*,*) isp, conc1d(isp)
                  end do
                  chm_error_l = .true.
                  return
               end if
!
            else
              kpp_solve = .false.
            end if
         end do
! Internal solver timestep
         hstart(i, k) = tstep_kpp * 60.0d0

         if (chm_soa_s == 'JIANG') then
     ! The concentrations of gas species involved in the SOA yield calculations
     ! are expected to be in the unit of ppmv
     ! voc_diff(1,:,:) -> change in ALKA mixing ratio in ppmv
     ! voc_diff(2,:,:) -> change in AROM mixing ratio in ppmv
     ! voc_diff(3,:,:) -> change in TOLU mixing ratio in ppmv
     ! voc_diff(4,:,:) -> change in ALKE mixing ratio in ppmv
     ! voc_diff(5,:,:) -> change in ISOP mixing ratio in ppmv
            do ii = 1, nsp_soa_gases
               isp = soa_gas_idx(ii)
               voc_diff(ii, i, k) = max(dble(tppmgs(isp, i, k)) - conc1d(isp), &
                                        0.0d0)
            end do
         end if
!
!  gas species:  convert units back to ppmv
         do isp = 1, nvar
            tppmgs(isp, i, k) = real(conc1d(isp))
         end do
!
! Update CH4 concentration
         if (chm_active_ch4_l) then
            conc_oh = 0.5d0 * (conc1d(ind_oh) + conc(ind_oh))
            conc1d(ind_ch4) = conc1d(ind_ch4) * exp(-1.0d0 * conc_oh * &
                                                rga(nn, 67) / 60.0d0 * &
                                                dble(chm_timestep))
            tppmgs(ind_ch4, i, k) = real(conc1d(ind_ch4))
         end if
!
      end do
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_adom2kppb_main [END]')
   if (chm_timings_L) call timing_stop_omp(334)
   !-----------------------------------------------------------------
   return
end subroutine mach_adom2kppb_main
