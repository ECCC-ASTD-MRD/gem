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
! Fichier/File   : mach_adom2yb_main.ftn90
! Creation       : P. Makar, B. Pabla, S. Menard Feb 2007 for GEM-MACH
!                  Janusz Pudykiewicz for CHRONOS 1995
! Description    : Entry point for ADOM2 Gas-phase chemistry process describing NOX/VOC system
!
! Extra info     : - Modifications to do cloud photolysis correction (P. Makar, Sept 1998)
!                  - Unified units for all gas species to be ug/kg when in dynamic
!                    bus. For lumped radicals, estimated molecular weight is set
!                    to 100 g/mol per Craig's and Sylvie's suggestion. (Verica,
!                    Feb 2017)
!
! Arguments:  IN
!               p2d      -> pressure (pa) on thermodynamic levels
!               tplus    -> temperature (deg k)
!               csza     -> Cosine of solar zenith angle
!               rjval_lookup -> photolysis rates from look-up table
!               step     -> Flag for first chem. step in current
!
!             OUT
!               voc_diff -> Changes in VOC species' concentrations for SOA calc
!
!             IN/OUT
!               tppmgs   -> Dynamic gas species' concentration in ppmv
!               rjval    -> MESSy photolysis rates
!
!============================================================================
!
!!if_on
subroutine mach_adom2yb_main(p2d, tplus, hu_vmr, sigt, rjval, tppmgs, &
                             rjval_lookup, csza, voc_diff, step, gni, gnk)
   use mach_pkg_gas_mod,     only: nspec, nsp_soa_gases, njrxs
   use mach_adom2_rates_mod, only: ntype2, ntype4
!!if_off
   use chm_utils_mod,        only: chm_lun_out, global_debug, chm_error_l,  &
                                   CHM_MSG_DEBUG, chm_timestep
   use chm_consphychm_mod,   only: consth, rgasi, avno, pi
   use mach_pkg_gas_mod,     only: soa_gas_idx, njidx, nvar, o2, aircon
   use mach_pkg_adom2_mod,   only: nreac, nprcf, ind_c2h6, ind_ch4, ind_h2o, &
                                   ind_o2, ind_no, ind_m, ind_oh,            &
                                   layerboundary1, layerboundary2,           &
                                   layerconcentration1, layerconcentration2, &
                                   layerconcentration3
   use chm_nml_mod,          only: chm_soa_s, nk_start, chm_strato_s,      &
                                   chm_step_factor, chm_active_ch4_l,      &
                                   chm_timings_l, chm_messy_jval_l
   use mach_gas_headers_mod, only: mach_adom2_drive
   use mach_adom2_rates_mod, only: mach_adom2_uprate
   use linoz_param,          only: hu_linoz_tropo     !From RPNPHY
   implicit none
!
! Declare subroutine arguments
!
!!if_on
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: gni, gnk
   real(kind=4),    intent   (in) :: p2d     (gni, gnk)
   real(kind=4),    intent   (in) :: hu_vmr  (gni, gnk)
   real(kind=4),    intent   (in) :: sigt    (gni, gnk)
   real(kind=4),    intent   (in) :: tplus   (gni, gnk)
   real(kind=4),    intent   (in) :: csza    (gni)
   real(kind=8),    intent(inout) :: rjval   (gni * gnk, njrxs)
   real(kind=4),    intent(inout) :: tppmgs  (nspec, gni, gnk)
   real(kind=8),    intent  (out) :: voc_diff(nsp_soa_gases, gni, gnk)
   real(kind=8),    intent   (in) :: rjval_lookup(gni * gnk, ntype2+ntype4)
!!if_off
!
!  Declare local variables
!
   integer(kind=4)             :: i, k, npt, ii, isp, indx
!
   real(kind=8)                :: rga(gni * gnk, nreac)
   real(kind=8)                :: bg(gni * gnk, nprcf)
!
   real(kind=4)                :: conc(nspec), conc_oh
   real(kind=4)                :: bgs(nprcf), rgs(nreac)

   real(kind=4)                :: cno(gni * gnk), zen(gni * gnk)

!  index remapping 2-d arrays to 1-d
   integer(kind=4)             :: nn

   logical(kind=4)             :: local_dbg
!
!  Declare external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp
!
!  CODE BEGINS HERE
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_adom2yb_main [BEGIN]')
   if (chm_timings_L) call timing_start_omp(334, 'mach_adom2yb_main', 330)

!  number of gridpoints in one slice of data
   npt = gni * gnk

   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         nn = nn + 1
         zen(nn) = acos(csza(i)) * 180.0 / pi ! convert radians to degrees
         cno(nn) = max(1.0e-12, tppmgs(ind_no, i, k))
      end do
   end do

!     Calculate reaction rate constants for current slice:
   call mach_adom2_uprate(tplus, p2d, cno, rjval_lookup, rga, bg, npt, gni, gnk)

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
! Return the look-up table base j-values in to 'rjval' for possible output
! Convert back to the original units of  s-1 (from min-1) for output
         rjval(:, i) = rga(:, njidx(i)) / 60.0d0
      end do
   end if
!
! Some constant values
   conc(ind_o2) = o2     !  should relative amount of o2 change? upper atm?
   conc(ind_m)  = aircon

!  Timestep in minutes
!
   voc_diff = 0.0d0

   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         nn = nn + 1

         do isp = 1, nvar
            conc(isp) = tppmgs(isp, i, k)
         end do
! 'Fixed' species in ppmv
         conc(ind_ch4)  = dble(tppmgs(ind_ch4, i, k))
         conc(ind_c2h6) = layerconcentration1 ! constant concentration of ethane  as a function of height.
         if (sigt(i, k) > layerboundary1) conc(ind_c2h6) = layerconcentration2
         if (sigt(i, k) > layerboundary2) conc(ind_c2h6) = layerconcentration3
!
!  determine value of water vapour from the meteorology:
!    (a) converting humidity variable to ppmv, and
!    (b) taking the average of the met value to determine the value for the given time step.
!        for humidity in kg(h2o)/kg(air) to ppmv, the conversion (consth) is * 1.e+06/(m(h2o)/m(air));
!        28.96 used for molecular mass of air, 18.015 for molecular mass of h2o
         conc(ind_h2o) = max(consth * 1.0e-10, hu_vmr(i, k))

         if ((chm_strato_s == 'NIL') .or. (hu_vmr(i, k) >= hu_linoz_tropo)) then
            do ii = 1, nreac
               rgs(ii)   = real(rga(nn, ii))
            end do

            do ii = 1, nprcf
               bgs(ii)   = real(bg(nn, ii))
            end do
!
!  set up for gas-phase solver, called in scalar mode
            if (step - chm_step_factor /= 0) then
!  for "initialization" steps, use the high accuracy, small
!  internal time step option in the solver (indx = 1)
               indx = 0
            else
               indx = 1
            end if

!  call solver for current grid-point's gas-phase chemistry:
            if (local_dbg) then
               do isp = 1, nspec
                  if(conc(isp) <= 0.0) then
                     write (chm_lun_out, *) 'initial concentration less than or equal to zero for species ', &
                                            isp, conc(isp)
                  end if
               end do
               do ii = 1, nreac
                  if(rgs(ii) <= 0.0 .and. csza(i) > 1.0e-5) then
                     write (chm_lun_out, *) 'zero rate constant found for reaction ', ii, k
                  end if
               end do
            end if

            call mach_adom2_drive(conc, indx, zen(nn), i, k, rgs, bgs)
            if (indx < 0) then
               write(0,*) '### Error in mach_adom2yb_main ###'
               write(0,*) '# Return code from mach_adom2_drive: ', indx
               select case(indx)
                  case(-1)
                     write(0,*) '# Integration failed'
                  case(-10)
                     write(0,*) '# Integration took too many steps'
                  case(-20, -21, -22)
                     write(0,*) '# Integration aborted because it was about'
                     write(0,*) '# to make a division by zero'
               end select
               write(0,*) '###          ABORT         ###'
               chm_error_l = .true.
               return
            end if
!
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
                  voc_diff(ii, i, k) = max(dble(tppmgs(isp, i, k) - conc(isp)), &
                                           0.0d0)
               end do
            end if
!
!  gas species:  convert units back to ppmv
!  Apply ADOM tendencies below hu_linoz (10 ppmv specific humidity) used with
!  nk_start (36 or aproximately 100mb)
!
            do isp = 1, nvar
               tppmgs(isp, i, k) = conc(isp)
            end do

! Update CH4 concentration
            if (chm_active_ch4_l) then
               conc_oh = 0.5 * (tppmgs(ind_oh, i, k) + conc(ind_oh))
               conc(ind_ch4) = conc(ind_ch4) * exp(-1.0 * conc_oh &
                                 * real(rga(nn, 67)) / 60.0 * chm_timestep)
               tppmgs(ind_ch4, i, k) = conc(ind_ch4)
            end if
!
         end if
!
      end do
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_adom2yb_main [END]')
   if (chm_timings_L) call timing_stop_omp(334)
   !-----------------------------------------------------------------
   return
end subroutine mach_adom2yb_main
