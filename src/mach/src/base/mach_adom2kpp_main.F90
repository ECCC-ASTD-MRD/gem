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
! Fichier/File   : mach_adom2kpp_main.ftn90
! Creation       : Jack Chen Dec 2017
! Description    : Entry point for KPP version for ADOM2 mechanism
!
! Extra info     :
!               - this version has nreact=112 instead of 114 as Rx 36 and 42 are merged to Rx 35 and 41, respetively.
!               - this version has nfix=6 instead of 5 due to the need of a 'dummy' reactoin product
!               - this version requires MESSy to produce Jx  (note different Rx index)
!               - gas-reaction coneff is now calculated from KPP code instead of "mach_gas_uprate_adom2.ftn90"
!               - conc unit is now mole/cm3 while rate unit is molec/cm3.sec instead of ppmv and ppmv/min in Y&B solver.
!               - update to variable index ordering to be similar to those used in KPP code
!                  - Modification to calculate NON-mass conserving condensable
!                    organic mass for later SOA calculations. Method of Odum
!                    et al., and the two product yield method from Seinfeld's
!                    lab, along with corrections to account for lumping, have
!                    been used. (P. Makar, July 1999).
!
! Arguments:  IN
!               p2d      -> pressure (pa) on thermodynamic levels
!               tplus    -> Air temperature (K)
!               hu_vmr   -> Atmospheric humidity (in ppmv)
!               sigt     -> Thermodynamic levels sigma coordinates
!               rjval    -> MeSSy photolysis rates
!
!             OUT
!               voc_diff -> Changes in VOC species' concentrations for SOA calc
!
!             IN/OUT
!               tppmgs   -> Dynamic gas species' concentration in ppmv
!               hstart   -> Solver internal timestep
!
!============================================================================
!
!!if_on
subroutine mach_adom2kpp_main(p2d, tplus, hu_vmr, sigt, rjval, tppmgs, hstart, &
                              voc_diff, gni, gnk)

   use mach_pkg_gas_mod,     only: nspec, nsp_soa_gases, njrxs
!!if_off
   use chm_utils_mod,        only: chm_lun_out, global_debug, chm_error_l, &
                                   CHM_MSG_DEBUG, chm_timestep
   use chm_consphychm_mod,   only: consth, rgasi, avno
   use chm_nml_mod,          only: chm_soa_s, nk_start, chm_strato_s,        &
                                   chm_timings_l, chm_active_ch4_l
   use mach_pkg_gas_mod,     only: soa_gas_idx, njidx, minRx, minconc,       &
                                   o2, aircon, nreact, nvar
   use mach_pkg_adom2_mod,   only: ind_c2h6, ind_ch4, ind_h2o, ind_o2,       &
                                   ind_no, ind_m, ind_dum, ind_oh,           &
                                   layerboundary1, layerboundary2,           &
                                   layerconcentration1, layerconcentration2, &
                                   layerconcentration3
   use mach_gas_headers_mod, only: mach_kpp_adom2_rates, mach_kpp_interface
   implicit none
!
! Declare subroutine arguments
!
!!if_on
   integer(kind=4), intent   (in) :: gni, gnk
   real(kind=4),    intent   (in) :: p2d     (gni, gnk)
   real(kind=4),    intent   (in) :: tplus   (gni, gnk)
   real(kind=4),    intent   (in) :: hu_vmr  (gni, gnk)
   real(kind=4),    intent   (in) :: sigt    (gni, gnk)
   real(kind=8),    intent   (in) :: rjval   (gni * gnk, njrxs)
   real(kind=4),    intent(inout) :: tppmgs  (nspec, gni, gnk)
   real(kind=8),    intent(inout) :: hstart  (gni, gnk)
   real(kind=8),    intent  (out) :: voc_diff(nsp_soa_gases, gni, gnk)
!!if_off
!
!  Declare local variables
!
   integer(kind=4) :: i, k, jj, isp
!  index remapping 2-d arrays to 1-d
   integer(kind=4) :: nn
!
!  no concentration array.  used in hcho + ho2 reaction rate (no. 53)
   real(kind=4)    :: cno(gni, gnk)
! conversion from ppm to molecules/cm3
   real(kind=4)    :: conx1(gni, gnk)
!
   real(kind=8)    :: rxrt2d(gni * gnk, nreact) ! reaction rate
   real(kind=8)    :: conc(nspec)     ! 1d conc. for all species
   real(kind=8)    :: rxrt(nreact)  ! 1d reaction rate for solver
!
! OH average concentration
   real(kind=8)    :: conc_oh
!
! ATOL - Absolute tolerance
   real(kind=8)    :: atol(nvar)
! RTOL - Relative tolerance
   real(kind=8)    :: rtol(nvar)
! Scalar/Vector tolerance control
   integer(kind=4) :: itolctr

   logical(kind=4) :: local_dbg
!
!  Declare external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp
!
!  CODE BEGINS HERE
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_adom2kpp_main [BEGIN]')
   if (chm_timings_L) call timing_start_omp(334, 'mach_adom2kpp_main', 330)

   do k = 1, gnk
      do i = 1, gni
! rgasi is gas constant, p2d is pressure in Pa ;  TPLUS in Kelvin
         conx1(i, k) = avno * p2d(i, k) * 1.E-12 / (rgasi  * tplus(i, k))
         cno(i, k)   = max(0.0, tppmgs(ind_no, i, k) * conx1(i, k))
      end do
   end do

!  Calculate reaction rate constants for current slice:
   call mach_kpp_adom2_rates(rxrt2d, tplus, p2d, conx1, cno, gni, gnk)
!
! update Rx rate array rxrt2d(:,:) with MESSY rates:
   do i = 1, njrxs
      rxrt2d(:, njidx(i)) = rjval(:, i) / 60.d0  ! convert to sec-1 for KPP rate
   end do
!
   voc_diff = 0.0d0
   conc(ind_dum) = minconc

!~~~>
   atol = 1.d-20 ! absolute tolerance
   rtol = 5.0d-2 ! relative tolerance
!   For Scalar tolerances (itolctr /= 0) the KPP solver code uses atol(1) and rtol(1)
!   For Vector tolerances (itolctr == 0) the code uses atol(1:nvar) and rtol(1:nvar)
   itolctr = 1

! Variable species convert ppmv to molec/cm3
   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         nn = nn + 1
         do isp = 1, nvar
            conc(isp) = dble(tppmgs(isp, i, k) * conx1(i, k))
         end do
! 'Fixed' species convert ppmv to molec/cm3
         conc(ind_ch4) = dble(tppmgs(ind_ch4, i, k) * conx1(i, k))
         if (chm_active_ch4_l) conc_oh = conc(ind_oh)
         conc(ind_o2)  = dble(o2  * conx1(i, k))  ! should relative amount of o2 change? upper atm?
         conc(ind_m)   = dble(aircon * conx1(i, k))

         conc(ind_c2h6) = dble(layerconcentration1 * conx1(i, k)) ! constant concentration of ethane  as a function of height.
         if (sigt(i, k) > layerboundary1) conc(ind_c2h6) = dble(layerconcentration2*conx1(i, k))
         if (sigt(i, k) > layerboundary2) conc(ind_c2h6) = dble(layerconcentration3*conx1(i, k))

!  determine value of water vapour from the meteorology:
!    (a) converting humidity variable to ppmv, and
!    (b) taking the average of the met value to determine the value for the given time step.
!        for humidity in kg(h2o)/kg(air) to ppmv, the conversion (consth) is * 1.e+06/(m(h2o)/m(air));
!        28.96 used for molecular mass of air, 18.015 for molecular mass of h2o
         conc(ind_h2o) = dble(max(consth * 1.0e-10, hu_vmr(i, k)) * conx1(i, k))

! update rx rate for this grid
         rxrt(1:nreact) = max(minRx, rxrt2d(nn, 1:nreact))

! call solver and return new conc1d and hstart
         call mach_kpp_interface(conc, rxrt, atol, rtol, hstart(i, k), itolctr)
         if (chm_error_l) return

         if (chm_soa_s == 'JIANG') then
     ! The concentrations of gas species involved in the SOA yield calculations
     ! are expected to be in the unit of ppmv
     ! voc_diff(1,:,:) -> change in ALKA mixing ratio in ppmv
     ! voc_diff(2,:,:) -> change in AROM mixing ratio in ppmv
     ! voc_diff(3,:,:) -> change in TOLU mixing ratio in ppmv
     ! voc_diff(4,:,:) -> change in ALKE mixing ratio in ppmv
     ! voc_diff(5,:,:) -> change in ISOP mixing ratio in ppmv
            do jj = 1, nsp_soa_gases
               isp = soa_gas_idx(jj)
               voc_diff(jj, i, k) = max(dble(tppmgs(isp, i, k)) - &
                                    (conc(isp) / dble(conx1(i, k))), 0.0d0)
            end do
         end if
!
!  gas species:  convert units back to ppmv
         do isp = 1, nvar
            tppmgs(isp, i, k) = real(conc(isp)) / conx1(i, k)
         end do
!
! Update CH4 concentration
         if (chm_active_ch4_l) then
            conc_oh = 0.5d0 * (conc_oh + conc(ind_oh))
            conc(ind_ch4) = conc(ind_ch4) * exp(-1.0d0 * conc_oh * &
                                          rxrt2d(nn, 65) * dble(chm_timestep))
            tppmgs(ind_ch4, i, k) = real(conc(ind_ch4)) / conx1(i, k)
         end if
!
      end do
   end do
!
   call msg_toall(CHM_MSG_DEBUG, 'mach_adom2kpp_main [END]')
   if (chm_timings_L) call timing_stop_omp(334)
   !-----------------------------------------------------------------
   return
end subroutine mach_adom2kpp_main
