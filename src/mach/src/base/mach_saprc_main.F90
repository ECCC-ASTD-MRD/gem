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
! Fichier/File   : mach_saprc_main.ftn90
! Creation       : Jack Chen modified from mach_adom2_main.ftn90 by '
!                  P. Makar, B. Pabla, S. Menard Feb 2007 for GEM-MACH
!
! Description    : Entry point for Gas-phase chemistry process
!                  Prepares species concentration and reaction rates and call
!                   chem. solver (KPP interface routine) and update results to
!                   bus. Also calls Linoz/SOA processes if invoked.
!
! Extra info     : Photolysis Rx Rate comes from the MESSy module
!                  Other gas-phase Rx Rate Coeff. from KPP generated routine
!                  need to be careful of the variable order specified in module
!
! Arguments:  IN
!               p2d    -> pressure field
!               tplus  -> temperature field
!               hu_ppm -> water vapor field (units ppmv)
!               rjval  -> photolysis rates
!                moi   -> total condensible organics (ug/m3)
!
!             OUT
!                dsoa  -> dSOA formed from input dVOC (ug/m3)
!
!             INOUT
!               tppmgs -> concentrations (ppmv)
!               hstart ->
!
!============================================================================
!
!!if_on
subroutine mach_saprc_main(p2d, tplus, hu_ppm, rjval, tppmgs, hstart, &
                           moi, dsoa, gni, gnk)
   use mach_pkg_gas_mod,     only: nspec, njrxs
!!if_off
   use chm_utils_mod,        only: chm_error_l,  CHM_MSG_DEBUG
   use chm_consphychm_mod,   only: consth, rgasi, avno
   use mach_pkg_gas_mod,     only: njidx, airppmv, n2ppmv, co2ppmv, o2ppmv, &
                                   ch4ppmv, minconc, minRX, nreact, nvar,   &
                                   soa_gas_idx, nsp_soa_gases, nvsoa
   use chm_nml_mod,          only: nk_start, chm_timings_l, chm_pkg_gas_s,  &
                                   chm_soa_s
   use mach_pkg_saprc_mod,   only: ind_M, ind_N2, ind_O2, ind_CH4, ind_CO2, &
                                   ind_H2O
   use mach_gas_headers_mod, only: mach_kpp_saprc07_rates, mach_kpp_interface, &
                                   mach_soa_odum
   implicit none
!
! Declare subroutine arguments
!
!!if_on
   integer(kind=4), intent(in) :: gni, gnk
   real(kind=4), intent   (in) :: p2d   (gni, gnk)
   real(kind=4), intent   (in) :: tplus (gni, gnk)
   real(kind=4), intent   (in) :: hu_ppm(gni, gnk)
   real(kind=8), intent   (in) :: rjval (gni * gnk, njrxs)
   real(kind=8), intent   (in) :: moi   (gni, gnk)
   real(kind=4), intent(inout) :: tppmgs(nspec, gni, gnk)
   real(kind=8), intent(inout) :: hstart(gni, gnk)
   real(kind=4), intent  (out) :: dsoa  (gni, gnk)
!!if_off
!
!  Declare local variables
!
   integer(kind=4) :: i, k, nn, ns, jj, isp

!  gas-phase species in ppmv units (nvar=transported,nfix=fixed,nspec=nvar+nfix)
   real(kind=8)    :: rxrt2d(gni * gnk, nreact) ! reaction rate
   real(kind=8)    :: conc1d(nspec)  ! 1d conc. for all species for solver
   real(kind=8)    :: rxrt1d(nreact) ! 1d reaction rate for solver
   real(kind=8)    :: atol(nvar), rtol(nvar) ! Absolute & relative tolerances
   real(kind=8)    :: conx1(gni * gnk) ! conversion from ppm to molecules/cm3
   real(kind=4)    :: minh2o
! Scalar/Vector tolerance control
   integer(kind=4) :: itolctr
!
!  SOA related fields
   real(kind=8) :: voc_diff (nsp_soa_gases, gni, gnk)
   real(kind=8) :: ksoa     (3, nvsoa, gni, gnk)
   real(kind=8) :: kHO2_r2o2(gni, gnk), kNO_r2o2 (gni, gnk)
!
!  Declare external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

!  CODE BEGINS HERE
!
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_saprc_main [BEGIN]')
   if (chm_timings_L) call timing_start_omp(334, 'mach_saprc_main', 330)

!! unit conversion factor
!  conx1 is for conversion from ppm to molecules/cm3, p2d in pa ; rho in kg/m3 ; TPLUS in Kelvin
   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
      do i = 1, gni
         nn = nn + 1
         conx1(nn) = dble(avno * p2d(i,k) / (rgasi  * tplus(i, k))) * 1.D-12
      end do
   end do
!
!---------------------
!  2D reaction rate array
!---------------------
   call mach_kpp_saprc07_rates(rxrt2d, ksoa, kNO_r2o2, kHO2_r2o2, tplus, &
                               conx1, gni, gnk)
!
!  complement rxrt2d rate array using previously calculated photolysis rates
   do i = 1, njrxs
      rxrt2d(:, njidx(i)) = rjval(:, i)
   end do

   atol(1:nvar) = 1.0d-3 ! absolute tolerance
!   atol(1:nvar) = 1.0d-5
   rtol(1:nvar) = 5.0d-2 ! relative tolerance
!   For Scalar tolerances (itolctr /= 0)  the KPP solver code uses atol(1) and rtol(1)
!   For Vector tolerances (itolctr == 0) the code uses atol(1:nvar) and rtol(1:nvar)
   itolctr = 1

   voc_diff = 0.0d0

   nn = (nk_start - 1) * gni
   do k = nk_start, gnk
!---------------------
!  concentration array
!---------------------
      conc1d = minconc ! Initialize

      do i = 1, gni
         nn = nn + 1

         do ns = 1, nvar
            conc1d(ns) = dble(tppmgs(ns, i, k)) * conx1(nn) !ppmv to molec/cm3
         end do

! Specify fixed species convert ppmv to molec/cm3
!! severial fixed species for now, ideally species such as CO2, N, H2 should come from climatology
         conc1d( ind_M   ) = airppmv * conx1(nn)
         if (chm_pkg_gas_s(1:9) == 'SAPRC07CS') then
            conc1d( ind_N2  ) = n2ppmv * conx1(nn)
         end if
         conc1d( ind_O2  ) = o2ppmv  * conx1(nn)
         conc1d( ind_CH4 ) = ch4ppmv * conx1(nn)
         conc1d( ind_CO2 ) = co2ppmv * conx1(nn)

!  determine value of water vapour from the meteorology:
!    (a) converting humidity variable to molec/cm3
!    (b) taking the average of the met value to determine the value for the given time step.  for humidity in
         minh2o = max(consth * 1.0e-10, hu_ppm(i,k))
         conc1d( ind_H2O ) = dble(minh2o) * conx1(nn)

! update rx rate for this grid
         rxrt1d(1:nreact) = max(MinRX, rxrt2d(nn,1:nreact))

         call mach_kpp_interface(conc1d, rxrt1d, atol, rtol, hstart(i, k), itolctr)
         if (chm_error_l) return

!! get dVOC for SOA in ppmv
         if (chm_soa_s /= 'NIL') then
            do jj = 1, nsp_soa_gases
               isp = soa_gas_idx(jj)
               voc_diff(jj, i, k) = max(dble(tppmgs(isp, i, k)) - &
                                    (conc1d(isp) / conx1(nn)), 0.0d0)
            end do
         end if

!! convert molec/cm3 to ppmv
         do ns = 1, nvar
            tppmgs(ns, i, k) = real(max(MinConc, conc1d(ns) / conx1(nn)))
         end do
      end do
   end do

   if (chm_soa_s(1:4) == 'ODUM') then
      call mach_soa_odum(moi, dsoa, voc_diff, ksoa, kNO_r2o2, kHO2_r2o2, &
                         tplus, conx1, tppmgs, gni, gnk)
   else
      dsoa = 0.0
   end if

   call msg_toall(CHM_MSG_DEBUG, 'mach_saprc_main [END]')
   if (chm_timings_L) call timing_stop_omp(334)

   !-----------------------------------------------------------------
   return
end subroutine mach_saprc_main
