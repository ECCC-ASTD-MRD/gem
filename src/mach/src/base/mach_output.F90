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
! Fichier/File   : mach_output.ftn90
! Creation       : S. Menard  ,  GEM-MACH, April 2008.
! Description    : 1) Units conversion for gases          (ug/kg) --> (ppb)
!                  2) Units conversion for speciated PM   (ug/kg) --> (ug/m3)
!                  3) Compute PM2.5 (AF) and PM10 (AC)    (ug/m3)
!                  4) Compute AQHI25 and AQHI10           (no unit)
!                  5) Compute dry deposition for NOy
!                  6) Add density to diagnostic level
!
! Extra info     : AQHI formulations
! -----------
!                  aqhi(pm2.5) = 10/10.4 * ( 100 * (exp(0.000871 * no2) - 1 + exp(0.000537 * o3) - 1 +
!                              exp(0.000487 * pm2.5) - 1) )
!                  aqhi(pm10)  = 10/11.7 * ( 100 * (exp(0.000871 * no2) - 1 + exp(0.000537 * o3) - 1 +
!                              exp(0.000297 * pm10)   -1) )
! Arguments:
!            IN
!
!               chem_tr  -> Chemical species' concentrations (ug/kg)
!
!            OUT
!               busvol   -> Chemistry volatile bus
!
!==============================================================================================
!
!!if_on
subroutine mach_output(busper, busvol, chem_tr, metvar2d, metvar3d, landuse)

   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_metvar_mod,       only: MV2D_PPLUS, MV2D_TDIAG, MV2D_QDIAG, MV3D_RHO
   use chm_utils_mod,        only: ik, CHM_MSG_DEBUG
   use chm_consphychm_mod,   only: mwt_air, rgasd, delta
   use chm_nml_mod,          only: aerosize, chm_aqhi_l, chm_pkg_pm_s, &
                                   chm_diag_aerosols_l, chm_ss_ao_l, chm_timings_L
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_O3, sp_NO2, sp_NOy, sp_AERO, sp_LU15,    &
                                   sp_AQ25, sp_AQ10, sp_AC, sp_AF, sp_density, &
                                   sp_NAF, sp_NAC, sp_NAT, sp_NUMC
   use mach_cam_utils_mod,   only: isize, icom, iae_SS
   use mach_cam_utils_mod,   only: rhop0, pvol, pm_agrege_nbin
   use mach_pkg_gas_mod,     only: nsp_gas_ppb_out, gas_ppb_out_list, &
                                   NOy_species_list, nb_NOy_species
   use mach_drydep_mod,      only: chm_lu15_out_l

   implicit none
!!if_on
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: landuse (chm_ni, lucprm)
!!if_off
!
!  Local variables
!
   integer(kind=4)            :: i, k, ibin, iae, isp, nsp, this_ik
   integer(kind=4)            :: busid, busid2
   real(kind=4),    parameter :: maxradius_pm25 = 1.280, maxradius_pm10 = 5.120, &
                                 maxradius_pm1 = 0.640, N_mwt = 14.
   real(kind=4)               :: dryflx(chm_ni)
   real(kind=4)               :: factor, o3_part, no2_part, pm25_part, pm10_part
   real(kind=4)               :: rhod(chm_ni)
   real(kind=4)               :: volratpm1, volratpm2p5, volratpm10, aerovol
   real(kind=4)               :: aeronum(chm_ni, chm_nk + 1, isize)
   ! Aerosol species (components) concentration in ug/m3
   real(kind=4), dimension(chm_ni, chm_nk + 1, icom, isize)  :: pm_outv
   real(kind=4), dimension(chm_ni, chm_nk + 1, icom), target :: &
                                               spm1,  & ! Speciated PM1
                                               spm25, & ! Speciated PM2.5
                                               spm10, & ! Speciated PM10
                                               spmt     ! Speciated total PM
   real(kind=4), dimension(:, :, :), pointer :: spm
   character(len=2), dimension(pm_agrege_nbin) :: pm_agrege_dx
!
!  External subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_output [BEGIN]')
   if (chm_timings_L) call timing_start_omp(380, 'mach_output', 480)
!
! Evaluate air density at the model surface (diagnostic level)
   do i = 1, chm_ni
      rhod(i) = metvar2d(i, MV2D_PPLUS) / (rgasd * metvar2d(i, MV2D_TDIAG) * &
                                      (1.0 + delta * metvar2d(i, MV2D_QDIAG)))
   end do

! 1) Units conversion for gases     (ug/kg) --> (ppb)
! =======================================================================================================
! NOy:
   if (sp_NOy > 0) then
      if (sm(sp_NOy) % out_offset > 0) then
         busid = sm(sp_NOy)% out_offset - 1
         do k = 1, chm_nk + 1
            do i = 1, chm_ni
               busid = busid + 1
               busvol(busid) = 0.0
            end do
         end do
!
         do nsp = 1, nb_NOy_species
            isp = NOy_species_list(nsp)
            busid = sm(sp_NOy) % out_offset
            factor = mwt_air / sm(isp) % mol_wt

            do k = 1, chm_nk + 1
               do i = 1, chm_ni
                  this_ik = ik(i, k, chm_ni)
                  busvol(busid + this_ik) = busvol(busid + this_ik) + &
                                            chem_tr(i, k, isp) * factor
               end do
            end do
         end do
      end if
   end if
!
! =============================================================================
! Default: O3, NO, NO2, SO2:
   do nsp = 1, nsp_gas_ppb_out
      isp = gas_ppb_out_list(nsp)
      busid = sm(isp) % out_offset
      factor = mwt_air / sm(isp) % mol_wt

      do k = 1, chm_nk + 1
         do i = 1, chm_ni
            this_ik = ik(i, k, chm_ni)
            busvol(busid + this_ik) = chem_tr(i, k, isp) * factor
         end do
      end do
   end do

! 2) Units conversion for particles (ug/kg) --> (ug/m3)
! =============================================================================
   if (chm_pkg_pm_s(1:3) == 'CAM') then
!
      do ibin = 1, isize
         do iae = 1, icom
            nsp = (iae - 1) * isize + ibin + sp_AERO - 1
            do k = 1, chm_nk
               do i = 1, chm_ni
                  pm_outv(i, k, iae, ibin) = metvar3d(i, k, MV3D_RHO) * &
                                             chem_tr(i, k, nsp)
               end do
            end do

            do i = 1, chm_ni
               pm_outv(i, chm_nk+1, iae, ibin) = pm_outv(i, chm_nk, iae, ibin)
!               pm_outv(i, chm_nk+1, iae, ibin) = rhod(i) * &
!                                                 chem_tr(i, chm_nk + 1, nsp)
            end do
         end do
      end do

! 3)
! a. Compute speciated PM1 (optional), PM2.5, PM10, and PM total (optional)
! ==============================================================================
! aerosize = 0.005, 0.010, 0.020, 0.040, 0.080, 0.160, 0.320, 0.640, 1.280, 2.560,
!            5.120, 10.240, 20.480
! aerosize = 0.005, 1.280, 5.12

! Initializations;
      spm1 = 0.0
      spm25 = 0.0
      spm10 = 0.0
      spmt  = 0.0
      busid = -1
      do k = 1, chm_nk + 1
         do i = 1, chm_ni
            busid = busid + 1
            busvol(sm(sp_AF) % out_offset + busid) = 0.0
            busvol(sm(sp_AC) % out_offset + busid) = 0.0
         end do
      end do
!
!  Value of volrat* is the ratio of volume of spherical shells for the PM:
      do ibin = 1, isize
! PM1 accumulator
         if (aerosize(ibin + 1) < maxradius_pm1) then
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spm1(i, k, iae) = spm1(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin)
                  end do
               end do
            end do
!
! Partial PM1 and PM2.5 accumulators
         else if (aerosize(ibin + 1) == maxradius_pm1) then
            volratpm1 = (0.5**3 - aerosize(ibin)**3) / &
                        (aerosize(ibin+1)**3 - aerosize(ibin)**3)
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spm1(i, k, iae) = spm1(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin) * volratpm1

!  Add PM1 bins to PM2.5 accumulators:
                     spm25(i, k, iae) = spm25(i, k, iae) + spm1(i, k, iae) + &
                                   pm_outv(i, k, iae, ibin) * (1.0 - volratpm1)
                  end do
               end do
            end do
!
! PM2.5 accumulators
         else if ((aerosize(ibin + 1) > maxradius_pm1) .and. &
                  (aerosize(ibin + 1) < maxradius_pm25)) then
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spm25(i, k, iae) = spm25(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin)
                  end do
               end do
            end do
!
! Partial PM2.5 and PM10 accumulators
         else if (aerosize(ibin + 1) == maxradius_pm25) then
!            volratpm2p5 = (1.25**3 - aerosize(ibin)**3) / &
!                           (aerosize(ibin+1)**3 - aerosize(ibin)**3)
            volratpm2p5 = 1.0
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spm25(i, k, iae) = spm25(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin) * volratpm2p5

!  Add PM2.5 bins to PM10 accumulators:
                     spm10(i, k, iae) = spm10(i, k, iae) + spm25(i, k, iae) + &
                                 pm_outv(i, k, iae, ibin) * (1.0 - volratpm2p5)
                  end do
               end do
            end do
!
! PM10 accumulators
         else if ((aerosize(ibin + 1) > maxradius_pm25) .and. &
                  (aerosize(ibin + 1) < maxradius_pm10)) then
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spm10(i, k, iae) = spm10(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin)
                  end do
               end do
            end do
!
! Partial PM10 and PMtot accumulators
         else if (aerosize(ibin + 1) == maxradius_pm10) then
!            volratpm10 = (5.0**3 - aerosize(ibin)**3) / &
!                         (aerosize(ibin+1)**3 - aerosize(ibin)**3)
            volratpm10 = 1.0
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spm10(i, k, iae) = spm10(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin) * volratpm10

!  Add PM10 bins to PMtot accumulators:
                     spmt(i, k, iae) = spmt(i, k, iae) + spm10(i, k, iae) + &
                                 pm_outv(i, k, iae, ibin) * (1.0 - volratpm10)
                  end do
               end do
            end do
!
! PMtot accumulators:
         else
            do iae = 1, icom
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     spmt(i, k, iae) = spmt(i, k, iae) + &
                                        pm_outv(i, k, iae, ibin)
                  end do
               end do
            end do
         end if
      end do ! End ibin loop
!
! 3b) Evaluate PM2.5(AF) and PM10(AC) (by aggregating spm25 and spm10 across species)
!
      do iae = 1, icom
         if (iae == iae_SS .and. (.not. chm_ss_ao_l)) cycle ! SEA-SALT excluded by default setup
         do k = 1, chm_nk + 1
            do i = 1, chm_ni
               this_ik = ik(i, k, chm_ni)
               busvol(sm(sp_AF) % out_offset + this_ik) = &
                     busvol(sm(sp_AF) % out_offset + this_ik) + spm25(i, k, iae)
!
               busvol(sm(sp_AC) % out_offset + this_ik) = &
                     busvol(sm(sp_AC) % out_offset + this_ik) + spm10(i, k, iae)
            end do
         end do
      end do
!
! 3c). Save spm1, spm25, spm10, and spmt to the bus for output, if requested
! Partial PM2.5 and PM10 accumulators
      if (chm_diag_aerosols_l) then
!
         if (isize <= 2) then
            pm_agrege_dx = (/'AF', 'AC'/)
         else
            pm_agrege_dx = (/'AS', 'AF', 'AC', 'AT'/)
         end if
         do ibin = 1, pm_agrege_nbin
            if (pm_agrege_dx(ibin) == 'AS') then
               spm => spm1
            else if (pm_agrege_dx(ibin) == 'AF') then
               spm => spm25
            else if (pm_agrege_dx(ibin) == 'AC') then
               spm => spm10
            else if (pm_agrege_dx(ibin) == 'AT') then
               spm => spmt
            end if
            do iae = 1, icom
               nsp = (iae - 1) * pm_agrege_nbin + ibin + sp_AERO - 1
               busid = sm(nsp) % out_offset
               if (busid > 0) then
                  do k = 1, chm_nk + 1
                     do i = 1, chm_ni
                       busvol(busid + ik(i, k, chm_ni)) = spm(i, k, iae)
                     end do
                  end do
               end if
            end do
         end do
!
! Total aerosol number concentration per bin
         do ibin = 1, isize
            busid = sm(sp_NUMC + ibin - 1) % out_offset
            do k = 1, chm_nk + 1
               do i = 1, chm_ni
                  aerovol = 0.0
                  do iae = 1, icom
                     aerovol = aerovol + pm_outv(i, k, iae, ibin) / &
                                         rhop0(iae) * 1.E-09
                  end do
               ! Calculate the particle number density (# per m3 of air):
                  aeronum(i, k, ibin) = aerovol / pvol(ibin)
                  if (busid > 0) &
                     busvol(busid + ik(i, k, chm_ni)) = aeronum(i, k, ibin)
               end do
            end do
         end do
!
! 3d). FINE, COARSE and TOTAL PM number concentrations
! PM2.5 accumulators
         busid = -1
         do k = 1, chm_nk + 1
            do i = 1, chm_ni
               busid = busid + 1
               busvol(sm(sp_NAF) % out_offset + busid) = 0.0
               busvol(sm(sp_NAC) % out_offset + busid) = 0.0
               busvol(sm(sp_NAT) % out_offset + busid) = 0.0
            end do
         end do
!
         do ibin = 1, isize
            if (aerosize(ibin + 1) < maxradius_pm25) then
               busid = sm(sp_NAF) % out_offset
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     this_ik = ik(i, k, chm_ni)
                     busvol(busid + this_ik) = busvol(busid + this_ik) + &
                                               aeronum(i, k, ibin)
                  end do
               end do
! Partial PM2.5 and PM10 accumulators
            else if (aerosize(ibin + 1) == maxradius_pm25) then
               volratpm2p5 = (1.25**3 - aerosize(ibin)**3) / &
                             (aerosize(ibin+1)**3 - aerosize(ibin)**3)
               busid = sm(sp_NAF) % out_offset
               busid2 = sm(sp_NAC) % out_offset
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     this_ik = ik(i, k, chm_ni)
                     busvol(busid + this_ik) = busvol(busid + this_ik) +   &
                                               aeronum(i, k, ibin) * volratpm2p5
                     busvol(busid2 + this_ik) = busvol(busid2 + this_ik) + &
                                                busvol(busid + this_ik) +  &
                                        aeronum(i, k, ibin) * (1.0 - volratpm2p5)
                  end do
               end do
! PM10 accumulators
            else if ((aerosize(ibin + 1) > maxradius_pm25) .and. &
                     (aerosize(ibin + 1) < maxradius_pm10)) then
               busid = sm(sp_NAC) % out_offset
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     this_ik = ik(i, k, chm_ni)
                     busvol(busid + this_ik) = busvol(busid + this_ik) + &
                                               aeronum(i, k, ibin)
                  end do
               end do
! Partial PM10 and PMtot accumulators
            else if (aerosize(ibin + 1) == maxradius_pm10) then
               volratpm10 = (5.0**3 - aerosize(ibin)**3) / &
                            (aerosize(ibin+1)**3 - aerosize(ibin)**3)
               busid = sm(sp_NAC) % out_offset
               busid2 = sm(sp_NAT) % out_offset
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     this_ik = ik(i, k, chm_ni)
                     busvol(busid + this_ik) = busvol(busid + this_ik) + &
                                           aeronum(i, k, ibin) * volratpm10
                     busvol(busid2 + this_ik) = busvol(busid2 + this_ik) + &
                                               busvol(busid + this_ik) +   &
                                        aeronum(i, k, ibin) * (1.0 - volratpm10)
                  end do
               end do
! PMtot accumulators:
            else
               busid = sm(sp_NAT) % out_offset
               do k = 1, chm_nk + 1
                  do i = 1, chm_ni
                     this_ik = ik(i, k, chm_ni)
                     busvol(busid + this_ik) = busvol(busid + this_ik) + &
                                               aeronum(i, k, ibin)
                  end do
               end do
            end if
         end do ! End ibin loop

      end if !(chm_diag_aerosols_l)
!
! 4) Calculation of AQHI25 and AQHI10 (with O3 and NO2 concentrations in ppb)
! ===============================================================================
      if (chm_aqhi_l) then
         do k = 1, chm_nk + 1
            do i = 1, chm_ni
               this_ik = ik(i, k, chm_ni)
               no2_part = exp(0.000871 * chem_tr(i, k, sp_NO2) * mwt_air / &
                                        sm(sp_NO2) % mol_wt) - 1.0
               o3_part  = exp(0.000537 * chem_tr(i, k, sp_O3) * mwt_air / &
                                        sm(sp_O3) % mol_wt) - 1.0
               pm25_part = exp(0.000487 * busvol(sm(sp_AF) % out_offset + this_ik)) - 1.0
               pm10_part = exp(0.000297 * busvol(sm(sp_AC) % out_offset + this_ik)) - 1.0

               busvol(sm(sp_AQ25) % out_offset + this_ik) = (10. / 10.4) * &
                                      100. * (no2_part + o3_part + pm25_part)

               busvol(sm(sp_AQ10) % out_offset + this_ik) = (10. / 11.7) * &
                                      100. * (no2_part + o3_part + pm10_part)
            end do
         end do
      end if !chm_aqhi_l
!
   end if !chm_pkg_pm_s
!
! 5) Calculation of dry deposition of NOy
! =============================================================================

   if (sp_NOy > 0) then
      if (sm(sp_NOy) % dd_offset > 0) then
         dryflx = 0.0
         do nsp = 1, nb_NOy_species
            isp = NOy_species_list(nsp)
      !Only include NOy species with defined deposition field
            if (sm(isp) % dd_offset <= 0) cycle
            do i = 1, chm_ni
               dryflx(i) = dryflx(i) + busper(sm(isp) % dd_offset + i -1)
            end do
         end do
         do i = 1, chm_ni
            busper(sm(sp_NOy) % dd_offset + i -1) = dryflx(i)
         end do
      end if
   end if

! 6) Diagnostic level for density
! ==============================================================================
   if (sp_density > 0) then
      do k = 1, chm_nk
         do i = 1, chm_ni
            busvol(sm(sp_density) % out_offset + ik(i, k, chm_ni)) = &
                                                     metvar3d(i, k, MV3D_RHO)
         end do
      end do
      do i = 1, chm_ni
         busvol(sm(sp_density) % out_offset + ik(i, chm_nk + 1, chm_ni)) = &
                                                                        rhod(i)
      end do
   end if

! 6) Landuse output
! ==============================================================================
   if (chm_lu15_out_l) then
      do k = 1, lucprm
         do i = 1, chm_ni
            busper(sm(sp_LU15) % per_offset + ik(i, k, chm_ni)) = landuse(i, k)
         end do
      end do
   end if

   call msg_toall(CHM_MSG_DEBUG, 'mach_output [END]')
   if (chm_timings_L) call timing_stop_omp(380)
   !-----------------------------------------------------------------

   return
end subroutine mach_output
