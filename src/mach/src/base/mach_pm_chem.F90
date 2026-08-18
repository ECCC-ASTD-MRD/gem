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
! Fichier/File   : mach_pm_chem.ftn90
! Creation       : P. Huang, W. Gong, Mar. 2008 for adapting AURAMS version CAM
!                  P. Huang, Dec. 2007 for GEM-MACH
!
! Description    : Entry point for aerosol processes
!
! Extra info     : Including calculations of Sea-salt surface flux
!
! Arguments:  IN
!               istep     -> Flag for first chem. step in current
!               f_j       -> Slice number
!               iseasn    -> Seasonal categories for calculation of PM dry deposition
!               busvol    -> Volatile bus for chemistry
!               metvar3d(MV3D_FTOT)   -> Total cloud fraction at each layer (0.0-1.0)
!               metvar3d(MV3D_QCPLUS) -> Cloud water/ice (kg/kg) air
!               metvar3d(MV3D_TPLUS)  -> Temperature (K)
!               metvar3d(MV3D_HUPLUS) -> Specific humidity (kg H2O/kg air)
!               metvar3d(MV3D_RNFLX)  -> Liquid precipitation flux
!               metvar3d(MV3D_SNOFLX) -> Solid precipitation flux
!               metvar3d(MV3D_PPRO)   -> CLOUD TO RAIN COLL TEND (CONSUN + KFC)
!               metvar3d(MV3D_PEVP)   -> EVAP. OF PRECIP (CONSUN + KFC)
!               metvar2d(MV2D_WSDIAG) -> Screen level wind speed

!
!             IN/OUT
!               chem_tr               -> Chemical tracers concentrations (ug/kg)
!               busper                -> Permanent bus for chemistry
!
!             LOCAL
!               pres  -> 3-D pressure (Pa)
!               rho   -> 3-D air density array (kg/m^3)
!               sig   -> Local sigma values
!               fland -> Landuse
!               ra    -> Aerodynamic resistance
!
!============================================================================
!
!!if_on
subroutine mach_pm_chem(busvol, busper, chem_tr, metvar2d, pmetvar3d, fland, &
                        oldso4, iseasn, istep, f_j, pni, pnk, pnka, nmod, &
                        kmod, tracers_can, kcan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkc
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_metvar_mod,       only: mv3d_qcplus, mv3d_tplus, mv3d_huplus,  &
                                   mv3d_sigt,   mv3d_ftot,  mv3d_rnflx,   &
                                   mv3d_snoflx, mv3d_ppro,  mv3d_pevp,    &
                                   mv2d_pplus,  mv2d_ue,    mv2d_tdiag,   &
                                   mv3d_rho,    mv2d_wsdiag, mv3d_ncplus
   use chm_utils_mod,        only: chm_timestep, chm_lun_out, ik,         &
                                   global_debug, chm_error_l, CHM_MSG_DEBUG
   use mach_cam_headers_mod, only: mach_cam_main, mach_cam_sfss, mach_mie_opt
   use chm_consphychm_mod,   only: grav, consth
   use chm_nml_mod,          only: chm_soa_s, chm_seaflux_s, chm_strato_s, &
                                   chm_diag_wetdep_l, nk_start_pm, chm_timings_l, &
                                   chm_diag_aero_opt_l
   use chm_species_info_mod, only: sm, unassigned
   use chm_species_idx_mod,  only: sp_AERO, sp_CAT, sp_H2O2, sp_HCO3, sp_HION, &
                                   sp_WH2O, sp_WSO3, sp_HNO3, sp_LU15, sp_NH3, &
                                   sp_NWOC, sp_O3, sp_SO2, sp_SO4
   use mach_cam_utils_mod,   only: igs_SO2, igs_SO4, igs_O3, igs_H2O2, igs_HNO3, &
                                   igs_ROOH, igs_NH3, sp_gOOH, iae_SS, ip_wflx
   use mach_cam_utils_mod,   only: isize, icom, ntr, nswdep, iae1, iae2, &
                                   sigmathck, chk_trc
   use linoz_param,          only: hu_linoz_tropo
   implicit none

!!if_on
   integer(kind=4), intent   (in) :: istep
   integer(kind=4), intent   (in) :: f_j
   integer(kind=4), intent   (in) :: pnk, pnka
   integer(kind=4), intent   (in) :: pni
   integer(kind=4), intent   (in) :: nmod     (pni)
   integer(kind=4), intent   (in) :: kmod     (pni, chm_nk)
   integer(kind=4), intent   (in) :: iseasn   (chm_ni)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr  (chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: oldso4   (pni, pnka)
   real(kind=4),    intent   (in) :: metvar2d (chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: pmetvar3d(pni, pnka, SIZE_MV3D)
   real(kind=4),    intent   (in) :: fland    (chm_ni, lucprm)
   real(kind=4),    intent(inout), optional :: tracers_can(pni, nkc, nb_dyn_tracers)
   integer(kind=4), intent   (in), optional :: kcan    (pni, nkc)
!!if_off
!
!  Declare local variables
!
   real(kind=4),    parameter :: rso2oh2so4 = -64.6 / 98.1
   real(kind=4),    parameter :: smf = 1.0e-15
   real(kind=4)               :: hu_linoz_consth
   integer(kind=4)            :: l, isz, nsp, iae, iwf
   integer(kind=4)            :: ii, kk, kc, k2, this_ik, ix, kz
   real(kind=4)               :: c_so4
   integer(kind=4)            :: iseasnx(pni)
   real(kind=4)               :: surfwd (pni)
   real(kind=4)               :: tdiag  (pni)
   real(kind=4)               :: pfland (pni, lucprm)
   real(kind=4)               :: ra     (pni, lucprm)
   real(kind=4)               :: usi    (pni)
   real(kind=4)               :: psurf  (pni)
   real(kind=4)               :: rho    (pni, pnk)
   real(kind=4)               :: pres   (pni, pnk)
   real(kind=4)               :: sig    (pni, pnk)
   real(kind=4)               :: tcld3  (pni, pnk)
   real(kind=4)               :: fctr   (pni, pnk)
   real(kind=4)               :: frevp  (pni, pnk)
   real(kind=4)               :: dsig   (pni, pnk)
   real(kind=4)               :: thlev  (pni, pnk)
   real(kind=4)               :: soa    (pni, pnk)
   real(kind=4)               :: throw  (pni, pnk)
   real(kind=4)               :: rhrow  (pni, pnk)
   real(kind=4)               :: hurow  (pni, pnk)
   real(kind=4)               :: pretrow(pni, pnk, 2)
   real(kind=4)               :: zmlwc  (pni, pnk)
   real(kind=4)               :: rtso2  (pni, pnk)
   real(kind=4)               :: ccn    (pni, pnk)
   real(kind=4)               :: wetflx (pni, nswdep)
!
   real(kind=4)               :: dryflx (pni, icom)
   real(kind=4)               :: rsfrow (pni, isize)

   real(kind=4)               :: trwtrow(pni, pnk, isize)
   real(kind=4)               :: aero_tr(pni, pnk, ntr)

!  for unit convertion
   real(kind=4),    parameter ::  ug2kg = 1.E-9
   real(kind=4),    parameter ::  kg2ug = 1.E+9
!
   logical(kind=4)            :: local_dbg
   logical(kind=4)            :: canopy_columns

   real, external             :: sfohra
!
!  Declare external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_pm_chem [BEGIN]')
   if (chm_timings_L) call timing_start_omp(340, 'mach_pm_chem', 480)
!
! CODE BEGINS HERE
!
   hu_linoz_consth = hu_linoz_tropo / consth

   canopy_columns = .false.
   if (present(kcan)) canopy_columns = .true.

! 2D fields
   do ix = 1, pni
      ii = nmod(ix)
      surfwd(ix) = metvar2d(ii, MV2D_WSDIAG)
      tdiag(ix)  = metvar2d(ii, MV2D_TDIAG)
      iseasnx(ix)= iseasn(ii)
      usi(ix)    = max(metvar2d(ii, MV2D_UE), 0.001)
      psurf(ix)  = metvar2d(ii, MV2D_PPLUS)
!
      do l = 1, lucprm
         pfland(ix, l) = fland(ii, l)
         ra(ix, l) = busvol(sm(sp_LU15) % ra_offset + ik(ii, l, chm_ni))
      end do
   end do

! 3D Fields
   do kk = 1, pnk
      kz = kk + nk_start_pm - 1
      do ix = 1, pni
         rho(ix, kk)  = pmetvar3d(ix, kz, mv3d_rho)
         sig(ix, kk)  = pmetvar3d(ix, kz, mv3d_sigt)
         pres(ix, kk) = sig(ix, kk) * psurf(ix)
      end do
   end do

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   call sigmathck(dsig, sig, pni, pnk)

!   call sigmacal(SHJ, SHTJ, DSHJ, sig)
!>> sg:
!>> after call to sigmacal SHJ contains sigma levels where chemicals are located
!>> SHTJ contains the half-way values, except for SHTJ(1) = SHJ(1)
!>> and SHTJ(pnk + 1) = 1
!>> DSHJ are the thickness in sigma of the layers that have SHJ as their centers.

   if (chm_soa_s /= 'NIL')  then
      do kz = nk_start_pm, chm_nk
         do ix = 1, pni
            kk = kmod(ix, kz) - nk_start_pm + 1
            this_ik = ik(nmod(ix), kz, chm_ni)
            soa(ix, kk) = busvol(sm(sp_NWOC) % out_offset + this_ik) * &
                          ug2kg / chm_timestep     ! ug/kg -> kg/kg/s
         end do
      end do
      if (canopy_columns) then
         do kc = 1, nkc
            do ix = 1, pni
               kk = kcan(ix, kc) - nk_start_pm + 1
               this_ik = ik(nmod(ix), kc, chm_ni)
               soa(ix, kk) = busvol(sm(sp_NWOC) % vd_offset + this_ik) * &
                             ug2kg / chm_timestep !ug/kg->kg/kg/s
            end do
         end do
      end if
   else
      soa = 0.0
   end if
!
!  Calculate the average production rate of SO4 for use in CAM
!  in ug/kg/s converted to SO2 oxidation rate in kg of so2/kg/s for CAM
   c_so4 =  rso2oh2so4 * ug2kg / chm_timestep
!
!  load needed gas and aerosol tracers from busdyn (!ug/kg -> kg/kg)
   do kz = nk_start_pm, chm_nk
      do ix = 1, pni
         k2 = kmod(ix, kz)
         kk = k2 - nk_start_pm + 1
         ii = nmod(ix)
         aero_tr(ix, kk, igs_SO2)  = chem_tr(ii, kz, sp_SO2) * ug2kg
!        aero_tr(ix, kk, igs_SO4)  = chem_tr(ii, kz, sp_SO4) * ug2kg
         !using [SO4]g before gas chemistry for use in CAM
         aero_tr(ix, kk, igs_SO4)  = oldso4(ix, k2) * ug2kg
         aero_tr(ix, kk, igs_O3)   = chem_tr(ii, kz, sp_O3)   * ug2kg
         aero_tr(ix, kk, igs_H2O2) = chem_tr(ii, kz, sp_H2O2) * ug2kg
         aero_tr(ix, kk, igs_HNO3) = chem_tr(ii, kz, sp_HNO3) * ug2kg
         aero_tr(ix, kk, igs_ROOH) = chem_tr(ii, kz, sp_gOOH) * ug2kg
         aero_tr(ix, kk, igs_NH3)  = chem_tr(ii, kz, sp_NH3)  * ug2kg

         do nsp = 1, icom * isize
            aero_tr(ix, kk, (nsp + iae1 - 1)) = chem_tr(ii, kz, (nsp + sp_AERO - 1)) * ug2kg
         end do

         rtso2(ix, kk) = (chem_tr(ii, kz, sp_SO4) - oldso4(ix, k2)) * c_so4
      end do
   end do

   if (canopy_columns) then
      do kc = 1, nkc
         do ix = 1, pni
            k2 = kcan(ix, kc)
            kk = k2 - nk_start_pm + 1
            aero_tr(ix, kk, igs_SO2)  = tracers_can(ix, kc, sp_SO2) * ug2kg
            aero_tr(ix, kk, igs_SO4)  = oldso4(ix, k2) * ug2kg
            aero_tr(ix, kk, igs_O3)   = tracers_can(ix, kc, sp_O3)   * ug2kg
            aero_tr(ix, kk, igs_H2O2) = tracers_can(ix, kc, sp_H2O2) * ug2kg
            aero_tr(ix, kk, igs_HNO3) = tracers_can(ix, kc, sp_HNO3) * ug2kg
            aero_tr(ix, kk, igs_ROOH) = tracers_can(ix, kc, sp_gOOH) * ug2kg
            aero_tr(ix, kk, igs_NH3)  = tracers_can(ix, kc, sp_NH3)  * ug2kg
            do nsp = 1, icom * isize
               aero_tr(ix, kk, (nsp + iae1 - 1)) = tracers_can(ix, kc, (nsp + sp_AERO - 1)) * ug2kg
            end do
!
            rtso2(ix, kk) = (tracers_can(ix, kc, sp_SO4) - oldso4(ix, k2)) * c_so4
         end do
      end do
   end if

   do kk = 1, pnk
      kz = kk + nk_start_pm - 1
      do ix = 1, pni
         thlev(ix, kk)   = dsig(ix, kk) * psurf(ix) / (grav * rho(ix, kk))
         pretrow(ix, kk, 1) = pmetvar3d(ix, kz, MV3D_RNFLX)  !liquid precip.
         pretrow(ix, kk, 2) = pmetvar3d(ix, kz, MV3D_SNOFLX) !solid precip.
         frevp(ix, kk)   = pmetvar3d(ix, kz, MV3D_PEVP)
         throw(ix, kk)   = pmetvar3d(ix, kz, MV3D_TPLUS)
         hurow(ix, kk)   = pmetvar3d(ix, kz, MV3D_HUPLUS)
         rhrow(ix, kk)   = sfohra(hurow(ix, kk), throw(ix, kk), pres(ix, kk))
         rhrow(ix, kk)   = min(max(rhrow(ix, kk), 0.0), 1.0)
         tcld3(ix, kk)   = pmetvar3d(ix, kz, MV3D_FTOT)
          ! zmlwc contains both the stratiform and convective cloud water content
         zmlwc(ix, kk)   = pmetvar3d(ix, kz, MV3D_QCPLUS)
         fctr(ix, kk)    = pmetvar3d(ix, kz, MV3D_PPRO) / (zmlwc(ix, kk) + smf)
         fctr(ix, kk)    = min(1.0, fctr(ix, kk))
          ! This CCN is only used if chm_indirect_l is enabled
         ccn(ix, kk)     = pmetvar3d(ix, kz, MV3D_NCPLUS)
      end do
   end do

   where (zmlwc < 1.0e-7)
      zmlwc = 0.0
      tcld3 = 0.0
   end where

!  call aerosol physics
   call chk_trc(aero_tr, iae1, iae2, f_j, 'bef_camV', istep, pni, pnk)
   if (chm_error_l) return

   call mach_cam_main(f_j, throw, rhrow, aero_tr, istep, zmlwc, pretrow,      &
                      pres, rho, thlev, rtso2, pfland, soa, iseasnx, ra, usi, &
                      tcld3, fctr, frevp, wetflx, dryflx, ccn, trwtrow,       &
                      pni, pnk)
   if (chm_error_l) return

   call chk_trc(aero_tr, iae1, iae2, f_j, 'aft_camV', istep, pni, pnk)
   if (chm_error_l) return

!  surface flux

   select case (chm_seaflux_s)
      case ('GONG_MONAHAN')
         if (local_dbg) then
            write (chm_lun_out, *) 'Compute sea-salt surface flux by cam scheme: ', chm_seaflux_s
         end if

!  sea-salt surface flux
         call mach_cam_sfss(tdiag, surfwd, rsfrow, pfland, pni)
!  adding surfing production flux
         do ix = 1, pni
            do isz = 1, isize
               nsp = (iae_SS - 1) * isize + isz + (iae1 - 1)
               aero_tr(ix, pnk, nsp) =                                  &
                                aero_tr(ix, pnk, nsp) + rsfrow(ix, isz) &
                       / (rho(ix, pnk) * thlev(ix, pnk)) * chm_timestep
            end do
         end do
!
      case ('GONG_MONAHAN_F')
         if (local_dbg) then
            write (chm_lun_out, *) '> Warning '
            write (chm_lun_out, *) '> sea-salt flux treated as emission : ', chm_seaflux_s
         end if
      case default
         if (local_dbg) then
            write (chm_lun_out, *) '> Warning '
            write (chm_lun_out, *) '> No surface flux scheme for sea-salt selected: ', chm_seaflux_s
         end if
   end select
!
   do kz = nk_start_pm, chm_nk
      do ix = 1, pni
         kk = kmod(ix, kz) - nk_start_pm + 1
         ii = nmod(ix)
!
! Apply lid for CAM2 tendencies of 10 ppmv specific humidity
!
         if ((chm_strato_s == 'NIL') .or. (hurow(ix, kk) >= hu_linoz_consth)) then
!           !kg/kg -> ug/kg
            chem_tr(ii, kz, sp_SO2 ) = aero_tr(ix, kk, igs_SO2)  * kg2ug
            chem_tr(ii, kz, sp_SO4 ) = aero_tr(ix, kk, igs_SO4)  * kg2ug
            chem_tr(ii, kz, sp_O3  ) = aero_tr(ix, kk, igs_O3)   * kg2ug
            chem_tr(ii, kz, sp_H2O2) = aero_tr(ix, kk, igs_H2O2) * kg2ug
            chem_tr(ii, kz, sp_HNO3) = aero_tr(ix, kk, igs_HNO3) * kg2ug
            chem_tr(ii, kz, sp_gOOH) = aero_tr(ix, kk, igs_ROOH) * kg2ug
            chem_tr(ii, kz, sp_NH3 ) = aero_tr(ix, kk, igs_NH3)  * kg2ug

            do nsp = 1, icom * isize
               chem_tr(ii, kz, (nsp + sp_AERO - 1)) = aero_tr(ix, kk, (nsp + iae1 - 1)) * kg2ug
            end do
!
         end if      !  if ((chm_strato_s == 'NIL') .or. (hu_vmr(i,k) >= hu_linoz ))
!
      end do
   end do

   if (canopy_columns) then
      do kc = 1, nkc
         do ix = 1, pni
            k2 = kcan(ix, kc)
            kk = k2 - nk_start_pm + 1
            tracers_can(ix, kc, sp_SO2 ) = aero_tr(ix, kk, igs_SO2)  * kg2ug
            tracers_can(ix, kc, sp_SO4 ) = aero_tr(ix, kk, igs_SO4)  * kg2ug
            tracers_can(ix, kc, sp_O3  ) = aero_tr(ix, kk, igs_O3)   * kg2ug
            tracers_can(ix, kc, sp_H2O2) = aero_tr(ix, kk, igs_H2O2) * kg2ug
            tracers_can(ix, kc, sp_HNO3) = aero_tr(ix, kk, igs_HNO3) * kg2ug
            tracers_can(ix, kc, sp_gOOH) = aero_tr(ix, kk, igs_ROOH) * kg2ug
            tracers_can(ix, kc, sp_NH3 ) = aero_tr(ix, kk, igs_NH3)  * kg2ug

            do nsp = 1, icom * isize
               tracers_can(ix, kc, (nsp + sp_AERO - 1)) = aero_tr(ix, kk, nsp + iae1 - 1) * kg2ug
            end do
         end do
      end do
   end if

! Wet deposition diagnostic
   if (chm_diag_wetdep_L) then
      do ix = 1, pni
         ii = nmod(ix) - 1
         busper(sm(sp_WSO3) % wd_offset + ii) = -wetflx(ix, 1) + &
                               busper(sm(sp_WSO3) % wd_offset + ii)
         busper(sm(sp_H2O2) % wd_offset + ii) = -wetflx(ix, 2) + &
                               busper(sm(sp_H2O2) % wd_offset + ii)
         busper(sm(sp_gOOH) % wd_offset + ii) = -wetflx(ix, 3) + &
                               busper(sm(sp_gOOH) % wd_offset + ii)
         busper(sm(sp_CAT)  % wd_offset + ii) = -wetflx(ix, 7) + &
                               busper(sm(sp_CAT)  % wd_offset + ii)
         busper(sm(sp_HCO3) % wd_offset + ii) = -wetflx(ix, 8) + &
                               busper(sm(sp_HCO3) % wd_offset + ii)
         busper(sm(sp_HION) % wd_offset + ii) = -wetflx(ix, 9) + &
                               busper(sm(sp_HION) % wd_offset + ii)
         busper(sm(sp_WH2O) % wd_offset + ii) = -wetflx(ix, 10) + &
                               busper(sm(sp_WH2O) % wd_offset + ii)
      end do
      do iae = 1, icom
         nsp = iae * isize + sp_AERO - 1  ! Use the species last bin's index
         if (sm(nsp) % wd_offset > 0) then
            iwf = ip_wflx(iae)               ! wet dep index
            do ix = 1, pni
               ii = nmod(ix) - 1
               busper(sm(nsp) % wd_offset + ii) = - wetflx(ix, iwf) + &
                                     busper(sm(nsp) % wd_offset + ii)
            end do
         end if
      end do
   end if

! Dry deposition diagnostic for aerosols
   do iae = 1, icom
      nsp = iae * isize + sp_AERO - 1 ! Use the aerosol's last bin index
      if (sm(nsp) % dd_offset > 0) then
         do ix = 1, pni
            ii = nmod(ix) - 1
            busper(sm(nsp) % dd_offset + ii) = busper(sm(nsp) % dd_offset + ii) + &
                                                   dryflx(ix, iae)
         end do
      end if
   end do

! Diagnostic aerosol optical properties (for estimating AOD, etc)
   if (chm_diag_aero_opt_l) then
      call mach_mie_opt(busper, busvol, aero_tr, trwtrow, thlev, rho, pni, pnk, &
                        nmod, kmod)
   end if

   call msg_toall(CHM_MSG_DEBUG, 'mach_pm_chem [END]')
   if (chm_timings_L) call timing_stop_omp(340)
   !-----------------------------------------------------------------

   return
end subroutine mach_pm_chem
