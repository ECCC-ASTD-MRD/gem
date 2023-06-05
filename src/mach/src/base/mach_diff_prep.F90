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
! Fichier/File   : mach_diff_prep.ftn90
! Creation       : Paul Makar, Deji Akingunola, Feb 2016
! Description    : Interface for calling the GEM-MACH vertical diffusion process
!                  with (optional) multiple diffusion coefficient profiles
!
! Extra info     :
!
! Arguments:
!           IN
!              metvar2d  -> Selected 2D meteorological fields
!              metvar3d  -> Selected 3D meteorological fields
!
!           IN/OUT
!              busper      -> Permanent bus
!              busvol      -> Volatile bus
!              chem_tr     -> Tracers concentrations
!              tracers_can -> Tracers concnetration within canopy layers
!
!=============================================================================
!
!!if_on
subroutine mach_diff_prep(busper, busvol, chem_tr, metvar2d, metvar3d, &
                          metvar3dnocan, metvar3dcan, emisbio_can, tracers_can, &
                          kmod, kcan, imod, imod2, ni_can, ni_nocan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkc, nkt
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_pkg_gas_mod,     only: num_be_sp
!!if_off
   use chm_metvar_mod,       only: MV2D_PPLUS, MV2D_DXDY, MV3D_HUPLUS, &
                                   MV3D_ZPLUS, MV3D_ZMOM, MV3D_SIGM,   &
                                   MV3D_SIGT,  MV3D_RHO,  MV3D_TPLUS, MV3D_KT
   use chm_utils_mod,        only: ik, chm_timestep, CHM_MSG_DEBUG
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_KTN
   use chm_nml_mod,          only: chm_timings_L, chm_vit_l
   use mach_headers_mod,     only: mach_diffusion, mach_vit_diffusivity

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ni_can, ni_nocan
   integer(kind=4), intent   (in) :: kmod(ni_can, chm_nk), kcan(ni_can, nkc)
   integer(kind=4), intent   (in) :: imod(ni_can), imod2(ni_nocan)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in), target :: metvar3dnocan(ni_nocan, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in), target :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   real(kind=4),    intent   (in) :: emisbio_can(ni_can, nkt, num_be_sp)
   real(kind=4),    intent(inout) :: tracers_can(ni_can, nkc, nb_dyn_tracers)
!!if_off
!
!  Declaration of local variables
!
   integer(kind=4) :: i, k, this_ik, sp, ii, ic, kc, kk
   integer(kind=4) :: echoice
   real(kind=4)    :: flux
   real(kind=4)    :: conc(ni_nocan, chm_nk, nb_dyn_tracers)
   real(kind=4)    :: vd  (ni_nocan, nb_dyn_tracers)
   real(kind=4)    :: kt  (ni_nocan, chm_nk)
   real(kind=4)    :: kvkt(ni_nocan, chm_nk)
   real(kind=4)    :: vt  (ni_nocan, chm_nk)
   real(kind=4)    :: dxdy(ni_nocan), psurf(ni_nocan)
   real(kind=4), dimension(ni_nocan, chm_nk, nb_dyn_tracers) :: conc0,    &
                                                                conc_mob, &
                                                                conc_vit
   real(kind=4)    :: conc_can(ni_can, nkt, nb_dyn_tracers)
   real(kind=4)    :: vd_can  (ni_can, nb_dyn_tracers)
   real(kind=4)    :: kt_can  (ni_can, nkt)
   real(kind=4)    :: vt_can  (ni_can, nkt)
   real(kind=4)    :: kvkt_can(ni_can, nkt)
   real(kind=4)    :: kvkt_can0(ni_can, chm_nk), zmomcan(ni_can, chm_nk)
   real(kind=4)    :: dxdy_can(ni_can), psurf_can(ni_can)
   real(kind=4)    :: zcan(ni_can, nkc)
   real(kind=4), dimension(ni_can, nkt, nb_dyn_tracers) :: conc0_can,    &
                                                           conc_ant_can, &
                                                           conc_mob_can, &
                                                           conc_vit_can, &
                                                           conc_bio_can
   logical(kind=4) :: calculate(ni_can,nb_dyn_tracers)
   real(kind=4), dimension(:, :), pointer :: tt, hu, rho, sigm, sigt
   real(kind=4), dimension(:, :), pointer :: zmom, zplus

   real(kind=4), parameter :: min_value = tiny(1.0) !1.0E-12 !tiny(1.0)
!
!  Declaration of external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp, mfotvt
!
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_diffusion [BEGIN]')
   if (chm_timings_L) call timing_start_omp(325, 'mach_diffusion', 480)
!
   if (ni_nocan > 0) then
!
      nullify(tt, hu, rho, sigm, sigt, zplus, zmom)
      tt   (1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_TPLUS)
      hu   (1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_HUPLUS)
      rho  (1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_RHO)
      sigt (1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_SIGT)
      sigm (1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_SIGM)
      zmom (1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_ZMOM)
      zplus(1:ni_nocan, 1:chm_nk) => metvar3dnocan(:, :, MV3D_ZPLUS)
!
      call mfotvt(vt, tt, hu, ni_nocan, chm_nk, ni_nocan)
!
      do i = 1, ni_nocan
         ii = imod2(i)
         psurf(i) = metvar2d(ii, MV2D_PPLUS)
         dxdy (i) = metvar2d(ii, MV2D_DXDY)
      end do

      do k = 1, chm_nk
         do i = 1, ni_nocan
            ii = imod2(i)
            kt(i, k) = busvol(sm(sp_KTN) % out_offset + ik(ii, k, chm_ni))
!
!  Construct the array of species concentrations and the vertical diffusivities
            do sp = 1, nb_dyn_tracers
               conc(i, k, sp) = chem_tr(ii, k, sp)
               conc_mob(i, k, sp) = conc(i, k, sp)
               conc0(i, k, sp) = conc(i, k, sp)
            end do
         end do
      end do
!
!  Deposition velocities placed in deposition velocity array. Deposition
!  velocities for PMs are handled by the CAM package.
      vd = 0.0
      do sp = 1, nb_dyn_tracers
         if (sm(sp) % vd_offset > 0) then
            do i = 1, ni_nocan
               vd(i, sp) = busvol(sm(sp) % vd_offset + imod2(i) - 1)
            end do
!
!  Depostion flux diagnostic for gas
            if (sm(sp) % dd_offset > 0) then
               do i = 1, ni_nocan
                  ii = imod2(i) - 1
                  flux = chm_timestep * rho(i, chm_nk) * conc(i, chm_nk, sp) * &
                         vd(i, sp) * 1.e-6 / sm(sp) % mol_wt
                  busper(sm(sp) % dd_offset + ii) = - flux + &
                                         busper(sm(sp) % dd_offset + ii)
               end do
            end if
         end if
!
      end do
!
!  Non on-road mobile anthropogenic/biogenic area sources.
      echoice = 3
      if (chm_vit_l) echoice = 1
      call mach_diffusion(busper, busvol, conc, vd, psurf, dxdy, rho, kt, vt, &
           sigm, sigt, zplus, zmom, echoice, imod2, ni_nocan, chm_nk)
!
!
!  on-road mobile anth area sources with VIT KT:
!  If VIT is enabled, apply the VIT KT only to mobile emissions taking place
!  on roadways.
      if (chm_vit_l) then
         conc_vit = 0.0
!
         echoice = 0
!  (1) Zero emissions solution with original diffusivities.
         call mach_diffusion(busper, busvol, conc0, vd, psurf, dxdy, rho, kt, &
              vt, sigm, sigt, zplus, zmom, echoice, imod2, ni_nocan, chm_nk)

! Evaluate the vehicle induced turbulence diffusivities (kvkt)
         call mach_vit_diffusivity(busper, dxdy, zmom, kvkt, imod2, ni_nocan)
! Add up the VIT and meteorological induced KT
         kt = kt + kvkt
         echoice = 2
         call mach_diffusion(busper, busvol, conc_mob, vd, psurf, dxdy, rho, kt,&
              vt, sigm, sigt, zplus, zmom, echoice, imod2, ni_nocan, chm_nk)
!
!  Calculate contribution due to on-road mobile anth area sources:
         do sp = 1, nb_dyn_tracers
            if (sm(sp) % mae_offset > 0) then
               do i = 1, ni_nocan
!  Calculate increment only if species had non-zero emissions.
                  if (busper(sm(sp) % mae_offset + imod2(i) - 1) > 0.0) then
                     do k = 1, chm_nk
                        conc_vit(i, k, sp) = conc_mob(i, k, sp) - conc0(i, k, sp)
                     end do
                  end if
               end do
            end if
         end do
!  Thus, the conc_vit array at this point contains the change in concentration
!  due to vehicle emissions, using the VIT KT
!
!  Add the changes in column mass mixing ratio due to VIT and onroad emissions,
!  to the solution which includes biogenic emissions and other anthropogenic
!  emissions.
!
         conc = conc_vit + conc
      end if
!
!  Return the new chemical concentrations back into the Dynamic bus
!
      do sp = 1, nb_dyn_tracers
         do k = 1, chm_nk
            do i = 1, ni_nocan
               chem_tr(imod2(i), k, sp) = max(conc(i, k, sp), min_value)
            end do
         end do
      end do
   end if
!
!********************************************************************
!  Canopy columns cases (canopy model switched on):
!********************************************************************
   if (ni_can > 0) then
!
      nullify(tt, hu, rho, sigm, sigt, zplus, zmom)
      tt   (1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_TPLUS)
      hu   (1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_HUPLUS)
      rho  (1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_RHO)
      sigt (1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_SIGT)
      sigm (1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_SIGM)
      zmom (1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_ZMOM)
      zplus(1:ni_can, 1:nkt) => metvar3dcan(:, :, MV3D_ZPLUS)
!
      call mfotvt(vt_can, tt, hu, ni_can, nkt, ni_can)
!
      do ic = 1, ni_can
         ii = imod(ic)
         psurf_can(ic) = metvar2d(ii, MV2D_PPLUS)
         dxdy_can (ic) = metvar2d(ii, MV2D_DXDY)
!
         do k = 1, chm_nk
            kk = kmod(ic, k)
         ! zmomcan on resolved model levels for evaluating vkt diffusivities
            zmomcan(ic, k) = metvar3d(ic, k, MV3D_ZMOM)
         ! original diffusivities and species concentration on canopy levels
            kt_can(ic, kk) = busvol(sm(sp_KTN) % out_offset + ik(ii, k, chm_ni))
            do sp = 1, nb_dyn_tracers
               conc_can(ic, kk, sp)     = chem_tr(ii, k, sp)
               conc0_can(ic, kk, sp)    = conc_can(ic, kk, sp)
               conc_mob_can(ic, kk, sp) = conc_can(ic, kk, sp)
               conc_bio_can(ic, kk, sp) = conc_can(ic, kk, sp)
            end do
         end do
      end do

!  Evaluate the vehicle induced turbulence diffusivities (kvkt)
      kvkt_can0 = 0.0
      kvkt_can  = 0.0
      if (chm_vit_l) then
         call mach_vit_diffusivity(busper, dxdy_can, zmomcan, kvkt_can0, imod, &
                                   ni_can)
      end if
!
      do ic = 1, ni_can
         ii = imod(ic)
         do kc = 1, nkc
            kk = kcan(ic, kc)
            zcan(ic, kc) = zplus(ic, kk)
            do sp = 1, nb_dyn_tracers
               conc_can(ic, kk, sp)     = tracers_can(ic, kc, sp)
               conc0_can(ic, kk, sp)    = conc_can(ic, kk, sp)
               conc_mob_can(ic, kk, sp) = conc_can(ic, kk, sp)
               conc_bio_can(ic, kk, sp) = conc_can(ic, kk, sp)
            end do
!
!  For anthropogenic KT, the canopy layer value is assumed to be the same as
!  the resolved layer in which it resides.
            do k = chm_nk - 5, chm_nk
!  If the canopy layer kc falls within a given resolved scale layer:
               if (zcan(ic, kc) < zmomcan(ic, k-1) .and. &
                   zcan(ic, kc) >= zmomcan(ic, k)) then
                   this_ik = ik(ii, k-1, chm_ni)
                   kt_can(ic, kk) = busvol(sm(sp_KTN) % out_offset + this_ik)
                   kvkt_can(ic, kk) = kvkt_can0(ic, k-1)
               end if
            end do
!  If the canopy layer kc falls below the lowest resolved scale layer...
            if (zcan(ic,kc) < zmomcan(ic, chm_nk)) then
               this_ik = ik(ii, chm_nk, chm_ni)
               kt_can(ic, kk) = busvol(sm(sp_KTN) % out_offset + this_ik)
               kvkt_can(ic, kk) = kvkt_can0(ic, chm_nk)
            end if
!
         end do
      end do
!
!  Deposition velocities placed in deposition velocity array.
      vd_can = 0.0
      do sp = 1, nb_dyn_tracers
         if (sm(sp) % vd_offset > 0) then
            do ic = 1, ni_can
               vd_can(ic, sp) = busvol(sm(sp) % vd_offset + imod(ic) - 1)
            end do
!
!  Depostion flux diagnostic for gas
            if (sm(sp) % dd_offset > 0) then
               do ic = 1, ni_can
                  ii = imod(ic) - 1
                  flux = chm_timestep * rho(ic, nkt) * conc_can(ic, nkt, sp) * &
                         vd_can(ic, sp) * 1.e-6 / sm(sp) % mol_wt
                  busper(sm(sp) % dd_offset + ii) = -flux + &
                                         busper(sm(sp) % dd_offset + ii)
               end do
            end if
         end if
!
      end do
!
!  (1) Zero emissions solution with original diffusivities.
      echoice = 0
      call mach_diffusion(busper, busvol, conc0_can, vd_can, psurf_can, &
           dxdy_can, rho, kt_can, vt_can, sigm, sigt, zplus, zmom, echoice, &
           imod, ni_can, nkt)

!  Non on-road mobile anth area sources.  Here, the assumption has been made that
!  these sources are being emitted in small towns and clearings, not directly
!  under the forest canopy, with the original resolved scale diffusivities.
!
      echoice = 3
      if (chm_vit_l) echoice = 1
      call mach_diffusion(busper, busvol, conc_can, vd_can, psurf_can, dxdy_can,&
           rho, kt_can, vt_can, sigm, sigt, zplus, zmom, echoice, imod, ni_can, &
           nkt)
! Calculate contribution due to anth (non-mobile) area sources:
      conc_ant_can = 0.0
!
!  The following loop goes through the different options for which the subsequent increment
!  in concentration should be calculated, which is then used to calculate conc_ant_can; The
!  increment in concentration due to ae_offset, fae_offset, and be_offset emissions if VIT
!  is "on", and ae_offset, fae_offset, be_offset, and mae_offset if VIT is "off".
!
      calculate = .false.
      do sp = 1, nb_dyn_tracers
!  Check to see if ae_offset, fae_offset or be_offset species are present and have
!  non-zero emissions.  These are needed for both vit and non-vit simulations.
        if (sm(sp) % ae_offset > 0) then
          do ic = 1, ni_can
            if(busper(sm(sp) % ae_offset + imod(ic)  -  1) > 0.0) calculate(ic,sp) = .true.
          end do
        end if
        if(sm(sp) % fae_offset > 0) then
          do ic = 1, ni_can
            if(busper(sm(sp) % fae_offset + imod(ic) - 1) > 0.0) calculate(ic,sp) = .true.
          end do
        end if
        if(sm(sp) % be_offset > 0) then
          do ic = 1, ni_can
            if(busper(sm(sp) % be_offset + imod(ic) - 1) > 0.0) calculate(ic,sp) = .true.
          end do
        end if
!  Check for mae_offset emissions being present ? only needed if VIT option is not used,
!  and mobile area source is being added as normal surface flux with no KT changes:
        if((.not. chm_vit_l) .and. (sm(sp) % mae_offset > 0)) then
          do ic = 1, ni_can
            if(busper(sm(sp) % mae_offset + imod(ic)  -  1) > 0.0) calculate(ic,sp) = .true.
          end do
         end if
      end do

      do sp = 1,nb_dyn_tracers
!  Calculate increment only if species had non-zero emissions for the cases noted above.
         do ic = 1,ni_can
           if(calculate(ic,sp)) then
                  do kk = 1, nkt
                     conc_ant_can(ic, kk, sp) = conc_can(ic, kk, sp) - &
                                               conc0_can(ic, kk, sp)
                  end do
           end if
         end do
      end do
!
!  on-road mobile anth area sources with VIT KT:
!  If VIT is enabled, apply the VIT KT only to mobile emissions taking place
!  on roadways.
      conc_vit_can = 0.0
      if (chm_vit_l) then
!
!  Add up the VIT and meteorological induced KT
         kvkt_can = kt_can + kvkt_can
         echoice = 2
         call mach_diffusion(busper, busvol, conc_mob_can, vd_can, psurf_can, &
              dxdy_can, rho, kvkt_can, vt_can, sigm, sigt, zplus, zmom, echoice,&
              imod, ni_can, nkt)
!
!  Calculate contribution due to on-road mobile anth area sources:
         do sp = 1, nb_dyn_tracers
            if (sm(sp) % mae_offset > 0) then
               do ic = 1, ni_can
!  Calculate increment only if species had non-zero emissions.
                  if (busper(sm(sp) % mae_offset + imod(ic) - 1) > 0.0) then
                     do kk = 1, nkt
                        conc_vit_can(ic, kk, sp) = conc_mob_can(ic, kk, sp) - &
                                                   conc0_can(ic, kk, sp)
                     end do
                  end if
               end do
            end if
         end do
!  Thus, the conc_vit array at this point contains the change in concentration
!  due to vehicle emissions, using the VIT KT
      end if
!
!  Last diffusion call:  biogenic emissions only , with modified canopy diffusivities
!
      do kk = 1, nkt    ! For canopy grids
         do ic = 1, ni_can
            kt_can(ic, kk) = metvar3dcan(ic, kk, MV3D_KT)
         end do
      end do
!  The canopy emissions are included as diffusion boundary conditions into the
!  individual model layers:
      echoice = -1
      call mach_diffusion(busper, busvol, conc_bio_can, vd_can, psurf_can, &
           dxdy_can, rho, kt_can, vt_can, sigm, sigt, zplus, zmom, echoice,&
           imod, ni_can, nkt, emisbio_can)
!
      conc_can = conc_ant_can + conc_bio_can + conc_vit_can
!
!  Return the new chemical concentrations back into the bus(es)
      do k = 1, chm_nk
         do ic = 1, ni_can
            kk = kmod(ic, k)
            ii = imod(ic)
            do sp = 1, nb_dyn_tracers
               chem_tr(ii, k, sp) = max(conc_can(ic, kk, sp), min_value)
            end do
         end do
      end do
!
      do kc = 1, nkc
         do ic = 1, ni_can
            kk = kcan(ic, kc)
            do sp = 1, nb_dyn_tracers
               tracers_can(ic, kc, sp) = max(conc_can(ic, kk, sp), min_value)
            end do
         end do
      end do

   end if
!
   call msg_toall(CHM_MSG_DEBUG, 'mach_diffusion [END]')
   if (chm_timings_L) call timing_stop_omp(325)
   !-----------------------------------------------------------------
   return

end subroutine mach_diff_prep
