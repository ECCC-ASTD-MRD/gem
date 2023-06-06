!begin trap head
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
! Fichier/File   : mach_headers_mod.ftn90
! Creation       : H. Landry, S. Menard, Mai 2008
! Description    : Modules defining explicit interfaces for mach subroutines
!
! Extra info     :
!
!============================================================================

module mach_headers_mod
   interface
!end trap head

subroutine mach_biog_beis(emissbio, emissbio_can, summer_std_p, winter_std_p, &
                          sesn_coef, air_temp, zmomcan, light_correction,     &
                          laifrac, canopy_hgt, canopy_column, ni_can)
   use chm_ptopo_grid_mod,   only: chm_ni, nkt, nkc
   use mach_pkg_gas_mod,     only: num_be_std, num_be_sp
   integer(kind=4), intent (in) :: ni_can
   real(kind=4),    intent(out) :: emissbio        (num_be_sp, chm_ni)
   real(kind=4),    intent(out) :: emissbio_can    (ni_can, nkt, num_be_sp)
   real(kind=4),    intent (in) :: winter_std_p    (num_be_std, chm_ni)
   real(kind=4),    intent (in) :: summer_std_p    (num_be_std, chm_ni)
   real(kind=4),    intent (in) :: light_correction(nkc + 1, chm_ni)
   real(kind=4),    intent (in) :: air_temp        (nkc + 1, chm_ni)
   real(kind=4),    intent (in) :: laifrac         (nkc + 1, chm_ni)
   real(kind=4),    intent (in) :: canopy_hgt      (ni_can)
   real(kind=4),    intent (in) :: zmomcan         (ni_can, nkt)
   real(kind=4),    intent (in) :: sesn_coef       (chm_ni)
   logical(kind=4), intent (in) :: canopy_column   (chm_ni)
end subroutine mach_biog_beis

subroutine mach_biog_main(busper, metvar2d, metvar3d, iseasons,  metvar3dcan, &
                          emisbio_can, kcan, ni_can)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use mach_pkg_gas_mod,     only: num_be_sp
   integer(kind=4), intent   (in) :: ni_can
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: metvar2d   (chm_ni,SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d   (chm_ni,chm_nk,SIZE_MV3D)
   real(kind=4),    intent   (in) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   integer(kind=4), intent   (in) :: iseasons   (chm_ni)
   integer(kind=4), intent   (in) :: kcan       (ni_can, nkc)
   real(kind=4),    intent  (out) :: emisbio_can(ni_can, nkt, num_be_sp)
end subroutine mach_biog_main

subroutine mach_biog_parcalc(cosine_solar_p, PARDB, PARDIF, metvar2d)
   use chm_ptopo_grid_mod, only: chm_ni
   use chm_metvar_mod,     only: SIZE_MV2D
   real(kind=4), intent (in) :: cosine_solar_p(chm_ni)           ! cosine of solar angle
   real(kind=4), intent(out) :: pardb         (chm_ni)           ! PAR direct beam
   real(kind=4), intent(out) :: pardif        (chm_ni)           ! PAR diffuse!
   real(kind=4), intent (in) :: metvar2d      (chm_ni,SIZE_MV2D) ! Met. variables
end subroutine mach_biog_parcalc

subroutine mach_calc_diag(busvol, busper, chem_tr, metvar2d, metvar3d, step)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   integer(kind=4), intent(in)    :: step
   real   (kind=4), dimension(:), pointer, contiguous :: busper
   real   (kind=4), dimension(:), pointer, contiguous :: busvol
   real   (kind=4), intent   (in) :: chem_tr(chm_ni, chm_nk+1, nb_dyn_tracers)
   real   (kind=4), intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real   (kind=4), intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
end subroutine mach_calc_diag

subroutine mach_calc_season(metvar2d, iseasn)
   use chm_ptopo_grid_mod,   only: chm_ni
   use chm_metvar_mod,       only: SIZE_MV2D
   integer(kind=4), intent(out) :: iseasn  (chm_ni)
   real(kind=4),    intent (in) :: metvar2d(chm_ni, SIZE_MV2D)
end subroutine mach_calc_season

subroutine mach_diffusion(busper, busvol, conc, vd, psurf, dxdy, rho, kt, vt, &
                          sigm, sigt, zplus, zmom, echoice, nmod, dni, dnk, &
                          emisbio2D)
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_pkg_gas_mod,     only: num_be_sp
   integer(kind=4), intent   (in) :: dni, dnk
   integer(kind=4), intent   (in) :: nmod(dni)
   integer(kind=4), intent   (in) :: echoice
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: conc     (dni, dnk, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: vd       (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: kt       (dni, dnk)
   real(kind=4),    intent   (in) :: rho      (dni, dnk)
   real(kind=4),    intent   (in) :: vt       (dni, dnk)
   real(kind=4),    intent   (in) :: sigm     (dni, dnk)
   real(kind=4),    intent   (in) :: sigt     (dni, dnk)
   real(kind=4),    intent   (in) :: zplus    (dni, dnk)
   real(kind=4),    intent   (in) :: zmom     (dni, dnk)
   real(kind=4),    intent   (in) :: psurf    (dni)
   real(kind=4),    intent   (in) :: dxdy     (dni)
   real(kind=4),    intent   (in), optional :: emisbio2D(dni, dnk, num_be_sp)
end subroutine mach_diffusion

subroutine mach_diff_boundary(conc, emissions, vd, ubf, kt, rho, zf, zh, &
                              dxdy, dni, dnk)
   use chm_species_info_mod, only: nb_dyn_tracers
   integer(kind=4), intent   (in) :: dni, dnk
   real(kind=4),    intent(inout) :: conc     (dni, dnk, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: emissions(dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: vd       (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: ubf      (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: kt       (dni, dnk)
   real(kind=4),    intent   (in) :: rho      (dni, dnk)
   real(kind=4),    intent   (in) :: zf       (dni, dnk)
   real(kind=4),    intent   (in) :: zh       (dni, dnk)
   real(kind=4),    intent   (in) :: dxdy     (dni)
end subroutine mach_diff_boundary

subroutine mach_diff_flux(conc, emissions, vd, ubf, kt, rho, zf, zh, dxdy, &
                          dni, dnk)
   use chm_species_info_mod, only: nb_dyn_tracers
   integer(kind=4), intent   (in) :: dni, dnk
   real(kind=4),    intent(inout) :: conc     (dni, dnk, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: emissions(dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: vd       (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: ubf      (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: kt       (dni, dnk)
   real(kind=4),    intent   (in) :: rho      (dni, dnk)
   real(kind=4),    intent   (in) :: zf       (dni, dnk)
   real(kind=4),    intent   (in) :: zh       (dni, dnk)
   real(kind=4),    intent   (in) :: dxdy     (dni)
end subroutine mach_diff_flux

subroutine mach_diff_rpnphy(conc, emissions, vd, kt, vt, sigt, sigm, &
                            dxdy, psurf, emissions2d, dni, dnk)
   use chm_species_info_mod, only: nb_dyn_tracers
   integer(kind=4), intent   (in) :: dni, dnk
   real(kind=4),    intent(inout) :: conc       (dni, dnk, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: emissions  (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: vd         (dni, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: kt         (dni, dnk)
   real(kind=4),    intent   (in) :: vt         (dni, dnk)
   real(kind=4),    intent   (in) :: sigt       (dni, dnk)
   real(kind=4),    intent   (in) :: sigm       (dni, dnk)
   real(kind=4),    intent   (in) :: dxdy       (dni)
   real(kind=4),    intent   (in) :: psurf      (dni)
   real(kind=4),    intent   (in) :: emissions2d(dni, dnk, nb_dyn_tracers)
end subroutine mach_diff_rpnphy

subroutine mach_diff_prep(busper, busvol, chem_tr, metvar2d, metvar3d, &
                          metvar3dnocan, metvar3dcan, emisbio_can, tracers_can, &
                          kmod, kcan, imod, imod2, ni_can, ni_nocan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkc, nkt
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_pkg_gas_mod,     only: num_be_sp
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
end subroutine mach_diff_prep

subroutine mach_input_check(busvol, metvar2d, metvar3d, landuse)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use mach_drydep_mod,      only: lucprm
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni,SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni,chm_nk,SIZE_MV3D)
   real(kind=4),    intent   (in) :: landuse (chm_ni, lucprm)
end subroutine mach_input_check

subroutine mach_landuse(busper, metvar2d, landuse_out)
   use chm_ptopo_grid_mod,  only: chm_ni
   use chm_metvar_mod,      only: SIZE_MV2D
   use mach_drydep_mod,     only: lucprm
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: metvar2d   (chm_ni, SIZE_MV2D)
   real(kind=4),    intent  (out) :: landuse_out(chm_ni, lucprm)
end subroutine mach_landuse

subroutine mach_main(busper, busvol, chem_tr, metvar2d, metvar3d, &
                     slab_index, step, ni_can, ni_nocan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   integer(kind=4), intent   (in) :: slab_index
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: ni_can, ni_nocan
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr(chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
end subroutine mach_main

subroutine mach_output(busper, busvol, chem_tr, metvar2d, metvar3d, landuse)

   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_drydep_mod,      only: lucprm
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: landuse (chm_ni, lucprm)
end subroutine mach_output

subroutine mach_plumerise(busvol, chem_tr, metvar2d, metvar3d, landuse, &
                          slab_index)
   use chm_metvar_mod,          only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,      only: chm_ni, chm_nk
   use chm_species_info_mod,    only: nb_dyn_tracers
   use mach_drydep_mod,         only: lucprm
   integer(kind=4), intent   (in) :: slab_index
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: landuse (chm_ni, lucprm)
end subroutine mach_plumerise

subroutine mach_plumerise_weight(cur_source, z_magl, index_above_pbl, &
                                 safe_inv_mo_length, weight, pbl_hgt, &
                                 ustar, z_temp, wnd_spd, stb_func)
   use chm_ptopo_grid_mod,      only: chm_nk
   integer(kind=4), intent (in) :: cur_source
   real(kind=4),    intent (in) :: z_magl (chm_nk + 1)
   real(kind=4),    intent(out) :: weight (chm_nk)
   integer(kind=4), intent (in) :: index_above_pbl
   real(kind=4),    intent (in) :: pbl_hgt
   real(kind=4),    intent (in) :: ustar
   real(kind=4),    intent (in) :: safe_inv_mo_length
   real(kind=4),    intent (in) :: z_temp (chm_nk + 1)
   real(kind=4),    intent (in) :: wnd_spd(chm_nk + 1)
   real(kind=4),    intent (in) :: stb_func(chm_nk)
end subroutine mach_plumerise_weight

subroutine mach_plumerise_weight4fire(cur_source, z_magl, pbl_indx, pbl_hgt, &
                                      safe_inv_mo_length, weight, F_bio, rho)
   use chm_ptopo_grid_mod,      only: chm_nk
   integer(kind=4), intent (in) :: cur_source
   real(kind=4),    intent (in) :: z_magl(chm_nk+1)
   real(kind=4),    intent(out) :: weight(chm_nk)
   real(kind=4),    intent (in) :: pbl_hgt
   real(kind=4),    intent (in) :: F_bio (4)
   integer(kind=4), intent (in) :: pbl_indx
   real(kind=4),    intent (in) :: safe_inv_mo_length
   real(kind=4),    intent (in) :: rho   (chm_nk)
end subroutine mach_plumerise_weight4fire

subroutine mach_tridiag (a, b, c, r, u, dni, dnk)
   integer(kind=4), intent (in) :: dni, dnk
   real(kind=4),    intent (in) :: a(dni, dnk)
   real(kind=4),    intent (in) :: b(dni, dnk)
   real(kind=4),    intent (in) :: c(dni, dnk)
   real(kind=4),    intent (in) :: r(dni, dnk)
   real(kind=4),    intent(out) :: u(dni, dnk)
end subroutine mach_tridiag

subroutine mach_vit_diffusivity(busper, dxdy, zmom, kvkt, nmod, ni)
   use chm_ptopo_grid_mod,   only: chm_nk
   integer(kind=4), intent (in) :: ni
   integer(kind=4), intent (in) :: nmod(ni)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent (in) :: dxdy(ni)
   real(kind=4),    intent (in) :: zmom(ni, chm_nk)
   real(kind=4),    intent(out) :: kvkt(ni, chm_nk)
end subroutine mach_vit_diffusivity

subroutine mach_canopy_levels(busper, metvar3dcan, metvar3dnocan, metvar3d, &
                            metvar2d, imod, imod2, kmod, kcan, ni_can, ni_nocan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   integer(kind=4), intent   (in) :: ni_can, ni_nocan
   integer(kind=4), intent  (out) :: imod(ni_can)
   integer(kind=4), intent  (out) :: imod2(ni_nocan)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent  (out) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   real(kind=4),    intent  (out) :: metvar3dnocan(ni_nocan, chm_nk, SIZE_MV3D)
   integer(kind=4), intent  (out) :: kmod(ni_can, chm_nk)
   integer(kind=4), intent  (out) :: kcan(ni_can, nkc)

end subroutine mach_canopy_levels

subroutine mach_canopy_transfer(chem_tr, tracers_can, metvar2d, metvar3d, &
                                metvar3dcan, kmod, kcan, imod, ni_can, flag)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D 
   use chm_species_info_mod, only: nb_dyn_tracers
   integer(kind=4),   intent(in)    :: ni_can
   real(kind=4),      intent(inout) :: chem_tr(chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),      intent(inout) :: tracers_can(ni_can, nkc, nb_dyn_tracers)
   real(kind=4),      intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),      intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),      intent   (in) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   integer(kind=4),   intent   (in) :: imod(ni_can)
   integer(kind=4),   intent   (in) :: kmod(ni_can, chm_nk)
   integer(kind=4),   intent   (in) :: kcan(ni_can, nkc)
   integer(kind=4),   intent   (in) :: flag
end subroutine mach_canopy_transfer

subroutine mach_stepinit(busper, step, trnch, ni_can)
   integer(kind=4), intent   (in) :: step, trnch
   integer(kind=4), intent  (out) :: ni_can
   real(kind=4),    dimension(:), pointer, contiguous :: busper
end subroutine mach_stepinit

 subroutine mach_lai_adjust(busper, landuse, trnch)
   use chm_ptopo_grid_mod,   only: chm_ni
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: trnch
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: landuse(chm_ni, lucprm)
end subroutine mach_lai_adjust

subroutine mach_cffeps_main(chem_tr, metvar2d, metvar3d, slab_index, istep)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   integer(kind=4), intent   (in) :: slab_index, istep
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
end subroutine mach_cffeps_main

 subroutine DIFUVDFj_2D(TU, U, KU, GU, JNG, R, ALFA, BETA, S, SK, &
                        BND, TAU, type, F, NU, NR, N, NK)
   integer, intent(in)  :: NU, NR, N, NK
   real,    intent(out) :: TU(NU, NK)
   real,    intent(in)  :: U(NU, NK), KU(NR, NK), GU(NR, NK), R(NR,NK)
   real,    intent(in)  :: JNG(NR, NK), BND(NR, NK)
   real,    intent(in)  :: ALFA(N), BETA(N), S(n,NK), SK(n,NK), TAU, F
   integer, intent(in)  :: type
end subroutine difuvdfj_2d

 subroutine mach_cffeps_calc(feps, emissions, r_smoke, dz, met, mws, ross, cmcDj, &
                             dcmcUTC, ihour, istep)
   use mach_cffeps_mod,           only: feps_type, met_type, emissions_type, &
                                        max_timesteps
   type(feps_type),      intent(inout) :: feps
   type(emissions_type), intent(inout) :: emissions
   type(met_type),       intent(in)    :: met
   integer(kind=4),      intent(in)    :: cmcDj, istep, ihour
   real(kind=4),         intent(in)    :: dcmcUTC
   real(kind=4),         intent(out)   :: dz      ! plume height in metres
   real(kind=4),         intent(out)   :: r_smoke ! mixing ratio of smoke to clear air (g/kg)
   real(kind=4),         intent(inout), dimension(MAX_TIMESTEPS) :: mws
   real(kind=4),         intent(inout), dimension(MAX_TIMESTEPS) :: ross
end subroutine mach_cffeps_calc

 subroutine mach_cffeps_procalc(ifuel, fbp, t1, t2, a, p, d)
   use mach_cffeps_mod,  only: fbp_type, cffeps_fbp_accel
   integer(kind=4), intent   (in) :: ifuel
   type(fbp_type),  intent   (in) :: fbp
   real(kind=4),    intent(inout) :: t1, t2, a, p, d
end subroutine mach_cffeps_procalc

 subroutine mach_cffeps_energy(fbp, ifuel, qplume, tsurf, area, perimeter, &
                               growth, fmc, mcL, mcF, mcH)
   use mach_cffeps_mod,           only: fbp_type
   type(fbp_type),  intent(in)  :: fbp
   real(kind=4),    intent(out) :: qplume
   real(kind=4),    intent(in)  :: tsurf
   real(kind=4),    intent(in)  :: mcL, mcF, mcH
   real(kind=4),    intent(in)  :: area, perimeter, growth, fmc
   integer(kind=4), intent(in)  :: ifuel
end subroutine mach_cffeps_energy

 SUBROUTINE MACH_SUNCOS(SCOS,LMX,XLAT,XLON,HZ,DATE)
          INTEGER(kind=4), intent(in):: LMX
          REAL(kind=4), intent(out)  :: SCOS(LMX)
          REAL(kind=4), intent(in)   :: XLAT(LMX),XLON(LMX),HZ,DATE
end subroutine mach_suncos

end interface
end module mach_headers_mod
