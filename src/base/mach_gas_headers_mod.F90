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
! Fichier/File   : mach_gas_headers_mod.ftn90
! Creation       : H. Landry, Sept 2008
! Description    : Modules defining explicit interfaces for mach subroutines
!                  related to gas phase chemistry
!
! Extra info     :
!
!============================================================================

module mach_gas_headers_mod
   interface
!end trap head

subroutine mach_gas_drydep_main(busper, busvol, metvar2d, lfu, iseasn)
   use chm_ptopo_grid_mod,   only: chm_ni
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: iseasn  (chm_ni)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
end subroutine mach_gas_drydep_main

subroutine mach_gas_drydep_solver(vd, aero_resist, diff_resist, surf_resist, &
                                  iseason, lfu, lai_2d, metvar2d, vdg)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm, nb_gas_depo
   use chm_ptopo_grid_mod,   only: chm_ni
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent(out) :: vd         (nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: aero_resist(chm_ni, lucprm)
   real(kind=4),    intent(out) :: diff_resist(nb_gas_depo, chm_ni)
   real(kind=4),    intent(out) :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent (in) :: lai_2d     (chm_ni)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
   real(kind=4), optional, intent(out) :: vdg (lucprm, chm_ni)
end subroutine mach_gas_drydep_solver

subroutine mach_gas_drydep_solver2(vd, aero_resist, diff_resist, surf_resist, &
                                   iseason, lfu, metvar2d, vdg)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm, nb_gas_depo
   use chm_ptopo_grid_mod,   only: chm_ni
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent(out) :: vd         (nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: aero_resist(chm_ni, lucprm)
   real(kind=4),    intent(out) :: diff_resist(nb_gas_depo, chm_ni)
   real(kind=4),    intent(out) :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
   real(kind=4), optional, intent(out) :: vdg (lucprm, chm_ni)
end subroutine mach_gas_drydep_solver2

subroutine mach_gas_bidi(vdg, metvar2d, busper, busvol)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm
   use chm_ptopo_grid_mod,   only: chm_ni
   real(kind=4),    intent   (in) :: vdg      (lucprm, chm_ni)
   real(kind=4),    intent   (in) :: metvar2d (chm_ni, SIZE_MV2D)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
end subroutine mach_gas_bidi

subroutine mach_gas_drydep_stat(vd, aero_resist, diff_resist, surf_resist, &
                                lfu, metvar2d)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm, nb_gas_depo
   use chm_ptopo_grid_mod,   only: chm_ni
   real(kind=4),    intent(in) :: vd         (nb_gas_depo, chm_ni)
   real(kind=4),    intent(in) :: aero_resist(chm_ni, lucprm)
   real(kind=4),    intent(in) :: diff_resist(nb_gas_depo, chm_ni)
   real(kind=4),    intent(in) :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   real(kind=4),    intent(in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent(in) :: metvar2d   (chm_ni, SIZE_MV2D)
end subroutine mach_gas_drydep_stat

subroutine mach_gas_drydep_ra(aero_resist, iseason, lfu, metvar2d)
   use mach_drydep_mod,    only: lucprm
   use chm_metvar_mod,     only: SIZE_MV2D
   use chm_ptopo_grid_mod, only: chm_ni
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
!  Land use and land use-dependant aerodynamic resistance - implicit conversion from the 1D vector
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent(out) :: aero_resist(chm_ni, lucprm)
end subroutine mach_gas_drydep_ra

subroutine mach_gas_drydep_ra2(aero_resist, iseason, lfu, metvar2d)
   use chm_metvar_mod,     only: SIZE_MV2D
   use mach_drydep_mod,    only: lucprm
   use chm_ptopo_grid_mod, only: chm_ni
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
!  Land use and land use-dependant aerodynamic resistance - implicit conversion from the 1D vector
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent(out) :: aero_resist(chm_ni, lucprm)
end subroutine mach_gas_drydep_ra2

subroutine mach_gas_main(busper, busvol, chem_tr, metvar2d, metvar3dnocan, &
                         step, ni_nocan, imod2, lfu)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: ni_nocan
   integer(kind=4), intent   (in) :: imod2   (ni_nocan)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3dnocan(ni_nocan, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
end subroutine mach_gas_main

subroutine mach_gas_canopy(busper, busvol, chem_tr, tracers_can, metvar2d, &
                           metvar3dcan, step, ni_can, imod, kcan, kmod, lfu)
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: step
   integer(kind=4), intent   (in) :: ni_can
   integer(kind=4), intent   (in) :: imod(ni_can)
   integer(kind=4), intent   (in) :: kcan(ni_can, nkc), kmod(ni_can, chm_nk)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),    intent(inout) :: tracers_can(ni_can, nkc, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
end subroutine mach_gas_canopy

subroutine mach_gas_strato(busper, hu_ppm, ztplus, sigt, o3, zpplus, lino3_new, &
                           gni, gnk, nmod)
   integer(kind=4), intent (in) :: gni, gnk
   integer(kind=4), intent (in) :: nmod     (gni)
   real   (kind=4), dimension(:), pointer, contiguous :: busper
   real   (kind=4), intent (in) :: zpplus   (gni)
   real   (kind=4), intent (in) :: o3       (gni, gnk)  !ug/kg
   real   (kind=4), intent (in) :: hu_ppm   (gni, gnk)  !ppmv
   real   (kind=4), intent (in) :: ztplus   (gni, gnk)  !K
   real   (kind=4), intent (in) :: sigt     (gni, gnk)  !
   real   (kind=4), intent(out) :: lino3_new(gni, gnk)  !ug /kg
end subroutine mach_gas_strato

subroutine mach_gas_messy(gmetvar3d, metvar2d, o3vmr, rj, nmod, gni, gnk)
   use chm_ptopo_grid_mod,  only: chm_ni
   use chm_metvar_mod,      only: SIZE_MV2D, SIZE_MV3D
   use mach_pkg_gas_mod,    only: njrxs
!  IO variable
   integer(kind=4),                                 intent(in)  :: gni, gnk
   integer(kind=4), dimension(gni),                 intent(in)  :: nmod
   real(kind=4),    dimension(chm_ni,   SIZE_MV2D), intent(in)  :: metvar2d
   real(kind=4),    dimension(gni, gnk, SIZE_MV3D), intent(in)  :: gmetvar3d
   real(kind=4),    dimension(gni, gnk)           , intent(in)  :: o3vmr
   real(kind=8),    dimension(gni* gnk, njrxs)    , intent(out) :: rj
end subroutine mach_gas_messy

subroutine mach_adom2_drive(yg, indx, zen, ig, jg, rgs, bgs)
   use mach_pkg_adom2_mod, only: nspec, nreac, nprcf
   integer(kind=4), intent(inout) :: indx
   real(kind=4),    intent   (in) :: zen
   integer(kind=4), intent   (in) :: ig
   integer(kind=4), intent   (in) :: jg
   real(kind=4),    intent(inout) :: yg  (nspec)
   real(kind=4),    intent   (in) :: rgs (nreac)
   real(kind=4),    intent   (in) :: bgs (nprcf)
end subroutine mach_adom2_drive

subroutine mach_adom2yb_main(p2d, tplus, hu_vmr, sigt, rjval, tppmgs, &
                             rjval_lookup, csza, voc_diff, step, gni, gnk)
   use mach_pkg_gas_mod,     only: nspec, nsp_soa_gases, njrxs
   use mach_adom2_rates_mod, only: ntype2, ntype4
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
end subroutine mach_adom2yb_main

subroutine mach_adom2kpp_main(p2d, tplus, hu_vmr, sigt, rjval, tppmgs, hstart, &
                              voc_diff, gni, gnk)

   use mach_pkg_gas_mod,     only: nspec, nsp_soa_gases, njrxs
   integer(kind=4), intent   (in) :: gni, gnk
   real(kind=4),    intent   (in) :: p2d     (gni, gnk)
   real(kind=4),    intent   (in) :: tplus   (gni, gnk)
   real(kind=4),    intent   (in) :: hu_vmr  (gni, gnk)
   real(kind=4),    intent   (in) :: sigt    (gni, gnk)
   real(kind=8),    intent   (in) :: rjval   (gni * gnk, njrxs)
   real(kind=4),    intent(inout) :: tppmgs  (nspec, gni, gnk)
   real(kind=8),    intent(inout) :: hstart  (gni, gnk)
   real(kind=8),    intent  (out) :: voc_diff(nsp_soa_gases, gni, gnk)
end subroutine mach_adom2kpp_main

subroutine mach_adom2kppb_main(p2d, tplus, hu_vmr, sigt, rjval, tppmgs, hstart, &
                               rjval_lookup, voc_diff, gni, gnk)
   use mach_pkg_gas_mod,     only: nspec, nsp_soa_gases, njrxs
   use mach_adom2_rates_mod, only: ntype2, ntype4
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
end subroutine mach_adom2kppb_main

SUBROUTINE mach_kppb_adom2day_rodas3(N, Y, nfix, fix, nreact, rconst, nprcf, bgs, &
                     delt, H, RelTol, AbsTol, IERR)
   integer(kind=4), intent(in) :: N, nfix, nreact, nprcf
   REAL(kind=8), INTENT(IN) ::  fix(nfix)
   REAL(kind=4), INTENT(IN) ::  delt
   REAL(kind=8), INTENT(IN) ::  rconst(nreact), bgs(nprcf)
   REAL(kind=8), INTENT(IN) ::  AbsTol(N), RelTol(N)
   INTEGER(kind=4), INTENT(OUT) :: IERR
   REAL(kind=8), INTENT(INOUT) :: H, Y(N)
end subroutine mach_kppb_adom2day_rodas3

SUBROUTINE mach_kppb_adom2night_rodas3 (N, Y, nfix, fix, nreact, rconst, nprcf, bgs, &
                     delt, H, RelTol, AbsTol, IERR)
  integer(kind=4), intent(in) :: N, nfix, nreact, nprcf
   REAL(kind=8), INTENT(IN) ::  fix(nfix)
   REAL(kind=4), INTENT(IN) ::  delt
   REAL(kind=8), INTENT(IN) ::  rconst(nreact), bgs(nprcf)
   REAL(kind=8), INTENT(IN) ::  AbsTol(N), RelTol(N)
   INTEGER(kind=4), INTENT(OUT) :: IERR
   REAL(kind=8), INTENT(INOUT) :: H, Y(N)
end subroutine mach_kppb_adom2night_rodas3

subroutine mach_kpp_adom2_rates(rct, temp2d, p2d, cfactor2d, cno2d, gni, gnk)
  use mach_pkg_adom2_mod,   only: nreact, nprcf, nreact_kpp
  integer(kind=4), intent(in) :: gni, gnk
  real(kind=8),    intent(out) :: rct(gni * gnk, nreact)
  real(kind=4),    intent(in)  :: temp2d(gni, gnk) ! Deg. Kelvin
  real(kind=4),    intent(in)  :: cfactor2d(gni, gnk) ! conversion from ppm to molecules/cm3
  real(kind=4),    intent(in)  :: cno2d(gni, gnk)  ! NO conc in molec/cm3
  real(kind=4),    intent(in)  :: p2d(gni, gnk)    ! Pressure profile in Pa
end subroutine mach_kpp_adom2_rates

 subroutine mach_kpp_interface(conc, rxt1d, atol, rtol, shstart, itolctr)
   use mach_pkg_gas_mod,    only: nspec, nreact, nvar
   real(kind=8),    intent(inout) :: conc (nspec) ! concentration for all species
   real(kind=8),    intent(in)    :: rxt1d(nreact) ! reaction rate
   real(kind=8),    intent(in)    :: atol (nvar)
   real(kind=8),    intent(in)    :: rtol (nvar)
   real(kind=8),    intent(inout) :: shstart
   integer(kind=4), intent(in)    :: itolctr
end subroutine mach_kpp_interface

subroutine mach_kppb_adom2_interface(conc, rgs, bgs, atol_in, rtol_in, hs, &
                                     day_night_flag)
  use mach_pkg_adom2_mod, only: nreac, nprcf, nspec
  real(kind=8),    intent(inout) :: hs
  real(kind=8),    intent(in)    :: atol_in, rtol_in
  real(kind=8),    intent(in)    :: rgs(nreac), bgs(nprcf)
  real(kind=8),    intent(inout) :: conc(nspec)
  integer(kind=4), intent(in)    :: day_night_flag
end subroutine mach_kppb_adom2_interface

subroutine mach_kpp_integrator(N,Y,Tstart,Tend, AbsTol,RelTol, FIX, RCONST, &
                               RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
  use chm_utils_mod,    only: dp
  use mach_pkg_gas_mod, only: nfix, nreact
   INTEGER,       INTENT(IN)    :: N
   REAL(kind=dp), INTENT(INOUT) :: Y(N)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCONST(NREACT)  ! JC add for OMP>1
   REAL(kind=dp), INTENT(IN)    :: FIX(NFIX)       ! JC add for OMP>1
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
end subroutine mach_kpp_integrator

subroutine mach_marine_halo(ozone, pressure, lfu, imod, ni, nk)
   use chm_ptopo_grid_mod,   only: chm_ni
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: ni, nk
   integer(kind=4), intent   (in) :: imod(ni)
   real(kind=4),    intent(inout) :: ozone(ni, nk)
   real(kind=4),    intent   (in) :: lfu(chm_ni, lucprm)
   real(kind=4),    intent   (in) :: pressure(ni, nk)
end subroutine mach_marine_halo

subroutine mach_kpp_saprc07_rates(rconst_out, ksoa, kNO_r2o2, kHO2_r2o2, tp, &
                                  conv, gni, gnk)
  use mach_pkg_gas_mod,   only: nreact, nvsoa
  integer(kind=4), intent(in) :: gni, gnk
  real(kind=8), intent(out)   :: rconst_out(gni*gnk, nreact)
  real(kind=8), intent(out)   :: ksoa      (3, nvsoa, gni, gnk)
  real(kind=8), intent(out)   :: kNO_r2o2  (gni, gnk)
  real(kind=8), intent(out)   :: kHO2_r2o2 (gni, gnk)
  real(kind=4), intent(in )   :: tp        (gni, gnk) ! Deg. Kelvin
  real(kind=8), intent(in )   :: conv      (gni*gnk)  ! concentration of air in molecules/cm3
end subroutine mach_kpp_saprc07_rates

subroutine mach_soa_jiang(moi, dsoa, soa_gas_diff, tp, p2d, gni, gnk)
   use mach_pkg_gas_mod,   only: nsp_soa_gases
   integer(kind=4), intent (in) :: gni, gnk
   real(kind=8),    intent (in) :: soa_gas_diff(nsp_soa_gases, gni, gnk)
   real(kind=8),    intent (in) :: moi         (gni, gnk)
   real(kind=4),    intent(out) :: dsoa        (gni, gnk)
   real(kind=4),    intent (in) :: p2d         (gni, gnk)
   real(kind=4),    intent (in) :: tp          (gni, gnk)
end subroutine mach_soa_jiang

subroutine mach_soa_odum(moi, dsoa, dvoc, ksoa, kNO_r2o2, kHO2_r2o2, tplus, &
                         conx1, trppm, gni, gnk)
   use chm_utils_mod,      only: dp, sp, i4
   use mach_pkg_gas_mod,   only: nsp_soa_gases, nspec, nvsoa, minconc
   integer(i4), intent (in) :: gni, gnk
   real(dp),    intent (in) :: dvoc (nsp_soa_gases, gni, gnk)
   real(dp),    intent (in) :: moi  (gni, gnk)
   real(dp),    intent (in) :: ksoa (3, nvsoa, gni, gnk)
   real(dp),    intent (in) :: kHO2_r2o2(gni, gnk)
   real(dp),    intent (in) :: kNO_r2o2 (gni, gnk)
   real(sp),    intent (in) :: tplus(gni, gnk)
   real(dp),    intent (in) :: conx1(gni * gnk)
   real(sp),    intent (in) :: trppm(nspec, gni, gnk)
   real(sp),    intent(out) :: dsoa (gni, gnk)
end subroutine mach_soa_odum

subroutine mach_saprc_main(p2d, tplus, hu_ppm, rjval, tppmgs, hstart, &
                           moi, dsoa, gni, gnk)
   use mach_pkg_gas_mod,     only: nspec, njrxs
   integer(kind=4), intent(in) :: gni, gnk
   real(kind=4), intent   (in) :: p2d   (gni, gnk)
   real(kind=4), intent   (in) :: tplus (gni, gnk)
   real(kind=4), intent   (in) :: hu_ppm(gni, gnk)
   real(kind=8), intent   (in) :: rjval (gni * gnk, njrxs)
   real(kind=8), intent   (in) :: moi   (gni, gnk)
   real(kind=4), intent(inout) :: tppmgs(nspec, gni, gnk)
   real(kind=8), intent(inout) :: hstart(gni, gnk)
   real(kind=4), intent  (out) :: dsoa  (gni, gnk)
end subroutine mach_saprc_main

end interface
end module mach_gas_headers_mod
