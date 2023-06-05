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
! Fichier/File   : mach_diffusion.ftn90
! Creation       : Paul Makar, Sylvie Gravel, Nov 2007
! Description    : Master subroutine for GEM-MACH vertical diffusion
!                  Finite differences are used to solve the diffusion equation, one-step backward
!                  difference in time (i.e. fully implicit), and 1st to second order in space,
!                  depending on the method chosen.  The interior rows (k = 2..NK-1) of the
!                  matrix system use second order finite differences, resulting in a tridiagonal
!                  system of equations.  There are 4 approached implemented, triggered by bc_type
!                  (see bc_type's description below)
!
! Extra info     : The method is described in detail in
!                  http://ulysse.cmc.ec.gc.ca/svn/gem-mach/doc/Vertical_diffusion_in_CHRONOS_AURAMS_and_GEM-MACH.doc
!
! Arguments:
!           IN
!              DXDY   -> Area of grid square
!              PSURF  -> Surface pressure
!              RHO    -> Air density
!              KT     -> Thermal vertical diffusion coefficient
!              VT     -> Virtual temperature profile
!              ZPLUS  -> Concentration layer height above sea-level
!              ZMOM   -> Diffusion constant height above sea level
!              busper -> Permanent bus
!              VD     -> Deposition velocities
!              emisbio-> 2-D biogenic emissions within forest canopy
!
!           IN/OUT
!              busvol -> Volatile bus
!              CONC   -> Chemical species' concentrations
!
!=============================================================================
!
!!if_on
subroutine mach_diffusion(busper, busvol, conc, vd, psurf, dxdy, rho, kt, vt, &
                          sigm, sigt, zplus, zmom, echoice, nmod, dni, dnk, &
                          emisbio2D)
   use chm_species_info_mod, only: nb_dyn_tracers
   use mach_pkg_gas_mod,     only: num_be_sp
!!if_off
   use chm_utils_mod,        only: chm_timestep, chm_lun_out, global_debug
   use chm_nml_mod,          only: chm_vert_diff_s, chm_met_modulation_s, &
                                   chm_vit_l, chm_debug_2d_i, dbg_itr
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_FMET, dbg_2D
   use chm_consphychm_mod,   only: grav
   use mach_headers_mod,     only: mach_diff_flux, mach_diff_boundary, mach_diff_rpnphy
   use mach_pkg_gas_mod,     only: be_species

   implicit none
!
!  Declaration of Subroutine arguments
!
!!if_on
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
!!if_off
!
!  Declaration of local variables
!
   integer(kind=4)             :: spx
   integer(kind=4)             :: i, k, sp
   real(kind=4)                :: sig_m(dnk + 1)
   real(kind=4)                :: mass_before, mass_after, mass_emiss, &
                                  mass_fix, mass_sum1, mass_sum2, p_scale
   real(kind=4), parameter     :: g2ug = 1.0e6
   real(kind=4), parameter     :: min_value = 1.0e-12
   real(kind=4)                :: emissions2d(dni, dnk, nb_dyn_tracers)
   real(kind=4)                :: emissions(dni, nb_dyn_tracers)
   real(kind=4)                :: ubf(dni, nb_dyn_tracers)
   real(kind=4)                :: fmet(dni), fmet_sp(dni)
   logical(kind=4)             :: local_dbg, local_dbg_2d
!
! -----------------------------------------------------------------------------------------------------------------
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
   local_dbg_2d = ((.false. .or. global_debug) .and. (chm_debug_2d_i > 3))
!
!  Mass balance diagnostics
!  build mass before in dbg_2d(1)
   if (local_dbg_2d .and. (dbg_itr > 0)) then
      mass_before = 0.0
      spx = dbg_itr
      sig_m(dnk + 1) = 1.0
      do i = 1, dni
         mass_sum1 = 0.0
         do k = 1, dnk
            sig_m(k) = sigm(i, k)
            p_scale = (sig_m(k + 1) - sig_m(k)) * psurf(i) * dxdy(i) / .980616e+1
            mass_sum1 = mass_sum1 + p_scale * conc(i, k, spx)
         end do
         busvol(sm(dbg_2d(1)) % out_offset + nmod(i) - 1) = mass_sum1
         mass_before = mass_before + mass_sum1
      end do
   end if
!
!  Construct arrays of emission fluxes; original units are g/s, convert to ug/s.
!
   emissions2d = 0.0
   emissions   = 0.0
   ubf         = 0.0
!
!  Values of echoice:
!  0:   Zero emissions solution with original diffusivities.
!  1:   Original solution:  all emissions (anthropogenic + biogenic) are in a surface flux
!       array which will be used as a diffusion boundary conditions.
!  2:   On-road emissions as flux boundary conditions to diffusion,
!       with additional vehicle induced turbulent diffusivities
   if (echoice == 1 .or. (.not. chm_vit_l)) then
!
!     Meteorological modulation of aerosol emissions
      if (chm_met_modulation_s /= 'NIL') then
         do i = 1, dni
            fmet(i) = busvol(sm(sp_FMET) % out_offset + nmod(i) - 1)
         end do
      end if

      do sp = 1, nb_dyn_tracers
!
! Meteorological modulation of crustal material emissions
         fmet_sp = 1.0
         if ((chm_met_modulation_s == 'CM_ONLY') .and. &
             (sm(sp) % ae_name(1:3) == 'ECM')) then
             fmet_sp = fmet
         end if
!
!  Area emissions added to emissions array:
!  original units are g/s:  convert to ug/s
         if (sm(sp) % ae_offset > 0) then
            do i = 1, dni
               emissions(i, sp) = emissions(i, sp) + &
                    (busper(sm(sp) % ae_offset + nmod(i) - 1) * g2ug * fmet_sp(i))
            end do
         end if

!  Fugitive area emissions added to emissions array, modulated by fmet:
!  original units are g/s:  convert to ug/s
         if (sm(sp) % fae_offset > 0) then
            do i = 1, dni
               emissions(i, sp) = emissions(i, sp) + &
                    (busper(sm(sp) % fae_offset + nmod(i) - 1) * g2ug * fmet(i))
            end do
         end if
!
!  Biogenic emissions added to emissions array:
!  original units are g/s:  convert to ug/s
         if (sm(sp) % be_offset > 0) then
            do i = 1, dni
               emissions(i, sp) = emissions(i, sp) + &
                                  (busper(sm(sp) % be_offset + nmod(i) - 1) * g2ug)
            end do
         end if
!
!  Bidirectional flux emissions added to emissions array:
!  original units are g/s:  convert to ug/s
         if (sm(sp) % bd_offset > 0) then
            do i = 1, dni
               emissions(i, sp) = emissions(i, sp) + &
                                  (busvol(sm(sp) % bd_offset + nmod(i) - 1) * g2ug)
            end do
         end if
      end do

   end if

   if (echoice > 1) then
      do sp = 1, nb_dyn_tracers
         if (sm(sp) % mae_offset > 0) then
            do i = 1, dni
               emissions(i, sp) = emissions(i, sp) + &
                      busper(sm(sp) % mae_offset + nmod(i) - 1) * g2ug
            end do
         end if
      end do

   end if
!
!  Biogenics as fluxes into model layers:
   if (present(emisbio2D)) then
      do sp = 1, num_be_sp
         do k = 1, dnk
            do i = 1, dni
               emissions2d(i, k, be_species(sp)) = emisbio2D(i, k, sp) * g2ug
            end do
         end do
      end do
   end if
!
!  Calculate new concentrations resulting from diffusion, all species in conc:
!
!  Each of the four bc_type options results in a different subroutine being
!  called; the subroutines are similar in nature in that the interior rows
!  of the matrix that is constructed within each are identical; differences
!  arise because of the boundary condition finite differences options.
!  Each subroutine in turn calls a tridiagonal solver, which solves the system
!  of equations and returns the new set of concentrations to the next level up
!  in the program.
!    chm_vert_diff_s  -> BOUNDARY
!                        FLUX
!                        RPNPHY : using the same approach as physics library
   select case (chm_vert_diff_s)
      case ('FLUX')
         call mach_diff_flux(conc, emissions, vd, ubf, kt, rho, zplus, zmom, &
                             dxdy, dni, dnk)
      case ('BOUNDARY')
         call mach_diff_boundary(conc, emissions, vd, ubf, kt, rho, &
                                 zplus, zmom, dxdy, dni, dnk)
      case ('RPNPHY', 'RPNPHY_I', 'RPNPHY_U')
         call mach_diff_rpnphy(conc, emissions, vd, kt, vt, sigt, sigm, &
                               dxdy, psurf, emissions2d, dni, dnk)
   end select
!
!  Mass balance diagnostics
!  build mass after in dbg_2d(2)
!  build mass after with min check in dbg_2d(3)
!  build mass of emissions in dbg_2d(4)
   if (local_dbg_2d .and. (dbg_itr > 0)) then
      mass_after = 0.
      mass_fix = 0.
      mass_emiss = 0.
      sig_m(dnk) = 1.0
      do i = 1, dni
         mass_sum1 = 0.0
         mass_sum2 = 0.0
         do k = 1, dnk
            sig_m(k) = sigm(i, k)
            p_scale = (sig_m(k + 1) - sig_m(k)) * psurf(i) * dxdy(i) / .980616e+1
            mass_sum1 = mass_sum1 + p_scale * conc(i, k, spx)
            mass_sum2 = mass_sum2 + p_scale * max(conc(i, k, spx), min_value)
         end do
         busvol(sm(dbg_2d(2)) % out_offset + nmod(i) - 1) = mass_sum1
         busvol(sm(dbg_2d(3)) % out_offset + nmod(i) - 1) = mass_sum2
         busvol(sm(dbg_2d(4)) % out_offset + nmod(i) - 1) = chm_timestep * emissions(i, spx)
         mass_emiss = mass_emiss + chm_timestep * emissions(i, spx)
         mass_after = mass_after + mass_sum1
         mass_fix   = mass_fix + mass_sum2
      end do
      if (local_dbg) write(chm_lun_out,1000) sm(spx) % dyn_name, mass_before, mass_emiss, mass_before+mass_emiss, mass_after, &
        mass_after-mass_before-mass_emiss, mass_fix
   end if
 1000    format('In mach_diffusion, mass changes for species ',a,' : ',/, &
         'mass before | mass of emissions | mass before + emissions | mass after | difference  | mass fixed ', /,&
             1x,e10.3,          10x,e10.3,          10x,e10.3,          9x,e10.3,    3x,e10.3,     4x,e10.3)
!

   return
!
end subroutine mach_diffusion
