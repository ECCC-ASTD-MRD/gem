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
! Fichier/File   : mach_diff_rpnphy.ftn90
! Creation       : Sylvie Gravel, May 2012
! Description    : Performs vertical diffusion using operator difuvdfj from physics library
!                  Adapted to staggered model
!
! Arguments:
!           IN
!              emissions  -> The sum of area and biogenic emissions
!              vd         -> Dry deposition velocities for gas
!              kt         -> Thermal vertical diffusion coefficient
!              dxdy       -> Area of grid square
!              psurf      -> Surface pressure
!              vt         -> Virtual temperature profile
!              sigt       -> Sigma coordinate on thermodynamic levels
!              sigm       -> Sigma coordinate on momentum levels
!
!           INOUT
!              conc       ->   species concentrations
!
!=============================================================================
!
!!if_on
subroutine mach_diff_rpnphy(conc, emissions, vd, kt, vt, sigt, sigm, &
                            dxdy, psurf, emissions2d, dni, dnk)
   use chm_species_info_mod, only: nb_dyn_tracers
!!if_off
   use chm_utils_mod,        only: chm_timestep, chm_lun_out, global_debug
   use chm_consphychm_mod,   only: grav, rgasd
   use chm_nml_mod,          only: chm_vert_diff_s, chm_ae_spread_l
   use phy_options,          only: tlift
   use mach_headers_mod,     only: difuvdfj_2d
   implicit none
!!if_on
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
!!if_off
!
!  local workspace arrays and variables:
!
   integer(kind=4)  :: i, k, ix, sp, itype, ni_net
   real(kind=4)     :: gam0, rsg
   real(kind=4), dimension(dni) :: gam1
   real(kind=4), dimension(dni * nb_dyn_tracers) :: alfa, beta
   real(kind=4), dimension(dni * nb_dyn_tracers, dnk) :: emis2D, ktsg, &
                                                         sigm_sg, sigt_sg
   real(kind=4), dimension(dni * nb_dyn_tracers, dnk) :: tconc, zero,conc_2d
   real(kind=4), dimension(dnk) :: e_weight
   real(kind=4)     :: denom, weight
   logical(kind=4)  :: local_dbg

!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   if (local_dbg) then
      write (chm_lun_out, *) "Entering mach_diff_rpnphy"
   end if

   e_weight = 0.0
   if (chm_ae_spread_l) then
     e_weight(dnk) = 0.5
     e_weight(dnk-1) = 0.5
   else
     e_weight(dnk) = 1.0
   end if
!
!  If emissions are provided as g/s, the weight factor is (g2ug*grav/(psurf*Surface*(sigm(k)-sigm(k-1))).
!  If emissions are provided as ug/s, the weight factor is (grav/(psurf*Surface*(sigm(k)-sigm(k-1)))).
!  Prepare weight factors to modify concentrations with 50% of emissions for
!  one layer above lowest prognostic level
!  remaining 50% will be added as part of surface flux
   if (chm_ae_spread_l) then
      k = dnk-1
      do i = 1, dni
         denom= psurf(i) * dxdy(i) * (sigm(i, k+1) - sigm(i, k))
         weight = grav * chm_timestep / denom
         do sp = 1, nb_dyn_tracers
            conc(i, k, sp) = conc(i, k, sp) + weight * emissions(i, sp) * &
                             e_weight(k)
         end do
      end do
   end if

   zero = 0.

   ni_net = dni * nb_dyn_tracers
!  Normalization factors for vertical diffusion in sigma coordinate factors are applied to diffusion coefficients
!  and to surface fluxes
   do i = 1, dni
      gam1(i) = grav / (psurf(i) * dxdy(i))
   end do

   rsg = (grav / rgasd)
   do k = 1, dnk
      ix = 0
      do i = 1, dni
         gam0  = rsg * sigt(i, k) / vt(i, k)
         do sp = 1, nb_dyn_tracers
            ix = ix + 1
            ktsg(ix, k) = kt(i, k) * gam0**2
            sigm_sg(ix, k) = sigm(i, k)
            sigt_sg(ix, k) = sigt(i, k)

            conc_2d(ix, k) = conc(i, k, sp)
            emis2D(ix, k) = gam1(i) * emissions2d(i, k, sp)
         end do
      end do
   end do

   if (tlift == 1) then
      itype = 6
   else
      itype = 5
   end if

   ix = 0
   if (chm_vert_diff_s == 'RPNPHY_I') then

      do i = 1, dni
         gam0 = -rsg * sigt(i, dnk) / vt(i, dnk)
         do sp = 1, nb_dyn_tracers
            ix = ix + 1
            alfa(ix) = gam1(i) * e_weight(dnk) * emissions(i, sp)
            beta(ix) = gam0 * vd(i, sp)
         end do
      end do

   else
      do i = 1, dni
         gam0 = -rsg * sigt(i, dnk) / vt(i, dnk)
         do sp = 1, nb_dyn_tracers
            ix = ix + 1
            alfa(ix) = gam1(i) * e_weight(dnk) * emissions(i, sp) + &
                       gam0 * vd(i, sp) * conc(i, dnk, sp)
         end do
      end do
      beta = 0.
   end if

!  Calculate tendancies due to vertical diffusion. Store in tconc.

   call difuvdfj_2D(tconc, conc_2d, ktsg, zero, zero, zero, alfa, beta,  &
                     sigm_sg, sigt_sg, emis2D, chm_timestep, itype, 1., &
                     ni_net, ni_net, ni_net, dnk)

   do k = 1, dnk
      ix = 0
      do i = 1, dni
         do sp = 1, nb_dyn_tracers
            ix = ix + 1
            conc(i, k, sp)  = conc(i, k, sp) + chm_timestep * tconc(ix, k)
         end do
      end do
   end do

   return

 end subroutine mach_diff_rpnphy
