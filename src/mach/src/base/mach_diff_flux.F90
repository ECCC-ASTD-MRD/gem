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
! Fichier/File   : mach_diff_flux.ftn90
! Creation       : P. Makar, W. Gong, S. Menard, May 2010. 
! Description    : Vertical diffusion for chemical species
!
! Extra info     : The method is described in "Vertical_diffusion_in_CHRONOS_AURAMS_and_GEM-MACH.doc"
!                  The specific variation calculated here the use of a 1/2 layer thick shadow zone
!                  above and below the region in order to define the derivatives at the model top
!                  and model bottom.  i.e. similar to the method used in CHRONOS.
!                  In order to reduce processing time, portions of the matrix system that are
!                  species invariant are calculated first, and the species dependant terms are
!                  added in where appropriate in the processing.
!
! Arguments:
!           IN
!              emissions  -> The sum of area and biogenic emissions
!              vd         -> Dry deposition velocities for gas
!              ubf        -> Species upper lid flux boundary condition index values
!              kt         -> Thermal vertical diffusion coefficient
!              dxdy       -> Area of grid square
!              rho        -> Air density
!              zf         -> Concentration-layer height above sea-level
!              zh         -> Diffusion constant height above sea-level
!
!           INOUT
!              CONC       ->   species concentrations
!
!=============================================================================
!
!!if_on
subroutine mach_diff_flux(conc, emissions, vd, ubf, kt, rho, zf, zh, dxdy, &
                          dni, dnk)
   use chm_species_info_mod, only: nb_dyn_tracers
!!if_off
   use chm_utils_mod,        only: chm_timestep, chm_lun_out, global_debug
   use mach_headers_mod,     only: mach_tridiag
   implicit none
!!if_on
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
!!if_off
!
!  local workspace arrays and variables:
!
   real(kind=4),   dimension(dni)          :: s_top, s_bot, s_dep
   real(kind=4),   dimension(dni, dnk)  :: a, b, c
   real(kind=4),   dimension(dni, dnk)  :: beta, tconc
   integer(kind=4)                            :: i, k, sp
   real(kind=4)                               :: b1
   logical(kind=4)                            :: local_dbg
!
!  diagnostic arrays
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   if (local_dbg) then
      write (chm_lun_out, *) "Entering mach_diff_flux"
   end if


!  First, calculate the coefficients for the desired matrix-vector system, for the
!  current time-step:
!  Interior rows of tridiagonal matrix:
   do k = 2, dnk - 1
      do i = 1, dni
         b1 = (zh(i, dnk - k) - zh(i, dnk - k + 1))
         c(i, k) = - chm_timestep * kt(i, dnk - k) /     &
                     ((zf(i, dnk - k) - zf(i, dnk - k + 1)) * b1)
         a(i, k) = - chm_timestep * kt(i, dnk - k + 1) / &
                     ((zf(i, dnk - k + 1) - zf(i, dnk - k + 2)) * b1)
         b(i, k) = 1. - c(i, k) - a(i, k)
      end do
   end do

   do i = 1, dni

!  Top and bottom (boundary condition rows)
!  bottom layer equations:
      b1      = (zf(i, dnk - 1) - zf(i, dnk)) **2
      c(i, 1) = - chm_timestep * kt(i, dnk - 1) / b1
      b(i, 1) = 1. - c(i, 1)

!  top layer equations:
      b1       = (zf(i, 1) - zf(i, 2)) **2
      a(i, dnk) = - chm_timestep * kt(i, 1) / b1
      b(i, dnk) =  1. - a(i, dnk)
      
      s_dep(i) = chm_timestep / (zf(i, dnk - 1) - zf(i, dnk))
      s_bot(i) = chm_timestep / &
             (rho(i, dnk) * dxdy(i) * (zf(i, dnk - 1) - zf(i, dnk)))
      s_top(i) = chm_timestep / (rho(i, 1) * dxdy(i) * (zf(i, 1) - zf(i, 2)))
   end do

!  Loop through all species:
   do sp = 1, nb_dyn_tracers

!  Prepare rhs of equation
      do k = 1, dnk
         do i = 1, dni
            beta(i, k) = conc(i, dnk - k + 1, sp)
         end do
      end do

      do i = 1, dni
!!!  Deposition and emissions as flux terms:  commented out in this version:
!  add deposition term:
         b(i, 1)    =  (1. - c(i, 1)) + 2. * vd(i, sp) * s_dep(i)
!  add emission term on rhs
         beta(i, 1) =  beta(i, 1) +  emissions(i, sp) * s_bot(i)
      end do

! upper lid source term added on rhs:
      do i = 1, dni
         beta(i, dnk) =  beta(i, dnk) - ubf(i, sp) * s_top(i)
      end do

      call mach_tridiag(a, b, c, beta, tconc, dni, dnk)

      do k = 1, dnk
         do i = 1, dni
            conc(i, k, sp) = tconc(i, dnk - k + 1)
         end do
      end do
!      
   end do

!  conc now holds the updated species concentrations
   if (local_dbg) then
      write (chm_lun_out, *) "Leaving mach_diff_flux"
   end if

   return

end subroutine mach_diff_flux
