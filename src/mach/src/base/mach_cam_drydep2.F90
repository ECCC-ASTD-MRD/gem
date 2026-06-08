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
! Fichier/File   : mach_cam_drydep2.ftn90
! Creation       : P.A. Makar (2011), Deji Akingunola (2018)
! Description    : Aerosol mass balance calculation for each layer.
!
! Extra info     : - Algorithm completely rewritten from original mach_tenddist
!                    August 2011) to track top and bottom layer interface
!                    fall velocities and remap the mass to the original
!                    layer structure
!
! Arguments:  IN
!
!                pdepv   -> gravitational settling velocity
!                thlev   -> Layer thickness [m]
!                rho     -> Layer air density [kg/m3]
!
!             IN/OUT
!
!                XROW    -> species concentrations
!
!             OUT
!
!                DRYFLX  -> Dry deposition flux
!
!============================================================================
!
!!if_on
subroutine mach_cam_drydep2(xrow, dryflx, pdepv, thlev, rho, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr, icom
!!if_off
   use chm_utils_mod,      only: chm_timestep
   use mach_cam_utils_mod, only: iae1, iae2, mwt_aero
   implicit none
!!if_on

   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent(out)   :: dryflx(pni, icom)
   real(kind=4),    intent(in)    :: pdepv (pni, pnk, isize)
   real(kind=4),    intent(in)    :: thlev (pni, pnk)
   real(kind=4),    intent(in)    :: rho   (pni, pnk)
!!if_off
!
! local variables
   integer(kind=4) :: nt, n, np, l, i, np2, k
   real(kind=8), dimension(pnk, icom * isize) :: xrownew, xrowold
   real(kind=8), dimension(pni, pnk)          :: oldtop, oldbot, zold
   real(kind=8), dimension(pni, isize)        :: minstep, dt
   integer(kind=4), dimension(pni, isize)     :: nsteps
   real(kind=8), dimension(pnk, isize)        :: znew, l1
   real(kind=8), dimension(pnk)               :: z1, z2
   real(kind=8), dimension(isize)             :: fz
   real(kind=4), dimension(icom)              :: fluxout
   real(kind=8), parameter                    :: small = 1.D-30   
   real(kind=8) :: colold, colnew, ratio, fluxout_step
   real(kind=8) :: dtmin
!  This version interpolates in mass/unit volume units, not mass mixing ratio
!
   xrownew  = 0.0d0
   dryflx = 0.0
!
   do i = 1, pni
      oldbot(i, pnk) = 0.0d0
      oldtop(i, pnk) = dble(thlev(i, pnk))
   end do
   do l = pnk-1, 1, -1
      do i = 1, pni
         oldbot(i, l) = oldtop(i, l+1)
         oldtop(i, l) = oldbot(i, l) + dble(thlev(i, l))
      end do
   end do

   do i = 1, pni
      zold(i, pnk) = oldbot(i, pnk)
      zold(i, 1) = oldtop(i, 1)
   end do
   do l = 2, pnk-1
      do i = 1, pni
         zold(i, l) = 0.5d0 * (oldtop(i, l) + oldbot(i, l))
      end do
   end do
!
!  Determine tentative minimum step size:
!
   minstep = 1.0D6
!  Find minimum timestep for the column at each particle size 
   do l = 1, pnk
      do n = 1,isize
         do i = 1, pni
            minstep(i, n) = min(minstep(i, n), 0.5D0 * &
                           (oldtop(i, l) - oldbot(i, l)) / dble(pdepv(i, l, n)))
            minstep(i, n) = min(minstep(i, n), dble(chm_timestep))
         end do
      end do
   end do
!
!  Determine final minimum step size:
!
   dtmin = 1.0D6
   do n = 1, isize
      do i = 1, pni
!  Substeps per model step:
         nsteps(i,n) = int(chm_timestep / real(minstep(i, n))) + 1
!  size of substep in seconds:
         dt(i,n) = dble(chm_timestep) / dble(nsteps(i, n))
         dtmin = min(dt(i, n), dtmin)
      end do
   end do
!
! Note that 1/2 the time required to transit the layer for the fastest moving
! particles has been chosen as the internal timestep. All particles of a given
! size will therefore remain in the same layer that they started in, which
! makes the backtrajectory calculation trivial for any given particle size.
! The backtrajectories all result in origin points higher than the destinations,
! than the destinations, but within the same layer as the destinations.
!
! For each column:
!
   do i = 1, pni 
      fluxout = 0.0
! Assign species invariant layer heights:
      do l = 2, pnk
         z2(l) = zold(i, l-1)
         z1(l) = zold(i, l)
      end do
! Set back-trajectory location for each particle size:
!  Note that zold at the model top will defacto be above the model top in this
!  configuration.
      do n = 1, isize
         do l = 1, pnk
            znew(l, n) = zold(i, l) + dt(i, n) * pdepv(i, l, n)
         end do
      end do

! Set size-dependant parameters for back-trajectories (Linear Interpolation):
!
      do n = 1, isize
         do l = 2, pnk
            l1(l,n) = (znew(l, n) - z1(l)) / (z2(l) - z1(l))
         end do
      end do
!
!  Diagnostics, for checking:
! Set initial concentrations in column:
      do np = iae1, iae2
         np2 = np - iae1 + 1
         do l = 1, pnk
            xrowold(l, np2) = xrow(i, l, np) * rho(i, l)
         end do
      end do
!
! Main advection loop
      do n = 1, isize
!
!  fz is the fraction of the top layer replaced by incoming mass from above, 
!  assumed small.
         fz(n) = min(max(dt(i, n) * dble(pdepv(i, 1, n)) / &
                 (zold(i, 1) - zold(i, 2)), 0.0D0), 1.0D0)
!
         do k = 1, nsteps(i, n)
! For each species:
            do nt = 1, icom
               np2 = isize * (nt - 1) + n
!
               do l = 2, pnk 
                  xrownew(l, np2) = max(small,                             &
                                   (xrowold(l-1, np2) - xrowold(l, np2)) * &
                                   l1(l, n) + xrowold(l, np2) )
               end do
!
! Model top: rather than extrapolate above the model top to determine particles
! coming down from above, instead use a weighted combination of the mass
! coming from above (assumed "small") and the mass remaining in the layer:
!
               xrownew(1, np2) = fz(n) * small + max(min(1.0D0 - fz(n), 1.0D0), 0.0D0) * xrowold(1, np2)
!
!
!  Mass conservation check in the column:
               colold = 0.D0
               colnew = 0.D0
               fluxout_step = dt(i, n) * dble(pdepv(i, pnk, n)) * 0.5d0 * &
                              (xrownew(pnk, np2) + xrowold(pnk, np2))
               do l = 2, pnk
                  colold = colold + ( zold(i, l-1) - zold(i, l)) * 0.5d0 *   &
                           (xrowold(l-1, np2) + xrowold(l, np2))
                  colnew = colnew + ( zold(i, l-1) - zold(i, l)) * 0.5d0 *   &
                           (xrownew(l-1, np2) + xrownew(l, np2))
               end do
               colnew = colnew + fluxout_step
               if (colnew > 0.0) then
                  ratio = colold / colnew
               else
                  ratio = 0.0D0
               end if
!
               do l = 1, pnk
                  xrownew(l, np2) = xrownew(l, np2) * ratio
! Update the concentration array:
                  xrowold(l, np2) = xrownew(l, np2)
               end do
!
               fluxout(nt) = fluxout(nt) + fluxout_step
!
! Next particle species
            end do  ! nt
! Do the next sub-step:
         end do !k
         
!  Evaluate the net dry deposition flux.
         do nt = 1, icom
            dryflx(i, nt) = dryflx(i, nt) + real(fluxout(nt))
         end do
!
!  Next particle bin size:
      end do ! n
!
!  Assign new concentrations:
      do np = iae1, iae2
         np2 = np - iae1 + 1
         do l = 1, pnk
            xrow(i, l, np) = real(xrownew(l, np2)) / rho(i, l)
         end do
      end do
!
!  Convert the units of the dry particle flux from kg/m2 to moles/m2:
! (Add -ve sign to align with previous GEMMACH convention of reporting deposition
!  as negative fluxes).
      do nt = 1, icom
         dryflx(i, nt) = -dryflx(i, nt) * 1.E3 / mwt_aero(nt)
      end do
!
   end do ! i

   return
end subroutine mach_cam_drydep2
