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
! Fichier/File   : mach_aurams_tsysms.ftn90
! Creation       : S. Gong, W. Gong, S.Gravel V. Bouchet, B. Pabla, GEM-MACH, July 2008
! Description    : Prints out row values (between il1 and il2) of array fld,
!                  for each of nlev levels at "jlat" latitude row inside routine.
!
! Extra info     :
!
! Arguments  IN
!              xrow(icdn) --> Tracer array (kg/kg) but at index xrow(icdn)=no cloud droplet/m3
!                             i.e. "th" for throw.
!              loc        --> location where printout is occurring, i.e.
!                             "after vrtdfs" being displayed.
!              pm_nk     --> number of vertical levels for variable.
!              chm_ni     --> number of equally-spaced longitudes to be printed out.
!              il1        --> starting longitude number for printout.
!              il2        --> finishing longitude number for printout.
!              nt         --> This is izon = 3
!
!            OUT
!              TMASS      --> Total Mass at site
!=============================================================================
!
!!if_on
subroutine mach_aurams_tsysms (xrow, TMASS, il1, il2, nt, loc)
   use chm_ptopo_grid_mod, only: chm_ni, pm_nk
   use mach_cam_utils_mod, only: ntr
!!if_off
   use chm_utils_mod,      only: chm_lun_out, global_debug
   use mach_cam_utils_mod, only: isize, iae1
   implicit none
!!if_on
   real(kind=4),    intent(out) :: tmass
   integer(kind=4), intent (in) :: il1
   integer(kind=4), intent (in) :: il2
   integer(kind=4), intent (in) :: nt
   integer(kind=4), intent (in) :: loc
   real(kind=4),    intent (in) :: xrow(chm_ni,pm_nk,ntr)
!!if_off
!
!  local variables
!
   integer(kind=4) :: n, l, il
   logical(kind=4) :: local_dbg

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   tmass = 0.0
   do n = (nt - 1) * isize + 1 + (iae1 - 1), nt * isize + (iae1 - 1)
      do l = 1, pm_nk
         do il = il1, il2
            tmass = tmass + xrow(il, l, n)
         end do
      end do
   end do

   if (local_dbg) then
      write (chm_lun_out, *)' totmas at site ', loc, ' is ', tmass
   end if
   return
end
