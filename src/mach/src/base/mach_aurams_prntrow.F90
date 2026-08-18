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
! Fichier/File   : mach_aurams_prntrow.ftn90
! Creation       : M. Lazare, S.Gravel, GEM-MACH, July 2008
! Description    : Prints out row values (between il1 and il2) of array fld,
!                  for each of nlev levels at "jlat" latitude row inside routine.
!
! Extra info     :
!
! Arguments  IN
!              fld      --> array to be printed out.
!              name     --> name of array to be included in printout,
!                           i.e. "th" for throw.
!              loc      --> location where printout is occurring, i.e.
!                           "after vrtdfs".
!                           being displayed.
!              nlev     --> number of levels for variable.
!              loff     --> number of levels (including moon layer) before
!                           variable has values, i.e. loff=msg+1 for estg,
!                           loff=0 for throw, loff=1 for urow. for single-
!                           level fields, some level indicator (i.e. could be
!                           done in a vertical level loop).
!                           moon layer referenced to as level "0" in printout.
!              il1      --> starting longitude number for printout.
!              il2      --> finishing longitude number for printout.
!              jlat     --> current latitude row for printout.
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_aurams_prntrow (fld, name, loc, nlev, loff, il1, il2, jlat)
   use chm_ptopo_grid_mod, only: chm_ni
!!if_off
   use chm_utils_mod,      only: chm_lun_out, global_debug

   implicit none
!!if_on
   integer(kind=4),   intent(in) :: nlev
   real(kind=4),      intent(in) :: fld(chm_ni,nlev)
   character(len=7),  intent(in) :: name
   character(len=14), intent(in) :: loc
   integer(kind=4),   intent(in) :: loff
   integer(kind=4),   intent(in) :: il1
   integer(kind=4),   intent(in) :: il2
   integer(kind=4),   intent(in) :: jlat
!!if_off
!
! Local varialbes
!
   integer(kind=4)  :: il, lind, l
   logical(kind=4)  :: local_dbg

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

!  distinguish between single-level and multi-level fields.
   if (nlev == 1) then
      if(local_dbg) then
         write(chm_lun_out, 6100)loc, loff, name, jlat, (fld(il, 1), il = il1, il2)
      end if
 6100 format ('0', a14, ', level', i3, ' for variable ', a7, ' at latitude row', i3, ' :', /, (5(1x, 1pe24.16)))
   else
      lind = loff - 1
      if(local_dbg) then
         write(chm_lun_out, 6200)
      end if
 6200 format ('1')
      do l = 1, nlev
         if(local_dbg) then
            write(chm_lun_out, 6300)loc, l + lind, name, jlat, (fld(il, l), il = il1, il2)
         end if
 6300    format ('0', a14, ', level', i3, ' for variable ', a7, ' at latitude row', i3, ' :', /, (5(1x, 1pe24.16)))
  200   continue
      end do
   end if

   return
end
