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
! Projet/Project  : GEM-MACH
! Fichier/File    : chm_mjrpts_get_size.ftn90
! Creation        : 
! Description     : Retreive major point vector length
!
! Extra info      :
!
! Arguments :
!
!          IN
!             file_unit      -->  file number for major point source file
!
!          OUT
!             numsrcs        -->  number of sources
!
!==============================================================================
!
!!if_on
subroutine chm_mjrpts_get_size(file_unit, numsrcs, istatus)
!!if_off
   use chm_utils_mod, only: global_debug, chm_lun_out
   implicit  none
!!if_on
   integer(kind=4), intent (in) :: file_unit
   integer(kind=4), intent(out) :: numsrcs
   integer(kind=4), intent(out) :: istatus
!!if_off
!
! Local variables
!
   integer(kind=4)   :: j1, j2, key
   logical(kind=4)   :: local_dbg
!
!  External Functions
   integer, external :: fstinf

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
   key = fstinf(file_unit, numsrcs, j1, j2, -1, ' ', -1, -1, -1, ' ', 'LAT')
   if (numsrcs > 0) then
      istatus = 0
   else
      istatus = -1
      if (local_dbg) then
         write(chm_lun_out, *) 'ERROR: expecting number of point sources > 0'
      end if
   end if

   return
end subroutine chm_mjrpts_get_size