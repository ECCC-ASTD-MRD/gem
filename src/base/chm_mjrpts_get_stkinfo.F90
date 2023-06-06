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
! Fichier/File    : chm_mjrpts_get_stkinfo.ftn90
! Creation        : A. Kallaur, H. Landry, S. Menard - janvier/fevrier 2007
! Description     : Retreive major point source stack information
!
! Extra info      :
!
! Arguments :
!
!          IN
!             file_unit      -->  file number for major point source file
!
!          OUT
!             all characteristics of the stacks (lat, lon, etc.)
!
!==============================================================================
!
!!if_on
subroutine chm_mjrpts_get_stkinfo(file_unit, numsrcs, istatus, &
                                  lat, lon, tem, hgt, dia, vel, typ)
!!if_off
   use chm_utils_mod, only: global_debug, chm_lun_out
   implicit  none
!!if_on
   integer(kind=4),        intent (in) :: numsrcs
   integer(kind=4),        intent (in) :: file_unit
   integer(kind=4),        intent(out) :: istatus
   real(kind=4),           intent(out) :: lat(numsrcs)
   real(kind=4),           intent(out) :: lon(numsrcs)
   real(kind=4), optional, intent(out) :: tem(numsrcs)
   real(kind=4), optional, intent(out) :: hgt(numsrcs)
   real(kind=4), optional, intent(out) :: dia(numsrcs)
   real(kind=4), optional, intent(out) :: vel(numsrcs)
   real(kind=4), optional, intent(out) :: typ(numsrcs)
!!if_off
!
! Local variables
!
   integer(kind=4)   :: j0, j1, j2, ier
!
!  External Functions
!
   integer, external :: fstlir
!
   istatus = 0
   ier = fstlir(lat, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'LAT')
   if (ier < 0) then
      write(chm_lun_out, *) 'ERROR: missing stack info for variable LAT'
      istatus = min(ier, istatus)
   end if
   ier = fstlir(lon, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'LON')
   if (ier < 0) then
      write(chm_lun_out, *) 'ERROR: missing stack info for variable LON'
      istatus = min(ier, istatus)
   end if
   if (present(tem)) then
      ier = fstlir(tem, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'TEM')
      if (ier < 0) then
         write(chm_lun_out, *) 'ERROR: missing stack info for variable TEM'
         istatus = min(ier, istatus)
      end if
   end if
   if (present(hgt)) then
      ier = fstlir(hgt, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'HGT')
      if (ier < 0) then
         write(chm_lun_out, *) 'ERROR: missing stack info for variable HGT'
         istatus = min(ier, istatus)
      end if
   end if
   if (present(dia)) then
      ier = fstlir(dia, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'DIA')
      if (ier < 0) then
         write(chm_lun_out, *) 'ERROR: missing stack info for variable DIA'
         istatus = min(ier, istatus)
      end if
   end if
   if (present(vel)) then
      ier = fstlir(vel, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'VEL')
      if (ier < 0) then
         write(chm_lun_out, *) 'ERROR: missing stack info for variable VEL'
         istatus = min(ier, istatus)
      end if
   end if
   if (present(typ)) then
      ier = fstlir(typ, file_unit, j0, j1, j2, -1, ' ', -1, -1, -1, ' ', 'TYP')
      if (ier < 0) then
         write(chm_lun_out, *) 'WARNING: no TYP record found, assume anthropogenic emissions'
         typ = 100.0
      end if
   end if

   return
end subroutine chm_mjrpts_get_stkinfo
