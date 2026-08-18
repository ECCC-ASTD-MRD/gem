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
! Fichier/File   : chm_fst_openfile.ftn90
! Creation       : A. Kallaur & H. Landry - Mars 2007
! Description    : Open an fst file
!
! Extra info     : There is a duplicate version of this code in the
!                  EMIS_PREP code library called em_open_fstfile.ftn90
!                  whose functionality is practically identical.
!
! Arguments:
!            IN
!              file_name -->  speaks for itself
!              file_type -->  RPN type value for fnom (usually something like 'RND+OLD')
!              options   -->  RPN options for fstouv  (usually something like 'RND')
!              lout      --> (File unit for disgnostic messages)
!
!            IN/OUT
!              file_unit --> the file unit to open (an initial value of 0 will
!                            assign a unit automatically - Recommended)
!
!==============================================================================
!
!!if_on
subroutine chm_fst_openfile(FILE_UNIT, file_name, file_type, options, istatus)
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent(inout) :: file_unit
   integer(kind=4), intent  (out) :: istatus
   character(*),    intent   (in) :: file_name
   character(*),    intent   (in) :: file_type
   character(*),    intent   (in) :: options
!!if_off
   integer  fnom, fstouv
   external fnom, fstouv

   istatus = 0
   istatus=fnom (file_unit, trim(file_name), file_type, 0)
   if (istatus >= 0) then
      istatus=fstouv (file_unit, options )
   endif

end subroutine chm_fst_openfile
