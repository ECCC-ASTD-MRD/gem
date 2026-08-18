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
! Fichier/File   : chm_fst_closefile.ftn90
! Creation       : A. Kallaur & H. Landry - Mars 2007
! Description    : Close an opened fst file
!
! Extra info     :
!
! Arguments:
!            IN
!              file_unit --> File unit to close
!              lout      --> File unit for diagnostic outputs
!
!==============================================================================
!!if_on
subroutine chm_fst_closefile(file_unit)
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent(in) :: file_unit
!!if_off
!
! Local Variables
!
   integer(kind=4) :: err
   integer  fclos, fstfrm
   external fclos, fstfrm

   err = fstfrm(file_unit)
   err = fclos(file_unit)

end subroutine chm_fst_closefile
