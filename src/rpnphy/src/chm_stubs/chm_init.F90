!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2020 - Air Quality Research Division &
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
! Fichier/File   : chm_init.ftn90
! Creation       : A. Akingunola (Winter 2020)
! Description    : Initialization of the sundry chemistry parameters at the
!                  beginning of each execution of the model
!        ------------> STUB VERSION  STUB VERSION  STUB VERSION <------------ !
!
!==============================================================================
function chm_init(F_path_S) result(F_istat)
   implicit none
   character(len=*), intent(in) :: F_path_S !# data/tables dir
   integer ::  F_istat

   F_istat = 1
   !#Note: The next line is only to prevent compiler complain of unused var
   if (F_path_S == '__SCRAP__') print *,'(chm_init) path='//F_path_S

   return
end function chm_init
