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
! Fichier/File   : mach_hetv_mod.ftn90
! Description    : Modules defining common variables for mach_hetv subroutines
!
! Extra info     :
!
!============================================================================

module mach_hetv_mod

   implicit none
   save
   integer(kind=4), parameter     :: maxhet = 3, maxghet = 3
!
   integer(kind=4), parameter     :: iter = 100
   integer(kind=4), parameter     :: itero = 4
   integer(kind=4), parameter     :: ndiv = 5
   real(kind=8), parameter        :: eps = 1.0d-09
   real(kind=8), parameter        :: eps2 = 1.0d-04
   real(kind=8), parameter        :: lwmin = 1.0d-20
   real(kind=8), parameter        :: tstd = 298.15d0
   real(kind=8), parameter        :: smrt = 1.0d-05
   real(kind=8), parameter        :: small = 1.0d-30
   real(kind=8), parameter        :: lolimit = 1.0d-24
  !gas constant, atm->Pa
   real(kind=8), parameter        :: rg = 8.3144D0 / 1.01325D+05 ! gas constant, atm->pa
   character(len = 20), parameter :: method = 'KM'

   integer(kind=4)                :: itermax

end module mach_hetv_mod
