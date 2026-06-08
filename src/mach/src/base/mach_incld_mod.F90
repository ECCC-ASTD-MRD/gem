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
! Fichier/File   : mach_incld_mod.ftn90
! Creation       : H. Landry, Sept 2008
! Description    : Modules defining common variables for mach_incld subroutines
!
! Extra info     :
!
!============================================================================

module mach_incld_mod

   save

   real(kind=4),               parameter :: gmin  = 1.e-20
   real(kind=4),               parameter :: aqmin = 1.e-20
   real(kind=4),               parameter :: dtmax = 75.
   real(kind=4), dimension(2), parameter :: dtsrt = (/ 0.001, 22.5/)
   logical(kind=4),            parameter :: domnfe = .false.
   logical(kind=4),            parameter :: doh2o2 = .true.
   logical(kind=4),            parameter :: dorooh = .true.
   real(kind=4)                          :: incld_dtmin

end module mach_incld_mod
