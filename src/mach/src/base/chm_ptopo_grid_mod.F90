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
! Fichier/File   : chm_ptopo_grid_mod.ftn90
! Creation       : H. Landry, Mai 2008
! Description    : Modules defining domain topology and local
!                  grid sizes.
!
! Extra info     :
!
!============================================================================
module chm_ptopo_grid_mod

   save

   integer(kind=4) :: chm_ni       ! GEM equiv. to p_ni (local x size - from physics)
   integer(kind=4) :: chm_nk       ! GEM equiv. to p_nk (local z size - from physics)
   integer(kind=4) :: pm_nk        ! number of vertical layers to which aerosol chemistry is applied 
                                   !        (pm_nk=chm_nk-nk_start_pm+1)
   integer(kind=4), parameter :: nkc = 3 ! number of canopy layers for shading effects 
   integer(kind=4) :: nkt, pm_nkc  ! Resolved model layers, plus canopy (nkc) layers
   
end module chm_ptopo_grid_mod
