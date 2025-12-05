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
! Fichier/File   : chm_mjrpts_sortinfo_mod.ftn90
! Creation       : A. Kallaur, H. Landry, Mai 2008
! Description    : Modules defining variables and arrays containing pertinent
!                  information to the sorting (by i and j) of the
!                  major point source emissions fields.
!
! Extra info     :
!
!============================================================================

module chm_mjrpts_sortinfo_mod

!     Name              Description
!======================================================================
!
!   ndata               Number of characteristics of each stack
!
   save

   integer(kind=4), parameter  :: ndata = 7
   integer(kind=4), parameter  :: i_gilc = 1, i_gjlc = 2, i_hgt = 3, &
                                  i_dia = 4,  i_tem = 5,  i_vel = 6, i_typ = 7
   integer(kind=4)             :: lnb_sources 
   integer(kind=4)             :: nb_me_species
   real(kind=4),    dimension(:,:), pointer     :: Lstack_Info, Lstack_emis
   integer(kind=4), dimension(:),   allocatable :: me_species_index
   
end module chm_mjrpts_sortinfo_mod
