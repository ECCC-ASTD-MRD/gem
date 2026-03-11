!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

module outusrdir

   implicit none
   public
   save
!
!revision
!______________________________________________________________________
!                                                                      |
!  VARIABLES FOR DEFINITION OF THE OUTPUT GRIDS (set_grid)             |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! OutUsrdir_MAX      | maximum number of usrdir that can be defined    |
! OutUsrdir_sets     | total number of sets of defined output grids    |
!  The following variables carry values for each defined output usrdir |
! OutUsrdir_id       | OutUsrdir_id(i) are the id of each defined grid |
! OutUsrdir_name_S   | OutUsrdir_name_S(i) are the name of each usrdir |
!----------------------------------------------------------------------
!
   integer, parameter :: OutUsrdir_MAX = 9, OutUsrdir_name_lenght = 20
   integer, dimension(OutUsrdir_MAX) :: OutUsrdir_id
   integer :: OutUsrdir_sets
   
   character(len=OutUsrdir_name_lenght), dimension(OutUsrdir_MAX) :: OutUsrdir_name_S

 end module outusrdir
