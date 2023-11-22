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

module ldnh
   implicit none
   public
   save

!
!**comdeck Ldnh.cdk
!
!______________________________________________________________________
!                                                                      |
!  LOCAL DIMENSIONS AND DATA TOPOLOGY WITHOUT HALO (set_transpose)     |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Ldnh_ni            | same as l_ni but without halo                   |
! Ldnh_nj            | same as l_nj but without halo                   |
! Ldnh_minx          | same as l_minx but without halo                 |
! Ldnh_maxx          | same as l_maxx but without halo                 |
! Ldnh_miny          | same as l_miny but without halo                 |
! Ldnh_maxy          | same as l_maxy but without halo                 |
!----------------------------------------------------------------------
!
   integer :: Ldnh_ni, Ldnh_nj,  Ldnh_i0, Ldnh_j0
   integer :: Ldnh_minx, Ldnh_maxx, Ldnh_miny, Ldnh_maxy

end module ldnh
