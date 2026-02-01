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

module MiMd
   implicit none
   public
   save

      type :: component
         character(len=256) name_S
         integer :: color, wpe0, rank
      end type component
      type (component) :: MiMd_world(100)
      
      character(len=256), dimension(:), allocatable :: names_S
      integer :: MiMd_Wmyproc, MiMd_Wnumproc, MiMd_ncolors
      integer, dimension(:,:), allocatable :: partners
      
end module MiMd
