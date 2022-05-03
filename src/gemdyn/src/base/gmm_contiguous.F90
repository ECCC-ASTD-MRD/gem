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
      module gmm_contiguous
      implicit none
      public
      save

      type :: gmm_contiguous_pntrs
         real, dimension(:,:  ), pointer :: pntr_2d
         real, dimension(:,:,:), pointer :: pntr_3d
         end type gmm_contiguous_pntrs
      integer gmm_nbplans
      real, pointer, contiguous, dimension (:) :: dynt2 => null()
      real, pointer, contiguous, dimension (:) :: dynt1 => null()
      real, pointer, contiguous, dimension (:) :: dynt0 => null()
      type(gmm_contiguous_pntrs), allocatable :: timlvl2(:)
      type(gmm_contiguous_pntrs), allocatable :: timlvl1(:)
      type(gmm_contiguous_pntrs), allocatable :: timlvl0(:)

      end module gmm_contiguous
