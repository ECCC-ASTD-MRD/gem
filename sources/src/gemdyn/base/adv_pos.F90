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
module adv_pos
   implicit none
   public
   save

   real, pointer, contiguous, dimension (:,:,:) ::  pxt => null()
   real, pointer, contiguous, dimension (:,:,:) ::  pyt => null()
   real, pointer, contiguous, dimension (:,:,:) ::  pzt => null()

   real, allocatable, dimension (:,:,:) ::  pxmu
   real, allocatable, dimension (:,:,:) ::  pymu
   real, allocatable, dimension (:,:,:) ::  pzmu

   real, allocatable, dimension (:,:,:) ::  pxmv
   real, allocatable, dimension (:,:,:) ::  pymv
   real, allocatable, dimension (:,:,:) ::  pzmv

   real, pointer, contiguous, dimension (:,:,:) ::  pxm_s => null()
   real, pointer, contiguous, dimension (:,:,:) ::  pym_s => null()
   real, pointer, contiguous, dimension (:,:,:) ::  pzm_s => null()

   real, allocatable, dimension (:,:)    :: pxtn
   real, allocatable, dimension (:,:)    :: pytn
   real, allocatable, dimension (:,:)    :: pztn

end module adv_pos
