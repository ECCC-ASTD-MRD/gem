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
module gmm_tracers
   implicit none
   public
   save
!
   real, pointer, contiguous, dimension (:,:,:) :: fld_mono => null()
   real, pointer, contiguous, dimension (:,:,:) :: fld_cub  => null()
   real, pointer, contiguous, dimension (:,:,:) :: fld_lin  => null()
   real, pointer, contiguous, dimension (:,:,:) :: fld_min  => null()
   real, pointer, contiguous, dimension (:,:,:) :: fld_max  => null()
   real, pointer, contiguous, dimension (:,:,:) :: fld_cub_o=> null()
   real, pointer, contiguous, dimension (:,:,:) :: fld_cub_i=> null()
   real, pointer, contiguous, dimension (:,:,:) :: airm1    => null()
   real, pointer, contiguous, dimension (:,:,:) :: airm0    => null()
   real, pointer, contiguous, dimension (:,:,:) :: pkps     => null()

   integer, parameter :: MAXNAMELENGTH = 32

   character(len=MAXNAMELENGTH) :: gmmk_mono_s, gmmk_cub_s, gmmk_lin_s, gmmk_min_s, gmmk_max_s, &
                                   gmmk_cub_o_s, gmmk_cub_i_s, gmmk_airm1_s, gmmk_airm0_s, gmmk_pkps_s

end module gmm_tracers
