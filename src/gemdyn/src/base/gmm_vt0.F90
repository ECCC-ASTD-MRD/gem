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
module gmm_vt0
   implicit none
   public
   save

      real, pointer, contiguous, dimension (:,:,:) ::  ut0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  vt0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  tt0 => null()
      real, pointer, contiguous, dimension (:,:)   ::  st0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  wt0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  qt0 => null()
      real, pointer, contiguous, dimension (:,:,:) :: zdt0 => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH), parameter:: gmmk_ut0_s   =  'URT0'
      character(len=MAXNAMELENGTH), parameter:: gmmk_vt0_s   =  'VRT0'
      character(len=MAXNAMELENGTH), parameter:: gmmk_tt0_s   =  'TT0'
      character(len=MAXNAMELENGTH), parameter:: gmmk_st0_s   =  'ST0'
      character(len=MAXNAMELENGTH), parameter:: gmmk_wt0_s   =  'WT0'
      character(len=MAXNAMELENGTH), parameter:: gmmk_qt0_s   =  'QT0'
      character(len=MAXNAMELENGTH), parameter:: gmmk_zdt0_s  =  'ZDT0'

end module gmm_vt0
