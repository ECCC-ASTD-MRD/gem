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
module gmm_vt2
   implicit none
   public
   save

      real, pointer, contiguous, dimension (:,:,:) ::  ut2 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  vt2 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  tt2 => null()
      real, pointer, contiguous, dimension (:,:)   ::  st2 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  wt2 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  qt2 => null()
      real, pointer, contiguous, dimension (:,:,:) :: zdt2 => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH), parameter:: gmmk_ut2_s   =  'URT2'
      character(len=MAXNAMELENGTH), parameter:: gmmk_vt2_s   =  'VRT2'
      character(len=MAXNAMELENGTH), parameter:: gmmk_tt2_s   =  'TT2'
      character(len=MAXNAMELENGTH), parameter:: gmmk_st2_s   =  'ST2'
      character(len=MAXNAMELENGTH), parameter:: gmmk_wt2_s   =  'WT2'
      character(len=MAXNAMELENGTH), parameter:: gmmk_qt2_s   =  'QT2'
      character(len=MAXNAMELENGTH), parameter:: gmmk_zdt2_s  =  'ZDT2'

end module gmm_vt2
