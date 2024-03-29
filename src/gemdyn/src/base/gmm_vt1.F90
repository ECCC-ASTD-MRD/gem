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
module gmm_vt1
   implicit none
   public
   save

      real, pointer, contiguous, dimension (:,:,:) ::  ut1 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  vt1 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  tt1 => null()
      real, pointer, contiguous, dimension (:,:)   ::  st1 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  wt1 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  qt1 => null()
      real, pointer, contiguous, dimension (:,:,:) :: zdt1 => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH), parameter:: gmmk_ut1_s   =  'URT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_vt1_s   =  'VRT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_tt1_s   =  'TT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_st1_s   =  'ST1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_wt1_s   =  'WT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_qt1_s   =  'QT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_zdt1_s  =  'ZDT1'

end module gmm_vt1
