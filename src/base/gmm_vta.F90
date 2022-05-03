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
module gmm_vta
   implicit none
   public
   save
!
!______________________________________________________________________
!                                                                      |
!  GMM variables that are digital filtered                             |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
!  uta               | U: x component of velocity                      |
!  vta               | V: y component of velocity                      |
!  tta               | T: temperature                                  |
! zdta               | Zdot: generalized vertical velocity             |
!  sta               | s: log of surface pressure minus a constant     |
!--------------------|-------------------------------------------------|
!  wta               | w: z component of velocity                      |
!  qta               | q: log of non-hydrostatic perturbation pressure |
!--------------------|-------------------------------------------------|
! trta               | passive tracer(s)                               |
!--------------------|-------------------------------------------------|
!

   real, pointer, contiguous, dimension (:,:,:)   :: uta  => null()
   real, pointer, contiguous, dimension (:,:,:)   :: vta  => null()
   real, pointer, contiguous, dimension (:,:,:)   :: tta  => null()
   real, pointer, contiguous, dimension (:,:,:)   :: qta  => null()
   real, pointer, contiguous, dimension (:,:,:)   :: wta  => null()
   real, pointer, contiguous, dimension (:,:,:)   :: zdta => null()
   real, pointer, contiguous, dimension (:,:)     :: sta  => null()
   real, pointer, contiguous, dimension (:,:,:,:) :: trta => null()

   integer, parameter :: MAXNAMELENGTH = 32

   character(len=MAXNAMELENGTH) :: gmmk_vta_s, gmmk_tta_s, gmmk_qta_s ,&
                 gmmk_zdta_s, gmmk_sta_s, gmmk_wta_s, gmmk_uta_s, gmmk_trta_s

end module gmm_vta
