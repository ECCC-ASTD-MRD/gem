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
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

!
!______________________________________________________________________
!                                                                      |
!  GMM variables at TIME t0                                            |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
!  ut0               | U: x component of velocity                      |
!  vt0               | V: y component of velocity                      |
!  tt0               | T: temperature                                  |
!  st0               | s: log of surface pressure minus a constant     |
!  wt0               | w:      dz/dt                                   |
!  qt0               | q: log of non-hydrostatic perturbation pressure |
!--------------------|-------------------------------------------------|
! zdt0               | Zetadot:dZeta/dt: generalized vertical velocity |
!--------------------|-------------------------------------------------|
! trt0               | passive tracer(s)                               |
!----------------------------------------------------------------------
!

      real, pointer, contiguous, dimension (:,:,:) ::  ut0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  vt0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  tt0 => null()
      real, pointer, contiguous, dimension (:,:)   ::  st0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  wt0 => null()
      real, pointer, contiguous, dimension (:,:,:) ::  qt0 => null()
      real, pointer, contiguous, dimension (:,:,:) :: zdt0 => null()
      real(kind=REAL64),pointer,dimension (:)     ::  p00 => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH) :: gmmk_wt0_s, gmmk_tt0_s, gmmk_zdt0_s, &
                                      gmmk_st0_s, gmmk_qt0_s, gmmk_ut0_s, gmmk_vt0_s, gmmk_p00_s

end module gmm_vt0
