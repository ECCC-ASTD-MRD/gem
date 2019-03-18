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
module ens_gmm_var
   implicit none
   public
   save
!
!______________________________________________________________________
!                                                                      |
!  GMM variables for Markov chains                                     |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! mcsph1             | 3D Markov chain                                 |
! difut1             | horizontal diffusion u-tendency                 |
! difvt1             | horizontal diffusion v-tendency                 |
! ugwdt1             | gravity wave         u-tendency                 |
! vgwdt1             | gravity wave         v-tendency                 |
!----------------------------------------------------------------------
!
!
! dimension (l_ni, l_nj ,(l_nk+2))
      real, pointer, dimension (:,:,:)     :: mcsph1  => null()
! dimension LDIST_SHAPE,  l_nk+1)
      real,    pointer, dimension (:,:,:)     :: difut1  => null()
      real,    pointer, dimension (:,:,:)     :: difvt1  => null()
      real,    pointer, dimension (:,:,:)     :: diout1  => null()
      real,    pointer, dimension (:,:,:)     :: diovt1  => null()
      real,    pointer, dimension (:,:,:)     :: ugwdt1  => null()
      real,    pointer, dimension (:,:,:)     :: vgwdt1  => null()
      real,    pointer, dimension (:,:,:)     :: ensdiv  => null()
      real,    pointer, dimension (:,:,:)     :: ensvor  => null()
! special shapes for spectral coefficients
      real, pointer, dimension (:,:,:)       :: ar_p     => null()
      real, pointer, dimension (:,:,:)       :: ai_p     => null()
      real, pointer, dimension (:,:)         :: ar_s     => null()
      real, pointer, dimension (:,:)         :: ai_s     => null()
      real, pointer, dimension (:,:,:)       :: br_p     => null()
      real, pointer, dimension (:,:,:)       :: bi_p     => null()
      real, pointer, dimension (:,:)         :: br_s     => null()
      real, pointer, dimension (:,:)         :: bi_s     => null()

! special shape for random number generator
      integer, pointer, dimension (:,:)       :: dumdum  => null()

! Legendre polynomial for SKEB Markov chain
      real*8, pointer, dimension (:,:,:)       ::   plg     => null()

      integer, parameter :: MAXNAMELENGTH    =  32

      character(len=MAXNAMELENGTH) :: gmmk_mcsph1_s
      character(len=MAXNAMELENGTH) :: gmmk_difut1_s,gmmk_difvt1_s
      character(len=MAXNAMELENGTH) :: gmmk_diout1_s,gmmk_diovt1_s
      character(len=MAXNAMELENGTH) :: gmmk_ugwdt1_s,gmmk_vgwdt1_s
      character(len=MAXNAMELENGTH) :: gmmk_ensdiv_s,gmmk_ensvor_s
      character(len=MAXNAMELENGTH) :: gmmk_ar_s,gmmk_ai_s,gmmk_br_s,gmmk_bi_s
      character(len=MAXNAMELENGTH) :: gmmk_ar_p,gmmk_ai_p,gmmk_br_p,gmmk_bi_p
      character(len=MAXNAMELENGTH) :: gmmk_dumdum_s,gmmk_plg_s

end module ens_gmm_var
