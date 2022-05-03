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
   use, intrinsic :: iso_fortran_env
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
! mcrhsint           | 2D Markov chain: rhs interpolation perturbation |
! mcsph1             | 3D Markov chain                                 |
! difut1             | horizontal diffusion u-tendency                 |
! difvt1             | horizontal diffusion v-tendency                 |
! ugwdt1             | gravity wave         u-tendency                 |
! vgwdt1             | gravity wave         v-tendency                 |
!----------------------------------------------------------------------
!
!
! dimension (LDIST_SHAPE)
      real, pointer, dimension (:,:)         :: mcrhsint=> null()
! dimension (l_ni, l_nj ,(l_nk+2))
      real, pointer, dimension (:,:,:)       :: mcsph1  => null()
! dimension LDIST_SHAPE,  l_nk+1)
      real, pointer, dimension (:,:,:)       :: difut1  => null()
      real, pointer, dimension (:,:,:)       :: difvt1  => null()
      real, pointer, dimension (:,:,:)       :: diout1  => null()
      real, pointer, dimension (:,:,:)       :: diovt1  => null()
      real, pointer, dimension (:,:,:)       :: ugwdt1  => null()
      real, pointer, dimension (:,:,:)       :: vgwdt1  => null()
      real, pointer, dimension (:,:,:)       :: ensdiv  => null()
      real, pointer, dimension (:,:,:)       :: ensvor  => null()
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
      real(kind=REAL64), pointer, dimension (:,:,:)         ::   pls      => null()
      real(kind=REAL64), pointer, dimension (:,:,:,:)       ::   plp     => null()

      integer, parameter :: MAXNAMELENGTH    =  32

      character(len=MAXNAMELENGTH), parameter:: gmmk_mcrhsint_s= 'MCRHSINT'
      character(len=MAXNAMELENGTH), parameter:: gmmk_mcsph1_s= 'MCSPH1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_difut1_s= 'DIFUT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_difvt1_s= 'DIFVT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_diout1_s= 'DIOUT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_diovt1_s= 'DIOVT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ugwdt1_s= 'UGWDT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_vgwdt1_s= 'VGWDT1'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ensdiv_s= 'ENSDIV'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ensvor_s= 'ENSVOR'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ar_s   = 'ARENS_S'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ai_s   = 'AIENS_S'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ar_p   = 'ARENS_P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_ai_p   = 'AIENS_P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_br_s   = 'BRENS_S'
      character(len=MAXNAMELENGTH), parameter:: gmmk_bi_s   = 'BIENS_S'
      character(len=MAXNAMELENGTH), parameter:: gmmk_br_p   = 'BRENS_P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_bi_p   = 'BIENS_P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_dumdum_s= 'DUMDUM'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pls_s  = 'P_LEGENS'
      character(len=MAXNAMELENGTH), parameter:: gmmk_plp_s  = 'P_LEGENP'
      
end module ens_gmm_var
