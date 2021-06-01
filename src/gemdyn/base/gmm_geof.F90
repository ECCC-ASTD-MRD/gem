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
module gmm_geof
   implicit none
   public
   save

!______________________________________________________________________
!                                                                      |
!  GMM VARIABLES ASSOCIATED WITH GEOPHYSICAL FIELDS (set_geof)         |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! fis0               | Phi srf at current timestep                     |
! sls                | large scale s (SLEVE)                           |
! topo_low           | Low resolution analysis orography before growth |
! topo_high          | High resolution target orography after growth   |
!-----------------------------------------------------------------------
!
!
      real, pointer, contiguous, dimension (:,:) :: fis0      => null()
      real, pointer, contiguous, dimension (:,:) :: orols     => null()
      real, pointer, contiguous, dimension (:,:) :: sls       => null()
      real, pointer, contiguous, dimension (:,:,:) :: topo_low  => null()
      real, pointer, contiguous, dimension (:,:,:) :: topo_high => null()
      real, pointer, contiguous, dimension (:,:) :: me_full => null()
      real, pointer, contiguous, dimension (:,:) :: me_large => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH), parameter :: gmmk_fis0_s      = 'FIS0'
      character(len=MAXNAMELENGTH), parameter :: gmmk_orols_s     = 'OROLS'
      character(len=MAXNAMELENGTH), parameter :: gmmk_sls_s       = 'SLS'
      character(len=MAXNAMELENGTH), parameter :: gmmk_topo_low_s  = 'TOPOLOW'
      character(len=MAXNAMELENGTH), parameter :: gmmk_topo_high_s = 'TOPOHIGH'
      character(len=MAXNAMELENGTH), parameter :: gmmk_me_full_s   = 'MEFULL'
      character(len=MAXNAMELENGTH), parameter :: gmmk_me_large_s  = 'MELARGE'

end module gmm_geof
