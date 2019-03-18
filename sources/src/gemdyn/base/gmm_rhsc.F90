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
module gmm_rhsc
   implicit none
   public
   save

!
!______________________________________________________________________
!                                                                      |
!  GMM VARIABLES ASSOCIATED WITH RHS (right-hand side equation)        |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
!   Contributions of time t1 to rhs of different equations             |
!--------------------|-------------------------------------------------|
! rhsu               | momentum equation along x                       |
! rhsv               | momentum equation along y                       |
! rhst               | thermodynamic equation                          |
! rhsc               | continuity equation                             |
! rhsw               | momentum equation along z                       |
! rhsf               | vertical velocity equation                      |
!--------------------|-------------------------------------------------|
! rhsp               |        Helmholtz equation                       |
! rhsb               |   boundary condition in opentop                 |
!--------------------|-------------------------------------------------|
! ruw1               | Ru interpolated from U_grid to G_grid           |
! ruw2               | advective contributions for U                   |
! rvw1               | Rv interpolated from V_grid to G_grid           |
! rvw2               | advective contributions for V                   |
!----------------------------------------------------------------------
!

      real, pointer, contiguous, dimension (:,:,:) :: rhsu => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsv => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhst => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsc => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsw => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsf => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsp => null()
      real, pointer, contiguous, dimension (:,:  ) :: rhsb => null()
      real, pointer, contiguous, dimension (:,:,:) :: ruw1 => null()
      real, pointer, contiguous, dimension (:,:,:) :: rvw1 => null()
      real, pointer, contiguous, dimension (:,:,:) :: ruw2 => null()
      real, pointer, contiguous, dimension (:,:,:) :: rvw2 => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsx => null()
      real, pointer, contiguous, dimension (:,:,:) :: rhsq => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH) :: &
               gmmk_rhsu_s,gmmk_rhsv_s,gmmk_rhst_s,gmmk_rhsc_s,&
               gmmk_rhsw_s,gmmk_rhsf_s,gmmk_rhsp_s,gmmk_rhsb_s,&
               gmmk_ruw1_s,gmmk_rvw1_s,gmmk_ruw2_s,gmmk_rvw2_s,&
               gmmk_rhsx_s,gmmk_rhsq_s

end module gmm_rhsc
