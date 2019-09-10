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
module gmm_orh
   implicit none
   public
   save
!
!______________________________________________________________________
!                                                                      |
! GMM variables associated with the Crank-Nicholson procedure          |
! and used to compute the RHS equation (set_crni)                      |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
!   Storing space for contributions of time t1                         |
!   to rhs of different equations                                      |
!  orhsu             | momentum equation along x                       |
!  orhsv             | momentum equation along y                       |
!  orhst             | thermodynamic equation                          |
!  orhsc             | continuity equation                             |
!  orhsw             | momentum equation along z                       |
!  orhsf             | vertical velocity equation                      |
!--------------------|-------------------------------------------------|
! Orh_icn            | current iteration of C-N procedure              |
! Orh_crank_L        | switch: .true. if current time step is C-N      |
!----------------------------------------------------------------------
!

      logical Orh_crank_L
      integer Orh_icn

      real, pointer, contiguous, dimension (:,:,:) :: orhsu => null()
      real, pointer, contiguous, dimension (:,:,:) :: orhsv => null()
      real, pointer, contiguous, dimension (:,:,:) :: orhst => null()
      real, pointer, contiguous, dimension (:,:,:) :: orhsc => null()
      real, pointer, contiguous, dimension (:,:,:) :: orhsw => null()
      real, pointer, contiguous, dimension (:,:,:) :: orhsf => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH) ::  &
           gmmk_orhsu_s, gmmk_orhsv_s, gmmk_orhst_s,&
           gmmk_orhsc_s, gmmk_orhsw_s, gmmk_orhsf_s

end module gmm_orh
