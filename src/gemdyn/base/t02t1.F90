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

!**s/r t02t1 -  Rename time level t0 -> t1
!
      subroutine t02t1()
      use dynkernel_options
      use glb_ld
      use gmm_itf_mod
      use gmm_vt1
      use gmm_vt0
      use tr3d
      implicit none
#include <arch_specific.hf>
!
!object
!     Associate the variables at time t1 to the space on disk and memory
!     associated with the variables at time t0


      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: ut_list  = [ 'URT0', 'URT1' ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: vt_list  = [ 'VRT0', 'VRT1' ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: tt_list  = [ 'TT0' , 'TT1'  ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: st_list  = [ 'ST0' , 'ST1'  ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: wt_list  = [ 'WT0' , 'WT1'  ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: qt_list  = [ 'QT0' , 'QT1'  ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2), parameter :: zdt_list = [ 'ZDT0', 'ZDT1' ]
      character(len=GMM_MAXNAMELENGTH) , dimension(2) :: tr_list
      integer :: i, istat
!
!     ---------------------------------------------------------------
!
      istat = gmm_shuffle( ut_list)
      istat = gmm_shuffle( vt_list)
      istat = gmm_shuffle( tt_list)
      istat = gmm_shuffle( st_list)
      istat = gmm_shuffle(zdt_list)
      istat = gmm_shuffle( wt_list)

      if (.not. Dynamics_hydro_L .or. trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H') then
         istat = gmm_shuffle(qt_list)
      end if

      do i=1,Tr3d_ntr
         tr_list(1) = 'TR/'//trim(Tr3d_name_S(i))//':M'
         tr_list(2) = 'TR/'//trim(Tr3d_name_S(i))//':P'
         istat = gmm_shuffle(tr_list)
      end do

      istat = gmm_get (gmmk_ut0_s , ut0)
      istat = gmm_get (gmmk_vt0_s , vt0)
      istat = gmm_get (gmmk_tt0_s , tt0)
      istat = gmm_get (gmmk_st0_s , st0)
      istat = gmm_get (gmmk_wt0_s , wt0)
      istat = gmm_get (gmmk_qt0_s , qt0)
      istat = gmm_get (gmmk_zdt0_s, zdt0)

      istat = gmm_get (gmmk_ut1_s , ut1)
      istat = gmm_get (gmmk_vt1_s , vt1)
      istat = gmm_get (gmmk_tt1_s , tt1)
      istat = gmm_get (gmmk_st1_s , st1)
      istat = gmm_get (gmmk_wt1_s , wt1)
      istat = gmm_get (gmmk_qt1_s , qt1)
      istat = gmm_get (gmmk_zdt1_s, zdt1)
!
!     ---------------------------------------------------------------
!
      return
      end
