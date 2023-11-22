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
      use rmn_gmm
      use gmm_contiguous
      use gmm_vt1
      use gmm_vt0
      use gmm_pw
      use mem_tracers
      use tr3d
      implicit none
#include <arch_specific.hf>
!
!object
!     Associate the variables at time t1 to the space on disk and memory
!     associated with the variables at time t0

      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: ut_list  = [ 'URT0', 'URT1' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: vt_list  = [ 'VRT0', 'VRT1' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: tt_list  = [ 'TT0' , 'TT1'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: st_list  = [ 'ST0' , 'ST1'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: wt_list  = [ 'WT0' , 'WT1'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: qt_list  = [ 'QT0' , 'QT1'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: zdt_list = [ 'ZDT0', 'ZDT1' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: dyn_list = [ 'DYNT0', 'DYNT1' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_uulist = ['PW_UU:P','PW_UU:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_vvlist = ['PW_VV:P','PW_VV:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_ttlist = ['PW_TT:P','PW_TT:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_pmlist = ['PW_PM:P','PW_PM:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_ptlist = ['PW_PT:P','PW_PT:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_gzlist = ['PW_GZ:P','PW_GZ:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_melist = ['PW_ME:P','PW_ME:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_p0list = ['PW_P0:P','PW_P0:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_pm8list = ['PW_PM8:P','PW_PM8:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_p08list = ['PW_P08:P','PW_P08:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: tracers = ['TRACERS:t0','TRACERS:t1']
      character(len=GMM_MAXNAMELENGTH) , dimension(2) :: tr_list
      integer :: i, istat, dim
      real, pointer, dimension (:) :: tr_tmp
!
!     ---------------------------------------------------------------
!
!$omp single
      istat = gmm_shuffle(  ut_list)
      istat = gmm_shuffle(  vt_list)
      istat = gmm_shuffle(  tt_list)
      istat = gmm_shuffle(  st_list)
      istat = gmm_shuffle( zdt_list)
      istat = gmm_shuffle(  wt_list)
      istat = gmm_shuffle( dyn_list)
      istat = gmm_shuffle(pw_uulist)
      istat = gmm_shuffle(pw_vvlist)
      istat = gmm_shuffle(pw_ttlist)
      istat = gmm_shuffle(pw_pmlist)
      istat = gmm_shuffle(pw_pm8list)
      istat = gmm_shuffle(pw_ptlist)
      istat = gmm_shuffle(pw_gzlist)
      istat = gmm_shuffle(pw_melist)
      istat = gmm_shuffle(pw_p0list)
      istat = gmm_shuffle(pw_p08list)
      istat = gmm_shuffle(tracers)

      if ((.not. Dynamics_hydro_L) .or. Dynamics_hauteur_L) then
         istat = gmm_shuffle(qt_list)
      end if

      tr_tmp => trt0
      trt0   => trt1
      trt1   => tr_tmp

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
      do i=1,Tr3d_ntr
         tracers_P(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt1((i-1)*dim+1:)
         tracers_M(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt0((i-1)*dim+1:)
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
      istat = gmm_get ('DYNT0', dynt0)

      istat = gmm_get (gmmk_ut1_s , ut1)
      istat = gmm_get (gmmk_vt1_s , vt1)
      istat = gmm_get (gmmk_tt1_s , tt1)
      istat = gmm_get (gmmk_st1_s , st1)
      istat = gmm_get (gmmk_wt1_s , wt1)
      istat = gmm_get (gmmk_qt1_s , qt1)
      istat = gmm_get (gmmk_zdt1_s, zdt1)
      istat = gmm_get ('DYNT1', dynt1)

      istat = gmm_get(gmmk_pw_uu_plus_s  ,pw_uu_plus)
      istat = gmm_get(gmmk_pw_vv_plus_s  ,pw_vv_plus)
      istat = gmm_get(gmmk_pw_tt_plus_s  ,pw_tt_plus)
      istat = gmm_get(gmmk_pw_pm_plus_s  ,pw_pm_plus)
      istat = gmm_get(gmmk_pw_pt_plus_s  ,pw_pt_plus)
      istat = gmm_get(gmmk_pw_gz_plus_s  ,pw_gz_plus)
      istat = gmm_get(gmmk_pw_me_plus_s  ,pw_me_plus)
      istat = gmm_get(gmmk_pw_p0_plus_s  ,pw_p0_plus)
      istat = gmm_get(gmmk_pw_pm_plus_8_s  ,pw_pm_plus_8)
      istat = gmm_get(gmmk_pw_p0_plus_8_s  ,pw_p0_plus_8)

      istat = gmm_get(gmmk_pw_uu_moins_s  ,pw_uu_moins)
      istat = gmm_get(gmmk_pw_vv_moins_s  ,pw_vv_moins)
      istat = gmm_get(gmmk_pw_tt_moins_s  ,pw_tt_moins)
      istat = gmm_get(gmmk_pw_pt_moins_s  ,pw_pt_moins)
      istat = gmm_get(gmmk_pw_gz_moins_s  ,pw_gz_moins)
      istat = gmm_get(gmmk_pw_pm_moins_s  ,pw_pm_moins)
      istat = gmm_get(gmmk_pw_me_moins_s  ,pw_me_moins)
      istat = gmm_get(gmmk_pw_p0_moins_s  ,pw_p0_moins)
      istat = gmm_get(gmmk_pw_pm_moins_8_s ,pw_pm_moins_8)
      istat = gmm_get(gmmk_pw_p0_moins_8_s ,pw_p0_moins_8)
!$omp end single
!
!     ---------------------------------------------------------------
!
      return
      end
