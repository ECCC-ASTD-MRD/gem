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

!**s/r t02t1 -  Rename time level t0 -> t2
!
      subroutine t02t2 ()
      use dynkernel_options
      use glb_ld
      use gmm_itf_mod
      use gmm_contiguous
      use gmm_vt2
      use gmm_vt0
      use mem_tracers
      use tr3d
      implicit none

!object
!     Associate the variables at time t2 to the space on disk and memory
!     associated with the variables at time t0

      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: ut_list  = [ 'URT0', 'URT2' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: vt_list  = [ 'VRT0', 'VRT2' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: tt_list  = [ 'TT0' , 'TT2'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: st_list  = [ 'ST0' , 'ST2'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: wt_list  = [ 'WT0' , 'WT2'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: qt_list  = [ 'QT0' , 'QT2'  ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: zdt_list = [ 'ZDT0', 'ZDT2' ]
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: dyn_list = [ 'DYNT0', 'DYNT2' ]

      character(len=GMM_MAXNAMELENGTH) , dimension(2) :: tr_list
      integer :: i, istat, dim
      real, pointer, dimension (:) :: tr_tmp
!
!     ---------------------------------------------------------------
!
      istat = gmm_shuffle(  ut_list)
      istat = gmm_shuffle(  vt_list)
      istat = gmm_shuffle(  tt_list)
      istat = gmm_shuffle(  st_list)
      istat = gmm_shuffle( zdt_list)
      istat = gmm_shuffle(  wt_list)
      istat = gmm_shuffle( dyn_list)

      if ((.not. Dynamics_hydro_L) .or. Dynamics_hauteur_L) then
         istat = gmm_shuffle(qt_list)
      end if

      tr_tmp => trt0
      trt0   => trt2
      trt2   => tr_tmp

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
      do i=1,Tr3d_ntr
         tracers_t2(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt2((i-1)*dim+1:)
         tracers_M(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt0((i-1)*dim+1:)
         tr_list(1) = 'TR/'//trim(Tr3d_name_S(i))//':M'
         tr_list(2) = 'TR/'//trim(Tr3d_name_S(i))//':t2'
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

      istat = gmm_get (gmmk_ut2_s , ut2)
      istat = gmm_get (gmmk_vt2_s , vt2)
      istat = gmm_get (gmmk_tt2_s , tt2)
      istat = gmm_get (gmmk_st2_s , st2)
      istat = gmm_get (gmmk_wt2_s , wt2)
      istat = gmm_get (gmmk_qt2_s , qt2)
      istat = gmm_get (gmmk_zdt2_s, zdt2)
      istat = gmm_get ('DYNT2', dynt2)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine t02t2
