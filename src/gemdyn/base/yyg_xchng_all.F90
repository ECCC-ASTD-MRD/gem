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

!**s/r yyg_xchng_all - Exchanges all Yin-Yang boundary conditions

      subroutine yyg_xchng_all ()
      use dynkernel_options
      use gmm_vt1
      use gmm_pw
      use glb_ld
      use tr3d
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

      character(len=GMM_MAXNAMELENGTH) :: tr_name
      logical :: using_qt1
      integer :: istat,n
      real, pointer, contiguous, dimension(:,:,:) :: tr1
!
!----------------------------------------------------------------------
!
      using_qt1 = ( .not. Dynamics_hydro_L ) .or. &
                  ( trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H')

      do n= 1, Tr3d_ntr
         tr_name = 'TR/'//trim(Tr3d_name_S(n))//':P'
         istat = gmm_get(tr_name, tr1)
         call yyg_xchng (tr1, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .true., 'CUBIC', .false.)
      end do

      istat = gmm_get (gmmk_wt1_s , wt1)
      istat = gmm_get (gmmk_zdt1_s,zdt1)
      istat = gmm_get (gmmk_st1_s , st1)
      istat = gmm_get (gmmk_pw_uu_plus_s,pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s,pw_vv_plus)
      istat = gmm_get (gmmk_pw_tt_plus_s,pw_tt_plus)

      call yyg_xchng_vec_q2q ( pw_uu_plus,pw_vv_plus, &
                               l_minx,l_maxx,l_miny,l_maxy, G_nk )

      call yyg_xchng ( pw_tt_plus, l_minx,l_maxx,l_miny,l_maxy,&
                           l_ni, l_nj, G_nk, .false., 'CUBIC', .false. )

      call yyg_xchng ( wt1 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                           G_nk, .false., 'CUBIC', .false. )

      call yyg_xchng ( zdt1, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                           G_nk, .false., 'CUBIC', .false. )

      call yyg_xchng ( st1 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                           1,    .false., 'CUBIC', .false. )

      if (using_qt1) then
         istat = gmm_get(gmmk_qt1_s,  qt1)
         call yyg_xchng ( qt1, l_minx,l_maxx,l_miny,l_maxy,&
                          l_ni,l_nj, G_nk+1, .false., 'CUBIC', .false. )
      end if

!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_xchng_all
