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
      use mem_tracers
      implicit none
#include <arch_specific.hf>

      logical :: using_qt1
      integer :: n
!
!----------------------------------------------------------------------
!
      using_qt1 = ( .not. Dynamics_hydro_L ) .or. Dynamics_hauteur_L

      do n= 1, Tr3d_ntr
         call yyg_int_xch_scal (tracers_P(n)%pntr, G_nk, .true., 'CUBIC', .false.)
      end do

      call yyg_xchng_vec_q2q ( pw_uu_plus,pw_vv_plus, &
                               l_minx,l_maxx,l_miny,l_maxy, G_nk )

      call yyg_int_xch_scal ( pw_tt_plus, G_nk, .false., 'CUBIC', .false. )

      call yyg_int_xch_scal ( wt1 , G_nk, .false., 'CUBIC', .false. )

      call yyg_int_xch_scal ( zdt1, G_nk, .false., 'CUBIC', .false. )

      call yyg_int_xch_scal ( st1 , 1   , .false., 'CUBIC', .false. )

      if (using_qt1) then
         call yyg_int_xch_scal ( qt1, G_nk+1, .false., 'CUBIC', .false. )
      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_xchng_all
