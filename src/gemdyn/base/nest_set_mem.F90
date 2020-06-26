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

!**s/r nest_set_mem - Allocate memory and set pointers
!		                for nesting variables

      subroutine nest_set_mem
      use mem_nest
      use lam_options
      use glb_ld
      use lun
      use tr3d
      implicit none
#include <arch_specific.hf>

      integer dimTot,dimHor,dim3d
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,1000)

      dimHor = (l_maxx-l_minx+1) * (l_maxy-l_miny+1)
      dim3d  = dimHor * l_nk
      dimTot = (6+Tr3d_ntr)*dim3d + 3*dimHor

      allocate (nest_now(dimTot))

      nest_u (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_now(        1:)
      nest_v (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_now(dim3d+  1:)
      nest_t (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_now(dim3d*2+1:)
      nest_w (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_now(dim3d*3+1:)
      nest_zd(l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_now(dim3d*4+1:)
      nest_q (l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> nest_now(dim3d*5+1:)
      nest_s     (l_minx:l_maxx,l_miny:l_maxy     )=> nest_now(dim3d*6+  dimHor+1:)
      nest_fullme(l_minx:l_maxx,l_miny:l_maxy     )=> nest_now(dim3d*6+2*dimHor+1:)
      nest_tr(l_minx:l_maxx,l_miny:l_maxy,1:Tr3d_ntr*l_nk)=> &
                                                      nest_now(dim3d*6+3*dimHor+1:)
      if (.not. Lam_ctebcs_L) then

      allocate (nest_deb(dimTot),nest_fin(dimTot))

      nest_u_deb (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_deb(         1:)
      nest_v_deb (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_deb(dim3d+  1:)
      nest_t_deb (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_deb(dim3d*2+1:)
      nest_w_deb (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_deb(dim3d*3+1:)
      nest_zd_deb(l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_deb(dim3d*4+1:)
      nest_q_deb (l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> nest_deb(dim3d*5+1:)
      nest_s_deb     (l_minx:l_maxx,l_miny:l_maxy     )=> nest_deb(dim3d*6+  dimHor+1:)
      nest_fullme_deb(l_minx:l_maxx,l_miny:l_maxy     )=> nest_deb(dim3d*6+2*dimHor+1:)
      nest_tr_deb(l_minx:l_maxx,l_miny:l_maxy,1:Tr3d_ntr*l_nk)=> &
                                                      nest_deb(dim3d*6+3*dimHor+1:)
      nest_u_fin (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_fin(         1:)
      nest_v_fin (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_fin(dim3d+  1:)
      nest_t_fin (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_fin(dim3d*2+1:)
      nest_w_fin (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_fin(dim3d*3+1:)
      nest_zd_fin(l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> nest_fin(dim3d*4+1:)
      nest_q_fin (l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> nest_fin(dim3d*5+1:)
      nest_s_fin (l_minx:l_maxx,l_miny:l_maxy         )=> nest_fin(dim3d*6+  dimHor+1:)
      nest_fullme_fin(l_minx:l_maxx,l_miny:l_maxy     )=> nest_fin(dim3d*6+2*dimHor+1:)
      nest_tr_fin(l_minx:l_maxx,l_miny:l_maxy,1:Tr3d_ntr*l_nk)=> &
                                                      nest_fin(dim3d*6+3*dimHor+1:)
      endif

      allocate (nest_weightm(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1),&
                nest_weightq(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1),&
                nest_weightu(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1),&
                nest_weightv(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1))

      call nest_init_weight(nest_weightm,0,0,-1,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)
      call nest_init_weight(nest_weightu,-1,0,0,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)
      call nest_init_weight(nest_weightv,0,-1,0,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)
      call nest_init_weight(nest_weightq,0,0,-1,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)

 1000 format( &
      /,'INITIALIZATION OF NESTING VARIABLE MEMORY', &
      /,'==============================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
