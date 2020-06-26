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

!**s/r nest_indata - Read and process nesting data during LAM
!                    integration for LBC.

      subroutine nest_indata (F_datev_S)
      use dynkernel_options
      use mem_nest
      use gmm_geof
      use inp_mod
      use glb_ld
      use lun
      use tr3d
      use vGrid_Descriptors, only : vgrid_descriptor
      use gem_timing
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(IN) :: F_datev_S
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      call gemtime_start ( 26, 'NEST_input', 10 )

      call inp_data ( nest_u_fin , nest_v_fin, nest_w_fin, nest_t_fin,&
                      nest_zd_fin, nest_s_fin, nest_tr_fin,&
                      nest_fullme_fin,.true.,F_datev_S    ,&
                      l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr )

      call derivate_data ( nest_zd_fin, nest_w_fin, nest_u_fin       ,&
                      nest_v_fin, nest_t_fin , nest_s_fin, nest_q_fin,&
                              l_minx,l_maxx,l_miny,l_maxy, G_nk,&
                              .not.Inp_zd_L, .not.Inp_w_L )

      call gemtime_stop (26)
!
!     ---------------------------------------------------------------
!
 1000 format(3X,'GETTING DATA FROM NEST TO BCS: (S/R NEST_INDATA)')

      return
      end

