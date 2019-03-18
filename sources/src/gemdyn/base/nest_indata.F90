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
!
      subroutine nest_indata (F_datev_S)
      use gmm_nest
      use gmm_geof
      use inp_mod
      use glb_ld
      use lun
      use gmm_itf_mod
      use gem_timing
      implicit none
#include <arch_specific.hf>

      character(len=*) F_datev_S

!author
!     Michel Desgagne   - Spring 2002
!
!revision
! v3_01 - Desgagne M.     - initial version
! v3_03 - Tanguay M.      - Adjoint Lam configuration
! v3_30 - Lee V.          - Hollow cubes and acid test for LAM
! v4_03 - Lee/Desgagne    - ISST
! v4_05 - Plante A.       - Top nesting
! v4_05 - Lepine M.       - VMM replacement with GMM
! v4_10 - Lee V.          - Remove TRNES on tracers,zd,w
! v4_4  - Plante A.       - Add computation of wt1
!                           Catch error on gmm


      integer istat
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      istat = gmm_get(gmmk_nest_u_fin_s ,nest_u_fin )
      istat = gmm_get(gmmk_nest_v_fin_s ,nest_v_fin )
      istat = gmm_get(gmmk_nest_w_fin_s ,nest_w_fin )
      istat = gmm_get(gmmk_nest_t_fin_s ,nest_t_fin )
      istat = gmm_get(gmmk_nest_zd_fin_s,nest_zd_fin)
      istat = gmm_get(gmmk_nest_s_fin_s ,nest_s_fin )
      istat = gmm_get(gmmk_nest_q_fin_s ,nest_q_fin )
      istat = gmm_get(gmmk_nest_fullme_fin_s,nest_fullme_fin)

      nest_zd_fin=0. ; nest_w_fin=0. ; nest_q_fin= 0.

      call gemtime_start ( 26, 'NEST_input', 10 )
      call inp_data ( nest_u_fin , nest_v_fin, nest_w_fin, nest_t_fin,&
                      nest_zd_fin, nest_s_fin, nest_q_fin            ,&
                      nest_fullme_fin, l_minx,l_maxx,l_miny,l_maxy   ,&
                      G_nk, .true., 'NEST/', ':F', F_datev_S)
      call gemtime_stop (26)

      call diag_zd_w ( nest_zd_fin, nest_w_fin                       ,&
                       nest_u_fin, nest_v_fin, nest_t_fin, nest_s_fin,&
                       l_minx,l_maxx,l_miny,l_maxy, G_nk             ,&
                       .not.Inp_zd_L, .not.Inp_w_L )
!
!     ---------------------------------------------------------------
!
 1000 format(3X,'GETTING DATA FROM NEST TO BCS: (S/R NEST_INDATA)')

      return
      end

