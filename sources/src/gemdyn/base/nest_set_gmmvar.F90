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

!**s/r nest_set_gmmvar - initialization of the commons for nesting variables
!		         within the Virtual Memory manager (VMM)
!
!#define SPY_VMM_CREATE spy_vmm_create

      subroutine nest_set_gmmvar
      use gmm_nest
      use lam_options
      use glb_ld
      use lun
      use tr3d
      use gmm_itf_mod
      use var_gmm
      implicit none
#include <arch_specific.hf>

!object
!	This subroutine initializes the commons containing the
!	keys used by the Virtual Memory Manager to identify the
!	nesting variables


      character(len=GMM_MAXNAMELENGTH) :: tr_name
      integer i,istat
      real, pointer, dimension (:,:,:) :: tr
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,1000)

      gmmk_nest_u_deb_s  = 'NEST_UA'
      gmmk_nest_v_deb_s  = 'NEST_VA'
      gmmk_nest_t_deb_s  = 'NEST_TA'
      gmmk_nest_s_deb_s  = 'NEST_SA'
      gmmk_nest_w_deb_s  = 'NEST_WA'
      gmmk_nest_q_deb_s  = 'NEST_QA'
      gmmk_nest_zd_deb_s = 'NEST_ZDA'
      gmmk_nest_fullme_deb_s = 'NEST_MEA'

      gmmk_nest_u_s      = 'NEST_U'
      gmmk_nest_v_s      = 'NEST_V'
      gmmk_nest_t_s      = 'NEST_T'
      gmmk_nest_s_s      = 'NEST_S'
      gmmk_nest_w_s      = 'NEST_W'
      gmmk_nest_q_s      = 'NEST_Q'
      gmmk_nest_zd_s     = 'NEST_ZD'
      gmmk_nest_fullme_s     = 'NEST_ME'

      gmmk_nest_u_fin_s  = 'NEST_UF'
      gmmk_nest_v_fin_s  = 'NEST_VF'
      gmmk_nest_t_fin_s  = 'NEST_TF'
      gmmk_nest_s_fin_s  = 'NEST_SF'
      gmmk_nest_w_fin_s  = 'NEST_WF'
      gmmk_nest_q_fin_s  = 'NEST_QF'
      gmmk_nest_zd_fin_s = 'NEST_ZDF'
      gmmk_nest_fullme_fin_s = 'NEST_MEF'

      istat = gmm_create(gmmk_nest_u_s , nest_u , meta3d_nk)
      istat = gmm_create(gmmk_nest_v_s , nest_v , meta3d_nk)
      istat = gmm_create(gmmk_nest_t_s , nest_t , meta3d_nk)
      istat = gmm_create(gmmk_nest_s_s , nest_s , meta2d   )
      istat = gmm_create(gmmk_nest_w_s , nest_w , meta3d_nk)
      istat = gmm_create(gmmk_nest_q_s , nest_q , meta3d_nk1)
      istat = gmm_create(gmmk_nest_zd_s, nest_zd, meta3d_nk)
      istat = gmm_create(gmmk_nest_fullme_s, nest_fullme, meta2d   )

      if (.not. Lam_ctebcs_L) then
         istat = gmm_create(gmmk_nest_u_deb_s , nest_u_deb , meta3d_nk)
         istat = gmm_create(gmmk_nest_v_deb_s , nest_v_deb , meta3d_nk)
         istat = gmm_create(gmmk_nest_t_deb_s , nest_t_deb , meta3d_nk)
         istat = gmm_create(gmmk_nest_s_deb_s , nest_s_deb , meta2d   )
         istat = gmm_create(gmmk_nest_w_deb_s , nest_w_deb , meta3d_nk)
         istat = gmm_create(gmmk_nest_q_deb_s , nest_q_deb , meta3d_nk1)
         istat = gmm_create(gmmk_nest_zd_deb_s, nest_zd_deb, meta3d_nk)
         istat = gmm_create(gmmk_nest_fullme_deb_s, nest_fullme_deb, meta2d   )

         istat = gmm_create(gmmk_nest_u_fin_s , nest_u_fin , meta3d_nk)
         istat = gmm_create(gmmk_nest_v_fin_s , nest_v_fin , meta3d_nk)
         istat = gmm_create(gmmk_nest_t_fin_s , nest_t_fin , meta3d_nk)
         istat = gmm_create(gmmk_nest_s_fin_s , nest_s_fin , meta2d   )
         istat = gmm_create(gmmk_nest_w_fin_s , nest_w_fin , meta3d_nk)
         istat = gmm_create(gmmk_nest_q_fin_s , nest_q_fin , meta3d_nk1)
         istat = gmm_create(gmmk_nest_zd_fin_s, nest_zd_fin, meta3d_nk)
         istat = gmm_create(gmmk_nest_fullme_fin_s, nest_fullme_fin, meta2d   )
      end if

      do i=1,Tr3d_ntr
         tr_name = Tr3d_name_S(i)
         nullify(tr)
         istat = gmm_create('NEST/'//trim(tr_name)//':C',tr,meta3d_nk)
         if (.not. Lam_ctebcs_L) then
            nullify(tr)
            istat = gmm_create('NEST/'//trim(tr_name)//':A',tr,meta3d_nk)
            nullify(tr)
            istat = gmm_create('NEST/'//trim(tr_name)//':F',tr,meta3d_nk)
         end if
      end do

      call nest_init_weight('M')
      call nest_init_weight('U')
      call nest_init_weight('V')
      call nest_init_weight('Q')

 1000 format( &
      /,'INITIALIZATION OF NESTING VARIABLE COMDECKS (S/R nest_set_gmmvar)', &
      /,'=================================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
