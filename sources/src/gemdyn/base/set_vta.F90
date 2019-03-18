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

!**s/r set_vta
!

!
      subroutine set_vta
      use gmm_vta
      use glb_ld
      use lun
      use gmm_itf_mod
      use var_gmm
      implicit none
#include <arch_specific.hf>
!
!author
!     alain patoine - march 1994
!
!revision
! v2_00 - Desgagne/Lee   - initial MPI version (from setvta v1_03)
! v2_21 - J. P. Toviessi - rename some model output variables
! v2_30 - Edouard S.     - remove pi' at the top
! v2_31 - Desgagne M.    - re-introduce 3D tracers*
! v4_05 - Lepine M.      - VMM replacement with GMM
!


      integer :: istat
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1000)

      nullify(uta,vta,wta,tta,qta,zdta,sta)

      gmmk_uta_s  = 'DIGF_UU'
      gmmk_vta_s  = 'DIGF_VV'
      gmmk_wta_s  = 'DIGF_WW'
      gmmk_tta_s  = 'DIGF_TT'
      gmmk_qta_s  = 'DIGF_Q '
      gmmk_zdta_s = 'DIGF_ZD'
      gmmk_sta_s  = 'DIGF_S '

      istat = gmm_create(gmmk_uta_s,uta,meta3d_nk,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(uta)'

      istat = gmm_create(gmmk_vta_s,vta,meta3d_nk,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(vta)'

      istat = gmm_create(gmmk_wta_s,wta,meta3d_nk,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(wta)'

      istat = gmm_create(gmmk_tta_s,tta,meta3d_nk,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(tta)'

      istat = gmm_create(gmmk_qta_s,qta,meta3d_nk1,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(qta)'

      istat = gmm_create(gmmk_zdta_s,zdta,meta3d_nk,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(zdta)'

      istat = gmm_create(gmmk_sta_s,sta,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_vta ERROR at gmm_create(sta)'
!
 1000 format( &
      /,'INITIALIZATION OF DIGITAL FILTER COMDECKS (S/R SET_VTA)', &
      /,'=======================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
