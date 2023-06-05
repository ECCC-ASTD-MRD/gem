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
      subroutine set_vta
      use gmm_vta
      use glb_ld
      use lun
      use gmm_table
      use var_gmm
      implicit none
#include <arch_specific.hf>

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

      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_uta_s  ; GMM_tbl%ara(gmm_cnt)='UU' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)='DGFU'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_vta_s  ; GMM_tbl%ara(gmm_cnt)='VV' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)='DGFV'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_tta_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)='DGFT'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_qta_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)='DGFQ'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_wta_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)='DGFW'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_zdta_s ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)='DGFZ'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_sta_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='DGFS'
!
 1000 format( &
      /,'INITIALIZATION OF DIGITAL FILTER COMDECKS (S/R SET_VTA)', &
      /,'=======================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
