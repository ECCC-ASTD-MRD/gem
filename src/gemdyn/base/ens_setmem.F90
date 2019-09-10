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

!**s/r ens_setmem - initialize ensemble prevision system
!
      subroutine ens_setmem (l_ni, l_nj, l_nk, Lun_out)
      use ens_param
      use ens_gmm_dim
      use ens_gmm_var
      use ens_options
      use gmm_itf_mod
      use var_gmm
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: l_ni, l_nj, l_nk, Lun_out
!
!     author
!     Lubos Spacek - February 2010
!
!     revision
! v4_12 - Spacek L.        - Initial version
! v_4.1.3 - N. Gagnon      - Change name of most parameters in the NAMELIST
!

!#include "ens_gmm_dim.cdk"
!#include "ens_param.cdk"

      integer :: istat
!-------------------------------------------------------------------
!
      if (.not.Ens_conf) return

      if (Lun_out > 0) write(Lun_out,6010)

      call gmm_build_meta3D(meta3d_sh2,          &
                            1,l_ni,0,0,l_ni,     &
                            1,l_nj,0,0,l_nj,     &
                            1,l_nk,0,0,l_nk,     &
                            0,GMM_NULL_FLAGS)
       call gmm_build_meta3D(meta3d_ar_p,                      &
                            1,Ens_ptp_lmax,0,0,Ens_ptp_lmax,  &
                            1,Ens_ptp_mmax,0,0,Ens_ptp_mmax,  &
                            1,Ens_ptp_ncha,0,0,Ens_ptp_ncha,&
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta3D(meta3d_br_p,                       &
                            1,Ens_ptp_lmax,0,0,Ens_ptp_lmax,  &
                            1,Ens_ptp_mmax,0,0,Ens_ptp_mmax,  &
                            1,Ens_ptp_ncha,0,0,Ens_ptp_ncha,&
                            0,GMM_NULL_FLAGS)
     call gmm_build_meta3D(meta3d_ai_p,                      &
                            1,Ens_ptp_lmax,0,0,Ens_ptp_lmax,  &
                            1,Ens_ptp_mmax,0,0,Ens_ptp_mmax,  &
                            1,Ens_ptp_ncha,0,0,Ens_ptp_ncha,&
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta3D(meta3d_bi_p,                       &
                            1,Ens_ptp_lmax,0,0,Ens_ptp_lmax,  &
                            1,Ens_ptp_mmax,0,0,Ens_ptp_mmax,  &
                            1,Ens_ptp_ncha,0,0,Ens_ptp_ncha,&
                            0,GMM_NULL_FLAGS)

      call gmm_build_meta2D(meta2d_ar_s,                      &
                            1,Ens_skeb_l,0,0,Ens_skeb_l,  &
                            1,Ens_skeb_m,0,0,Ens_skeb_m,  &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta2d_br_s,                       &
                            1,Ens_skeb_l,0,0,Ens_skeb_l,  &
                            1,Ens_skeb_m,0,0,Ens_skeb_m,  &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta2d_ai_s,                      &
                            1,Ens_skeb_l,0,0,Ens_skeb_l,  &
                            1,Ens_skeb_m,0,0,Ens_skeb_m,  &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta2d_bi_s,                       &
                            1,Ens_skeb_l,0,0,Ens_skeb_l,  &
                            1,Ens_skeb_m,0,0,Ens_skeb_m,  &
                            0,GMM_NULL_FLAGS)

      call gmm_build_meta3D(meta3d_plg,                       &
                            1,Ens_skeb_nlat,0,0,Ens_skeb_nlat,  &
                            1,Ens_skeb_l,0,0,Ens_skeb_l,  &
                            1,Ens_skeb_m,0,0,Ens_skeb_m,&
                            0,GMM_NULL_FLAGS)

      call gmm_build_meta2D(meta2d_dum,                       &
                            1,36,0,0,36,                      &
                            1,2*MAX2DC,0,0,2*MAX2DC,&
                            0,GMM_NULL_FLAGS)
      gmmk_mcsph1_s= 'MCSPH1'
      gmmk_difut1_s= 'DIFUT1'
      gmmk_difvt1_s= 'DIFVT1'
      gmmk_diout1_s= 'DIOUT1'
      gmmk_diovt1_s= 'DIOVT1'
      gmmk_ugwdt1_s= 'UGWDT1'
      gmmk_vgwdt1_s= 'VGWDT1'
      gmmk_ensdiv_s= 'ENSDIV'
      gmmk_ensvor_s= 'ENSVOR'
      gmmk_ar_s   = 'ARENS_S'
      gmmk_ai_s   = 'AIENS_S'
      gmmk_ar_p   = 'ARENS_P'
      gmmk_ai_p   = 'AIENS_P'
      gmmk_br_s   = 'BRENS_S'
      gmmk_bi_s   = 'BIENS_S'
      gmmk_br_p   = 'BRENS_P'
      gmmk_bi_p   = 'BIENS_P'
      gmmk_dumdum_s= 'DUMDUM'
      gmmk_plg_s  = 'P_LEGEN'

      istat = gmm_create(gmmk_mcsph1_s,mcsph1,meta3d_sh2,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'mcsph1'
      istat = gmm_create(gmmk_difut1_s,difut1,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'difut1'
      istat = gmm_create(gmmk_difvt1_s,difvt1,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'difvt1'
      istat = gmm_create(gmmk_diout1_s,diout1,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'diout1'
      istat = gmm_create(gmmk_diovt1_s,diovt1,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'diovt1'
      istat = gmm_create(gmmk_ugwdt1_s,ugwdt1,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ugwdt1'
      istat = gmm_create(gmmk_vgwdt1_s,vgwdt1,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'vgwdt1'
      istat = gmm_create(gmmk_ensdiv_s,ensdiv,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ensdiv'
      istat = gmm_create(gmmk_ensvor_s,ensvor,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ensvor'

      istat = gmm_create(gmmk_ar_s   ,ar_s   ,meta2d_ar_s,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ar_s'
      istat = gmm_create(gmmk_br_s   ,br_s   ,meta2d_br_s,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'br_s'
      istat = gmm_create(gmmk_ai_s   ,ai_s   ,meta2d_ai_s,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ai_s'
      istat = gmm_create(gmmk_bi_s   ,bi_s   ,meta2d_bi_s,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'bi_s'
      istat = gmm_create(gmmk_ar_p   ,ar_p   ,meta3d_ar_p,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ar_p'
      istat = gmm_create(gmmk_br_p   ,br_p   ,meta3d_br_p,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'br_p'
      istat = gmm_create(gmmk_ai_p   ,ai_p   ,meta3d_ai_p,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'ai_p'
      istat = gmm_create(gmmk_bi_p   ,bi_p   ,meta3d_bi_p,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'bi_p'

      istat = gmm_create(gmmk_plg_s,plg,meta3d_plg,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'plg'

      istat = gmm_create(gmmk_dumdum_s,dumdum,meta2d_dum,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat))write(*,6000)'dum'

 6000 format('ens_set_mem at gmm_create(',A,')')
 6010 format(/,'INITIALIZATION OF MEMORY FOR ENSEMBLES (S/R ENS_SETMEM)' &
             /(55('=')))
!
!-------------------------------------------------------------------
!
      return
      end
