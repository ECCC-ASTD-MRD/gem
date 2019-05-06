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

!  s/r fislh_hybrid   - Generates A and B of the hybrid coordinate
!                    Also sets Z and other related vertical parameters.
!
      subroutine fislh_hybrid ( F_hybuser, Nk )
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_put,&
           vgd_levels,VGD_OK,VGD_ERROR,vgd_print, VGD_LEN_NAME, VGD_NO_REF_NOMVAR, &
           vgd_stda76
      use vgrid_wb, only: vgrid_wb_put
      use gmm_pw
      use HORgrid_options
      use VERgrid_options
      use dyn_fisl_options
      use glb_ld
      use ctrl
      use cstv
      use lun
      use dimout
      use out_mod
      use levels
      use ver
      use wb_itf_mod
      implicit none
#include <arch_specific.hf>

      integer Nk
      real, dimension(Nk) :: F_hybuser        !user-specified hybrid coordinate values
!
!Authors: Claude Girard & Andre Plante, July 2017
!

      character(len=32), parameter  :: VGRID_M_S  = 'ref-m'
      character(len=32), parameter  :: VGRID_T_S  = 'ref-t'

      character(len=32) :: REFP0_S, REFP0_LS_S, dumc
      integer k,istat,pnip1,options_readwrite,options_readonly,ier
      integer, dimension(:), pointer :: wkpti
      real, dimension(:), pointer :: std_p_prof=>null(),wkpt
      real*8, parameter :: zero=0.d0, one=1.d0, half=0.5d0
      real*8, dimension(:), pointer :: wkpt8
      character(len=VGD_LEN_NAME) :: rfls_S
!     __________________________________________________________________
!
      allocate(   Ver_hyb%m(G_nk+1),       Ver_hyb%t(G_nk+1), &
           Ver_std_p_prof%m(G_nk+1),Ver_std_p_prof%t(G_nk+1), &
                  Ver_ip1%m(G_nk+1),       Ver_ip1%t(G_nk+1), &
                  Ver_a_8%m(G_nk+1),       Ver_a_8%t(G_nk+1), &
                  Ver_b_8%m(G_nk+1),       Ver_b_8%t(G_nk+1), &
                  Ver_c_8%m(G_nk+1),       Ver_c_8%t(G_nk+1), &
                Ver_z_8%m(0:G_nk+1),     Ver_z_8%t(0:G_nk+1), &
                                         Ver_z_8%x(0:G_nk  ), &
                 Ver_dz_8%m(G_nk  ),      Ver_dz_8%t(G_nk  ), &
                Ver_idz_8%m(G_nk  ),     Ver_idz_8%t(G_nk  ), &
               Ver_dbdz_8%m(G_nk  ),    Ver_dbdz_8%t(G_nk  ), &
               Ver_dcdz_8%m(G_nk  ),    Ver_dcdz_8%t(G_nk  ), &
                 Ver_wp_8%m(G_nk  ),      Ver_wm_8%m(G_nk  ), &
                Ver_onezero(G_nk+1),      Ver_zeronk(G_nk  ), &
               Ver_wpstar_8(G_nk  ),    Ver_wmstar_8(G_nk  ), &
              Ver_Tstar_8%m(G_nk+1),   Ver_Tstar_8%t(G_nk  ), &
                  Ver_wpA_8(G_nk  ),       Ver_wmA_8(G_nk  ) )

      Ver_code = 6

      ! Construct vertical coordinate
      Level_kind_ip1 = 21
      Level_version  = 1

      istat = vgd_new ( Ver_vgdobj, kind=Level_kind_ip1,&
                   version=Level_version, hyb=F_hybuser,&
               rcoef1=Hyb_rcoef(1), rcoef2=Hyb_rcoef(2),&
               rcoef3=Hyb_rcoef(3), rcoef4=Hyb_rcoef(4),&
                                         dhm=0., dht=0. )

      if (Lun_debug_L) istat = vgd_print(Ver_vgdobj)

      call handle_error_l(istat==VGD_OK,'fislh_hybrid','coordinate construction failed')

      ier = vgd_get(Ver_vgdobj,'RFLS large scale reference field',rfls_S,quiet=.true.)
      Schm_sleve_L = trim(rfls_S) /= VGD_NO_REF_NOMVAR

      ! Retrieve information required to fill model arrays
      nullify(wkpt,wkpti,wkpt8)
      if (vgd_get(Ver_vgdobj,'CA_M - vertical A coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_a_8%m = wkpt8(1:size(Ver_a_8%m)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(Ver_vgdobj,'CB_M - vertical B coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_b_8%m = wkpt8(1:size(Ver_b_8%m)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(Ver_vgdobj,'CA_T - vertical A coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_a_8%t = wkpt8(1:size(Ver_a_8%t)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(Ver_vgdobj,'CB_T - vertical B coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_b_8%t = wkpt8(1:size(Ver_b_8%t)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(Ver_vgdobj,'CC_M - vertical C coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_c_8%m = wkpt8(1:size(Ver_c_8%m)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(Ver_vgdobj,'CC_T - vertical C coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_c_8%t = wkpt8(1:size(Ver_c_8%t)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(Ver_vgdobj,'VCDM - vertical coordinate (m)'   ,wkpt) /= VGD_OK) istat = VGD_ERROR
      Ver_hyb%m = wkpt(1:size(Ver_hyb%m)); deallocate(wkpt); nullify(wkpt)
      if (vgd_get(Ver_vgdobj,'VCDT - vertical coordinate (t)'   ,wkpt) /= VGD_OK) istat = VGD_ERROR
      Ver_hyb%t = wkpt(1:size(Ver_hyb%t)); deallocate(wkpt); nullify(wkpt)
      if (vgd_get(Ver_vgdobj,'VIPM - level ip1 list (m)'        ,wkpti) /= VGD_OK) istat = VGD_ERROR
      Ver_ip1%m = wkpti(1:size(Ver_ip1%m)); deallocate(wkpti); nullify(wkpti)
      if (vgd_get(Ver_vgdobj,'VIPT - level ip1 list (t)'        ,wkpti) /= VGD_OK) istat = VGD_ERROR
      Ver_ip1%t = wkpti(1:size(Ver_ip1%t)); deallocate(wkpti); nullify(wkpti)

!      print*,'istat, VGD_ERROR',istat,VGD_ERROR

      print*,'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
      if( vgd_stda76(Ver_vgdobj, Ver_ip1%m, std_p_prof, 'PRESSURE') /= VGD_OK) istat = VGD_ERROR
      Ver_std_p_prof%m=std_p_prof
      deallocate(std_p_prof)
      if( vgd_stda76(Ver_vgdobj, Ver_ip1%t, std_p_prof, 'PRESSURE') /= VGD_OK) istat = VGD_ERROR
      Ver_std_p_prof%t=std_p_prof
      deallocate(std_p_prof)
      call handle_error_l(istat==VGD_OK, &
                'fislh_hybrid','retrieving coordinate info')

      print*,"TODO in fislh_hybrid compute Cstv_ptop_8?"
      Cstv_ptop_8=-1.

      !-------------------------------
      ! Define z(m/t/x) from A(m/t):
      !-------------------------------
      ! Level top is at level 1 plus half of delta between level 1 and 2.
      Ver_z_8%m(0) = 0.5d0 * ( 3.d0*Ver_a_8%m(1) - Ver_a_8%m(2) )
      do k = 1, G_nk+1
         Ver_z_8%m(k) = Ver_a_8%m(k)
      end do

     !Define the positions of true thermo levels
      Ver_z_8%t(0) = Ver_z_8%m(0)
      do k = 1, G_nk+1
         Ver_z_8%t(k) = Ver_a_8%t(k)
      end do

     !Define the positions of zeta_dot
      Ver_z_8%x(0) = Ver_z_8%m(0)
      do k = 1, G_nk-1
         Ver_z_8%x(k) = Ver_z_8%t(k)
      end do
      Ver_z_8%x(G_nk)=zero

      Ver_zmin_8 = Ver_z_8%m(G_nk+1)
      Ver_zmax_8 = Ver_z_8%m(0)

      if ( Ctrl_canonical_dcmip_L ) then
         Cstv_pref_8 = 100000.d0
         Cstv_ptop_8 = Cstv_pref_8 * exp(-Ver_z_8%m(0)/8780.2)
      end if

!     ----------------------
!     Compute dz, 1/dz, dBdz
!     ----------------------

      do k=1,G_nk
           Ver_dz_8%m(k) = Ver_z_8%x(k) - Ver_z_8%t(k-1)
          Ver_idz_8%m(k) = one/Ver_dz_8%m(k)
           Ver_dz_8%t(k) = Ver_z_8%m(k+1) - Ver_z_8%m(k)
          Ver_idz_8%t(k) = one/Ver_dz_8%t(k)
      end do

!     -------------------------------------------------------
!     Compute AVERGING WEIGHTS FROM THERMO TO MOMENTUM LEVELS
!     -------------------------------------------------------

      do k=1,G_nk
         Ver_wm_8%m(k) = (Ver_z_8%t(k)-Ver_z_8%m(k))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
         Ver_wmA_8(k)  = Ver_wm_8%m(k)
         Ver_wp_8%m(k) = one-Ver_wm_8%m(k)
         Ver_wpA_8(k)  = one-Ver_wm_8%m(k)
      end do
!
!     SPECIAL WEIGHTS due to last thermo level
!
      Ver_wmstar_8 = zero
      Ver_wpstar_8 = one
      Ver_wmstar_8(G_nk)=(Ver_z_8%x(G_nk)-Ver_z_8%t(G_nk)) &
                        /(Ver_z_8%x(G_nk)-Ver_z_8%x(G_nk-1))
      Ver_wpstar_8(G_nk)=one-Ver_wmstar_8(G_nk)
      Ver_wmA_8(G_nk)  = Ver_wp_8%m(G_nk)*Ver_wmstar_8(G_nk)+Ver_wm_8%m(G_nk)
      Ver_wpA_8(G_nk)  = one-Ver_wmA_8(G_nk)

!     -------------------------------------------------------
!     Initialize Ver_onezero
!     -------------------------------------------------------

      Ver_onezero=1.
      Ver_onezero(1)=0.

      Ver_zeronk=1.
      Ver_zeronk(G_nk)=0.

!     ----------------------------------------------------------
!     Save Ver_vgdobj and ip1m/t for output
!     ----------------------------------------------------------
      REFP0_S = 'PW_P0:P'  !# gmmk_pw_p0_plus_s !NOTE: could gmmk_* be defined as parameters in a .cdk, this way it could be used here and would be more consistent
      REFP0_LS_S = ' '
      if (Schm_sleve_L) REFP0_LS_S = 'PW_P0_LS'  !# gmmk_pw_p0_ls_s
      istat = vgrid_wb_put(VGRID_M_S, Ver_vgdobj, Ver_ip1%m,  &
           REFP0_S, REFP0_LS_S, F_overwrite_L=.true.)
      istat = vgrid_wb_put(VGRID_T_S, Ver_vgdobj, Ver_ip1%t,  &
           REFP0_S, REFP0_LS_S, F_overwrite_L=.true.)

      options_readwrite = WB_IS_LOCAL
      options_readonly = options_readwrite + WB_REWRITE_NONE

      istat= wb_put('model/Vgrid/size-hybm',size(Ver_hyb%m),options_readonly)
      istat= wb_put('model/Vgrid/size-hybt',size(Ver_hyb%t),options_readonly)
      istat= wb_put('model/Vgrid/hybm'     ,Ver_hyb%m      ,options_readwrite)
      istat= wb_put('model/Vgrid/hybt'     ,Ver_hyb%t      ,options_readwrite)
      istat= wb_put('model/Vgrid/am'       ,Ver_a_8%m      ,options_readonly)
      istat= wb_put('model/Vgrid/bm'       ,Ver_b_8%m      ,options_readonly)
      istat= wb_put('model/Vgrid/at'       ,Ver_a_8%t      ,options_readonly)
      istat= wb_put('model/Vgrid/bt'       ,Ver_b_8%t      ,options_readonly)
      istat= wb_put('model/Vgrid/rcoef'    ,Hyb_rcoef      ,options_readonly)
      istat= wb_put('model/Vgrid/vcode'    ,Ver_code       ,options_readonly)

      if (Lun_out > 0) then
         write (Lun_out,1005) G_nk,Hyb_rcoef
            do k=1,G_nk
               call convip(pnip1,Ver_hyb%m(k),3,1,dumc,.false.)
               write (Lun_out,1006) k,Ver_hyb%m(k),Ver_hyb%m(k), &
                                  Ver_hyb%m(k)-Ver_hyb%m(k+1),pnip1
            end do
      end if

 1005 format (/'STAGGERED VERTICAL LAYERING ON',I4,' MOMENTUM HEIGHT LEVELS WITH ', &
               'Grd_rcoef= ',4f7.2,':'/ &
               2x,'level',10x,'HYB',8x,'~HEIGHTS',5x,'~DELTA_Z',7x,'IP1')
 1006 format (1x,i4,3x,es15.5,2(6x,f6.0),4x,i10)
      return
      end
