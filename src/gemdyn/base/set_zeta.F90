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

!  s/r set_zeta    - Generates A and B of the hybrid coordinate
!                    Also sets Z and other related vertical parameters.
!
      subroutine set_zeta (F_hybuser, Nk)
      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use levels
      use ver
      use VERgrid_options
      use vGrid_Descriptors
      use vgrid_wb
      use wb_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Nk
      real, dimension(Nk), intent(in) :: F_hybuser        !user-specified hybrid coordinate values
!
! authors
!      A. Plante & C. Girard - CMC - janvier 2008
!
! object
!    To return A, B parameters for momentum and thermodynamic levels
!    These levels are used for DYNAMICAL CALCULATIONS in the models !
!
!    Also to return other parameters related to the vertical discretization
!
!         Z, dZ, 1/dZ, dBdZ, etc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! n staggered MOMENTUM & THERMO **LEVELS** !!!!!!!!!!!!!!!!!!!!!!!!!
!
!               ttttttttttttttttt Cstv_ztop_8=Ver_z_8%m(0)=Ver_z_8%t(0)=Ver_z_8%x(0)
!
!               - - - - - - - - - Ver_z_8%m(1)
!
!               ================= Ver_z_8%t(1)=Ver_z_8%x(1) = ( Ver_z_8%m(2) + Ver_z_8%m(1) ) / 2
!
!               - - - - - - - - - Ver_z_8%m(2)
!
!                      ...
!
!               - - - - - - - - - Ver_z_8%m(k)
!
!               ================= Ver_x_8%t(k)=Ver_z_8%x(k) = ( Ver_z_8%m(k+1) + Ver_z_8%m(k) )/2
!
!               - - - - - - - - - Ver_z_8%m(k+1)
!
!                      ...
!
!               - - - - - - - - - Ver_z_8%m(n)
!
!               ================= Ver_z_8%t(n)
!
!               sssssssssssssssss Cstv_zsrf_8=Ver_z_8%m(n+1)=Ver_z_8%x(n)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! n staggered MOMENTUM & THERMO **LAYERS** !!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!               Cstv_ztop_8 ttttttttttttttttttttttttttttttttttttttttttt Ver_z_8%x(0)
!                                                     \
!               Ver_z_8%m(1)- - - - - - - - - - - - - - Ver_dz_8%m(1)=Ver_z_8%x(1)-Cstv_ztop_8
!                                          /          /
!   Ver_z_8%m(2)-Ver_z_8%m(1)=Ver_dz_8%t(1)    ======================== Ver_z_8%x(1)
!                                          \          \
!               Ver_z_8%m(2)- - - - - - - - - - - - - - Ver_dz_8%m(2)=Ver_z_8%x(2)-Ver_z_8%x(1)
!                                                     /
!                           =========================================== Ver_z_8%x(2)
!
!                                                     ...
!
!                           =========================================== Ver_z_8%x(k-1)
!                                                     \
!               Ver_z_8%m(k)- - - - - - - - - - - - - - Ver_dz_8%m(k)=Ver_z_8%x(k)-Ver_z_8%x(k-1)
!                                          /          /
! Ver_z_8%m(k+1)-Ver_z_8%m(k)=Ver_dz_8%t(k)    ======================== Ver_z_8%x(k)
!                                          \
!             Ver_z_8%m(k+1)- - - - - - - - - - - - - - - - - - - - - -
!
!                                                     ...
!
!                           =========================================== Ver_z_8%x(n-1)
!                                                   \
!                                                    \
!                                                     \
!                           - - - - - - - - - - - - - - Ver_dz_8%m(n)=Ver_z_8%x(n)-Ver_z_8%x(n-1)
!                                          /          /
!    Cstv_zsrf_8-Ver_z_8%m(n)=Ver_dz_8%t(n)          /
!                                          \        /
!               Cstv_zsrf_8 sssssssssssssssssssssssssssssssssssssssssss Ver_z_8%x(n)
!
! arguments
! none
!

      character(len=32) :: REFP0_S, REFP0_LS_S, dumc
      integer k,istat,pnip1,err,options_readwrite,options_readonly
      integer, dimension(:), pointer :: wkpti
      real, dimension(:), pointer :: std_p_prof=>null(),wkpt
      real    height,heightp1
      real(kind=REAL64)  wk_8
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
      real(kind=REAL64), dimension(:), pointer :: wkpt8
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
                 Ver_wp_8%m(G_nk  ),      Ver_wp_8%t(G_nk  ), &
                 Ver_wm_8%m(G_nk  ),      Ver_wm_8%t(G_nk  ), &
                  Ver_bzz_8(G_nk  ),     Ver_onezero(G_nk+1), &
               Ver_wpstar_8(G_nk  ),    Ver_wmstar_8(G_nk  ), &
                  Ver_wpA_8(G_nk  ),       Ver_wmA_8(G_nk  ), &
                  Ver_wpM_8(G_nk  ),       Ver_wmM_8(G_nk  ), &
                  Ver_wpC_8(G_nk  ),       Ver_wmC_8(G_nk  ), &
                 Ver_czz_8(G_nk  ))

      Cstv_pref_8 = 100000.d0
      Ver_code    = 6

      ! Construct vertical coordinate
      Level_kind_ip1 = 5
      Level_version  = 5
      istat = 0

      Schm_sleve_L= .false. ; err= 0
      if (   Hyb_rcoef(3) >= 0. .or. Hyb_rcoef(4) >= 0. ) then
         if( Hyb_rcoef(3) <  0. .or. Hyb_rcoef(4) <  0. ) err= -1
         Schm_sleve_L= .true. ; Level_version  = 100
      endif
      if (err == 0) then
         if(Schm_sleve_L)then
            istat = vgd_new ( Ver_vgdobj, kind=Level_kind_ip1,&
                         version=Level_version, hyb=F_hybuser,&
                     rcoef1=Hyb_rcoef(1), rcoef2=Hyb_rcoef(2),&
                     rcoef3=Hyb_rcoef(3), rcoef4=Hyb_rcoef(4),&
                         ptop_out_8=wk_8, pref_8=Cstv_pref_8 ,&
                          dhm=0., dht=0., avg_L=.true. )
         else
            istat = vgd_new ( Ver_vgdobj, kind=Level_kind_ip1,&
                         version=Level_version, hyb=F_hybuser,&
                     rcoef1=Hyb_rcoef(1), rcoef2=Hyb_rcoef(2),&
                         ptop_out_8=wk_8, pref_8=Cstv_pref_8 ,&
                          dhm=0., dht=0., avg_L=.true. )
         endif
      endif
      call gem_error (min(err,istat),'SET_ZETA', &
                  'Incorrect vertical construct, check Hyb, Hyb_rcoef')

      Cstv_ptop_8= wk_8

      if (Lun_debug_L) istat = vgd_print(Ver_vgdobj)

      Cstv_Zsrf_8 = log(Cstv_pref_8)
      Cstv_Ztop_8 = log(Cstv_ptop_8)
      Ver_zmin_8 = Cstv_Ztop_8
      Ver_zmax_8 = Cstv_Zsrf_8

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
      if(Schm_sleve_L)then
         if (vgd_get(Ver_vgdobj,'CC_M - vertical C coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
         Ver_c_8%m = wkpt8(1:size(Ver_c_8%m)); deallocate(wkpt8); nullify(wkpt8)
         if (vgd_get(Ver_vgdobj,'CC_T - vertical C coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
         Ver_c_8%t = wkpt8(1:size(Ver_c_8%t)); deallocate(wkpt8); nullify(wkpt8)
      else
         Ver_c_8%m=zero; Ver_c_8%t=zero
      end if
      if (vgd_get(Ver_vgdobj,'VCDM - vertical coordinate (m)'   ,wkpt) /= VGD_OK) istat = VGD_ERROR
      Ver_hyb%m = wkpt(1:size(Ver_hyb%m)); deallocate(wkpt); nullify(wkpt)
      if (vgd_get(Ver_vgdobj,'VCDT - vertical coordinate (t)'   ,wkpt) /= VGD_OK) istat = VGD_ERROR
      Ver_hyb%t = wkpt(1:size(Ver_hyb%t)); deallocate(wkpt); nullify(wkpt)
      if (vgd_get(Ver_vgdobj,'VIPM - level ip1 list (m)'        ,wkpti) /= VGD_OK) istat = VGD_ERROR
      Ver_ip1%m = wkpti(1:size(Ver_ip1%m)); deallocate(wkpti); nullify(wkpti)
      if (vgd_get(Ver_vgdobj,'VIPT - level ip1 list (t)'        ,wkpti) /= VGD_OK) istat = VGD_ERROR
      Ver_ip1%t = wkpti(1:size(Ver_ip1%t)); deallocate(wkpti); nullify(wkpti)

      call handle_error_l(istat==VGD_OK,'set_zeta','retrieving coordinate info')
      if(Schm_sleve_l)then
         istat = vgd_levels(Ver_vgdobj,Ver_ip1%m,std_p_prof,sfc_field=100000.,sfc_field_ls=100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for m levels')
         Ver_std_p_prof%m=std_p_prof
         deallocate(std_p_prof)
         istat = vgd_levels(Ver_vgdobj,Ver_ip1%t,std_p_prof,sfc_field=100000.,sfc_field_ls=100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for t levels')
         Ver_std_p_prof%t=std_p_prof
      else
         istat = vgd_levels(Ver_vgdobj,Ver_ip1%m,std_p_prof,100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for m levels')
         Ver_std_p_prof%m=std_p_prof
         deallocate(std_p_prof)
         istat = vgd_levels(Ver_vgdobj,Ver_ip1%t,std_p_prof,100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for t levels')
         Ver_std_p_prof%t=std_p_prof
      end if

!     -------------------------------
!     Define zeta(m/t/x) from A(m/t):
!     -------------------------------

!     Ver_a_8%m(1:G_nk+1):       G_nk   momentum levels + surface
!     Ver_a_8%t(1:G_nk+1):       G_nk   thermo   levels + surface
!     Ver_z_8%m(0:G_nk+1): top + G_nk   momentum levels + surface
!     Ver_z_8%t(0:G_nk  ): top + G_nk   thermo   levels
!     Ver_z_8%x(0:G_nk  ): top + G_nk-1 thermo   levels + surface

      Ver_z_8%m(0) = Cstv_Ztop_8
      Ver_z_8%t(0) = Cstv_Ztop_8
      do k = 1, G_nk+1
         Ver_z_8%m(k) = Ver_a_8%m(k)
         Ver_z_8%t(k) = Ver_a_8%t(k)
      end do

      if ( Schm_autobar_L ) Ver_z_8%t(G_nk)=Cstv_Zsrf_8

     !Define the positions of zeta_dot
      Ver_z_8%x(0) = Cstv_Ztop_8
      do k = 1, G_nk-1
         Ver_z_8%x(k) = Ver_a_8%t(k)
      end do
      Ver_z_8%x(G_nk)=Cstv_Zsrf_8

!     ----------------------
!     Compute dZ, 1/dZ, dBdZ
!     ----------------------

      do k=1,G_nk
           Ver_dz_8%m(k) = Ver_z_8%x(k) - Ver_z_8%x(k-1)
      end do

      do k=1,G_nk
          Ver_idz_8%m(k) = one/Ver_dz_8%m(k)
           Ver_dz_8%t(k) = Ver_z_8%m(k+1) - Ver_z_8%m(k)
          Ver_idz_8%t(k) = one/Ver_dz_8%t(k)
         Ver_dbdz_8%t(k) = (Ver_b_8%m(k+1)-Ver_b_8%m(k))*Ver_idz_8%t(k)
         Ver_dcdz_8%t(k) = (Ver_c_8%m(k+1)-Ver_c_8%m(k))*Ver_idz_8%t(k)
         if(k == 1) then
            Ver_dbdz_8%m(k) = (Ver_b_8%t(k)-zero)*Ver_idz_8%m(k)
            Ver_dcdz_8%m(k) = (Ver_c_8%t(k)-zero)*Ver_idz_8%m(k)
         elseif(k == G_nk) then
            Ver_dbdz_8%m(k) = (one -Ver_b_8%t(k-1))*Ver_idz_8%m(k)
            Ver_dcdz_8%m(k) = (zero-Ver_c_8%t(k-1))*Ver_idz_8%m(k)
         else
            Ver_dbdz_8%m(k) = (Ver_b_8%t(k)-Ver_b_8%t(k-1))*Ver_idz_8%m(k)
            Ver_dcdz_8%m(k) = (Ver_c_8%t(k)-Ver_c_8%t(k-1))*Ver_idz_8%m(k)
         end if
      end do

      if(Schm_autobar_L) then
         do k=1,G_nk
            Ver_dbdz_8%m(k) = one/(Cstv_Zsrf_8-Ver_z_8%m(1))
            Ver_dcdz_8%m(k) = one/(Cstv_Zsrf_8-Ver_z_8%m(1))
         end do
      end if

!     -------------------------------------------------------
!     Compute AVERGING WEIGHTS FROM THERMO TO MOMENTUM LEVELS
!     -------------------------------------------------------

      do k=1,G_nk
         Ver_wmM_8(k) = (Ver_z_8%t(k)-Ver_z_8%m(k))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
         Ver_wpM_8(k) = one - Ver_wmM_8(k)

         Ver_wmC_8(k) = (Ver_z_8%x(k)-Ver_z_8%m(k))/(Ver_z_8%x(k)-Ver_z_8%x(k-1))
         Ver_wpC_8(k) = one - Ver_wmC_8(k)

         Ver_wmA_8(k) = (Ver_z_8%m(k)-Ver_z_8%x(k-1))/(Ver_z_8%x(k)-Ver_z_8%x(k-1))
         Ver_wpA_8(k) = one - Ver_wmA_8(k)

         Ver_wm_8%m(k)= Ver_wmA_8(k)
         Ver_wp_8%m(k)= Ver_wpA_8(k)
      end do

!     -------------------------------------------------------
!     Compute AVERGING WEIGHTS FROM MOMENTUM TO THERMO LEVELS
!     -------------------------------------------------------

      do k=1,G_nk-1
         Ver_wp_8%t(k) = half
         Ver_wm_8%t(k) = half
      end do
      Ver_wp_8%t(G_nk) = one
      Ver_wm_8%t(G_nk) = zero
!
!     SPECIAL WEIGHTS due to last thermo level
!
      Ver_wmstar_8 = zero
      Ver_wpstar_8 = one

      if ( .not.Schm_autobar_L ) then
         Ver_wmstar_8(G_nk)=half*Ver_dz_8%t(G_nk)/Ver_dz_8%m(G_nk)
         Ver_wpstar_8(G_nk)=one-Ver_wmstar_8(G_nk)
         Ver_wp_8%m(G_nk) = Ver_wpstar_8(G_nk) * Ver_wpA_8(G_nk)
         Ver_wm_8%m(G_nk) = one - Ver_wp_8%m(G_nk)
      end if
!     -------------------------------------------------------
!     Compute VERTICAL AVERGING OF B from thermo to momentum
!     -------------------------------------------------------

      do k=1,G_nk
         if(k == 1) then
            Ver_bzz_8(k) = Ver_wp_8%m(k)*Ver_b_8%t(k) &
                         + Ver_wm_8%m(k)*zero
            Ver_czz_8(k) = Ver_wp_8%m(k)*Ver_c_8%t(k) &
                         + Ver_wm_8%m(k)*zero
         elseif(k == G_nk) then
            Ver_bzz_8(k) = Ver_wp_8%m(k)*one &
                         + Ver_wm_8%m(k)*Ver_b_8%t(k-1)
            Ver_czz_8(k) = Ver_wp_8%m(k)*zero &
                         + Ver_wm_8%m(k)*Ver_c_8%t(k-1)
         else
            Ver_bzz_8(k) = Ver_wp_8%m(k)*Ver_b_8%t(k) &
                         + Ver_wm_8%m(k)*Ver_b_8%t(k-1)
            Ver_czz_8(k) = Ver_wp_8%m(k)*Ver_c_8%t(k) &
                         + Ver_wm_8%m(k)*Ver_c_8%t(k-1)
         end if
      end do

!     -------------------------------------------------------
!     Initialize Ver_onezero
!     -------------------------------------------------------

      Ver_onezero=1.
      Ver_onezero(1)=0.

      REFP0_S = 'PW_P0:P'  !# gmmk_pw_p0_plus_s !NOTE: could gmmk_* be defined as parameters in a .cdk, this way it could be used here and would be more consistent
      REFP0_LS_S = ' '
      if (Schm_sleve_L) REFP0_LS_S = 'PW_P0_LS'  !# gmmk_pw_p0_ls_s
      istat = vgrid_wb_put('ref-m', Ver_vgdobj, Ver_ip1%m,  &
                           REFP0_S, REFP0_LS_S, F_overwrite_L=.true.)
      istat = vgrid_wb_put('ref-t', Ver_vgdobj, Ver_ip1%t,  &
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
      istat= wb_put('model/Vgrid/ptop'     ,Cstv_ptop_8    ,options_readonly)
      istat= wb_put('model/Vgrid/pref'     ,Cstv_pref_8    ,options_readonly)
      istat= wb_put('model/Vgrid/vcode'    ,Ver_code       ,options_readonly)

      if (Lun_out > 0) then
         write (Lun_out,1005) G_nk,Hyb_rcoef
         do k=1,G_nk
            height = -16000./log(10.)*log(Ver_hyb%m(k))

            if (k < G_nk) then
               heightp1 = -16000./log(10.)*log(Ver_hyb%m(k+1))
            end if

            if (k == G_nk) then
               heightp1 = 0.
            end if

            call convip(pnip1,Ver_hyb%m(k),5,1,dumc,.false.)
            write (Lun_out,1006) k,Ver_hyb%m(k),height, &
                                 height-heightp1,pnip1
         end do
      end if

      if (Lun_debug_L) call prgenab()

 1005 format (/'STAGGERED VERTICAL LAYERING ON',I4,' MOMENTUM HYBRID LEVELS WITH ',&
               'Grd_rcoef= ',4f7.2/2x,'level',10x,'HYB',8x,'~HEIGHTS',5x,'~DELTA_Z',7x,'IP1')
 1006 format (1x,i4,3x,es15.5,2(6x,f6.0),4x,i10)
!
!     __________________________________________________________________
!
      return
      end
