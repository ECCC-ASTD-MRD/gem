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

      subroutine adz_set
      use HORgrid_options
      use ctrl
      use glb_ld
      use adz_options
      use adz_mem
      use lam_options
      use ver
      use gmm_pw
      use gmm_itf_mod
      use rstr
      implicit none

#include "gmm_gem_flags.hf"
#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      integer  j, k, k0, BCS_BASE, pnz, ext
      real*8, parameter :: EPS_8= 1.D-5
      real*8 :: ra,rb,rc,rd
      real*16 :: smallest_dz,posz

      type(gmm_metadata) :: meta, mymeta
      integer*8 :: flag_m_f
      integer :: kp1, flag_r_n, istat,istatu,istatv
!
!---------------------------------------------------------------------
!
      Adz_maxcfl= max(1,Grd_maxcfl)
      Adz_halox = Adz_maxcfl + 1
      Adz_haloy = Adz_halox

      Adz_lminx = 1    - Adz_halox
      Adz_lmaxx = l_ni + Adz_halox
      Adz_lminy = 1    - Adz_haloy
      Adz_lmaxy = l_nj + Adz_haloy

      Adz_nit = Adz_lmaxx - Adz_lminx + 1
      Adz_njt = Adz_lmaxy - Adz_lminy + 1
      Adz_nij = Adz_nit * Adz_njt
      Adz_2dnh = l_ni * l_nj
      Adz_3dnh = Adz_2dnh * l_nk

      Adz_ioff = l_i0 - 2
      Adz_joff = l_j0 - 2

      ext=1
      if (Grd_yinyang_L) ext=2

      Adz_i0=   1  + pil_w - (ext+1)*west
      Adz_in= l_ni - pil_e +  ext   *east
      Adz_j0=   1  + pil_s - (ext+1)*south
      Adz_jn= l_nj - pil_n +  ext   *north

      Adz_i0u = Adz_i0 ;  Adz_inu = Adz_in
      Adz_j0v = Adz_j0 ;  Adz_jnv = Adz_jn
      if (.not. Grd_yinyang_L) then
         Adz_i0u= 1 + pil_w     ; Adz_j0v= 1 + pil_s
         Adz_inu= l_niu - pil_e ; Adz_jnv= l_njv - pil_n
      end if

      Adz_k0 = Lam_gbpil_t+1
      Adz_k0t= Adz_k0
      if(Lam_gbpil_t > 0) Adz_k0t= Adz_k0 - 1
      Adz_k0m= max(Adz_k0t-2,1)

      Adz_kkmax= l_nk-2

      BCS_BASE= 4
      if (Grd_yinyang_L) BCS_BASE = 3

      Adz_iminposx = l_i0+adz_lminx   + EPS_8
      Adz_imaxposx = l_i0+adz_lmaxx-2 - EPS_8
      Adz_iminposy = l_j0+adz_lminy   + EPS_8
      Adz_imaxposy = l_j0+adz_lmaxy-2 - EPS_8
      if (l_west ) Adz_iminposx = l_i0+     BCS_BASE   + EPS_8
      if (l_east ) Adz_imaxposx = l_i0+l_ni-BCS_BASE-1 - EPS_8
      if (l_south) Adz_iminposy = l_j0+     BCS_BASE   + EPS_8
      if (l_north) Adz_imaxposy = l_j0+l_nj-BCS_BASE-1 - EPS_8

      allocate (  Adz_delz_m(0:l_nk),  Adz_delz_t(0:l_nk), &
                 Adz_odelz_m(0:l_nk), Adz_odelz_t(0:l_nk) )
      do k = 0,l_nk
         kp1 = k+1
         Adz_delz_m (k) = Ver_z_8%m(kp1)-Ver_z_8%m(k)
         Adz_delz_t (k) = Ver_z_8%t(kp1)-Ver_z_8%t(k)
         Adz_odelz_m(k) = 1.0d0 / Adz_delz_m(k)
         Adz_odelz_t(k) = 1.0d0 / Adz_delz_t(k)
      end do

      smallest_dz= minval(Adz_delz_m(1:l_nk))
      adz_ovdzm_8 = 1.0d0/smallest_dz
      pnz = nint(1.0+(Ver_z_8%m(l_nk+1)-Ver_z_8%m(0))*adz_ovdzm_8)
      allocate ( Adz_search_m(pnz) )
      k0 = 0
      do k= 1, pnz
         posz = Ver_z_8%m(0) + dble(k-1) * smallest_dz
         if (posz > Ver_z_8%m(k0+1)) k0 = min(l_nk+1, k0+1)
         Adz_search_m(k) = k0
      end do

      smallest_dz= minval(Adz_delz_t(1:l_nk))
      adz_ovdzt_8 = 1.0d0/smallest_dz
      pnz = nint(1.0+(Ver_z_8%t(l_nk)-Ver_z_8%t(0))*adz_ovdzt_8)
      allocate ( Adz_search_t(pnz) )
      k0 = 0
      do k= 1, pnz
         posz = Ver_z_8%t(0) + dble(k-1) * smallest_dz
         if (posz > Ver_z_8%t(k0+1)) k0 = min(l_nk+1, k0+1)
         Adz_search_t(k) = k0
      end do

      allocate ( Adz_cy_8(Adz_lminy:Adz_lmaxy) )

      do j = Adz_lminy, Adz_lmaxy
         Adz_cy_8(j) = 1.d0 / cos(G_yg_8(l_j0+j-1))
      end do

      allocate( Adz_zabcd_8%t(l_nk), Adz_zabcd_8%m(l_nk),&
                Adz_zbacd_8%t(l_nk), Adz_zbacd_8%m(l_nk),&
                Adz_zcabd_8%t(l_nk), Adz_zcabd_8%m(l_nk),&
                Adz_zdabc_8%t(l_nk), Adz_zdabc_8%m(l_nk) )

      do k= 2, l_nk-2
         ra = Ver_z_8%m(k-1)
         rb = Ver_z_8%m(k)
         rc = Ver_z_8%m(k+1)
         rd = Ver_z_8%m(k+2)
         Adz_zabcd_8%m(k) = 1.0/triprod(ra,rb,rc,rd)
         Adz_zbacd_8%m(k) = 1.0/triprod(rb,ra,rc,rd)
         Adz_zcabd_8%m(k) = 1.0/triprod(rc,ra,rb,rd)
         Adz_zdabc_8%m(k) = 1.0/triprod(rd,ra,rb,rc)

         ra = Ver_z_8%t(k-1)
         rb = Ver_z_8%t(k)
         rc = Ver_z_8%t(k+1)
         rd = Ver_z_8%t(k+2)
         Adz_zabcd_8%t(k) = 1.0/triprod(ra,rb,rc,rd)
         Adz_zbacd_8%t(k) = 1.0/triprod(rb,ra,rc,rd)
         Adz_zcabd_8%t(k) = 1.0/triprod(rc,ra,rb,rd)
         Adz_zdabc_8%t(k) = 1.0/triprod(rd,ra,rb,rc)
      end do

      allocate (&
         Adz_uu_ext( Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_vv_ext( Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_ww_ext( Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_uvw_d(3,Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_uu_arr(l_ni,l_nj,l_nk), Adz_uu_dep(l_ni,l_nj,l_nk),&
         Adz_vv_arr(l_ni,l_nj,l_nk), Adz_vv_dep(l_ni,l_nj,l_nk),&
         Adz_ww_arr(l_ni,l_nj,l_nk), Adz_ww_dep(l_ni,l_nj,l_nk) )

      Adz_uu_ext=0. ; Adz_vv_ext=0. ; Adz_ww_ext=0.

      allocate (Adz_pxyzmu (3,l_ni,l_nj,l_nk), &
                Adz_pxyzmv (3,l_ni,l_nj,l_nk), &
                Adz_pxyzt  (3,l_ni,l_nj,l_nk))
      Adz_pxyzmu= 0. ; Adz_pxyzmv= 0. ; Adz_pxyzt= 0.

      allocate ( Adz_wpxyz(-1:l_ni+2,-1:l_nj+2,l_nk,3) )
      Adz_wpxyz(-1:0,:,:,:)=0. ; Adz_wpxyz(l_ni+1:l_ni+2,:,:,:)=0.
      Adz_wpxyz(:,-1:0,:,:)=0. ; Adz_wpxyz(:,l_nj+1:l_nj+2,:,:)=0.

      flag_m_f = FLAG_LVL_M
      flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

      call gmm_build_meta4D (meta,  1,3   ,0,0,   3, &
                                    1,l_ni,0,0,l_ni, &
                                    1,l_nj,0,0,l_nj, &
                                    1,l_nk,0,0,l_nk, &
                                    0,GMM_NULL_FLAGS)
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)

      istat= gmm_create('ADZ_PXYZM',Adz_pxyzm, mymeta, flag_r_n)
      istat= gmm_get('ADZ_PXYZM',Adz_pxyzm)

!     Wind information at the lowest thermodynamic level
!     from the physics surface layer scheme
      nullify(Adz_uslt,Adz_vslt)
      if ( .not. Ctrl_phyms_L ) Adz_slt_winds= .false.
      if (Adz_slt_winds) then
         istatu= gmm_get(gmmk_pw_uslt_s, Adz_uslt)
         istatv= gmm_get(gmmk_pw_vslt_s, Adz_vslt)
         if ((istatu/=0).or.(istatv/=0)) Adz_slt_winds= .false.
      end if

      Adz_niter = Adz_itraj
      if ( .not. Rstri_rstn_L ) call adz_inittraj
!
!---------------------------------------------------------------------
!
      return
contains

      real*8 function triprod(za,zb,zc,zd)
         real*8, intent(in) :: za,zb,zc,zd
         triprod = ((za-zb)*(za-zc)*(za-zd))
      end function triprod

      end subroutine adz_set
