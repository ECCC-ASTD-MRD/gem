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

      subroutine adz_interp_traj
      use adz_mem
      use adz_options
      use adv_pos
      use glb_ld
      use cstv
      use dcst
      use gmm_vt0
      use ver
      use geomh
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer :: i,j,k,k00,ii1,jj1
      real(kind=REAL64), dimension(-1:l_ni+2,-1:l_nj+2,l_nk) :: xm,ym,zm
      real(kind=REAL64) :: aa,bb,cc, xt,yt,zt, wdt,ww,wp,wm
      real(kind=REAL64) :: hh, x1, x2, x3, x4
      real(kind=REAL64), dimension(2:l_nk-2) :: w1, w2, w3, w4
      real(kind=REAL64) :: zmin_bound, zmax_bound
      real(kind=REAL64), parameter :: half=0.5d0
!
!---------------------------------------------------------------------
!
      zmin_bound = dble(0)
      zmax_bound = dble(l_nk+1)

      aa = -0.0625d0
      bb = +0.5625d0
      cc = geomh_hx_8 * 0.5d0

      call rpn_comm_xch_halo_8 (Adz_wpxyz, -1,l_ni+2, -1,l_nj+2,&
                 l_ni,l_nj, 3*l_nk, 2,2, .false.,.false., l_ni,0)

      do k= Adz_k0, l_nk
         do j= Adz_j0, Adz_jn
            do i= Adz_i0u, Adz_inu
               xm(i,j,k) =  aa * (Adz_wpxyz(i-1,j,k,1) + Adz_wpxyz(i+2,j,k,1)) &
                          + bb * (Adz_wpxyz(i  ,j,k,1) + Adz_wpxyz(i+1,j,k,1)) - 0.5d0
               ym(i,j,k) =  aa * (Adz_wpxyz(i-1,j,k,2) + Adz_wpxyz(i+2,j,k,2)) &
                          + bb * (Adz_wpxyz(i  ,j,k,2) + Adz_wpxyz(i+1,j,k,2))
               zt        =  aa * (Adz_wpxyz(i-1,j,k,3) + Adz_wpxyz(i+2,j,k,3)) &
                          + bb * (Adz_wpxyz(i  ,j,k,3) + Adz_wpxyz(i+1,j,k,3))
               Adz_pmu(1,i,j,k)= min(max(xm(i,j,k),Adz_iminposx),Adz_imaxposx)
               Adz_pmu(2,i,j,k)= min(max(ym(i,j,k),Adz_iminposy),Adz_imaxposy)
               Adz_pmu(3,i,j,k)= min(max(zt       ,zmin_bound  ),zmax_bound  )
               zm(i,j,k) = Adz_pmu(3,i,j,k)
            end do
         end do
      end do

      do k= Adz_k0, l_nk
         do j= Adz_j0v, Adz_jnv
            do i= Adz_i0, Adz_in
               xm(i,j,k) =  aa * (Adz_wpxyz(i,j-1,k,1) + Adz_wpxyz(i,j+2,k,1)) &
                          + bb * (Adz_wpxyz(i,j  ,k,1) + Adz_wpxyz(i,j+1,k,1))
               ym(i,j,k) =  aa * (Adz_wpxyz(i,j-1,k,2) + Adz_wpxyz(i,j+2,k,2)) &
                          + bb * (Adz_wpxyz(i,j  ,k,2) + Adz_wpxyz(i,j+1,k,2)) - 0.5d0
               zt =         aa * (Adz_wpxyz(i,j-1,k,3) + Adz_wpxyz(i,j+2,k,3)) &
                          + bb * (Adz_wpxyz(i,j  ,k,3) + Adz_wpxyz(i,j+1,k,3))
               Adz_pmv(1,i,j,k)= min(max(xm(i,j,k),Adz_iminposx),Adz_imaxposx)
               Adz_pmv(2,i,j,k)= min(max(ym(i,j,k),Adz_iminposy),Adz_imaxposy)
               Adz_pmv(3,i,j,k)= min(max(zt       ,zmin_bound  ),zmax_bound  )
               zm(i,j,k) = Adz_pmv(3,i,j,k)
           end do
         end do
      end do

! Prepare parameters for cubic vertical intepolation
      do k=2,l_nk-2
         hh = Ver_z_8%t(k)
         x1 = Ver_z_8%m(k-1)
         x2 = Ver_z_8%m(k  )
         x3 = Ver_z_8%m(k+1)
         x4 = Ver_z_8%m(k+2)
         w1(k) = lagcoef( hh, x1, x2, x3, x4 )
         w2(k) = lagcoef( hh, x2, x1, x3, x4 )
         w3(k) = lagcoef( hh, x3, x1, x2, x4 )
         w4(k) = lagcoef( hh, x4, x1, x2, x3 )
      end do

      k00=max(Adz_k0t-1,1)
      if (.not.associated(pxt)) then
         allocate (pxt(l_ni,l_nj,l_nk),pyt(l_ni,l_nj,l_nk),pzt(l_ni,l_nj,l_nk))
      endif

      do k= max(k00,2), l_nk-2
         do j= 1, l_nj
!!DIR$ IVDEP
            do i= 1, l_ni
               xt = w1(k)*Adz_wpxyz(i,j,k-1,1)+ &
                    w2(k)*Adz_wpxyz(i,j,k  ,1)+ &
                    w3(k)*Adz_wpxyz(i,j,k+1,1)+ &
                    w4(k)*Adz_wpxyz(i,j,k+2,1)
               yt = w1(k)*Adz_wpxyz(i,j,k-1,2)+ &
                    w2(k)*Adz_wpxyz(i,j,k  ,2)+ &
                    w3(k)*Adz_wpxyz(i,j,k+1,2)+ &
                    w4(k)*Adz_wpxyz(i,j,k+2,2)
               wdt= w1(k)*Adz_uvw_dep(3,i,j,k-1)+ &
                    w2(k)*Adz_uvw_dep(3,i,j,k  )+ &
                    w3(k)*Adz_uvw_dep(3,i,j,k+1)+ &
                    w4(k)*Adz_uvw_dep(3,i,j,k+2)
               zt = Ver_z_8%t(k) - Cstv_dtzD_8*  wdt &
                                 - Cstv_dtzA_8* zdt0(i,j,k)
               xm(i,j,k)=xt; ym(i,j,k)=yt
               xt = min(max(xt,Adz_iminposx),Adz_imaxposx)
               yt = min(max(yt,Adz_iminposy),Adz_imaxposy)
               Adz_pxyzt(1,i,j,k)= xt
               Adz_pxyzt(2,i,j,k)= yt
               Adz_pxyzt(3,i,j,k)= zindx(zt)
               zm (i,j,k)= Adz_pxyzt(3,i,j,k)
               pzt(i,j,k)= zt
               ii1 = xt ; jj1 = yt
               pxt(i,j,k) = G_xg_8(ii1) + (xt-dble(ii1))*geomh_hx_8
               pyt(i,j,k) = G_yg_8(jj1) + (yt-dble(jj1))*geomh_hy_8
            end do
         end do
      end do

      ww=Ver_wmstar_8(l_nk)
      wp=(Ver_z_8%t(l_nk  )-Ver_z_8%m(l_nk-1))*Ver_idz_8%t(l_nk-1)
      wm=1.d0-wp

      k= l_nk-1
      do j= 1, l_nj
!!DIR$ IVDEP
         do i= 1, l_ni
            xt = (Adz_wpxyz(i,j,k,1)+Adz_wpxyz(i,j,k+1,1))*half
            yt = (Adz_wpxyz(i,j,k,2)+Adz_wpxyz(i,j,k+1,2))*half
            wdt= (Adz_uvw_dep(3,i,j,k)+Adz_uvw_dep(3,i,j,k+1))*half
            zt =Ver_z_8%t(k) - Cstv_dtzD_8*  wdt &
                             - Cstv_dtzA_8*zdt0(i,j,k)
            xm(i,j,k)=xt; ym(i,j,k)=yt
            xt = min(max(xt,Adz_iminposx),Adz_imaxposx)
            yt = min(max(yt,Adz_iminposy),Adz_imaxposy)
            Adz_pxyzt(1,i,j,k)= xt
            Adz_pxyzt(2,i,j,k)= yt
            Adz_pxyzt(3,i,j,k)= zindx(zt)
            zm (i,j,k)= Adz_pxyzt(3,i,j,k)
            pzt(i,j,k)= zt
            ii1 = xt ; jj1 = yt
            pxt(i,j,k) = G_xg_8(ii1) + (xt-dble(ii1))*geomh_hx_8
            pyt(i,j,k) = G_yg_8(jj1) + (yt-dble(jj1))*geomh_hy_8
         end do
      end do

      k= l_nk
      if (Adz_slt_winds) then
         do j= 1, l_nj
!!DIR$ IVDEP
         do i= 1, l_ni
            xt= (geomh_x_8(i) - Dcst_inv_rayt_8 * Cstv_dt_8 * &
                     Adz_uslt(i,j)*Adz_cy_8(j) - G_xg_8(1)) * &
                     geomh_inv_hx_8 + 1.d0
            yt= (geomh_y_8(j) - Dcst_inv_rayt_8 * Cstv_dt_8 * &
                 Adz_vslt(i,j) - G_yg_8(1)) * geomh_inv_hy_8 + 1.d0
            xm(i,j,k)=xt; ym(i,j,k)=yt
            xt = min(max(xt,Adz_iminposx),Adz_imaxposx)
            yt = min(max(yt,Adz_iminposy),Adz_imaxposy)
            Adz_pxyzt(1,i,j,k)= xt
            Adz_pxyzt(2,i,j,k)= yt
            ii1 = xt ; jj1 = yt
            pxt(i,j,k) = G_xg_8(ii1) + (xt-dble(ii1))*geomh_hx_8
            pyt(i,j,k) = G_yg_8(jj1) + (yt-dble(jj1))*geomh_hy_8
         end do
         end do
      else
         do j= 1, l_nj
!!DIR$ IVDEP
         do i= 1, l_ni
            xt = wp*Adz_wpxyz(i,j,k,1)+wm*Adz_wpxyz(i,j,k-1,1)
            yt = wp*Adz_wpxyz(i,j,k,2)+wm*Adz_wpxyz(i,j,k-1,2)
            xm(i,j,k)=xt; ym(i,j,k)=yt
            xt = min(max(xt,Adz_iminposx),Adz_imaxposx)
            yt = min(max(yt,Adz_iminposy),Adz_imaxposy)
            Adz_pxyzt(1,i,j,k)= xt
            Adz_pxyzt(2,i,j,k)= yt
            ii1 = xt ; jj1 = yt
            pxt(i,j,k) = G_xg_8(ii1) + (xt-dble(ii1))*geomh_hx_8
            pyt(i,j,k) = G_yg_8(jj1) + (yt-dble(jj1))*geomh_hy_8
         end do
         end do
      end if

      do j= 1, l_nj
!!DIR$ IVDEP
      do i= 1, l_ni
         wdt= (Adz_uvw_dep(3,i,j,k-1)+Adz_uvw_dep(3,i,j,k))*half
         zt = Ver_z_8%t(k)-ww*(Cstv_dtD_8*  wdt &
             +Cstv_dtA_8*zdt0(i,j,k-1))
         Adz_pxyzt(3,i,j,k)= zindx(zt)
         pzt(i,j,k)= zt
         zm (i,j,k)= Adz_pxyzt(3,i,j,k)
      end do
      end do

      if (k00==1) then
         do j= 1, l_nj
!!DIR$ IVDEP
            do i= 1, l_ni
               xt = (Adz_wpxyz(i,j,1,1 )+Adz_wpxyz (i,j,2,1))*half
               yt = (Adz_wpxyz(i,j,1,2 )+Adz_wpxyz (i,j,2,2))*half
               wdt= (Adz_uvw_dep(3,i,j,1)+Adz_uvw_dep(3,i,j,2))*half
               zt = Ver_z_8%t(1) - Cstv_dtzD_8*wdt &
                                 - Cstv_dtzA_8*zdt0(i,j,1)
               xm(i,j,1)=xt; ym(i,j,1)=yt
               xt = min(max(xt,Adz_iminposx),Adz_imaxposx)
               yt = min(max(yt,Adz_iminposy),Adz_imaxposy)
               Adz_pxyzt(1,i,j,1)= xt
               Adz_pxyzt(2,i,j,1)= yt
               Adz_pxyzt(3,i,j,1)= zindx(zt)
               zm (i,j,1)= Adz_pxyzt(3,i,j,1)
               pzt(i,j,1)= zt
               ii1 = xt ; jj1 = yt
               pxt(i,j,1) = G_xg_8(ii1) + (xt-dble(ii1))*geomh_hx_8
               pyt(i,j,1) = G_yg_8(jj1) + (yt-dble(jj1))*geomh_hy_8
            end do
         end do
      end if

      Adz_pt   (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, Adz_k0t:l_nk)=&
      Adz_pxyzt(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, Adz_k0t:l_nk)

      Adz_pb   (:,Adz_i0b:Adz_inb, Adz_j0b:Adz_jnb, Adz_k0t:l_nk)=&
      Adz_pxyzt(:,Adz_i0b:Adz_inb, Adz_j0b:Adz_jnb, Adz_k0t:l_nk)
!
!---------------------------------------------------------------------
!
      return

contains

      real(kind=REAL64) function lagcoef( x, x1, x2, x3, x4 )
      implicit none
      real(kind=REAL64), intent(in) :: x, x1, x2, x3, x4
      lagcoef = ( ( x  - x2 ) * ( x  - x3 ) * ( x  - x4 ) )/ &
                ( ( x1 - x2 ) * ( x1 - x3 ) * ( x1 - x4 ) )
      return
      end function lagcoef

      real(kind=REAL64) function zindx (posz)
      use, intrinsic :: iso_fortran_env
      implicit none
      real(kind=REAL64), intent(inout) :: posz
      integer :: kk1,nb
      posz = min(max(posz,Ver_zmin_8),Ver_zmax_8)
      kk1 = (posz - ver_z_8%t(0)  ) * adz_ovdzt_8 + 1.d0
      kk1 = min(max(1,kk1),ubound(adz_search_t,1))
      kk1 = adz_search_t(kk1)
      if ( sig*posz >  sig*ver_z_8%t(min(kk1+1,l_nk+1))) kk1= kk1 + 1
      if ( sig*posz <  sig*ver_z_8%t(kk1)              ) kk1= kk1 - 1
      kk1 = min(l_nk+1,max(0,kk1))
      nb  = max(min(kk1,G_nk-1),1)
      zindx= (posz-ver_z_8%t(nb))*Adz_odelz_t(nb) + dble(nb)
      return
      end function zindx

      end subroutine adz_interp_traj

