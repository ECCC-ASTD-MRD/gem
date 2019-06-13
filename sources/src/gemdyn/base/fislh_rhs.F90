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
!
!*s/r rhs_H - compute the right-hand sides: Ru, Rv, Rc, Rt, Rw,
!             and save the results for next iteration in the o's
!           - Height-type vertical coordinate
!
      subroutine fislh_rhs ( F_oru, F_orv, F_ort, F_orw, F_orc, F_orf, &
                             F_u  ,   F_v,   F_t,   F_w,  F_zd,   F_q, &
                             Minx,Maxx,Miny,Maxy, Nk )
      use HORgrid_options
      use gem_options
      use dyn_fisl_options
      use coriolis
      use geomh
      use tdpack
      use gmm_itf_mod
      use gmm_phy
      use glb_ld
      use cstv
      use ver
      use lun
      use metric
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_oru,F_orv,F_ort
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_orw,F_orc,F_orf
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_u,F_v,F_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)    :: F_w,F_zd
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(inout) :: F_q

!Author: Claude Girard, July 2017

      integer :: i0,  in,  j0,  jn
      integer :: i0u, inu, i0v, inv
      integer :: j0u, jnu, j0v, jnv

      integer :: i, j, k, km, kp, nij, jext, istat
      real(kind=REAL64)  :: div, barz, barzp, u_interp, v_interp, t_interp, w1
      real(kind=REAL64)  :: phy_bA_m_8, phy_bA_t_8
      real(kind=REAL64), dimension(l_ni,l_nj) :: xtmp_8, ytmp_8
      real(kind=REAL64), parameter :: one=1.d0, zero=0.d0, half=0.5d0
      real, dimension(Minx:Maxx,Miny:Maxy,l_nk), target :: zero_array
      zero_array=0.
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      if (Schm_phycpl_S == 'RHS') then
         istat = gmm_get(gmmk_phy_uu_tend_s,phy_uu_tend)
         istat = gmm_get(gmmk_phy_vv_tend_s,phy_vv_tend)
         istat = gmm_get(gmmk_phy_tv_tend_s,phy_tv_tend)
         phy_bA_m_8 = 0.d0
         phy_bA_t_8 = 0.d0
      else if (Schm_phycpl_S == 'AVG') then
         istat = gmm_get(gmmk_phy_uu_tend_s,phy_uu_tend)
         istat = gmm_get(gmmk_phy_vv_tend_s,phy_vv_tend)
         istat = gmm_get(gmmk_phy_tv_tend_s,phy_tv_tend)
         phy_bA_m_8 = Cstv_bA_m_8
         phy_bA_t_8 = Cstv_bA_8
      else
         phy_uu_tend => zero_array
         phy_vv_tend => zero_array
         phy_tv_tend => zero_array
         phy_bA_m_8 = 0.d0
         phy_bA_t_8 = 0.d0
      end if

!     Common coefficients

      jext = Grd_bsc_ext1 + 1

!     Exchanging halos for derivatives & interpolation
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call rpn_comm_xch_halo( F_u , l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj ,G_nk  , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( F_v , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_njv,G_nk  , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( F_t , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,G_nk  , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( F_q,  l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,G_nk+1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

      nij = l_ni*l_nj

!     Indices to compute Rc, Rt, Rw
      i0 = 1
      j0 = 1
      in = l_ni
      jn = l_nj
      if (l_west)  i0 = 4+jext
      if (l_east)  in = l_niu-2-jext
      if (l_south) j0 = 4+jext
      if (l_north) jn = l_njv-2-jext

!     Additional indices to compute Ru, Rv
      i0u = 1
      i0v = 1
      inu = l_niu
      inv = l_ni
      j0u = 1
      j0v = 1
      jnu = l_nj
      jnv = l_njv

      if (l_west ) i0u = 2+jext
      if (l_west ) i0v = 3+jext
      if (l_east ) inu = l_niu-1-jext
      if (l_east ) inv = l_niu-1-jext
      if (l_south) j0u = 3+jext
      if (l_south) j0v = 2+jext
      if (l_north) jnu = l_njv-1-jext
      if (l_north) jnv = l_njv-1-jext

!$omp parallel private(km,kp,barz,barzp, &
!$omp     u_interp,v_interp,t_interp, &
!$omp     div,xtmp_8,ytmp_8,w1)

!$omp do
      do k = 1,l_nk
         km=max(k-1,1)
         kp=min(k+1,l_nk)

         !********************************
         ! Compute Ru: RHS of U equation *
         !********************************
         do j= j0u, jnu
            do i= i0u, inu

               !assuming T at momentum level 1 equal to T at thermo level 3/2
               barz  = Ver_wp_8%m(k)*F_t(i  ,j,k)+Ver_wm_8%m(k)*F_t(i  ,j,km)
               barzp = Ver_wp_8%m(k)*F_t(i+1,j,k)+Ver_wm_8%m(k)*F_t(i+1,j,km)
               t_interp = (barz + barzp)*half/Cstv_Tstr_8

               v_interp = 0.25d0*(F_v(i,j,k)+F_v(i,j-1,k)+F_v(i+1,j,k)+F_v(i+1,j-1,k))

               F_oru(i,j,k) = Cstv_invT_m_8  * F_u(i,j,k) - Cstv_Beta_m_8 * ( &
                              t_interp * ( F_q(i+1,j,k) - F_q(i,j,k) ) * geomh_invDX_8(j)          &
                            - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * F_u(i,j,k) ) * v_interp )  &
                            + Cstv_Beta_m_8 * t_interp * mc_Jx(i,j,k) * ( &
                              Ver_wp_8%m(k)*half*( (F_q(i+1,j,k+1)-F_q(i+1,j,k ))*mc_iJz(i+1,j,k )   &
                                                 + (F_q(i  ,j,k+1)-F_q(i  ,j,k ))*mc_iJz(i  ,j,k ) ) &
                            + Ver_wm_8%m(k)*half*( (F_q(i+1,j,k  )-F_q(i+1,j,km))*mc_iJz(i+1,j,km)   &
                                                 + (F_q(i  ,j,k  )-F_q(i  ,j,km))*mc_iJz(i  ,j,km) ) )
               F_oru(i,j,k) = F_oru(i,j,k) + ((1d0-phy_bA_m_8)/Cstv_bA_m_8) * phy_uu_tend(i,j,k)
            end do
         end do

         !********************************
         ! Compute Rv: RHS of V equation *
         !********************************
         do j = j0v, jnv
            do i = i0v, inv

               barz  = Ver_wp_8%m(k)*F_t(i,j  ,k)+Ver_wm_8%m(k)*F_t(i,j  ,km)
               barzp = Ver_wp_8%m(k)*F_t(i,j+1,k)+Ver_wm_8%m(k)*F_t(i,j+1,km)
               t_interp = ( barz + barzp)*half/Cstv_Tstr_8

               u_interp = 0.25d0*(F_u(i,j,k)+F_u(i-1,j,k)+F_u(i,j+1,k)+F_u(i-1,j+1,k))

               F_orv(i,j,k) = Cstv_invT_m_8  * F_v(i,j,k) - Cstv_Beta_m_8 * ( &
                              t_interp * ( F_q(i,j+1,k) - F_q(i,j,k) ) * geomh_invDY_8             &
                            + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp )   &
                            + Cstv_Beta_m_8 * t_interp * mc_Jy(i,j,k) * ( &
                              Ver_wp_8%m(k)*half*( (F_q(i,j+1,k+1)-F_q(i,j+1,k ))*mc_iJz(i,j+1,k )   &
                                                 + (F_q(i,j  ,k+1)-F_q(i,j  ,k ))*mc_iJz(i,j  ,k ) ) &
                            + Ver_wm_8%m(k)*half*( (F_q(i,j+1,k  )-F_q(i,j+1,km))*mc_iJz(i,j+1,km)   &
                                                 + (F_q(i,j  ,k  )-F_q(i,j  ,km))*mc_iJz(i,j  ,km) ) )
               F_orv(i,j,k) = F_orv(i,j,k) + ((1d0-phy_bA_m_8)/Cstv_bA_m_8) * phy_vv_tend(i,j,k)
            end do
         end do

         !******************************************
         ! Compute Rw & Rt: RHS of w & T equations *
         !******************************************
         xtmp_8(:,:) = one
         do j = j0, jn
            do i = i0, in
               xtmp_8(i,j) = F_t(i,j,k) / real(Cstv_Tstr_8)
            end do
         end do
         call vlog( ytmp_8, xtmp_8, nij )
         w1=half/(cpd_8*Cstv_Tstr_8)
         do j= j0, jn
            do i= i0, in

               F_orw(i,j,k) = Cstv_invT_nh_8 * F_w(i,j,k) - Cstv_Beta_nh_8 *xtmp_8(i,j)* &
                              ( (F_q(i,j,k+1)-F_q(i,j,k))*mc_iJz(i,j,k)                 &
                              -grav_8*( one-one/xtmp_8(i,j) ) )

               F_ort(i,j,k) = Cstv_invT_8 * ( ytmp_8(i,j) - w1*(F_q(i,j,k+1)+F_q(i,j,k)) ) &
                             - Cstv_Beta_8 * mu_8 * F_w(i,j,k)
               F_ort(i,j,k) = F_ort(i,j,k) + ((1d0-phy_bA_t_8)/Cstv_bA_8) * 1./F_t(i,j,k) * phy_tv_tend(i,j,k)

               F_orf(i,j,k) = Cstv_invT_nh_8 * (ztht(i,j,k)-Ver_z_8%t(k)) &
                            - Cstv_Beta_nh_8 * ( Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km) - F_w(i,j,k) )
            end do
         end do

         !*****************************************
         ! Compute Rc: RHS of continuity equation *
         !*****************************************
         w1=epsi_8/grav_8
         do j= j0, jn
            do i= i0, in
               div = (F_u (i,j,k)-F_u (i-1,j,k))*geomh_invDXM_8(j)                                 &
                   + (F_v (i,j,k)*geomh_cyM_8(j)-F_v (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                   + (F_zd(i,j,k)-Ver_onezero(k)*F_zd(i,j,km))*Ver_idz_8%m(k)*Ver_wpstar_8(k)      &
                   + half * ( mc_Ix(i,j,k)*(F_u(i,j,k)+F_u(i-1,j,k))                               &
                   + mc_Iy(i,j,k)*(F_v(i,j,k)+F_v(i,j-1,k)) )                             &
                   + mc_Iz(i,j,k)*(Ver_wpA_8(k)*F_zd(i,j,k)+Ver_wmA_8(k)*Ver_onezero(k)*F_zd(i,j,km))
               F_orc (i,j,k) = Cstv_invT_8 *  w1 * F_q(i,j,k) +   Cstv_invT_8 * mc_logJz(i,j,k)  &
                             - Cstv_Beta_8 * ( div-epsi_8*(Ver_wp_8%m(k)*F_w(i,j,k)+Ver_onezero(k)*Ver_wm_8%m(k)*F_w(i,j,km)) )
               F_orc(i,j,k) = F_orc(i,j,k) + ((1d0-phy_bA_t_8)/Cstv_bA_8) * &
                             (              Ver_wp_8%m(k)*phy_tv_tend(i,j,k )/F_t(i,j,k ) + &
                             Ver_onezero(k)*Ver_wm_8%m(k)*phy_tv_tend(i,j,km)/F_t(i,j,km) )
            end do
         end do

      end do
!$omp end do
!$omp  end parallel

1000  format(3X,'COMPUTE THE RIGHT-HAND-SIDES: (S/R RHS_H)')
!
!     ---------------------------------------------------------------
!
      return
      end
