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
!---------------------------------- LICENCE END --------------------------------

!**s/r  bacp - backsubstitution: obtain new values of the variables:
!                                                   u,v,w,t,q,s,zd
!                from new P , the right-hand sides (Ru,Rv,Rt,Rw,Rf)
!                             and non-linear terms (Nu,Nv,Nt,Nw,Nf)

      subroutine bac ( F_lhs_sol       , F_sl ,F_fis      , &
                       F_u   , F_v     , F_w  , F_t       , &
                       F_s   , F_zd    , F_q  , F_nest_q  , &
                       F_ru  , F_rv    , F_rt , F_rw , F_rf , F_rb, &
                       F_nu  , F_nv    , F_nt , F_nw , F_nf , F_nb, &
                       Minx,Maxx,Miny,Maxy, ni,nj,Nk, i0, j0, k0, in, jn )

      use cstv
      use dynkernel_options
      use dyn_fisl_options
      use ctrl
      use gem_options
      use geomh
      use glb_ld
      use glb_pil
      use HORgrid_options
      use lun
      use tdpack
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, ni, nj, Nk, i0, j0, k0, in, jn
      real(kind=REAL64),  intent(in) :: F_lhs_sol (ni,nj,Nk)
      real, dimension(Minx:Maxx,Miny:Maxy),      intent(in) :: F_fis, F_sl, F_rb, F_nb
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(in) :: F_ru, F_rv, F_rt, F_rw,   &
                                                               F_rf, F_nu, F_nv, F_nt, F_nw, F_nf
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk+1), intent(in) :: F_nest_q
      real, dimension(Minx:Maxx,Miny:Maxy),        intent(out) :: F_s
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk),   intent(out) :: F_u, F_v, F_w, F_zd, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk+1), intent(out) :: F_q


      integer :: i, j, k, km, nij, k0t
      real(kind=REAL64)  :: w1, w2, w3, w4, Pbar, qbar
      real(kind=REAL64), dimension(i0:in,j0:jn) :: xtmp_8, ytmp_8
      real  , dimension(Minx:Maxx,Miny:Maxy,Nk+1) :: GP
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=.5d0
!     __________________________________________________________________
!
      if (Ctrl_testcases_adv_L) then
         call canonical_cases ("BAC")
         return
      end if

      if (Lun_debug_L) write(Lun_out,1000)

      nij = (in - i0 + 1)*(jn - j0 + 1)

      k0t=k0
      if(Schm_opentop_L) k0t=k0-1

      do k=k0,l_nk
         do j= j0, jn
            do i= i0, in
               GP(i,j,k) = sngl(F_lhs_sol(i,j,k))
            end do
         end do
      end do

!
!     Compute P at top and bottom
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      if ( Schm_opentop_L ) then
         GP(:,:,1:k0-2) = 0.0

         do j= j0, jn
            do i= i0, in
               GP(i,j,k0-1) = Ver_alfat_8 * F_lhs_sol(i,j,k0) &
                            + Ver_cst_8*(F_rb(i,j)-F_nb(i,j))
            end do
         end do

      end if

      do j= j0, jn
         do i= i0, in
            GP(i,j,l_nk+1)  = Ver_alfas_8 * F_lhs_sol(i,j,l_nk)   &
                            + Ver_betas_8 * F_lhs_sol(i,j,l_nk-1) &
                            - Ver_css_8*(F_rt(i,j,l_nk)-F_nt(i,j,l_nk))
         end do
      end do

      call rpn_comm_xch_halo(GP,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1, &
                  G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

!
!     Compute w
!     ~~~~~~~~~

      do k=k0t,l_nk
         km = max(k-1,k0t)
         w1 = Cstv_tau_m_8*Rgasd_8*Cstv_Tstr_8/grav_8
         do j= j0, jn
            do i= i0, in
               Pbar = Ver_wpstar_8(k)*GP(i,j,k+1) &
                    + Ver_wmstar_8(k)*half*(GP(i,j,k)+GP(i,j,km))
               Pbar = Ver_wp_8%t(k)*Pbar+Ver_wm_8%t(k)*GP(i,j,k)
               F_w(i,j,k) = w1 * ( F_rf(i,j,k) - F_nf(i,j,k) &
               + gama_8 * ( (GP(i,j,k+1)-GP(i,j,k))*Ver_idz_8%t(k) &
                                        + cappa_8 * Pbar ) )
            end do
         end do
      end do

      if(.not.Dynamics_hydro_L) then

!        Compute q
!        ~~~~~~~~~
!
!        N.B.  Top Boundary condition:
!                 Closed Top(k0 == 1):  F_q(i,j,k0) = 0
!                   Open Top(k0 /= 1):  F_q(i,j,k0t) is externally specified

         if (Schm_opentop_L) then

            do j= j0, jn
               do i= i0, in
                  F_q(i,j,k0t)=F_nest_q(i,j,k0t)
               end do
            end do

         end if
!
!        Note : we cannot use omp on loop k
!               due to vertical dependency F_q(i,j,k)
         do k=k0t,l_nk
            km=max(k-1,1)
            w4 = one/(one+Ver_wp_8%t(k)*Ver_wpstar_8(k)*Ver_dz_8%t(k))
            w3 = half*Ver_wp_8%t(k)*Ver_wmstar_8(k)*Ver_dz_8%t(k)*w4
            w2 = (one-(Ver_wm_8%t(k)+half*Ver_wp_8%t(k)*Ver_wmstar_8(k))*Ver_dz_8%t(k))*w4
            w1 = Ver_dz_8%t(k)/grav_8*w4

            do j= j0, jn
               do i= i0, in
                  F_q(i,j,k+1) = - w3 * F_q(i,j,km)   &
                                 + w2 * F_q(i,j,k)    &
                                 - w1 * ( F_rw(i,j,k) - F_nw(i,j,k)  &
                                    - Cstv_invT_nh_8 * F_w(i,j,k)  )
               end do
            end do

         end do

      end if

!     Compute U & V
!     ~~~~~~~~~~~~~

      do k=k0,l_nk
         do j= j0, jn
            do i= i0, l_niu-pil_e
               F_u(i,j,k) = Cstv_tau_m_8*(F_ru(i,j,k)-F_nu(i,j,k) - (GP(i+1,j,k)-GP(i,j,k))*geomh_invDXMu_8(j))
            end do
         end do

         do j= j0, l_njv-pil_n
            do i= i0, in
               F_v(i,j,k) = Cstv_tau_m_8*(F_rv(i,j,k)-F_nv(i,j,k) - (GP(i,j+1,k)-GP(i,j,k))*geomh_invDYMv_8(j))
            end do
         end do
      end do

!     Compute s
!     ~~~~~~~~~
      w1 = one/(Rgasd_8*Cstv_Tstr_8)

      do j= j0, jn
         do i= i0, in
            F_s(i,j) = w1*(GP(i,j,l_nk+1)-F_fis(i,j))-F_q(i,j,l_nk+1)
         end do
      end do

!     Compute zdot
!     ~~~~~~~~~~~~

!     N.B.  Top Boundary condition:
!              Closed Top(k0t == 1):  F_zd(i,j,k0t) = 0
!                Open Top(k0t /= 1):  F_zd(i,j,k0t) is computed

      do k=k0t,l_nk-1
         w1=gama_8*Ver_idz_8%t(k)
         w2=gama_8*epsi_8
         w3=Cstv_invT_8*Cstv_bar1_8
         do j= j0, jn
            do i= i0, in
               Pbar= Ver_wp_8%t(k)*GP(i,j,k+1)+Ver_wm_8%t(k)*GP(i,j,k)
               qbar=(Ver_wp_8%t(k)*F_q(i,j,k+1)+Ver_wm_8%t(k)*F_q(i,j,k))
               F_zd(i,j,k)=-Cstv_tau_m_8*( F_rt(i,j,k)- F_nt(i,j,k) &
                          + w1 * ( GP(i,j,k+1)-GP(i,j,k) ) - w2 * Pbar ) &
                          - w3 * ( Ver_b_8%t(k)*F_s(i,j) &
                                  +Ver_c_8%t(k)*F_sl(i,j)+qbar)
            end do
         end do
      end do

!     Compute FI' (into GP)
!     ~~~~~~~~~~~

      do k=k0t,l_nk
         km=max(k-1,1)
         w1=Rgasd_8*Cstv_Tstr_8
         do j= j0, jn
            do i= i0, in
               GP(i,j,k)=GP(i,j,k)-w1*(Ver_b_8%m(k)*F_s(i,j)  &
                                      +Ver_c_8%m(k)*F_sl(i,j) &
                                      +F_q(i,j,k) )
            end do
         end do
      end do


      do j= j0, jn
         do i= i0, in
            GP(i,j,l_nk+1)=F_fis(i,j)
         end do
      end do

!     Compute T
!     ~~~~~~~~~

      do k=k0t,l_nk
         km=max(k-1,1)
         do j= j0, jn
            do i= i0, in
               qbar=Ver_wpstar_8(k)*F_q(i,j,k+1)+Ver_wmstar_8(k)*half*(F_q(i,j,k)+F_q(i,j,km))
               qbar=Ver_wp_8%t(k)*qbar+Ver_wm_8%t(k)*F_q(i,j,k)
               ytmp_8(i,j)=-qbar
            end do
         end do
         call vexp( xtmp_8, ytmp_8, nij )
         do j= j0, jn
            do i= i0, in
               xtmp_8(i,j)=xtmp_8(i,j)*(one+Ver_dbdz_8%t(k)*F_s(i,j)  &
                                           +Ver_dcdz_8%t(k)*F_sl(i,j))
            end do
         end do
         call vrec ( ytmp_8, xtmp_8, nij )
         w1=Ver_idz_8%t(k)/rgasd_8
         do j= j0, jn
            do i= i0, in
               F_t(i,j,k)=ytmp_8(i,j)*(Cstv_Tstr_8-w1*(GP(i,j,k+1)-GP(i,j,k)))
            end do
         end do
      end do

      if(Schm_autobar_L) then
         F_t=Cstv_Tstr_8 ; F_zd=0. ! not necessary but safer
      end if
!
1000  format (5X,'BACK SUBSTITUTION: (S/R BAC)')
!     __________________________________________________________________
!
      return
      end

