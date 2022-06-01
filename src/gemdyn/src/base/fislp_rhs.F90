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
!*s/r fislp_rhs - compute the right-hand sides: Ru, Rv, Rc, Rt, Rw, Rf,
!                 save the results for next iteration in the o's
!
      subroutine fislp_rhs ( F_oru, F_orv, F_orc, F_ort, F_orw, F_orf, &
                             F_u, F_v, F_w, F_t, F_s, F_zd, F_q, &
                             F_sl, F_fis, Minx, Maxx, Miny, Maxy, Nk )
      use coriolis
      use cstv
      use dynkernel_options
      use dyn_fisl_options
      use gem_options
      use hvdif_options
      use geomh
      use glb_ld
      use gmm_phy, only: phy_uu_tend, phy_vv_tend, phy_tv_tend, phy_cplm, phy_cplt
      use HORgrid_options
      use lun
      use tdpack
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy), intent(inout) :: F_s
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_fis, F_sl
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(in)    :: F_w, F_zd
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(inout) :: F_u, F_v, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(out)   :: F_oru, F_orv, F_orc, F_ort, F_orw, F_orf
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk+1), intent(inout) :: F_q

      integer :: i0,  in,  j0,  jn
      integer :: i0u, inu, i0v, inv
      integer :: j0u, jnu, j0v, jnv

      integer :: i, j, k, km, kq, nij, jext
      real(kind=REAL64) :: tdiv, BsPqbarz, barz, barzp,  &
                u_interp, v_interp, t_interp, mu_interp, &
                w1, w2
      real(kind=REAL64)  xtmp_8(l_ni,l_nj), ytmp_8(l_ni,l_nj)
      real, dimension(Minx:Maxx,Miny:Maxy,l_nk+1) :: BsPq, FI
      real, dimension(Minx:Maxx,Miny:Maxy,l_nk) :: MU
      real(kind=REAL64), parameter :: one=1.d0, zero=0.d0, half=0.5d0
      logical :: smago_in_rhs_L
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      smago_in_rhs_L=.false.

!     Common coefficients

      jext = Grd_bsc_ext1 + 1

      nij = l_ni*l_nj

!     Indices to compute Rc, Rt, Rf, Rw
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

      call diag_fip(FI, F_s, F_sl, F_t, F_q, F_fis, l_minx,l_maxx,l_miny,l_maxy, l_nk, &
                    i0u,inu+1,j0v,jnv+1)

      if (Dynamics_hydro_L) then
         MU = 0.
      else
         call diag_mu(MU, F_q, F_s, F_sl, l_minx,l_maxx,l_miny,l_maxy, l_nk, &
                      i0u,inu+1,j0v,jnv+1)
      end if

      do k=1,l_nk+1
         kq=max(2,k)
         do j=j0v,jnv+1
            do i=i0u,inu+1
               BsPq(i,j,k) = Ver_b_8%m(k)*F_s(i,j) &
                            +Ver_c_8%m(k)*F_sl(i,j) + F_q(i,j,k)
            end do
         end do
      end do

      do k = 1,l_nk
         km=max(k-1,1)

!********************************
! Compute Ru: RHS of U equation *
!********************************


!     Compute V barX in wk2
!     ~~~~~~~~~~~~~~~~~~~~~
         do j= j0u, jnu
            do i= i0u, inu

               barz  = Ver_wpM_8(k)*MU(i  ,j,k)+Ver_wmM_8(k)*MU(i  ,j,km)
               barzp = Ver_wpM_8(k)*MU(i+1,j,k)+Ver_wmM_8(k)*MU(i+1,j,km)
               mu_interp = ( barz + barzp)*half

               barz  = Ver_wpM_8(k)*F_t(i  ,j,k)+Ver_wmM_8(k)*F_t(i  ,j,km)
               barzp = Ver_wpM_8(k)*F_t(i+1,j,k)+Ver_wmM_8(k)*F_t(i+1,j,km)
               t_interp = (barz + barzp)*half

               v_interp = 0.25d0*(F_v(i,j,k)+F_v(i,j-1,k)+F_v(i+1,j,k)+F_v(i+1,j-1,k))

               F_oru(i,j,k) = Cstv_invT_m_8  * F_u(i,j,k) - Cstv_Beta_m_8 * ( &
                              rgasd_8 * t_interp * ( BsPq(i+1,j,k) - BsPq(i,j,k) ) * geomh_invDXMu_8(j)  &
                            + ( one + mu_interp)* (   FI(i+1,j,k) -   FI(i,j,k) ) * geomh_invDXMu_8(j)  &
                            - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * F_u(i,j,k) ) * v_interp ) + &
                              ((1d0-phy_cplm(i,j))/Cstv_bA_m_8) * phy_uu_tend(i,j,k)
            end do
         end do

!********************************
! Compute Rv: RHS of V equation *
!********************************

!     Compute U barY in wk2
!     ~~~~~~~~~~~~~~~~~~~~~
         do j = j0v, jnv
            do i = i0v, inv

               barz  = Ver_wpM_8(k)*MU(i,j  ,k)+Ver_wmM_8(k)*MU(i,j  ,km)
               barzp = Ver_wpM_8(k)*MU(i,j+1,k)+Ver_wmM_8(k)*MU(i,j+1,km)
               mu_interp = ( barz + barzp )*half

               barz  = Ver_wpM_8(k)*F_t(i,j  ,k)+Ver_wmM_8(k)*F_t(i,j  ,km)
               barzp = Ver_wpM_8(k)*F_t(i,j+1,k)+Ver_wmM_8(k)*F_t(i,j+1,km)
               t_interp = ( barz + barzp)*half

               u_interp = 0.25d0*(F_u(i,j,k)+F_u(i-1,j,k)+F_u(i,j+1,k)+F_u(i-1,j+1,k))

               F_orv(i,j,k) = Cstv_invT_m_8  * F_v(i,j,k) - Cstv_Beta_m_8 * ( &
                              rgasd_8 * t_interp * ( BsPq(i,j+1,k) - BsPq(i,j,k) ) * geomh_invDYMv_8(j) &
                                  + ( one + mu_interp)* (   FI(i,j+1,k) -   FI(i,j,k) ) * geomh_invDYMv_8(j) &
                                          + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp ) + &
                                          ((1d0-phy_cplm(i,j))/Cstv_bA_m_8) * phy_vv_tend(i,j,k)
            end do
         end do

!********************************************
! Compute Rt: RHS of thermodynamic equation *
! Compute Rx: RHS of Ksidot equation        *
! Compute Rc: RHS of Continuity equation    *
!********************************************

         xtmp_8(:,:) = one
         do j = j0, jn
            do i = i0, in
               xtmp_8(i,j) = F_t(i,j,k) / Cstv_Tstr_8
            end do
         end do
         call gem_vlog( ytmp_8, xtmp_8, nij )
         do j= j0, jn
            do i= i0, in
               w1=Ver_wpstar_8(k)*BsPq(i,j,k+1)+Ver_wmstar_8(k)*half*(BsPq(i,j,k)+BsPq(i,j,km))
               BsPqbarz = Ver_wp_8%t(k)*w1+Ver_wm_8%t(k)*BsPq(i,j,k)
               F_ort(i,j,k) = Cstv_invT_8 * ( ytmp_8(i,j) - Cstv_bar1_8 * cappa_8 * BsPqbarz ) &
                            + Cstv_Beta_8 * cappa_8 *(Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km)) + &
                            ((1d0-phy_cplt(i,j))/Cstv_bA_8) * 1./F_t(i,j,k) * phy_tv_tend(i,j,k)
            end do
         end do

         xtmp_8(:,:) = one
         do j = j0, jn
            do i = i0, in
               xtmp_8(i,j) = one + Ver_dbdz_8%m(k) * F_s(i,j)  &
                                 + Ver_dcdz_8%m(k) * F_sl(i,j)
            end do
         end do
         call gem_vlog( ytmp_8, xtmp_8, nij)
         do j = j0, jn
            do i = i0, in
               tdiv = (F_u (i,j,k)-F_u (i-1,j,k))*geomh_invDXM_8(j) &
                    + (F_v (i,j,k)*geomh_cyM_8(j)-F_v (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + (F_zd(i,j,k)-Ver_onezero(k)*F_zd(i,j,km))*Ver_idz_8%m(k) &
                    + Ver_wpC_8(k) * F_zd(i,j,k) + Ver_wmC_8(k) * Ver_onezero(k) * F_zd(i,j,km)
               F_orc (i,j,k) = Cstv_invT_8 * ( Cstv_bar1_8*(Ver_b_8%m(k)*F_s(i,j)  &
                                                           +Ver_c_8%m(k)*F_sl(i,j)) + ytmp_8(i,j) ) &
                             - Cstv_Beta_8 * tdiv
            end do
         end do

   !********************************************
   ! Compute Rw: RHS of  w equation            *
   ! Compute Rf: RHS of  FI equation           *
   ! Compute Rq: RHS of  qdot equation         *
   !********************************************

         if(.not.Dynamics_hydro_L) then
            do j= j0, jn
               do i= i0, in
                  F_orw(i,j,k) = Cstv_invT_nh_8 * F_w(i,j,k) &
                               + Cstv_Beta_nh_8 * grav_8 * MU(i,j,k)
               end do
            end do
         end if
         do j= j0, jn
            do i= i0, in
               w1=Ver_wpstar_8(k)*FI(i,j,k+1)+Ver_wmstar_8(k)*half*(FI(i,j,k)+FI(i,j,km))
               w2=Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km)
               F_orf(i,j,k) = Cstv_invT_8 * ( Ver_wp_8%t(k)*w1+Ver_wm_8%t(k)*FI(i,j,k) )  &
                            + Cstv_Beta_8 * Rgasd_8 * Cstv_Tstr_8 * w2 &
                            + Cstv_Beta_8 * grav_8 * F_w(i,j,k)
            end do
         end do

      end do

      if (smago_in_rhs_L) then
         call smago_in_rhs ( F_oru, F_orv, F_orw, F_ort, F_u, F_v, F_w, F_t, F_s, &
                          F_sl,l_minx,l_maxx,l_miny,l_maxy,G_nk )
      end if

      !Reinforce F_oru/F_orv constant in vertical when AUTOBAR
      !-------------------------------------------------------
      if (Schm_autobar_L) then

         do k=2,l_nk
            do j=j0u,jnu
               do i=i0u,inu
                  F_oru(i,j,k) = F_oru(i,j,1)
               end do
            end do
         end do

         do k=2,l_nk
            do j=j0v,jnv
               do i=i0v,inv
                  F_orv(i,j,k) = F_orv(i,j,1)
               end do
            end do
         end do

      end if

1000  format(3X,'COMPUTE THE RIGHT-HAND-SIDES: (S/R RHS)')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine fislp_rhs
