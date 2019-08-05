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

!**s/r  bac_H- backsubstitution: obtain new values of the variables:u,v,w,t,q,zd
!                from new q , the right-hand sides (Ru,Rv,Rw,Rt,Rf)
!                             and non-linear terms (Nu,Nv,Nw,Nt)
!             - Height-type vertical coordinate
!

      subroutine fislh_bac ( F_lhs_sol, &
                             F_u   , F_v  , F_t  , F_w , F_zd  , F_q , &
                             F_ru  , F_rv , F_rt , F_rw , F_rf , F_rb, &
                             F_nu  , F_nv , F_nt , F_nw , F_nb,          &
                             Minx, Maxx, Miny, Maxy, ni, nj, Nk, i0, j0,k0, in, jn )
      use fislh_sol
      use gem_options
      use dynkernel_options
      use dyn_fisl_options
      use geomh
      use HORgrid_options
      use tdpack
      use glb_pil
      use glb_ld
      use lun
      use cstv
      use ver
      use metric
      use ctrl
      use, intrinsic :: iso_fortran_env
      implicit none
!
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, ni,nj,Nk , i0, j0, k0,in, jn

      real(kind=REAL64), dimension(ni,nj,Nk),                intent(in)  :: F_lhs_sol
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out) :: F_u,F_v,F_t,F_w,F_zd
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(out) :: F_q
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)  :: F_ru,F_rv,F_rt,F_rw,F_rf
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)  :: F_nu,F_nv,F_nt,F_nw
      real, dimension(Minx:Maxx,Miny:Maxy),       intent(in) ::  F_rb ,F_nb
!
!Author: Claude Girard, July 2017
!        Syed Husain, June 2019 (major revision)
!        Abdessamad Qaddouri, July 2019 (opentop)


#include <arch_specific.hf>

!
      integer :: i, j, k, km, k0t
      real(kind=REAL64)  :: Buoy
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!     __________________________________________________________________
!
      if (Ctrl_testcases_adv_L) then
         call canonical_cases ("BAC")
         return
      end if

      if (Lun_debug_L) write(Lun_out,1000)
      k0t=k0
      if ( Schm_opentop_L ) k0t=k0-1


      do k=k0, l_nk
         do j= j0, jn
            do i= i0, in
               F_q(i,j,k) = sngl(F_lhs_sol(i,j,k))
            end do
         end do
      end do
!
!$omp do
      do j= j0, jn
         do i= i0, in
            F_q(i,j,l_nk+1) = mc_alfas_H_8(i,j) * F_q(i,j,l_nk)   &
                         - mc_betas_H_8(i,j) * F_q(i,j,l_nk-1) &
                         + mc_css_H_8(i,j)   * (F_rt(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rt(i,j,l_nk-1) &
                            + Cstv_invT_m_8*(F_rf(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rf(i,j,l_nk-1))) &
                         - mc_css_H_8(i,j)   * (F_nt(i,j,l_nk ) - Ver_wmstar_8(G_nk)*F_nt(i,j,l_nk-1))
         end do
      end do
!$omp end do
!
!
      if ( Schm_opentop_L ) then
         F_q(:,:,1:k0-2) = 0.0
!$omp do
         do j= j0, jn
            do i= i0, in
               F_q(i,j,k0-1) = Ver_alfat_8 * F_lhs_sol(i,j,k0) &
                         - Ver_cst_8*(F_rb(i,j)-F_nb(i,j))
            end do
         end do
!$omp enddo
      end if
!

!$omp single
      call rpn_comm_xch_halo(F_q,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!$omp end single

!$omp parallel private(km,Buoy)
!
!$omp do
      do k=k0, l_nk
         km=max(k-1,1)

         do j= j0, jn
            do i= i0, l_niu-pil_e

   !           Compute U
   !           ~~~~~~~~~
               F_u(i,j,k) = Cstv_tau_m_8*(F_ru(i,j,k)-F_nu(i,j,k)  &
                          - (F_q(i+1,j,k)-F_q(i,j,k))*geomh_invDX_8(j) &
                          + isol_i*mc_Jx(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (F_q(i+1,j,k+1)-F_q(i+1,j,k ))*mc_iJz(i+1,j,k )   &
                                                + (F_q(i  ,j,k+1)-F_q(i  ,j,k ))*mc_iJz(i  ,j,k ) ) &
                          + Ver_wm_8%m(k)*half*( (F_q(i+1,j,k  )-F_q(i+1,j,km))*mc_iJz(i+1,j,km)  &
                                               + (F_q(i  ,j,k  )-F_q(i  ,j,km))*mc_iJz(i  ,j,km) ) ))
            end do
         end do

         do j= j0, l_njv-pil_n
            do i= i0, in

   !           Compute V
   !           ~~~~~~~~~
               F_v(i,j,k) = Cstv_tau_m_8*(F_rv(i,j,k)-F_nv(i,j,k) &
                          - (F_q(i,j+1,k)-F_q(i,j,k))*geomh_invDYMv_8(j) &
                          + isol_i*mc_Jy(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (F_q(i,j+1,k+1)-F_q(i,j+1,k ))*mc_iJz(i,j+1,k )   &
                                                + (F_q(i,j  ,k+1)-F_q(i,j  ,k ))*mc_iJz(i,j  ,k ) ) &
                            +Ver_wm_8%m(k)*half*( (F_q(i,j+1,k  )-F_q(i,j+1,km))*mc_iJz(i,j+1,km)  &
                                                + (F_q(i,j  ,k  )-F_q(i,j  ,km))*mc_iJz(i,j  ,km) ) ))
            end do
         end do
      end do
!$omp end do

!$omp do
      do k=k0t, l_nk
         do j= j0, jn
            do i= i0, in

   !           Compute w
   !           ~~~~~~~~~
               F_w(i,j,k) = Cstv_tau_m_8 * ( F_rt(i,j,k) - F_nt(i,j,k) &
                          - gama_8*((isol_i*mc_iJz(i,j,k) + isol_d*Ver_idz_8%t(k))&
                            *(F_q(i,j,k+1)-F_q(i,j,k)) &
                          - Cstv_bar1_8*half*mu_8*(F_q(i,j,k+1)+F_q(i,j,k))))

   !           Compute zdot
   !           ~~~~~~~~~~~~
               F_zd(i,j,k) = (F_rf(i,j,k) + F_w(i,j,k))*Ver_zeronk(k)


   !           Compute T
   !           ~~~~~~~~~
               Buoy = (F_q(i,j,k+1)-F_q(i,j,k))*      &
                      (isol_i*mc_iJz(i,j,k)+isol_d*Ver_idz_8%t(k))  &
                    + F_w(i,j,k)*Cstv_invT_nh_8 - F_rw(i,j,k) + F_nw(i,j,k)
               F_t(i,j,k) = Cstv_Tstr_8 * (Buoy / grav_8+one )
            end do
         end do
      end do
!$omp enddo

      if (Schm_autobar_L) then
         F_t=Cstv_Tstr_8 ; F_zd=0. ! not necessary but safer
      end if

!$omp end parallel
!
1000  format (5X,'BACK SUBSTITUTION: (S/R BAC_H)')
!     __________________________________________________________________
!
      return
      end
