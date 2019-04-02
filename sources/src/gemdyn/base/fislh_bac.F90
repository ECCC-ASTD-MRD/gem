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
                             F_ru  , F_rv , F_rt , F_rw , F_rf ,       &
                             F_nu  , F_nv , F_nt , F_nw ,              &
                             Minx, Maxx, Miny, Maxy, ni, nj, Nk, i0, j0, in, jn )
      use fislh_sol
      use gem_options
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
      implicit none
!
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, ni,nj,Nk , i0, j0, in, jn

      real*8, dimension(ni,nj,Nk),                intent(in)  :: F_lhs_sol
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out) :: F_u,F_v,F_t,F_w,F_zd
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(out) :: F_q
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)  :: F_ru,F_rv,F_rt,F_rw,F_rf
      real,   dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)  :: F_nu,F_nv,F_nt,F_nw
!
!Author: Claude Girard, July 2017

#include <arch_specific.hf>

!
      integer :: i, j, k, km
      real*8  :: Buoy
      real*8, parameter :: zero=0.d0, one=1.d0, half=0.5d0
!     __________________________________________________________________
!
      if (Lun_debug_L) write(Lun_out,1000)

      do k=1,l_nk
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
                            + mc_css_H_8(i,j)   * (F_rt(i,j,l_nk)-F_nt(i,j,l_nk))
         end do
      end do
!$omp end do
!
      call rpn_comm_xch_halo(F_q,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0)


!$omp parallel private(km,Buoy)
!
!$omp do
      do k=1,l_nk
         km=max(k-1,1)

         do j= j0, jn
            do i= i0, l_niu-pil_e

   !           Compute U
   !           ~~~~~~~~~
               F_u(i,j,k) = Cstv_tau_m_8*(F_ru(i,j,k)-F_nu(i,j,k)  &
                          - (F_q(i+1,j,k)-F_q(i,j,k))*geomh_invDX_8(j) &
                          + isol_i*mc_Jx(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (F_q(i+1,j,k+1)-F_q(i+1,j,k ))*mc_iJz(i+1,j,k )   &
                                                 +(F_q(i  ,j,k+1)-F_q(i  ,j,k ))*mc_iJz(i  ,j,k ) ) &
                            +Ver_wm_8%m(k)*half*( (F_q(i+1,j,k  )-F_q(i+1,j,km))*mc_iJz(i+1,j,km)  &
                                                 +(F_q(i  ,j,k  )-F_q(i  ,j,km))*mc_iJz(i  ,j,km) ) ))
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
                                                 +(F_q(i,j  ,k+1)-F_q(i,j  ,k ))*mc_iJz(i,j  ,k ) ) &
                            +Ver_wm_8%m(k)*half*( (F_q(i,j+1,k  )-F_q(i,j+1,km))*mc_iJz(i,j+1,km)  &
                                                 +(F_q(i,j  ,k  )-F_q(i,j  ,km))*mc_iJz(i,j  ,km) ) ))
            end do
         end do
      end do
!$omp end do

!$omp do
      do k=1,l_nk-1
         do j= j0, jn
            do i= i0, in

   !           Compute zdot
   !           ~~~~~~~~~~~~
               F_zd(i,j,k) = Cstv_tau_m_8 * ( F_rt(i,j,k) - F_nt(i,j,k) &
                        - gama_8*((isol_i*mc_iJz(i,j,k) + isol_d*Ver_idz_8%t(k))*&
                                 (F_q(i,j,k+1)-F_q(i,j,k)) &
                            - half*mu_8*(F_q(i,j,k+1)+F_q(i,j,k))) )

   !           Compute w
   !           ~~~~~~~~~
               F_w(i,j,k) = - F_rf(i,j,k) + F_zd(i,j,k)

   !           Compute T
   !           ~~~~~~~~~
               Buoy = (F_q(i,j,k+1)-F_q(i,j,k))*      &
                      (isol_i*mc_iJz(i,j,k)+isol_d*Ver_idz_8%t(k))  &
                    + F_zd(i,j,k)*Cstv_invT_nh_8 - F_rw(i,j,k) + F_nw(i,j,k)
               F_t(i,j,k) = Cstv_Tstr_8 / (one - Buoy / grav_8 )
            end do
         end do
      end do
!$omp end do

!$omp do
      do j= j0, jn
         do i= i0, in

   !        Compute w(N) since zdot(N)=0
   !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            F_w(i,j,G_nk) = - F_rf(i,j,G_nk) + Ver_wmstar_8(G_nk)*F_zd(i,j,G_nk-1)

            F_zd(i,j,G_nk) = zero

   !        Compute T(N)
   !        ~~~~~~~~~~~~~
            Buoy = (F_q(i,j,G_nk+1)-F_q(i,j,G_nk))*     &
                   (isol_i*mc_iJz(i,j,G_nk) + isol_d*Ver_idz_8%t(G_nk))   &
                 + F_zd(i,j,G_nk-1)*Ver_wmstar_8(G_nk)*Cstv_invT_nh_8     &
                 - F_rw(i,j,G_nk) + F_nw(i,j,G_nk)
            F_t(i,j,G_nk) = Cstv_Tstr_8 / (one - Buoy / grav_8 )

         end do
      end do
!$omp end do

!$omp end parallel
!
1000  format (5X,'BACK SUBSTITUTION: (S/R BAC_H)')
!     __________________________________________________________________
!
      return
      end

