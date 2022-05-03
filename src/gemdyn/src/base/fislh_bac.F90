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

      subroutine fislh_bac ( F_dt_8, i0, j0, k0, in, jn ,k0t )
      use fislh_sol
      use gem_options
      use dynkernel_options
      use dyn_fisl_options
      use geomh
      use HORgrid_options
      use tdpack
      use gmm_vt0
      use mem_tstp
      use glb_pil
      use glb_ld
      use cstv
      use ver
      use metric
      use ctrl
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

!Author: Claude Girard, July 2017
!        Syed Husain, June 2019 (major revision)
!        Abdessamad Qaddouri, July 2019 (opentop)

      integer :: i, j, k, km
      real :: w5,w6
      real(kind=REAL64) :: tau_m_8, tau_nh_8, invT_m_8, invT_nh_8, Buoy
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      if (Ctrl_testcases_adv_L) then
!$omp single
         call canonical_cases ("BAC")
!$omp end single
         return
      end if

      tau_m_8  = F_dt_8 * Cstv_bA_m_8
      tau_nh_8 = F_dt_8 * Cstv_bA_nh_8
      invT_m_8 = one/tau_m_8
      invT_nh_8= one/tau_nh_8

!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
            do i= i0, in
               qt0(i,j,k) = sngl(lhs_sol(i,j,k))
            end do
         end do
      end do
!$omp enddo

!$omp do
      do j= j0, jn
         do i= i0, in
            qt0(i,j,l_nk+1) = GVM%mc_alfas_H_8(i,j) * qt0(i,j,l_nk)   &
                         - GVM%mc_betas_H_8(i,j) * qt0(i,j,l_nk-1) &
                         + GVM%mc_css_H_8(i,j)   * (rhst(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhst(i,j,l_nk-1) &
                            + invT_m_8*(rhsf(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhsf(i,j,l_nk-1))) &
                         - GVM%mc_css_H_8(i,j)   * (nl_t(i,j,l_nk ) - Ver_wmstar_8(G_nk)*nl_t(i,j,l_nk-1))
         end do
      end do
!$omp enddo

      if ( Schm_opentop_L ) then
!$omp do
         do k=1, k0-2
            qt0(:,:,k) = 0.0
         end do
!$omp enddo
!$omp do
         do j= j0, jn
            do i= i0, in
               w5= GVM%mc_alfat_8(i,j) * lhs_sol(i,j,k0)
               w6= GVM%mc_cst_8  (i,j) * (rhsb(i,j)-nl_b(i,j))
               qt0(i,j,k0-1) = w5 - w6
            end do
         end do
!$omp enddo
      end if

!$omp single
      call rpn_comm_xch_halo(qt0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!$omp end single

!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
            km=max(k-1,1)
!DIR$ SIMD
            do i= i0, l_niu-pil_e

   !           Compute U
   !           ~~~~~~~~~
               ut0(i,j,k) = tau_m_8*(rhsu(i,j,k)-nl_u(i,j,k)  &
                          - (qt0(i+1,j,k)-qt0(i,j,k))*geomh_invDX_8(j) &
                          + isol_i*GVM%mc_Jx_8(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (qt0(i+1,j,k+1)-qt0(i+1,j,k ))*GVM%mc_iJz_8(i+1,j,k )   &
                                                + (qt0(i  ,j,k+1)-qt0(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k ) ) &
                          + Ver_wm_8%m(k)*half*( (qt0(i+1,j,k  )-qt0(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)  &
                                               + (qt0(i  ,j,k  )-qt0(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km) ) ))
            end do
         end do
      end do
!$omp enddo nowait
      
!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, l_njv-pil_n
            km=max(k-1,1)
!DIR$ SIMD
            do i= i0, in

   !           Compute V
   !           ~~~~~~~~~
               vt0(i,j,k) = tau_m_8*(rhsv(i,j,k)-nl_v(i,j,k) &
                          - (qt0(i,j+1,k)-qt0(i,j,k))*geomh_invDYMv_8(j) &
                          + isol_i*GVM%mc_Jy_8(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (qt0(i,j+1,k+1)-qt0(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )   &
                                                + (qt0(i,j  ,k+1)-qt0(i,j  ,k ))*GVM%mc_iJz_8(i,j  ,k ) ) &
                            +Ver_wm_8%m(k)*half*( (qt0(i,j+1,k  )-qt0(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)  &
                                                + (qt0(i,j  ,k  )-qt0(i,j  ,km))*GVM%mc_iJz_8(i,j  ,km) ) ))
            end do
         end do
      end do
!$omp enddo nowait

!$omp do collapse(2)
      do k=k0t, l_nk
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in

   !           Compute w
   !           ~~~~~~~~~
               wt0(i,j,k) = tau_m_8 * ( (rhst(i,j,k) - nl_t(i,j,k)) * Cstv_bar1_8 &
                          - gama_8*((isol_i*GVM%mc_iJz_8(i,j,k) + isol_d*Ver_idz_8%t(k))&
                            *(qt0(i,j,k+1)-qt0(i,j,k)) &
                          - Cstv_bar1_8*half*mu_8*(qt0(i,j,k+1)+qt0(i,j,k))))

   !           Compute zdot
   !           ~~~~~~~~~~~~
               zdt0(i,j,k) = (rhsf(i,j,k) + wt0(i,j,k))*Ver_zeronk(k)


   !           Compute T
   !           ~~~~~~~~~
               Buoy = (qt0(i,j,k+1)-qt0(i,j,k))*      &
                      (isol_i*GVM%mc_iJz_8(i,j,k)+isol_d*Ver_idz_8%t(k))  &
                    + wt0(i,j,k)*invT_nh_8 - (rhsw(i,j,k) - nl_w(i,j,k)) * Cstv_bar1_8
               tt0(i,j,k) = Cstv_Tstr_8 / (one - Buoy / grav_8 )

            end do
         end do
      end do
!$omp enddo

      if (Schm_autobar_L) then
!$omp single
         tt0=Cstv_Tstr_8 ; zdt0=0. ! not necessary but safer
!$omp end single
      end if
!
!     ---------------------------------------------------------------
!
      return
      end subroutine fislh_bac
