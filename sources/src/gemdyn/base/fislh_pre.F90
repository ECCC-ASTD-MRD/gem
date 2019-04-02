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
!**s/r pre_H -  Combine some rhs obtaining Rt' and Rc", the linear
!               contributions to the rhs of Helmholtz equation
!            -  Height-type vertical coordinate

!
      subroutine fislh_pre ( F_ru ,F_rv ,F_rt ,F_rw ,F_rc, F_rf, &
                            Minx, Maxx, Miny, Maxy, i0, j0, in, jn, Nk )
      use HORgrid_options
      use gem_options
      use geomh
      use tdpack
      use gmm_itf_mod
      use glb_ld
      use lun
      use cstv
      use ver
      use metric
      use fislh_sol

      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy, i0, j0, in, jn, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_ru,F_rv,F_rt
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_rw,F_rc,F_rf
!
!Author: Claude Girard, July 2017
!
      integer :: i, j, k, km
      real    :: w0
      real*8  :: div, c1, w1, w2
      real*8, parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      c1 = grav_8 * Cstv_tau_8

      call rpn_comm_xch_halo( F_ru, l_minx,l_maxx,l_miny,l_maxy, l_niu, l_nj,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( F_rv, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_njv,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

!*************************************
! Combination of governing equations *
!*************************************

!$omp parallel private(km,div,w1,w2)

!$omp do
      do k=1, l_nk

         do j= j0, jn
         do i= i0, in

            F_rc(i,j,k) = F_rc(i,j,k) - Cstv_invT_8 * mc_logJz(i,j,k)

!           Compute Rc' combining Ru, Rv and Rc
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            div  = (F_ru(i,j,k)-F_ru(i-1,j,k))*geomh_invDXM_8(j) &
                 + (F_rv(i,j,k)*geomh_cyM_8(j)-F_rv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                 + isol_i*half*(mc_Ix(i,j,k)*(F_ru(i,j,k)+F_ru(i-1,j,k)) &
                               +mc_Iy(i,j,k)*(F_rv(i,j,k)+F_rv(i,j-1,k)) )

            F_rc(i,j,k) = div - Cstv_invT_m_8 * F_rc(i,j,k)

!           Compute Rf'
!           ~~~~~~~~~~~
            w0 = Cstv_invT_nh_8*(ztht(i,j,k)-Ver_z_8%t(k))
            F_rf(i,j,k) = F_rf(i,j,k) - w0

!           Compute Rw'
!           ~~~~~~~~~~~
            F_rw(i,j,k) = F_rw(i,j,k) + Cstv_invT_nh_8*F_rf(i,j,k)

!           Compute Rt'
!           ~~~~~~~~~~~
            F_rt(i,j,k) = F_rt(i,j,k) + mu_8*F_rf(i,j,k)

!           Compute Rt" combining Rt'and Rw'
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            F_rt(i,j,k) = gama_8 * ( c1 * F_rt(i,j,k) + F_rw(i,j,k) )

         end do
         end do

      end do
!$omp end do

!$omp do
      do k=1,l_nk
         km=max(k-1,1)
         do j= j0, jn
         do i= i0, in
            w1= (Ver_idz_8%m(k) + (isol_i*mc_Iz(i,j,k) - epsi_8)*Ver_wpA_8(k))*Ver_zeronk(k)
            w2= (Ver_idz_8%m(k) - (isol_i*mc_Iz(i,j,k) - epsi_8)*Ver_wmA_8(k))*Ver_onezero(k)

!           Compute Rc"
!           ~~~~~~~~~~~
            F_rc(i,j,k) = F_rc(i,j,k) + Cstv_invT_m_8*epsi_8*&
                 (Ver_wp_8%m(k)*F_rf(i,j,k)+Ver_wm_8%m(k)*Ver_onezero(k)*F_rf(i,j,km))

!           Compute Rc"', combining Rc" and Rt"
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            F_rc(i,j,k) = F_rc(i,j,k) + w1 * F_rt(i,j,k) - w2 * F_rt(i,j,km)

         end do
         end do
      end do
!$omp end do

!     Apply lower boundary conditions
!
!$omp do
      do j= j0, jn
      do i= i0, in
         F_rt(i,j,l_nk) = F_rt(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rt(i,j,l_nk-1)
      end do
      end do
!$omp end do

!$omp end parallel

1000  format(3X,'UPDATE  THE RIGHT-HAND-SIDES: (S/R PRE_H)')
!
!     ---------------------------------------------------------------
!
      return
      end
