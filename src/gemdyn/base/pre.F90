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
!**s/r prep - Add metric corrections to r.h.s. of momentum equations.
!               Compute advective contributions on geopotential grid.
!               Interpolate advection contribution from geopotential
!               grid to wind grids. Update r.h.s with advective
!               contributions.
!               Combine some rhs obtaining Rt", Rf" and Rc", the linear
!               contributions to the rhs of Helmholtz equation
!
      subroutine pre ( F_ru, F_rv, F_fis, F_rc, F_rt, &
                       F_rw, F_rf, F_rb, F_nest_t, Minx, Maxx, Miny, Maxy, &
                       i0, j0, in, jn, k0, nk )
      use cstv
      use gem_options
      use dyn_fisl_options
      use geomh
      use glb_ld
      use gmm_itf_mod
      use HORgrid_options
      use lun
      use tdpack
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, i0, j0, in, jn, k0, nk
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_fis
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_rb
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_ru, F_rv, F_rt, F_rf
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: F_rc
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) ::  F_rw, F_nest_t


      integer :: i, j, k, km, k0t
      real(kind=REAL64)  :: rdiv, w1, w2, w3, w4
      real    :: w_rt
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0 , &
                           alpha1=-1.d0/16.d0 , alpha2=9.d0/16.d0
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      k0t = k0
      if(Schm_opentop_L) k0t= k0-1

      call rpn_comm_xch_halo( F_ru, l_minx,l_maxx,l_miny,l_maxy, l_niu, l_nj,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( F_rv, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_njv,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

!*************************************
! Combination of governing equations *
!*************************************

      do k=k0, l_nk
         km=max(k-1,1)
         w1= Ver_igt_8*Ver_wpA_8(k)
         w2= Ver_igt_8*Ver_wmA_8(k)*Ver_onezero(k)
         do j= j0, jn
            do i= i0, in

   !           Combine continuity & w equations
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               F_rc(i,j,k) = F_rc(i,j,k) - w1*F_rw(i,j,k ) - w2*F_rw(i,j,km)

   !           Compute the divergence of the RHS of momentum equations
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               rdiv  = (F_ru(i,j,k)-F_ru(i-1,j,k))*geomh_invDXM_8(j) &
                     + (F_rv(i,j,k)*geomh_cyM_8(j)-F_rv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j)

   !           Combine divergence & continuity equations : Rc"
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               F_rc(i,j,k) = rdiv - F_rc(i,j,k) / Cstv_tau_m_8

            end do
         end do
      end do

      if (Schm_opentop_L) then

         do j= j0, jn
            do i= i0, in
               F_rb(i,j) = F_rt(i,j,k0t)
            end do
         end do

      end if


      do k=k0t,l_nk

!        Compute Rt" & Rf"
!        ~~~~~~~~~~~~~~~~~

         w1 = cappa_8 /( Rgasd_8 * Cstv_Tstr_8 )
         w2 = Cstv_invT_m_8 / ( cappa_8 + epsi_8 )
         do j= j0, jn
            do i= i0, in
   !           Combine Rt and Rw
   !           ~~~~~~~~~~~~~~~~~
               w_rt = F_rt(i,j,k) + Ver_igt_8 * F_rw(i,j,k)

   !           Compute Rt"
   !           ~~~~~~~~~~~
               F_rt(i,j,k) = w2 * ( w_rt + Ver_igt2_8 * F_rf(i,j,k) )

   !           Compute Rf"
   !           ~~~~~~~~~~~
               F_rf(i,j,k) = w2 * ( w_rt - w1 * F_rf(i,j,k) )
            end do
         end do

      end do

      do j= j0, jn
         do i= i0, in
   !        Adjust Rt" at last level : Rt"'
   !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            F_rt(i,j,l_nk) = (F_rt(i,j,l_nk)-Ver_wmstar_8(l_nk)*F_rt(i,j,l_nk-1)) &
                             /Ver_wpstar_8(l_nk)
         end do
      end do

!************************************************************
! The linear contributions to the RHS of Helmholtz equation *
!************************************************************

!     Finish computations of RP(in Rc), combining Rc", Rt", Rf"
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Avoid mixing double*real additions to real, multiply by double
!     and then putting it back into real... (to avoid roundoff errors)

      do k=k0,l_nk
         km=max(k-1,1)
         w1= (Ver_idz_8%m(k) + Ver_wp_8%m(k))* Cstv_bar1_8
         w2= ((Ver_idz_8%m(k) - Ver_wm_8%m(k))*Ver_onezero(k))* Cstv_bar1_8
         w3=Ver_wpA_8(k)*epsi_8*Cstv_bar1_8
         w4=Ver_wmA_8(k)*epsi_8*Ver_onezero(k)*Cstv_bar1_8
         do j= j0, jn
            do i= i0, in
               F_rc(i,j,k) = F_rc(i,j,k) - Cstv_bar0_8 * F_fis(i,j)  &
                             - w1 * F_rt(i,j,k) + w2 * F_rt(i,j,km)  &
                             - w3 * F_rf(i,j,k) - w4 * F_rf(i,j,km)
            end do
         end do
      end do

      if (Schm_opentop_L) then

!        Apply opentop boundary conditions

         w1=Cstv_invT_8/Cstv_Tstr_8

         do j= j0, jn
            do i= i0, in
               F_rb(i,j) = F_rt(i,j,k0t) - Ver_ikt_8*(F_rb(i,j) - w1*(F_nest_t(i,j,k0t)-Cstv_Tstr_8))
               F_rc(i,j,k0  ) = F_rc(i,j,k0  ) - Ver_cstp_8 * F_rb(i,j)
            end do
         end do

      end if

!     Apply lower boundary conditions
!
      w1 = Cstv_invT_8*Cstv_invT_m_8 / ( Rgasd_8 * Cstv_Tstr_8 )

      do j= j0, jn
         do i= i0, in
            F_rt(i,j,l_nk) = F_rt(i,j,l_nk) - w1 * F_fis(i,j)
            F_rt(i,j,l_nk) = Ver_wpstar_8(l_nk) * F_rt(i,j,l_nk)
            F_rc(i,j,l_nk) = F_rc(i,j,l_nk) + Ver_cssp_8 * F_rt(i,j,l_nk)
         end do
      end do

1000  format(3X,'UPDATE  THE RIGHT-HAND-SIDES: (S/R PRE)')
!
!     ---------------------------------------------------------------
!
      return
      end
