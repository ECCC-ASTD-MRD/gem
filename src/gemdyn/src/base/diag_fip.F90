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

!**s/r diag_fip - Computes geopotential perturbation

      subroutine diag_fip( F_fip, F_s, F_sl, F_t, F_q, F_fis, &
                           Minx,Maxx,Miny,Maxy, Nk, i0,in,j0,jn )
      use dynkernel_options
      use tdpack
      use glb_ld
      use dyn_fisl_options
      use cstv
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk, i0, in, j0, jn
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1), intent(out) :: F_fip
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1), intent(in) :: F_q
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_s, F_fis, F_sl
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_t

!author
!
! Claude Girard
!
!arguments
!  Name                         Description
!-------------------------------------------------------------
! F_fip          - geopotential perturbation
! F_s            - log(pis/pref)
! F_t            - temperature
! F_fis          - surface geopotential


      integer i,j,k,km
      real(kind=REAL64), parameter :: one = 1.d0, half = .5d0
      real(kind=REAL64)  qbar,w1
!
!     ---------------------------------------------------------------
!
      do j=j0,jn
         do i=i0,in
            F_fip(i,j,G_nk+1)= F_fis(i,j)
         end do
      end do

      if (Dynamics_hydro_L) then

         do j=j0,jn
            do k= G_nk,1,-1
               w1= rgasd_8*Ver_dz_8%t(k)
               do i= i0,in
                  F_fip(i,j,k)= F_fip(i,j,k+1)+w1*(F_t(i,j,k)*(one &
                               +Ver_dbdz_8%t(k)*F_s(i,j)  &
                               +Ver_dcdz_8%t(k)*F_sl(i,j))-Cstv_Tstr_8)
               end do
            end do
         end do

      else

         do j=j0,jn
            do k= G_nk,1,-1
               km= max(1,k-1)
               w1= rgasd_8*Ver_dz_8%t(k)
               do i= i0,in
                  qbar=Ver_wpstar_8(k)*F_q(i,j,k+1)+Ver_wmstar_8(k)*half*(F_q(i,j,k)+F_q(i,j,km))
                  qbar=Ver_wp_8%t(k)*qbar+Ver_wm_8%t(k)*F_q(i,j,k)*Ver_onezero(k)
                  F_fip(i,j,k)= F_fip(i,j,k+1)+w1*(F_t(i,j,k)*exp(-qbar)*(one &
                               +Ver_dbdz_8%t(k)*F_s(i,j)  &
                               +Ver_dcdz_8%t(k)*F_sl(i,j))-Cstv_Tstr_8)
               end do
            end do
         end do

      end if
!
!     ---------------------------------------------------------------
!
      return
      end
