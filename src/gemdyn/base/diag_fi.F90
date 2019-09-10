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

!**s/r diag_fi - Computes geopotential by vertically integrating the temperature

      subroutine diag_fi ( F_fi, F_s, F_t, F_q, &
                           Minx,Maxx,Miny,Maxy, Nk, i0,in,j0,jn )
      use gmm_geof
      use dynkernel_options
      use tdpack
      use glb_ld
      use cstv
      use ver
      use gmm_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy,Nk,i0,in,j0,jn
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1), intent(out) ::  F_fi
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1), intent(in) :: F_q
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_s
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_t


!author: Andre Plante july 2006.

!arguments
!  Name                        Description
!------------------------------------------------------------
! F_fi          - geopotential
! F_s           - log(pis/pistars)
! F_t           - temperature
! F_fis         - surface geopotential

      integer :: i,j,k,km,istat
      real(kind=REAL64), parameter :: one = 1.d0, half = 0.5d0
      real(kind=REAL64) :: qbar, w1
!
!     ---------------------------------------------------------------

      istat = gmm_get (gmmk_fis0_s,fis0)
      istat = gmm_get (gmmk_sls_s ,sls )

!$omp parallel private(qbar,w1,km)

!$omp do
      do j=j0,jn
         do i=i0,in
            F_fi(i,j,G_nk+1)= fis0(i,j)
         end do
      end do
!$omp enddo

      if (Dynamics_hydro_L) then

!$omp do
      do j=j0,jn
         do k= G_nk,1,-1
            w1= rgasd_8*Ver_dz_8%t(k)
            do i= i0,in
               F_fi(i,j,k)= F_fi(i,j,k+1)+w1*F_t(i,j,k)* &
                   (one+Ver_dbdz_8%t(k)*(F_s(i,j)+Cstv_Sstar_8)+Ver_dcdz_8%t(k)*(sls(i,j)+Cstv_Sstar_8))
            end do
         end do
      end do
!$omp enddo

      else

!$omp do
      do j=j0,jn
         do k= G_nk,1,-1
            km= max(1,k-1)
            w1= rgasd_8*Ver_dz_8%t(k)
            do i= i0,in
               qbar=Ver_wpstar_8(k)*F_q(i,j,k+1)+Ver_wmstar_8(k)*half*(F_q(i,j,k)+F_q(i,j,km))
               qbar=Ver_wp_8%t(k)*qbar+Ver_wm_8%t(k)*F_q(i,j,k)*Ver_onezero(k)
               F_fi(i,j,k)= F_fi(i,j,k+1)+w1*F_t(i,j,k)*exp(-qbar)* &
                   (one+Ver_dbdz_8%t(k)*(F_s(i,j)+Cstv_Sstar_8)+Ver_dcdz_8%t(k)*(sls(i,j)+Cstv_Sstar_8))
            end do
         end do
      end do
!$omp enddo

      end if
!$omp end parallel
!
!     ---------------------------------------------------------------
!
      return
      end
