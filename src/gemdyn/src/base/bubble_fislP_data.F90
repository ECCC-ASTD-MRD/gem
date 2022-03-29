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
!**s/r bubble_fislP_data - generates initial condition for Robert's bubble
!                    experiment (Robert 1993 JAS) - FISL pressure coord.
!
      subroutine bubble_fislP_data (F_t, Mminx, Mmaxx, Mminy, Mmaxy, nk)
      use bubble_options
      use gem_options
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Mminx,Mmaxx,Mminy,Mmaxy,nk
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,nk), intent(out) :: F_t

      integer :: i, j, k, ii
      real(kind=REAL64) :: pp, ex, theta, r,rad
!
!     ---------------------------------------------------------------
!
      if(.not.bubble_gaus_L) then

         do k=1,g_nk
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  ii=i+l_i0-1
                  theta = bubble_theta
                  if ( (((ii)-bubble_ictr)**2 +((k)-bubble_kctr)**2) < bubble_rad**2 ) then
                     theta = theta + 0.5d0
                  end if
                  pp = exp(Ver_a_8%t(k))
                  ex = (pp/Cstv_pref_8)**cappa_8
                  F_t(i,j,k) = theta * ex
               end do
            end do
         end do

      else

         ! Gaussian

         do k=1,g_nk
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  ii=i+l_i0-1
                  r = sqrt( (dble((ii)-bubble_ictr)**2   &
                           + dble((k )-bubble_kctr)**2) )&
                           * dble(bubble_dx)
                  rad = bubble_rad*dble(bubble_dx)
                  if ( r <= rad ) then
                     theta = bubble_theta + 0.5d0
                  else
                     theta = bubble_theta + 0.5d0 * &
                             exp( -(r - rad)**2 / 100.d0**2 );
                  end if
                  pp = exp(Ver_a_8%t(k))
                  ex = (pp/Cstv_pref_8)**cappa_8
                  F_t(i,j,k) = theta * ex
               end do
            end do
         end do

      end if
!
!     -----------------------------------------------------------------
!
      return
      end
