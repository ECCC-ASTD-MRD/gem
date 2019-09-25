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

!**s/r diag_mu calculate mu: ratio of vertical to gravitational accelerations
!
      subroutine diag_mu( F_mu, F_q, F_s, F_sl, &
                          Minx,Maxx,Miny,Maxy, Nk, i0,in,j0,jn )
      use cstv
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy,i0,in,j0,jn,Nk
      real F_mu(Minx:Maxx,Miny:Maxy,Nk), F_q(Minx:Maxx,Miny:Maxy,Nk+1)
      real  F_s(Minx:Maxx,Miny:Maxy), F_sl(Minx:Maxx,Miny:Maxy)

!author
!
! claude girard 2013
!
!object
!       see id section
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_mu         O    -
! F_s          I    - log(pi_s/pref)
! F_q          I    - log(p/pi)
! i0,in,j0,jn  I    - index over which computation will be made.


      integer i,j,k,km
      real(kind=REAL64)  w1, qbar
      real(kind=REAL64), parameter :: one = 1.d0, half = .5d0
!
!     ---------------------------------------------------------------
!
      do k=1,Nk
         km=max(1,k-1)
         do j= j0, jn
            do i= i0, in
               w1 = one + Ver_dbdz_8%t(k)*(F_s(i,j)  + Cstv_Sstar_8) &
                        + Ver_dcdz_8%t(k)*(F_sl(i,j) + Cstv_Sstar_8)
               F_mu(i,j,k) = Ver_idz_8%t(k)*(F_q(i,j,k+1)-F_q(i,j,k))/w1
               qbar = Ver_wpstar_8(k)*F_q(i,j,k+1)+Ver_wmstar_8(k)*half*(F_q(i,j,k)+F_q(i,j,km))
               qbar = Ver_wp_8%t(k)*qbar+Ver_wm_8%t(k)*F_q(i,j,k)
               F_mu(i,j,k) = exp(qbar)*(F_mu(i,j,k)+one)-one
            end do
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end
