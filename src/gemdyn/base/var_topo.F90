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

!**s/r var_topo - varies topography incrementally from analysis
!                 topography to target topography

      subroutine var_topo (F_topo, F_step, Minx,Maxx,Miny,Maxy)
      use gmm_geof
      use gem_options
      use tdpack
      use glb_ld
      use gmm_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy
      real F_topo (Minx:Maxx,Miny:Maxy), F_step


      integer i,j, gmmstat
      real(kind=REAL64), parameter :: one = 1.0d0
      real(kind=REAL64), parameter :: two = 2.0d0
      real(kind=REAL64)  lt, pio2, f, b
!     __________________________________________________________________
!
      gmmstat = gmm_get(gmmk_topo_low_s , topo_low )
      gmmstat = gmm_get(gmmk_topo_high_s, topo_high)

      if ( .not. Vtopo_L ) then
         F_topo (1:l_ni,1:l_nj) =  topo_high(1:l_ni,1:l_nj)
         return
      end if

      lt  = Vtopo_ndt
      pio2= pi_8 / two

      f= max(0.,min(F_step-Vtopo_start,real(Vtopo_ndt)))

      if ( Vtopo_ndt <= 0 ) then
         b= 1
         if (F_step < Vtopo_start) b= 0
      else
         b = f / lt
      end if

      b = (cos(pio2 * (one-b) ))**2

      do j= 1, l_nj
      do i= 1, l_ni
         F_topo (i,j) = (one - b)*topo_low(i,j) + b*topo_high(i,j)
      end do
      end do
!     __________________________________________________________________
!
      return
      end
