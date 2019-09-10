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

!**s/r hz_div_damp

      subroutine hz_div_damp ( F_du,F_dv, F_u, F_v, &
                              i0u,inu,j0u,jnu,i0v,inv,j0v,jnv, &
                              Minx,Maxx,Miny,Maxy,Nk )
      use cstv
      use dcst
      use dyn_fisl_options
      use hvdif_options
      use geomh
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,NK), intent (inout) :: F_du, F_dv
      real, dimension(Minx:Maxx,Miny:Maxy,NK), intent (in)    :: F_u, F_v
      integer, intent(in) :: i0u,inu,j0u,jnu,i0v,inv,j0v,jnv
!author
!   Claude Girard
!

      integer :: i,j,k
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: div
      real(kind=REAL64) :: kdiv_damp,kdiv_damp_max
!
!     ---------------------------------------------------------------
!
      kdiv_damp_max=0.25d0*(Dcst_rayt_8*geomh_hx_8)**2/Cstv_dt_8
      kdiv_damp=Hzd_div_damp*kdiv_damp_max/Cstv_bA_m_8

      do k=1,Nk
         do j=j0v,jnv+1
         do i=i0u,inu+1
            div(i,j,k) = (F_u (i,j,k)-F_u (i-1,j,k))*geomh_invDXM_8(j) &
                       + (F_v (i,j,k)*geomh_cyM_8(j)-F_v (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j)
         end do
         end do
      end do

      do k =1, Nk
         do j=j0u,jnu
         do i=i0u,inu
            F_du(i,j,k) = F_du(i,j,k)+kdiv_damp*(div(i+1,j,k)-div(i,j,k))*geomh_invDXMu_8(j)
         end do
         end do
         do j=j0v,jnv
         do i=i0v,inv
            F_dv(i,j,k) = F_dv(i,j,k)+kdiv_damp*(div(i,j+1,k)-div(i,j,k))*geomh_invDYMv_8(j)
         end do
         end do
      end do

!     ---------------------------------------------------------------

!
      return
      end subroutine hz_div_damp
