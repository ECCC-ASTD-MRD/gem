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

!**s/r hwnd_stag - Staggering/unstaggering of horizontal wind components

      subroutine hwnd_stag2 ( F_du,F_dv, F_su,F_sv, Minx,Maxx,Miny,Maxy,&
                              NK, F_i0,F_in,F_j0,F_jn, F_stag_L )
      use gem_options
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical,intent(IN) :: F_stag_L
      integer,intent(IN) :: Minx,Maxx,Miny,Maxy,NK,F_i0,F_in,F_j0,F_jn
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(IN) :: F_su, F_sv
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(OUT):: F_du, F_dv

      integer i,j,k, i0u,j0u,inu,jnu, i0v,j0v,inv,jnv
      real(kind=REAL64), parameter :: half=0.5d0 , hpo=1.5d0, &
                      alpha1=-1.d0/16.d0 , alpha2=9.d0/16.d0
!
!     ---------------------------------------------------------------
!
! When F_stag_L=.true.  it is assumed that F_su,F_sv are unstaggered wind
!                       components both defined (1:l_ni,1:l_nj)
! When F_stag_L=.false. it is assumed that F_su(1:l_niu,1:l_lnj and
!                       F_sv(1:l_ni,1:l_njv) are staggered wind components
! on output F_du,F_dv will always be defined (1:l_ni,1:l_nj)

      if (F_stag_L) then

         i0u = F_i0 + west
         inu = F_in - east
         j0u = 1   -G_haloy
         jnu = l_nj+G_haloy
         i0v = 1   -G_halox
         inv = l_ni+G_halox
         j0v = F_j0 + south
         jnv = F_jn - north

         call rpn_comm_xch_halo (F_su,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,Nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo (F_sv,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,Nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         do k = 1,Nk
            do j = j0u, jnu
               do i = i0u, inu
                  F_du(i,j,k) =  ( F_su(i-1,j,k) + F_su(i+2,j,k) )*alpha1 &
                               + ( F_su(i  ,j,k) + F_su(i+1,j,k) )*alpha2
               end do
            end do
            if (l_west) then
               do j = j0u, jnu
                  F_du(F_i0,j,k) = (F_su(F_i0,j,k) + F_su(F_i0+1,j,k)) *half
               end do
            end if
            if (l_east) then
               do j = j0u, jnu
                  F_du(F_in  ,j,k) = (F_su(F_in,j,k) + F_su(F_in+1,j,k))  *half
                  F_du(F_in+1,j,k) = hpo*F_su(F_in+1,j,k) - F_su(F_in,j,k)*half
               end do
            end if

            do j = j0v, jnv
               do i = i0v, inv
                  F_dv(i,j,k) =  ( F_sv(i,j-1,k) + F_sv(i,j+2,k) )*alpha1 &
                               + ( F_sv(i,j  ,k) + F_sv(i,j+1,k) )*alpha2
               end do
            end do
            if (l_south) then
               do i = i0v, inv
                  F_dv(i,F_j0,k) = (F_sv(i,F_j0,k) + F_sv(i,F_j0+1,k)) * half
               end do
            end if
            if (l_north) then
               do i = i0v, inv
                  F_dv(i,F_jn  ,k) = (F_sv(i,F_jn,k) + F_sv(i,F_jn+1,k))  *half
                  F_dv(i,F_jn+1,k) = hpo*F_sv(i,F_jn+1,k) - F_sv(i,F_jn,k)*half
               end do
            end if
         end do

      else

         i0u = F_i0 + 2*west
         inu = F_in -   east
         j0u = 1   -G_haloy
         jnu = l_nj+G_haloy
         i0v = 1   -G_halox
         inv = l_ni+G_halox
         j0v = F_j0 + 2*south
         jnv = F_jn -   north
         
         call rpn_comm_xch_halo (F_su,l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj,Nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo (F_sv,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_njv,Nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         do k = 1,Nk
            do j = j0u, jnu
               do i = i0u, inu
                  F_du(i,j,k) =  ( F_su(i-2,j,k) + F_su(i+1,j,k) )*alpha1 &
                               + ( F_su(i-1,j,k) + F_su(i  ,j,k) )*alpha2
               end do
            end do
            if (l_west) then
               do j = j0u, jnu
                  F_du(F_i0+1,j,k) = (F_su(F_i0,j,k) + F_su(F_i0+1,j,k))   *half
                  F_du(F_i0  ,j,k) = hpo*F_su(F_i0,j,k) - F_su(F_i0+1,j,k) *half
               end do
            end if
            if (l_east) then
               do j = j0u, jnu
                  F_du(F_in  ,j,k) = (F_su(F_in,j,k) + F_su(F_in-1,j,k))  *half
                  F_du(F_in+1,j,k) = hpo*F_su(F_in,j,k) - F_su(F_in-1,j,k)*half
               end do
            end if

            do j = j0v, jnv
               do i = i0v, inv
                  F_dv(i,j,k) =  ( F_sv(i,j-2,k) + F_sv(i,j+1,k) )*alpha1 &
                               + ( F_sv(i,j-1,k) + F_sv(i  ,j,k) )*alpha2
               end do
            end do
            if (l_south) then
               do i = i0v, inv
                  F_dv(i,F_j0+1,k) = (F_sv(i,F_j0,k) + F_sv(i,F_j0+1,k))   *half
                  F_dv(i,F_j0  ,k) = hpo*F_sv(i,F_j0,k) - F_sv(i,F_j0+1,k) *half
               end do
            end if
            if (l_north) then
               do i = i0v, inv
                  F_dv(i,F_jn  ,k) = (F_sv(i,F_jn,k) + F_sv(i,F_jn+1,k))  *half
                  F_dv(i,F_jn+1,k) = hpo*F_sv(i,F_jn,k) - F_sv(i,F_jn-1,k)*half
               end do
            end if
         end do

      end if
!
!     ---------------------------------------------------------------
!
      return
      end subroutine hwnd_stag2
