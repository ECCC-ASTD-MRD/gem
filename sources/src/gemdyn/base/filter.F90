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

!**s/r filter - 2-delta x filter

      subroutine filter2 (F_fd, F_lx, F_coef, Minx,Maxx,Miny,Maxy,Nk)
      use gem_options
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy,Nk, F_lx
      real F_fd(Minx:Maxx,Miny:Maxy,Nk),  F_coef


      integer i, j, k, i0, in, j0, jn, n
      real w1(l_minx:l_maxx,l_miny:l_maxy,nk),coefficient
      real*8, parameter :: half=0.5d0
!
!-------------------------------------------------------------------
!
      if ( ( F_coef >= 0.0001 ) .and. ( F_coef <= 0.5 ) ) then

         coefficient = 1. - F_coef
         i0 = 1 + pil_w
         if (l_west ) i0 = i0+1
         in = l_niu - pil_e
         j0 = 1 + pil_s
         if (l_south) j0 = j0+1
         jn = l_njv - pil_n

         do n=1,F_lx

            call rpn_comm_xch_halo (F_fd,l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

            do k=1,Nk
               do j= 1+pil_s, l_nj-pil_n
               do i=i0,in
                  w1(i,j,k)= F_coef* (F_fd(i-1,j,k)+F_fd(i+1,j,k))*half &
                             + coefficient* F_fd (i,j,k)
               end do
               end do
            end do
            if (l_west) then
               do k=1,Nk
                  do j= 1+pil_s, l_nj-pil_n
                     w1(i0-1,j,k) = F_fd(i0-1,j,k)
                  end do
               end do
            end if
            if (l_east) then
               do k=1,Nk
                  do j= 1+pil_s, l_nj-pil_n
                     w1(in+1,j,k) = F_fd(in+1,j,k)
                  end do
               end do
            end if

            call rpn_comm_xch_halo (w1,  l_minx,l_maxx,l_miny,l_maxy, &
               l_ni,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

            do k=1,Nk
               do j=j0,jn
               do i= 1+pil_w, l_ni-pil_e
                  F_fd(i,j,k)= F_coef * (w1(i,j-1,k)+w1(i,j+1,k))*half &
                               + coefficient * w1(i,j,k)
                  end do
               end do
            end do

            if (l_south) then
               do k=1,Nk
                  do i= 1+pil_w, l_ni-pil_e
                     F_fd(i,j0-1,k) = w1(i,j0-1,k)
                  end do
               end do
            end if
            if (l_north) then
               do k=1,Nk
                  do i= 1+pil_w, l_ni-pil_e
                     F_fd(i,jn+1,k) = w1(i,jn+1,k)
                  end do
               end do
            end if

         end do

      end if
!
!-------------------------------------------------------------------
!
      return
      end
