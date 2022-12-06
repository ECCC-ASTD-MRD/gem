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

!**   s/r vzpng_zmean - Calculate add/substract mean zonal wind component
!
      subroutine vspng_zmean( F_pert, F_field, tsum_8, Minx, Maxx, Miny, Maxy, Nk, lmean )
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer Minx,Maxx,Miny,Maxy, Nk
      real F_field(Minx:Maxx,Miny:Maxy,Nk),F_pert(Minx:Maxx,Miny:Maxy,Nk)
      real(kind=REAL64) tsum_8(l_nj,Nk)
      logical lmean
!
!object
!     The subroutine calculates the zonal mean of the zonal
!     of the wind component and adds/substracts it to/from
!     the corresponding component. It is supposed to be used only for
!     non rotated grid.
!
      integer i,j,k
      integer i0,in,j0,jn,k0,kn,njk,err
      real(kind=REAL64) sum_8(l_nj,Nk),gsum_8(l_nj,Nk)
!
!--------------------------------------------------------------------
!
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj
      k0 = 1
      kn = Nk
      njk= l_nj*Nk

      if(lmean) then
      sum_8(:,:)=0.

         do k=k0,kn
            do j=j0,jn
               sum_8(j,k) = 0.
               do i=i0,in
                  sum_8(j,k) = F_field(i,j,k)+sum_8(j,k)
               end do
            end do
         end do

         call rpn_comm_ALLREDUCE ( sum_8, gsum_8, l_nj*Nk, &
                      "MPI_DOUBLE_PRECISION","MPI_SUM","EW",err )

         tsum_8=gsum_8/real(G_ni)

         do k=k0,kn
            do j=j0,jn
               do i=i0,in
                  F_pert(i,j,k) = F_field(i,j,k) - tsum_8(j,k)
               end do
            end do
         end do

      else

         do k=k0,kn
            do j=j0,jn
               do i=i0,in
                  F_pert(i,j,k) = F_field(i,j,k) + tsum_8(j,k)
               end do
            end do
         end do

      end if
!
!----------------------------------------------------------------
!
      return
      end subroutine vspng_zmean
