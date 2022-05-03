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

!**s/r glbsum8 - Do a mpi-bit reproducible global 2D sum on real(kind=REAL64) field

      subroutine glbsum8 (F_sum_8,F_field_8,Minx,Maxx,Miny,Maxy,Nk,&
                             Gminx,Gmaxx,Gminy,Gmaxy)
      use gem_options
      use glb_ld
      use glb_pil
      use HORgrid_options
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
      integer Minx,Maxx,Miny,Maxy,NK,Gminx,Gmaxx,Gminy,Gmaxy
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk),intent(in) :: F_field_8
      real(kind=REAL64), dimension(Nk)                    ,intent(out):: F_sum_8
      real(kind=REAL64), dimension (:,:,:), allocatable :: wk1_8
      real(kind=REAL64), dimension(Nk) :: p_sum_8

      integer err,i,j,k
!
!     ---------------------------------------------------------------
!
         p_sum_8 = 0.0d0
         if (Ptopo_myproc == 0) then
            allocate (wk1_8(G_ni,G_nj,NK))
         else
            allocate (wk1_8(1,1,1))
         end if

         call RPN_COMM_coll (wk1_8,1,G_ni,1,G_nj,G_ni,G_nj,NK,0,0,2, &
                   F_field_8,minx,maxx,miny,maxy,0,0,err)

         if (Ptopo_myproc == 0) then
            p_sum_8 = 0.0d0
            do k=1,NK
            do j=Gminy,Gmaxy
            do i=Gminx,Gmaxx
               p_sum_8(k)=p_sum_8(k) + wk1_8(i,j,k)
            end do
            end do
            end do
           !print *,'p_sum=',p_sum_8
           ! F_sum_8 = p_sum_8
         end if
         F_sum_8 = p_sum_8
         call RPN_COMM_bcast (F_sum_8,NK,"MPI_DOUBLE_PRECISION",0,"grid",err)
         deallocate(wk1_8)

      return
      end
