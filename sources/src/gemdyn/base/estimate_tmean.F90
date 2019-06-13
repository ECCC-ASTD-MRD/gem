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
!
!**s/p estimate_tmean - Estimate mean temperature on Momentum levels

      subroutine estimate_tmean (F_tmean,F_t,Minx,Maxx,Miny,Maxy,F_nk)
      use HORgrid_options
      use geomh
      use glb_ld
      use cstv
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

      !Arguments
      !---------
      integer,                                   intent(in) :: Minx,Maxx,Miny,Maxy !I, Dimension H
      integer,                                   intent(in) :: F_nk                !I, Number of vertical levels
      real, dimension(F_nk),                     intent(out):: F_tmean             !I: Mean temperature levels
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(in) :: F_t                 !I: Temperature

      !----------------------------------------------------------
      integer i,j,k,err
      real(kind=REAL64) avg_8(F_nk+1),g_avg_8(F_nk+1)
      character(len= 9) communicate_S
      real(kind=REAL64), parameter :: ZERO_8 = 0.0d0
      !----------------------------------------------------------

      avg_8 = ZERO_8

      call grid_area_mask (geomh_area_8,geomh_mask_8,l_ni,l_nj)

      !Estimate Local mean on Thermo levels
      !------------------------------------
      do k=1,F_nk+1

         if (k/=F_nk+1) then

            do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e

               avg_8(k) = avg_8(k) + F_t(i,j,k) * geomh_area_8(i,j) * geomh_mask_8(i,j)

            end do
            end do

         else

            do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e

               avg_8(k) = avg_8(k) + geomh_area_8(i,j) * geomh_mask_8(i,j)

            end do
            end do

         end if

      end do

      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"

      !Estimate Global mean on Thermo levels using MPI_ALLREDUCE
      !---------------------------------------------------------
      call rpn_comm_ALLREDUCE (avg_8,g_avg_8,F_nk+1,"MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      do k=1,F_nk
         g_avg_8(k) = g_avg_8(k)/g_avg_8(F_nk+1)
      end do

      !Estimate Global mean on Momentum levels
      !---------------------------------------
      do k=1,F_nk

         F_tmean(k) = g_avg_8(k)

      end do

      return
      end
