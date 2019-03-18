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

!**s/r mass_tr - Evaluate Mass of Tracer (assuming in Mixing Ratio)

      subroutine mass_tr (F_mass_tracer_8,F_tracer,F_air_mass,F_minx,F_maxx,F_miny,F_maxy,F_nk, &
                          F_i0,F_in,F_j0,F_jn,F_k0)

      use HORgrid_options
      use geomh

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      real*8,            intent(out):: F_mass_tracer_8                            !Mass of Tracer
      integer,           intent(in) :: F_minx,F_maxx,F_miny,F_maxy                !Dimension H
      integer,           intent(in) :: F_nk                                       !Number of vertical levels
      integer,           intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0                   !Scope of operator
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_tracer   !Current Tracer (Mixing Ratio)
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_air_mass !Air mass

      !object
      !=======================================================
      !     Evaluate Mass of Tracer (assuming in Mixing Ratio)
      !=======================================================

      integer :: i,j,k,err
      real*8  :: c_mass_8,c_level_8(F_nk),gc_mass_8
      character(len= 9) :: communicate_S

      !-------------------------------------------------------

      !Evaluate Local Mass
      !-------------------
!$omp parallel do private(k,i,j) shared(c_level_8)
      do k=F_k0,F_nk
         c_level_8(k) = 0.0d0
         do j=F_j0,F_jn
         do i=F_i0,F_in
            c_level_8(k) = c_level_8(k) + F_tracer(i,j,k) * F_air_mass(i,j,k) * geomh_mask_8(i,j)
         end do
         end do
      end do
!$omp end parallel do

      c_mass_8 = 0.0d0

      do k=F_k0,F_nk
         c_mass_8 = c_mass_8 + c_level_8(k)
      end do

      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"

      !Evaluate Global Mass using MPI_ALLREDUCE
      !----------------------------------------
      call rpn_comm_ALLREDUCE (c_mass_8,gc_mass_8,1,"MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      F_mass_tracer_8 = gc_mass_8

      return

      end subroutine mass_tr
