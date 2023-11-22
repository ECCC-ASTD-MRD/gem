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

      subroutine mass_tr_hlt (F_mass_tracer_8,F_tracer,F_air_mass,F_minx,F_maxx,F_miny,F_maxy,F_nk, &
                          F_i0,F_in,F_j0,F_jn,F_k0)

      use adz_mem
      use dynkernel_options
      use geomh
      use masshlt
      use omp_lib
      use HORgrid_options
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      real(kind=REAL64), intent(out):: F_mass_tracer_8                            !Mass of Tracer
      integer,           intent(in) :: F_minx,F_maxx,F_miny,F_maxy                !Dimension H
      integer,           intent(in) :: F_nk                                       !Number of vertical levels
      integer,           intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0                   !Scope of operator
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_tracer   !Current Tracer (Mixing Ratio)
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_air_mass !Air mass

      !object
      !=======================================================
      !     Evaluate Mass of Tracer (assuming in Mixing Ratio)
      !=======================================================

      include 'mpif.h'
      include 'rpn_comm.inc'
      integer :: i,j,k,err,comm
      real(kind=REAL64) :: c_mass_8, c_avg_8
      real(kind=REAL64) :: gathS(Ptopo_numproc*Ptopo_ncolors)
!
!---------------------------------------------------------------------
!     
      c_mass_8 = 0.0d0

      !Evaluate Local Mass
      !-------------------
      if (Schm_autobar_L) then
      
!$omp do
         do j=F_j0,F_jn
            do i=F_i0,F_in
               c_mass_8 = c_mass_8 + F_tracer(i,j,1) * geomh_area_mask_8(i,j)
            end do
         end do
!$omp enddo nowait

      else

!$omp do collapse(2)
         do k=F_k0,F_nk
            do j=F_j0,F_jn
               do i=F_i0,F_in
                  c_mass_8 = c_mass_8 + F_tracer(i,j,k) * F_air_mass(i,j,k) * geomh_mask_8(i,j)
               end do
            end do
         end do
!$omp enddo nowait 

      end if
      thread_sum(1,OMP_get_thread_num()) = c_mass_8

      comm = RPN_COMM_comm ('MULTIGRID')
!$OMP BARRIER

      !Evaluate Global Mass
      !----------------------------------------
!$omp single
      c_avg_8=sum(thread_sum(1,:))
      call MPI_Allgather(c_avg_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)

      g_avg_8(1) = sum(gathS) 
!$omp end single

      F_mass_tracer_8 = g_avg_8(1)
!
!---------------------------------------------------------------------
!
      return
      end subroutine mass_tr_hlt
