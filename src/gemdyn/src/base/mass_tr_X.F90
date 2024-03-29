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

!**s/r mass_tr_X - Evaluate Mass of Component (%w or %m) for all tracers using Bermejo-Conde
!                  (assuming in Mixing Ratio)

      subroutine mass_tr_X (F_mass_X_8,F_bc,F_air_mass,F_minx,F_maxx,F_miny,F_maxy,F_nk, &
                            F_i0,F_in,F_j0,F_jn,F_k0,F_ntr_bc,F_kind_X)

      use adz_mem
      use dynkernel_options
      use omp_timing
      use geomh
      use HORgrid_options
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,           intent(in)                       :: F_ntr_bc                    !Number of tracers using B-C
      type(bc),          intent(in),  dimension(F_ntr_bc) :: F_bc                        !Pointers for tracers using B-C
      integer,           intent(in)                       :: F_kind_X                    !1=Component %w / 2=Component %m
      real(kind=REAL64), intent(out), dimension(F_ntr_bc) :: F_mass_X_8                  !Mass of Component
      integer,           intent(in)                       :: F_minx,F_maxx,F_miny,F_maxy !Dimension H
      integer,           intent(in)                       :: F_nk                        !Number of vertical levels
      integer,           intent(in)                       :: F_i0,F_in,F_j0,F_jn,F_k0    !Scope of operator
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_air_mass        !Air mass

      !object
      !===============================================================================
      !     Evaluate Mass of Component (%w or %m) for all tracers using Bermejo-Conde
      !     (assuming in Mixing Ratio)
      !===============================================================================

      include 'mpif.h'
      include 'rpn_comm.inc'
      integer :: n,i,j,k,err,comm
      real(kind=REAL64), dimension(F_ntr_bc) :: c_mass_8,gc_mass_8
      real, pointer, dimension(:,:,:) :: F_tr_X
      real(kind=REAL64) :: gathV(F_ntr_bc,Ptopo_numproc*Ptopo_ncolors)
!
!---------------------------------------------------------------------
!
!      call gtmg_start (15, 'MASS__', 74)

      c_mass_8 = 0.0d0

!      call gtmg_start (18, 'SOMME_', 15)

      do n=1,F_ntr_bc

         if (F_kind_X==1) F_tr_X => F_bc(n)%w !Component %w
         if (F_kind_X==2) F_tr_X => F_bc(n)%m !Component %m

         !Evaluate Local Mass
         !-------------------
         if (Schm_autobar_L) then

            do j=F_j0,F_jn
               do i=F_i0,F_in
                  c_mass_8(n) = c_mass_8(n) + F_tr_X(i,j,1) * geomh_area_mask_8(i,j)
               end do
            end do

         else

            do k=F_k0,F_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     c_mass_8(n) = c_mass_8(n) + F_tr_X(i,j,k) * F_air_mass(i,j,k) * geomh_mask_8(i,j)
                  end do
               end do
            end do

         end if

      end do

!      call gtmg_stop  (18)

!      call gtmg_start (19, 'REDUCE', 15)

      comm = RPN_COMM_comm ('MULTIGRID')

      !Evaluate Global Mass
      !----------------------------------------
      call MPI_Allgather(c_mass_8,F_ntr_bc,MPI_DOUBLE_PRECISION,gathV,F_ntr_bc,MPI_DOUBLE_PRECISION,comm,err)

      do i=1,F_ntr_bc

         gc_mass_8(i) = sum(gathV(i,:))

      end do

      F_mass_X_8 = gc_mass_8

!      call gtmg_stop (19)

!      call gtmg_stop (15)
!
!---------------------------------------------------------------------
!
      return
      end subroutine mass_tr_X
