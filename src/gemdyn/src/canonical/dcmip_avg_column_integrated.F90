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

!**s/r dcmip_avg_column_integrated - Evaluate Average Column Integrated of F_tr

      subroutine dcmip_avg_column_integrated (F_avg_tr,F_tr,Minx,Maxx,Miny,Maxy,F_nk)

      use glb_ld
      use gmm_pw

      use, intrinsic :: iso_fortran_env
      implicit none

      !Arguments
      !---------
      integer,                                   intent(in) :: Minx,Maxx,Miny,Maxy  !Dimension H
      integer,                                   intent(in) :: F_nk                 !Number of vertical levels
      real, dimension(Minx:Maxx,Miny:Maxy)     , intent(out):: F_avg_tr             !Average Column Integrated of F_tr
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(in) :: F_tr                 !Tracer

      !object
      !===============================================
      !     Evaluate Average Column Integrated of F_tr
      !===============================================

      integer :: i,j,k
      real, pointer, dimension(:,:,:) :: pm
      real(kind=REAL64) :: avg_8,dp_8
!
!---------------------------------------------------------------------
!
      !Obtain Pressure Momentum at TIME M
      !----------------------------------
      pm => pw_pm_plus

      F_avg_tr(:,:) = 0.

      do j=1,l_nj
      do i=1,l_ni
         avg_8 = 0.0d0
         do k=1,F_nk
            dp_8 = (pm(i,j,k+1) - pm(i,j,k))
            avg_8 = avg_8 + dp_8
            F_avg_tr(i,j) = F_avg_tr(i,j) + F_tr(i,j,k) * dp_8
         end do
         F_avg_tr(i,j) = F_avg_tr(i,j)/avg_8
      end do
      end do
!
!---------------------------------------------------------------------
!
      return
      end
