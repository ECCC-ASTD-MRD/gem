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

!**s/r dry_sfc_pressure_8 - Compute dry air surface pressure REAL64
!
      subroutine dry_sfc_pressure_8 (F_drysfcp0_8, presT_8, p0T_8, &
                               Minx,Maxx,Miny,Maxy,Nk,F_k0,F_timelevel_S)
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=1), intent(IN) :: F_timelevel_S
      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,Nk,F_k0
      real(kind=REAL64) F_drysfcp0_8(Minx:Maxx,Miny:Maxy   ),&
                             presT_8(Minx:Maxx,Miny:Maxy,Nk),&
                               p0T_8(Minx:Maxx,Miny:Maxy)
!
!---------------------------------------------------------------------
!
!$omp parallel
      call dry_sfc_pressure_hlt_8 (F_drysfcp0_8, presT_8, p0T_8,&
                                   F_k0,F_timelevel_S)
!$omp end parallel
!
!---------------------------------------------------------------------
!
      return
      end
