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

!**s/r get_air_mass - Evaluate Air Mass at appropriate time

      subroutine get_air_mass (F_air_mass,F_time,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0)

      implicit none

      !arguments
      !---------
      integer,          intent(in) :: F_time                                       !Time 1/0
      integer,          intent(in) :: F_minx,F_maxx,F_miny,F_maxy                  !Dimension H
      integer,          intent(in) :: F_k0                                         !Scope of operator
      integer,          intent(in) :: F_nk                                         !Number of vertical levels
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_air_mass !Air_mass
!
!---------------------------------------------------------------------
!
!$omp parallel
      call get_air_mass_hlt (F_air_mass, F_time,F_minx,F_maxx,&
                             F_miny,F_maxy,F_nk,F_k0)
!$omp end parallel
!
!---------------------------------------------------------------------
!
      return
      end subroutine get_air_mass
