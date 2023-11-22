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

!**s/r wil_galewski_geo_8 - Initialize geopotential in balanced state of Galewski's case

      function wil_galewski_geo_8 (evaluated_at_lat_8)

      use tdpack

      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64) wil_galewski_geo_8,evaluated_at_lat_8

      !object
      !=================================================================
      !     Initialize geopotential in balanced state of Galewski's case
      !     Tellus, 56A, 429-440, 2004
      !=================================================================

      !---------------------------------------------------------------

      integer jj,j_partial
      real(kind=REAL64)   wil_galewski_wind_8,intr_8,wind_8,fonc_8,snlat_8,cslat_8,lat_8
      external wil_galewski_wind_8
      integer, parameter :: J_MAX = 100
      real bound1,bound2,lat_4, &
           gauss_latitude(J_MAX),gauss_weight(J_MAX)

      !---------------------------------------------------------------

      !------------------------------------------------------
      !Thini is the height in gradient wind balance with uini
      !- ( af + u*tan(theta) )*u = dthini/d(theta)
      !------------------------------------------------------

      j_partial = J_MAX

      bound1 = -pi_8/2.
      bound2 = evaluated_at_lat_8

      !Evaluate parameters for Gaussian integral
      !-----------------------------------------
      call wil_gauleg (bound1,bound2,gauss_latitude,gauss_weight,j_partial)

      intr_8 = 0.

      !Do Gaussian integral w.r.t theta
      !--------------------------------
      do jj = 1,j_partial

         lat_4 =      gauss_latitude(jj)
         lat_8 = dble(gauss_latitude(jj))

         snlat_8 = sin(lat_8)
         cslat_8 = cos(lat_8)

         wind_8 = wil_galewski_wind_8 (lat_8,1)

         fonc_8 = 0.

         if (cslat_8/=0.) fonc_8 = - ( 2.*omega_8*snlat_8*rayt_8 +  &
                                   wind_8*snlat_8/cslat_8 ) * wind_8

          intr_8 = intr_8 + gauss_weight(jj) * fonc_8

      enddo

      wil_galewski_geo_8 = intr_8

      !---------------------------------------------------------------

      return
      end
