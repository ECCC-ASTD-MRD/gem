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

!**s/r wil_galewski_wind_8 - Initialize wind in balanced state of Galewski's case

      function wil_galewski_wind_8 (evaluated_at_lat_8,kind)

      use tdpack

      use, intrinsic :: iso_fortran_env
      implicit none

      integer kind
      real(kind=REAL64) wil_galewski_wind_8,evaluated_at_lat_8

      !object
      !=========================================================
      !     Initialize wind in balanced state of Galewski's case
      !     Tellus, 56A, 429-440, 2004
      !     ----------------------------------------------------
      !     kind = 1 = u wind
      !     kind = 2 = vorticity
      !=========================================================

      !---------------------------------------------------------------

      real(kind=REAL64) umax_8,lat0_8,lat1_8,emag_8,ratio_8,wind_wk_8,vort_wk_8, &
             part1_8,part2_8,part3_8,lat_8

      !---------------------------------------------------------------

      umax_8  = 80.d0
      lat0_8  = pi_8/7.
      lat1_8  = pi_8/2. - lat0_8
      emag_8  = exp(-4./(lat1_8-lat0_8)**2)
       lat_8  = evaluated_at_lat_8

      wind_wk_8 = 0.
      vort_wk_8 = 0.

      if (lat_8 > lat0_8.and.lat_8 < lat1_8) then

          !Uini is a barotropically unstable state (with jet at 45N)
          !----------------------------------------------------------
          ratio_8   = (lat_8-lat0_8)*(lat_8-lat1_8)
          ratio_8   = 1./ratio_8
          wind_wk_8 = (umax_8/emag_8)*exp(ratio_8)

          !Vorini is the corresponding vorticity = -1/a cos(theta)*d(cos(theta)*uini)/d(theta)
          !-----------------------------------------------------------------------------------
          part1_8   = tan(lat_8)
          part3_8   = (lat_8-lat0_8)*(lat_8-lat1_8)
          part3_8   = part3_8**2
          part2_8   = (lat_8-lat0_8)+(lat_8-lat1_8)
          vort_wk_8 = (wind_wk_8/rayt_8) * ( part1_8 + part2_8/part3_8 )

      endif

      if (kind == 1) wil_galewski_wind_8 = wind_wk_8
      if (kind == 2) wil_galewski_wind_8 = vort_wk_8

      !---------------------------------------------------------------

      return
      end
