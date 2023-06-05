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

!**s/r wil_eval_fld - Evaluates the analytic solutions of either the zonal or meridional velocity
!                     or the geopotential height on at given latitude, longitude and time (MATSUNO)

      function wil_eval_fld (lat,lon,time,k,n,amp,field) result (solution)

      use wil_options
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      real(kind=REAL64)            :: solution !Analytic Solution for required field

      real(kind=REAL64),intent(in) :: lat, &   !Latitude (radians)
                                      lon, &   !Longitude(radians)
                                      time,&   !Time (sec)
                                      amp      !Wave amplitude (m/sec)

      integer,          intent(in) :: k, &     !Spherical wave-number (dimensionless)
                                      n        !Wave-mode (dimensionless)

      character(len=*), intent(in) :: field    !'phi' for geopotential height
                                               !  'u' for zonal velocity
                                               !  'v' for meridional velocity

      !object
      !==================================================================================
      !     Evaluates the analytic solutions of either the zonal or meridional velocity
      !     or the geopotential height on at given latitude, longitude and time (MATSUNO)
      !     -----------------------------------------------------------------------------
      !     Shamir et al.,2019,GMD,12,2181-2193
      !     Extracted from https://zenodo.org/record/2605203#.XPe8_EO1uEJ
      !==================================================================================

      complex(kind=REAL64), parameter :: IMAG = cmplx(0.d0,1.d0)

      complex(kind=REAL64) :: wil_eval_fld_amp,f_hat

      external :: wil_eval_fld_amp
!
!---------------------------------------------------------------------
!
      !Evaluate latitude-dependent amplitude of required field
      !-------------------------------------------------------
      f_hat = wil_eval_fld_amp (lat,k,n,amp,field)

      !Adding time and longitude dependence
      !------------------------------------
      solution = real( exp( IMAG * (k * lon - Williamson_omega_8 * time) ) * f_hat )
!
!---------------------------------------------------------------------
!
      return
      end function
