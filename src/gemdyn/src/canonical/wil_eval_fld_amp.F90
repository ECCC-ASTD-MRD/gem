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

!**s/r wil_eval_fld_amp - Evaluates the latitude dependent amplitudes at a given latitude point (MATSUNO)

      function wil_eval_fld_amp (lat,k,n,amp,field) result (fld_amp)

      use dcst
      use wil_options
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      complex(kind=REAL64)            :: fld_amp !Amplitude for the required field

      real(kind=REAL64),   intent(in) :: lat     !Latitude (radians)

      integer,             intent(in) :: k, &    !Spherical wave-number (dimensionless)
                                         n       !Wave-mode (dimensionless)

      real(kind=REAL64),   intent(in) :: amp     !Wave amplitude (m/sec)

      character(len=*),    intent(in) :: field   !'phi' for geopotential height
                                                 !  'u' for zonal velocity
                                                 !  'v' for meridional velocity

      !object
      !====================================================================================
      !     Evaluates the latitude dependent amplitudes at a given latitude point (MATSUNO)
      !     Shamir et al.,2019,GMD,12,2181-2193
      !     -------------------------------------------------------------------------------
      !     Extracted from https://zenodo.org/record/2605203#.XPe8_EO1uEJ
      !====================================================================================

      real(kind=REAL64) :: Lamb

      complex(kind=REAL64), parameter :: IMAG = cmplx(0.d0,1.d0)

      complex(kind=REAL64) :: wil_eval_mer_vel, &
                              v_hat,v_hat_plus_1,v_hat_minus_1,u_hat,p_hat

      external :: wil_eval_mer_vel

      common /Earth_dict/  OMEGA_,G,A,H0
      real(kind=REAL64) :: OMEGA_,G,A,H0
!
!---------------------------------------------------------------------
!
      !Lamb's parameter
      !----------------
      Lamb = (2. * OMEGA_ * A)**2 / (G * H0)

      !Evaluate the meridional velocity amp first
      !------------------------------------------
      v_hat = wil_eval_mer_vel (lat,Lamb,n,amp)

      !Evaluate functions for u and phi
      !--------------------------------
      v_hat_plus_1  = wil_eval_mer_vel (lat,Lamb,n + 1,amp)
      v_hat_minus_1 = wil_eval_mer_vel (lat,Lamb,n - 1,amp)

      !Eq. (6a) in the text
      !--------------------
      if (field == 'v') then

         fld_amp = v_hat

      !Eq. (6b) in the text
      !--------------------
      elseif (field == 'u') then

         u_hat = (- ((n + 1) / 2.0)**0.5 * (Williamson_omega_8 / (G * H0)**0.5 + k / A) * &
                 v_hat_plus_1 - ((n) / 2.0)**0.5 * (Williamson_omega_8 / (G * H0)**0.5 -  &
                 k / A) * v_hat_minus_1)

         !Pre-factors
         !-----------
         u_hat = G * H0 * Lamb**0.25 / &
                 (IMAG * A * (Williamson_omega_8**2 - G * H0 * (k / A)**2)) * u_hat

         fld_amp = u_hat

      !Eq. (6c) in the text
      !--------------------
      elseif (field == 'phi') then

         p_hat = (- ((n + 1) / 2.0)**0.5 * (Williamson_omega_8 + (G * H0)**0.5 * k / A) * &
                 v_hat_plus_1 + ((n) / 2.0)**0.5 * (Williamson_omega_8 - (G * H0)**0.5 *  &
                 k / A) * v_hat_minus_1)

         p_hat = G * H0 * Lamb**0.25 / &
                 (IMAG * A * (Williamson_omega_8**2 - G * H0 * (k / A)**2)) * p_hat

         fld_amp = p_hat

      else

         call gem_error (-1,'WIL_EVAL_FLD_AMP','MATSUNO field amplitudes NOT VALID')

      end if
!
!---------------------------------------------------------------------
!
      return
      end function
