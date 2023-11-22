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

!**s/r wil_eval_omega - Unpack Earth dictionary and Evaluates the wave frequencies for a given wave-number
!                       and wave-mode (MATSUNO)
!                       For further details see Eqs. (2)-(5) in Shamir et al.,2019,GMD,12,2181-2193

      function wil_eval_omega (k,n) result (omega)

      use dcst
      use tdpack, only: pi_8, grav_8
      use wil_options

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      real(kind=REAL64)   :: omega !Wave frequency in rad/sec for each wave-type (ROSSBY, EIG, WIG)

      integer, intent(in) :: k, &  !Spherical wave-number (dimensionless)
                             n     !Wave-mode (dimensionless)

      !object
      !=======================================================================================
      !     Unpack Earth dictionary and Evaluates the wave frequencies for a given wave-number
      !     and wave-mode (MATSUNO)
      !     ----------------------------------------------------------------------------------
      !     Extracted from https://zenodo.org/record/2605203#.XPe8_EO1uEJ
      !=======================================================================================

      integer :: j

      real(kind=REAL64) :: delta0,delta4,omegaj(3)

      complex(kind=REAL64), parameter :: IMAG = cmplx(0.d0,1.d0)

      complex(kind=REAL64) :: deltaj

      common /Earth_dict/  OMEGA_,G,A,H0
      real(kind=REAL64) :: OMEGA_,G,A,H0
!
!---------------------------------------------------------------------
!
      !Unpack Earth dictionary
      !-----------------------
      OMEGA_ = Dcst_omega_8
      G      = grav_8
      A      = Dcst_rayt_8
      H0     = Williamson_mean_depth_8

      !Evaluate Eq. (4) in the text
      !----------------------------
      omegaj(:) = 0.

      delta0 = 3. * (G * H0 * (k / A)**2 + 2. * OMEGA_ * (G * H0)**0.5 / A * (2 * n + 1))

      delta4 = -54. * OMEGA_ * G * H0 * k / A**2

      do j = 1,3

         deltaj = (delta4**2 - 4. * delta0**3 + 0. * IMAG)**0.5
         deltaj = (0.5 * (delta4 + deltaj))**(1. / 3.)
         deltaj = deltaj * exp(2. * pi_8 * IMAG * j / 3.)

         !Evaluate Eq. (3) in the text
         !----------------------------
         omegaj(j) = real(-1. / 3. * (deltaj + delta0 / deltaj))

      end do

      !As Eq. (5) in the text
      !----------------------
      if (Williamson_wave_type == 'ROSSBY') then

         omega = - minval(abs(omegaj))

      elseif (Williamson_wave_type == 'WIG') then

         omega = minval(omegaj)

      elseif (Williamson_wave_type == 'EIG') then

         omega = maxval(omegaj)

      end if
!
!---------------------------------------------------------------------
!
      return
      end function
