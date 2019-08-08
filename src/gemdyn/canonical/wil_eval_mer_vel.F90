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

!**s/r wil_eval_mer_vel - Evaluates the meridional velocity amplitude at a given latitude point and
!                         a given wave-amplitude (MATSUNO)
!                         See Eq.(6a) in Shamir et al.,2019,GMD,12,2181-2193

      function wil_eval_mer_vel (lat,Lamb,n,amp) result (velocity)

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      complex(kind=REAL64)           :: velocity !Meridional velocity amplitude

      real(kind=REAL64), intent(in)  :: lat,  &  !Latitude (radians)
                                        Lamb, &  !Lamb's parameter ~ 2935 for Earth's parameters
                                        amp      !Wave amplitude (m/sec)

      integer,           intent (in) :: n        !Wave-mode (dimensionless)

      !object
      !===============================================================================
      !     Evaluates the meridional velocity amplitude at a given latitude point and
      !     a given wave-amplitude (MATSUNO)
      !     See Eq.(6a) in Shamir et al.,2019,GMD,12,2181-2193
      !     --------------------------------------------------------------------------
      !     Extracted from https://zenodo.org/record/2605203#.XPe8_EO1uEJ
      !===============================================================================

      real(kind=REAL64) :: wil_eval_hermite_pol,y,ex

      external :: wil_eval_hermite_pol
!
!---------------------------------------------------------------------
!
      !Re-scale latitude
      !-----------------
      y = Lamb**0.25 * lat

      !Gaussian envelope
      !-----------------
      ex = exp(-0.5 * y**2)

      velocity = amp * ex * wil_eval_hermite_pol (y,n)
!
!---------------------------------------------------------------------
!
      return
      end function
