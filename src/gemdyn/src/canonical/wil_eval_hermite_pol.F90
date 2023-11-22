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

!**s/r wil_eval_hermite_pol - Evaluates the normalized Hermite polynomial of degree n at point/s x
!                             using the three-term recurrence relation
!                             For further details see Eq. (7) in Shamir et al.,2019,GMD,12,2181-2193

      recursive function wil_eval_hermite_pol (x,n) result (hermite)

      use tdpack, only: pi_8
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      real(kind=REAL64)            :: hermite !Evaluation of the normalized Hermite polynomial

      real(kind=REAL64),intent(in) :: x       !Point where the evaluation takes place

      integer,          intent(in) :: n       !Degree of the Hermite polynomial

      !object
      !===========================================================================
      !     Evaluates the normalized Hermite polynomial of degree n at point/s x
      !     using the three-term recurrence relation
      !     For further details see Eq. (7) in Shamir et al.,2019,GMD,12,2181-2193
      !     ----------------------------------------------------------------------
      !     Extracted from https://zenodo.org/record/2605203#.XPe8_EO1uEJ
      !===========================================================================
!
!---------------------------------------------------------------------
!
      !Main evaluation
      !---------------
      if (n < 0) then

         hermite = 0.0d0

      elseif (n == 0) then

         hermite = 1.0d0 / pi_8**0.25

      elseif (n == 1) then

         hermite = (4.0 / pi_8)**0.25 * x

      elseif (n >= 2) then

         hermite = ((2.0 / n)**0.5 * x   * wil_eval_hermite_pol(x, n - 1) - &
                    ((n - 1.0) / n)**0.5 * wil_eval_hermite_pol(x, n - 2))

      end if
!
!---------------------------------------------------------------------
!
      return
      end function
