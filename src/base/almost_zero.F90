!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

logical function almost_zero(A) result(retval)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   !@object Check if value A is close to zero, up to machine precision
   !@author Stéphane Gaudreault -- June 2014
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"

   real(REAL64), intent(in) :: A

   integer(INT64) :: aBit
   integer(INT64), parameter :: two_complement = int(Z'80000000', kind=IDOUBLE)
   aBit = 0

   aBit = transfer(A, aBit)
   if (aBit < 0) then
      aBit = two_complement - aBit
   end if
   ! lexicographic order test with a tolerance of 1 adjacent float
   retval = (abs(aBit) <= 1)
end function almost_zero
