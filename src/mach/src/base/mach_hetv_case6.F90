!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_hetv_case6.ftn90
! Creation       : P. Makar,  V. Bouchet,  A. Nenes,  S. Gravel,  B. Pabla,  S. Menard
! Description    : Heterogeneous chemistry solver for case6.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver,  recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne6 (inclusive)
!                  has been pre-sorted,  and contains the gridpoint data that must
!                  be solved using case 6's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne6.
!
!                  Units on input are molality; moles/kg H2O.
!
! Extra info     : Athanasios Nenes,  Mail Code 210-41,  Dept. of Chemical Engineering,
!                  California Institute of Technology,  Pasadena,  Ca.,  91125,  USA
!                  nenes@its.caltech.edu
!
!                 The conditions for which this case is called are as follows:
!                     1.5 <= tats <  2.0,   rh <  mdrh_leto_ambis  -> calcb1aa_v
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages,  as numbered):
!                     None; solid phase only
!                  The system of equations to be solved in this case are as follows:
!                     TA = 3 (NH4)3H(SO4)2 + 2 (NH4)2SO4,
!                     TS = 2 (NH4)3H(SO4)2 + (NH4)2SO4.
!                  The solution to the system of equations:
!                     (NH4)3H(SO4)2 = 2 TS - TA,
!                     (NH4)2SO4  = 2 TA - 3 TS
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case6(npts, so4_i, no3_i, nh4_i, hso4_i,   &
                           hno3_i, h_i, nh3_i, amsul_i, leto_i, &
                           lwn_i, ts_i, ta_i, tn_i, ncas)
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: amsul_i (npts)
   real(kind=8),    intent(inout) :: leto_i  (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
!!if_off
!
!   Dry particle, with letovicite and ammonium sulphate being the only solids possible.  Solve the system of
!   two equations in two unknowns, place all nitrate into HNO3, zero the remaining terms, and exit.

   where (ncas == 6)
      leto_i  = max((2.0d0 * ts_i - ta_i), 0.0d0)
      amsul_i = max((2.0d0 * ta_i - 3.0d0 * ts_i), 0.d0)
      hno3_i  = tn_i
      nh3_i   = 0.0d0 
      so4_i   = 0.0d0
      hso4_i  = 0.0d0
      no3_i   = 0.0d0
      h_i     = 0.0d0
      nh4_i   = 0.0d0
      lwn_i   = 0.0d0
   end where

   return
end subroutine mach_hetv_case6
