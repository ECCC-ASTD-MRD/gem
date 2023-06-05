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
! Fichier/File   : mach_hetv_case7.ftn90
! Creation       : P. Makar,  V. Bouchet,  A. Nenes,  S. Gravel,  B. Pabla,  S. Menard
! Description    : Heterogeneous chemistry solver for case7.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne7 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case7's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne7.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                        1.5 <= TA/TS < 2.0, mdrh_leto_amsul <= RH < drhleto
!                  This case uses linear interpolation between the solutions for the
!                  systems of equations for "dry" and "wet" conditions.
!                  The dry case and conditions are:
!                     case 6, 1.5 <= TA/TS < 2.0,  rh < mdrh_leto_amsul
!                  The wet case and wet case conditions are:
!                     case 8, 1.5 <= TA/TS < 2.0, mdrh_leto_amsul <= RH < drhamsul
!                  Interpolation is based on the parameter WF1; the linear interpolant
!                  between wet and dry systems.  When WF1 = 1, the dry system is the
!                  solution, when WF1 = 0, the wet system is the solution.  Values of
!                  WF between these limits indicate an interpolation between the solutions
!                  for the wet and dry cases.
!                  For the conditions noted above, the value of WF is given by:
!                  wf = ( drhleto - rh ) / ( drhleto - mdrh_leto_amsu
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case7(npts, nr, ne7, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, amsul_i, leto_i, mdrh_leto_amsul_i,     &
                           drh_leto_i, lwn_i, t_i, rh_i, ts_i, ta_i, tn_i,     &
                           k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_case8
   use chm_utils_mod,         only: chm_error_l
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne7
   integer(kind=4), intent   (in) :: ncas             (npts)
   real(kind=8),    intent   (in) :: k0               (nr)
   real(kind=8),    intent   (in) :: p1               (nr)
   real(kind=8),    intent   (in) :: p2               (nr)
   real(kind=8),    intent(inout) :: so4_i            (npts)
   real(kind=8),    intent(inout) :: no3_i            (npts)
   real(kind=8),    intent(inout) :: nh4_i            (npts)
   real(kind=8),    intent(inout) :: hso4_i           (npts)
   real(kind=8),    intent(inout) :: hno3_i           (npts)
   real(kind=8),    intent(inout) :: h_i              (npts)
   real(kind=8),    intent(inout) :: nh3_i            (npts)
   real(kind=8),    intent(inout) :: amsul_i          (npts)
   real(kind=8),    intent(inout) :: leto_i           (npts)
   real(kind=8),    intent(inout) :: lwn_i            (npts)
   real(kind=8),    intent   (in) :: t_i              (npts)
   real(kind=8),    intent   (in) :: rh_i             (npts)
   real(kind=8),    intent   (in) :: ta_i             (npts)
   real(kind=8),    intent   (in) :: ts_i             (npts)
   real(kind=8),    intent   (in) :: tn_i             (npts)
   real(kind=8),    intent   (in) :: mdrh_leto_amsul_i(npts)
   real(kind=8),    intent   (in) :: drh_leto_i       (npts)
!!if_off
!
!  Local variables:
!
   real(kind=8), dimension(ne7) :: so4, no3, nh4, hso4, hno3, h
   real(kind=8), dimension(ne7) :: nh3, amsul, leto
   real(kind=8), dimension(ne7) :: rh
   real(kind=8), dimension(ne7) :: ta, ts, tn
   real(kind=8), dimension(ne7) :: mdrh_leto_amsul, drh_leto
!
   integer(kind=4)              :: i
   real(kind=8)                 :: wfmin, wfmax
   real(kind=8), dimension(ne7) :: amsuld, letod, hno3d, nh3d
   real(kind=8), dimension(ne7) :: wf, wf1

!  Calculate wet particle gases and components using Case 8.
   call mach_hetv_case8(npts, nr, ne7, 7, so4_i, no3_i, nh4_i, hso4_i, &
                        hno3_i, h_i, nh3_i, amsul_i, lwn_i, t_i, rh_i, &
                        ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   if (chm_error_l) return

   so4   = pack(so4_i,  ncas == 7)
   nh4   = pack(nh4_i,  ncas == 7)
   no3   = pack(no3_i,  ncas == 7)
   hso4  = pack(hso4_i, ncas == 7)
   hno3  = pack(hno3_i, ncas == 7)
   h     = pack(h_i,    ncas == 7)
   nh3   = pack(nh3_i,  ncas == 7)
   amsul = pack(amsul_i, ncas == 7)
   leto  = pack(leto_i, ncas == 7)
   mdrh_leto_amsul = pack(mdrh_leto_amsul_i, ncas == 7)
   drh_leto        = pack(drh_leto_i,        ncas == 7)
   ta    = pack(ta_i,   ncas == 7)
   ts    = pack(ts_i,   ncas == 7)
   tn    = pack(tn_i,   ncas == 7)
   rh    = pack(rh_i,   ncas == 7)

!  Calculate value of WF.
   wfmin = 2.0d0
   wfmax = -2.0d0
   do i = 1, ne7
      wf(i) = (drh_leto(i) - rh(i)) / (drh_leto(i) - mdrh_leto_amsul(i))
      wfmin = min(wfmin, wf(i))
      wfmax = max(wfmax, wf(i))
   end do
   if (wfmax > 1.0d0) then
      write(0, *) '### Error in mach_hetv_case7 ###'
      write(0, *) '# interpolant out of bounds ( > 1) in case7.'
      write(0, *) '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if
   if (wfmin < 0.0d0) then
      write(0, *) '### Error in mach_hetv_case7 ###'
      write(0, *) '# interpolant out of bounds ( < 0) in case7.'
      write(0, *) '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if
   do i = 1, ne7
      wf1(i) = min(max(1.0d0 - wf(i),0.D0),1.D0)
   end do

!  Calculate dry particle gases and components using (inlined) Case 6.
   do i = 1, ne7
      letod(i)  = max(2.0d0 * ts(i) - ta(i),0.D0)
      amsuld(i) = max(2.0d0 * ta(i) - 3.0d0 * ts(i),0.D0)
      hno3d(i)  = tn(i)
      nh3d(i)   = 0.0d0
   end do

!  Interpolate between the two cases to get the final result (which
!  replaces the array variables used for the wet case):
   do i = 1, ne7
      h(i)    = wf1(i) * h(i)
      nh4(i)  = wf1(i) * nh4(i)
      no3(i)  = wf1(i) * no3(i)
      hso4(i) = wf1(i) * hso4(i)
      so4(i)  = wf1(i) * so4(i)
!  Solid phases, gases:
      amsul(i) = wf(i) * amsuld(i) + wf1(i) * amsul(i)
      leto(i)  = wf(i) * letod(i)  + wf1(i) * leto(i)
      nh3(i)   = wf(i) * nh3d(i)   + wf1(i) * nh3(i)
      hno3(i)  = wf(i) * hno3d(i)  + wf1(i) * hno3(i)
   end do

   so4_i   = unpack(so4,   ncas == 7, so4_i)
   no3_i   = unpack(no3,   ncas == 7, no3_i)
   nh4_i   = unpack(nh4,   ncas == 7, nh4_i)
   hso4_i  = unpack(hso4,  ncas == 7, hso4_i)
   hno3_i  = unpack(hno3,  ncas == 7, hno3_i)
   h_i     = unpack(h,     ncas == 7, h_i)
   nh3_i   = unpack(nh3,   ncas == 7, nh3_i)
   amsul_i = unpack(amsul, ncas == 7, amsul_i)
   leto_i  = unpack(leto,  ncas == 7, leto_i)

   return
end subroutine mach_hetv_case7
