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
! Fichier/File   : mach_hetv_case4.ftn90
! Creation       : P. Makar, V. Bouchet, A. Nenes, S. Gravel, B. Pabla, S. Menard
! Description    : Heterogeneous chemistry solver for case4.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne4 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case4's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne4.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                        1 <= TA/TS < 1.5, mdrh_leto_ambis <= RH < drh_ambis
!                  This case uses linear interpolation between the solutions for the
!                  systems of equations for "dry" and "wet" conditions.
!                  The dry case and conditions are:
!                     case 3, 1 <= TA/TS < 1.5, rn < mdrh_leto_ambis
!                  The wet case and wet case conditions are:
!                     case 5, 1 <= TA/TS < 1.5, drh_ambis <= RH < drh_leto
!                  Interpolation is based on the parameter WF1; the linear interpolant
!                  between wet and dry systems.  When WF1 = 1, the dry system is the
!                  solution, when WF1 = 0, the wet system is the solution.  Values of
!                  WF between these limits indicate an interpolation between the solutions
!                  for the wet and dry cases. For the conditions noted above,
!                  the value of WF is given by:
!                  wf = ( drh_ambis - rh ) / ( drh_ambis - mdrh_leto_ambis)
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case4(npts, nr, ne4, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, ambis_i, leto_i, mdrh_leto_ambis_i,     &
                           drh_ambis_i, lwn_i, t_i, rh_i, ts_i, ta_i, tn_i,    &
                           k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_case5
   use chm_utils_mod,         only: chm_error_l
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne4
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
   real(kind=8),    intent(inout) :: ambis_i          (npts)
   real(kind=8),    intent(inout) :: leto_i           (npts)
   real(kind=8),    intent(inout) :: lwn_i            (npts)
   real(kind=8),    intent   (in) :: t_i              (npts)
   real(kind=8),    intent   (in) :: rh_i             (npts)
   real(kind=8),    intent   (in) :: ta_i             (npts)
   real(kind=8),    intent   (in) :: ts_i             (npts)
   real(kind=8),    intent   (in) :: tn_i             (npts)
   real(kind=8),    intent   (in) :: mdrh_leto_ambis_i(npts)
   real(kind=8),    intent   (in) :: drh_ambis_i      (npts)
!!if_off
!
!  Local variables:
!
   real(kind=8), dimension(ne4) :: so4, no3, nh4, hso4, hno3, h
   real(kind=8), dimension(ne4) :: nh3, ambis, leto
   real(kind=8), dimension(ne4) :: rh
   real(kind=8), dimension(ne4) :: ta, ts, tn
   real(kind=8), dimension(ne4) :: mdrh_leto_ambis, drh_ambis
!
   integer(kind=4)              :: i
   real(kind=8)                 :: wfmin, wfmax
   real(kind=8), dimension(ne4) :: ambisd, letod, hno3d, nh3d
   real(kind=8), dimension(ne4) :: wf, wf1

!  Calculate wet particle gases and components
!  using Case .
   call mach_hetv_case5(npts, nr, ne4, 4, so4_i, no3_i, nh4_i, hso4_i, &
                        hno3_i, h_i, nh3_i, leto_i, lwn_i, t_i, rh_i,  &
                        ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   if (chm_error_l) return

   so4   = pack(so4_i,  ncas == 4)
   nh4   = pack(nh4_i,  ncas == 4)
   no3   = pack(no3_i,  ncas == 4)
   hso4  = pack(hso4_i, ncas == 4)
   hno3  = pack(hno3_i, ncas == 4)
   h     = pack(h_i,    ncas == 4)
   nh3   = pack(nh3_i,  ncas == 4)
   ambis = pack(ambis_i, ncas == 4)
   leto  = pack(leto_i, ncas == 4)
   mdrh_leto_ambis  = pack(mdrh_leto_ambis_i, ncas == 4)
   drh_ambis        = pack(drh_ambis_i,       ncas == 4)
   ta    = pack(ta_i,   ncas == 4)
   ts    = pack(ts_i,   ncas == 4)
   tn    = pack(tn_i,   ncas == 4)
   rh    = pack(rh_i,   ncas == 4)

!  Calculate dry particle gases and components using (inlined) Case 3 .
   do i = 1, ne4
      ambisd(i) = max(3.0d0 * ts(i) - 2.0d0 * ta(i), 0.0d0)
      letod(i)  = max(ta(i) - ts(i), 0.0d0)
      hno3d(i)  = tn(i)
      nh3d(i)   = 0.0d0
   end do

!  Calculate value of WF.
   wfmin = 2.0d0
   wfmax = -2.0d0
   do i = 1, ne4
      wf(i) = (drh_ambis(i) - rh(i)) / (drh_ambis(i) - mdrh_leto_ambis(i))
      wfmin = min(wfmin, wf(i))
      wfmax = max(wfmax, wf(i))
   end do
   if (wfmax > 1.0d0) then
      write(0, *) '### Error in mach_hetv_case4 ###'
      write(0, *) '# interpolant out of bounds ( > 1) in case4.'
      write(0, *) '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if
   if (wfmin < 0.0d0) then
      write(0, *) '### Error in mach_hetv_case4 ###'
      write(0, *) '# interpolant out of bounds ( < 0) in case4.'
      write(0, *) '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if
   do i = 1, ne4
      wf1(i) = min(max(1.0d0 - wf(i),0.D0),1.D0)
   end do

!  Interpolate between the two cases to get the final result (which
!  replaces the array variables used for the wet case):
   do i = 1, ne4
      h(i)    = wf1(i) * h(i)
      nh4(i)  = wf1(i) * nh4(i)
      no3(i)  = wf1(i) * no3(i)
      hso4(i) = wf1(i) * hso4(i)
      so4(i)  = wf1(i) * so4(i)

!  Solid phases, gases:
      ambis(i) = wf(i) * ambisd(i) + wf1(i) * ambis(i)
      leto(i)  = wf(i) * letod(i)  + wf1(i) * leto(i)
      nh3(i)   = wf(i) * nh3d(i)   + wf1(i) * nh3(i)
      hno3(i)  = wf(i) * hno3d(i)  + wf1(i) * hno3(i)
   end do

   so4_i   = unpack(so4,   ncas == 4, so4_i)
   no3_i   = unpack(no3,   ncas == 4, no3_i)
   nh4_i   = unpack(nh4,   ncas == 4, nh4_i)
   hso4_i  = unpack(hso4,  ncas == 4, hso4_i)
   hno3_i  = unpack(hno3,  ncas == 4, hno3_i)
   h_i     = unpack(h,     ncas == 4, h_i)
   nh3_i   = unpack(nh3,   ncas == 4, nh3_i)
   ambis_i = unpack(ambis, ncas == 4, ambis_i)
   leto_i  = unpack(leto,  ncas == 4, leto_i)

   return
end subroutine mach_hetv_case4
