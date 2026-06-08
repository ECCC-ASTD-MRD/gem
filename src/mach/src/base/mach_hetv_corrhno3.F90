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
! Fichier/File   : mach_hetv_corrhno3.ftn90
! Creation       : P. Makar,  V. Bouchet,  A. Nenes,  S. Gravel,  B. Pabla,  S. Menard
! Description    : Heterogeneous chemistry solver for hno3 correction.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length ne (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using the HNO3 algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                        tats <  2.0, All relative humidities
!
!                  This subroutine is called as a correction on earlier calculations,
!                  and determines the extent to which HNO3g may dissolve into the
!                  acid aerosols (if they have sufficient liquid water).
!
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!
!                     HNO3eq <=> Heq + NO3eq , kNO3
!                  The system of equations to be solved in this case are as follows:
!
!                     kHNO3 (HNO3eq) - (Heq)(SO4eq) = 0,
!                     Heq = Hbefore + del
!                     NO3eq = NO3before + del
!                     HNO3eq = HNO3before - del
!
!                  ...where "del" represents the moles/kg H2O of HNO3 that dissolves to the
!                 aerosol phase, "eq" represents the equilibrium solution, and "before"
!                  represents the starting values before the correction is made.  Note
!                  that NO3before is assumed to be zero in the actual implementation.
!
!                  The solution to the system of equations:
!
!                     del = 0.5 * ( - (Hbefore + kHNO3) +
!                           sqrt( (Hbefore + khno3)^2 + 4 (kHNO3 * HNO3before) ) )
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_corrhno3(npts, nr, ne, nc, so4_i, no3_i, nh4_i, hso4_i, hno3_i,  &
                              h_i, lwn_i, t_i, rh_i, rho_i, k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_water, mach_hetv_activity
   use mach_hetv_mod,         only: lwmin, tstd, rg
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne
   integer(kind=4), intent   (in) :: ncas(npts)
   real(kind=8),    intent   (in) :: k0(nr)
   real(kind=8),    intent   (in) :: p1(nr)
   real(kind=8),    intent   (in) :: p2(nr)
   real(kind=8),    intent(inout) :: so4_i (npts)
   real(kind=8),    intent(inout) :: no3_i (npts)
   real(kind=8),    intent(inout) :: nh4_i (npts)
   real(kind=8),    intent(inout) :: hso4_i(npts)
   real(kind=8),    intent(inout) :: hno3_i(npts)
   real(kind=8),    intent(inout) :: h_i   (npts)
   real(kind=8),    intent(inout) :: lwn_i (npts)
   real(kind=8),    intent   (in) :: t_i   (npts)
   real(kind=8),    intent   (in) :: rh_i  (npts)
   real(kind=8),    intent   (in) :: rho_i (npts)
!!if_off
!
!  Local variables:
!
   integer(kind=4) :: i

   real(kind=8), dimension(ne) :: so4, no3, nh4, hso4, hno3, h, lwn
   real(kind=8), dimension(ne) :: t, aw, rho
!
!  Workspace arrays:
   real(kind=8), dimension(ne) :: del
   real(kind=8), dimension(ne) :: g_h_no3, g_h_hso4
   real(kind=8), dimension(ne) :: g_nh4_no3, g_h2_so4
   real(kind=8), dimension(ne) :: lwo, awu, law
   real(kind=8), dimension(ne) :: so4b, no3b, nh4b, hso4b, hb
!  arrays for reaction rates:
   real(kind=8), dimension(ne) :: khno3
!  Local variables that may be used in the solution:
   real(kind=8) :: work1, work2, lwo_lwn

   so4  = pack(so4_i,  ncas == nc)
   nh4  = pack(nh4_i,  ncas == nc)
   no3  = pack(no3_i,  ncas == nc)
   hso4 = pack(hso4_i, ncas == nc)
   hno3 = pack(hno3_i, ncas == nc)
   h    = pack(h_i,    ncas == nc)
!   lwn  = pack(lwn_i,  ncas == nc)
   aw   = pack(rh_i,   ncas == nc)
   t    = pack(t_i,    ncas == nc)
   rho  = pack(rho_i,  ncas == nc)

!   Initial guess for liquid water content (will be replaced in calculations to follow:
   do i = 1, ne
      lwo(i) = 1.d0
!
!  polynomials for electrolytes only good between aw=0.002 and aw=0.98
      awu(i) = min(max(aw(i), 2.0d-03), 0.98d0)
      law(i) = dble(log(real(awu(i))))    !double precision log unnec.
   end do

!  No iteration required; single calculation.
!  Calculate liquid water, based on current ion concentration.
   call mach_hetv_water(ne, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)

!  (2)  Calculate new concentrations of ions and
!  other phases; multiply old concentration by
!  ratio of old to new liquid water contents
   do i = 1, ne
      lwo_lwn = lwo(i) / max(lwn(i), lwmin)
      hno3(i) = max(hno3(i) * lwo_lwn, lwmin)
      so4(i)  = max(so4(i)  * lwo_lwn, lwmin)
      h(i)    = max(h(i)    * lwo_lwn, lwmin)
      no3(i)  = max(no3(i)  * lwo_lwn, lwmin)
      nh4(i)  = max(nh4(i)  * lwo_lwn, lwmin)
      hso4(i) = max(hso4(i) * lwo_lwn, lwmin)
!
!  Update liquid water content
      lwo(i) = lwn(i)
   end do

!  Prior to calculating activity coefficients, check to see if the aerosol is
!  dry (in which case lwn will be <= the minimum value of 1D-9).
!  For the activity coefficient calculation which follows, these "low water"
!  cases have their molar values replaced by a single bogus set of ion
!  concentrations.  This is to allow the activity coefficient calculation to
!  proceed without division by zero problems.  Later in this routine, the
!  HNO3 calculation is only performed for those gridpoints with liquid water contents
!  greater than lwmin; the use of bogus values here does not affect later code.

   do i = 1, ne
      if(lwn(i) <= lwmin) then
         so4b(i)  = 2.D0
         nh4b(i)  = 3.D0
         hb(i)    = 0.D0
         hso4b(i) = 0.D0
         no3b(i)  = 0.D0
      else
         so4b(i)  = so4(i)
         nh4b(i)  = nh4(i)
         hb(i)    = h(i)
         hso4b(i) = hso4(i)
         no3b(i)  = no3(i)
      end if
   end do

!  Calculate activity coefficients
   call mach_hetv_activity(ne, t, hb, hso4b, nh4b, no3b, so4b, g_h_hso4, &
                           g_h2_so4, g_h_no3, g_nh4_no3)

!  Factor activity coefficients into equilibrium constants
!  to simplify calculations.
!  *  Conversion factors:
!  *
!  *        atmospheres -> moles/kg H2O:
!  *
!  *   atmospheres *1.01325E5(kg m^-1 s^-2 /atm) /
!  *    { 8.3144 (kg m^2 s^-2 mol^-1 K^-1) * T (K) * Lw (kg H2O/kg air)
!  *          * rho (kg air/m3 air)  }    =  moles/(kg H2O)
!  *
!  *        atmospheres^-1 -> (moles/kg H2O)^-1 :
!  *
!  *         1          8.3144 * T * Lw * rho     kg H2O
!  *   ------------- * ----------------------  = --------
!  *   atmospheres      1.01325E5                 moles

!  Calculate rates for Equilibrium constants
   do i = 1, ne
      work1 = tstd / t(i) - 1.D0
      work2 = 1.d0 + log(tstd / t(i)) - tstd / t(i)
!  1:  HNO3(g_eq) <=> HNO3(aq)
!  2:  HNO3(aq) <=> H+ + NO3-
!  Net  HNO3(g_eq) <=> H+ + NO3- : khno3=k1*k2
      khno3(i) = (k0(1) * exp(p1(1) * work1 + p2(1) * work2)) *  &
                 (k0(2) * exp(p1(2) * work1 + p2(2) * work2))
! 1, 2:  HNO3(g_eq) <=> HNO3(aq) <=> H+ + NO3-
! k1: Original units (moles HNO3(aq)/kg H2O)/(atm HNO3(g)), but
! gas is already in moles/(kg H2O); multiply by (R T Lwn rho).
! k2: Divide by the square of the HNO3 activity coefficient
      khno3(i) = khno3(i) * rg * t(i) * lwn(i) * rho(i) / (g_h_no3(i) * g_h_no3(i))
   end do

!   Solve system of equations.
   do i = 1, ne
      del(i) = 0.0d0
      if (lwn(i) > lwmin) then
         del(i) = 0.5d0 * ( - (h(i) + khno3(i)) + sqrt((h(i) + khno3(i)) * (h(i) + khno3(i)) &
                 + 4.d0 * khno3(i) * hno3(i)) )
         del(i) = min(hno3(i), del(i))
         del(i) = max(del(i), 0.0d0)
      end if
   end do
   do i = 1, ne
      hno3(i) = hno3(i) - del(i)
      h(i)    = h(i) + del(i)
      no3(i)  = del(i)
   end do

!  Convert values back to moles/kg air from moles/kg water by multiplying by the liquid water content:
   do i = 1, ne
      hno3(i)  = hno3(i)  * lwn(i)
      so4(i)   = so4(i)   * lwn(i)
      h(i)     = h(i)     * lwn(i)
      no3(i)   = no3(i)   * lwn(i)
      nh4(i)   = nh4(i)   * lwn(i)
      hso4(i)  = hso4(i)  * lwn(i)
   end do

   so4_i  = unpack(so4,  ncas == nc, so4_i)
   nh4_i  = unpack(nh4,  ncas == nc, nh4_i)
   no3_i  = unpack(no3,  ncas == nc, no3_i)
   hso4_i = unpack(hso4, ncas == nc, hso4_i)
   hno3_i = unpack(hno3, ncas == nc, hno3_i)
   h_i    = unpack(h,    ncas == nc, h_i)
   lwn_i  = unpack(lwn,  ncas == nc, lwn_i)

   return
end subroutine mach_hetv_corrhno3
