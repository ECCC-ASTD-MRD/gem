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
! Fichier/File   : mach_soa_jiang.ftn90
! Creation       : P. Makar, C. Stroud, Feb 2007 for GEM-MACH
! Description    : Calculates SOA yields for the ADOM speciation using
!                  the Jiang approach.
!
! Extra info     : Based on AURAMS chemf.ftn90 code
!                  SOA calculations, gas-phase chemistry solver.
!
! Modified       : A. Akingunola, Nov 2017 - to optimize the code, isolated comonly used calculation
!                  structures that involve reading values from the bus into the individual variables,
!                  and updated constants to be of double precission
!                  V. Savic-Jovcic, June 2018 - to increase readability of the code,
!                  introduced iay_functions and a three-step do loop for Newton's iteration
!                  C. Stroud, October 2018 - added isoprene as a precursor to SOA, based
!                  on yield in Barsanti, K. et al. (2013)
!                  change name from "mach_gas_soa_yield" to "mach_soa_jiang"
!
! Arguments:  IN
!                p2d            -> pressure (pa)
!                tp             -> temperature (k)
!                moi            -> total condensible organics in ug/m3
!                soa_gas_diff   -> dvoc (ppmv)
!
!	      OUT
!	         dsoa           -> dsoa formed from input dvoc (ug/m3)
!
!============================================================================
!
!!if_on
subroutine mach_soa_jiang(moi, dsoa, soa_gas_diff, tp, p2d, gni, gnk)
   use mach_pkg_gas_mod,   only: nsp_soa_gases
!!if_off
   use chm_utils_mod,      only: CHM_MSG_DEBUG, ik !, chm_lun_out
   use chm_consphychm_mod, only: rgasi
   use chm_nml_mod,        only: nk_start, chm_timings_l
   implicit none
!!if_on
   integer(kind=4), intent (in) :: gni, gnk
   real(kind=8),    intent (in) :: soa_gas_diff(nsp_soa_gases, gni, gnk)
   real(kind=8),    intent (in) :: moi         (gni, gnk)
   real(kind=4),    intent(out) :: dsoa        (gni, gnk)
   real(kind=4),    intent (in) :: p2d         (gni, gnk)
   real(kind=4),    intent (in) :: tp          (gni, gnk)
!!if_off
!
!  Declaration of local variables
!
   integer(kind=4)  :: i, k, mm
   real(kind=8)     :: tt_fac1, tt_fac2, tt_fac3
   real(kind=8)     :: alka_var, arom_var, tolu_var, alke_var1, alke_var2, isop_var
   real(kind=8)     :: x, f, df
   real(kind=8)     :: m0, p
   real(kind=8)     :: falka, farom, ftolu, falke1, falke2, falke3, fisop
   real(kind=8)     :: kom1a, kom1b, kom1c, kom1d, kom1e, kom1f, kom1g
   real(kind=8)     :: kom2a, kom2b, kom2c, kom2d, kom2e, kom2f, kom2g
   real(kind=8)     :: dena, denb, denc, dend, dene, denf, deng
   real(kind=8)     :: y1, y2, y3, y4, y5, y6, y7
   real(kind=8)     :: dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7
!  Temperature-independent SOA product parameters (Jiang, W. 2003):
   real(kind=8), parameter :: alpha1a = 0.01d0,   alpha2a = 0.3d0    ! For Alkane + OH (Anthro)
   real(kind=8), parameter :: alpha1b = 0.038d0,  alpha2b = 0.167d0  ! For Aromatic + OH
   real(kind=8), parameter :: alpha1c = 0.071d0,  alpha2c = 0.138d0  ! For Toluene + OH
   real(kind=8), parameter :: alpha1d = 0.01d0,   alpha2d = 0.3d0    ! For Alkene1 + OH (Anthro)
   real(kind=8), parameter :: alpha1e = 0.038d0,  alpha2e = 0.326d0  ! For Alkene2 + OH (Alpha-Pinene)
   real(kind=8), parameter :: alpha1f = 0.13d0,   alpha2f = 0.406d0  ! For Alkene3 + OH (Beta-Pinene)
   real(kind=8), parameter :: alpha1g = 0.0288d0, alpha2g = 0.232d0  ! For Isoprene
!
!  Declaration of external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   call msg_toall(CHM_MSG_DEBUG, 'mach_soa_jiang [BEGIN]')
   if (chm_timings_L) call timing_start_omp(336, 'mach_soa_jiang', 330)

!  Initializations
   dsoa = 0.0
!
!  Calculate the change in concentration for the SOA producing
!  species, multiply by appropriate conversion factors and yields
!  to calculate the total SOA mass added.  Note that in this initial
!  version, the condensable gas produced is assumed to go entirely to
!  aerosol; the initial condensable gas mixing ratio is assumed to
!  be the same as the initial SOA mass mixing ratio.  When CAM is
!  coupled, the initial SOA mas will have to be taken from the
!  sum of the SOA masses in each bin.
!  Note that ALL of the condensable gas is removed in CAM; the addition to the
!  total condensable gas calculated here is that from only the current operator
!  step.
!
   !  Secondary organic aerosol yield loop:
   do k = nk_start, gnk
      do i = 1, gni
!    Add primary and secondary to give m0. Assign poa to p locally in loop.
!    convert existing condensable mass (assumed to be entirely aerosol, here)
!    back to ug / m3 from ug / kg
!    added p=0 for absorption of SOA into entire OC phase
         m0 = moi(i, k)

!  The p=0 line changes the IAY formula so that the SOA can partition into both
!  the SOA and POA. If you want to partition the SOA only into existing SOA
!  then we use the first line commented out.
!         p = poa(i, k) * dble(rho(i, k))
         p = 0.0d0
         if (m0 /= 0.0d0) then
!  Calculate temperature-dependent factors for K dependence, based on the work of
!  Sheehan and Bowman (ES&T 35, no 11, 2129-2135, 2001.  Sources for the enthalpy of
!  vapourization of the different SOA products are given  for each factor.
!  Measurement t* temperatures were 35C (308.15K) for all species but isoprene. For
!  isoprene, measurement t* temperatures were 21.85C (295 K).
!
!  First, evaluate commonly used variables
            tt_fac1 = dble(tp(i, k)) / 308.15d0
            tt_fac2 = 1.0d0 / dble(tp(i, k)) - 1.0d0 / 308.15d0
            tt_fac3 = dble(p2d(i, k)) / dble(tp(i, k))
            alka_var  = max(0.0d0, soa_gas_diff(1, i, k)) * 7.226D-01 * tt_fac3
            arom_var  = max(0.0d0, soa_gas_diff(2, i, k)) * 1.6898D+01 * tt_fac3
            tolu_var  = max(0.0d0, soa_gas_diff(3, i, k)) * 1.14315D+01 * tt_fac3
            alke_var1 = max(0.0d0, soa_gas_diff(4, i, k)) * 5.167D-02 * tt_fac3
            alke_var2 = max(0.0d0, soa_gas_diff(4, i, k)) * 3.810D+00 * tt_fac3
            isop_var  = max(0.0d0, soa_gas_diff(5, i, k)) * 8.193D+00 * tt_fac3
!
!  H for Alkanes:  using C10 monocarboxylic acid extrapolated from Tao and McMurry,
!  ES&T 23, 1519-1523, 1989's formula for even numbered monocarboxylic acids.
            falka = tt_fac1 * exp(136.7d+03 / dble(rgasi) * tt_fac2)
!  H for Aromatics:  using Sheehan and Bowman's value of 17.5Kcal / mol = 73.3 KJ / mol
            farom = tt_fac1 * exp(73.3d+03 / dble(rgasi) * tt_fac2)
!  H for Toluene:  using Sheehan and Bowman's value of 17.5Kcal / mol = 73.3 KJ / mol
            ftolu = farom
!  H for alkenes (anth.) : Tao and McMurry for C10 monocarboxylic acids, again:
            falke1 = falka
!  H for alkenes (alpha-pinene): using Sheehan and Bowman's pinic acid value of
!  22.3 Kcal / mol = 93.35 KJ / mol
            falke2 = tt_fac1 * exp(93.35d+03 / dble(rgasi) * tt_fac2)
!  H for alkenes (other monoterpenes):  use same value as alpha-pinene:
            falke3 = falke2
!  H for isoprene
            fisop = dble(tp(i, k)) / 295.00d0 * exp( 42.0d+03 / dble(rgasi) * (1.0d0 / dble(tp(i, k)) - 1.0d0 / 295.00d0) )
!  Calculate condensable gas mass mixing ratio in ug / kg, CREATED during the current
!  time step.  Note that this replaces the total organic mass.  A three-stage
!  Newton's iteration of the Jiang et al equation is used:

!  Define SOA product parameters that are temperature dependent (Jiang, W. 2003)
!  For Alkane + OH (Anthro)
            kom1a   = 0.02d0 * falka
            kom2a   = 5.0D-04 * falka
!  For Aromatic + OH
            kom1b   = 0.042d0 * farom
            kom2b   = 1.4D-03 * farom
!  For Toluene + OH
            kom1c   = 0.053d0 * ftolu
            kom2c   = 1.9D-03 * ftolu
!  For Alkene1 + OH (Anthro)
            kom1d   = 0.02d0 * falke1
            kom2d   = 5.0D-04 * falke1
!  For Alkene2 + OH (Alpha-Pinene)
            kom1e   = 0.171d0 * falke2
            kom2e   = 4.0D-03 * falke2
!  For Alkene3 + OH (Beta-Pinene)
            kom1f   = 0.044d0 * falke3
            kom2f   = 4.9D-03 * falke3
!  For Isoprene
            kom1g   = 1.620d0 * fisop
            kom2g   = 8.62D-03 * fisop
!  Initial guess for new value of m0 = old value:
            x = m0
!  A three-step do loop for the three-stage Newton's iteration:
            do mm = 1, 3
! Switch for jiang (explicit x and p) or odum soa yield (p = x) formula
!!! Odom switch now in seprate subroutine
!              if (chm_soa_s == 'ODUM') then
!                 p = x   ! this line will convert jiang formula to odum formula
!              endif

!  Define IAY functions for each ROG
               dena = iay_den(alpha1a, kom1a, alpha2a, kom2a, x, p)
               y1 = iay_y(alpha1a, kom1a, alpha2a, kom2a, x, dena)
!
               denb = iay_den(alpha1b, kom1b, alpha2b, kom2b, x, p)
               y2 = iay_y(alpha1b, kom1b, alpha2b, kom2b, x, denb)
!
               denc = iay_den(alpha1c, kom1c, alpha2c, kom2c, x, p)
               y3 = iay_y(alpha1c, kom1c, alpha2c, kom2c, x, denc)
!
               dend = iay_den(alpha1d, kom1d, alpha2d, kom2d, x, p)
               y4 = iay_y(alpha1d, kom1d, alpha2d, kom2d, x, dend)
!
               dene = iay_den(alpha1e, kom1e, alpha2e, kom2e, x, p)
               y5 = iay_y(alpha1e, kom1e, alpha2e, kom2e, x, dene)
!
               denf = iay_den(alpha1f, kom1f, alpha2f, kom2f, x, p)
               y6 = iay_y(alpha1f, kom1f, alpha2f, kom2f, x, denf)
!
               deng = iay_den(alpha1g, kom1g, alpha2g, kom2g, x, p)
               y7 = iay_y(alpha1g, kom1g, alpha2g, kom2g, x, deng)

!  Define derivatives of IAY with respect to x for each ROG
               dydx1 = iay_dydx(alpha1a, kom1a, alpha2a, kom2a, x, p, dena)
!
               dydx2 = iay_dydx(alpha1b, kom1b, alpha2b, kom2b, x, p, denb)
!
               dydx3 = iay_dydx(alpha1c, kom1c, alpha2c, kom2c, x, p, denc)
!
               dydx4 = iay_dydx(alpha1d, kom1d, alpha2d, kom2d, x, p, dend)
!
               dydx5 = iay_dydx(alpha1e, kom1e, alpha2e, kom2e, x, p, dene)
!
               dydx6 = iay_dydx(alpha1f, kom1f, alpha2f, kom2f, x, p, denf)
!
               dydx7 = iay_dydx(alpha1g, kom1g, alpha2g, kom2g, x, p, deng)

!  Define function:
               f = x - m0 - (alka_var  * y1 + arom_var  * y2 + tolu_var  * y3 + &
                             alke_var1 * y4 + alke_var2 * y5 + alke_var2 * y6 + &
                             isop_var  * y7)

!  Define derivative of function:
               df = 1.0d0 - (alka_var  * dydx1 + arom_var  * dydx2 + tolu_var  * dydx3 + &
                             alke_var1 * dydx4 + alke_var2 * dydx5 + alke_var2 * dydx6 + &
                             isop_var  * dydx7)

!  Next approximation of variable:
               x = x - f / df
            end do
!        if(chm_lun_out > 0) write(chm_lun_out, *) '3rd iteration f df x: ', f, df, x
!  'x' is the new value of M0 at the current time.  What's needed for later size apportionment
!  in CAM is the change in M0 values:  Assign this to "dsoa" (ug/m3)
!
            dsoa(i, k) = real(dmax1(x - m0, 0.0d0))

         end if
      end do
   end do

   if (chm_timings_L) call timing_stop_omp(336)

   contains

! Common denominator in IAY function:
   real(kind=8) function iay_den(a1,k1,a2,k2,xx,pp)
   real(kind=8), intent(in) :: a1, k1, a2, k2, xx, pp
        iay_den = a1 * k1 * (k1 * xx**2 + pp) / (1.0d0 + k1 * xx)**2 &
                + a2 * k2 * (k2 * xx**2 + pp) / (1.0d0 + k2 * xx)**2
   end function iay_den

! IAY function:
   real(kind=8) function iay_y(a1,k1,a2,k2,xx,den)
   real(kind=8), intent(in) :: a1, k1, a2, k2, xx, den
        iay_y = (xx * (a1 * k1 / (1.0d0 + k1 * xx) + a2 * k2 / (1.0d0 + k2 * xx)))**2 / den
   end function iay_y

! Derivative of IAY function:
   real(kind=8) function iay_dydx(a1,k1,a2,k2,xx,pp,den)
   real(kind=8), intent(in) :: a1, k1, a2, k2, xx, pp, den
        iay_dydx = 2.0d0 * xx * (a1 * k1 / (1.0d0 + k1 * xx) + a2 * k2 / (1.0d0 + k2 * xx))**2 / den &
                 + 2.0d0 * xx**2 * (a1 * k1 / (1.0d0 + k1 * xx) + a2 * k2 / (1.0d0 + k2 * xx)) * &
                   (-1.0d0 * a1 * k1**2 / (1.0d0 + k1 * xx)**2 - a2 * k2**2 / (1.0d0 + k2 * xx)**2) / den &
                 - xx**2 * (a1 * k1 / (1.0d0 + k1 * xx) + a2 * k2 / (1.0d0 + k2 * xx))**2 * &
                   (2.0d0 * a1 * k1**2 * xx / (1.0d0 + k1 * xx)**2 - 2.0d0 * a1 * k1**2 * (k1 * xx**2 + pp) / (1.0d0 + k1 * xx)**3 &
                    + 2.0d0 * a2 * k2**2 * xx / (1.0d0 + k2 * xx)**2 - 2.0d0 * a2 * k2**2 * (k2 * xx**2 + pp) / (1.0d0 + k2 * xx)**3) / den**2
   end function iay_dydx

end subroutine mach_soa_jiang
