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
! Fichier/File   : mach_hetv_water.ftn90
! Creation       : P. Makar, V. Bouchet, A. Nenes, S. Gravel, B. Pabla, S. Menard
! Description    : Uses the ZSR approach to estimate the water content of an
!                  aerosol.  Simply put, the ZSR approach requires that the
!                  user determine a set of electrolyte moles/kg air, which, when
!                  dissociated, would give rise to the same moles/kg air of
!                  ions as is present in the aerosol.  The ZSR formula then states:
!                    Lw = Sum (i)  { [Moles/(kg air) of the i'th electrolyte] / [Molality
!                  of the i'th electrolyte at the same water activity as the
!                  aerosol solution ] }.  Lw is the liquid water content,
!                  in kg water / kg air.  Unfortunately, more than one combination of
!                  electrolytes is possible, and more than one possible set of ions from each
!                  electrolyte is possible.  Information from Seinfeld's group (Seinfeld,
!                  personal communication, 2000) suggests that the resulting variations in
!                  aerosol water content are usually not significant.  This routine uses
!                  the following methodology for determining the moles electrolyte/ kg air from
!                  the moles ion/kg air:
!                    (1)  Number of moles of letovicite (NH4)3H(SO4)2 is taken from
!                    the smallest of the number of moles of SO4=, 1/3 of the NH4-,
!                    or half of the HSO4-.  The last is to allow some partitioning between
!                    ammonium bisulphate and letovicite for the HSO4- ions.
!
!                    (2)  Number of moles of ammonium bisulphate (NH4HSO4) is taken from
!                    the smallest of the number of moles of the remaining HSO4- and NH4+.
!
!                    (3)  Number of moles of ammonium sulphate (NH4)2SO4 is taken from the
!                    smallest of the number of half of the remaining moles of NH4+, and
!                    the remaining moles of SO4=.
!
!                    (4)  Number of moles of ammonium nitrate NH4NO3 is taken from the
!                    smallest of the number of remaining moles of NH4+ and NO3-.
!
!                    (5)  Number of moles of H2SO4 is taken from the minimum of the
!                    remaining number of moles of SO4= and half of the remaining moles of
!                    H+.
!
!                    (6)  Number of moles of H2SO4 is increased by the minimum of any
!                    remaining HSO4- and the remaining H+.
!
!                    (7)  Number of moles of HNO3 is taken from the minimum of the remaining
!                    number of moles of H+ and NO3-.
!
!                   Note:  the above priority was taken due to the following reasons:
!                    (a)  Ammonium nitrate and HNO3 are not found in solution unless there
!                   is sufficient NH4+ to neutralize the aerosol.  Therefore, the attempt is
!                   made to use up the ammonium on sulphate electrolytes before attempting
!                   to create ammonium nitrate.
!                    (b)  Letovicite and ammonium bisulphate are assumed to be both low in
!                   magnitude compared to ammonium sulphate, and (here, for the purposes
!                   of partitioning the ion moles) the bisulphate ion is not assumed
!                   to partion further into SO4= and H+.  The equilibrium between HSO4- and
!                   SO4=, maintained in the rest of the code, will ensure that HSO4- is
!                   relatively low and hence the water content attributed to letovicite and
!                   ammonium bisulphate will be relatively low.
!                   (c)  H2SO4 is given priority over HNO3 formation, since the latter is
!                   unlikely to form in significant quantities in the presence of the former
!                   (based on gamma values for HNO3 entering aerosols, in DeMore et al., 1997).
!                  Important note:  The ZSR technique, lacking a unique solution to the water
!                  content problem, should be compared to measurements at the first
!                  available opportunity.
!
! Extra info     :  Set minimum possible liquid water content of aerosols:
!                   Note that this value corresponds to 1.5E-10 ppbv H2O on a molar
!                   basis, while the US standard atmosphere's water profile
!                   has lower limit of 0.1 ppmv up into the stratosphere.
!                   A lower limit is necessary to prevent division by zero errors
!                   in the case of very low water concentrations.
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_water(ne, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)
!!if_off

   use mach_hetv_mod,   only: lwmin
   implicit none
!!if_on
   
   integer(kind=4), intent (in) :: ne
   real(kind=8),    intent (in) :: so4 (ne)
   real(kind=8),    intent (in) :: h   (ne)
   real(kind=8),    intent (in) :: no3 (ne)
   real(kind=8),    intent (in) :: nh4 (ne)
   real(kind=8),    intent (in) :: hso4(ne)
   real(kind=8),    intent (in) :: awu (ne)
   real(kind=8),    intent (in) :: law (ne)
   real(kind=8),    intent (in) :: lwo (ne)
   real(kind=8),    intent(out) :: lwn (ne)
!!if_off
!
! Local variables
!
   integer(kind=4)           :: i
   real(kind=8)              :: wleto,  wnitac, wsulac,  wamsul
   real(kind=8)              :: wambis, wamnit, wsulac1, wsulac2
   real(kind=8)              :: wso4w,  whw,    wno3w,   wnh4w
   real(kind=8)              :: whso4w

!----------------------------------------------------------------------
!Name      |Description                           |T  |Dimensions |In/|
!          |                                      |Y  |           |Out|
!          |                                      |P  |           |Loc|
!          |                                      |E  |           |cal|
!----------------------------------------------------------------------
! ambis    | Estimated molality of ammonium       | R | ne        | L |
!          |bisulphate [NH4HSO4], which, when     |   |           |   |
!          |dissolved with the other              |   |           |   |
!          |electrolytes, will help create a      |   |           |   |
!          |solution with the same ionic          |   |           |   |
!          |concentrations as the aerosol         |   |           |   |
!          |solution.                             |   |           |   |
! amnit    | Estimated molality of ammonium       | R | ne        | L |
!          |nitrate [NH4NO3], which, when         |   |           |   |
!          |dissolved with the other              |   |           |   |
!          |electrolytes, will help create a      |   |           |   |
!          |solution with the same ionic          |   |           |   |
!          |concentrations as the aerosol         |   |           |   |
!          |solution.                             |   |           |   |
! amsul    | Estimated molality of ammonium       | R | ne        | L |
!          |sulphate [(NH4)2SO4], which, when     |   |           |   |
!          |dissolved with the other              |   |           |   |
!          |electrolytes, will help create a      |   |           |   |
!          |solution with the same ionic          |   |           |   |
!          |concentrations as the aerosol         |   |           |   |
!          |solution.                             |   |           |   |
! aw       | Water activity (relative humidity on | R | ne        | I |
!          |0 - 1 scale).                         |   |           |   |
! awu      | Water activity used in the           | R | ne        | L |
!          |calculations.  The polynomials for    |   |           |   |
!          |single electrolyte solution           |   |           |   |
!          |molalities as a function of water     |   |           |   |
!          |activity are only accurate between aw |   |           |   |
!          |= 0.002 and aw = 0.98.  The aw values |   |           |   |
!          |read in to the subroutine are clipped |   |           |   |
!          |at these limits and placed into awu,  |   |           |   |
!          |which is then used in water           |   |           |   |
!          |calculations.                         |   |           |   |
! h        | Array of initial H+ ion molalities   | R | ne        | I |
!          |(moles/kg)                            |   |           |   |
! hso4     | Array of initial HSO4- ion           | R | ne        | I |
!          |molalities (moles/kg)                 |   |           |   |
! hso4w    | Working array for HSO4-.  Number of  | R | ne        | L |
!          |moles of HSO4- remaining after each   |   |           |   |
!          |partition.                            |   |           |   |
! hw       | Working array for H+.  Number of     | R | ne        | L |
!          |moles of H+ remaining after each      |   |           |   |
!          |partition.                            |   |           |   |
! i        | Loop index.                          | I |           | L |
! law      | Natural logarithm of aw.  Some of    | R | ne        | L |
!          |the polynomials are in ln(aw) rather  |   |           |   |
!          |than aw, for greater accuracy.        |   |           |   |
! leto     | Estimated molality of letovicite     | R | ne        | L |
!          |[NH4)3H(SO4)2], which, when dissolved |   |           |   |
!          |with the other electrolytes, will     |   |           |   |
!          |help create a solution with the same  |   |           |   |
!          |ionic concentrations as the aerosol   |   |           |   |
!          |solution.                             |   |           |   |
! lwn      | New water content of aerosol (kg H2O | R | ne        | O |
!          |/ kg air)                             |   |           |   |
! lwo      | Old water content of aerosol (kg H2O | R | ne        | I |
!          |/ kg air).                            |   |           |   |
! ne       | Number of points for which water     | I |           | I |
!          |content is required.                  |   |           |   |
! nh4      | Array of initial NH4+ ion molalities | R | ne        | I |
!          |(moles/kg)                            |   |           |   |
! nh4w     | Working array for NH4+.  Number of   | R | ne        | L |
!          |moles of NH4+ remaining after each    |   |           |   |
!          |partition.                            |   |           |   |
! nitac    | Estimated molality of nitric acid    | R | ne        | L |
!          |[HNO3], which, when dissolved with    |   |           |   |
!          |the other electrolytes, will help     |   |           |   |
!          |create a solution with the same ionic |   |           |   |
!          |concentrations as the aerosol         |   |           |   |
!          |solution.                             |   |           |   |
! no3      | Array of initial NO3- ion molalities | R | ne        | I |
!          |(moles/kg)                            |   |           |   |
! no3w     | Working array for NO3-.  Number of   | R | ne        | L |
!          |moles of NO3- remaining after each    |   |           |   |
!          |partition.                            |   |           |   |
! so4      | Array of initial SO4= ion molalities | R | ne        | I |
!          |(moles/kg)                            |   |           |   |
! so4w     | Working array for SO4=.  Number of   | R | ne        | L |
!          |moles of SO4= remaining after each    |   |           |   |
!          |partition.                            |   |           |   |
! sulac    | Estimated molality of sulphuric acid | R | ne        | L |
!          |[H2SO4], which, when dissolved with   |   |           |   |
!          |the other electrolytes, will help     |   |           |   |
!          |create a solution with the same ionic |   |           |   |
!          |concentrations as the aerosol         |   |           |   |
!          |solution.                             |   |           |   |
! sulac1   | Sulphuric acid creating SO4=.        | R | ne        | L |
! sulac2   | Sulphuric acid creating HSO4-.       | R | ne        | L |
!----------------------------------------------------------------------

!  Assign values to working arrays for number of moles of each ion:
!  Note:  on entry, so4, h, no3, nh4, and hso4 are assumed to be in units
!  of moles/kg H2O, with the water content taken from the previous (lwo)
!  value.  Convert this back to moles/(kg air):
   do i = 1, ne
      wso4w   = so4(i)  * lwo(i)
      whw     = h(i)    * lwo(i)
      wno3w   = no3(i)  * lwo(i)
      wnh4w   = nh4(i)  * lwo(i)
      whso4w  = hso4(i) * lwo(i)
!  Letovicite is formed from the available NH4, SO4, H, and half of the HSO4
      wleto = min(wnh4w / 3.0D0, wso4w) !3 moles of NH4 per letovic
      wleto = min(wleto, 0.5D0 * whso4w)   !1.2 of avail. HSO4 all
!  Subtract the number of moles of letovicite (NH4)3H(SO4)2 from the remaining
!  moles of electrolyte components:
      whso4w = whso4w - wleto
      wnh4w = wnh4w - 3.0D0 * wleto
      wso4w = wso4w - wleto
!  Remaining bisulphate enters ammonium bisulphate NH4HSO4:
      wambis = min(whso4w, wnh4w)
      whso4w = whso4w - wambis
      wnh4w = wnh4w - wambis
!  Remaining sulphate used to make ammonium sulphate (NH4)2SO4:
      wamsul = min(wnh4w * 0.5D0, wso4w)  !2 moles of NH4 per ammo
      wnh4w = wnh4w - 2.0D0 * wamsul
      wso4w = wso4w - wamsul
!  Remaining ammonia used to make ammonium nitrate NH4NO3:
      wamnit = min(wnh4w, wno3w)
      wno3w = wno3w - wamnit
!  Remaining sulphate used to make H2SO4:
      wsulac1 = min(wso4w, 0.5D0 * whw)
      whw = whw - 2.0D0 * wsulac1
!  Remaining bisulphate used to make H2SO4:
      wsulac2 = min(whso4w, whw)
      whw = whw - wsulac2
      wsulac = wsulac1 + wsulac2
!  Remaining nitrate used to make HNO3:
      wnitac = min(wno3w, whw)


!  Note that the working (w) ion arrays above may have leftover mass,
!  which could not be combined into electrolytes.  These arrays may
!  be written out for statistical purposes, if desired.
!  The moles/kg of electrolytes can now be used to calculate water content,
!  using the ZSR formula.  Lw = sum (i) of M(i)/m0(i, aw(i)), where
!  M(i) is the number of moles/(kg air) of the i'th electrolyte in the aerosol
!  solution, and m0(i, aw(i)) is the number of moles /(kg water) in a solution
!  having the same water activity as the mixture but composed solely of
!  electrolyte i.
!  The formulae for the terms in the denominator resulted from using
!  Excel to create polynomials for m as a function of aw from formulae
!  for aw as a function of m.  The latter formula came from the following
!  sources:
!  Species      Chemical     Variable     Source
!               Formula        name
!  Ammonium
!  nitrate      NH4NO3         amnit    Chan et al., 1992
!
!  Letovicite   (NH4)3H(SO4)2  leto     Tang and Munkelwitz, 1994
!
!  Ammonium
!  Bisulphate   NH4HSO4        ambis    Tang and Munkelwitz, 1994
!
!  Ammonium
!  Sulphate     (NH4)2SO4      amsul    Tang and Munkelwitz, 1994
!
!  Sulphuric
!  Acid         H2SO4          sulac    Pilinis and Seinfeld, 1987
!
!  Nitric acid  HNO3           nitac    Pilinis and Seinfeld, 1987
!
      lwn(i) = wamnit / (((((-2.057930d+03 * awu(i) + 8.019187d+03) *          &
                             awu(i) - 1.240833d+04) * awu(i) + 9.620334d+03) * &
                             awu(i) - 3.846213d+03) * awu(i) + 6.727864d+02) + &
               wleto / ((((5.84250d+01 * awu(i) - 2.17648d+02) *               &
                           awu(i) + 3.08113d+02) * awu(i) - 2.10129d+02) *     &
                           awu(i) + 6.13002d+01) +                             &
               wambis / (((((-2.55060d-01 * law(i) - 4.05279d+00) *            &
                             law(i) - 1.87386d+01) * law(i) - 1.51834d+01) *   &
                             law(i) - 3.09739d+01) * law(i)) +                 &
               wamsul / ((((1.25661d+02 * awu(i) - 4.81538d+02) *              &
                            awu(i) + 6.77193d+02) * awu(i) - 4.44043d+02) *    &
                            awu(i) + 1.22771d+02) +                            &
               wsulac / ((((((-2.61973d-02 * law(i) - 5.15658d-01) *           &
                              law(i) - 3.79080d+00) * law(i) - 1.26892d+01) *  &
                              law(i) - 1.95522d+01) * law(i) - 1.96898d+01) *  &
                              law(i) + 3.24434d-01) +                          &
               wnitac / (((((-9.82824d-01 * law(i) - 7.22682d+00) *            &
                             law(i) - 2.10233d+01) * law(i) - 2.76882d+01) *   &
                             law(i) - 2.59928d+01) * law(i))
!
      lwn(i) = max(lwn(i), lwmin)
   end do

   return
end subroutine mach_hetv_water
