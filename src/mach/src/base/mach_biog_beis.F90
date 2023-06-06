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
! Fichier/File   : mach_biog_beis.ftn90
! Creation       : H. Landry, S. Menard, A. Kallaur - mai 2006
! Description    : Compute the biogenics emissions
!
! Extra Info     : From bioe3.ftn90, AURAMS v1.3.1
!                  Geron et. al.
!                  An Improved Model for Estimating Emissions of Volatile Organic
!                  Compounds from Forests in the Eastern United States, JGR Vol. 99
!                  No. D6, pp 12, 773-12, 791, June 20 1994
!
!                  Vose J. M. and W. T. Swank
!                  Assessing Seasonal Leaf Area Dynamics and Vertical Leaf Area
!                  Distributions in Eastern White Pine (Pinus Strobus L.) with a
!                  Portable Light Meter, Tree Physiol., 7, pp 125-134 1990.
!
!                  Vose et al (1995)
!                  Vertical leaf area distribution, light transmittance, and
!                  application of the Beer-Lambert Law in four mature hardwood stands
!                  in the Southern Appalachians, Canadian Journal of Forest Research
!
!                  E.J. Williams, A. Guenther and F.C. Fehsenfeld
!                  An Inventory of Nitric Oxide Emissions from Soils in the United
!                  States, JGR Vol. 97 pp 7511-7519, May 20 1992
!
!                : Renamed to mach_biog_beis.ftn90 and rewritten to enable emissions
!                  into canopy shaded layer (Deji Akingunola and Paul Makar, Nov. 2020)
!
! Arguments:  IN
!               summer_std_p     -->  Summer biogenic standard basal rate
!               winter_std_p     -->  Winter biogenic standard basal rate
!               air_temp         -->  Air temperature
!               light_correction -->  cosine of solar angle
!               canopy_column    -->  Canopy presence logical
!               canopy_hgt       -->  Height of canopy within model grid
!               laifrac          -->  Fraction of LAI at heights of 1.0, 0.5 and 0.2 canopy height
!               sesn_coef        -->  seasonal coefficients
!
!             OUT
!               emissbio         -->  gas-phase species' biogenic emissions
!
!==============================================================================
!
!!if_on
subroutine mach_biog_beis(emissbio, emissbio_can, summer_std_p, winter_std_p, &
                          sesn_coef, air_temp, zmomcan, light_correction,     &
                          laifrac, canopy_hgt, canopy_column, ni_can)
   use chm_ptopo_grid_mod,   only: chm_ni, nkt, nkc
   use mach_pkg_gas_mod,     only: num_be_std, num_be_sp
!!if_off
   use chm_nml_mod,          only: chm_pkg_gas_s
   use chm_utils_mod,        only: chm_error_l, global_debug
   use chm_consphychm_mod,   only: rgasi
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_ALKA, sp_ALKE, sp_ISOP, sp_OLE2, &
                                   sp_ALK3 !, sp_NO
   use mach_pkg_gas_mod,     only: be_std_name, be_species

   implicit none
!!if_on
   integer(kind=4), intent (in) :: ni_can
   real(kind=4),    intent(out) :: emissbio        (num_be_sp, chm_ni)
   real(kind=4),    intent(out) :: emissbio_can    (ni_can, nkt, num_be_sp)
   real(kind=4),    intent (in) :: winter_std_p    (num_be_std, chm_ni)
   real(kind=4),    intent (in) :: summer_std_p    (num_be_std, chm_ni)
   real(kind=4),    intent (in) :: light_correction(nkc + 1, chm_ni)
   real(kind=4),    intent (in) :: air_temp        (nkc + 1, chm_ni)
   real(kind=4),    intent (in) :: laifrac         (nkc + 1, chm_ni)
   real(kind=4),    intent (in) :: canopy_hgt      (ni_can)
   real(kind=4),    intent (in) :: zmomcan         (ni_can, nkt)
   real(kind=4),    intent (in) :: sesn_coef       (chm_ni)
   logical(kind=4), intent (in) :: canopy_column   (chm_ni)
!!if_off
!
!  Local (internal) variables
!
   real(kind=4),    dimension(num_be_std) :: basal_rate

   integer(kind=4), dimension(nkt, ni_can)    :: nfrcbio
   integer(kind=4), dimension(2, nkt, ni_can) :: ifrcbio
   real(kind=4),    dimension(2, nkt, ni_can) :: frcbio
   real(kind=4),    dimension(num_be_sp, nkc + 1) :: added_voc
!
!     * Formula variables
   real(kind=4) :: tt, tfac
   real(kind=4) :: no_corr_fact, iso_corr_fact, voc_mono_corr_fact
   real(kind=4) :: ovoc_alka_ratio, ovoc_alke_ratio, mono_alke_ratio, vocsum
   real(kind=4), dimension(num_be_std, nkc + 1) :: be_std_emiss

   real(kind=4), parameter :: small = tiny(1.0)

   integer(kind=4) :: i_isop, i_no, i_alka, i_alke, i_mono, i_ovoc
!     * Loop counters
   integer(kind=4) :: i, ic, kc, k, kk, nsp, nlev

! Estimate biogenic emissions level fractions within canopy
   if (any(canopy_column)) then
      call mach_canopy_bio_levels(zmomcan, canopy_hgt, frcbio, nfrcbio, ifrcbio, ni_can, nkt)
      if (chm_error_l) return
   end if

   i_no = 1   ! findloc(be_species, sp_NO, dim = 1)
   i_isop = findloc(be_species, sp_ISOP, dim = 1)
   i_mono = findloc(be_std_name, "MONO", dim = 1)
   i_ovoc = findloc(be_std_name, "OVOC", dim = 1)
   if (chm_pkg_gas_s(1:5) == 'ADOM2') then
      i_alka = findloc(be_species, sp_ALKA, dim = 1)
      i_alke = findloc(be_species, sp_ALKE, dim = 1)
      mono_alke_ratio = sm(be_species(i_alke)) % mol_wt / 136.0
   else if (chm_pkg_gas_s(1:7) == 'SAPRC07') then
      i_alka = findloc(be_species, sp_ALK3, dim = 1)
      i_alke = findloc(be_species, sp_OLE2, dim = 1)
   end if
   ovoc_alka_ratio = sm(be_species(i_alka)) % mol_wt / 148.0
   ovoc_alke_ratio = sm(be_species(i_alke)) % mol_wt / 148.0

!  Initializations
   emissbio = 0.0
   emissbio_can = 0.0

! ----- Main loop
   ic = 0
   do i = 1, chm_ni
!
!  Applying the seasonal coefficients to the emissions
!  Set the basal emission rates (which would apply for the whole column and its LAI):
      do nsp = 1, num_be_std
         basal_rate(nsp) = (summer_std_p(nsp, i) - winter_std_p(nsp, i)) * &
                           sesn_coef(i) + winter_std_p(nsp, i)
      end do
!
      if (canopy_column(i)) then
         ic = ic + 1
         nlev = nkc + 1
      else
         nlev = 1
      end if
!
!* -----------------------------------------------------------------------------
!*    Temperature correction for NO emissions
!*       For reference, beis 2 had 4 formulas for forest, grassland, wetland an agricultural use
!*       Only the agricultural part has been kept and applies to all. As far as I can tell
!*            nof(i, j, 4) = exp(alfa*(0.72*tt(i, j)+52.28-tb))
!*       (Ref Geron et. al. and Williams et. al.)
!*
!*       NOTE : under -4 C, the soil is frozen and the correction is consistent
!* -----------------------------------------------------------------------------
!
      tt = min(303.0, air_temp(nlev, i))
      if (tt > 268.890) then
         no_corr_fact = exp(0.05112 * tt - 15.68248)  ! agriculture
      else
         no_corr_fact = 0.0
      end if
!
!* -----------------------------------------------------------------------------
!*    Calculate the NO emission rate (surface only)
!* -----------------------------------------------------------------------------
      emissbio(i_no, i) = basal_rate(i_no) * no_corr_fact
!
!  For the other VOC and monoterpernes species, emissions are assumed to be proportional
!  to the amount of foliage within the layer, in turn multiplied by the appropriate scaling factors:
      do kc = 1, nlev
         tt = min(air_temp(kc, i), 315.0)
!
!* -----------------------------------------------------------------------------
!*    Temperature correction for isoprene
!*       The formula used is (ref. Geron et al) :
!*
!*                     exp(ct1*(T-Ts)/R*Ts*t)
!*         ct(i, j) = --------------------------
!*                   1 + exp(ct2*(T-Tm)/R*Ts*t)
!*
!*       NOTE : over 40 C leaf stomata are closed and the correction is constant
!* -----------------------------------------------------------------------------
         tfac = 1.0 / (rgasi * tt * 303.0)
         iso_corr_fact = exp(9.5E4 * (tt - 303.0) * tfac) / &
                             (1.0 + exp(2.3E5 * (tt - 314.0) * tfac))
!
         be_std_emiss(i_isop, kc) = basal_rate(i_isop) * iso_corr_fact *  &
                                    laifrac(kc, i) * light_correction(kc, i)   ! g/s
         added_voc(i_isop, kc) = be_std_emiss(i_isop, kc)  ! g/s added to layer
!
!* -----------------------------------------------------------------------------
!*    Temperature correction for mono and other VOC
!*       Calculate Temperature correction factors for Monoterpenes(M) and
!*       OVOC (Other VOC). The correction function is ref. (Geron et. al) :
!*
!*         OV(or M) = (OVs or Ms)*exp[B(T-Ts)]
!*
!*         NOTE1: The same formula is used for OVOC as Monoterpenes.
!*         NOTE2: See
!* -----------------------------------------------------------------------------
         voc_mono_corr_fact = exp(0.09 * (tt - 303.0))
!
! All the other standard biogenic emission 'species' are assumed to also utilize the OVOC
! temperature correction factor
         do nsp = i_isop + 1, num_be_std
            be_std_emiss(nsp, kc) = basal_rate(nsp) * laifrac(kc, i) * voc_mono_corr_fact ! g/s
         end do
      end do
!
      if (chm_pkg_gas_s(1:5) == 'ADOM2') then
!     Convert mono and voc to alke and alka
!
!     ALKE:  sum of 50% of the VOC mass and all of the monoterpene,
!     both weighted by the ratio of ALKE mass to (monoterpene or OVOC) mass.
!     The weighting is to convert moles of the emitted species, monoterpene
!     or OVOC to moles of the model species, ALKE.
!     OVOC molecular mass of 147.26056 is from assuming 50% monoterpene (c10H16)
!     and 50% C10 alcohol (CH3(CH2)8CH2OH).
!     To be consistent with values in BEIS, mono mass is 136.0 and ovoc is 148.0
!
!     ALKA:  50% of the OVOC, molecular mass weighted to get moles of ALKA
!
         do kc = 1, nlev
            added_voc(i_alka, kc) = 0.5 * ovoc_alka_ratio * be_std_emiss(i_ovoc, kc)
            added_voc(i_alke, kc) = 0.5 * ovoc_alke_ratio * be_std_emiss(i_ovoc, kc) + &
                                    mono_alke_ratio * be_std_emiss(i_mono, kc)
         end do
!
      else if (chm_pkg_gas_s(1:7) == 'SAPRC07') then
!    As for ADOM, split OVOC equally between the alkanes and alkenes
         do kc = 1, nlev
            added_voc(i_mono, kc) = be_std_emiss(i_mono, kc)
            added_voc(i_alka, kc) = 0.5 * ovoc_alka_ratio * be_std_emiss(i_ovoc, kc)
            added_voc(i_alke, kc) = 0.5 * ovoc_alke_ratio * be_std_emiss(i_ovoc, kc)
         end do
      end if
!
      if (canopy_column(i)) then

!  At this point, added_vocs are the amounts of each chemical, in g/s,
!  added to each of the the following canopy layers:
!  where canopy height = hc:
!  kc = 1:  hc to 0.75 hc
!  kc = 2:  0.75 hc to 0.35 hc
!  kc = 3:  0.35 hc to the surface
!  kc = nlev: surface
!
! Note that the layer additions noted above make use of light attenuation
! factors which are not for the entire layer (while the emission factors have
! been designed for an integrated total). Therefore, the mass of biogenic
! emissions will be conserved by ratioing to the totals from the full layer
! calculation:
!
         do nsp = i_no + 1, num_be_sp
            vocsum = max(small, sum(added_voc(nsp, 1:nkc)))  !  totals in g/s:
            do kc = 1, nkc
               added_voc(nsp, kc) = added_voc(nsp, kc) / vocsum * added_voc(nsp, nlev)
            end do
!
! At this point, the total mass of biogenic VOCs summed across canopy
! layers will be the same as added as a surface flux term in the original
! approach.  It will be distributed according to the response of the
! original layers to the variations in light, LAI, etc.
!
!  Add to appropriate resolved model or canopy model layers
!  Need to
!  (1) Determine the g/s addition to each layer
            do k = 1, nkt
               if (nfrcbio(k, ic) > 0) then
                  do kk = 1, nfrcbio(k, ic)
                     kc = ifrcbio(kk, k, ic)
                     emissbio_can(ic, k, nsp) = emissbio_can(ic, k, nsp) + &
                                              added_voc(nsp, kc) * frcbio(kk, k, ic)
                  end do
               end if
            end do
         end do
         emissbio_can(ic, nkt, i_no) = emissbio(i_no, i)
!
! For non-canopy columns;
      else
         do nsp = i_no + 1, num_be_sp
            emissbio(nsp, i) = added_voc(nsp, nlev)
         end do
      end if
!
   end do ! end of main loop over chm_ni

   return
!
   contains
!
   subroutine mach_canopy_bio_levels(zmomcan, hc, frcbio, nfrcbio, ifrcbio, ni_can, nkt)
      integer(kind=4), intent(in)  :: nkt, ni_can
      real(kind=4),    intent(in)  :: hc(ni_can)
      real(kind=4),    intent(in)  :: zmomcan(ni_can, nkt)
      integer(kind=4), intent(out) :: nfrcbio(nkt, ni_can)
      integer(kind=4), intent(out) :: ifrcbio(2, nkt, ni_can)
      real(kind=4),    intent(out) :: frcbio(2, nkt, ni_can)
!
      real(kind=4), parameter :: dzmin = 0.0005
      integer(kind=4)         :: ic, k, k2, kk, kcan_top
      real(kind=4)            :: zmombio_diff
      real(kind=4), dimension(4)     :: zmombio
      real(kind=4), dimension(nkt+1) :: zmom_1d
      real(kind=4), dimension(3)     :: sumc2
      logical(kind=4)                :: local_dbg
!
!  Within the revised biogenic emissions module, the emissions are worked out
!  for three layers (where hc is the height of the canopy):
!   (1) hc to 0.75 hc, thickness 0.25 hc
!   (2) 0.75 hc to 0.35 hc, thickness 0.40 hc
!   (3) 0.35 hc to the surface, thickness 0.35 hc
!  An equivalent mapping is required to take ->these<- three layers to
!  the extended set of nkt levels in the canopy column.  Each of the three layers
!  will be greater than or equal in size to the corresponding combined layers
!  in the column, with the exception of the top-most layer, so the same algorithm
!  as above is followed, with an additional case.
!
!  nfrcbio(k,ic)   : the number of biogenic emissions model levels
!                    contributing to canopy level k
!  ifrcbio(n,k,ic) : the index of the biogenic emissions level
!                    contributing to canopy level k (n is at most 2)
!  frcbio(n,k,ic)  : the fractional contribution of the biogenic
!                    emissions level ifrcbio which goes to canopy level k
!
      nfrcbio = 0
      ifrcbio = 0
      frcbio = 0.0
!
!  Check for coincident layers first:
!
      do ic = 1, ni_can
!
! Array zmombio generates the upper and lower boundaries of the three regions within the canopy
! which are used to generate by-layer emissions of biogenic VOC species (biogenic NO is assumed
! to be emitted from the decaying leaf material on the ground).
         zmombio(1) = hc(ic)
         zmombio(2) = 0.75 * hc(ic)
         zmombio(3) = 0.35 * hc(ic)
         zmombio(4) = 0.0
!
         zmom_1d(nkt + 1) = 0.0
         zmom_1d(1:nkt)   = zmomcan(ic, 1:nkt)
! The model level above the tallest canopy in grid
         kcan_top = 2
         do k = nkt, 3, -1
            if (zmom_1d(k) > hc(ic)) then
               kcan_top = k - 1
               exit
            end if
         end do
         do k = kcan_top, nkt
            do kk = 1, nkc
!(0)  Upper and lower boundaries coincide; trivial case
! Note the use of differences and a limit of dzmin to detect "same layer", required
! due to round-off issues in the numbers
               if ((abs(zmom_1d(k) - zmombio(kk)) < dzmin .and. abs(zmom_1d(k+1) - zmombio(kk+1)) < dzmin) .or. &
! (1a) Upper boundaries of combined and biogenic emissions layers coincide,
! lower boundary of biogenic emissions layer is above combined model layer:
                   (abs(zmom_1d(k) - zmombio(kk)) < dzmin .and. (zmombio(kk+1) - zmom_1d(k+1)) > dzmin) .or. &
!  (2b) Lower boundaries coincide, upper boundary of biogenic layer is within canopy layer:
                   (abs(zmom_1d(k+1) - zmombio(kk+1)) < dzmin .and. (zmom_1d(k) - zmombio(kk)) > dzmin) .or. &
!  (3b) Biogenic emissions layer exists inside a canopy model layer
                   ((zmom_1d(k) - zmombio(kk)) > dzmin .and. (zmombio(kk+1) - zmom_1d(k+1)) > dzmin)) then
                  nfrcbio(k, ic) = 1
                  ifrcbio(1, k, ic) = kk
                  frcbio(1, k, ic) = 1.0
                  cycle
               end if
!  (1) Lower boundaries coincide, upper boundary of combined layer is within biogenic emissions layer:
               zmombio_diff = 1.0  / (zmombio(kk) - zmombio(kk+1))
               if ((abs(zmom_1d(k+1) - zmombio(kk+1)) < dzmin .and. (zmombio(kk) - zmom_1d(k)) > dzmin) .or. &
!  (2) Both combined model layer boundaries exist inside a biogenic emissions model layer:
                   ((zmombio(kk) - zmom_1d(k)) >  dzmin .and. (zmom_1d(k+1) - zmombio(kk+1)) > dzmin)) then
                  nfrcbio(k, ic) = 1
                  ifrcbio(1, k, ic) = kk
                  frcbio(1, k, ic) = (zmom_1d(k) - zmom_1d(k+1)) * zmombio_diff
               end if
! (3) Upper boundaries of combined and biogenic emissions layers coincide,
! lower boundary of combined layer is within biogenic emissions layer:
!   (4) Biogenic emissions layer boundary splits a combined canopy layer:
               if ((zmom_1d(k) >= zmombio(kk)) .and. (zmom_1d(k+1) < zmombio(kk)) .and. &
                   (zmom_1d(k+1) > zmombio(kk+1))) then
! Upper contribution from emissions layer kk, to combined model layer k:
                  nfrcbio(k, ic) = nfrcbio(k, ic) + 1
                  ifrcbio(nfrcbio(k, ic), k, ic) = kk
                  frcbio(nfrcbio(k, ic), k, ic) = (zmombio(kk) - zmom_1d(k+1)) * zmombio_diff
! Lower contribution from emissions layer kk, to combined model layer k+1:
                  nfrcbio(k+1, ic) = nfrcbio(k+1, ic) + 1
                  ifrcbio(nfrcbio(k+1, ic), k+1, ic) = kk
                  frcbio(nfrcbio(k+1, ic), k+1, ic) = (zmom_1d(k+1) - zmombio(kk+1)) * zmombio_diff
               end if
            end do
         end do
!
         local_dbg = (.false. .or. global_debug)
         if (local_dbg) then
!  Check on the values of the fractions:  they should sum to unity across the number
!  of original model levels!
            sumc2 = 0.
            do k = nkt, 1, -1
               if (nfrcbio(k, ic) > 0) then
                  do kk = 1, nfrcbio(k, ic)
                     kc = ifrcbio(kk, k, ic)
                     sumc2(kc) = sumc2(kc) + frcbio(kk, k, ic)
                  end do
               end if
            end do
            do kk = nkc, 1, -1
               if (sumc2(kk) < 0.999 .or. sumc2(kk) > 1.001 ) then
               write(6,*) 'Layer mismatch in biogenic level to canopy level matching - stopping code'
               write(6,*) 'sum of non-zero contributions from 3 canopy layer in column ',ic, &
               ' layer ',kk,' is ',sumc2(kk),' (should be unity)'
               write(6,*) 'values of zmombio: ',(zmombio(k), k = 1,nkc+1)
               write(6,*) 'values of zmomcan: '
               do k = kcan_top, nkt+1
                  write(*,*) k, zmom_1d(k)
               end do
               write(6,*) 'values of nfrcbio: ',(nfrcbio(k,ic), k = 1,nkt)
               write(6,*) 'values of ifrcbio: ',((ifrcbio(k2,k,ic), k2 = 1, nfrcbio(k,ic)), k = 1,nkt)
               write(6,*) 'values of frcbio: ',((frcbio(k2,k,ic), k2 = 1, nfrcbio(k,ic)), k = 1,nkt)
               write(6,*) 'values of sumc2: ',(sumc2(k), k=1,nkc)
               chm_error_l = .true.
               return
               end if
            end do
         end if
!
      end do
22 format(a100,1x,i4,1x,2(1x,i3),1x,4(1x,1pe10.3,1x),1pe11.4)
!
      return
   end subroutine mach_canopy_bio_levels

end subroutine mach_biog_beis
