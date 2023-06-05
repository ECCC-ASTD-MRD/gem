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
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_cam_utils_mod.ftn90
! Creation       : P. Huang, Apr. 2008
! Description    : Declare parameters used for the PM dry deposition scheme
!                  in CAM
!                  Declare variables related to the indices of gas and PM
! Modification   : V. Savic-Jovcic, March 2016: Cleaned up unused variables
!                  A. Akingunola and V. Savic-Jovcic, April 2016:
!                   - Renamed mach_cam_pre subroutine to mach_cam_const, and
!                     merged it with mach_cam_pre_mod to this module.
!                   - mach_cam_const defines gaseous and particle species that
!                     are treated by CAM, as well as coagulation parameters
!                  A. Akingunola, June 2017:
!                   - Centralize the evaluation of the scaling factors (dvn),
!                     and other sub-bin parameters that are used to redistribute the
!                     bin mass to sub-bins in the cam_intrsec and cam_drydepo routines.
!
!============================================================================

module mach_cam_utils_mod
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, pm_nk
   use chm_utils_mod,        only: post_increment, chm_stop
   use chm_species_info_mod, only: species_master
   use chm_species_idx_mod,  only: sp_SO2, sp_H2O2, sp_ROOH, sp_O3, sp_HNO3, &
                                   sp_NH3, sp_SO4, sp_XOOH
   use chm_nml_mod,          only: aerosize, chm_pkg_gas_s ! , nk_start_pm

   implicit none

!   integer(kind=4)  pm_nk     ! number of vertical layers on which aerosol chemistry is applied (pm_nk=chm_nk-nk_start_pm+1)

   integer(kind=4)  icob      ! size bin number with which coagulation will apply

   integer(kind=4)  isize     ! number of size bins
   integer(kind=4)  nbnd      ! number of size sub-divisions
   integer(kind=4)  ntr       ! total number of tracers used in cam (isize * icom + maxns + camextr)
   integer(kind=4)  nsb       ! total number of sub-bins of particles used in cam (icom * nbnd)

   integer(kind=4), parameter :: maxns   = 7   ! No of gas-phase chemistry species
   integer(kind=4), parameter :: maxnsg  = 12  ! No. of gas/part species
   integer(kind=4), parameter :: maxnsaq = 12  ! No. of Aqueous species (excluding fixed H2OA)
   integer(kind=4), parameter :: icom    = 8   ! number of dry aerosol species
   integer(kind=4), parameter :: nswdep  = 15  ! number of species in wet dep flux output

   integer(kind=4), parameter :: maxsize = 12  ! maximum number of size bins

   real(kind=4),    parameter :: tmin    = 0.0 ! minimum thresold value

!  volume and area of unit-radius sphere, commonly used in CAM subroutines
   real(kind=4),    parameter :: ursv = 4.0 * acos(-1.0) / 3. ! Unity-Radius Sphere Volume
   real(kind=4),    parameter :: ursa = 4.0 * acos(-1.0)      ! Unity-Radius Sphere Area

!  dry aerosol names and density [kg/m^3]
!  aerosol type data are in the following order:
!  1. (nh4)2so4
!  2. sea salt
!  3. secondary organic
!  4. nh4no3
!  5. ammonium (modelled as (nh4)2so4)
!  6. soil dust
!  7. black carbon
!  8. primary organic
   character(len=8), dimension(icom) :: aeroname = (/"SULPHATE", "SEA-SALT", "OMCARBON", "NITRATES", &
                                                     "AMMONIUM", "SOILDUST", "BLCARBON", "PMCARBON"/)
! Aliases for 'aeroname', follows the aeroname ordering
   character(len=2), dimension(icom) :: aero_sname = (/"SU", "SS", "OC", "NI", "AM", "CM", "EC", "PC"/)
   real, dimension(icom) :: rhop0 = (/ 1769.0, 2170.0, 1300.0, 1725.0, 1754.3, 2650.0, 1500.0, 1300.0 /)
   real, dimension(icom) :: mwt_aero = (/ 96.06360000, 67.18000000, 132.13420000, 62.00500000, &
                                          18.03850000, 60.08000000, 12.01100000, 132.13420000 /)
!
! SOLUBLE AEROSOL PROPERTIES
!  namsol are the permitted soluble aerosol names
!   character(len=8), dimension(numsol) :: namsol = &
!                (/'SULPHATE', 'SEA-SALT', 'OMCARBON', 'NITRATES', 'AMMONIUM'/)
!  namins are the permitted insoluble aerosol names
!   character(len=8), dimension(numins) :: namins = &
!                (/'SOILDUST', 'BLCARBON', 'PMCARBON'/)
!  soot, dust and primary organic are insoluble
!
!  numsol = total number of soluble aerosol types possible
!  numins = total number of insoluble aerosol types possible
   integer(kind=4), parameter :: numsol = 5, numins = 3
! soluble aerosol types:
!     1 = (nh4)2so4
!     2 = sea salt
!     3 = secondary organic
!     4 = nh4no3
!     5 = ammonium (modelled as ammonium sulfate)   !wanmin
! nu_sol are the ions per solute molecule
   real(kind=4), dimension(numsol), parameter :: nu_sol =  &
                                 (/3.0, 2.165, 1.0, 2.0, 3.0/)
! mw_sol are the molecular weights of the dry aerosol components of soluble
! aerosol types
   real(kind=4), dimension(numsol), parameter :: mw_sol =  &
                                 (/132.1342, 67.180, 400.0, 80.0435, 132.1342/)
! deliq is the deliquescence relative humidity of soluble aerosol types
   real(kind=4), dimension(numsol), parameter :: deliq =   &
                                 (/0.70, 0.74, 0.70, 0.70, 0.70/)
! recry is the recrystallization relative humidity of soluble aerosol types
   real(kind=4), dimension(numsol), parameter :: recry =   &
                                 (/0.50, 0.45, 0.50, 0.50, 0.50/)
! phit(j, i) is the phi(aw) table of coefficients
! i is the solute index:
! j is the data :
!  1 = minimum aw for polynomial fit
!  2 = aw break point between two polynomials
!  3 = x^3 coef for aw > awbreak
!  4 = x^2 coef
!  5 = x^1 coef
!  6 = x^0 coef
!  7 = x^4 coef for aw <= awbreak
!  8 = x^3 coef
!  9 = x^2 coef
!  10= x^1 coef
!  11= x^0 coef
   real(kind=4), dimension(11, numsol), parameter :: phit =                                    &
      reshape ( (/0.39, 0.92, 457.060777, -1280.47495, 1194.81750, -370.739425,     &
                 -1.62440470, 4.07342346, -5.61205075, 3.873682106, -0.216021389,   &
                  0.44, 0.92, 410.74729, -1138.2693, 1049.2792, -320.74562,         &
                 -5.79690208, 17.7685336, -22.5253540, 11.8087027, -0.48210984,     &
                  0.20, 0.85, 83.445, -225.41, 203.45, -60.477,                     &
                  8.9501, -20.468, 16.738, -5.137, 1.0115,                          &
                  0.275, 0.81, 7.6174049, -19.354181, 17.103802, -4.5561686,        &
                 -1.1108526, 3.7035588, -5.1408203, 4.0788267, -0.77326108,         &
                  0.39, 0.92, 457.060777, -1280.47495, 1194.81750, -370.739425,     &
                -1.62440470, 4.07342346, -5.61205075, 3.873682106, -0.216021389/),  &
      (/11, numsol/) )
!
!  condensation/nucleation time intervals
   real(kind=4), dimension(15) :: condnu

!  indices for model gases
   integer(kind=4) igs_SO2
   integer(kind=4) igs_SO4
   integer(kind=4) igs_O3
   integer(kind=4) igs_H2O2
   integer(kind=4) igs_HNO3
   integer(kind=4) igs_ROOH
   integer(kind=4) igs_NH3
! Proxy for the gas-phase mechanism specific organic peroxide specie index
   integer(kind=4) sp_gooH

!  species used by CAM
!   integer(kind=4) igs_H2S, igs_DMS
   integer(kind=4) iae1, iae2

!  indices of aerosol components
   integer(kind=4) iae_SU
   integer(kind=4) iae_SS
   integer(kind=4) iae_OC
   integer(kind=4) iae_NI
   integer(kind=4) iae_AM
   integer(kind=4) iae_SD
   integer(kind=4) iae_EC
   integer(kind=4) iae_PC

!  Wet deposition fields are in the following order
!  1. SO2(g)[igs_SO2];   2. H2O2(g) [igs_H2O2];   3. ROOH(g)[igs_ROOH];
!  4. SO4(=)[iae_SU] ;   5. NO3(-)  [iae_NI];     6. NH4(+) [iae_AM];
!  7. CAT1(s)        ;   8. HCO3(g)         ;     9. H+(ion)[1+4+5+6+8-7];
!  10.H2O            ;  11. SS [iae_SS]     ;    12. OC [iae_OC]
!  13.CM [iae_SD]    ;  14. EC [iae_EC]     ;    15. PC [iae_PC]
!
!  indirect addressing for wet flux (wetflx)
   integer(kind=4), dimension(icom) :: ip_wflx

!  additional parameters used in CAM
   real(kind=4),    dimension(maxns)                       :: mwt_igs
   integer(kind=4), dimension(maxns)                       :: ipos_g

   integer(kind=4), dimension(maxsize)                     :: igf
   integer(kind=4), dimension(maxsize, maxsize*maxsize, 2) :: igfij
   real(kind=4),    dimension(maxsize)                     :: pvol
   real(kind=4),    dimension(2, maxsize)                  :: binrange
   real(kind=4),    dimension(maxsize, maxsize, maxsize)   :: coagfr
!
   real(kind=4),    dimension(maxsize)                     :: sub_rn
   real(kind=4),    dimension(maxsize)                     :: sub_a_n
   real(kind=4),    dimension(maxsize)                     :: sub_pvol
   real(kind=4),    dimension(maxsize)                     :: dvn

!  Diagnostic aerosol fields output parameter
! Total number of diagnostic aggregate PM bins; In 2-bin mode, this is set to
! 2 bins represented in the model, representing (PM2.5 and PM10); however for
! in 12-bin mode it represents (PM1, PM2.5, PM10, and PM_total)
   integer(kind=4) :: pm_agrege_nbin

   save

   contains


!============================================================================
! Name           : sigmacal
!
! Creation       : P. Huang, Mar. 2008 for adapting AURAMS version CAM
!                  P. Huang, Dec. 2007 for GEM-MACH
!
! Description    : 1. Calculate interface values of sigma for CAM model
!                  2. Check tracer mass for CAM model
!
! Extra info     : Varify size bin configration
!
! Arguments:  IN
!               sig            -> Local sigma values
!               ni             -> No. x-direction (W-E) gridpoints
!               nk             -> No. z-direct. vertical levels
!
!             OUT
!               SHTJ          -> Local interface (top) sigma value
!               DSHJ          -> Sigma difference of bottom and top of a layer
!
!============================================================================

   subroutine sigmacal(SHTJ, DSHJ, sig)

      implicit none
      real(kind=4),    intent(in)  :: sig (chm_ni, chm_nk)
      real(kind=4),    intent(out) :: shtj(chm_ni, chm_nk+1), dshj(chm_ni, chm_nk)
      integer(kind=4)              :: i, k

      do i = 1, chm_ni
         shtj(i, 1) = sig(i, 1)
         do  k = 2, chm_nk
            shtj(i, k) = (sig(i, k) + sig(i, k-1)) * 0.5
         end do
!>>      shtj(i, chm_nk + 1) = sig(i, chm_nk)    sg: I assume that this was done because sig(chm_nk) was
!>>                                                  expected to contain the surface (ie sig=1) therefore
!>>                                                  I replace it by the following line
         shtj(i, chm_nk + 1) = 1.0
         do k = 1, chm_nk
            dshj(i, k) = shtj(i, k + 1) - shtj(i, k)
         end do
      end do

      return
   end subroutine sigmacal


!============================================================================
! Name           : sigmathck
!
! Creation       : V. Savic-Jovcic, Mar. 2016 adapted sigmacal into sigmadep
!
! Description    : 1. Calculate depth of model layer as a difference of sigma levels
!
! Arguments:  IN
!               sig       -> local sigma values
!
!             OUT
!               DSIG      -> thickness of the layers that have sig as the centers
!
!============================================================================

   subroutine sigmathck(dsig, sig, pni, pnk)

      implicit none
      integer(kind=4), intent(in)  :: pni, pnk
      real(kind=4),    intent(in)  :: sig (pni, pnk)
      real(kind=4),    intent(out) :: dsig(pni, pnk)
      integer(kind=4)              :: i, k
      real(kind=4)                 :: shtj(pni, pnk+1)

      do i = 1, pni
         shtj(i, 1) = sig(i, 1)
         do  k = 2, pnk
            shtj(i, k) = (sig(i, k) + sig(i, k-1)) * 0.5
         end do
!>>      shtj(i, pm_nk + 1) = sig(i, pm_nk)    sg: P. Huang assumed that this was done because sig(pm_nk) was
!>>                                                  expected to contain the surface (ie sig=1) therefore
!>>                                                  it was replaced by the following line
         shtj(i, pnk + 1) = 1.0
         do k = 1, pnk
            dsig(i, k) = shtj(i, k + 1) - shtj(i, k)
         end do
      end do

      return
   end subroutine sigmathck


!============================================================================
! Name           : mach_cam_const
!
! Description    : 1. Set indices of tracer for gas and PM applied in CAM model
!                  2. Pre-calculations for CAM model
!
! Extra info: CALCULATED INDICES AND PROPERTIES:
!        BINRANGE       -> Bin sizes of dry aerosol (m)
!        COAGFR         -> Coagulation transformation ratio
!        PVOL           -> Bin volume of dry aerosol
!        IGFIJ          -> Indices of coagulation transformation
!        IGF            -> Total indices of coagulation
!
! Arguments:  None
!
! Extra info: Naming convention for bins:
!              - The 2 first characters are used to identify the species
!              - The other 2 are be used for the bin number
!                 - Use numbers for variables in the dynamic bus and
!                 - Use letters for the emissions _C for coarse, _F for fine
!
!============================================================================
   subroutine mach_cam_const()
      use chm_utils_mod,     only: global_debug, chm_lun_out, chm_timestep
      use chm_nml_mod,       only: chm_intrsec_ver, chm_intrsec_ndiv

      implicit none
!
!  Declaration of local variables
!
      integer(kind=4) :: camextr  ! extra species after paticle species

      real(kind=4)    :: pop, pop1

      integer(kind=4) :: i, j, k, n, itot, idx, jsub
      real(kind=4)    :: voij, vokp1, vok, vokm1
      real(kind=4)    :: ri(maxsize)
      real(kind=4)    :: rn(2, nbnd)
      real(kind=4)    :: rndiv, dr_orig, dr_new
      real(kind=4)    :: dp1, dm1, dnp1_dlogd, dnm1_dlogd, &
                         dvp1_dlogd, dvm1_dlogd, vsum
! Parameters for trimodal distribution:
      real(kind=4),               parameter :: sqrt2pi  = sqrt(2.0 * acos(-1.0))
!  Parameters for assumed particle size distribution (Seinfeld and Pandis, 2006,
!  after Jaenicke, 1993, urban aerosol distribution):
      real(kind=4), dimension(3), parameter :: nnn  = (/ 9.93E4, 1.11E3, 3.64E4 /)
      real(kind=4), dimension(3), parameter :: dpm  = (/ 0.013, 0.014, 0.050 /)
      real(kind=4), dimension(3), parameter :: lsig = (/ 0.245, 0.666, 0.337 /)
      logical(kind=4) :: local_dbg

!  Check that the runtime isize is within the assumes maximum bin size
      if (isize > maxsize .or. nbnd > maxsize .or. icob > maxsize) then
         write(*, *)  'CAM module parameters binrange, pvol, coagfr .etc.', &
                    & ' are computed for aerosol bin size ', maxsize
         write(*, *)  'The parameters can be adjusted for isize > 12 in ', &
                    & 'mach_cam_utils_mod.ftn90'
         call chm_stop('Problem with maxsize in mach_cam_utils_mod', -1)
         return
      end if

      local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

! Determine total number of subbins of aerosols treated in CAM
      nsb = icom * nbnd

! Index of gas species in CAM master array
      idx = 1
      igs_SO4  = post_increment(idx)
      igs_SO2  = post_increment(idx)
      igs_H2O2 = post_increment(idx)
      igs_ROOH = post_increment(idx)
      igs_HNO3 = post_increment(idx)
      igs_NH3  = post_increment(idx)
      igs_O3   = post_increment(idx)

! Generic pointer to the organic peroxide gas specie
      if (chm_pkg_gas_s(1:5) == 'ADOM2') then
         sp_gooH = sp_ROOH
      else if (chm_pkg_gas_s(1:7) == 'SAPRC07') then
         sp_gooH = sp_XOOH
      end if
      mwt_igs(igs_SO4)  = species_master(sp_SO4) % mol_wt
      mwt_igs(igs_SO2)  = species_master(sp_SO2) % mol_wt
      mwt_igs(igs_H2O2) = species_master(sp_H2O2) % mol_wt
      mwt_igs(igs_ROOH) = species_master(sp_gOOH) % mol_wt
      mwt_igs(igs_HNO3) = species_master(sp_HNO3) % mol_wt
      mwt_igs(igs_NH3)  = species_master(sp_NH3) % mol_wt
      mwt_igs(igs_O3)   = species_master(sp_O3) % mol_wt

     ! Position of the gas species in the aqueous phase solver
      ipos_g(igs_SO4)  = 0
      ipos_g(igs_SO2)  = 1
      ipos_g(igs_H2O2) = 2
      ipos_g(igs_ROOH) = 3
      ipos_g(igs_HNO3) = 7
      ipos_g(igs_NH3)  = 8
      ipos_g(igs_O3)   = 11

!  Sequence of aerosol species
!     aerosol species ID: 1- SULPHATE,    2- SEA-SALT
!                         3- OMCARBON,    4- NITRATES
!                         5- AMMONIUM,    6- CRUSTAL MATERIAL/SOILDUST
!                         7- BLCARBON,    8- PMCARBON
      do idx = 1, icom
         if (aeroname(idx) == 'SULPHATE') iae_SU = idx
         if (aeroname(idx) == 'SEA-SALT') iae_SS = idx
         if (aeroname(idx) == 'OMCARBON') iae_OC = idx
         if (aeroname(idx) == 'NITRATES') iae_NI = idx
         if (aeroname(idx) == 'AMMONIUM') iae_AM = idx
         if (aeroname(idx) == 'SOILDUST') iae_SD = idx
         if (aeroname(idx) == 'BLCARBON') iae_EC = idx
         if (aeroname(idx) == 'PMCARBON') iae_PC = idx
      end do

!  indirect addressing for wet flux (wetflx)
      ip_wflx(iae_SU)  = 4
      ip_wflx(iae_SS)  = 11
      ip_wflx(iae_OC)  = 12
      ip_wflx(iae_NI)  = 5
      ip_wflx(iae_AM)  = 6
      ip_wflx(iae_SD)  = 13
      ip_wflx(iae_EC)  = 14
      ip_wflx(iae_PC)  = 15

      iae1 = maxns + 1
      iae2 = (icom * isize) + iae1 - 1

      idx = iae2 + 1
! Add 'camextr' fields here
!      igs_H2S = post_increment(idx)
!      igs_DMS = post_increment(idx)
      camextr = idx - iae2 - 1

! Determine total number of tracers treated in CAM
      ntr = isize * icom + maxns + camextr

!  bin size of dry aerosol
      binrange = -999.0
      do n = 1, isize
         binrange(1, n) = aerosize(n)   * 1.e-6
         binrange(2, n) = aerosize(n+1) * 1.e-6
      end do

!  bin volume of dry aerosol
      pvol = -999.0
      do i = 1, isize
         ri(i)   = (binrange(1, i) + binrange(2, i)) * 0.5
         pvol(i) = ursv * ri(i) * ri(i) * ri(i)
      end do

!  coagulation transformation ratio
      coagfr = 0.0
      igf    = 0
      igfij  = 0
      do i = 1, icob
         do j = i, icob
!  new volume of coagulated particle [i, j]
            voij = pvol(i) + pvol(j)
            do k = j, icob
               vok = pvol(k)
               if (k == icob .and. voij >= vok) then
                  coagfr(i, j, k) = 1.0
               end if
               if (k < icob) then
                  vokp1 = pvol(k + 1)
                  if( voij >= vok .and. voij < vokp1) then
                     coagfr(i, j, k) = vok / voij * (vokp1 - voij) / (vokp1 - vok)
                  end if
               end if
               if (k > 1) then
                  vokm1 = pvol(k - 1)
                  if( voij > vokm1 .and. voij < vok) then
                     coagfr(i, j, k) = 1.0 - coagfr(i, j, k - 1)
                  end if
               end if
               coagfr(j, i, k) = coagfr(i, j, k)
            end do
         end do
      end do

      do k = 1, icob
         itot = 0
         do j = 1, icob
            do i = 1, icob
               if (coagfr(i, j, k) > 0.0) then
                  itot = itot + 1
                  igfij(k, itot, 1) = i
                  igfij(k, itot, 2) = j
               end if
            end do
         end do

         igf(k) = itot
         i = igfij(k, itot, 1)
         j = igfij(k, itot, 2)

         if (local_dbg) then
            write (0, *)'^^^^^^ ', k, igf(k), i, j
         end if
      end do
!
!  compute the condensation/nucleation time intervals
      do i = 1, 15
         pop1 = chm_timestep / exp(real(15 - i) / 1.7)
         if (i >= 2) then
            condnu(i) = pop1 - pop
            pop = pop1
         else
            condnu(i) = pop1
            pop = pop1
         end if
      end do
!
! Evaluate the scaling factors (dvn) that can be used to redistribute the
! bin mass to sub-bins in the cam_intrsec and cam_drydepo routines.
! sub_rn, sub_a_n, and sub_pvol contains the subdivided radius, area,
! and volume arrays respectively
!
      if (isize == 12) then
!  skip sub-division when the 12-bin structure is used
         do k = 1, nbnd
            sub_rn(k)   = ri(k)
            sub_a_n(k)  = ursa * sub_rn(k)**2
            sub_pvol(k) = pvol(k)
            dvn(k) = 1.0
         end do
!
      else
!
         sub_rn   = -999.0
         sub_a_n  = -999.0
         sub_pvol = -999.0
         dvn = -999.0
!  Commence subdivision of original bins:
         rndiv = real(chm_intrsec_ndiv)
!
!  Determine bin boundaries of the subdivided distribution:
!
         if (chm_intrsec_ver == 1 .or. chm_intrsec_ver == 4) then
            do k = 1, isize
               do j = 1, chm_intrsec_ndiv
                  jsub = (k - 1) * chm_intrsec_ndiv + j
                  rn(1, jsub) = binrange(1, k) + float(j - 1) *              &
                                (binrange(2, k) - binrange(1, k)) / rndiv
               end do
            end do
!
         else if (chm_intrsec_ver == 2 .or. chm_intrsec_ver == 3 .or. chm_intrsec_ver == 5 .or. chm_intrsec_ver == 6) then
            do k = 1, isize
               do j = 1, chm_intrsec_ndiv
                  jsub = (k - 1) * chm_intrsec_ndiv + j
                  rn(1, jsub) = exp(alog(binrange(1, k)) + float(j - 1) *    &
                                (alog(binrange(2, k)) - alog(binrange(1, k))) / rndiv)
               end do
            end do
         end if
!
         do k = 2, nbnd
            rn(2, k - 1) = rn(1, k)
         end do
         rn(2, nbnd) = binrange(2, isize)
         do k = 1, nbnd
            sub_rn(k)   = (rn(1, k) + rn(2, k)) * 0.5
            sub_a_n(k)  = ursa * sub_rn(k)**2
            sub_pvol(k) = sub_a_n(k) * sub_rn(k) / 3.0
         end do
!
         if (chm_intrsec_ver == 1 .or. chm_intrsec_ver == 4) then
!  version 1 of possible distribution of particle concentration and number density:
!  Use delta(quantity) / delta(radius) and the leftmost boundary value of
!  the quantity being zero to divide up quantities:
            do k = 1, isize
               dr_orig = (binrange(2, k) - binrange(1, k))
               do j = 1, chm_intrsec_ndiv
                  jsub = (k - 1) * chm_intrsec_ndiv + j
                  dr_new = rn(2, jsub) - rn(1, jsub)
                  dvn(jsub) = dr_new / dr_orig
               end do
            end do
!
         else if (chm_intrsec_ver == 2 .or. chm_intrsec_ver == 5) then
!  version 2 of possible distribution of particle concentration and number density:
!  Use delta(quantity) / delta(ln radius) and the leftmost boundary value of
!  the quantity being zero to divide up quantities:
!
            do k = 1, isize
               dr_orig = alog(binrange(2, k) / binrange(1, k))
               do j = 1, chm_intrsec_ndiv
                  jsub = (k - 1) * chm_intrsec_ndiv + j
                  dr_new = alog(rn(2, jsub) / rn(1, jsub))
                  dvn(jsub) = dr_new / dr_orig
               end do
            end do
!
         else if (chm_intrsec_ver == 3 .or. chm_intrsec_ver == 6) then
!  particle concentration is redistributed according to urban volume distribution
!  of Jaenicke, 1993, quoted in Seinfeld and Pandis 2006 pages 370-371.
!  Use delta(quantity) / delta(ln radius) and the leftmost boundary value of
!  the quantity being zero to divide up quantities:
!
!  First, define particle number and volume distribution according
!  to Jaenicke, 1993:
!
            do i = 1, nbnd
               dp1 = rn(2, i) * 2.0 * 1.E6  ! diameter in um
               dm1 = rn(1, i) * 2.0 * 1.E6  ! diameter in um
               dnp1_dlogd = 0.
               dnm1_dlogd = 0.
! construct trimodal dN / dlogD value at boundary radii of sub-bin:
               do j = 1, 3
                  dnp1_dlogd = dnp1_dlogd + nnn(j) / (sqrt2pi * lsig(j))      &
                               * exp( - 1.0 * (alog10(dp1)-alog10(dpm(j)))**2 &
                               / (2. * lsig(j)**2) )
                  dnm1_dlogd = dnm1_dlogd + nnn(j) / (sqrt2pi * lsig(j))      &
                               * exp( - 1.0 * (alog10(dm1)-alog10(dpm(j)))**2 &
                               / (2. * lsig(j)**2) )
               end do
! construct corresponding dV / dlogD value at boundary radii of sub-bin:
               dvp1_dlogd = dnp1_dlogd * ursv * (dp1*0.5)**3
               dvm1_dlogd = dnm1_dlogd * ursv * (dm1*0.5)**3
! integrate between these values :  volume = area under curve. Rectangle rule used:
               dvn(i) = (dvp1_dlogd + dvm1_dlogd) * (alog10(dp1)-alog10(dm1))
            end do
!!
!  The dvn array contains a "typical urban" trimodal volume (i.e. mass) distribution.
!  The sub-bin boundaries rn(2, i) to rn(1, i) have been constructed to fall within the
!  original bin boundaries.  The next step determines scaling factors to use
!  within the bin, based on the relative amount of sub-bin volume within the bin:
            do k = 1, isize
               vsum = 0.
               do j = 1, chm_intrsec_ndiv
                  jsub = (k - 1) * chm_intrsec_ndiv + j
                  vsum = vsum + dvn(jsub)
               end do
               do j = 1, chm_intrsec_ndiv
                  jsub = (k - 1) * chm_intrsec_ndiv + j
                  dvn(jsub) = dvn(jsub) / vsum
               end do
            end do
!
! At this point, the array dvn contains scaling factors that can be used to
! redistribute the bin mass to sub-bins. Within the bin, the trimodal
! distribution will be used, but discontinuities may occur at the
! sub-bin borders.
!
         end if
!
      end if
!
! For PM diagnostics
      if (isize == 2) then
         pm_agrege_nbin = 2
      else
         pm_agrege_nbin = 4
      end if
!
      ! To print details of PM name and size bin structure
      if (local_dbg) then
         write (chm_lun_out, *) 'cam -> isize, icob, icom'
         write (chm_lun_out, 333)  isize, icob, icom
         write (chm_lun_out, 1111) (aeroname(n), n = 1, icom)
         write (chm_lun_out, 2222) (rhop0(n), n = 1, icom)
         write (chm_lun_out, *) '       cam size bin configurations '
         do n = 1, isize
            write (chm_lun_out, *) binrange(1, n), '   ', binrange(2, n)
         end do
      end if
 1111 format (' aerosol type(s) -> ', 8('| ', a8, ' |', 1x))
 2222 format ('     density(kg/m3) ', 9('  ', f8.2, '  ', 1x))
 333  format (8x, i4, 1x, i4, 1x, i3)

      return
   end subroutine mach_cam_const


!============================================================================
! Name           : chk_trc
!
! Creation       : P. Huang, Mar. 2008 for CAM in GEM-MACH
!
! Description    : Check tracer value to avoid any unreasonable tracer value crashing run of CAM
!
! Arguments:  IN
!               trac            -> Tracer array to be checked
!               ntr             -> Number of total tracers for trac
!               strt_trci       -> Start index of the tracers to be checked
!               end_trci        -> End index of the tracers to be checked
!               JLAT            -> J Slice number
!               where_str       -> Identifier shows location of this ckecking process applied
!               f_step          -> Time step in current
!               ipass           -> Identifier to print location of each call when ipass = 9
!
!             OUT
!
!============================================================================

   subroutine chk_trc(trac, strt_trci, end_trci, JLAT, where_str, f_step, &
                      pni, pnk)
      use chm_utils_mod, only: chm_error_l

      implicit none
      integer(kind=4),  intent(in)  :: pni, pnk
      integer(kind=4),  intent(in)  :: strt_trci, end_trci, jlat, f_step
      real(kind=4),     intent(in)  :: trac(pni, pnk, ntr)
      character(len=8), intent(in)  :: where_str
      integer(kind=4)               :: n, i, l

      do n = strt_trci, end_trci
         do l = 1, pnk
            do i = 1, pni
               if (trac(i, l, n) > 0.1 .or. trac(i, l, n) < 0.0) then
                  write (0, *)' tracer concentration is too high or negative'
                  write (0, *)' overflow is likely. program stops at:', where_str
                  write (0, *)' tracer ', n, '= ', trac(i, l, n)
                  write (0, *)' level    grid pt     slice       step# '
                  write (0, *)   l,    i,     jlat,        f_step
                  chm_error_l = .true.
                  return
               end if
            end do
         end do
      end do

      return
   end subroutine chk_trc

end module mach_cam_utils_mod
