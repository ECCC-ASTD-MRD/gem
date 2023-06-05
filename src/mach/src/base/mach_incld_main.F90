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
! Fichier/File   : mach_incld_main.ftn90
! Creation       : Sylvain Menard, Wanmin Gong, Ashu Dastoor, Ping Huang,
!                  Craig Stroud, Sylvie Gravel, B. Pabla for GEM-MACH, June 2008
!
! Description    : This is a vectorize version of the ADOM stratus aqueous phase chemistry solver.
!                  Modified for hosting aq chem and tracer-cloud interaction.
!                  It also splits aerosol OC into primary and secondary components.
!                  It compute cloud-to-rain conversion and wet flux
!
! Extra info     : A note on units:
!                  On input it is assumed that gaseous species concentrations and
!                  aerosol species concentrations are in kg/kg-of-air, cloud
!                  water content q_bin are in kg/m^3-of-air. Internally for the aqueous-
!                  phase chemistry integration, gaseous species concentrations are
!                  converted to PPMV, aerosol (or aqueous) species concentration to
!                  moles per litre of water (MOLAR), and cloud water content (liquid
!                  and solid) to g/m^3-of air. Concentration units are all converted
!                  back before exiting this subroutine.
!
! Arguments:
!           IN
!
!            tempk       --> Atmospheric Temperature (Kelvin)
!            psacw       --> CW to snow collection rates(gm/m3 s)
!            rcrit       --> Critical radius expressed in terms of real bin index which is determine from aerosol activation
!            roarow      --> air density (kg/m3)
!            fctr        --> Cloud-to-rain conversion rate (s-1)
!            RAD1        --> CW droplet radii (m)
!
!           OUT
!            flux        --> Cloud-to-rain conversion rate (s-1)
!
!           IN/OUT
!            G           --> gas/part species conc (kg/kg_air)
!            Q_BIN       --> Liquid water conc of cloudwater, Ice/Snow, & rainwater in air (kg/m**3-air)
!            AEROCON     --> Aerosol concentration (aq and non aq species) kg/kg_air
!            AERONUM     --> Number conconcenration
!=============================================================================
!
!!if_on
subroutine mach_incld_main(gaz_conc, aerocon, q_bin, tempk, psacw, rad1, rcrit,&
                           roarow, ibulk, flux, fctr, aeronum, pni, pnk)
   use mach_cam_utils_mod,     only: icom, maxnsg, maxns, isize, nswdep
!!if_off
   use chm_consphychm_mod,     only: rgasd, rho_h2o, mwt_air
   use chm_utils_mod,          only: chm_timestep, chm_error_l, CHM_MSG_DEBUG
   use chm_nml_mod,            only: chm_bkgd_co2, chm_timings_L
   use mach_incld_headers_mod, only: mach_incld_soleq, mach_incld_upaqr, &
                                     mach_incld_dochem, mach_incld_steady
   use mach_incld_mod,         only: aqmin, gmin
   use mach_cam_utils_mod,     only: iae_SU, iae_NI, iae_AM, mwt_aero, & !iae_SD
                                     maxnsaq, maxnsg, mwt_igs, ipos_g, ip_wflx
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   integer(kind=4), intent   (in) :: ibulk
   real(kind=4),    intent(inout) :: gaz_conc(pni, pnk, maxns)
   real(kind=4),    intent(inout) :: aerocon (pni, pnk, icom, isize)
   real(kind=4),    intent(inout) :: q_bin   (pni, pnk, isize, 2)
   real(kind=4),    intent   (in) :: tempk   (pni, pnk)
   real(kind=4),    intent   (in) :: psacw   (pni, pnk)
   real(kind=4),    intent   (in) :: rad1    (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rcrit   (pni, pnk)
   real(kind=4),    intent   (in) :: roarow  (pni, pnk)
   real(kind=4),    intent  (out) :: flux    (pni, pnk, nswdep)
   real(kind=4),    intent   (in) :: fctr    (pni, pnk)
   real(kind=4),    intent(inout) :: aeronum (pni, pnk, isize)
!!if_off
!
! Local variables
!
   integer(kind=4)                 :: nptsnz
   integer(kind=4)                 :: ii, kk, jj, mm, ijk
   integer(kind=4)                 :: ibin, isize_0, nmax, nq, ncw
   integer(kind=4)                 :: idrf, ij, iwf, iae
   integer(kind=4), parameter      :: icom_aq = 7 ! Number  of aerosol components
                                                  ! taking part in aqueous chemistry

   real(kind=4)                    :: convert, ffr, lrt, resid, factor, kg2mol
   real(kind=4), parameter         :: smf = 1.0e-30
   real(kind=4), parameter         :: dtspmx = 225.0
   real(kind=4)                    :: q_bin_t, aq_bulk, aqnew_diff, st_ion, remfrc
   real(kind=4)                    :: qbin_bulk
   real(kind=4)                    :: time, tfinal, tout, tin, dtsplt
   real(kind=4)                    :: rndrop, bb, arg, aa, aaa

!  set mean free path in air (m), sticking coeff, etc.
   real(kind=4), parameter         :: mfpl = 1.0e-07, stik = 1.0e-3
   real(kind=4), parameter         :: bbconst = 0.7 + (1.333333333 * &
                                                       (1.0 - stik) / stik)
!  set molecular weights of gases
   real(kind=4), dimension(maxnsg), parameter :: gmwt = &
                                                (/64.0, 34.0,  62.0, 98.0, &
                                                 115.0, 132.0, 63.0, 17.0, &
                                                 80.0,  100.0, 48.0, 44.0/)
   real(kind=4), dimension(maxnsg) :: sq_gmwt

   real(kind=4), dimension(pni*pnk, maxnsg)             :: gnew, g0
   real(kind=4), dimension(pni*pnk, maxnsg)             :: delg
   real(kind=4), dimension(pni*pnk, maxnsaq, 2)         :: aqnew, aqnew_before
   real(kind=4), dimension(pni*pnk, maxnsaq, 2)         :: daq_chem, daq_flux
   real(kind=4), dimension(pni*pnk, 2)                  :: qnew
   real(kind=4), dimension(pni*pnk)                     :: tempknew, tempinew
   real(kind=4), dimension(pni*pnk)                     :: rtempnew, psacwnew
   real(kind=4), dimension(pni*pnk, maxnsg, 2)          :: pgscav
   real(kind=4), dimension(pni*pnk, 25, 2)              :: raq
   real(kind=4), dimension(pni*pnk, 5, 2)               :: baq
   real(kind=4), dimension(pni, pnk, maxnsg)            :: gtmp
   real(kind=4), dimension(pni, pnk)                    :: rhom, fract
   real(kind=4), dimension(pni, pnk, isize, 2, icom_aq) :: aq, aq_orig

   real(kind=4), dimension(pni, pnk, icom_aq)           :: aqfrc
   real(kind=4), dimension(pni, pnk, isize)             :: vfrac

   real(kind=4), dimension(pni, pnk, 2, icom_aq)        :: daq_chem_f, daq_flux_f

   real(kind=4), dimension(pni*pnk, nswdep)             :: fluxnew
   integer(kind=4), dimension(pni*pnk)                  :: iaq
!
!  External subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_incld [BEGIN]')
   if (chm_timings_L) call timing_start_omp(350, 'mach_incld', 340)

   sq_gmwt = sqrt(gmwt)

   q_bin_t     = 0.0

   nptsnz = pni * pnk

   dtsplt = dtspmx  ! initial dtsplt (225 s)

   gaz_conc = max(gaz_conc, gmin)

   aq = 1.0e-18
   aq(:, :, :, :, 6) = 1.0e-7
   aq(:, :, :, :, 7) = 1.0e-7
!
   do mm = 2, maxns
      ij = ipos_g(mm)

      do kk = 1, pnk
         do ii = 1, pni
!  gas perfect law .
!  here we transform gas from kg/kg to ppmv assuming v1/v2(volume ratio) = n1/n2 (molar ratio)
!  this is an approximation involving (p1=p2)
            gtmp(ii, kk, ij) = 1.0e+06 * gaz_conc(ii, kk, mm) * (mwt_air / mwt_igs(mm))
         end do
      end do

   end do

!  Initialization
!
!  (Note H2SO4 in the
!  gas/particle list is not the same as H2SO4, i.e. sulphuric acid in
!  gas-phase, in AURAMS.

   do kk = 1, pnk
      do ii = 1, pni
         gtmp(ii, kk, 4)  = 0.0 !h2so4
         gtmp(ii, kk, 5)  = 0.0 !nh4hso4
         gtmp(ii, kk, 6)  = 0.0 !nh42so4
         gtmp(ii, kk, 9)  = 0.0 !nh4no3
         gtmp(ii, kk, 10) = 0.0 !dust
         gtmp(ii, kk, 12) = chm_bkgd_co2
      end do
   end do

   if (ibulk == 1) then
!  initialization of vfrac
      vfrac = 0.0

!  computing qbin_bulk (kg/m3) from qbin (kg/m3) and also
!  Computing VFRAC ... After discussions, we decide to calculate VFRAC by using qbulk and q_bin

      do kk = 1, pnk
         do ii = 1, pni
            qbin_bulk = 0.0 ! water
            isize_0 = int(rcrit(ii, kk))
            do ibin = isize_0, isize
               qbin_bulk = qbin_bulk + q_bin(ii, kk, ibin, 1)
            end do

!  we must calculate vfrac here because information's about cwc
!  in each individual bin is lost a few line below when we
!  transfer qbulk into the first size bin of q_bin .
            do ibin = isize_0, isize
               if (qbin_bulk > 0.0) then
                  vfrac(ii, kk, ibin) = q_bin(ii, kk, ibin, 1) / qbin_bulk
               end if
            end do
!  transfer qbin_bulk into the first size bin of q_bin(overwrite the first bin!)
            q_bin(ii, kk, 1, 1) = qbin_bulk ! save qbin_bulk in first bin
            q_bin(ii, kk, 1, 2) = 0.0

!  assign qbulk to qbin(ii) for ii=isize_0, isize
            do ibin = isize_0, isize
               q_bin(ii, kk, ibin, 1) = qbin_bulk
               q_bin(ii, kk, ibin, 2) = 0.0
            end do
         end do
      end do
   end if

!  Assign aerocon to aq and unit conversion (kg/kg -> MOLAR)
!  and partitioning conc for partially activated bin

   do kk = 1, pnk
      do ii = 1, pni
!  Air density (in g/m^3)
         rhom(ii, kk) = roarow(ii, kk) * 1000.0

         isize_0 = int(rcrit(ii, kk))
         fract(ii, kk) = max(0.0, (1.0 - (rcrit(ii, kk) - real(isize_0))))

         do ibin = isize_0, isize
            q_bin_t = q_bin(ii, kk, ibin, 1) + q_bin(ii, kk, ibin, 2)
            if (q_bin_t <= 0) cycle

!  Convert = den_air * den_water / LWC / MWT
!  (kg/m^3a) (kg/m^3w)  (kg/m^3a)
!  Note => kmol/m^3 = mol/L, rhom in g/m^3 and den_water = 1000

            convert = rhom(ii, kk) / q_bin_t

            do jj = 1, 2
               ffr = q_bin(ii, kk, ibin, jj) / q_bin_t
               kg2mol = convert * ffr
               aq(ii, kk, ibin, jj, 1) = aerocon(ii, kk, iae_SU, ibin) / mwt_aero(iae_SU) * kg2mol ! so4=
               aq(ii, kk, ibin, jj, 2) = aerocon(ii, kk, iae_NI, ibin) / mwt_aero(iae_NI) * kg2mol ! no3-
               aq(ii, kk, ibin, jj, 3) = aerocon(ii, kk, iae_AM, ibin) / mwt_aero(iae_AM) * kg2mol ! nh4+
               ! Exclude CAT1 until speciation from dust (iae_CM) becomes available
               aq(ii, kk, ibin, jj, 4) = 0.0     ! CAT1+

               aq(ii, kk, ibin, jj, 5) = 0.0 ! HCO3-  initial from dust

!    evaluate H+ and OH-
               aq(ii, kk, ibin, jj, 6) = 2.0 * aq(ii, kk, ibin, jj, 1) + &
                                         aq(ii, kk, ibin, jj, 2) - aq(ii, kk, ibin, jj, 3)
!    added to avoid negative proton concentrations.
               if (aq(ii, kk, ibin, jj, 6) <= 0.0) then
!    add residual to oh-
                  resid = 1.0e-7 - aq(ii, kk, ibin, jj, 6)
                  aq(ii, kk, ibin, jj, 6) = 1.0e-7
               else
                  resid = 0.0
               end if
!    estimate oh- from h+
               aq(ii, kk, ibin, jj, 7) = 1.0e-14 / aq(ii, kk, ibin, jj, 6) + resid
            end do
         end do

         if (isize_0 <= isize) then
            do mm = 1, icom_aq
               do jj = 1, 2
                  aq(ii, kk, isize_0, jj, mm) = aq(ii, kk, isize_0, jj, mm) * fract(ii, kk)
               end do
            end do
         end if
      end do
   end do

!  End of the initialisation

!  new stuff for bulk chemistry

!     step a) in Wanmin Gong Documentation for Aqueous phase
!     Bulk Chemistry implementation
!
!     sum initial binned aqueous-phase tracer concentrations to
!     produced initial bulk aqueous-phase tracer (eg. AQ_BULK)
!     to go into the chemistry/wet-flux cycle
!     W.G. June 30. 2000.
!
   if (ibulk == 1) then
!  We save a copy of aq_bin in aq_orig (to be used later)
      aq_orig = aq

!  computing aq_bulk from aq_bin
      do mm = 1, icom_aq
         do jj = 1, 2
            do kk = 1, pnk
               do ii = 1, pni
                  aq_bulk = 0.0
                  do ibin = int(rcrit(ii, kk)), isize   ! all activated bins are
                     aq_bulk = aq_bulk + aq(ii, kk, ibin, jj, mm)
                  end do
                  aq(ii, kk, 1, jj, mm) = aq_bulk ! save aq_bulk in first bin
               end do
            end do
         end do
      end do
   end if
!
!  Reduce dimension from (pni, pnk) to (pni*pnk)
!  similar operation for size-dependant arrays are done later

   do mm = 1, maxnsg          ! maxns is replace by maxnsg
      do kk = 1, pnk
         do ii = 1, pni
            ijk = (kk - 1) * pni + ii
            gnew(ijk, mm) = gtmp(ii, kk, mm)
            g0(ijk, mm) = gtmp(ii, kk, mm)
         end do
      end do
   end do
!
   do kk = 1, pnk
      do ii = 1, pni
         ijk = (kk - 1) * pni + ii
         tempknew(ijk) = tempk(ii, kk)
         psacwnew(ijk) = psacw(ii, kk)
!
!  compute tempi and rtemp (rtw)
         tempinew(ijk) = 1.0 / tempk(ii, kk)
         rtempnew(ijk) = mwt_air / rhom(ii, kk) * 1.0e+09 / rho_h2o
      end do
   end do
!
!  end of buffer dimension reduction
!
!  The following computations are carried out by size bins (Loop over size bins)
   if (ibulk == 1) then
      nmax = 1
   else if (ibulk == 0) then
      nmax = isize
   end if
!
   fluxnew = 0.0
!
!initialize delg
   delg = 0.0

   do ibin = 1, nmax
!
!  Bulk Chemistry implementation
!  At the end of each chemistry (DOCHEM) and Wet fluxes calculation
!  steps within the chemistry /wet/fluxes loops, compute the
!  changes in bulk aerosol tracer concentration , due to
!  each processes and accumulate them, separatley, in two
!  additionnal arrays (DAQ_CHEM, and DAQ_FLX)

      if (ibulk == 1) then
!  Initilalization of DAQ_CHEM and DAQ_FLUX
         daq_flux = 0.0
         daq_chem = 0.0
      end if

!  calculate ionic strength and assign zero liquid water to the grids with
!  ionic strength greater than 0.5 or where cloud water radii (rad1) is zero
      iaq = 0
      nq  = 0
      do kk = 1, pnk
         do ii = 1, pni
            st_ion = 0.5 * (aq(ii, kk, ibin, 1, 1) * 4.0 +                    &
                            aq(ii, kk, ibin, 1, 2) + aq(ii, kk, ibin, 1, 3) + &
                            aq(ii, kk, ibin, 1, 4) + aq(ii, kk, ibin, 1, 5))
!
            if (st_ion > 0.01 .or. rad1(ii, kk, ibin) <= 0.0) then
               q_bin(ii, kk, ibin, 1) = 0.0
               q_bin(ii, kk, ibin, 2) = 0.0
            end if
            if ((q_bin(ii, kk, ibin, 1) > 0.0) .or. &
                (q_bin(ii, kk, ibin, 2) > 0.0)) nq = nq + 1
         end do
      end do

!  If total q_bin in the current 2-D slice for the ibin is zero,
!  then cycle to the next size bin.
      if (nq == 0) cycle

      aqnew = aqmin    ! initialisation
!  Reduce dimension from (ii, kk) to (ijk) for aq, q_bin and rad1
      do jj = 1, 2
         do kk = 1, pnk
            do ii = 1, pni
               ijk = (kk - 1) * pni + ii

!  re-initialise aq for non cloudy grids this is added in for output cloud ion conc
!  and pH, and should not afftect the cloud chemistry integration. (WG, Mar 2001)
               if (q_bin(ii, kk, ibin, jj) <= 0.0) then
                  aq(ii, kk, ibin, jj, :) = aqmin
                  aq(ii, kk, ibin, jj, 6) = 1.0e-07
                  aq(ii, kk, ibin, jj, 7) = 1.0e-07
               end if

!  kg/m^3 to g/m^3
               qnew(ijk, jj) = q_bin(ii, kk, ibin, jj) * 1000.0

               aqnew(ijk, 4, jj) = aq(ii, kk, ibin, jj, 1)   ! so4=
               aqnew(ijk, 5, jj) = aq(ii, kk, ibin, jj, 2)   ! no3-
               aqnew(ijk, 6, jj) = aq(ii, kk, ibin, jj, 3)   ! nh4+
               aqnew(ijk, 7, jj) = aq(ii, kk, ibin, jj, 4)   ! CAT1+
               aqnew(ijk, 8, jj) = aq(ii, kk, ibin, jj, 5)   ! HCO3-
               aqnew(ijk, 9, jj) = aq(ii, kk, ibin, jj, 6)   ! H+
               aqnew(ijk, 10, jj) = aq(ii, kk, ibin, jj, 7)  ! OH-
            end do
         end do
      end do
      if (ibulk == 1) aqnew_before = aqnew
!
      pgscav = 0.0
      do kk = 1, pnk
         do ii = 1, pni
            ijk = (kk - 1) * pni + ii
            if (qnew(ijk, 1) <= 0.0) cycle
!
!  Assign continuous scavenging rates for gas/part.
!  number of cloudwater drops/m**3 air
            rndrop = 2.387324146e-7 * qnew(ijk, 1) / &
                     (rad1(ii, kk, ibin) * rad1(ii, kk, ibin)**2)
             bb = bbconst * mfpl / rad1(ii, kk, ibin) + 1.0
            arg = 21171.42715 * tempknew(ijk)
             aa = 4.188790204 * mfpl * sqrt(arg)
            aaa = rndrop * aa * rad1(ii, kk, ibin) / bb
            pgscav(ijk, 1, 1)  = aaa / sq_gmwt(1)
            pgscav(ijk, 2, 1)  = aaa / sq_gmwt(2)
            pgscav(ijk, 3, 1)  = aaa / sq_gmwt(3)
            pgscav(ijk, 7, 1)  = aaa / sq_gmwt(7)
            pgscav(ijk, 8, 1)  = aaa / sq_gmwt(8)
            pgscav(ijk, 11, 1) = aaa / sq_gmwt(11)
            pgscav(ijk, 12, 1) = aaa / sq_gmwt(12)
!
!  hno3
            pgscav(ijk, 7, 2) = psacwnew(ijk) * 1.5 * 2.0
!  nh3
            pgscav(ijk, 8, 2) = psacwnew(ijk) * 1.5 * 2.0
         end do
      end do

!  sorting out grids for re-packing into snow(sn) and cloud water(cw) only grids
      ncw = 0
      do ijk = 1, nptsnz
! find the points that contains cloud water part cw
         if (qnew(ijk, 1) > 0.0) then
            ncw = ncw + 1
            iaq(ijk) = 2
         end if
      end do

!  Compute variable coefficients for cloud and snow
      do jj = 1, 2
         do ijk = 1, nptsnz
            lrt = qnew(ijk, jj) * rtempnew(ijk)
            baq(ijk, 1, jj) = 1.0e+15
            if (lrt > 0.0) then
               baq(ijk, 1, jj) = 1.0 / (lrt + smf)
            end if
            baq(ijk, 2, jj) = lrt
            baq(ijk, 3, jj) = 2.0  * baq(ijk, 1, jj)
            baq(ijk, 4, jj) = 0.25 * baq(ijk, 1, jj)
            baq(ijk, 5, jj) = 0.01 * baq(ijk, 1, jj)
         end do
      end do

!  Compute rate constants / variable coefficients
      call mach_incld_upaqr(aqnew, tempinew, raq, qnew, nptsnz, &
                            pgscav, rtempnew)
!
!  Initial equilibrium to initialize the system

      call mach_incld_soleq(gnew, aqnew, baq, raq, nptsnz, iaq, ncw)

!  step b) Compute changes in BULK aerosol tracer concentration
      if (ibulk == 1) then
         do jj = 1, maxnsaq
            do ijk = 1, nptsnz
               daq_chem(ijk, jj, 1) = daq_chem(ijk, jj, 1) +          &
                           (aqnew(ijk, jj, 1) - aqnew_before(ijk, jj, 1))
               aqnew_before(ijk, jj, 1) = aqnew(ijk, jj, 1)
            end do
         end do
      end if
!
!  Update rate constants that depends on concentration
      call mach_incld_upaqr(aqnew, tempinew, raq, qnew, nptsnz)
!  set ozone concentration to equilibrium value
      do ijk = 1, nptsnz
         if (qnew(ijk, 1) > 0.0) then
            aqnew(ijk, 12, 1) = baq(ijk, 1, 1) * raq(ijk, 4, 1) * &
                                gnew(ijk, 11) / raq(ijk, 5, 1)
         end if
      end do
!
!  Now entering in dochem -

      time = 0.0        ! this is a scalar (initial time: 0 s)
      tfinal = chm_timestep    ! advection time step passed through arguement

   !  Set out times
      tout = min(0.6, 0.1 * dtsplt)
      tin = 0.0       ! scalar parameter
      idrf = 1
!
!  Initial concentration

      do while ((time + 1.0e-02) < tfinal)

!  integrate mass transfer / chemistry for tout
!  to get chemistry started before advection

         call mach_incld_dochem(gnew, aqnew, qnew, tempknew, tempinew, baq, &
                                raq, psacwnew, tin, tout, idrf, nptsnz)
         if (chm_error_l) return

!  step b) Compute changes in BULK aerosol tracer concentration

         if (ibulk == 1) then
            do jj = 1, 2
               do mm = 1, maxnsaq
                  do ijk = 1, nptsnz
                     daq_chem(ijk, mm, jj) = daq_chem(ijk, mm, jj) +     &
                           (aqnew(ijk, mm, jj) - aqnew_before(ijk, mm, jj))

                     aqnew_before(ijk, mm, jj) = aqnew(ijk, mm, jj)
                  end do
               end do
            end do
         end if

         if ((time + dtsplt) > tfinal) dtsplt = tfinal - time
         tout = dtsplt

!  advance time
         if (idrf > 1) time = time + dtsplt
         idrf = 2
      end do

!  compute delg

      do mm = 1, maxnsg
         do ijk = 1, nptsnz
            delg(ijk, mm) = delg(ijk, mm) + (gnew(ijk, mm) - g0(ijk, mm))

            ! and reset gnew
            gnew(ijk, mm) = g0(ijk, mm)
         end do
      end do

!  Cloud water to rain conversion
!  fluxnew in moles per m3.

      if (ibulk == 1) then
! Wet flux indices (iwf) defined in mach_cam_utils_mod
         do iwf = 1, 8
            do kk = 1, pnk
               do ii = 1, pni
                  ijk = (kk - 1) * pni + ii
                  factor = 1.0e+09 / (rtempnew(ijk) * rho_h2o) * 1.0e-06
                  aqnew(ijk, iwf, 1) = max(aqmin, aqnew(ijk, iwf, 1) * &
                                           exp(-fctr(ii, kk) * chm_timestep))
                  aqnew_diff = aqnew(ijk, iwf, 1) - aqnew_before(ijk, iwf, 1)
                  fluxnew(ijk, iwf)  = fluxnew(ijk, iwf) - &
                                       (aqnew_diff * baq(ijk, 2, 1) * factor)
                  daq_flux(ijk, iwf, 1)  = daq_flux(ijk, iwf, 1) + aqnew_diff
               end do
            end do
         end do
!
!  steady state solution for H+:
         call mach_incld_steady(aqnew, nptsnz, 1)

         do ijk = 1, nptsnz
            fluxnew(ijk, 9) = fluxnew(ijk, 1) + 2.0 * fluxnew(ijk, 4) +  &
                              fluxnew(ijk, 5) - fluxnew(ijk, 6) +        &
                              fluxnew(ijk, 8) - fluxnew(ijk, 7)
            fluxnew(ijk, 9) = max(0.0, fluxnew(ijk, 9))
            daq_flux(ijk, 9, 1) = daq_flux(ijk, 9, 1) +  &
                                  (aqnew(ijk, 9, 1) - aqnew_before(ijk, 9, 1))
         end do
      end if

!  Return HSO3-, H2O2(aq) and ROOH(aq) to SO2, H2O2(g) and ROOH(g)

      do jj = 1, 2
         do mm = 1, 3
            do ijk = 1, nptsnz
               delg(ijk, mm) = delg(ijk, mm) + aqnew(ijk, mm, jj) * baq(ijk, 2, jj)
            end do
         end do
      end do

!  update aq (with

      do jj = 1, 2
         do kk = 1, pnk
            do ii = 1, pni
               ijk = (kk - 1) * pni + ii
               aq(ii, kk, ibin, jj, 1) = max(aqmin, aqnew(ijk, 4, jj))
               aq(ii, kk, ibin, jj, 2) = max(aqmin, aqnew(ijk, 5, jj))
               aq(ii, kk, ibin, jj, 3) = max(aqmin, aqnew(ijk, 6, jj))
            end do
         end do
      end do

   end do

!  end of the size-bin loop

   if (ibulk == 1) then
      nmax = 1
   else if (ibulk == 0) then
      nmax = isize
   end if

!  update gnew after integration over all bins
   do mm = 1, maxnsg
      do ijk = 1, nptsnz
         gnew(ijk, mm) = gnew(ijk, mm) + delg(ijk, mm)
      end do
   end do

! bring buffers back to original format
!      (pni*pnk) -> (pni, pnk)

   do mm = 1, maxnsg
      do kk = 1, pnk
         do ii = 1, pni
            ijk = (kk - 1) * pni + ii
            gtmp(ii, kk, mm) = gnew(ijk, mm)
         end do
      end do
   end do

! initialise the (output) wet flux
   flux = 0.0

   do iwf = 1, 9
      do kk = 1, pnk
         do ii = 1, pni
            ijk = (kk - 1) * pni + ii
            flux(ii, kk, iwf) = fluxnew(ijk, iwf)
         end do
      end do
   end do

!  step c) in Wanmin Gong Documentation for Aqueous phase
!  Bulk Chemistry implementation
!  Update the binned auqueous array with DAQ_CHEM distributed
!  over the bins according to VFRAC

   if (ibulk == 1) then
      do jj = 1, 2
         do kk = 1, pnk
            do ii = 1, pni
               ijk = (kk - 1) * pni + ii
               daq_chem_f(ii, kk, jj, 1) = daq_chem(ijk, 4, jj)
               daq_flux_f(ii, kk, jj, 1) = daq_flux(ijk, 4, jj)
               daq_chem_f(ii, kk, jj, 2) = daq_chem(ijk, 5, jj)
               daq_flux_f(ii, kk, jj, 2) = daq_flux(ijk, 5, jj)
               daq_chem_f(ii, kk, jj, 3) = daq_chem(ijk, 6, jj)
               daq_flux_f(ii, kk, jj, 3) = daq_flux(ijk, 6, jj)
            end do
         end do
      end do

      do mm = 1, 3
         do jj = 1, 2
            do ibin = 1, isize
               do kk = 1, pnk
                  do ii = 1, pni
                     aq_orig(ii, kk, ibin, jj, mm) = aq_orig(ii, kk, ibin, jj, mm) +  &
                                                     vfrac(ii, kk, ibin) * daq_chem_f(ii, kk, jj, mm)
                  end do
               end do
            end do
         end do
      end do
!
! *** distribute bulk cloud-to-rain over bins by concentration
! *** ratio instead of liquid water volume ratio
! *** (WG, Oct. 2001)
!
      do mm = 1, 3
         do kk = 1, pnk
            do ii = 1, pni
               isize_0 = int(rcrit(ii, kk))
               aq_bulk = 0.0
               do ibin = isize_0, isize
                  aq_bulk = aq_bulk + aq_orig(ii, kk, ibin, 1, mm)
               end do
               aqfrc(ii, kk, mm) = daq_flux_f(ii, kk, 1, mm) / (aq_bulk + 0.1 * aqmin)
               aqfrc(ii, kk, mm) = min(0.0, max(-1.0, aqfrc(ii, kk, mm)))
            end do
         end do
      end do
!
      do mm = 1, 3
         do jj = 1, 2
            do ibin = 1, isize
               do kk = 1, pnk
                  do ii = 1, pni
                     aq(ii, kk, ibin, jj, mm) = aq_orig(ii, kk, ibin, jj, mm) * (aqfrc(ii, kk, mm) + 1.0)
                     aq(ii, kk, ibin, jj, mm) = max(0.0, aq(ii, kk, ibin, jj, mm))
                  end do
               end do
            end do
         end do
      end do

      do kk = 1, pnk
         do ii = 1, pni
            qbin_bulk = q_bin(ii, kk, 1, 1)
!  calculate rain water (mol per m3) due to cloud-to-rain
!  Note: qbin_bulk is in kg/m3_air.
            flux(ii, kk, 10) = qbin_bulk * fctr(ii, kk) *  chm_timestep * 1000.0 / 18.0

            remfrc = min(0.0, max(-1.0, aqfrc(ii, kk, 1)))
            isize_0 = int(rcrit(ii, kk))
            do ibin = isize_0, isize
!
!  adjust number concentration for cloud-to-rain removal
!  and account for wet removal of non aqueous-phase aerosol
!  components (internally mixed), starting from seasalt according to
!  the odering defined/stated in mach_cam_utils_mod .
!  Note: aerocon is in kg/kg, and flux in moles per m3 air.
               aeronum(ii, kk, ibin) = (1.0 + remfrc) * aeronum(ii, kk, ibin)
               do iae = 1, icom
                  iwf = ip_wflx(iae)
                  if (iwf <= 10) cycle    ! skip SO4(=), NO3(-), and NH4(+)
                  flux(ii, kk, iwf) = flux(ii, kk, iwf) - remfrc * &
                                     aerocon(ii, kk, iae, ibin) * rhom(ii, kk) / mwt_aero(iae)
                  aerocon(ii, kk, iae, ibin) = (1.0 + remfrc) * aerocon(ii, kk, iae, ibin)
               end do
            end do
         end do
      end do
   end if

!  buffer back to original format completed-

!  Update aerosol concentration after aqueous-phase processes

   do kk = 1, pnk
      do ii = 1, pni
         isize_0 = int(rcrit(ii, kk))
         do ibin = isize_0 + 1, isize
            q_bin_t = q_bin(ii, kk, ibin, 1) + q_bin(ii, kk, ibin, 2)
!  convert = den_air * den_water / lwc / mwt
!  (kg/m^3a) (kg/m^3w)  (kg/m^3a)
!  note => kmol/m^3 = mol/l, rhom in g/m^3 and den_water = 1000.
            if (q_bin_t <= 0) cycle
            convert = rhom(ii, kk) / q_bin_t
            aerocon(ii, kk, iae_SU, ibin) = (aq(ii, kk, ibin, 1, 1) + aq(ii, kk, ibin, 2, 1)) / convert * mwt_aero(iae_SU)
            aerocon(ii, kk, iae_NI, ibin) = (aq(ii, kk, ibin, 1, 2) + aq(ii, kk, ibin, 2, 2)) / convert * mwt_aero(iae_NI)
            aerocon(ii, kk, iae_AM, ibin) = (aq(ii, kk, ibin, 1, 3) + aq(ii, kk, ibin, 2, 3)) / convert * mwt_aero(iae_AM)
! no updating for CM for now
!      aerocon(ii, kk, iae_SD, ibin) = (aq(ii, kk, ibin, 1, 4) + aq(ii, kk, ibin, 2, 4)) / &
!                                            convert * mwt_aero(iae_SD) * 4.0
         end do
         if (isize_0 > isize) cycle
         q_bin_t = q_bin(ii, kk, isize_0, 1) + q_bin(ii, kk, isize_0, 2)
         if (q_bin_t <= 0) cycle
         convert = rhom(ii, kk) / q_bin_t
         aerocon(ii, kk, iae_SU, isize_0) = aerocon(ii, kk, iae_SU, isize_0) * (1.0 - fract(ii, kk)) + (aq(ii, kk, isize_0, 1, 1) +  &
                                     aq(ii, kk, isize_0, 2, 1)) / convert * mwt_aero(iae_SU)
         aerocon(ii, kk, iae_NI, isize_0) = aerocon(ii, kk, iae_NI, isize_0) * (1.0 - fract(ii, kk)) + (aq(ii, kk, isize_0, 1, 2) +  &
                                     aq(ii, kk, isize_0, 2, 2)) / convert * mwt_aero(iae_NI)
         aerocon(ii, kk, iae_AM, isize_0) = aerocon(ii, kk, iae_AM, isize_0) * (1.0 - fract(ii, kk)) + (aq(ii, kk, isize_0, 1, 3) +  &
                                     aq(ii, kk, isize_0, 2, 3)) / convert * mwt_aero(iae_AM)

!      aerocon(ii, kk, iae_SD, isize_0) = aerocon(ii, kk, iae_SD, isize_0) *  &
!             (1.0 - fract(ii, kk)) + (aq(ii, kk, isize_0, 1, 4) +  &
!              aq(ii, kk, isize_0, 2, 4)) / convert * mwt_aero(iae_SD) * 4.0
      end do
   end do

!  Return gaz_conc(pni, maxnsg, pnk) to gaz_conc(pni, maxns, pnk) with update.
!  The original gaseous species concentrations were stored in gtmp.

   do mm = 2, maxns
      ij = ipos_g(mm)
      do kk = 1, pnk
         do ii = 1, pni
            gaz_conc(ii, kk, mm) = 1.0e-06 * gtmp(ii, kk, ij) * (mwt_igs(mm) / mwt_air)
         end do
      end do
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_incld [END]')
   if (chm_timings_L) call timing_stop_omp(350)
   !-----------------------------------------------------------------

   return

end subroutine mach_incld_main
