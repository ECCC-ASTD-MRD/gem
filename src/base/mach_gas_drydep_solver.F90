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
! Projet / Project : GEM-MACH
! Fichier / File   : mach_gas_drydep_solver.ftn90
! Creation         : A. Robichaud - Dec 2002
! Description      : new scheme for dry deposition of gas based on multiple
!                    resistance approach and single-layer "big leaf" approach
!                    (mostly based on work of WESELEY (1989,1996), ZHANG et al.
!                    2002 and ROBICHAUD, 1991 and 1994).
!
! Extra Info       : Modified by Bin He, Modify for thread safety - Oct 2003
!
!                    Modified by A. Kallaur, M. Moran, P.A. Beaulieu, A. Robichaud for GEM-MACH - Feb 2008
!
!                    Modified by S.R.Beagley - June 2013 (snow fraction usage)
!                                                         sea-ice and open water used where needed.
!
!                    Modified by Verica Savic-Jovcic per instructions of Paul Makar - Nov 2017
!                                (limitted soil resistance and
!                                 mesophyll resistance included in calculations of surface resistance)
!
!                    Added option to use satellite-derived LAI field
!                    - P. Makar & A. Akingunola, Feb. 2020
!
!                    Added option to use parameterization of marine halogen chemistry following Sarwar et al. (2015)
!                    - D. Pendelbery, K. Toyota & V. Savic-Jovcic, March 2020
!
!                    For a good review of dry deposition modeling, see
!                    JACOBSON, 1999, "Fundamentals of Atmospheric Modeling",
!                    WESELEY and HICKS, 2000, "Atmospheric Environment" or
!                    Seinfeld and Pandis, 1998 "Atmospheric Chemistry and Physics".
!
!
! Arguments:  IN
!
!    metvar2d(:, MV2D_ILMO)  -> Inverse of Monin-Obukhov length
!
!    metvar2d(:, MV2D_TSURF)   -> Surface temperature (K)
!                NOTE: Two choices here:
!                      1) Take TSURF; averaged skin Temp from PHY (ISBA)
!                      2) Take TT(1) 1st level model temp
!                      In discussion with B. Bilodeau Nov 22 2007, he mentioned
!                      that the most appropriate option is 1).
!
!    metvar2d(:, MV2D_UE) -> Surface friction velocity (m/s) taken from Phyvar(UE)
!
!    metvar2d(i, MV2D_QDIAG)        -> Specific humidity
!
!    metvar2d(:, MV2D_FLUSOLIS)-> Downwards visible solar flux
!                Phyvar(FLUSOLIS)  (W/m2)
!
!    metvar2d(:, MV2D_RAINRATE) -> Precipitation rate (m/s)
!                Phyvar (U1 or RAINRATE from volatile bus)
!
!     iseason -> Seasonal categories
!                  1  midsummer with lush vegetation
!                  2  autumn with cropland before harvest
!                  3  later autumn after frost, no snow
!                  4  winter, snow on ground and subfreezing
!                  5  transitional spring with partially green short annuals
!
!         lfu -> Land Form Use
!                CHEMVar(LAND_USE_15) from Chemical permanent bus
!                Derived from the CMC26 category data set into 15 category
!                set in subroutine "mach_landuse"
!
!    metvar2d(:, MV2D_SNODP)  -> Snow depth (m)
!                PhyVar(SD)
!
!    metvar2d(:,MV2D_PPLUS)   -> Surface presssure (Pa)
!                PhyVar(2p -> pplus from Dyn bus)
!
! AERO_RESIST -> Aerodynamic resistance for deposition gas species (s/m)
!
! Arguments:     OUT
!    VD       -> Deposition velocity for deposition gas species (m/s)
!    VDG      -> Deposition velocity for ammonia through ground (m/s) (optional; required when chm_ammonia_bidi_s != 'OFF')
!
! DIFF_RESIST -> Molecular diffusion resistance for deposition gas species (s/m)
!
! SURF_RESIST -> Total surface resistance for deposition gas species (s/m)
!
!============================================================================!
!
!!if_on
subroutine mach_gas_drydep_solver(vd, aero_resist, diff_resist, surf_resist, &
                                  iseason, lfu, lai_2d, metvar2d, vdg)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm, nb_gas_depo
   use chm_ptopo_grid_mod,   only: chm_ni
!!if_off
   use chm_metvar_mod,       only: MV2D_UE,    MV2D_ILMO,     MV2D_PPLUS,    &
                                   MV2D_TSURF, MV2D_FLUSOLIS, MV2D_RAINRATE, &
                                   MV2D_SNODP, MV2D_QDIAG
   use chm_utils_mod,        only: chm_lun_out, global_debug
   use chm_nml_mod,          only: chm_o3icedep, chm_bkgd_co2, chm_dep_lai2d_l, &
                                   chm_mar_halo_l, chm_gas_drydep_s
   use chm_consphychm_mod,   only: karman, tcdk, rgasd
   use mach_drydep_mod,      only: gas_depo, nsn, prandtl, b4, ao, bo, co, &
                                   dzero, isimple, inew, insz, laindex_sat,&
                                   rcutd, rgdso2, rgdo3, rexpo3, rexpso2,  &
                                   rcanp, rsmin, tmin, tmax, topt, laindex
   use chm_species_info_mod, only: species_master
   use chm_species_idx_mod,  only: sp_SO4, sp_O3, sp_SO2, sp_NH3
   use chm_datime_mod,       only: imonth
   implicit none
!!if_on
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent(out) :: vd         (nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: aero_resist(chm_ni, lucprm)
   real(kind=4),    intent(out) :: diff_resist(nb_gas_depo, chm_ni)
   real(kind=4),    intent(out) :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent (in) :: lai_2d     (chm_ni)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
   real(kind=4), optional, intent(out) :: vdg (lucprm, chm_ni)
!!if_off
!
!  Local variables
!
   integer(kind=4) :: i, sp, specie
   real(kind=4)    :: ilmob, ksqp, ksde, ksca, kst
   real(kind=4)    :: psurf, tsurf, humr
   real(kind=4)    :: sat_h2o_pressure, srad, b7, de, tl
   real(kind=4)    :: th, t0, rcut, rsx, rconv, rexp
   real(kind=4)    :: rsoil, rcan, rground, iuss, lain, laid, rnsinv
   real(kind=4)    :: rinvrcx, precip, rmx, vdf, rmx1, surf_resist_max
   real(kind=4)    :: wst, schmidt_num, expts
   real(kind=4)    :: ustar, tsurfk, tsurfm, scpr_w
   integer(kind=4) :: isnow, isea, nlus
   logical(kind=4) :: local_dbg
   real(kind=4)    :: dh2odx(nb_gas_depo), mol_diff_coef(nb_gas_depo)
   real(kind=4)    :: ksx(lucprm), lai(lucprm), sumlai, sumvegfrac
   real(kind=4)    :: lfu1(chm_ni, lucprm)
   real(kind=4)    :: alpha_sp, beta_sp, fzero_sp, hstar_sp
! Local constants
   real(kind=4), parameter :: max_ilmob = 100.

!  Code for dry deposition of gaseous species begins here

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

!  Molecular diffusivity of gas with respect to H2O and maximum surface resistance
   select case (trim(chm_gas_drydep_s))
      case ("ROBICHAUD")
         do sp = 1, nb_gas_depo
            dh2odx(sp) = sqrt(2.608 / (1.0 + 28.9644 / &
                         species_master(gas_depo(sp) % sp_id) % mol_wt))
            schmidt_num = 0.84 * dh2odx(sp)
            mol_diff_coef(sp) = 2.0 * exp(0.6667 * log(schmidt_num / prandtl))
         end do
         surf_resist_max=9999.0
      case ("ROBICHAUD3")
         do sp = 1, nb_gas_depo
            dh2odx(sp) = (0.66344 + 0.11921 * &
                         log(species_master(gas_depo(sp) % sp_id) % mol_wt))**4
            schmidt_num = 0.624 * dh2odx(sp)
            mol_diff_coef(sp) = 2.0 * exp(0.6667 * log(schmidt_num / prandtl))
         end do
         surf_resist_max=1.0e6
   end select

!  Initialize
   vd = 0.0
   diff_resist = 0.0
   surf_resist = 0.0
   if (present(vdg)) vdg = 0.0

   lfu1 = lfu

! Begin main loop over Domain

   do i = 1, chm_ni
      psurf = metvar2d(i, MV2D_PPLUS)

    ! Converting from Kelvin to Celsius for dry deposition equations.
    ! Also, set a lower limit of -60C for temperature to prevent model
    ! crashing over Antartica which happened in the global version

      tsurf = max(metvar2d(i, MV2D_TSURF) - tcdk, -60.)
      srad  = max(26.76, metvar2d(i, MV2D_FLUSOLIS))
      ustar = max(metvar2d(i, MV2D_UE), 0.001)

      isnow    = 0
      isea     = iseason(i)

    ! Check for snow depth at particular grid point, and
    ! override Season value if there is.
      if (metvar2d(i, MV2D_SNODP) > 0.005) then
         isnow = 1
         isea  = 4
         expts = 1000.0 * exp(-tsurf - 4.0)
      end if

    ! sat_h2o_pressure from GILL (in millibar)
      sat_h2o_pressure = max(10.0 ** ((0.7859 + 0.03477 * tsurf) / (1.0 + 0.00412 * tsurf)), 1.0) !mb

      humr = min(psurf * .01 * metvar2d(i, MV2D_QDIAG) / (0.622 * sat_h2o_pressure), 1.0) ! (pa-->mb)

      b7  = ao * log(log(srad)) - bo

      de  = (1.0 - humr) * sat_h2o_pressure ! de: in millibar

    ! Jarvis scheme for water vapor
    ! ksqp = correction for solar flux, kst: correction for temp.,
    ! ksde = correction for humidity (vapor pressure),
    ! ksca = correction for CO2

      ksqp = co * log(srad) - dzero
      if (ksqp < 0.0) ksqp = 0.0

    ! Jacobson  scheme for ksqp (isimple =1 for testing purposes. Default is 0)
      if (isimple == 1) ksqp = 1.0 / (1.0 + (200.0 / srad)**2)

      ksde = 1.0 - 0.03 * de
      if (ksde < 0.0) ksde = 0.0

      ksca = 1.0 - b7 * chm_bkgd_co2
      if (chm_bkgd_co2 < 100.0) ksca = 1.0
      if (chm_bkgd_co2 > 1000.0) ksca = 0.0

    ! resistance of gases to buoyant convection
      rconv = 100.0 * (1.0 + 1000.0 / (srad + 10.0))

    ! check if precipitation is present
      wst    = 0.0
      precip = 3.6E06 * metvar2d(i, MV2D_RAINRATE)
      if ((precip > 1.0) .or. (humr >= 0.95)) wst = 0.5

    ! loop over all LUC (Land Use Categories)
      ksx = 0.0

      do nlus = 1, lucprm

         if (lfu1(i, nlus) < 0.01) cycle

         tl  = tmin(nlus)
         th  = tmax(nlus)
         t0  = topt(nlus)

         if ((tsurf > tl) .and. (tsurf < th)) then
            kst = ((tsurf - tl) * (th - tsurf) / ((t0 - tl) * (th - t0))) ** b4
         else
            kst = 0.0
         end if

       ! total correction factor for stomatal conductance
         ksx(nlus) = max(ksqp * kst * ksde * ksca, 0.0001)

      end do

    ! sumlai is the area-weighted net LAI of the grid cell, according to the
    ! deposition alogarithm's assumed land-use dependent LAI values.
      sumlai = tiny(1.0)

      if (chm_dep_lai2d_l) then

       ! LAI relative contribution for the given land use type is calculated from the
       ! lookup table, but the satellite-derived LAI is used to give the net value of
       ! each land use type, and satellite-derived LAI fractions are also used.
       ! Here, lai(nlus) are a set of values which, if area-weighted, would give the
       ! same total grid-cell LAI as the satellite value, but preserve the same relative
       ! ratios of LAI between land use types as is assumed in the deposition code.
       !
       ! Determine the relative land use fraction of all potentially vegetated surfaces
       ! for the current grid cell:
         sumvegfrac = sum(lfu1(i, 1:7)) + sum(lfu1(i, 9:11)) + lfu1(i,15)

         do nlus = 1, lucprm
            lai(nlus) = max(laindex_sat(nlus, imonth), 1.0E-6)
            sumlai = sumlai + lai(nlus) * lfu1(i, nlus)
         end do

         do nlus = 1, lucprm
          ! Note that grid cells with LU <= 0.01 are cycled in the code which follows.
          ! Also, note that sumlai is based on a minimum of 1E-06
            if (lfu1(i, nlus) >= 0.01 .and. sumvegfrac >= 0.01) then
               lai(nlus) = max(lai(nlus) / sumlai * lai_2d(i), 1.0E-6)
            end if
         end do

       ! Note that lai(nlus) as generated above is:
       ! (1) Making use of satellite-derived monthly north american relative LAI values
       ! (2) Using the product of (1) and the local land use fraction to scale the
       !     grid-cell total LAI values,
       ! (3) Divided by the area fraction of the land use type, in order to get lai values
       !     that when multiplied by the land use fraction, would add to the satellite
       !     derived values.
         sumlai = tiny(1.0)
         do nlus = 1,lucprm
            sumlai = sumlai + lai(nlus) * lfu1(i,nlus)   !  == lai_2d
         end do
!
      else

       ! LAI from a look up table and add up very small value to avoid division by zero
         do nlus = 1,lucprm
            lai(nlus) = max(laindex(nlus, isea), 1.0E-6)
            sumlai = sumlai + lfu1(i, nlus) * lai(nlus)
         end do

      end if
!
    ! Store this area-weighted net lai back in lai_2d so that it can be transferred back to the main
    ! code on output:
    ! lai_2d(i) = sumlai

      sp_loop: do sp = 1, nb_gas_depo
         specie = gas_depo(sp) % sp_id
         hstar_sp = gas_depo(sp) % hstar
         fzero_sp = gas_depo(sp) % fzero
         alpha_sp = gas_depo(sp) % alpha
         beta_sp  = gas_depo(sp) % beta

       ! Before beginning loop over land-use categories, check if specie is SO4.
       ! If so, VD is calculated and and the loop is cycled.
         if (specie == sp_SO4) then

            ilmob = sign(min(abs(metvar2d(i, MV2D_ILMO)), max_ilmob), metvar2d(i, MV2D_ILMO))
          ! For sulfates use formulation parametrization
          ! from walcek et al. (1986)

          ! stable case
            if (ilmob >= 0.0) then
               vd(sp, i) = metvar2d(i, MV2D_UE) / 500.0
          ! unstable case
            else
               vd(sp, i) = (metvar2d(i, MV2D_UE) / 500.0) *        &
                           (1.0 + exp(-0.667 * log(-300.0 * ilmob)))
            end if

            cycle sp_loop

         end if

         diff_resist(sp, i) = mol_diff_coef(sp) / (karman * ustar)
         diff_resist(sp, i) =  min(9999.0, diff_resist(sp, i))

       ! leaf mesophyl resistance
         if (hstar_sp >= 0.0 .and. fzero_sp >= 0.0) then
            rmx1 = 1.0 / (hstar_sp / 3000.0 + 100.0 * fzero_sp)
         else
            rmx1 = 0.0
         end if

       ! loop over all LUC (Land Use Categories)
         do nlus = 1, lucprm

            if (lfu1(i, nlus) < 0.01) cycle

            if (specie == sp_O3 .and. local_dbg) then  ! Ozone case
                if ((ksx(nlus) <= 0.0) .or. (ksx(nlus) > 1.5)) then
                   write(chm_lun_out, *) 'Ozone ksx( ',i,' ) = ',ksx(nlus)
                end if
            end if

          ! final stomatal resistance for gas "sp"
            rsx = rsmin(nlus, isea) / ksx(nlus)

          ! correction rsx, rcut and rmx scaled by LAI (ValMartin et al, GRL, 2014)
          ! stomatal resistance
            rsx = min((rsx * dh2odx(sp) / lai(nlus)), 9999.0)
          ! mesophyl resistance
            rmx = rmx1 / lai(nlus)

          ! cuticular resistance and resistance on leaves
            if (hstar_sp >= 0.0 .and. fzero_sp >= 0.0) then
             ! cuticular resistance
               rcut = rcutd(nlus, isea) / (lai(nlus) *                 &
                      (1.0E-5 * hstar_sp + fzero_sp))
             ! resistance on leaves, twigs,bark and other exposed surfaces (WESELEY)
               rexp = 1.0 / (1.0E-5 * hstar_sp / rexpso2(nlus, isea) + &
                      fzero_sp / rexpo3(nlus, isea))
            else
               rcut = 9999.0
               rexp = 9999.0
            end if

          ! soil resistance from AURAMS
            if ((alpha_sp > 0.0) .or. (beta_sp > 0.0)) then
               rsoil = 1.0 / (alpha_sp / rgdso2(nlus, isea) + &
                       beta_sp / rgdo3(nlus, isea))
            else
               rsoil = 9999.0
            end if

          ! canopy resistance (sometimes called rac: areodynamic canopy resis.)
            rcan = rcanp(nlus, isea)

          ! total surface resistance (see JACOBSON 1999, WESELEY 1989)
          ! this part is only for ozone: new scheme for non-stomatal for ozone
          ! ref. ZHANG et al., Atmos. Env., 36, 2002.
          ! however if inew=1 treat all gases with this scheme
            if (chm_mar_halo_l .and. (specie == sp_O3) .and. (nlus == 14)) then

               rground = 0.0
             ! the following is from Sarwar et al. (2015)
               tsurfm = max(-2.0,tsurf)
               tsurfk = tsurfm+273.15
             ! from Kenjiro (Sc/Pr)^(-2/3) for water with temperature dependency
             ! temperature depencies for Heff, Ki and Dw from Pound et al. 2020, and for Ci from Sarwar et al. (2015)
               scpr_w = 4.2737e-7 * (tsurfm*tsurfm+155.*tsurfm+3700.) * 10.**(247.8/(tsurfk-140.)) * exp(1896./tsurfk)
               rinvrcx = 1.75 * karman * sqrt(psurf/rgasd/tsurfk/1027.)*ustar/scpr_w**(0.66667) +  &
                         10**(-0.25 - 0.013*tsurfm) *  &
                         sqrt( exp(-8772.2/tsurfk + 51.5) * 1.46e6*exp(-9134./tsurfk) * 1.1e-6*exp(-1896./tsurfk) )

            else if ((insz == 1) .and. ((specie == sp_O3) .or. (inew == 1))) then

             ! the following is from ZHANG. et AL. (2002), Atmospheric Env.
             ! if rain or dew is present, switch to another equation (see, ZHANG
             ! et al., 2002).
               if ((lai(nlus) > 0.001) .and. (lai(nlus) < 10.0)) then

                  iuss = 1.0 / (ustar * ustar)
                  lain = exp(0.25 * log(lai(nlus)))
                  laid = 1.0 / lain

                  rground = rcan * iuss * lain + rsoil

                  if (wst > 0.499) then
                   ! for wet canopies
                     rcut = rcut / 20
                     rnsinv = 1.0 / rground + sqrt(lai(nlus)) * ustar / rcut
                  else
                   ! for dry canopies
                     rnsinv = 1.0 / rground + 1.0 / (rcut * exp(-0.03 * 100 * humr) * laid / ustar)
                  end if

                  rinvrcx = (1.0 - wst) / (rsx + rmx) + rnsinv

               end if

          ! adjustment for snow (see Robichaud, 1991)
          ! Only applies to SO2
            else if (specie == sp_SO2 .and. isnow == 1 .and. nlus /= 13 .and. nlus /= 14) then

               rground = rcan + rsoil + expts

               rinvrcx = 1.0 / (rcut + expts) + 1.0 / (rconv + rexp + expts) + 1.0 / rground

            else

               rground = rcan + rsoil

               rinvrcx = (1.0 - wst) / (rsx + rmx) + 1.0 / rcut + 1.0 / (rconv + rexp) + 1.0 / rground

            end if

            surf_resist(nlus, sp, i) = 1.0 / rinvrcx
            surf_resist(nlus, sp, i) = max(0.0, surf_resist(nlus, sp, i))
            surf_resist(nlus, sp, i) = min(surf_resist_max, surf_resist(nlus, sp, i))

          ! dry deposition for gas "sp" for one LUC
            vdf = 1.0 / (aero_resist(i, nlus) + diff_resist(sp, i) + surf_resist(nlus, sp, i))

            if (vdf < 0.0) then
               if (local_dbg) then
                  write(chm_lun_out, *) '>Warning'
                  write(chm_lun_out, *) '>Dry deposition speed is negative at i = ', i, ' and nlus ', nlus, ' for species index', specie
               end if
               vdf = 0.0
            end if

          ! special cases (see: WESELEY and HICKS, 2000, Atmospheric Env.)
          ! where current scheme does not work well
            if (specie == sp_O3) then

               if (nlus == 13) then
                  vdf = 3.0E-4 ! original values was vdf = 1.0E-4
               else if (nlus == 14 .and. (.not. chm_mar_halo_l))  then
                  vdf = 3.0E-4
               else if ((nlus == 12) .and. (chm_o3icedep == 1))  then
                  vdf = 1.0E-4
               else if ((nlus == 9) .and. (isnow == 0)) then
                  vdf = 1.5E-3
               end if

            end if

            vd(sp, i) = vd(sp, i) + lfu1(i, nlus) * vdf

            vd(sp, i) = min(0.1, vd(sp, i))

          ! deposition velocity contribution from ground for ammonia
            if (specie==sp_NH3 .and. present(vdg)) then
               vdg(nlus, i) = lfu1(i, nlus) * vdf * min(1.0, max(0.0, surf_resist(nlus, sp, i) / rground))
            end if

         end do      ! loop over all land use categories

      end do sp_loop ! loop over all deposition species

   end do  ! loop over grid point row

   return

end subroutine mach_gas_drydep_solver
