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
! Projet / Project : GEM-MACH version 2
! Fichier / File   : mach_drydep_gas_solver.ftn90
! Creation         : A. Robichaud - Dec 2002, UPGRADE DEC 2008 (version 2)
! Description      : new scheme for dry deposition of gas based on multiple
!                    resistance approach and single-layer "big leaf" approach
!                    (originally based on work of WESELEY (1989,1996),ROBICHAUD, 1991 and 1994,
!                     JARVIS 1976 and ZHANG ET AL. 2002 (for non-stomatal resistance))
!                     For ocean and inland water, new parametrization are given for ozone
!                     that could be extended for all other slightly soluble gases
!                     ( see Seinfield and Pandis, 2006, Wanninkhof 1992)
!
! Extra Info       : Modified by Bin He, Modify for thread safety - Oct 2003
!
!                    Modified by A. Kallaur, M. Moran, P.A. Beaulieu, A. Robichaud for GEM-MACH - Feb 2008
!
!                    Upgrade dec 2008. Alain Robichaud for:
!                    1) ustar now calculated locally (for each land use "lus") and not from GEM
!                    2) new formulation for surface resistance for ozone over ocean for soluble and slightly soluble gases
!                    3) modified roughness length for ocean and water
!                    4) minimum value for surf resis and aero_resis imposed to avoid unrealistic high values vd
!                    5) temperature dependency for solubility into water introduced: ats
!                    6) surface resistance of H2O2 does not follow Weseley's model anymore (following Ganzeveld)
!                    8) include bug fix for the computation of water vapour deficit (ksde)    
!                    9) formula for dh2odx was not precise and is replaced
!
!                    Modified by A. Robichaud (May 2011). Ref. Science , Karl et al, 2010.
!                    Rate of metabolization of Oxygenated VOCs stronger than previously thought
!
! Arguments:  IN
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
subroutine mach_gas_drydep_solver2(vd, aero_resist, diff_resist, surf_resist, &
                                   iseason, lfu, metvar2d, vdg)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm, nb_gas_depo
   use chm_ptopo_grid_mod,   only: chm_ni
!!if_off
   use chm_metvar_mod,       only: MV2D_UE,    MV2D_ILMO,     MV2D_PPLUS,    &
                                   MV2D_TSURF, MV2D_FLUSOLIS, MV2D_RAINRATE, &
                                   MV2D_SNODP, MV2D_DLAT,     MV2D_DLON,     &
                                   MV2D_QDIAG, MV2D_WSDIAG
   use chm_utils_mod,        only: chm_lun_out, global_debug
   use chm_nml_mod,          only: chm_bkgd_co2
   use chm_consphychm_mod,   only: karman, pi, tcdk
   use mach_drydep_mod,      only: gas_depo, nsn, prandtl, b4, ao, bo, co,     &
                                   dzero, isimple, inew2, insz,                &
                                   rcutd, rgdso2, rgdo3, rexpo3, rexpso2, zz0, &
                                   rcanp, rsmin, tmin, tmax, topt, laindex
   use chm_species_info_mod, only: species_master
   use chm_species_idx_mod,  only: sp_SO4, sp_O3, sp_HNO3, sp_NO, sp_NO2, &
                                   sp_TOLU, sp_H2O2, sp_SO2, sp_NH3
   implicit none
!!if_on
   integer(kind=4), intent (in) :: iseason    (chm_ni)
   real(kind=4),    intent(out) :: vd         (nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: aero_resist(chm_ni, lucprm)
   real(kind=4),    intent(out) :: diff_resist(nb_gas_depo, chm_ni)
   real(kind=4),    intent(out) :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   real(kind=4),    intent (in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent (in) :: metvar2d   (chm_ni, SIZE_MV2D)
   real(kind=4), optional, intent(out) :: vdg (lucprm, chm_ni)
!!if_off
!
!  Local variables
!
   real(kind=4)    :: ilmob, kst, ksqp, ksde, ksca, fctdsol, rmwt
   real(kind=4)    :: psim, tsurf, humr, href1, vvent, schmidt_num
   real(kind=4)    :: sat_h2o_pressure, srad, b7, de, tl
   real(kind=4)    :: th, t0, rcut, rsx, rconv, rexp
   real(kind=4)    :: rsoil, rcan, rground, iuss, lain, laid, lai, rnsinv
   real(kind=4)    :: rinvrcx, precip, rmx, rmx1, vdf
   real(kind=4)    :: wst, psim1, psim2, expts
   real(kind=4)    :: roughlen, ref_hgt, ylat, ylon
   real(kind=4)    :: hatil, kinvw, dynvw, dmx, scxw, sco2, sa
   real(kind=4)    :: klm, rth1, rth2, rth4
   real(kind=4)    :: rad_to_deg 
   integer(kind=4) :: i, ii, sp, specie
   integer(kind=4) :: isnow, isea, nlus, iocnew
   logical(kind=4) :: local_dbg
   real(kind=4)    :: ustar(lucprm), ksx(lucprm)
   real(kind=4)    :: dh2odx(nb_gas_depo), mol_diff_coef(nb_gas_depo)
   real(kind=4)    :: alpha_sp, beta_sp, fzero_sp, hstar_sp, ats_sp
! Local constants
   real(kind=4), parameter :: max_ilmob = 100.

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

!  Code for dry deposition of gaseous species begins here

! Molecular diffusivity of gas  with respect to H2O
   do sp = 1, nb_gas_depo
      rmwt = species_master(gas_depo(sp) % sp_id) % mol_wt
      dh2odx(sp) = 0.472 + 0.033 * rmwt - 0.0002 * rmwt**2 + 748E-9 * rmwt**3
      schmidt_num = 0.68 * dh2odx(sp)
      mol_diff_coef(sp) = 2.0 * exp(0.6667 * log(schmidt_num / prandtl))
   end do

!  Initialize deposition velocity table
   vd = 0.0
   diff_resist = 0.0
   surf_resist = 0.0
   if (present(vdg)) vdg = 0.0

   rad_to_deg = 180.0 / pi

! Begin main loop over Domain
   do i = 1, chm_ni
      ylat  = metvar2d(i, MV2D_DLAT) * rad_to_deg
      ylon  = metvar2d(i, MV2D_DLON) * rad_to_deg

      tsurf = max(metvar2d(i, MV2D_TSURF) - tcdk, -60.0)
      srad  = max(26.76, metvar2d(i, MV2D_FLUSOLIS))
      vvent = metvar2d(i, MV2D_WSDIAG)

      isnow    = 0
      isea     = iseason(i)

      if (metvar2d(i, MV2D_SNODP) > 0.005) then
         isnow = 1
         isea  = 4
         expts   = 1000.0 * exp(-tsurf - 4.0)
      end if
      if ( metvar2d(i, MV2D_SNODP) < .001 .and. isea == 4 .and. tsurf > 5.) then
         isea  = 5    ! transition season
      end if
      if (ylat >= 60.0 .and. ylat <= 75.0) then
         if (isea == 4 .and. tsurf >= 5.0) isea = 5
      end if

      sat_h2o_pressure = max(10.0 **((0.7859 + 0.03477 * tsurf) / (1.0 + 0.00412 * tsurf)), 1.0) !mb
      humr = min(metvar2d(i, MV2D_PPLUS) * .01 * metvar2d(i, MV2D_QDIAG) / (0.622 * sat_h2o_pressure), 1.0) ! (pa-->mb) 

      b7  = ao * log(log(srad)) - bo

      de  = (1.0 - humr) * sat_h2o_pressure ! de: in millibar

   ! Jarvis scheme for water vapor
   ! ksqp = correction for solar flux, kst: correction for temp.,
   ! ksde = correction for humidity (vapor pressure),
   ! ksca = correction for CO2

      ksqp = co * log(srad) - dzero
      if (ksqp < 0.0) ksqp = 0.0
      if (isimple == 1) ksqp = 1.0 / (1.0 + (200.0 / srad)**2)

      ksde = 1.0 - 0.03 * de 
      if (ksde < 0.0) ksde = 0.0

      ksca = 1.0 - b7 * chm_bkgd_co2
      if (chm_bkgd_co2 < 100.0) ksca = 1.0
      if (chm_bkgd_co2 > 1000.0) ksca = 0.0

    ! resistance of gases to buoyant convection
      rconv = 100.0 * (1.0 + 1000.0 / (srad + 10.0))

      wst    = 0.0
      precip = 3.6E06 * metvar2d(i, MV2D_RAINRATE)
      if ((precip > 1.0) .or. (humr >= 0.95)) wst = 0.5

      ilmob = sign(min(abs(metvar2d(i, MV2D_ILMO)), max_ilmob), metvar2d(i, MV2D_ILMO))

! loop over all LUC (Land Use Categories)
      do nlus = 1, lucprm

         if (lfu(i, nlus) < 0.01) cycle
     
         ref_hgt = 10.0 * zz0(nlus, isea)
         ref_hgt = max(1.5, ref_hgt)
         if (nlus == 13 .or. nlus == 14) ref_hgt = 10.0

         roughlen = zz0(nlus, isea)

         psim = 0.0

         href1 = max((ref_hgt - roughlen), 1.5)
         if ((ref_hgt * ilmob) > 0.01) then
            psim = -5.0 * (href1 * ilmob)
         else if ((ref_hgt * ilmob) < -0.01) then
            psim1 = exp(0.032 + 0.448 * log(-ref_hgt * ilmob) - 0.264 * log(-1.0 * ref_hgt * ilmob))
            psim2 = exp(0.032 + 0.448 * log(-roughlen * ilmob) - 0.264 * log(-1.0 * roughlen * ilmob)) 
            psim  = psim1 - psim2
         end if

         ustar(nlus) = max(karman * vvent / (log(ref_hgt / roughlen) - psim), 0.001)
         if (nlus == 13 .or. nlus == 14) then
            do ii = 1,10
               if (nlus == 13) roughlen = 0.11 * 1.5E-5 / ustar(nlus)
               if (nlus == 14) roughlen = 0.0144 *ustar(nlus)**2 / 9.81 + 1.463e-5/(9.1 * ustar(nlus)) 
               if ((ref_hgt * ilmob) < -0.01) then
                  psim1 = exp(0.032 + 0.448 * log(-ref_hgt * ilmob) - 0.264 * log(-1.0 * ref_hgt * ilmob))
                  psim2 = exp(0.032 + 0.448 * log(-roughlen * ilmob) - 0.264 * log(-1.0 * roughlen * ilmob)) 
                  psim  = psim1 - psim2
               end if
               ustar(nlus) = karman * vvent / (log(ref_hgt / roughlen) - psim)
            end do
         end if

         tl  = tmin(nlus)
         th  = tmax(nlus)
         t0  = topt(nlus)
    
     ! kst: correction for temp.,
         if ((tsurf > tl) .and. (tsurf < th)) then
            kst = ((tsurf - tl) * (th - tsurf) / ((t0 - tl) * (th - t0))) ** b4
         else
            kst = 0.0
         end if

     ! total correction factor for stomatal conductance
         ksx(nlus) = max(ksqp * kst * ksde * ksca, 0.1)
      end do

      sp_loop: do sp = 1, nb_gas_depo
         specie = gas_depo(sp) % sp_id
         hstar_sp = gas_depo(sp) % hstar2
         fzero_sp = gas_depo(sp) % fzero2
         alpha_sp = gas_depo(sp) % alpha
         beta_sp  = gas_depo(sp) % beta
         ats_sp   = gas_depo(sp) % ats

!  Before beginning loop over land-use categories, check if specie is SO4.
!  If so, VD is calculated and and the loop is cycled.
         if (specie == sp_SO4) then

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

   ! leaf mesophyl resistance
         fctdsol = exp(ats_sp * (1.0 / (tsurf + tcdk) - 1.0 / 298.0))
         if (hstar_sp >= 0.0 .and. fzero_sp >= 0.0) then
            rmx1 = 1.0 / (hstar_sp * fctdsol / 3000.0 + &
                   100.0 * fzero_sp)
         else
            rmx1 = 9999.0
         end if

! loop over all LUC (Land Use Categories)
         do nlus = 1, lucprm

            if (lfu(i, nlus) < 0.01) cycle
     
     ! LAI from a look up table and add up very small value to avoid division by zero
            lai = min(laindex(nlus, isea), 1.0E-6)

    ! canopy resistance (sometimes called rac: areodynamic canopy resis.)
            rcan = rcanp(nlus, isea)

    ! soil resistance from AURAMS
            rsoil = 1.0 / (alpha_sp / rgdso2(nlus, isea) + &
                    beta_sp / rgdo3(nlus, isea))

            diff_resist(sp, i) = mol_diff_coef(sp) / &
                                 (karman * (ustar(nlus) + 0.001))
            diff_resist(sp, i) = min(9999.0, diff_resist(sp, i))

            if (local_dbg .and. specie == sp_O3) then  ! Ozone case
               if ((ksx(nlus) <= 0.0) .or. (ksx(nlus) > 1.5)) then
                  write(chm_lun_out, *) 'Ozone ksx( ',i,' )= ',ksx(nlus)
               end if
            end if

     ! final stomatal resistance for gas "sp_i"
            rsx = rsmin(nlus, isea) / ksx(nlus)
     ! correction rsx, rcut and rmx  scaled by LAI (ValMartin et al, GRL, 2014)
            rsx = min((rsx * dh2odx(sp) / lai), 9999.0)
            rmx = rmx1 / lai
     ! cuticular resistance and resistance on leaves
            if (hstar_sp >= 0.0 .and. fzero_sp >= 0.0) then
     ! cuticular resistance
               rcut = rcutd(nlus, isea) /(lai * &
                      (1.0E-5 * hstar_sp * fctdsol + fzero_sp))
    ! resistance on leaves, twigs,bark and other exposed surfaces (WESELEY)
               rexp = 1.0 /                     &
                      (1.0E-5 * hstar_sp * fctdsol / rexpso2(nlus, isea) + &
                      fzero_sp / rexpo3(nlus, isea))
            else
               rcut = 9999.0
               rexp = 9999.0
            end if

         ! total surface resistance (see JACOBSON 1999, WESELEY 1989)
         ! inew2=0,insz=1  : calculate non-stomatal resistance for ozone only
         ! inew2=1         : calculate non-stomatal resistance for all species
            if ((inew2 == 1) .or. ((specie == sp_O3) .and. (insz == 1))) then

                if ((lai > 0.001) .and. (lai < 10.0)) then

                   iuss = 1.0 / (ustar(nlus) * ustar(nlus))
                   lain = exp(0.25 * log(lai))
                   laid = 1.0 / lain

                   rground = rcan * iuss * lain + rsoil
                   
                   if (wst > 0.499) then
                  ! for wet canopies
                      rcut = rcut / 20.0
                      rnsinv = 1.0 / rground + sqrt(lai) * ustar(nlus) / rcut
                   else
                  ! for dry canopies
                      rnsinv = 1.0 / rground + 1.0 / (rcut * exp(-0.03 * 100.0 * humr) * laid / ustar(nlus))
                   end if
                   
                   rinvrcx = (1.0 - wst) / rsx + rnsinv

                end if

            ! adjustment for snow (see for ex. figure 3, Robichaud 1991 for SO2)
            else if (specie == sp_SO2 .and. isnow == 1 .and. (nlus /= 13 .and. nlus /= 14)) then

               rground = rcan + rsoil + expts
               
               rinvrcx = 1.0 / (rcut + expts) + 1.0 / (rconv + rexp + expts) + &
                         1.0 / rground

            else

               rground = rcan + rsoil
               
               rinvrcx = (1.0 - wst) / (rsx + rmx) + 1.0 / rcut + 1.0 / (rconv + rexp) + 1.0 / rground
               
            end if

            surf_resist(nlus, sp, i) = 1.0 / rinvrcx

!&&&&&&&&&&&&&&&&&&&   ocean block &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            iocnew = 0

            if ((nlus == 13 .or. nlus == 14) .and. (specie /= sp_HNO3)) then  ! boucle A
               if (iocnew == 1) then     ! boucle B

                   if (hstar_sp <= 1.0 .and. hstar_sp > 0.) then  ! boucle C

                      hatil = 3.0 * (hstar_sp / 1.3E-2)
                      kinvw = 1.0E-6 * (1.77849 - 0.05472 * tsurf + 9.0205E-4 * tsurf**2 - 5.93E-6 * tsurf**3)
         ! correction for averaged ocean salinity (based on table 11-8, Knauss, 1978, Introd. Phys. Ocean)
                      if (nlus == 14) kinvw = kinvw * 1.06
                      dynvw = kinvw * 1000.0
                      dmx = 5.95E-15 * tsurf / dynvw
                      if (specie == sp_NO .or. specie == sp_NO2) dmx = 0.75 * dmx
                      if (specie == sp_TOLU) dmx = 0.5 * dmx

                      scxw = kinvw / dmx
                      sco2 = 2073.1 - 125.62 * tsurf + 3.26276 * tsurf**2 - 0.043219 * tsurf**3
                      sa = sco2 / scxw

                      if (vvent <= 3.6)  then 
                         klm = 0.17 * vvent * sa**0.67 / 3600.0
                      else if (vvent > 13.) then
                         klm = (0.612 * sa**0.67 + (5.9 * vvent - 49.9) * sa**0.5) / 3600.0;
                      else
                         klm = (0.612 * sa**0.67 + (2.85 * vvent - 10.26) *sa**0.5) / 3600.0;
                      endif
        
                      if (((ylat >= 35.0 .and. ylat <= 65.0 ).and.(ylon >= 310.0 .and. ylon < 355.0))  .or.   &
                  ((ylat >= 35.0 .and. ylat <= 65.0) .and. (ylon >= 150.0 .and. ylon < 220.0))) then
                         if (isea == 4) then 
         ! W&M 1999 applies
                            klm = (0.0283 * vvent**3) / ((sa**0.5) * 3600.0)
                         else if (isea == 3 .or. isea == 5) then
         ! W92 applies
                            klm = 0.31 * vvent**2 * sa**0.5 / 3600.0
                         end if
                      else if (ylat <= -35.0) then
                         klm = (0.0283 * vvent**3) / ((sa**0.5) * 3600.0)
                      end if
         ! surface resistance over ocean (convert into s/m units)

                      surf_resist(nlus, sp, i) = 100.0 / (klm * hatil * fctdsol)
                      if (nlus == 14) then
                         rth1 = 5000.0 * 1.3E-2 / hstar_sp  ! for other gases than ozone scaled by henry solubility constant hstar2
                         if (surf_resist(nlus, sp, i) >= rth1) surf_resist(nlus, sp, i) = rth1
                      else
                         rth2 = 3300.0 * 1.3E-2 / hstar_sp
                         rth4 = 1000.0 * 1.3E-2 / hstar_sp
                         if (surf_resist(nlus, sp, i) >= rth2) surf_resist(nlus, sp, i) = rth2
                         if (surf_resist(nlus, sp, i) <= rth4) surf_resist(nlus, sp, i) = rth4
                      end if

                   else if (hstar_sp <= 100.0 ) then  ! moderately  soluble gases 
                      surf_resist(nlus, sp, i) = rsoil
                   else
                      surf_resist(nlus, sp, i) = 5.0
                   end if ! fin boucle C 

               else  ! case when a constant for rc is imposed 

                   surf_resist(nlus, sp, i) = 3300.0 * 1.3E-2 * (dh2odx(sp) / 1.63) / hstar_sp
                   
               end if ! fin boucle B

               rground = surf_resist(nlus, sp, i)
                
            end if ! fin boucle A
!&&&&&&&&&&&&&&&&&&&   end of ocean block &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

            if (specie == sp_H2O2) surf_resist(nlus, sp, i) = 5.0
                  
            if ((specie == sp_O3) .and. (nlus == 8) .and. (srad <= 75.0)) surf_resist(nlus, sp, i) = 9999.0
            surf_resist(nlus, sp, i) = min(9999.0, surf_resist(nlus, sp, i))
            surf_resist(nlus, sp, i) = max(5.0, surf_resist(nlus, sp, i))

            vdf = 1.0 / (aero_resist(i, nlus) + diff_resist(sp, i) + surf_resist(nlus, sp, i))

            if (vdf < 0.0) then
!           if(chm_lun_out > 0) write(chm_lun_out, *) '>Warning'
!           if(chm_lun_out > 0) write(chm_lun_out, *) '>Dry deposition speed is negative at i = ', i
               vdf = 0.0
            end if

            if (specie == sp_O3)  then
               if (((nlus == 9) .or. (nlus ==10)) .and. (isnow == 0)) then 
                  vdf = min(1.5E-3, vdf)
               end if
               if ((nlus == 8) .and. (srad < 75.0)) then
                  vdf = 0.1E-3;
               end if
!           if ((ylat .le. 40.0) .and. (srad < 75.0) .and. (ylon.lt.160.0) .and. (nlus .le.12)) vdf=1.0E-5 
            end if
             
            if (specie == sp_HNO3 .and. ((nlus == 13) .or. (nlus == 14))) then
               vdf = 1.0E-2
            end if

         ! computation for all landuses
            vd(sp, i) = vd(sp, i) + lfu(i, nlus) * vdf

            vd(sp, i) = min(0.1, vd(sp, i))

         ! deposition velocity contribution from ground for ammonia
            if (specie==sp_NH3 .and. present(vdg)) then
               vdg(nlus, i) = lfu(i, nlus) * vdf * min(1.0, max(0.0, surf_resist(nlus, sp, i) / rground))
            end if
            
         end do      ! loop over all land use categories
      end do sp_loop ! loop over all species

   end do  ! loop over grid point row
   
   return
end subroutine mach_gas_drydep_solver2
