!---------------------------------- licence begin -------------------------------
! gem-mach - atmospheric chemistry library for the gem numerical atmospheric model
! copyright (c) 2007-2013 - air quality research division &
!                           national prediction operations division
!                           environnement canada
! this library is free software; you can redistribute it and/or
! modify it under the terms of the gnu lesser general public
! license as published by the free software foundation; either
! version 2.1 of the license, or (at your option) any later version.
!
! this library is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the gnu
! lesser general public license for more details.
!
! you should have received a copy of the gnu lesser general public
! license along with this library; if not, write to the free software
! foundation, inc., 51 franklin street, fifth floor, boston, ma  02110-1301  usa
!---------------------------------- licence end ---------------------------------
!============================================================================!
!         environnement canada         |        environment canada           !
!                                      |                                     !
! - service meteorologique du canada   | - meteorological service of canada  !
! - direction generale des sciences    | - science and technology branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! projet/project : gem-mach
! fichier/file   : mach_gas_messy.ftn90
! creation       : Jack Chen - 2016-Feb
!
! description    : entry point for photolysis rate (J) calculation to MESSYJVAL module
!                  reaction mapping is HARDCODED between MESSY and gas mechanisms!
!                  output J from MESSY in unit 1/sec (converted to 1/min for ADOM2)
!                  MESSY module src (as Jun2016 renamed to ->) :
!                      messy_jval.ftn90   -> "mach_messy_jval_mod.ftn90"
!                      messy_jval_jvpp.inc-> "mach_messy_jval_inc.ftn90"
!                      messy_module.ftn90 -> "mach_pkg_messy_mod.ftn90"
!
!                  See "mach_messy_jval_mod.ftn90" for changes from original MESSY code
!                  See comments in messy_module.ftn90 for additional photolysis rates
!
! arguments:  in
!               metvar3d
!               metvar2d
!               p2d      ->
!               o3vmr    -> ozone in mol/mol
!
!             out
!               rj( gni*gnk, nrjxs )-> slab photolysis rates for njrxs [1/min]
!
!============================================================================
!
!!if_on
subroutine mach_gas_messy(gmetvar3d, metvar2d, o3vmr, rj, nmod, gni, gnk)
   use chm_ptopo_grid_mod,  only: chm_ni
   use chm_metvar_mod,      only: SIZE_MV2D, SIZE_MV3D
   use mach_pkg_gas_mod,    only: njrxs
!!if_off
   use chm_utils_mod,       only: ik, CHM_MSG_DEBUG, dp
   use chm_metvar_mod,      only: MV2D_CANG,  MV3D_TPLUS, MV2D_PPLUS,  MV3D_HUPLUS, &
                                  MV3D_LWC,   MV2D_MG,    MV3D_CLDRAD, MV2D_AL5,    &
                                  MV3D_SIGM,  MV3D_SIGT,  MV2D_H,     MV3D_ZMOM,    &
                                  MV2D_DLAT,   MV2D_DLON
   use chm_consphychm_mod,  only: du_o3, grav, avno, mwt_air
   use chm_nml_mod,         only: chm_timings_L, chm_pkg_gas_s
   use chm_datime_mod,      only: ijul_day, chm_dttim_s
   use mach_pkg_messy_mod,  only: u0lim, messy_njx, &
                                  jval_NO2, jval_NOO2,  jval_NO2O, jval_O3P,  &
                                  jval_O1D, jval_HONO,  jval_HNO3, jval_HNO4, &
                                  jval_H2O2,jval_CH3OOH,jval_CHOH, jval_COH2, &
                                  jval_CH3CHO,jval_CH3COCH3, jval_MGLYOX,     &
                                  jval_CH4, jval_N2O, jval_N2O5, jval_NO,     &
                                  jval_NO3NOO, jval_O2, jval_PAN
   use mach_messy_jval_mod, only: mach_messy_jvalues
   use mach_headers_mod,    only: mach_suncos

   implicit none
!
!!if_on
!  IO variable
   integer(kind=4),                                 intent(in)  :: gni, gnk
   integer(kind=4), dimension(gni),                 intent(in)  :: nmod
   real(kind=4),    dimension(chm_ni,   SIZE_MV2D), intent(in)  :: metvar2d
   real(kind=4),    dimension(gni, gnk, SIZE_MV3D), intent(in)  :: gmetvar3d
   real(kind=4),    dimension(gni, gnk)           , intent(in)  :: o3vmr
   real(kind=8),    dimension(gni* gnk, njrxs)    , intent(out) :: rj
!!if_off

!  local variables
!  inputs to jvalues, use full gnk, note on dummy nk+1 level for colo3
   real(kind=4)   , dimension(gni, gnk+1) :: relo3  ! O3 vmr dummy top lay. [mol/mol]
   real(kind=4)   , dimension(gni, gnk+1) :: colo3  ! O3 column with dummy top lay. [mol/mol]
   real(kind=4)   , dimension(gni, gnk)   :: temp   ! K
   real(kind=4)   , dimension(gni, gnk)   :: hum    ! Specific humidity (kg/kg)
   real(kind=4)   , dimension(gni, gnk)   :: p2d    ! Pressure in thermodynamic levels (Pa)
   real(kind=4)   , dimension(gni, gnk)   :: rhum   ! relative humidity 1-100%
   real(kind=4)   , dimension(gni, gnk)   :: clp    ! cloud liquid water path g/m2
   real(kind=4)   , dimension(gni, gnk)   :: aclc   ! cloud fraction (0-1)
   real(kind=4)   , dimension(gni     )   :: csza   ! cos solar zenith angle - recalculated
   real(kind=4)   , dimension(gni     )   :: slf    ! sea-land-frac. (0-1)
   real(kind=4)   , dimension(gni     )   :: albedo ! cumulative albedo (0-1)
   real(kind=4)   , dimension(gni     )   :: dlat   ! Latitudes
   real(kind=4)   , dimension(gni     )   :: dlon   ! Longitudes
   real(kind=4)   , dimension(gni, gnk)   :: mh     ! height, momentum level, grid-top (m)
   integer(kind=4), dimension(gni     )   :: pblidx ! pbl layer index (1=top, nk=surface)
   integer(kind=4), dimension(gni     )   :: iu0    ! filter for sorting by cosine zenith angle
! midatmo-consider solar Lyman-alpha abs line (60-70km) for photolysis of CO2, HO2 and O2 see: DOI:10.1029/97GL52690
   logical(kind=4)                        :: midatmo = .false.
   real(kind=4)   , dimension(gni, gnk)   :: dz2d    ! layer thickness in momentum level (m)
   real(kind=4)                           :: liqwcin ! cloud liquid water content

! outputs from jvalues (moved from MODULE to argument for OMP>1 fix).
   real(dp)   , dimension(messy_njx, gni, gnk) :: jval ! see messy_module for MESSY_NJX
! calheat-calculate of UV heating rates by O2 and O3 output in rh_o2, rh_o3
!   logical :: calheat = .false.
!   real(kind=8)   , dimension(gni, gnk)   :: rh_o2
!   real(kind=8)   , dimension(gni, gnk)   :: rh_o3
!
   ! molec/cm2 -> DU ; 1 DU = 2.687 E+19 molecules of O3 per square meter
   real(kind=4), parameter :: cst_o3 = (1.0E-4 / grav) * &
                              (avno / (1.0E-3 * mwt_air)) * 1.0E+03 / du_o3
   real(kind=4)            :: hrnow, njul_day
   real(kind=4)            :: o3_over, psurf
   real(kind=4)            :: sig_m(gnk + 1), dp_clp(gni, gnk)
   integer(kind=4)         :: gni_day
   integer(kind=4)         :: i, k, ii, ik0, ihour, iminute
! external subroutines
   external mfohr4, msg_toall, timing_start_omp, timing_stop_omp

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_gas_messy [BEGIN]')
   if (chm_timings_L) call timing_start_omp(331, 'mach_gas_messy', 330)

!! ozone column - using the same formulation in mach_linoz_xcol (deGrandpre)
!!                (DU->mole/cm2) - add dummy one layer top-layer
   do i = 1, gni
      ii = nmod(i)

      psurf = metvar2d(ii, MV2D_PPLUS)
      sig_m(1) = gmetvar3d(i, 1, MV3D_SIGM)
      relo3(i, 1) = o3vmr(i, 1)
      do k = 2, gnk
         sig_m(k) = gmetvar3d(i, k, MV3D_SIGM)
      end do
      sig_m(gnk + 1) = 1.0
!
      o3_over = 0.
      o3_over = cst_o3 * o3vmr(i, 1) * (sig_m(1) * psurf)
      colo3(i, 1) = o3_over * du_o3 * 1.E-3
      do k = 1, gnk
         dp_clp(i, k) = (sig_m(k + 1) - sig_m(k)) * psurf
         o3_over = o3_over +  cst_o3 * o3vmr(i, k) * dp_clp(i, k)
         colo3(i, k + 1) = o3_over * du_o3 * 1.E-3

!! ozone mixing ratio - add one dummy top-layer-relo3(:1)=relo3(:,2)=o3vmr(:1)
         relo3(i, k + 1) = o3vmr(i, k)

!! Pressure (Pa)
         p2d(i, k) = gmetvar3d(i, k, MV3D_SIGT) * psurf
      end do
   end do

!-- prepare input met. variables:
   do k = 1, gnk
      do i = 1, gni
!! temperature (K)
         temp(i, k) = gmetvar3d(i, k, MV3D_TPLUS)

!! Specific Humidity (Kg/kg)
         hum(i, k) = gmetvar3d(i, k, MV3D_HUPLUS)

!! momentum height (m)
         mh(i, k) = gmetvar3d(i, k, MV3D_ZMOM)

!! cloud fraction from GEM radiative transfer calculation - from GEM "prep_cw_rad.F90"
         aclc(i, k) = gmetvar3d(i, k, MV3D_CLDRAD)
!
!! cloud liquid path follow GEM "prep_cw_rad.F90",  kg/m2->g/m2
         liqwcin = max(gmetvar3d(i, k, MV3D_LWC), 0.0) / &
                   max(gmetvar3d(i, k, MV3D_CLDRAD), 0.05)
         clp(i,k) = max(liqwcin, 0.0) * dp_clp(i, k) / grav * 1000.
      end do
   end do

!! relative humidity (0-100) - from GEM "calcdiag.F90"
   call mfohr4(rhum, hum, temp, p2d, gni, gnk, gni, .false.)
   rhum = min(max(rhum * 100., 0.), 100.)

   do i = 1, gni
      ii = nmod(i)

!! sea-land fraction (1-0)
      slf(i) = metvar2d(ii, MV2D_MG)

!! cumulative albedo from 4 landuse types
      albedo(i) = metvar2d(ii, MV2D_AL5)

      dlat(i) = metvar2d(ii, MV2D_DLAT)
      dlon(i) = metvar2d(ii, MV2D_DLON)
   end do

!! cosine solar-zenith-angle
! recalculate cos(sza) since value from GEM is limited to 0-90deg.
! following GEM's ccc1_cccmarad, ccc2_cccmarad, newrad but with modified "mach_suncos"

   njul_day = real(ijul_day)
   read(chm_dttim_s(10:11), '(I2)') ihour
   read(chm_dttim_s(12:13), '(I2)') iminute
   hrnow = real(ihour) + (real(iminute) / 60.0)
   call mach_suncos(csza, gni, dlat, dlon, hrnow, njul_day)

   do i = 1, gni
!! get pbl layers in slab for simplified aerosol loading in MESSY
!! evaluate layer index of pbl from the ground up (surf. layer=gnk)
      do k = gnk, 1, -1
         pblidx(i) = k
         if (mh(i,k) >= metvar2d(nmod(i), MV2D_H)) exit
      end do
!
!! get layer thickness
      do k = 1, gnk-1
         dz2d(i, k) = mh(i, k) - mh(i, k+1)
      end do
      dz2d(i, gnk) = mh(i, gnk)
   end do

   !Indirect addressing for sunlit grid points based on czsa
   gni_day = 0
   do i = 1, gni
      if (csza(i) >= u0lim) then
         gni_day = gni_day + 1
         iu0(gni_day) = i
      end if
   end do

!! Call MESSY output jval(:,:,:) in sec-1
   if (gni_day > 0) then
      call mach_messy_jvalues(colo3, csza, p2d, rhum, temp, albedo, aclc, slf, &
                              clp, dz2d, midatmo, pblidx, iu0, gni_day, jval,  &
                              gni, gnk)
!
!! Transfer to output 1d bus
      if (chm_pkg_gas_s(1:5) == 'ADOM2') then
! note from Craig S.
! JNO2*0.005 is a reasonable estimate for DIAL
! (dicarbonyls - mostly for the aromatic ring opening products with a double bond in the middle of the molecule)
! Two photolysis rates that span the range possible for DIAL species, derived from TUV,
!  CH2=C(CH3)CHO -> Products               4.960E-06
!  CH3COCHO -> CH3CO + HCO                 1.344E-04
! Here is the corresponding J(NO2)
!   NO2 -> NO + O(3P)                      9.063E-03
! ==> JNO2*0.005 is a reasonable estimate falling in-between these two extreme.
         ik0 = 0
         do k = 1, gnk
            do i = 1, gni
               ik0 = ik0 + 1
               rj(ik0, 1 ) = jval(jval_NO2,      i, k) !NO2 + hv =  NO + O
               rj(ik0, 2 ) = jval(jval_NOO2,     i, k) !NO3 + hv =  NO
               rj(ik0, 3 ) = jval(jval_NO2O,     i, k) !NO3 + hv =  NO2 + O
               rj(ik0, 4 ) = jval(jval_O3P,      i, k) !O3 + hv = O
               rj(ik0, 5 ) = jval(jval_O1D,      i, k) !O3 + hv = O1D
               rj(ik0, 6 ) = jval(jval_HONO,     i, k) !HONO + hv =  OH + NO
               rj(ik0, 7 ) = jval(jval_HNO3,     i, k) !HNO3 + hv  =  NO2 + OH
               rj(ik0, 8 ) = jval(jval_HNO4,     i, k) !HNO4 + hv = NO2 + HO2
               rj(ik0, 9 ) = jval(jval_H2O2,     i, k) !H2O2 + hv = 2 OH
               rj(ik0, 10) = jval(jval_CH3OOH,   i, k) !ROOH + hv =  HO2 + OH
               rj(ik0, 11) = jval(jval_CHOH,     i, k) !HCHO + hv {+ 2 O2} = 2 HO2 + CO
               rj(ik0, 12) = jval(jval_COH2,     i, k) !HCHO +  hv = CO
               rj(ik0, 13) = jval(jval_CH3CHO,   i, k) !ALD2 + hv {+ 2 O2} = HCHO + RO2 + RO2R  + CO + HO2
               rj(ik0, 14) = jval(jval_CH3COCH3, i, k) !MEK + hv = ALD2 + MCO3 + RO2R + RO2
               rj(ik0, 15) = jval(jval_MGLYOX,   i, k) !MGLY + hv = MCO3 + CO + HO2
               rj(ik0, 16) = jval(jval_NO2,      i, k) * 5.0d-3 !DIAL + hv = HO2 + CO + MCO3
            end do
         end do
!
         rj = rj * 60.0d0  ! covert to from sec-1 to min-1 for ADOM
!
      else if (chm_pkg_gas_s(1:7) == 'SAPRC07') then  ! both SAPRC07C and SAPRC07CS
         ik0 = 0
         do k = 1, gnk
            do i = 1, gni
               ik0 = ik0 + 1
               rj(ik0, 1 ) = jval( jval_NO2     ,i,k) !NO2 + hv = NO + O3P
               rj(ik0, 2 ) = jval( jval_NOO2    ,i,k) !NO3 + hv = NO +O2
               rj(ik0, 3 ) = jval( jval_NO2O    ,i,k) !NO3 + hv = NO2 + O3P
               rj(ik0, 4 ) = jval( jval_O1D     ,i,k) !O3 + hv = O1D + O2
               rj(ik0, 5 ) = jval( jval_O3P     ,i,k) !O3 + hv = O3P + O2
               rj(ik0, 6 ) = jval( jval_HONO    ,i,k) !HONO + hv = OH + NO
               rj(ik0, 7 ) = jval( jval_HNO3    ,i,k) !HNO3 + hv = OH + NO2
               rj(ik0, 8 ) = jval( jval_HNO4    ,i,k) !HNO4 + hv = 0.61HO2 + 0.61NO2 + 0.39OH + 0.39NO3
               rj(ik0, 9 ) = jval( jval_H2O2    ,i,k) !H2O2 + hv = 2OH
               rj(ik0, 10) = jval( jval_PAN     ,i,k) !PAN + hv = 0.6MECO3 + 0.6NO2 + 0.4RO2R + 0.4HCHO + 0.4CO2 + 0.4NO3
               rj(ik0, 11) = jval( jval_PAN     ,i,k) !PAN2 + hv = 0.6RCO3 + 0.6NO2 + 0.4RO2R + 0.4CCHO + 0.4CO2 + 0.4NO3
               rj(ik0, 12) = jval( jval_CHOH    ,i,k) !HCHO + hv = 2HO2 + CO
               rj(ik0, 13) = jval( jval_COH2    ,i,k) !HCHO + hv = CO + H2
               rj(ik0, 14) = jval( jval_CH3CHO  ,i,k) !CCHO + hv = HO2 + RO2R + CO + HCHO
               rj(ik0, 15) = jval( jval_CH3CHO  ,i,k) * 1.403E-3 / 4.16E-4 !RCHO + hv = CCHO + RO2R + CO + HO2 -- by Craig S
               rj(ik0, 16) = jval( jval_H2O2    ,i,k) * 3.94E-4 / 5.64E-4  !XOOH + hv = HO2 + OH -- by Craig S
               rj(ik0, 17) = jval( jval_MGLYOX  ,i,k) !MGLY + hv = HO2 + CO + MECO3
               rj(ik0, 18) = jval( jval_MGLYOX  ,i,k) * 0.387 / 1.56E-2 !AFG1 + hv = 1.023HO2+0.173RO2R+0.305MECO3+0.5RCO3+0.695CO+0.173HCHO+0.418MGLY+0.003AFG1+0.753XC -- by Craig S
               rj(ik0, 19) = jval( jval_MGLYOX  ,i,k) * 0.387 / 1.56E-2 !AFG2 + hv = PROD2 - XC -- by Craig S
               rj(ik0, 20) = jval( jval_CH3COCH3,i,k) * 1.97e-4 / 2.48e-5 !IPRD + hv = 0.655HO2+0.141RO2R+0.084OH+0.354MECO3+0.343RCO3+0.084R2O2+0.867CO+0.429HCHO+0.184CCHO+0.246PROD2+0.123XC -- by Craig
               rj(ik0, 21) = jval( jval_CH3COCH3,i,k) * 4.69e-6 / 2.48e-5 !PROD2 + hv = 0.913RO2R+0.4MECO3+0.6RCO3+0.087RO2N+0.677R2O2+0.303HCHO+0.163CCHO+0.78RCHO-0.091XC -- by Craig
               rj(ik0, 22) = jval( jval_NO      ,i,k) * 2.35E-4 / 0.723 !RNO3 + hv = 0.344HO2+0.554RO2R+NO2+0.102RO2N+0.167R2O2+0.135HCHO+0.444CCHO+0.137RCHO+0.528PROD2+0.786XC -- by Craig
               rj(ik0, 23) = jval( jval_H2O2    ,i,k) * 3.94E-4 / 5.64E-4 !IOOH + hv = DUM -- Jack Proxy need check!!
            end do
         end do
         if (chm_pkg_gas_s(1:9) == 'SAPRC07CS') then
            ik0 = 0
            do k = 1, gnk
               do i = 1, gni
                  ik0 = ik0 + 1
                  rj(ik0, 24) = jval( jval_NO      ,i,k) !NO + hv = O1D + N
                  rj(ik0, 25) = jval( jval_N2O     ,i,k) !N2O + hv = N2 + O1D
                  rj(ik0, 26) = jval( jval_N2O5    ,i,k) !N2O5 + hv = NO2 + NO3
                  rj(ik0, 27) = jval( jval_NO3NOO  ,i,k) !N2O5 + hv = NO + NO3 + O3P
                  rj(ik0, 28) = jval( jval_O2      ,i,k) !O2 + hv = O3P
                  rj(ik0, 29) = jval( jval_CH4     ,i,k) !CH4 + hv = 2HO2 + CO + H2 {+ 2H}
               end do
            end do
         end if
      end if
!
   else
      jval = 0.0d0
      rj   = 0.0d0
   end if
!
   call msg_toall(CHM_MSG_DEBUG, 'mach_gas_messy [END]')
   if (chm_timings_L) call timing_stop_omp(331)
   !-----------------------------------------------------------------
   return

end subroutine mach_gas_messy
