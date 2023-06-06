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
! Fichier/File   : mach_pkg_saprc_mod.ftn90
! Creation       : Jack C. 2017
! Description    : Defines the indices initialisation s/r and the meta info of
!                  fields that are part of the SAPRC07{C,CS} gas package
!
! Apr 2017 - J. Chen - initialize species for gas-mechanism independent processes
!                    - called by "mach_main.ftn90"
!
!============================================================================

module mach_pkg_saprc_mod
   use chm_nml_mod,           only: chm_pkg_gas_s, chm_pkg_pm_s,        &
                                    chm_gas_drydep_s, chm_messy_jval_l
   use chm_utils_mod,         only: pre_increment
   use chm_species_idx_mod
   use chm_species_info_mod,  only: species_master
   use mach_drydep_mod,       only: depo, gas_depo, nb_gas_depo
   use mach_pkg_gas_mod,      only: gas_species, nvar, nfix, nspec, nreact,   &
                                    NOy_species_list, nb_NOy_species,         &
                                    num_be_sp, num_be_std, be_species,        &
                                    be_std_descr, be_std_name,                &
                                    njidx, njrxs, nonzero, lu_nonzero,        &
                                    ind_M, ind_N2, ind_O2, ind_H2O, ind_H2,   &
                                    ind_DUM, ind_CH4, ind_CO2, nsp_soa_gases, &
                                    soa_gas_idx, nvsoa, ind_O3, ind_NO,       &
                                    ind_NO3, ind_HO2
   use mach_pkg_messy_mod,    only: messy_njx, messy_jx_list

   implicit none

!! SAPRC mechanism generic and parameters
! NREACT - Number of reactions
   integer(kind=4), parameter :: nreact_saprc07 = 158

! for NOy deposition/output
   integer(kind=4), parameter :: nsp_NOy_saprc07 = 11
   integer(kind=4), target    :: NOy_species_list_saprc07(nsp_NOy_saprc07)

!! Photolysis reactions rates
! njrxs is number of photolysis reactions required (SAPRC07CS=21)
   integer(kind=4), parameter :: saprc07cs_njrxs = 29
   integer(kind=4), parameter :: saprc07c_njrxs = 23
! njidx is the photolysis reaction index associated with the host gas mechanism
!  index order needs to match that defined in "mach_gas_messy.ftn90"
   integer(kind=4), target    :: saprc07_njidx(saprc07cs_njrxs)

!! species list needed to undergo deposition
   integer(kind=4), parameter :: nsp_depo_saprc07 = 41
   type(depo), target         :: gas_depo_saprc07(nsp_depo_saprc07)

!! SAPRC07CS gas-chem specific variables and parameters

!! list of variable species that are advected
   integer(kind=4), parameter :: num_saprc07cs = 59
   integer(kind=4), parameter :: num_saprc07c  = 56
   integer(kind=4), target    :: saprc07cs_species(num_saprc07cs)
   integer(kind=4), target    :: saprc07c_species(num_saprc07c)

!! Biogenic standard emission rate
   integer(kind=4), parameter :: num_saprc07_be_std = 4
   integer(kind=4), parameter :: num_saprc07_be_sp = 5
   integer(kind=4), target    :: saprc07_be_species (num_saprc07_be_sp)
   character(len=4), target   :: saprc07_be_std_name(num_saprc07_be_std)
   character(len=2), target   :: saprc07_be_std_desc(num_saprc07_be_std)

!! list of soa species
   integer(kind=4), parameter :: num_saprc07_soa_gas = 8
   integer(kind=4), parameter :: num_saprc07_nvsoa = 12
   integer(kind=4), target    :: saprc07_soa_gas_idx(num_saprc07_soa_gas)

!! for MESSy phtolysis rate
! njx is the number of J-values to be calculated in MESSy
   integer(kind=4), parameter :: messy_njx_saprc07cs = 21
   integer(kind=4), parameter :: messy_njx_saprc07c = 15
   integer(kind=4), target    :: messy_jx_list_saprc07(messy_njx_saprc07cs)

   save

!============================================================================

  contains

!============================================================================
! Name           : pkg_saprc_idxinit
!
! Description    : Initialize the indices of the species used with various SAPRC mechanisms.
!                  Starts with the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer(kind=4) function pkg_saprc_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start
         integer(kind=4) start_index, nb_fields

!!! following to be done by one thread only

! Please keep alphabetical order
! (variable species)
         start_index = idx
         sp_AFG1   = pre_increment(idx)
         sp_AFG2   = pre_increment(idx)
         sp_ALK3   = pre_increment(idx)
         sp_ALK4   = pre_increment(idx)
         sp_ARO1   = pre_increment(idx)
         sp_ARO2   = pre_increment(idx)
         sp_BZO    = pre_increment(idx)
         sp_CCHO   = pre_increment(idx)
         sp_CH4    = pre_increment(idx)
         sp_CO     = pre_increment(idx)
         sp_CRES   = pre_increment(idx)
         sp_ETHE   = pre_increment(idx)
         sp_H2O2   = pre_increment(idx)
         sp_HCHO   = pre_increment(idx)
         sp_HNO3   = pre_increment(idx)
         sp_HNO4   = pre_increment(idx)
         sp_HO2    = pre_increment(idx)
         sp_HONO   = pre_increment(idx)
         sp_IEPX   = pre_increment(idx)
         sp_IOOH   = pre_increment(idx)
         sp_IPRD   = pre_increment(idx)
         sp_ISOP   = pre_increment(idx)
         sp_MCO3   = pre_increment(idx)
         sp_MGLY   = pre_increment(idx)
         if (chm_pkg_gas_s == 'SAPRC07CS') then
            sp_N      = pre_increment(idx)
            sp_N2O    = pre_increment(idx)
         end if
         sp_N2O5   = pre_increment(idx)
         sp_NO     = pre_increment(idx)
         sp_NO2    = pre_increment(idx)
         sp_NO3    = pre_increment(idx)
         sp_O1D    = pre_increment(idx)
         sp_O3     = pre_increment(idx)
         sp_O3P    = pre_increment(idx)
         sp_OH     = pre_increment(idx)  ! shldn't be transported, keep for now
         sp_OLE1   = pre_increment(idx)
         sp_OLE2   = pre_increment(idx)
         sp_PAN    = pre_increment(idx)
         sp_PAN2   = pre_increment(idx)
         sp_PRD2   = pre_increment(idx)
         sp_R2O2   = pre_increment(idx)
         sp_RCHO   = pre_increment(idx)
         sp_RCO3   = pre_increment(idx)
         sp_RN3a   = pre_increment(idx)
         sp_RO2N   = pre_increment(idx)
         sp_RO2R   = pre_increment(idx)
         sp_SO2    = pre_increment(idx)
         sp_SO4    = pre_increment(idx)
         sp_TERP   = pre_increment(idx)
         sp_XOOH   = pre_increment(idx)
         sp_yIOH   = pre_increment(idx)
! product only species
         sp_DUM    = pre_increment(idx) ! fixed
         sp_H2     = pre_increment(idx) ! fixed
         sp_XC     = pre_increment(idx) ! product only spc
         sp_XN     = pre_increment(idx) ! product only
! fixed species
         sp_CO2    = pre_increment(idx) ! fixed (product spc only)
         sp_H2O    = pre_increment(idx)  ! fixed
         sp_M      = pre_increment(idx)  ! fixed
         if (chm_pkg_gas_s == 'SAPRC07CS') then
            sp_N2  = pre_increment(idx)  ! fixed
         end if
         sp_O2     = pre_increment(idx)  ! fixed

         if (chm_pkg_pm_s /= 'NIL') then
            sp_NH3    = pre_increment(idx)
         end if

         NOy_species_list_saprc07(1:nsp_NOy_saprc07) = &
                     (/ sp_NO,   sp_NO2,  sp_HNO3, sp_HNO4, sp_N2O5, sp_NO3,  &
                        sp_PAN,  sp_PAN2, sp_RN3a, sp_RO2N, sp_HONO /)
         nb_NOy_species = nsp_NOy_saprc07
         NOy_species_list => NOy_species_list_saprc07

! Photolysis reactions
         saprc07_njidx = (/  1,  16,  17,  18,  19,  23,  28,  34,  41,  57, &
                            66,  85,  86,  90,  93,  96,  97, 104, 107, 111, &
                           113, 115, 145, 147, 148, 149, 150, 151, 152 /) ! entry 23-29 for SAPRC07CS
         njrxs = saprc07cs_njrxs   ! SAPRC07CS
         if (trim(chm_pkg_gas_s) == 'SAPRC07C') then
            njrxs = saprc07c_njrxs ! subset of saprc07cs_njrxs
         end if
         njidx => saprc07_njidx(1:njrxs)

! MESSy J-values to be calculated
         if (chm_messy_jval_l) then
            messy_jx_list_saprc07 = (/  2,  3,  4,  5,  6,  7,  9, 10, &
                                       11, 12, 14, 15, 17, 18, 19,  &
                                        1,  8, 56, 57, 66, 72  /) ! last 6 entry for SAPRC07CS
            messy_njx = messy_njx_saprc07cs

            if (trim(chm_pkg_gas_s) == 'SAPRC07C') then
               njrxs = saprc07c_njrxs ! subset of saprc07cs_njrxs
               messy_njx = messy_njx_saprc07c
            end if
            messy_jx_list => messy_jx_list_saprc07(1:messy_njx)
         end if

         if (chm_gas_drydep_s /= 'NIL') call define_deposited_species_saprc07()

! Set species for the gas-phase mechanism solver
         if (trim(chm_pkg_gas_s) == 'SAPRC07CS') then
            call set_kpp_parameters_saprc07cs()
!
         else if (trim(chm_pkg_gas_s) == 'SAPRC07C') then
            call set_kpp_parameters_saprc07c()
         end if

!  Set the biogenic emission field names and description
! IMPORTANT NOTE: 'NO' must be listed first for both the standard emission rate names
!                 and model be_species list. Similarly 'OVOC' needs to be listed lastly.
         saprc07_be_std_name = (/'NO  ', 'ISOP', 'MONO', 'OVOC'/)
         saprc07_be_std_desc = (/'NO', 'IO', 'MO', 'VO'/)
         saprc07_be_species = (/sp_NO, sp_ISOP, sp_TERP, sp_ALK3, sp_OLE2/)
         num_be_std = num_saprc07_be_std
         num_be_sp  = num_saprc07_be_sp
         be_std_descr => saprc07_be_std_desc
         be_std_name  => saprc07_be_std_name
         be_species   => saprc07_be_species

         nb_fields = idx - start_index
         pkg_saprc_idxinit = nb_fields

         return
      end function pkg_saprc_idxinit

!============================================================================
! defines selected KPP solver parameters for the set mechanism
!============================================================================
      subroutine set_kpp_parameters_saprc07cs
         implicit none

         nreact = 158
         nvar = 52
         nfix = 7
         nonzero = 604
         lu_nonzero = 653

         saprc07cs_species = &
            (/ sp_IEPX, sp_SO4 , sp_CO2,  sp_XN  , sp_XC  , sp_SO2 , &
               sp_IOOH, sp_yIOH, sp_H2O2, sp_HONO, sp_N2O , sp_XOOH, &
               sp_N   , sp_ARO1, sp_ARO2, sp_ALK4, sp_PAN , sp_PAN2, &
               sp_N2O5, sp_HNO4, sp_O1D , sp_BZO , sp_AFG2, sp_CRES, &
               sp_HNO3, sp_CO  , sp_AFG1, sp_ETHE, sp_MGLY, sp_R2O2, &
               sp_ISOP, sp_OLE1, sp_HCHO, sp_TERP, sp_RCHO, sp_OLE2, &
               sp_IPRD, sp_RN3a, sp_O3P , sp_CCHO, sp_O3  , sp_ALK3, &
               sp_RO2R, sp_HO2 , sp_RO2N, sp_MCO3, sp_OH  , sp_NO3 , &
               sp_PRD2, sp_NO  , sp_NO2 , sp_RCO3, sp_M   , sp_N2  , &
               sp_O2  , sp_H2O , sp_H2  , sp_CH4 , sp_DUM /)
         ind_CO2 = 3  ; ind_M   = 53 ; ind_N2  = 54 ; ind_O2  = 55
         ind_H2O = 56 ; ind_H2  = 57 ; ind_CH4 = 58 ; ind_DUM = 59
         nspec = num_saprc07cs
         gas_species => saprc07cs_species
!
! Set SOA scheme parameters
         nsp_soa_gases = num_saprc07_soa_gas
         nvsoa         = num_saprc07_nvsoa
         ! gas_species indices for sp_ALK4, sp_OLE1, sp_OLE2, sp_ARO1, &
         !                         sp_ARO2, sp_ISOP, sp_CRES, sp_TERP
         saprc07_soa_gas_idx = (/ 16, 32, 36, 14, 15, 31, 24, 34 /)
         soa_gas_idx => saprc07_soa_gas_idx
         ind_O3 = 41 ; ind_NO = 50 ; ind_NO3 = 48 ; ind_HO2 = 44

         return
      end subroutine set_kpp_parameters_saprc07cs

      subroutine set_kpp_parameters_saprc07c
         implicit none

         nreact = 146
         nvar = 50
         nfix = 6
         nonzero = 587
         lu_nonzero = 634

         saprc07c_species = &
            (/ sp_IEPX, sp_SO4 , sp_CO2,  sp_XN  , sp_XC  , sp_SO2 , &
               sp_O1D , sp_IOOH, sp_yIOH, sp_H2O2, sp_N2O5, sp_HONO, &
               sp_XOOH, sp_ARO1, sp_ARO2, sp_ALK4, sp_PAN , sp_PAN2, &
               sp_HNO4, sp_BZO , sp_AFG2, sp_CRES, sp_HNO3, sp_CO  , &
               sp_AFG1, sp_ETHE, sp_MGLY, sp_R2O2, sp_ISOP, sp_OLE1, &
               sp_HCHO, sp_TERP, sp_RCHO, sp_OLE2, sp_IPRD, sp_RN3a, &
               sp_O3P , sp_CCHO, sp_O3  , sp_ALK3, sp_RO2R, sp_HO2 , &
               sp_RO2N, sp_MCO3, sp_NO2 , sp_NO  , sp_NO3 , sp_PRD2, &
               sp_OH  , sp_RCO3, sp_M   , sp_O2  , sp_H2O , sp_H2  , &
               sp_CH4 , sp_DUM /)
         ind_CO2 = 3  ; ind_M = 51 ; ind_O2 = 52 ; ind_H2O = 53 ; ind_H2 = 54
         ind_CH4 = 55 ; ind_DUM = 56
         nspec = num_saprc07c
         gas_species => saprc07c_species
!
! Set SOA scheme parameters
         nsp_soa_gases = num_saprc07_soa_gas
         nvsoa         = num_saprc07_nvsoa
         ! gas_species indices for sp_ALK4, sp_OLE1, sp_OLE2, sp_ARO1, &
         !                         sp_ARO2, sp_ISOP, sp_CRES, sp_TERP
         saprc07_soa_gas_idx = (/ 16, 30, 34, 20, 21, 29, 22, 32 /)
         soa_gas_idx => saprc07_soa_gas_idx
         ind_O3 = 39 ; ind_NO = 46 ; ind_NO3 = 47 ; ind_HO2 = 42
!
         return
      end subroutine set_kpp_parameters_saprc07c

!============================================================================
! defines gas-phase mechanism (SAPRC07{C,CS}) specific deposition species and parameters
!============================================================================
      subroutine define_deposited_species_saprc07
         implicit none

!  Note: the order of the species for dry deposition is correlated to the order in
!  tables alpha, beta, hstar and fzero in mach_drydep_mod.ftn90. Do not change
!  the order without making sure that those tables have been modified as well.
         gas_depo_saprc07 % sp_id   = (/  &
                  sp_AFG1, sp_AFG2, sp_ALK3, sp_ALK4, sp_ARO1, sp_ARO2, &
                  sp_CCHO, sp_CRES, sp_ETHE, sp_HCHO, sp_N2O5, sp_HNO4, &
                  sp_HNO3, sp_HONO, sp_HO2,  sp_H2O2, sp_SO4,  sp_IPRD, &
                  sp_ISOP, sp_MCO3, sp_MGLY, sp_NO,   sp_NO2,  sp_OLE1, &
                  sp_OLE2, sp_O3,   sp_PAN,  sp_PAN2, sp_PRD2, sp_RCHO, &
                  sp_RCO3, sp_RN3a, sp_RO2N, sp_RO2R, sp_R2O2, sp_SO2,  &
                  sp_TERP, sp_XOOH, sp_IEPX, sp_yIOH, sp_NH3 /)

!  alpha and Beta are weights to be applied to resistance for SO2 and
!  Ozone for a given gas (from Table 1 of Zhang et al., 2002);
         gas_depo_saprc07 % alpha = (/ &
                   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
                   0.0,  0.0,  0.0,  0.8, 10.0,  5.0,  &
                  10.0,  2.0,  0.0,  1.0,  1.0,  0.0,  &
                   0.0,  0.0,  0.01, 0.0,  0.0,  0.0,  &
                   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
                   0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  &
                   0.0,  0.1,  1.0,  1.0,  1.0        /)

         gas_depo_saprc07 % beta  = (/ &
                  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  &
                  0.05,  0.05,  0.05,  0.2,  10.0,   5.0,   &
                 10.0,   2.0,   0.5,   1.0,   1.0,   0.05,  &
                  0.5,   0.5,   0.0,   0.01,  0.8,   0.05,  &
                  0.05,  1.0,   0.6,   0.6,   0.05,  0.05,  &
                  0.5,   0.5,   0.5,   0.5,   0.5,   0.0,   &
                  0.05,  0.8,   1.0,   1.0,   0.0          /)

   !  data for hstar (from Table 2 of Wesely, 1989): !  henry's constant: Solubility in water
   !  Those are used to evaluate how much gas is dissolved in mesophyll, cuticle, and other vegetation surfaces
         gas_depo_saprc07 % hstar =  (/ &
                  -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, &
                     15.0, -9999.0, -9999.0,  6000.0,  1.0E14,   1.0E7, &
                   1.0E14,   1.0E5, -9999.0,   1.0E5, -9999.0,    15.0, &
                  -9999.0, -9999.0, -9999.0,  3.0E-3,  1.0E-2, -9999.0, &
                  -9999.0,  1.0E-2,     3.6,     3.6, -9999.0,    15.0, &
                  -9999.0, -9999.0, -9999.0, -9999.0, -9999.0,   1.0E5, &
                  -9999.0, -9999.0,   1.0E5,   1.0E5,  2.0E04          /)

         gas_depo_saprc07 % fzero = (/  &
                  -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0,  &
                      0.0, -9999.0, -9999.0,     0.0,     0.0,     0.0,  &
                      0.0,     0.1, -9999.0,     1.0,     0.0,     0.0,  &
                  -9999.0, -9999.0, -9999.0,     0.0,     0.1, -9999.0,  &
                  -9999.0,     1.0,     0.1,     0.1, -9999.0,     0.0,  &
                  -9999.0, -9999.0, -9999.0, -9999.0, -9999.0,     0.0,  &
                  -9999.0, -9999.0,     1.0,     1.0,     0.0           /)

! Undefined parameters
         gas_depo_saprc07 % hstar2 = -9999.0
         gas_depo_saprc07 % fzero2 = 0.0
         gas_depo_saprc07 % ats    = -9999.0
         ! gas_depo_saprc07(nh3_index) % hstar2 = 61.0
         ! gas_depo_saprc07(nh3_index) % ats    = 3400.0

! Note: order for sname must correspond to that of sp_id
         gas_depo_saprc07 % sname =  &
                 (/ "AF1",  "AF2",  "AK3",  "AK4",  "AR1",  "AR2", &
                    "CCH",  "CRE",  "ETH",  "HCH",  "N25",  "HN4", &
                    "HN3",  "HNO",  "HO2",  "H22",  "SO4",  "IPR", &
                    "ISO",  "MC3",  "MGL",  "NO ",  "NO2",  "OL1", &
                    "OL2",  "O3 ",  "PAN",  "PN2",  "PRD",  "RCH", &
                    "RC3",  "RN3",  "R2N",  "R2R",  "R22",  "SO2", &
                    "TRP",  "XOO",  "IEP",  "YIH",  "NH3" /)

         gas_depo_saprc07 % descr = &
                 (/ "F1",  "F2",  "K3",  "K4",  "R1",  "R2", &
                    "CC",  "CR",  "ET",  "HC",  "N5",  "H4", &
                    "H3",  "HN",  "HO",  "H2",  "S4",  "IP", &
                    "IO",  "MC",  "MG",  "NO",  "N2",  "L1", &
                    "L2",  "O3",  "PA",  "P2",  "PR",  "RH", &
                    "RC",  "R3",  "RN",  "RR",  "22",  "S2", &
                    "TR",  "XO",  "IE",  "YI",  "NH" /)
!
         if (chm_pkg_pm_s /= 'NIL') then
            nb_gas_depo = nsp_depo_saprc07
         else
            nb_gas_depo = nsp_depo_saprc07 - 1  ! exclude NH3
         end if
         ! associate pointer and target
         gas_depo(1:nb_gas_depo) => gas_depo_saprc07(1:nb_gas_depo)

         return
      end subroutine define_deposited_species_saprc07

end module mach_pkg_saprc_mod
