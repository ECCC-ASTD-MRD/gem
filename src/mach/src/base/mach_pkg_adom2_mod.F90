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
! Fichier/File   : mach_pkg_adom2_mod.ftn90
! Creation       : H. Landry, Dec. 2007
! Description    : Defines the indices initialisation s/r and the meta info of
!                  fields that are part of the ADOM2 gas package
!
!
!============================================================================

module mach_pkg_adom2_mod
   use chm_nml_mod,          only: chm_gas_drydep_s, chm_pkg_gas_s, &
                                   chm_messy_jval_l
   use chm_utils_mod,        only: pre_increment
   use chm_species_idx_mod
   use chm_species_info_mod, only: species_master
   use mach_drydep_mod,      only: depo, gas_depo, nb_gas_depo
   use mach_pkg_gas_mod,     only: gas_species, nvar, nfix, nspec, nreact,  &
                                   NOy_species_list, nb_NOy_species,        &
                                   nsp_soa_gases, soa_gas_idx,              &
                                   be_species, num_be_sp, num_be_std,       &
                                   be_std_descr, be_std_name,               &
                                   njidx, njrxs, lu_nonzero, chm_noy_out_l, &
                                   ind_c2h6, ind_ch4, ind_dum, ind_h2o, &
                                   ind_m, ind_o2, ind_oh, ind_no
   use mach_pkg_messy_mod,   only: messy_njx, messy_jx_list

!! ADOM2 gas-chem specific variables and parameters
!! species list needed to undergo deposition
   integer(kind=4), parameter :: nsp_depo_adom2 = 28
   type(depo), target         :: gas_depo_adom2(nsp_depo_adom2)

!! list of variable species that are advected
   integer(kind=4), parameter :: num_adom2_species = 45
   integer(kind=4), target    :: adom2_species(num_adom2_species)
   integer(kind=4), parameter :: num_adom2kpp_species = 46
   integer(kind=4), target    :: adom2kpp_species(num_adom2kpp_species)

!! Biogenic standard emission rate
   integer(kind=4), parameter :: num_adom2_be_std = 4
   integer(kind=4), parameter :: num_adom2_be_sp  = 4
   integer(kind=4), target    :: adom2_be_species (num_adom2_be_sp)
   character(len=4), target   :: adom2_be_std_name(num_adom2_be_std)
   character(len=2), target   :: adom2_be_std_desc(num_adom2_be_std)

!! list of soa species
   integer(kind=4), parameter :: num_adom2_soa_gas = 5
   integer(kind=4), target    :: adom2_soa_gas_idx(num_adom2_soa_gas)

!! for mach_output
   integer(kind=4), parameter :: nsp_NOy_adom2 = 6
   integer(kind=4), target    :: NOy_species_list_adom2(nsp_NOy_adom2)

!! for MESSy phtolysis rate
! njrxs is number of photolysis reactions required (ADOM=16)
   integer(kind=4), parameter :: adom2_njrxs = 16
! njidx is the photolysis reaction index associated with the host gas mechanism
!  index order needs to match that defined in "mach_gas_messy.ftn90"
   integer(kind=4), target    :: adom2_njidx(adom2_njrxs)
! number of Jx to calculate in MESSy (ADOM2)
   integer(kind=4), parameter :: messy_njx_adom2 = 15
   integer(kind=4), target    :: messy_jx_list_adom2(messy_njx_adom2)


! Sundry ADOM2 mechanism parameters;
!  For vertical profile of ethane in ADOM2 solvers
   ! Height of the treshold values (sigma)
   real(kind=4), parameter :: layerboundary1 = 0.869, layerboundary2 = 0.924
   real(kind=4), parameter :: layerconcentration1 = 0.001, &
                              layerconcentration2 = 0.003, &
                              layerconcentration3 = 0.005

!  nreac is number of reactions in the ADOM2 mechanism, Y&B solver
   integer(kind=4), parameter :: nreac = 114
!  nreact_kpp is number of reactions in the ADOM2 mechanism, KPP solver
   integer(kind=4), parameter :: nreact_kpp = 112
!  nprcf is number of variable product coefficients (function of t, etc.)
!  used in the adom mechani
   integer(kind=4), parameter :: nprcf = 26
!  number of advected species - "nh3"
! (i.e. advected species which are chemically active)
   integer(kind=4), parameter :: nsi = 28

!  time step sizes for gas phase mechanism
   real(kind=4), parameter               :: hstart  = 1.0e-6
   real(kind=4), parameter, dimension(2) :: hsub  = (/ 0.2,  2.0 /)
   real(kind=4), parameter, dimension(2) :: hmax  = (/ 5.0,  15.0 /)
   real(kind=4), parameter               :: ymings = 1.0e-5
   real(kind=4), parameter               :: ymin2gs = 1.0e-20
! END adomYB solver
!
   save

   contains
!============================================================================
! Name           : pkg_adom2_idxinit
!
! Description    : Initialize the indices of ADOM2 package.  Starts with
!                  the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer(kind=4) function pkg_adom2_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start

         integer(kind=4) start_index, nb_fields


! Please keep alphabetical order

         start_index = idx

         sp_ALD2 = pre_increment(idx)
         sp_ALKA = pre_increment(idx)
         sp_ALKE = pre_increment(idx)
         sp_AROM = pre_increment(idx)
         sp_BZO  = pre_increment(idx)
         sp_C3H8 = pre_increment(idx)
         sp_CH4  = pre_increment(idx)
         sp_CO   = pre_increment(idx)
         sp_CRES = pre_increment(idx)
         sp_CRG1 = pre_increment(idx)
         sp_CRG2 = pre_increment(idx)
         sp_DIAL = pre_increment(idx)
         sp_ETHE = pre_increment(idx)
         sp_H2O2 = pre_increment(idx)
         sp_HCHO = pre_increment(idx)
         sp_HNO3 = pre_increment(idx)
         sp_HNO4 = pre_increment(idx)
         sp_HO2  = pre_increment(idx)
         sp_HONO = pre_increment(idx)
         sp_ISOP = pre_increment(idx)
         sp_MCO3 = pre_increment(idx)
         sp_MEK  = pre_increment(idx)
         sp_MGL0 = pre_increment(idx)
         sp_N2O5 = pre_increment(idx)
         sp_NH3  = pre_increment(idx)
         sp_NO   = pre_increment(idx)
         sp_NO2  = pre_increment(idx)
         sp_NO3  = pre_increment(idx)
         sp_O    = pre_increment(idx)
         sp_O3   = pre_increment(idx)
         sp_OH   = pre_increment(idx)
         sp_OSD  = pre_increment(idx)
         sp_PAN  = pre_increment(idx)
         sp_R2O2 = pre_increment(idx)
         sp_RNO3 = pre_increment(idx)
         sp_RO2  = pre_increment(idx)
         sp_RO2N = pre_increment(idx)
         sp_RO2R = pre_increment(idx)
         sp_ROOH = pre_increment(idx)
         sp_SO2  = pre_increment(idx)
         sp_SO4  = pre_increment(idx)
         sp_TOLU = pre_increment(idx)
!
! fixed species
         sp_C2H6 = pre_increment(idx)
         sp_H2O  = pre_increment(idx)
         sp_M    = pre_increment(idx)
         sp_O2   = pre_increment(idx)
         if (chm_pkg_gas_s == 'ADOM2KPP') then
            sp_DUM = pre_increment(idx)
         end if
!
         if (chm_noy_out_l) then
            NOy_species_list_adom2 = (/ sp_NO, sp_NO2, sp_HNO3, sp_PAN, &
                                        sp_HONO, sp_RNO3/)
            nb_NOy_species = nsp_NOy_adom2
            NOy_species_list => NOy_species_list_adom2
         end if
!
         ! Set the indices of the gas species in the solver
         ! NOTE: NH3 is not part of this list, as it is not active in the gas
         !       mechanism
! NVAR - Number of Variable species
         nvar = 40
! NREACT is number of reactions in the ADOM2 mechanism with KPP solver, plus the 'nprcf
         nreact = nreact_kpp + nprcf
         if (trim(chm_pkg_gas_s) == 'ADOM2KPP') then
!   These parameters need to be consistent in all adomKPP routines:
! NFIX - Number of Fixed species
            nfix = 6
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
            lu_nonzero = 435

            adom2kpp_species = &
                      (/ sp_SO4 , sp_OSD , sp_H2O2, sp_PAN , sp_C3H8, sp_HONO, &
                         sp_N2O5, sp_ALKA, sp_HNO4, sp_TOLU, sp_DIAL, sp_AROM, &
                         sp_BZO , sp_SO2 , sp_CO  , sp_CRES, sp_MGL0, sp_HNO3, &
                         sp_ROOH, sp_CRG1, sp_CRG2, sp_RNO3, sp_ALKE, sp_ISOP, &
                         sp_ETHE, sp_O   , sp_MEK , sp_O3  , sp_R2O2, sp_NO3 , &
                         sp_RO2R, sp_RO2 , sp_HO2 , sp_ALD2, sp_RO2N, sp_MCO3, &
                         sp_NO  , sp_HCHO, sp_OH  , sp_NO2 , sp_M   , sp_O2  , &
                         sp_H2O , sp_CH4 , sp_C2H6, sp_DUM /)
            ind_m     = 41 ; ind_o2  = 42 ; ind_h2o = 43 ; ind_ch4 = 44
            ind_c2h6  = 45 ; ind_dum = 46 ; ind_no  = 37 ; ind_oh  = 39
            if (chm_messy_jval_l) then
               adom2_njidx = (/ 1, 13, 14, 15, 16, 20, &
                               23, 30, 36, 46, 47, 48, &
                               53, 60, 62, 101 /)
            else
               adom2_njidx = (/ 1, 16, 13, 14, 15, 20, &
                               23, 30, 36, 46, 47, 48, &
                               53, 60, 62, 101 /)
            end if
            gas_species => adom2kpp_species
            nspec = num_adom2kpp_species
            ! gas_species indices for sp_alka, sp_arom, sp_tolu, sp_alke, sp_isop
            nsp_soa_gases = num_adom2_soa_gas
            adom2_soa_gas_idx = (/ 8, 12, 10, 23, 24 /)
            soa_gas_idx => adom2_soa_gas_idx
!
         else  !(chm_pkg_gas_s == 'ADOM2' .or. chm_pkg_gas_s == 'ADOM2KPPB')
            nfix = 5
            adom2_species = &
                    (/ sp_SO2,  sp_SO4,  sp_NO,   sp_NO2,  sp_O3,   sp_H2O2, &
                       sp_HNO3, sp_PAN,  sp_C3H8, sp_ALKA, sp_ETHE, sp_ALKE, &
                       sp_TOLU, sp_AROM, sp_HCHO, sp_ALD2, sp_MEK,  sp_MGL0, &
                       sp_DIAL, sp_ROOH, sp_CRES, sp_HONO, sp_RNO3, sp_ISOP, &
                       sp_HO2,  sp_RO2,  sp_MCO3, sp_CO,   sp_OSD,  sp_O   , &
                       sp_NO3,  sp_N2O5, sp_HNO4, sp_OH,   sp_RO2R, sp_R2O2, &
                       sp_RO2N, sp_BZO,  sp_CRG1, sp_CRG2, sp_CH4 , sp_C2H6, &
                       sp_H2O , sp_O2 ,  sp_M /)
            ind_ch4 = 41 ; ind_c2h6 = 42 ; ind_h2o = 43 ; ind_o2 = 44 ; ind_m = 45
            ind_no  = 3  ; ind_oh   = 34
            adom2_njidx = (/ 1, 13, 14, 15, 16, 20, &
                            23, 30, 37, 48, 49, 50, &
                            55, 62, 64, 103 /)
            ! gas_species indices for sp_alka, sp_arom, sp_tolu, sp_alke, sp_isop
            nsp_soa_gases = num_adom2_soa_gas
            adom2_soa_gas_idx = (/ 10, 14, 13, 12, 24 /)
            soa_gas_idx => adom2_soa_gas_idx
            gas_species => adom2_species
            nspec = num_adom2_species
         end if
!
         njrxs = adom2_njrxs
         njidx => adom2_njidx

!  select jval to calculate in MESSy:
         messy_jx_list_adom2(1:messy_njx_adom2) = &
                                         (/ 2,   3,   4,   5,   6, &
                                            7,   9,  10,  12,  13, &
                                           14,  15,  17,  18,  19 /)
         messy_jx_list => messy_jx_list_adom2
         messy_njx = messy_njx_adom2

         if (chm_gas_drydep_s /= 'NIL') call define_adom2_deposited_species()

!  Set the biogenic emission field names and description
         adom2_be_std_name = (/'NO  ', 'ISOP', 'OVOC', 'MONO'/)
         adom2_be_std_desc = (/'NO', 'IO', 'VO', 'MO'/)
         adom2_be_species  = (/sp_NO, sp_ISOP, sp_ALKA, sp_ALKE/)
         num_be_std = num_adom2_be_std
         num_be_sp  = num_adom2_be_sp
         be_std_descr => adom2_be_std_desc
         be_std_name  => adom2_be_std_name
         be_species => adom2_be_species

         nb_fields = idx - start_index
         pkg_adom2_idxinit = nb_fields

      end function pkg_adom2_idxinit


!============================================================================
! defines gas-phase mechanism (ADOM2) specific deposition species and parameters

      subroutine define_adom2_deposited_species()
         implicit none
!
!  Note: the order of the species for dry deposition is correlated to the order in
!  tables alpha, beta, hstar and fzero in mach_drydep_mod.ftn90. Do not change
!  the order without making sure that those tables have been modified as well.
         gas_depo_adom2 % sp_id = (/ &
             sp_SO2,  sp_SO4,  sp_NO,   sp_NO2,  sp_O3,   sp_H2O2, sp_HNO3, &
             sp_PAN,  sp_C3H8, sp_ALKA, sp_ETHE, sp_ALKE, sp_TOLU, sp_AROM, &
             sp_HCHO, sp_ALD2, sp_MEK,  sp_MGL0, sp_DIAL, sp_ROOH, sp_CRES, &
             sp_HONO, sp_RNO3, sp_ISOP, sp_HO2,  sp_RO2,  sp_MCO3, sp_NH3 /)

!  alpha and Beta are weights to be applied to resistance for SO2 and
!  Ozone for a given gas (from Table 1 of Zhang et al., 2002);
!  corrected NO beta value on 15 Aug. 2014
!  Corrected unoxygenated VOCs, HO2, RO2, MCO3 values in Nov. 2017,
!  based on Paul's request from May 2016.
         gas_depo_adom2 % alpha = (/ 1.0, 1.0, 0.0, 0.0 , 0.0, 1.0, 10.0,  &
                                     0.0, 0.0, 0.0, 0.0 , 0.0, 0.0,  0.0,  &
                                     0.8, 0.0, 0.0, 0.01, 0.0, 0.1,  0.01, &
                                     2.0, 0.0, 0.0, 0.0 , 0.0, 0.0,  1.0/)
         gas_depo_adom2 % beta  = (/ 0.0, 1.0,  0.01, 0.8, 1.0,  1.0, 10.0, &
                                     0.6, 0.0,  0.0,  0.0, 0.0,  0.0, 0.0,  &
                                     0.2, 0.05, 0.05, 0.0, 0.05, 0.8, 0.0,  &
                                     2.0, 0.0,  0.0,  0.0, 0.0,  0.0, 0.0/)
         !  hstar henry's constant: Solubility in water  (from Table 2 of Wesely, 1989)
         gas_depo_adom2 % hstar =     (/                               &
                1.0E05, -9999.0,  2.0E-3, 1.0E-02, 1.0E-02, 1.0E05,  1.0E14,  &
                   3.6, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, &
                6000.0,    15.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, &
                1.0E05, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0,  2.0E04/)
         ! fzero: Chemical reactivity
         gas_depo_adom2 % fzero = (/                                   &
                   0.0,     0.0,     0.0,     0.1,     1.0,     1.0,     0.0, &
                   0.1, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, &
                   0.0,     0.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0, &
                   0.1, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0,     0.0 /)

! as above but for use with alternative formulation Robichaud2

         gas_depo_adom2 % hstar2 = (/                                  &
                1.2,    -9999.0,  1.9E-3, 1.2E-02, 1.3E-02,  1.0E05,   2.6E6, &
                3.6,     1.5E-3, -9999.0,  4.7E-3, -9999.0,    0.15, -9999.0, &
             6000.0,       15.0,    20.0,     1.0,     1.0,    30.0, -9999.0, &
               49.0,    -9999.0,  1.3E-2, -9999.0, -9999.0, -9999.0,    61.0 /)

         gas_depo_adom2 % fzero2 = (/                                  &
                   0.0,     0.0,     0.0,     0.1,     1.0,     1.0,     1.0, &
                   0.1,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                   1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0, &
                   0.1,     0.0,     0.0,     1.0,     0.0,     0.1,     0.0 /)

! data for temperature dependency for solubility
! (Sander, 1999 (Handbook) and Seinfeld and Pandis, 2006

         gas_depo_adom2 % ats = (/                                            &
                3020.,     0.0,   1500.0,  2500.0,  2300.0,  6800.0,  8700.0, &
               5900.0,   2700.,  -9999.0,  1800.0, -9999.0,  4000.0, -9999.0, &
               6800.0,   5800.,    5000., -9999.0, -9999.0, -9999.0, -9999.0, &
               4800.0, -9999.0, -9999.0, -9999.0, -9999.0,  -9999.0,  3400.0 /)

         gas_depo_adom2 % sname = &
                 (/ "SO2", "SO4", "NO ", "NO2", "O3 ", "H22", "HN3", &
                    "PAN", "C38", "A2 ", "ETH", "A3 ", "TOL", "ARO", &
                    "HCH", "ALD", "MEK", "MGL", "DIA", "ROO", "CRE", &
                    "HNO", "RN3", "ISO", "HO2", "RO2", "MC3", "NH3" /)

         gas_depo_adom2 % descr = &
                 (/ "S2", "S4", "NO", "N2", "O3", "H2", "H3", &
                    "PA", "C3", "A2", "ET", "A3", "TL", "AR", &
                    "HC", "AD", "MK", "MG", "DI", "RH", "CR", &
                    "HN", "R3", "IO", "HO", "RO", "MC", "NH" /)

!
! associate the generic gas deposition pointer to the ADOM-II specific target
!
         nb_gas_depo = nsp_depo_adom2
         gas_depo => gas_depo_adom2 ! associate pointer and target

         return
      end subroutine define_adom2_deposited_species

end module mach_pkg_adom2_mod


