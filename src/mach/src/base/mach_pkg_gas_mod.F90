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
! Fichier / File   : mach_pkg_gas_mod.ftn90
! Creation         : Deji Akingunola - May 2019
! Description      : To hold common/generalize parameters for the various
!                    gas-phase mechanisms
!
!============================================================================

module mach_pkg_gas_mod
   use chm_nml_mod,  only: chm_bkgd_ch4, chm_bkgd_co2
   implicit none

   ! Total number of gas species defined in  the model. Note that the actual number
   ! used at runtime depends on the selected gas-phase mechanisms
   integer(kind=4), parameter  :: nb_gas_species = 67
!
   ! Dynamically active (advected) gas tracers
   integer(kind=4), pointer, dimension(:) :: gas_species => null()

   integer(kind=4), pointer, dimension(:) :: NOy_species_list => null()
   integer(kind=4)                        :: nb_NOy_species = 0
   logical(kind=4) :: chm_noy_out_l   ! Is Atmospheric NOY requested for output

!  SOA yield calculations
   integer(kind=4), pointer, dimension(:) :: soa_gas_idx => null()
   integer(kind=4)                        :: nsp_soa_gases = 0
   integer(kind=4)                        :: nvsoa = 0

! MESSy photolysis reactions
   integer(kind=4), pointer, dimension(:) :: njidx => null()
   integer(kind=4)                        :: njrxs = 0

! Biogenic standard emission rate
   character(len=2), pointer, dimension(:) :: be_std_descr => null()
   character(len=4), pointer, dimension(:) :: be_std_name  => null()
   integer(kind=4), pointer, dimension(:)  :: be_species => null()
   integer(kind=4)                         :: num_be_std = 0
   integer(kind=4)                         :: num_be_sp = 0

   ! concentrations of fix species
   real(kind=8), parameter :: o2ppmv  = 2.1d05
   real(kind=8), parameter :: airppmv = 1.0d06
   real(kind=8), parameter :: n2ppmv  = 7.8d05
   real(kind=8), parameter :: minconc = tiny(1.0d0) ! minimum concentration
   ! minimum KPP concentration in molecules/cc
   ! (roughly equal to min. concentration of 1.0e-20ppmv used by tthe Y&B solver for ADOM2 mechanism)
   real(kind=8), parameter :: minkppconc = 2.0d-7
   real(kind=8), parameter :: minRx   = tiny(1.0d0) ! minimum reaction rate for MESSY
   real(kind=4), parameter :: o2      = 2.1e05
   real(kind=4), parameter :: aircon  = 1.0e06

   real(kind=4) :: ch4
   real(kind=8) :: ch4ppmv
   real(kind=8) :: co2ppmv

   ! Initialize the number of gas species for PPB output
   integer(kind=4)         :: nsp_gas_ppb_out = 0
   integer(kind=4), dimension(nb_gas_species) :: gas_ppb_out_list = 0

! Common KPP number fields
! NVAR - Number of Variable species
   integer (kind=4) :: nvar = 0
! NFIX - Number of Fixed species
   integer (kind=4) :: nfix = 0
!  Total number of species in the gas mechanism (variable + fixed)
   integer (kind=4) :: nspec = 0
! NREACT - Number of reactions
   integer (kind=4) :: nreact = 0
! NONZERO - Number of nonzero entries in Jacobian
   integer (kind=4) :: nonzero = 0
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
   integer (kind=4) :: lu_nonzero = 0

! Index declaration for fixed species
   integer(kind=4)  :: ind_c2h6 = 0
   integer(kind=4)  :: ind_ch4 = 0
   integer(kind=4)  :: ind_co2 = 0
   integer(kind=4)  :: ind_dum = 0
   integer(kind=4)  :: ind_h2 = 0
   integer(kind=4)  :: ind_h2o = 0
   integer(kind=4)  :: ind_m = 0
   integer(kind=4)  :: ind_n2 = 0
   integer(kind=4)  :: ind_o2 = 0
   integer(kind=4)  :: ind_oh = 0
   integer(kind=4)  :: ind_no = 0
! Additional index declaration for ODUM SOA species
   integer(kind=4)  :: ind_o3 = 0
   integer(kind=4)  :: ind_ho2 = 0
   integer(kind=4)  :: ind_no3 = 0

   save

   contains

!============================================================================
! Name           : pkg_gas_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
      subroutine pkg_gas_metainit()
         use chm_nml_mod,          only: chm_get_ae_emis_l, chm_get_be_emis_l, &
                                         chm_get_mj_emis_l, chm_htap_emis_l,   &
                                         chm_mono_s, chm_mass_s,               &
                                         chm_active_ch4_l, gas_ppb_out_list_s, &
                                         chm_ammonia_bidi_s
         use chm_species_idx_mod
         use chm_species_info_mod, only: species_master, unassigned
         use mach_drydep_mod,      only: gas_depo, nb_gas_depo
         implicit none
!
!        local variable
         character(len=22) :: tr_prop_s
         character(len=40), dimension(nb_gas_species) :: descr
         character(len=4),  dimension(nb_gas_species) :: dynname, outname, &
                                                         emisname, bioname
         integer(kind=4),   dimension(nb_gas_species) :: sp_id
         real(kind=4),      dimension(nb_gas_species) :: molwgt
         character(len=31) :: lname
         character(len=4)  :: sname

         integer(kind=4)   :: idx, isp, busid

#define CHM_SPECIE(indx, spx, DNAME, ONAME, ENAME, BNAME, DESC, MWGT) idx = idx + 1 ; sp_id(idx) = spx ; \
 dynname(idx) = DNAME ; outname(idx) = ONAME ; emisname(idx) = ENAME ; \
 bioname(idx) = BNAME ; descr(idx) = DESC ; molwgt(idx) = MWGT

         tr_prop_s = ';min=0'
         select case (chm_mono_s)
           case ('CLIP')
             tr_prop_s = trim(tr_prop_s)//';monot=1'
           case ('ILMC')
             tr_prop_s = trim(tr_prop_s)//';monot=2'
           case ('NIL')
             continue
         end select
         select case (chm_mass_s)
           case ('BC')
             tr_prop_s = trim(tr_prop_s)//';massc=1'
           case ('NIL')
         end select

         idx = 0
!! 'Lumped photoreactive monounsaturated dicarbonyl aromatic fragmentation products that photolyze to form radicals.' - Aromatic unsaturated ring
         CHM_SPECIE(idx, sp_AFG1, 'TAF1', 'AFG1', 'EAF1', '    ','Photoreactive monounsaturated dicarbonyl  ', 98.09994)
!! 'Lumped photoreactive monounsaturated dicarbonyl aromatic fragmentation products that photolyze to form non-radical products' - Aromatic unsaturated ring
         CHM_SPECIE(idx, sp_AFG2, 'TAF2', 'AFG2', '    ', '    ','Lumped monounsaturated dicarbonyl aromatic', 98.09994)
         CHM_SPECIE(idx, sp_ALD2, 'TALD', 'ALD2', 'EALD', '    ','Acetaldehyde and higher aldehydes         ', 44.05)
         CHM_SPECIE(idx, sp_ALK3, 'TAK3', 'ALK3', 'EA3 ', 'BPAA','Mid-sized alkanes                         ', 58.61)
         CHM_SPECIE(idx, sp_ALK4, 'TAK4', 'ALK4', 'EA4 ', '    ','Mid to Large alkanes                      ', 77.60)
         CHM_SPECIE(idx, sp_ALKA, 'TA3 ', 'ALKA', 'EA3 ', 'BPAA','C4+ alkanes                               ', 93.43)
         CHM_SPECIE(idx, sp_ALKE, 'TA2 ', 'ALKE', 'EA2 ', 'BPAE','C3+ alkenes                               ', 57.30)
!! 'Aromatics with kOH < 2x104 ppm-1 min-1.' - Emitted Lumped model VOC
         CHM_SPECIE(idx, sp_ARO1, 'TAR1', 'ARO1', 'EAR1', '    ','Moderate reacting aromatics               ', 95.16)
         CHM_SPECIE(idx, sp_ARO2, 'TAR2', 'ARO2', 'EAR2', '    ','Fast reacting aromatics (kOH > 2.E4ppm/min', 118.72)
         CHM_SPECIE(idx, sp_AROM, 'TARO', 'AROM', 'EARO', '    ','Lumped higher aromatics                   ', 117.97)
         CHM_SPECIE(idx, sp_BZO,  'TBZO', 'BZO ', '    ', '    ','Phenoxy radical                           ', 93.00)
         CHM_SPECIE(idx, sp_C3H8, 'TC38', 'C3H8', 'EC38', '    ','Propane and other slowly reacting organics', 44.09)
!! 'Acetaldehyde' - Explicit and Lumped Molecule Reactive Organic
         CHM_SPECIE(idx, sp_CCHO, 'TCHC', 'CCHO', 'ECHC', '    ','Acetaldehyde                              ', 44.05)
         CHM_SPECIE(idx, sp_CH4,  'TCH4', 'CH4G', 'ECH4', '    ','Methane                                   ', 16.04)
         CHM_SPECIE(idx, sp_CO,   'TCO ', 'CO_P', 'ECO ', '    ','Carbon monoxide                           ', 28.00)
         CHM_SPECIE(idx, sp_CRG1, 'TCR1', 'CRG1', '    ', '    ','Criegee radical #1                        ', 46.00)
         CHM_SPECIE(idx, sp_CRG2, 'TCR2', 'CRG2', '    ', '    ','Criegee radical #2                        ', 60.00)
         CHM_SPECIE(idx, sp_CRES, 'TCRE', 'CRES', 'ECRE', '    ','Lumped creosols                           ', 108.13)
         CHM_SPECIE(idx, sp_DIAL, 'TDIA', 'DIAL', '    ', '    ','General dicarbonyl                        ', 84.00)
         CHM_SPECIE(idx, sp_ETHE, 'TETH', 'ETHE', 'EETH', '    ', 'Ethene                                   ', 28.05)
         CHM_SPECIE(idx, sp_H2O2, 'TH22', 'H2O2', '    ', '    ','Hydrogen peroxide                         ', 34.00)
         CHM_SPECIE(idx, sp_HCHO, 'THCH', 'HCHO', 'EHCH', '    ','Formaldehyde                              ', 30.03)
         CHM_SPECIE(idx, sp_HNO3, 'THN3', 'HNO3', '    ', '    ','Nitric acid                               ', 63.00)
         CHM_SPECIE(idx, sp_HNO4, 'THN4', 'HNO4', '    ', '    ','Pernitric acid                            ', 79.01)
         CHM_SPECIE(idx, sp_HO2,  'THO2', 'HO2 ', '    ', '    ','Hydroperoxyl radical                      ', 33.00)
         CHM_SPECIE(idx, sp_HONO, 'THON', 'HONO', 'EHON', '    ','Nitrous acid                              ', 47.00)
         CHM_SPECIE(idx, sp_IEPX, 'TIPE', 'IEPX', '    ', '    ','Isoprene epoxydiols                       ', 118.0)
         CHM_SPECIE(idx, sp_IPRD, 'TIPR', 'IPRD', 'EIPR', '    ','Lumped isoprene product species           ', 100.12)
         CHM_SPECIE(idx, sp_IOOH, 'TIOH', 'IOOH', '    ', '    ','IOOH + OH = OH + IEPX                     ', 1.00)
         CHM_SPECIE(idx, sp_ISOP, 'TISO', 'ISOP', 'EISO', 'BPIO','Isoprene                                  ', 68.11)
!! Peroxy Propionyl and higher peroxy acyl Radicals - active_radicals
         CHM_SPECIE(idx, sp_MCO3, 'TMC3', 'MCO3', '    ', '    ','Acetyl Peroxy (CH3CO3) Radical            ', 75.00)
         CHM_SPECIE(idx, sp_MEK,  'TMEK', 'MEK ', 'EMEK', '    ','Methyl ethyl ketone and higher ketones    ', 72.10)
         CHM_SPECIE(idx, sp_MGL0, 'TMGL', 'MGLY', '    ', '    ','Methyl glyoxal (ADOM2)                    ', 72.00)
         CHM_SPECIE(idx, sp_MGLY, 'TMGL', 'MGLY', 'EMGL', '    ','Methyl glyoxal                            ', 72.00)
         CHM_SPECIE(idx, sp_N,    'TAN ', 'N_A ', '    ', '    ','Atomic nitrogen                           ', 14.00)
         CHM_SPECIE(idx, sp_N2O,  'TN2O', 'N2OG', '    ', '    ','Nitrous oxide                             ', 44.01)
         CHM_SPECIE(idx, sp_N2O5, 'TN25', 'N2O5', '    ', '    ','Dinitrogen pentoxide                      ', 108.00)
         CHM_SPECIE(idx, sp_NH3,  'TNH3', 'NH3 ', 'ENH3', '    ','Ammonia                                   ', 17.03)
         CHM_SPECIE(idx, sp_NO,   'TNO ', 'NO  ', 'ENO ', 'BPNO','Nitric oxide                              ', 30.00)
         CHM_SPECIE(idx, sp_NO2,  'TNO2', 'N2  ', 'ENO2', '    ','Nitrogen dioxide                          ', 46.00)
         CHM_SPECIE(idx, sp_NO3,  'TNO3', 'NO3 ', '    ', '    ','Nitrate radical                           ', 62.00)
         CHM_SPECIE(idx, sp_O1D,  'TO1D', 'O1D ', '    ', '    ','Excited-state oxygen atom                 ', 16.00)
         CHM_SPECIE(idx, sp_O3P,  'TO3P', 'O3P ', '    ', '    ','Ground-state oxygen atom                  ', 16.00)
         CHM_SPECIE(idx, sp_O,    'TO  ', 'O3P ', '    ', '    ','Ground-state oxygen atom                  ', 16.00)
         CHM_SPECIE(idx, sp_OSD,  'TOSD', 'O1D ', '    ', '    ','Excited-state oxygen atom                 ', 16.00)
         CHM_SPECIE(idx, sp_O3,   'TO3 ', 'O3  ', '    ', '    ','Ozone                                     ', 48.00)
         CHM_SPECIE(idx, sp_OH,   'TOH ', 'OH  ', '    ', '    ','Hydroxyl radical                          ', 17.00)
!! Alkenes (other than ethene) with kOH < 7x104 ppm-1 min-1. - Emitted Lumped model VOC
         CHM_SPECIE(idx, sp_OLE1, 'TOL1', 'OLE1', 'EOL1', '    ','Moderate reacting alkenes                 ', 72.34)
         CHM_SPECIE(idx, sp_OLE2, 'TOL2', 'OLE2', 'EOL2', 'BPAE','Fast reacting alkenes(kOH > 7.E4ppm/min   ', 75.78)
         CHM_SPECIE(idx, sp_PAN,  'TPAN', 'PAN ', '    ', '    ','Peroxyacetyl nitrate (PAN)                ', 121.00)
         CHM_SPECIE(idx, sp_PAN2, 'TPN2', 'PAN2', '    ', '    ','PPN and other higher alkyl PAN analogues  ', 1.00)
         CHM_SPECIE(idx, sp_PRD2, 'TPRD', 'PRD2', 'EPR2', '    ','Ketones and other oxygenated non-aldehyde ', 116.16)
         CHM_SPECIE(idx, sp_R2O2, 'TR22', 'R2O2', '    ', '    ','General organic peroxyl radical #2        ', 100.00)
!! Lumped C3+ Aldehydes - Explicit and Lumped Molecule Reactive Organic
         CHM_SPECIE(idx, sp_RCHO, 'TRCH', 'RCHO', 'ERCH', '    ','Propionaldehyde and larger aldehydes      ', 58.08)
         CHM_SPECIE(idx, sp_RCO3, 'TRCO', 'RCO3', '    ', '    ','Acyl Peroxy Radical (C3+)                 ', 89.10)
         CHM_SPECIE(idx, sp_RO2N, 'TR2N', 'RO2N', '    ', '    ','Alkyl nitrate => organic peroxyl radical  ', 100.00)
         CHM_SPECIE(idx, sp_RO2R, 'TR2R', 'RO2R', '    ', '    ','General organic peroxyl radical #1        ', 100.00)
         CHM_SPECIE(idx, sp_RNO3, 'TRN3', 'RNO3', '    ', '    ','Organic nitrate (RNO3)                    ', 121.00)
         CHM_SPECIE(idx, sp_RN3a, 'TRN3', 'RNO3', '    ', '    ','Organic nitrate (RNO3 - SO7C)             ', 147.18)
         CHM_SPECIE(idx, sp_RO2,  'TRO2', 'RO2 ', '    ', '    ','Total RO2 radicals                        ', 61.00)
         CHM_SPECIE(idx, sp_ROOH, 'TROO', 'ROOH', '    ', '    ','Lumped Organic peroxide                   ', 62.00)
         CHM_SPECIE(idx, sp_SO2,  'TSO2', 'S2  ', 'ESO2', '    ','Sulphur dioxide                           ', 64.00)
         CHM_SPECIE(idx, sp_SO4,  'TSO4', 'S4  ', 'ESO4', '    ','Sulphuric acid gas                        ', 96.00)
         CHM_SPECIE(idx, sp_TERP, 'TTRP', 'TERP', 'ETRP', 'BTRP','Monoterpenes - Emitted Lumped model VOCs  ', 136.24)
         CHM_SPECIE(idx, sp_TOLU, 'TTOL', 'TOLU', 'ETOL', '    ','Toluene and other mono-alkylbenzenes      ', 92.13)
!! 'Lumped organic hydroperoxide structure, used to represent the effect of photolysis of the hydroperoxide group in the SAPRC-99 peroxy radical representation (Carter, 2000a)' - Explicit and Lumped Molecule Reactive Organic
         CHM_SPECIE(idx, sp_XOOH, 'TROO', 'XOOH', '    ', '    ','Lumped Organic peroxide (S07C)            ', 1.00)
         CHM_SPECIE(idx, sp_yIOH, 'TIYH', 'yIOH', '    ', '    ','yIOH                                      ', 1.0)

         do isp = 1, nb_gas_species

            if ((sp_id(isp) == sp_CH4) .and. (.not. chm_active_ch4_l)) cycle

            if (sp_id(isp) > 0) then
               if (dynname(isp)(1:1) /= ' ') then
                  species_master(sp_id(isp)) % dyn_name   = dynname(isp)
                  species_master(sp_id(isp)) % dyn_string = "VN=TR/"//trim(dynname(isp))//":P"//trim(tr_prop_s)// &
                                                       & "; ON="//dynname(isp)//"; VD="//trim(descr(isp))//    &
                                                       & "; VS=T;VB=D1"

               end if

               species_master(sp_id(isp)) % mol_wt    = molwgt(isp)
!  For gaseous species output in ppb
               if (gas_ppb_out_list_s(1) == 'TOUT' .or. gas_ppb_out_list_s(1) == 'ALL ' .or. &
                   any(gas_ppb_out_list_s == outname(isp))) then
                  nsp_gas_ppb_out = nsp_gas_ppb_out + 1
                  species_master(sp_id(isp)) % out_name   = outname(isp)
                  species_master(sp_id(isp)) % out_string = "VN=sp_"//trim(outname(isp))//"; ON="// &
                                                       & outname(isp)//"; VD=PPB output of "//trim(descr(isp))//    &
                                                       & "; VS=T;VB=V0"
         ! Set the list (indices) of gas to output in PPB (from input namelist list)
                  gas_ppb_out_list(nsp_gas_ppb_out) = sp_id(isp)
               end if
!  Case where biogenic emissions option is set
!  Biogenic modulated rates (change every timestep)
               if (chm_get_be_emis_l .and. bioname(isp)(1:1) == 'B') then
                  species_master(sp_id(isp)) % be_name   = bioname(isp)
                  species_master(sp_id(isp)) % be_string = "VN=BIO_"//trim(outname(isp))//"; ON="// &
                                                       & bioname(isp)//"; VD=Biogenic emission of "//trim(outname(isp))// &
                                                       & " (g); VS=A;VB=P0"
               end if
!  Case where area emissions option is set
               if (chm_get_ae_emis_l .and. emisname(isp)(1:1) == 'E') then
                  species_master(sp_id(isp)) % ae_name   = emisname(isp)
                  species_master(sp_id(isp)) % ae_string = "VN=EMIS_"//trim(outname(isp))//"; ON="// &
                                                       & emisname(isp)//"; VD=Area emission of "//trim(outname(isp))// &
                                                       & " (g/s); VS=A;VB=P1"
               end if
!  Case where major point emissions option is set
               if (chm_get_mj_emis_l .and. (.not. chm_htap_emis_l) .and. &
                   emisname(isp)(1:1) == 'E') then
                  species_master(sp_id(isp)) % me_name   = emisname(isp)
               end if

            end if
         end do

         if (chm_htap_emis_l) then
! only SO2 from HTAP volcano, no other mjr_pt
            species_master(sp_SO2) % me_name   = 'ESO2'
! exclude so4 for HTAP simulations
            if (sp_SO4 > 0) then
               species_master(sp_SO4) % ae_name   = unassigned
               species_master(sp_SO4) % ae_name   = unassigned
               species_master(sp_SO4) % ae_string = unassigned
            end if
         end if

! ammonia bidirectional flux
         if (trim(chm_ammonia_bidi_s) /= 'OFF' .and. sp_NH3 > 0) then

            species_master(sp_NH3) % bd_name   = "NHBD"
            species_master(sp_NH3) % bd_string = "VN=BIDI_NH3; ON=NHBD; VD=NH3 bidirectional flux (g/s) ; VS=A ; VB=V0"

            if (trim(chm_ammonia_bidi_s) == 'GEP2D') then
               species_master(sp_NH3) % gep_name   = "NHGP"
               species_master(sp_NH3) % gep_string = "VN=EP_NH3; ON=NHGP; VD=NH3 ground emissions potential ; VS=A*15 ; VB=P0"
            end if

         end if

! Gas species deposition velocity
         do isp = 1, nb_gas_depo
            busid = gas_depo(isp) % sp_id
            sname = "V" // gas_depo(isp) % sname
            lname = trim(gas_depo(isp) % sname)
            species_master(busid) % vd_name   = sname
            species_master(busid) % vd_string = "VN="//trim(lname)//"_vd; ON="//sname// "; VD="// &
                                           & trim(lname)//" deposition velocity (m/s); VS=A;VB=V0"
         end do

         ! Set fixed background concentration for CH4 and CO2
         ch4     = chm_bkgd_ch4
         ch4ppmv = dble(chm_bkgd_ch4)
         co2ppmv = dble(chm_bkgd_co2)

         return
      end subroutine pkg_gas_metainit

end module mach_pkg_gas_mod
