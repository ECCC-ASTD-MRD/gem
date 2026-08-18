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
! Fichier/File   : mach_pkg_diag_mod.ftn90
! Creation       : S. Menard, August. 2013
! Description    : Defines the indices initialisation s/r and the meta info of
!                  fields that are part of the dry and wet deposition package
!
!============================================================================

module mach_pkg_diag_mod
   use chm_utils_mod,        only: pre_increment
   use chm_nml_mod,          only: chm_pkg_pm_s
   use chm_species_idx_mod
   use chm_species_info_mod, only: species_master
   use mach_cam_utils_mod,   only: isize

   save

   contains
!============================================================================
! Name           : pkg_diag_wet_idxinit
!
! Description    : Initialize the indices of diagnostic fields.  Starts with
!                  the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer function pkg_diag_wet_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start

         integer(kind=4) start_index, nb_fields

         start_index = idx

! Note only these need to be incremented, other are already incremented in other module packages.
         sp_WSO3    = pre_increment(idx)
         sp_CAT     = pre_increment(idx)
         sp_HCO3    = pre_increment(idx)
         sp_HION    = pre_increment(idx)
         sp_WH2O    = pre_increment(idx)

         nb_fields = idx - start_index
         pkg_diag_wet_idxinit = nb_fields


      end function pkg_diag_wet_idxinit
!============================================================================
! Name           : pkg_diag_col_idxinit
!
! Description    : Initialize the indices of diagnostic fields.  Starts with
!                  the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer function pkg_diag_col_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start

         integer(kind=4) start_index, nb_fields

         start_index = idx

         sp_tend_O3  = pre_increment(idx)
         sp_diff_O3  = pre_increment(idx)
! 3D total columns
         sp_col_O3   = pre_increment(idx)
         sp_col_NO2  = pre_increment(idx)
         sp_col_NO   = pre_increment(idx)
         sp_col_HU   = pre_increment(idx)
! 2D tropo. columns
         sp_tcol_O3   = pre_increment(idx)
         sp_tcol_NO2  = pre_increment(idx)
         sp_tcol_NO   = pre_increment(idx)
         sp_tcol_PAN  = pre_increment(idx)
         sp_tcol_HNO3 = pre_increment(idx)
         sp_tcol_N2O5 = pre_increment(idx)
         sp_tcol_SO2  = pre_increment(idx)
         sp_tcol_SO4  = pre_increment(idx)
         sp_tcol_CO   = pre_increment(idx)
         sp_tcol_HCHO = pre_increment(idx)
         sp_tcol_ISOP = pre_increment(idx)
!
         nb_fields = idx - start_index
         pkg_diag_col_idxinit = nb_fields


      end function pkg_diag_col_idxinit
!
!============================================================================
! Name           : pkg_diag_pm_idxinit
!
! Description    : Initialize the indices of the aerosols' diagnostic fields.  Starts with
!                  the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer function pkg_diag_pm_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start

         integer(kind=4) start_index, nb_fields

         start_index = idx

         sp_NAF = pre_increment(idx)
         sp_NAC = pre_increment(idx)
         sp_NAT = pre_increment(idx)

         sp_NUMC = pre_increment(idx)
         idx     = isize + sp_NUMC - 1  ! sp_NUMC takes up 'isize' places

         nb_fields = idx - start_index
         pkg_diag_pm_idxinit = nb_fields

      end function pkg_diag_pm_idxinit
!
!============================================================================
! Name           : pkg_diag_col_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
!
      subroutine pkg_diag_col_metainit()
         implicit none
!
! 3D
         species_master(sp_tend_o3) % out_name    =  "TGO3"
         species_master(sp_tend_o3) % out_string  =  "VN=phytd_to3; ON=TGO3; VD=OZONE TENDENCY GAS-PHASE (ug /kg/ sec) ; VS=T;VB=v0"

         species_master(sp_diff_o3) % out_name    =  "TDO3"
         species_master(sp_diff_o3) % out_string  =  "VN=TDO3_out_std; ON=TDO3; VD=OZONE TENDENCY DIFFUSION (ug /kg/ sec) ; VS=T;VB=v0"

         species_master(sp_col_O3) % out_name    =  "YO3"
         species_master(sp_col_O3) % out_string  =  "VN=YO3_out_std; ON=YO3; VD=OZONE 3D TOTAL COLUMN (D.U.) ; VS=T;VB=v0"

         species_master(sp_col_NO2) % out_name    =  "YNO2"
         species_master(sp_col_NO2) % out_string  =  "VN=YNO2_out_std; ON=YNO2; VD=NO2 3D TOTAL COLUMN (E+15 molec. cm-2) ; VS=T;VB=v0"

         species_master(sp_col_NO ) % out_name    =  "YNO"
         species_master(sp_col_NO ) % out_string  =  "VN=YNO_out_std; ON=YNO; VD=NO 3D TOTAL COLUMN (E+15 molec. cm-2) ; VS=T;VB=v0"

         species_master(sp_col_HU ) % out_name    =  "YHU"
         species_master(sp_col_HU ) % out_string  =  "VN=YHU_out_std; ON=YHU; VD=H2O VAPOR 3D TOTAL COLUMN (E+15 molec. cm-2) ; VS=T;VB=v0"

         species_master(sp_tcol_O3)  % out_name   =  "ZO3"
         species_master(sp_tcol_O3)  % out_string =  "VN=ZO3_out_std; ON=ZO3; VD=TROPO. COLUMN O3 (molec /cm2) ; VS=A;VB=v0"

         species_master(sp_tcol_NO2) % out_name   =  "ZNO2"
         species_master(sp_tcol_NO2) % out_string =  "VN=ZNO2_out_std; ON=ZNO2; VD=TROPO. COLUMN NO2 (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_NO ) % out_name   =  "ZNO"
         species_master(sp_tcol_NO ) % out_string =  "VN=ZNO_out_std;  ON=ZNO;  VD=TROPO. COLUMN NO  (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_PAN) % out_name   =  "ZPAN"
         species_master(sp_tcol_PAN) % out_string =  "VN=ZPAN_out_std;  ON=ZPAN;  VD=TROPO. COLUMN PAN (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_HNO3) % out_name   =  "ZHN3"
         species_master(sp_tcol_HNO3) % out_string =  "VN=ZHN3_out_std;  ON=ZHN3;  VD=TROPO. COLUMN HNO3 (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_N2O5) % out_name   =  "ZN25"
         species_master(sp_tcol_N2O5) % out_string =  "VN=ZN25_out_std;  ON=ZN25;  VD=TROPO. COLUMN N2O5 (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_SO2) % out_name   =  "ZSO2"
         species_master(sp_tcol_SO2) % out_string =  "VN=ZSO2_out_std; ON=ZSO2; VD=TROPO. COLUMN SO2 (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_SO4) % out_name   =  "ZSO4"
         species_master(sp_tcol_SO4) % out_string =  "VN=ZSO4_out_std; ON=ZSO4; VD=TROPO. COLUMN SO4 (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_CO ) % out_name   =  "ZCO"
         species_master(sp_tcol_CO ) % out_string =  "VN=ZCO_out_std;  ON=ZCO;  VD=TROPO. COLUMN CO  (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_HCHO) % out_name   =  "ZHCH"
         species_master(sp_tcol_HCHO) % out_string =  "VN=ZHCH_out_std;  ON=ZHCH;  VD=TROPO. COLUMN HCHO (molec /cm2); VS=A;VB=v0"

         species_master(sp_tcol_ISOP) % out_name   =  "ZISO"
         species_master(sp_tcol_ISOP) % out_string =  "VN=ZISO_out_std;  ON=ZISO;  VD=TROPO. COLUMN ISOP (molec /cm2); VS=A;VB=v0"
! 2D end

      end subroutine pkg_diag_col_metainit
!
!============================================================================
! Name           : pkg_diag_acc_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
!
      subroutine pkg_diag_acc_metainit()
         implicit none
!
! 3D Accumulators (gas)
!
         species_master(sp_col_O3)  % per_name   =  "AYO3"
         species_master(sp_col_O3)  % per_string =  "VN=AYO3_a3d; ON=AYO3; VD=Accum. O3 total column D.U.; VS=T;VB=p0"

         species_master(sp_col_NO2) % per_name   =  "AYN2"
         species_master(sp_col_NO2) % per_string =  "VN=AYN2_a3d; ON=AYN2; VD=Accum. NO2 total column Peta molec. cm-2; VS=T;VB=p0"

         species_master(sp_col_NO)  % per_name   =  "AYNO"
         species_master(sp_col_NO)  % per_string =  "VN=AYNO_a3d; ON=AYNO; VD=Accum. NO total column Peta molec. cm-2; VS=T;VB=p0"

         species_master(sp_col_HU)  % per_name   =  "AYHU"
         species_master(sp_col_HU)  % per_string =  "VN=AYHU_a3d; ON=AYHU; VD=Accum. Water Vapor total column Peta molec. cm-2; VS=T;VB=p0"

         species_master(sp_tend_O3) % per_name   =  "AGO3"
         species_master(sp_tend_O3) % per_string =  "VN=AGO3_a3d; ON=AGO3; VD=Accum. O3 tendency gas-phase  ug /kg /sec; VS=T;VB=p0"

         species_master(sp_O3 ) % per_name   =  "AMO3"
         species_master(sp_O3 ) % per_string =  "VN=AMO3_a3d; ON=AMO3; VD=Accum. O3  mmr ug /kg; VS=T;VB=p0"

         species_master(sp_NO2) % per_name   =  "AMN2"
         species_master(sp_NO2) % per_string =  "VN=AMN2_a3d; ON=AMN2; VD=Accum. NO2 mmr ug /kg; VS=T;VB=p0"

         species_master(sp_NO ) % per_name   =  "AMNO"
         species_master(sp_NO ) % per_string =  "VN=AMNO_a3d; ON=AMNO; VD=Accum. NO  mmr ug /kg; VS=T;VB=p0"

         species_master(sp_PAN) % per_name   =  "AMPN"
         species_master(sp_PAN) % per_string =  "VN=AMPN_a3d; ON=AMPN; VD=Accum. PAN mmr ug /kg; VS=T;VB=p0"

         species_master(sp_HNO3) % per_name   =  "AMH3"
         species_master(sp_HNO3) % per_string =  "VN=AMH3_a3d; ON=AMH3; VD=Accum. HNO3 mmr ug /kg; VS=T;VB=p0"

         species_master(sp_N2O5) % per_name   =  "AM25"
         species_master(sp_N2O5) % per_string =  "VN=AM25_a3d; ON=AM25; VD=Accum. N2O5 mmr ug /kg; VS=T;VB=p0"

         species_master(sp_SO2) % per_name   =  "AMS2"
         species_master(sp_SO2) % per_string =  "VN=AMS2_a3d; ON=AMS2; VD=Accum. SO2 mmr ug /kg; VS=T;VB=p0"

         species_master(sp_SO4) % per_name   =  "AMS4"
         species_master(sp_SO4) % per_string =  "VN=AMS4_a3d; ON=AMS4; VD=Accum. SO4 mmr ug /kg; VS=T;VB=p0"

         species_master(sp_CO ) % per_name   =  "AMCO"
         species_master(sp_CO ) % per_string =  "VN=AMCO_a3d; ON=AMCO; VD=Accum. CO  mmr ug /kg; VS=T;VB=p0"

         species_master(sp_HCHO) % per_name   =  "AMHC"
         species_master(sp_HCHO) % per_string =  "VN=AMHC_a3d; ON=AMHC; VD=Accum. HCHO mmr ug /kg; VS=T;VB=p0"

         species_master(sp_ISOP) % per_name   =  "AMIS"
         species_master(sp_ISOP) % per_string =  "VN=AMIS_a3d; ON=AMIS; VD=Accum. ISOP mmr ug /kg; VS=T;VB=p0"
!
! 2D Accumulators (12 gas-phase)
!
         species_master(sp_tcol_O3 ) % per_name   =  "AZO3"
         species_master(sp_tcol_O3 ) % per_string =  "VN=AZO3_a2d; ON=AZO3; VD=Accum.  O3 tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_NO2) % per_name   =  "AZN2"
         species_master(sp_tcol_NO2) % per_string =  "VN=AZN2_a2d; ON=AZN2; VD=Accum. NO2 tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_NO ) % per_name   =  "AZNO"
         species_master(sp_tcol_NO ) % per_string =  "VN=AZNO_a2d; ON=AZNO; VD=Accum. NO  tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_PAN) % per_name   =  "AZPN"
         species_master(sp_tcol_PAN) % per_string =  "VN=AZPN_a2d; ON=AZPN; VD=Accum. PAN tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_HNO3) % per_name   =  "AZH3"
         species_master(sp_tcol_HNO3) % per_string =  "VN=AZH3_a2d; ON=AZH3; VD=Accum. HNO3 tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_N2O5) % per_name   =  "AZ25"
         species_master(sp_tcol_N2O5) % per_string =  "VN=AZ25_a2d; ON=AZ25; VD=Accum. N2O5 tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_SO2) % per_name   =  "AZS2"
         species_master(sp_tcol_SO2) % per_string =  "VN=AZS2_a2d; ON=AZS2; VD=Accum. SO2 tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_SO4) % per_name   =  "AZS4"
         species_master(sp_tcol_SO4) % per_string =  "VN=AZS4_a2d; ON=AZS4; VD=Accum. SO4 tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_CO ) % per_name   =  "AZCO"
         species_master(sp_tcol_CO ) % per_string =  "VN=AZCO_a2d; ON=AZCO; VD=Accum. CO tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_HCHO) % per_name   =  "AZHC"
         species_master(sp_tcol_HCHO) % per_string =  "VN=AZHC_a2d; ON=AZHC; VD=Accum. HCHO tropo. column molec /cm2; VS=A;VB=p0"

         species_master(sp_tcol_ISOP) % per_name   =  "AZIS"
         species_master(sp_tcol_ISOP) % per_string =  "VN=AZIS_a2d; ON=AZIS; VD=Accum. ISOP tropo. column molec /cm2; VS=A;VB=p0"

      end subroutine pkg_diag_acc_metainit
!
!============================================================================
! Name           : pkg_diag_dry_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
!
      subroutine pkg_diag_dry_metainit()
         use mach_cam_utils_mod, only: isize, icom, aeroname, aero_sname
         use mach_pkg_gas_mod,   only: NOy_species_list
         use mach_drydep_mod,    only: gas_depo, nb_gas_depo

         use phy_options,        only: phyoutlist_S
         use clib_itf_mod,       only: clib_tolower
         implicit none

         character(len=4) :: sname, sname_out, lname
         integer(kind=4)  :: nsp, busid, istat
         logical(kind=4)  :: noy_dep, noy_specie
!
         ! Deposition flux of NOy (if requested)
         noy_dep    = .false.
         noy_specie = .false.
         if (any(phyoutlist_S == 'dnoy')) then
            species_master(sp_NOy) % dd_name   = "DNOY"
            species_master(sp_NOy) % dd_string = "VN=NOY_dd; ON=DNOY; VD=Dry deposition of NOy; VS=A;VB=p0"
            noy_dep = .true.
         end if
!
! Dry deposition (gas)
!
         do nsp = 1, nb_gas_depo
            busid = gas_depo(nsp) % sp_id
            lname = species_master(busid) % dyn_name

            ! Deposition flux
            sname = "D" // trim(gas_depo(nsp) % sname)
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (noy_dep) noy_specie = any(NOy_species_list == busid)
            if (any(phyoutlist_S == sname_out) .or. noy_specie) then
               species_master(busid) % dd_name   = sname
               species_master(busid) % dd_string = "VN="//trim(lname)//"_dd; ON="//sname// &
                                                 & "; VD=Dry deposition of gaseous "//     &
                                                 & trim(lname)//"; VS=A;VB=p0"
            end if

            ! Quasi-laminar sublayer resistance (Wesely's Rb term)
            sname = trim(gas_depo(nsp) % descr)//"RB"
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (any(phyoutlist_S == sname_out)) then
               species_master(busid) % rb_name   = sname
               species_master(busid) % rb_string = "VN="//trim(lname)//"_rb; ON="//sname//"; VD="// &
                                           & trim(lname)//" molecular diffusion resistance (s/m); VS=A;VB=V0"
            end if

            ! Total surface resistance
            sname = trim(gas_depo(nsp) % descr)//"RC"
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (any(phyoutlist_S == sname_out)) then
               species_master(busid) % rc_name   = sname
               species_master(busid) % rc_string = "VN="//trim(lname)//"_rc; ON="//sname//"; VD="// &
                                           & trim(sname)//" total surface resistance (s/m); VS=A*15;VB=V0"
            end if

            ! Deposition velocity for ground surface pathway
            sname = trim(gas_depo(nsp) % descr)//"VG"
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (any(phyoutlist_S == sname_out)) then
               species_master(busid) % vdg_name   = sname
               species_master(busid) % vdg_string = "VN="//trim(lname)//"_vdg; ON="//sname//"; VD="// &
                                           & trim(sname)//" deposition velocity for ground surface pathway (m/s); VS=A*15;VB=V0"
            end if
            
         end do
!
! Aerosol species dry deposition
         do nsp = 1, icom
            busid = nsp * isize + sp_AERO - 1  ! Use each aerosol's last bin index
            sname = 'D' // aero_sname(nsp)
            if (aero_sname(nsp) == 'PC') sname = 'DPOM'
            if (aero_sname(nsp) == 'OC') sname = 'DSOM'
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (any(phyoutlist_S == sname_out) .and. busid > 0) then
               species_master(busid) % dd_name   = sname
               species_master(busid) % dd_string = "VN="//trim(sname)//"_dd; ON="//sname// &
                                                 & "; VD=Dry deposition of (P) "// &
                                                 & aeroname(nsp) //"(moles/m2); VS=A;VB=p0"
            end if
         end do

         return
      end subroutine pkg_diag_dry_metainit
!
!============================================================================
! Name           : pkg_diag_wet_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
!
      subroutine pkg_diag_wet_metainit()
         use mach_cam_utils_mod, only: icom, aeroname, aero_sname
         use phy_options,        only: phyoutlist_S
         use clib_itf_mod,       only: clib_tolower
         implicit none

         integer(kind=4)  :: busid, iae, istat
         character(len=4) :: sname, sname_out

! wet deposition
         if (sp_H2O2 > 0) then
            species_master(sp_H2O2) % wd_name   =  "WH22"
            species_master(sp_H2O2) % wd_string =  "VN=WH22_wd; ON=WH22; VD=Wet deposition of H2O2 (hydrogen peroxide); VS=A;VB=p0"
         end if

         if (sp_ROOH > 0) then
            species_master(sp_ROOH) % wd_name   =  "WROO"
            species_master(sp_ROOH) % wd_string =  "VN=WROO_wd; ON=WROO; VD=Wet deposition of ROOH (organic peroxide); VS=A;VB=p0"
         end if

         if (sp_XOOH > 0) then
            species_master(sp_XOOH) % wd_name   =  "WROO"
            species_master(sp_XOOH) % wd_string =  "VN=WROO_wd; ON=WROO; VD=Wet deposition of XOOH (organic peroxide); VS=A;VB=p0"
         end if

         species_master(sp_WSO3) % wd_name   =  "WSO3"
         species_master(sp_WSO3) % wd_string =  "VN=WSO3_wd; ON=WSO3; VD=Wet deposition of HSO3- (moles/m2); VS=A;VB=p0"

         species_master(sp_CAT) % wd_name   =  "CAT"
         species_master(sp_CAT) % wd_string =  "VN=CAT_wd; ON=CAT; VD=Wet deposition of base cations (moles/m2); VS=A;VB=p0"

         species_master(sp_HCO3) % wd_name   =  "HCO3"
         species_master(sp_HCO3) % wd_string =  "VN=HCO3_wd; ON=HCO3; VD=Wet deposition of HCO3(-) (moles/m2); VS=A;VB=p0"

         species_master(sp_HION) % wd_name   =  "HION"
         species_master(sp_HION) % wd_string =  "VN=HION_wd; ON=HION; VD=Wet deposition of H(+) (moles/m2); VS=A;VB=p0"

         species_master(sp_WH2O) % wd_name   =  "WH2O"
         species_master(sp_WH2O) % wd_string =  "VN=WH2O_wd; ON=WH2O; VD=Wet deposition of aerosol water; VS=A;VB=p0"

         ! Aerosol components wet deposition
         do iae = 1, icom
            busid = iae * isize + sp_AERO - 1  ! Use the species last bin's index
            sname = 'W' // aero_sname(iae)
            if (aero_sname(iae) == 'SU') sname = 'WSO4'
            if (aero_sname(iae) == 'NI') sname = 'WNO3'
            if (aero_sname(iae) == 'AM') sname = 'WNH4'
            if (aero_sname(iae) == 'PC') sname = 'WPOM'
            if (aero_sname(iae) == 'OC') sname = 'WSOM'
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (any(phyoutlist_S == sname_out) .and. busid > 0) then
               species_master(busid) % wd_name   = sname
               species_master(busid) % wd_string = "VN="//trim(sname)//"_wd; ON="//sname// &
                                                 & "; VD=wet deposition of (P) "// &
                                                 & aeroname(iae) //"(moles/m2); VS=A;VB=p0"
            end if

         end do

         return
      end subroutine pkg_diag_wet_metainit

!
!============================================================================
! Name           : pkg_diag_pm_metainit
!
! Description    : Initialize the meta information for diagnostic PM fields
!
! Arguments:  None
!============================================================================
!
      subroutine pkg_diag_pm_metainit()
         use mach_cam_utils_mod, only: icom, aeroname, aero_sname, pm_agrege_nbin
         use phy_options,        only: phyoutlist_S
         use clib_itf_mod,       only: clib_tolower
         implicit none

         integer(kind=4)   :: nbin, nsp, iae, istat
         character(len=1)  :: zbin
         character(len=2)  :: cbin, &
                              ! Strings to describe speciated PM25 and PM10
                              pm_agrege_nm(pm_agrege_nbin), &
                              ! Strings to name speciated PM25 and PM10 output
                              pm_agrege_dx(pm_agrege_nbin)
         character(len=4)  :: sname, sname_out

         species_master(sp_NAF) % out_name   = "NAFP"
         species_master(sp_NAF) % out_string = "VN=NAF_vol_std; ON=NAFP; VD=PM2.5 number conc. (/m3) ; VS=T;VB=V0"

         species_master(sp_NAC) % out_name   = "NACP"
         species_master(sp_NAC) % out_string = "VN=NAC_vol_std; ON=NACP; VD=PM10 number conc. (/m3) ; VS=T;VB=V0"

         species_master(sp_NAT) % out_name   = "NATP"
         species_master(sp_NAT) % out_string = "VN=NAT_vol_std; ON=NATP; VD=PM total number conc. (/m3) ; VS=T;VB=V0"

         do nbin = 1, isize
            write(zbin, '(Z1)') nbin
            sname = "NUM" // zbin
            sname_out = sname
            istat = clib_tolower(sname_out)
            if (any(phyoutlist_S == sname_out)) then
               nsp = sp_NUMC + nbin - 1
               species_master(nsp) % out_name   = sname
               species_master(nsp) % out_string = "VN="//sname//"_out_std; ON="//sname// &
                                             &  "; VD= Aerosol number conc. in bin "//zbin//" (/m3) ; VS=T;VB=V0"
            end if
         end do

         if (isize <= 2) then
            pm_agrege_dx = (/'25', '10'/)
            pm_agrege_nm = (/'AF', 'AC'/)
         else
            pm_agrege_dx = (/'1 ', '25', '10', 'T '/)
            pm_agrege_nm = (/'AS', 'AF', 'AC', 'AT'/)
         end if
         do nbin = 1, pm_agrege_nbin

            cbin = pm_agrege_dx(nbin)
            do iae = 1, icom
               sname = aero_sname(iae)//pm_agrege_nm(nbin)
               sname_out = sname
               istat = clib_tolower(sname_out)
               if (any(phyoutlist_S == sname_out)) then
                  nsp = (iae - 1) * pm_agrege_nbin + nbin + sp_AERO - 1
                  species_master(nsp) % out_name   = sname
                  species_master(nsp) % out_string = "VN="//sname//"_out_std; ON="//sname// &
                                                    & "; VD="//aeroname(iae)//" PM"//  &
                                                    & cbin//" component (ug/m3);VS=T;VB=v0"
               end if
            end do
         end do

         return
      end subroutine pkg_diag_pm_metainit

end module mach_pkg_diag_mod
