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

module mach_pkg_misc_mod
   use chm_nml_mod,          only: chm_get_be_emis_l, chm_aqhi_l, chm_soa_s,   &
                                   chm_pkg_pm_s, chm_get_on_emis_l, chm_vit_l, &
                                   chm_diag_aero_opt_l, aero_opt_wavel,        &
                                   chm_canopy_shading_l, chm_dep_lai2d_l,      &
                                   chm_sat_seasons_l, chm_strato_s, chm_biog_s,&
                                   chm_pkg_gas_s
   use chm_utils_mod,        only: pre_increment
   use chm_species_idx_mod,  only: sp_AQ10, sp_AQ25, sp_density, sp_FMET,    &
                                   ce_LAI, sp_LAI, sp_LU15, sp_KTN, sp_NOy,  &
                                   sp_AC, sp_AF, sp_NWOC,  sp_JNO2, sp_STSE, &
                                   sp_be_std, sp_CVKT, sp_MVKT, sp_TVKT,     &
                                   sp_OPTD, sp_SSCA, sp_ASYM, var_step,      &
                                   ce_FRT, ce_POP, sp_HC, sp_CLU, sp_FRT,    &
                                   sp_FRL, sp_CRL,   &
                                   lin_c2, lin_c4, lin_c5, lin_c6, lin_c7
   use chm_species_info_mod, only: species_master, nb_dyn_tracers, UNASSIGNED, &
                                   chem_ent_vars, ent_vars_num
   use mach_pkg_gas_mod,     only: chm_noy_out_l, num_be_std, be_std_descr, &
                                   be_std_name
   use mach_drydep_mod,      only: chm_lu15_out_l
   use mach_mie_data_mod,    only: mach_mie_optical_properties, nwl_aod
   use phy_options,          only: phyoutlist_S, llinoz

   save

   contains
!============================================================================
! Name           : pkg_misc_idxinit
!
! Description    : Initialize the indices of miscelaneous fields.  Starts with
!                  the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer function pkg_misc_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start

         integer(kind=4) start_index, nb_fields

         start_index = idx

         sp_LU15 = pre_increment(idx)

         if (chm_get_be_emis_l .or. chm_dep_lai2d_l) then
            sp_LAI   = pre_increment(idx)
         end if

         if (chm_biog_s(1:3) /= 'NIL') then
            sp_be_std  = pre_increment(idx)
            idx = sp_be_std + (num_be_std * 2) - 1
         end if

         if (chm_soa_s /= 'NIL') then
            sp_NWOC = pre_increment(idx)
         end if

         if (chm_pkg_pm_s(1:3) == 'CAM') then
            sp_AC = pre_increment(idx)
            sp_AF = pre_increment(idx)
         end if

         if (chm_aqhi_l) then
            sp_AQ25 = pre_increment(idx)
            sp_AQ10 = pre_increment(idx)
         end if

         if (chm_noy_out_l) then
            sp_NOy = pre_increment(idx)
         end if

         sp_KTN = pre_increment(idx)

         if (any(phyoutlist_S == 'rho')) then
            sp_density  = pre_increment(idx)
         end if

         if (any(phyoutlist_S == 'jno2')) then
            sp_JNO2  = pre_increment(idx)
         end if

         if (chm_pkg_pm_s /= 'NIL') then
            sp_FMET  = pre_increment(idx)
         endif

         if (chm_vit_l) then
            sp_CVKT = pre_increment(idx)
            sp_MVKT = pre_increment(idx)
            sp_TVKT = pre_increment(idx)
         end if

         if (chm_diag_aero_opt_l) then
            sp_OPTD = pre_increment(idx)
            idx = sp_OPTD + nwl_aod - 1
            sp_SSCA = pre_increment(idx)
            idx = sp_SSCA + nwl_aod - 1
            sp_ASYM = pre_increment(idx)
            idx = sp_ASYM + nwl_aod - 1
         end if

! Satellite LAI based seasonal mask
         if (chm_sat_seasons_l) then
            ce_LAI  = pre_increment(ent_vars_num)
            sp_STSE = pre_increment(idx)
         end if

         if (chm_canopy_shading_l) then
            ce_POP = pre_increment(ent_vars_num)
            ce_FRT = pre_increment(ent_vars_num)

            sp_HC  = pre_increment(idx)
            sp_CLU = pre_increment(idx)
            sp_FRT = pre_increment(idx)
            sp_FRL = pre_increment(idx)
            idx = sp_FRL + 2 ! 3-levels BELD3-based LAI fractions
            sp_CRL = pre_increment(idx)
            idx = sp_CRL + 3 ! 4-levels BELD3-based cumulative LAI fractions
         end if

         if (chm_strato_s == 'LINOZ' .and. (.not. llinoz)) then
            lin_c2 = pre_increment(idx)
            lin_c4 = pre_increment(idx)
            lin_c5 = pre_increment(idx)
            lin_c6 = pre_increment(idx)
            lin_c7 = pre_increment(idx)
         end if

! KPP solver internal time-step
         if (chm_pkg_gas_s(1:8) == 'ADOM2KPP' .or. chm_pkg_gas_s(1:5) == 'SAPRC') then
            var_step = pre_increment(idx)
         end if

         nb_fields = idx - start_index
         pkg_misc_idxinit = nb_fields

      end function pkg_misc_idxinit

!============================================================================
! Name           : pkg_misc_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
      subroutine pkg_misc_metainit()
       use clib_itf_mod,      only: clib_tolower
       implicit none

         integer           :: i, nsp, nn, istat, init_laik
         character(len=1)  :: zbin, init_c
         character(len=4)  :: sname, sname_out
         character(len=11) :: lamc
         character(len=12) :: lname

         if (chm_lu15_out_l) then
            species_master(sp_LU15) % per_name   =  "LU15"
            species_master(sp_LU15) % per_string =  "VN=land_use_15; ON=LU15; VD=15 category land use for dry dep; VS=A*15;VB=p0"
         end if
         species_master(sp_LU15) % ra_name    =  "LURA"
         species_master(sp_LU15) % ra_string  =  "VN=AREO_RESIST; ON=LURA; VD=aerodynamic resistance (s/m); VS=A*15;VB=v0"

         init_laik = 1
         if (chm_sat_seasons_l) then
            init_laik = 0
            chem_ent_vars(ce_LAI) % ent_name   =  "LAI_ENT"
            chem_ent_vars(ce_LAI) % ent_string =  "VN=LAI_ENT; ON=LAIE; VD=Monthly Leaf Area index input; VS=A*12;VB=E1"
         end if

         if (chm_get_be_emis_l .or. chm_dep_lai2d_l) then
            write(init_c, '(I1)') init_laik
            species_master(sp_LAI) % per_name     =  "LAIK"
            species_master(sp_LAI) % per_string   =  "VN=LAI_per_std; ON=LAIK; VD=Leaf Area index for BEIS3; VS=A;VB=p"//init_c
         end if

         if (sp_be_std > 0) then
            nn = sp_be_std
            do nsp = 1, num_be_std
               ! Firstly, summer standard biogenic base rate
               sname = 'ES' // be_std_descr(nsp)
               species_master(nn) % per_name     = sname
               species_master(nn) % per_string   = 'VN=BSE_' // sname // '; ON=' // sname // '; VD=' // be_std_name(nsp) // &
                                                  &   ' Biogenic std summer emission (g); VS=A;VB=p1'
               ! Secondly, winter standard biogenic base rate
               sname = 'EW' // be_std_descr(nsp)
               species_master(nn+1) % per_name   = sname
               species_master(nn+1) % per_string = 'VN=BWE_' // sname // '; ON=' // sname // '; VD=' // be_std_name(nsp) // &
                                                  &   ' Biogenic std winter emission (g); VS=A;VB=p1'
               nn = nn + 2
            end do
         end if

         if (chm_pkg_pm_s(1:3) == 'CAM') then
            species_master(sp_AF) % out_name    =  'AF'
            species_master(sp_AF) % out_string  =  'VN=AF_out_std; ON=AF; VD= PM2.5 ; VS=T;VB=v0'

            species_master(sp_AC) % out_name    =  'AC'
            species_master(sp_AC) % out_string  =  'VN=AC_out_std; ON=AC; VD= PM10 ; VS=T;VB=v0'
         end if

         species_master(sp_KTN) % out_name    =  "KTN"
         species_master(sp_KTN) % out_string  =  "VN=KTN_out_std; ON=KTN; VD=KT new for vert diff; VS=T;VB=v0"

         if (chm_aqhi_l) then
            species_master(sp_AQ25) % out_name    =  "AQ25"
            species_master(sp_AQ25) % out_string  =  "VN=AQ25_out_std; ON=AQ25; VD= Air Quality index 2.5 ; VS=T;VB=v0"

            species_master(sp_AQ10) % out_name    =  "AQ10"
            species_master(sp_AQ10) % out_string  =  "VN=AQ10_out_std; ON=AQ10; VD= Air Quality index 10.0; VS=T;VB=v0"
         endif

         if (chm_soa_s /= 'NIL') then
            species_master(sp_NWOC) % out_name   =  'NWOC'
            species_master(sp_NWOC) % out_string =  'VN=newsoa; ON=NWOC; VD=SOA CREATED IN TIMESTEP   ; VS=T;VB=V0'
         end if

         if (any(phyoutlist_S == 'noy')) then
            species_master(sp_NOy)  % out_name    =  'NOY'
            species_master(sp_NOy)  % out_string  =  'VN=NOY_out_std; ON=NOY; VD=NOy (ppb); VS=T;VB=V0'
         end if

         if (sp_density > 0) then
            species_master(sp_density) % out_name    =  "RHO"
            species_master(sp_density) % out_string  =  "VN=sp_density; ON=RHO; VD=chem air density; VS=T;VB=V0"
         end if

         if (sp_JNO2 > 0) then
            species_master(sp_JNO2) % out_name    =  "JNO2"
            species_master(sp_JNO2) % out_string  =  "VN=sp_JNO2; ON=JNO2; VD=NO2 J-value (s-1); VS=T;VB=v0"
         end if

         if (chm_pkg_pm_s /= 'NIL') then
            species_master(sp_FMET) % out_name    =  "FMET"
            species_master(sp_FMET) % out_string  =  "VN=FMET_out_std; ON=FMET; VD=Modulation factor for fugitive dust emissions; VS=A;VB=v0"
         end if

         if (chm_vit_l) then
            species_master(sp_CVKT) % per_name    =  "CVKT"
            species_master(sp_CVKT) % per_string  =  "VN=on_CVKT; ON=CVKT; VD=Onroad CVKT; VS=A;VB=p1"

            species_master(sp_MVKT) % per_name    =  "MVKT"
            species_master(sp_MVKT) % per_string  =  "VN=on_MVKT; ON=MVKT; VD=Onroad MVKT; VS=A;VB=p1"

            species_master(sp_TVKT) % per_name    =  "TVKT"
            species_master(sp_TVKT) % per_string  =  "VN=on_TVKT; ON=TVKT; VD=Onroad TVKT; VS=A;VB=p1"
         end if

! On-road area emissions
         if (chm_get_on_emis_l) then
            do nsp = 1, nb_dyn_tracers
              ! Skip if specie is sea-salt or area emissions is not defined
               if ((species_master(nsp) % ae_name == UNASSIGNED) .or. &
                   (species_master(nsp) % ae_name(1:3) == 'ESS')) cycle
               sname = "M" // species_master(nsp) % ae_name(2:4)
               species_master(nsp) % mae_name   = trim(sname)
               species_master(nsp) % mae_string = "VN=ON_"//trim(sname)//"; ON="//trim(sname)//  &
                                              & " ; VD=ONROAD MOBILE AREA EMISSION "//&
                                              & trim(sname)//" ; VS=A;VB=p1"
            end do

         end if

         if (chm_diag_aero_opt_l) then
            do nn = 1, nwl_aod
               write(zbin, '(Z1)') nn
               write(lamc, '(ES11.4)') aero_opt_wavel(nn)

               sname = "AOD"//zbin
               sname_out = sname
               istat = clib_tolower(sname_out)
               if (any(phyoutlist_S == sname_out)) then
                  species_master(nn + sp_OPTD - 1) % out_name   = sname
                  species_master(nn + sp_OPTD - 1) % out_string = "VN=COL_"//sname//"; ON="//sname// &
                                                                  & "; VD=Column AOD at "//lamc//" microns; VS=A;VB=v0"
               end if

               species_master(nn + sp_OPTD - 1) % per_name   = "OPD"//zbin
               species_master(nn + sp_OPTD - 1) % per_string = "VN=OPTD"//zbin//"; ON=OPD"//zbin//"; VD=Opt Depth at "//lamc//" microns; VS=T;VB=p0"

               species_master(nn + sp_ASYM - 1) % per_name   = "ASY"//zbin
               species_master(nn + sp_ASYM - 1) % per_string = "VN=ASYM"//zbin//"; ON=ASY"//zbin//"; VD=Asym Factor at "//lamc//" microns; VS=T;VB=p0"

               species_master(nn + sp_SSCA - 1) % per_name   = "SSC"//zbin
               species_master(nn + sp_SSCA - 1) % per_string = "VN=SSCA"//zbin//"; ON=SSC"//zbin//"; VD=Scatt. Albedo at "//lamc//" microns; VS=T;VB=p0"
            end do
!
!           Call Mie lookup table for direct feedback
            call mach_mie_optical_properties()
!
         end if

         if (chm_canopy_shading_l) then
!  Volatile bus variable for SOA formation under canopy
            if (chm_soa_s /= 'NIL') then
               lname = "cansoa"
               species_master(sp_NWOC) % vd_name   = "CNWC"
               species_master(sp_NWOC) % vd_string = "VN="//lname//"; ON=CNWC; VD=CANOPY SOA CREATED IN TIMESTEP ; VS=A*3;VB=V0"
            end if

            chem_ent_vars(ce_POP) % ent_name   = "POPU_ENT"
            chem_ent_vars(ce_POP) % ent_string = "VN=POPU_ENT; ON=PP_E; IN=POPU; VD=Population input; VS=A;VB=E1"

            chem_ent_vars(ce_FRT) % ent_name   = "FRT_BELD3"
            chem_ent_vars(ce_FRT) % ent_string = "VN=FRT_BELD3; ON=FR_E; IN=FRT; VS=A;VB=E1; " // &
                                              & "VD=BELD3-based forest fractional land use (E) "

            lname = "FRT_MASK"
            species_master(sp_FRT) % per_name   = "FRT"
            species_master(sp_FRT) % per_string = "VN="//lname//"; ON=FRT; VD=Forest canopy mask;VS=A;VB=p0"

            lname = "HC_satell"
            species_master(sp_HC) % per_name    = "FCH"
            species_master(sp_HC) % per_string  = "VN="//lname//"; ON=FCH; VD=Canopy height(m) from satellite retrievals ;VS=A;VB=p1"

            lname = "CLU_BELD3"
            species_master(sp_CLU) % per_name   = "CLU"
            species_master(sp_CLU) % per_string = "VN="//lname//"; ON=CLU; VD=BELD3-based Clumping Index;VS=A;VB=p1"

            do i = 1, 3
               write(zbin, '(Z1)') i
               sname = "F"//zbin//"R"
               lname = trim(sname)//"_BELD3"
               species_master(sp_FRL + i - 1) % per_name   = trim(sname)
               species_master(sp_FRL + i - 1) % per_string = "VN="//lname//"; ON="//sname// &
                                                    & "; VD=BELD3-based LAI fraction in #"  &
                                                    & //zbin//" forest layer;VS=A;VB=p1"
            end do

            do i = 1, 4
               write(zbin, '(Z1)') i
               sname = "C"//zbin//"R"
               lname = trim(sname)//"_BELD3"
               species_master(sp_CRL + i - 1) % per_name   = trim(sname)
               species_master(sp_CRL + i - 1) % per_string = "VN="//lname//"; ON="//sname// &
                                                    & "; VD=BELD3-based cum LAI fraction above #" &
                                                    & //zbin//" forest interface;VS=A;VB=p1"
            end do
         end if

         if (sp_STSE > 0) then
            species_master(sp_STSE) % per_name   =  "STMK"
            species_master(sp_STSE) % per_string =  "VN=Sat_seasn_msk_output; ON=STMK; VS=A;VB=p0; " // &
                                                &  "VD=Output Satellite LAI based Seasonal Mask"
         end if

         if (chm_strato_s == 'LINOZ' .and. (.not. llinoz)) then
            species_master(lin_c2) % per_name   = "TTCE"
            species_master(lin_c2) % per_string = "VN=LINOZ_c2; ON=TTCE; VD=Temperature climatology ERA5 2010-2016 (DegK); VS=T ; VB=p1; IN=TTCE"
            species_master(lin_c4) % per_name   = "PML"
            species_master(lin_c4) % per_string = "VN=LINOZ_c4; ON=PML ; VD=LINOZ Coeff. c4: P-L   Coeff     ; VS=T ; VB=p1; IN=PML "
            species_master(lin_c5) % per_name   = "DPDO"
            species_master(lin_c5) % per_string = "VN=LINOZ_c5; ON=DPDO; VD=LINOZ Coeff. c5: d(P-L)/dO3 Coeff; VS=T ; VB=p1; IN=DPDO"
            species_master(lin_c6) % per_name   = "DPDT"
            species_master(lin_c6) % per_string = "VN=LINOZ_c6; ON=DPDT; VD=LINOZ Coeff. c6: d(P-L)/dT       ; VS=T ; VB=p1; IN=DPDT"
            species_master(lin_c7) % per_name   = "DPDC"
            species_master(lin_c7) % per_string = "VN=LINOZ_c7; ON=DPDC; VD=LINOZ Coeff. c7: d(P-L)/dcolo3   ; VS=T ; VB=p1; IN=DPDC"
         end if

         if (var_step > 0) then
            species_master(var_step) % per_name   =  "HSTR"
            species_master(var_step) % per_string =  "VN=HSTART; ON=HSTR; VS=T; VB=p0; " // &
                                                  &  "VD=KPP solver internal timestep"
            if (chm_canopy_shading_l) then
! Use the dd_offset component to store hstart under canopy columns
               species_master(var_step) % dd_name   = "HSTC"
               species_master(var_step) % dd_string = "VN=HSTART_CAN; ON=HSTC; VS=A*3; VB=p0; " // &
                                                  &   "VD=CANOPY KPP solver internal timestep"
            end if
         end if

         return
      end subroutine pkg_misc_metainit

end module mach_pkg_misc_mod

