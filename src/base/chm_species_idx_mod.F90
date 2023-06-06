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

! This module defines all the species that will be eventually use in the model
! One integer per species/fields


module chm_species_idx_mod
#if defined(MACH_TENDENCIES)
   use chm_utils_mod,        only: MAX_DEBUG_VAR, MAX_TRACERS
   use chm_species_info_mod, only: tendency_idx
#else
   use chm_utils_mod,        only: MAX_DEBUG_VAR
#endif
   save
! These variables will be used to access the master array of fields (species or else)
! Please keep alphabetical order
! MACH specific variables
   integer(kind=4) :: sp_density = 0

! Gas indices - ADOM2-only variables
   integer(kind=4) :: sp_AFG1 = 0
   integer(kind=4) :: sp_AFG2 = 0
   integer(kind=4) :: sp_ALK3 = 0
   integer(kind=4) :: sp_ALK4 = 0
   integer(kind=4) :: sp_ALKA = 0
   integer(kind=4) :: sp_ALKE = 0
   integer(kind=4) :: sp_ALD2 = 0
   integer(kind=4) :: sp_ARO1 = 0
   integer(kind=4) :: sp_ARO2 = 0
   integer(kind=4) :: sp_AROM = 0
   integer(kind=4) :: sp_BZO = 0
   integer(kind=4) :: sp_C3H8 = 0
   integer(kind=4) :: sp_CCHO = 0
   integer(kind=4) :: sp_CH4 = 0
   integer(kind=4) :: sp_CO = 0
   integer(kind=4) :: sp_CRES = 0
   integer(kind=4) :: sp_CRG1 = 0
   integer(kind=4) :: sp_CRG2 = 0
   integer(kind=4) :: sp_DIAL = 0
   integer(kind=4) :: sp_ETHE = 0
   integer(kind=4) :: sp_H2O2 = 0
   integer(kind=4) :: sp_HCHO = 0
   integer(kind=4) :: sp_HNO3 = 0
   integer(kind=4) :: sp_HNO4 = 0
   integer(kind=4) :: sp_HO2 = 0
   integer(kind=4) :: sp_HONO = 0
   integer(kind=4) :: sp_IEPX = 0
   integer(kind=4) :: sp_IPRD = 0
   integer(kind=4) :: sp_IOOH = 0
   integer(kind=4) :: sp_ISOP = 0
   integer(kind=4) :: sp_MCO3 = 0
   integer(kind=4) :: sp_MEK = 0
   integer(kind=4) :: sp_MGL0 = 0
   integer(kind=4) :: sp_MGLY = 0
   integer(kind=4) :: sp_N = 0
   integer(kind=4) :: sp_N2 = 0
   integer(kind=4) :: sp_N2O = 0
   integer(kind=4) :: sp_N2O5 = 0
   integer(kind=4) :: sp_NH3 = 0
   integer(kind=4) :: sp_NO = 0
   integer(kind=4) :: sp_NO2 = 0
   integer(kind=4) :: sp_NO3 = 0
   integer(kind=4) :: sp_O = 0
   integer(kind=4) :: sp_O3 = 0
   integer(kind=4) :: sp_O1D = 0
   integer(kind=4) :: sp_O3P = 0
   integer(kind=4) :: sp_OH = 0
   integer(kind=4) :: sp_OLE1 = 0
   integer(kind=4) :: sp_OLE2 = 0
   integer(kind=4) :: sp_OSD = 0
   integer(kind=4) :: sp_PAN = 0
   integer(kind=4) :: sp_PAN2 = 0
   integer(kind=4) :: sp_PRD2 = 0
   integer(kind=4) :: sp_R2O2 = 0
   integer(kind=4) :: sp_RCHO = 0
   integer(kind=4) :: sp_RCO3 = 0
   integer(kind=4) :: sp_RNO3 = 0
   integer(kind=4) :: sp_RN3a = 0
   integer(kind=4) :: sp_RO2 = 0
   integer(kind=4) :: sp_RO2N = 0
   integer(kind=4) :: sp_RO2R = 0
   integer(kind=4) :: sp_ROOH = 0
   integer(kind=4) :: sp_SO4 = 0
   integer(kind=4) :: sp_TERP = 0
   integer(kind=4) :: sp_TOLU = 0
   integer(kind=4) :: sp_SO2 = 0
   integer(kind=4) :: sp_XOOH = 0
   integer(kind=4) :: sp_yIOH = 0

! Fixed (or steady state product-only) gas species
   integer(kind=4) :: sp_C2H6 = 0
   integer(kind=4) :: sp_CO2 = 0
   integer(kind=4) :: sp_DUM = 0
   integer(kind=4) :: sp_H2 = 0
   integer(kind=4) :: sp_H2O = 0
   integer(kind=4) :: sp_M = 0
   integer(kind=4) :: sp_O2 = 0
   integer(kind=4) :: sp_XC = 0
   integer(kind=4) :: sp_XN = 0

! Biogenic (BEIS) emission module variables
   integer(kind=4) :: sp_LAI = 0   ! Leaf Area Index
   integer(kind=4) :: sp_be_std = 0   ! Standard biogenic emissions starting index

   integer(kind=4) :: sp_LU15 = 0  ! 15-Categories Land use fractions
   integer(kind=4) :: sp_KTN = 0   ! Modified eddy diffusivity

! Fugitive dust emissions modulation factor
   integer(kind=4) :: sp_FMET = 0  ! Fugitive aerosol emissions modulation factor

! Organic aerosol related
   integer(kind=4) :: sp_NWOC = 0  ! SOA yield

! For output only
   integer(kind=4) :: sp_AQ25 = 0  ! PM2.5 based AQHI index
   integer(kind=4) :: sp_AQ10 = 0  ! PM10 based AQHI index
   integer(kind=4) :: sp_JNO2 = 0  ! NO2 photolysis rate (s-1) (at the surface)

! LINOZ Coefficients
   integer(kind=4) :: lin_c2 = 0 ! Clim. TT for LINOZ
   integer(kind=4) :: lin_c4 = 0 ! LINOZ Coeff.  c4
   integer(kind=4) :: lin_c5 = 0 ! LINOZ Coeff.  c5
   integer(kind=4) :: lin_c6 = 0 ! LINOZ Coeff.  c6
   integer(kind=4) :: lin_c7 = 0 ! LINOZ Coeff.  c7

!  Aerosol-module
   integer(kind=4) :: sp_AERO = 0 ! Aerosol tracers starting index
   integer(kind=4) :: sp_AC = 0   ! total coarse fraction
   integer(kind=4) :: sp_AF = 0   ! total fine fraction

! Diagnostic for Linoz
   integer(kind=4) :: sp_NOY = 0
! 3D Total columns
   integer(kind=4) :: sp_tend_o3 = 0
   integer(kind=4) :: sp_diff_o3 = 0
   integer(kind=4) :: sp_col_o3 = 0
   integer(kind=4) :: sp_col_no2 = 0
   integer(kind=4) :: sp_col_no = 0
   integer(kind=4) :: sp_col_hu = 0
! 2D Tropo. columns
   integer(kind=4) :: sp_tcol_o3 = 0
   integer(kind=4) :: sp_tcol_no2 = 0
   integer(kind=4) :: sp_tcol_no = 0
   integer(kind=4) :: sp_tcol_pan = 0
   integer(kind=4) :: sp_tcol_hno3 = 0
   integer(kind=4) :: sp_tcol_n2o5 = 0
   integer(kind=4) :: sp_tcol_so2 = 0
   integer(kind=4) :: sp_tcol_so4 = 0
   integer(kind=4) :: sp_tcol_co = 0
   integer(kind=4) :: sp_tcol_hcho = 0
   integer(kind=4) :: sp_tcol_isop = 0

! Onroad VKT (Vehicle-Kilometers-Traveled) fields
   integer(kind=4) :: sp_CVKT = 0  ! Cars VKT
   integer(kind=4) :: sp_MVKT = 0  ! Mid-sized VKT
   integer(kind=4) :: sp_TVKT = 0  ! Trucks VKT

! Diagnostic for wet deposition
   integer(kind=4) :: sp_WSO3 = 0
   integer(kind=4) :: sp_CAT = 0
   integer(kind=4) :: sp_HCO3 = 0
   integer(kind=4) :: sp_HION = 0
   integer(kind=4) :: sp_WH2O = 0

! Diagnostic for total aerosol
   integer(kind=4) :: sp_NAC = 0  ! Diagnostic PM10 number density
   integer(kind=4) :: sp_NAF = 0  ! Diagnostic PM2.5 number density
   integer(kind=4) :: sp_NAT = 0  ! Diagnostic PM Total number density
   integer(kind=4) :: sp_NUMC = 0 ! Aerosol number concentration per bin

! Diagnostic aerosol optical properties
   integer(kind=4) :: sp_OPTD = 0 ! Optical depth
   integer(kind=4) :: sp_ASYM = 0 ! Asymmetric factor
   integer(kind=4) :: sp_SSCA = 0 ! Single-scattering albedo

! Forest canopy shading
   integer(kind=4) :: sp_HC  = 0 ! Forest canopy height from satellite retrievals
   integer(kind=4) :: sp_CLU = 0 ! Clumping index from BELD3 data
   integer(kind=4) :: sp_FRT = 0 ! Forest canopy mask
   integer(kind=4) :: sp_FRL = 0 ! fraction of LAI in each 3 forest layer from BELD3 data
   integer(kind=4) :: sp_CRL = 0 ! cumulative fraction of LAI above each (4) forest layer interface
   integer(kind=4) :: ce_FRT = 0 ! Forest fraction from BELD3 data
   integer(kind=4) :: ce_POP = 0 ! Population per grid cell (Entry)
!
! LAI
   integer(kind=4) :: ce_LAI  = 0 ! Monthly satellite-derived LAI input (Entry bus)
   integer(kind=4) :: sp_STSE = 0 ! Satellite LAI based seasonal mask
!
! KPP internal time-step
   integer(kind=4) :: var_step = 0
!
! Debug output variable
   integer(kind=4) :: dbg_2D(MAX_DEBUG_VAR) = 0
   integer(kind=4) :: dbg_3D(MAX_DEBUG_VAR) = 0

#if defined(MACH_TENDENCIES)
! Tendencies fields
! First dimension = looping index (idx:idx+nb_tracers)
! Second dimension = corresponding index to species_master (to link this field to it's original tracer)
   type(tendency_idx) :: tendency_fields(MAX_TRACERS)
#endif

end module chm_species_idx_mod
