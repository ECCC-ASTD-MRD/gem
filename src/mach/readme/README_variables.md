Names of GEM-MACH Output Variables
========================================================================================

In GEM-MACH, chemical gas and aerosol species are advected as passive tracers. Therefore, their concentrations can be output from the dynamics part of the model (dm) in the units that are expected in the dynamical core of the model (ug/kg). Because these units are less commonly used in atmospheric chemistry, physics part of the model has a capability of outputing the concentrations in the more commonly used units (pm output in units of ppbV for gas concentrations and ug/m3 for aerosol concentrations).


Gas Species Concentrations and Emissions
----------------------

| dm Output (ug/kg) | pm Output (ppbV) (1) | Emissions (g/s) (2) | Biogenic Emissions (g/s) | Description | Mechanism / Note |
| :------ | :----- | :----- | :----- | :------: | :-----: |
| TAF1 | AFG1 | EAF1 |      | Photoreactive monounsaturated dicarbonyl   | SAPRC07C,S only |
| TAF2 | AFG2 |      |      | Lumped monounsaturated dicarbonyl aromatic | SAPRC07C,S only |
| TALD | ALD2 | EALD |      | Acetaldehyde and higher aldehydes          | ADOM2 only |
| TAK3 | ALK3 | EA3  | BPAA | Mid-sized alkanes                          | SAPRC07C,S only |
| TAK4 | ALK4 | EA4  |      | Mid to Large alkanes                       | SAPRC07C,S only |
| TA3  | ALKA | EA3  | BPAA | C4+ alkanes                                | ADOM2 only |
| TA2  | ALKE | EA2  | BPAE | C3+ alkenes                                | ADOM2 only |
| TAR1 | ARO1 | EAR1 |      | Moderate reacting aromatics                | SAPRC07C,S only |
| TAR2 | ARO2 | EAR2 |      | Fast reacting aromatics (kOH > 2.E4ppm/min | SAPRC07C,S only |
| TARO | AROM | EARO |      | Lumped higher aromatics                    | ADOM2 only |
| TBZO | BZO  |      |      | Phenoxy radical                            | |
| TC38 | C3H8 | EC38 |      | Propane and other slowly reacting organics | ADOM2 only |
| TCHC | CCHO | ECHC |      | Acetaldehyde                               | SAPRC07C,S only |
| TCH4 | CH4G | ECH4 |      | Methane                                    | Available only when chm\_active\_ch4\_l=.true. |
| TCO  | CO_P | ECO  |      | Carbon monoxide                            | |
| TCR1 | CRG1 |      |      | Criegee radical #1                         | ADOM2 only |
| TCR2 | CRG2 |      |      | Criegee radical #2                         | ADOM2 only |
| TCRE | CRES | ECRE |      | Lumped creosols                            | |
| TDIA | DIAL |      |      | General dicarbonyl                         | ADOM2 only |
| TETH | ETHE | EETH |      | Ethene                                     | |
| TH22 | H2O2 |      |      | Hydrogen peroxide                          | |
| THCH | HCHO | EHCH |      | Formaldehyde                               | |
| THN3 | HNO3 |      |      | Nitric acid                                | |
| THN4 | HNO4 |      |      | Pernitric acid                             | |
| THO2 | HO2  |      |      | Hydroperoxyl radical                       | |
| THON | HONO | EHON |      | Nitrous acid                               | |
| TIPE | IEPX |      |      | Isoprene epoxydiols                        | SAPRC07C,S only |
| TIPR | IPRD | EIPR |      | Lumped isoprene product species            | SAPRC07C,S only |
| TIOH | IOOH |      |      | Hydroxyhydroperoxide from Isoprene         | SAPRC07C,S only |
| TISO | ISOP | EISO | BPIO | Isoprene                                   | |
| TMC3 | MCO3 |      |      | Acetyl Peroxy (CH3CO3) Radical             | |
| TMEK | MEK  | EMEK |      | Methyl ethyl ketone and higher ketones     | |
| TMGL | MGLY |      |      | Methyl glyoxal (ADOM2)                     | |
| TMGL | MGLY | EMGL |      | Methyl glyoxal                             | SAPRC07C,S only |
| TAN  | N_A  |      |      | Atomic nitrogen                            | SAPRC07C,S only |
| TN2O | N2OG |      |      | Nitrous oxide                              | SAPRC07CS only |
| TN25 | N2O5 |      |      | Dinitrogen pentoxide                       | |
| TNH3 | NH3  | ENH3 |      | Ammonia                                    | |
| TNO  | NO   | ENO  | BPNO | Nitric oxide                               | |
| TNO2 | N2   | ENO2 |      | Nitrogen dioxide                           | |
| TNO3 | NO3  |      |      | Nitrate radical                            | |
| TO1D | O1D  |      |      | Excited-state oxygen atom                  | SAPRC07C,S only |
| TO3P | O3P  |      |      | Ground-state oxygen atom                   | SAPRC07C,S only |
| TO   | O3P  |      |      | Ground-state oxygen atom                   | |
| TOSD | O1D  |      |      | Excited-state oxygen atom                  | ADOM2 only |
| TO3  | O3   |      |      | Ozone                                      | |
| TOH  | OH   |      |      | Hydroxyl radical                           | |
| TOL1 | OLE1 | EOL1 |      | Moderate reacting alkenes                  | SAPRC07C,S only |
| TOL2 | OLE2 | EOL2 | BPAE | Fast reacting alkenes(kOH > 7.E4ppm/min    | SAPRC07C,S only |
| TPAN | PAN  |      |      | Peroxyacetyl nitrate (PAN)                 | |
| TPN2 | PAN2 |      |      | PPN and other higher alkyl PAN analogues   | SAPRC07C,S only |
| TPRD | PRD2 | EPR2 |      | Ketones and other oxygenated non-aldehyde  | SAPRC07C,S only |
| TR22 | R2O2 |      |      | General organic peroxyl radical #2         | |
| TRCH | RCHO | ERCH |      | Propionaldehyde and larger aldehydes       | SAPRC07C,S only |
| TRCO | RCO3 |      |      | Acyl Peroxy Radical (C3+)                  | SAPRC07C,S only |
| TR2N | RO2N |      |      | Alkyl nitrate => organic peroxyl radical   | |
| TR2R | RO2R |      |      | General organic peroxyl radical #1         | |
| TRN3 | RNO3 |      |      | Organic nitrate (RNO3)                     | |
| TRO2 | RO2  |      |      | Total RO2 radicals                         | |
| TROO | ROOH |      |      | Lumped Organic peroxide                    | |
| TSO2 | S2   | ESO2 |      | Sulphur dioxide                            | |
| TSO4 | S4   | ESO4 |      | Sulphuric acid gas                         | |
| TTRP | TERP | ETRP | BTRP | Monoterpenes - Emitted Lumped model VOCs   | SAPRC07C,S only |
| TTOL | TOLU | ETOL |      | Toluene and other mono-alkylbenzenes       | |
| TIYH | yIOH |      |      | yIOH                                       | SAPRC07C,S only |

Notes:
(1) By assigning values to 'gas\_ppb\_out\_list\_s namelist key, choose which gas species diagnostic variables (second column) will be made available to request in output.
(2) To get the names of on-road gas emissions fields, replace the first letter "E" of the name in the third, "Emissions", column with the letter "M".

Aerosol Species Mass Concentrations and Emissions
----------------------

| dm Output (ug/kg) (1) | pm Output (ug/m3) (2) | Area Emission (g/s) (3) | Fugitive Emission (g/s) (3) | On-road Emissions (g/s) | Description | Note |
| :------ | :----- | :----- | :----- | :----- | :------: | :-----: |
| TSUX | SUYY | ESUX | FSUX | MSUX | Sulphate |  |
| TSSX | SSYY | ESSX (4)|  |  | Sea Salt |  |
| TOCX | OCYY |  |  |  | Secondary Organic Carbon |  |
| TNIX | NIYY | ENTX | FNTX | MNTX | Nitrates |  |
| TAMX | AMYY | EAMX | FAMX | MAMX | Ammonium |  |
| TCMX | CMYY | ECMX | FCMX | MCMX | Crustal Material (Soil-dust) |  |
| TECX | ECYY | EECX | FECX | MECX | Black Carbon |  |
| TPCX | PCYY | EPCX | FPCX | MPCX | Primary Organic Carbon |  |

Notes:<br>
(1) 'X' in the first, third, fourth and fifth column stands for the bin number of the size distribution: 1 and 2 for 2-bin model version, and 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, and C for 12-bin model version.<br>
(2) 'YY' in the second column stands for the PM size: 'AS' for Sub-micron (<=1um), 'AF' for Fine (<=2.5um), 'AC' for Coarse (<=10um), and 'AT' for Total (all sizes). 2-bin model has only 'AF' and 'AC' defined, while 12-bin model has all four.<br>
This output is available only when chm\_diag\_aerosols\_l is set to .true.<br>
(3) The two largest bins of the 12-bin size distribution (11th and 12th, or B and C) have neither area nor major point emissions available, while the 8th bin lacks only major point emissions.<br>
(4) Sea Salt emissions are evaluated for all bin sizes within the model at the runtime.


Diagnostic Aerosol Concentrations
----------------------

| Field Name (1) | Description | Note |
| :------ | :------: | :-----: |
| AF | PM2.5 Mass Concentration (ug/m3) |  |
| AC | PM10 Mass Concentration (ug/m3) |  |
| NAFP | PM2.5 (Fine) Number concentration (/m3) | (2) |
| NACP | PM10 (Coarse) Number concentration (/m3) | (2) |
| NATP | PM Total Number Concentration (/m3) | (2) |
| NUMX (1) | Bin-Specific Aerosol Number Concentration (/m3) | (2) |

Notes:<br>
(1) As above, 'X' stands for the bin number of the size distribution: 1 and 2 for 2-bin model version, and 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, and C for 12-bin model version.<br>
(2) All but AF and AC are available when chm\_diag\_aerosols\_l is set to .true. AF and AC do not require a namelist key to be set.


Diagnostic Dry Deposition Fields for Gas Species
----------------------

| Deposition Flux (moles/m2) | Deposition Velocity (m/s) | Total surface resistance (s/m) | Quasi-laminar sublayer resistance (Wesely's Rb term) (s/m) | Description | Note |
| :------ | :----- | :----- | :----- | :------: | :-----: |
| DALD | VALD | ADRC | ADRB | Acetaldehyde and higher aldehydes |  |
| DA3 | VA3 | A3RC | A3RB | C4+ alkanes |  |
| DA2 | VA2 | A2RC | A2RB | C3+ alkenes |  |
| DARO | VARO | ARRC | ARRB | Lumped higher aromatics |  |
| DC38 | VC38 | C3RC | C3RB | Propane and other slowly reacting organics |  |
| DCRE | VCRE | CRRC | CRRB | Lumped creosols |  |
| DDIA | VDIA | DIRC | DIRB | General dicarbonyl |  |
| DETH | VETH | ETRC | ETRB | Ethene |  |
| DH22 | VH22 | H2RC | H2RB | Hydrogen peroxide |  |
| DHCH | VHCH | HCRC | HCRB | Formaldehyde |  |
| DHN3 | VHN3 | H3RC | H3RB | Nitric acid |  |
| DHO2 | VHO2 | HORC | HORB | Hydroperoxyl radical |  |
| DHNO | VHNO | HNRC | HNRB | Nitrous acid |  |
| DISO | VISO | IORC | IORB | Isoprene |  |
| DMC3 | VMC3 | MCRC | MCRB | Acetyl Peroxy (CH3CO3) Radical |  |
| DMEK | VMEK | MKRC | MKRB | Methyl ethyl ketone and higher ketones |  |
| DMGL | VMGL | MGRC | MGRB | Methyl glyoxal |  |
| DNH3 | VNH3 | NHRC | NHRB | Ammonia |  |
| DNO | VNO | NORC | NORB | Nitric oxide |  |
| DNO2 | VNO2 | N2RC | N2RB | Nitrogen dioxide |  |
| DO3 | VO3 | O3RC | O3RB | Ozone |  |
| DPAN | VPAN | PARC | PARB | Peroxyacetyl nitrate (PAN) |  |
| DRN3 | VRN3 | R3RC | R3RB | Organic nitrate (RNO3) |  |
| DRO2 | VRO2 | RORC | RORB | Total RO2 radicals |  |
| DROO | VROO | RHRC | RHRB | Lumped Organic peroxide (ROOH) |  |
| DSO2 | VSO2 | S2RC | S2RB | Sulphur dioxide |  |
| DSO4 | VSO4 | S4RC | S4RB | Sulphuric acid gas |  |
| DTOL | VTOL | TLRC | TLRB | Toluene and other mono-alkylbenzenes |  |
| DNOY |  |  |  | Total NOy deposition | Available only when chm\_noy\_out\_l is set to .true. |

Note: Available when chm\_diag\_drydep\_L is set to .true.


Diagnostic Dry Deposition Fields for Aerosol Species
----------------------

| Deposition Flux (moles/m2) | Description | Note |
| :------ | :------: | :-----: |
| DAM | Dry Deposition of Ammonium |  |
| DNI | Dry Deposition of Nitrate |  |
| DSU | Dry Deposition of Sulphate |  |
| DEC | Dry Deposition of Black Carbon |  |
| DPOM | Dry Deposition of Primary Organic Carbon | |
| DSOM | Dry Deposition of Secondary Organic Carbon | |
| DCM | Dry Deposition of Crustal Material (Soil-dust) | |
| DSS | Dry Deposition of Sea Salt |  |

Note: Available when chm\_diag\_drydep\_L is set to .true.


Diagnostic Wet Deposition Fields
----------------------

| Deposition Flux (moles/m2) | Description | Note |
| :------ | :------: | :-----: |
| WH22 | Wet deposition of H2O2 (hydrogen peroxide) |  |
| WROO | Wet deposition of ROOH (organic peroxide) |  |
| WSO3 | Wet deposition of HSO3- |  |
| CAT | Wet deposition of base cations |  |
| HCO3 | Wet deposition of HCO3(-) |  |
| HION | Wet deposition of H(+) |  |
| WH2O | Wet deposition of aerosol water |  |
| WNH4 | Wet deposition of Ammonium particles |  |
| WNO3 | Wet deposition of Nitrate particles |  |
| WSO4 | Wet deposition of Sulphate particles |  |
| WEC | Wet deposition of Black Carbon particles |  |
| WPOM | Wet deposition of Primary Organic Carbon particles |  |
| WSOM | Wet deposition of Secondary Organic Carbon particles |  |
| WCM | Wet deposition of Crustal Material (Soil-dust) particles |  |
| WSS | Wet deposition of Sea Salt particles |  |

Note: Available when chm\_diag\_wetdep\_L is set to .true.

Bidirectional Flux Fields (ammonia)
----------------------

| Field Name | Description | Note |
| :------ | :------: | :-----: |
| NHVG | Deposition velocity for ground pathway (s/m) | output |
| NHBD | NH3 emissions from bidirectional flux scheme (g/s) | output |
| NHGP | NH3 static ground emissions potential | input |
| NHGD | NH3 dynamic ground emissions potential | input |
| NHTA | NH3 bidirectional flux time scale (hours) | output |
| SLPH | Soil PH | input |

Miscellaneous Fields
----------------------

| Field Name | Description | Note |
| :------ | :------: | :-----: |
| AQ25 | Air Quality index 2.5 | Available when chm\_aqhi\_l is set to .true. |
| AQ10 | Air Quality index 10.0 | Available when chm\_aqhi\_l is set to .true. |
| LU15 | 15 category land use for dry dep |  |
| LURA | aerodynamic resistance (s/m) |  |
| LAIK | Leaf Area index for BEIS3 |  |
| LAIE | Monthly Leaf Area index | input |
| ESNO | Biogenic std summer emission NO(g/s) |  |
| ESIO | Biogenic std summer EMISS ISOP(g/s) |  |
| ESVO | Biogenic std summer EMISS VOC(g/s) |  |
| ESMO | Biogenic std summer EMISS MONO(g/s) |  |
| EWNO | Biogenic std winter emission NO(g/s) |  |
| EWIO | Biogenic std winter EMISS ISOP(g/s) |  |
| EWVO | Biogenic std winter EMISS VOC(g/s) |  |
| EWMO | Biogenic std winter EMISS MONO(g/s) |  |
| KTN | KT new for vertical diffusion |  |
| NWOC | SOA CREATED IN TIMESTEP |  |
| NOY | NOy (ppb) |  |
| RHO | chem air density |  |
| JNO2 | NO2 J-value (s-1) |  |
| FMET | Modulation factor for fugitive dust emissions |  |
| CVKT | Vehicle Kilometers Traveled for cars | Required as hourly input when vehicle induced turbulence activated (i.e. chm\_vit\_l=.true.)  |
| MVKT | Vehicle Kilometers Traveled for mid-sized vehicles | Required as hourly input when vehicle induced turbulence activated (i.e. chm\_vit\_l=.true.)  |
| TVKT | Vehicle Kilometers Traveled for trucks | Required as hourly input when vehicle induced turbulence activated (i.e. chm\_vit\_l=.true.)  |
| AODW | Column aerosol optical depth (AOD) at specified Walength | 'W' stands for a hexa-decimal formatted index of the specified mid-band wavelengths provided in aero\_opt\_wavel (can be up to 15 different values) |
| OPDW | Optical Depth at specified wavelength | 'W' stands for a hexa-decimal formatted index of the specified mid-band wavelengths provided in aero\_opt\_wavel (can be up to 15 different values) |
| ASYW | Asym Factor at specified wavelength | 'W' stands for a hexa-decimal formatted index of the specified mid-band wavelengths provided in aero\_opt\_wavel (can be up to 15 different values) |
| SSCW | Single-Scattering Albedo at specified wavelength | 'W' stands for a hexa-decimal formatted index of the specified mid-band wavelengths provided in aero\_opt\_wavel (can be up to 15 different values) |
| POPU | Grid-sized Human Population input | input |
| PP\_E | Grid-sized Human Population input | output |
| FCH | Forest canopy height (m) from satellite retrievals | input |
| CLU | Clumping Index from BELD3 data | input |
| FRT | Forest fractional land use from BELD3 data | input |
| FR\_E | Forest fractional land use from BELD3 data | output |
| FRT | Forest canopy mask | output |
| FXR | Fraction of LAI in each 3 forest layers from BELD3 data | input <br> X stands for the forest layer number, counting down from the top.  |
| CXR | Cumulative fraction of LAI above each of 4 forest layer interface from BELD3 data | input <br> X stands for the interface number, counting down from the top (i.e. at heights of 1.0, 0.5, 0.2 and 0. canopy height). |
| STMK | Satellite LAI-based Seasonal Mask |  |
| O3FK | Ozone 3D Climatology Fortuin-Kelder 1998 |  |
| TTCE | Temperature climatology ERA5 2010-2016 (DegK) |  |
| PML | LINOZ Coefficient c4: P-L |  |
| DPDO | LINOZ Coefficient c5: d(P-L)/dO3 |  |
| DPDT | LINOZ Coefficient c6: d(P-L)/dT |  |
| DPDC | LINOZ Coefficient c7: d(P-L)/dcolo3 |  |
| NHVG | Deposition velocity for ground pathway (s/m) | Used in NH3 bidirectional flux (output) |
| NHBD | NH3 emissions from bidirectional flux scheme (g/s) | Used in NH3 bidirectional flux (output) |
| NHGP | NH3 ground emissions potential | Used in NH3 bidirectional flux (input) |
| HSTR | KPP-solver internal time step on the model vertical levels ; Used in gas mechanisms using KPP rodas solver (output) |
| HSTC | KPP-solver internal time step on the additional canopy vertical levels ; Used in gas mechanisms using KPP rodas solver (output) |


Debug variables
---------------

Depending on the debugging needs, developers can assign values of internal model variables to the debug variables.

To output these debug variables, besides the request in `outcfg.out`, a number of the debug variables is required to be provided in the namelist.

Names of these debug variables are automatically created based on their requested number:
* The first 3 characters are `2DB` for 2D variables, and `3DB` for 3D variables
* The forth character depends on the number of variables requested.
 * For the maximum of 18 variables, the values are: `1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I`.
 * For any number less than 18, the values are a subset of the above list, starting with 1 and ending with the chacarcter corresponding to the requested number of variables.

Examples:

* 2DB1, 2DB2, ... 2DBA, 2DBB, ...
* 3DB1, 3DB2, ... 3DBA, 3DBB, ...

For more details on the debug variables, read [README_compilation](README_compilation.md).

