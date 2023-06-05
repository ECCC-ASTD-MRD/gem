GEM-MACH Chemistry Namelist Settings  (chemistry\_cfgs group) Plus Relevant GEM Namelist Settings
========================================================================================


#### GEM Keys Necessary When Running GEM-MACH

| Name | Group | Description | Default Value | Type |
| :------ | :----- | :----- | :------: | :-----: |
| tr3d\_anydate\_l | gem\_cfgs | Controls whether the tracers's validity time matches analysis <br> Has to be set to **.true.** in GEM-MACH | .false. | logical |
| kfcprod | convection\_cfgs | Compute production terms for Kain-Fritsch scheme <br> Has to be set to **.true.** in GEM-MACH   | .false. | logical |


#### GEM Keys Recommended When Running GEM-MACH

| Name | Group | Description | Default Value | Type |
| :------ | :----- | :----- | :------: | :-----: |
| indiag\_list\_s | physics\_cfgs | Controls initialization of the diagnostic-level values for the dynamical fields in Physics <br> Should be set to **'UU','VV','TT','HU','QC',** in GEM-MACH | all dynamic fields<br>(i.e., 'UU', 'VV', 'TT', 'HU', all of the hydrometeors, and all of the chemical tracers) | character(len=32) |


#### Chemistry Keys Defining GEM's Chemistry-Related Behavior

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| chm\_master | Master key to activate chemistry (.true.) or to run GEM without chemistry as if in chm\_stub mode (.false.)| .true. | logical | Usually set to .true. <br>When this key is set to .false., keys for chemistry feedback onto physics (chm\_direct\_l and chm\_indirect\_l) are overwritten to .false. too.|
| chm\_step\_factor | Chemistry to physics timestep factor (an integer multiple, 3 is often used) | -1 | integer | Set to 3 in RAQDPS022 |
| chm\_mass\_s | Name of mass conservation scheme applied to advection of **all** tracers.<br>Options:<br>  - 'BC' (Bermejo-Conde)<br>  - 'NIL' | 'BC' | character | Set to default in RAQDPS022.<br>To set for individual tracers, modify tr3d\_list\_s in gem\_cfgs group, as discussed in [GEM-MACH's wiki page](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH\_v2\_with\_GEM\_4.8#Monotonicity.2C\_Mass\_Conservation.2C\_and\_Initialization). |
| chm\_mono\_s | Name of monotonicity scheme applied to advection of all tracers.<br>Options:<br> - 'CLIP' <br>  - 'ILMC'  (Iterative Locally Mass Conserving)<br>  - 'NIL' | 'CLIP' | character | Set to 'ILMC' in RAQDPS022. <br> For individual tracers, modify tr3d\_list\_s in gem\_cfgs group, as discussed in [GEM-MACH's wiki page](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH\_v2\_with\_GEM\_4.8#Monotonicity.2C\_Mass\_Conservation.2C\_and\_Initialization).  See [deGrandpre et al. (2016)] (https://wiki.cmc.ec.gc.ca/images/c/c0/MWR-D-15-0142.pdf) for more information.|


#### Chemistry Vertical Boundary Settings

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| chm\_vert\_diff\_s | Name of vertical diffusion algorithm<br>Options:<br> - 'RPNPHY' <br>  - 'RPNPHY\_I' <br> - 'FLUX' <br>  - 'BOUNDARY' <br>  - 'NIL' | undefined | character | Set to 'RPNPHY' in RAQDPS022. |
| nk\_start | Top model level at which tropospheric gas-phase chemistry scheme (cf. chm\_pkg\_gas\_s) is activated | 1 | integer | Set to 33 (on L84) in RAQDPS022. <br>In combination with LINOZ, it is recommended to set to 33. |
| nk\_start\_pm | Top model level at which aerosol chemistry is activated | 1 | integer | Set to 33 in RAQDPS022. <br>In combination with LINOZ, it is recommended to set to 33. |


#### Modification of Meteorological Fields Used in Chemistry

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| chm\_kt\_minmax (3) | Min, max, urban-min (respectively; i.e., vector of length 3) limits for vertical diffusion coefficient | 0.1, 1500.0, 1.95 | real | Set to default values in RAQDPS022. <br>|
| chm\_pblh\_min\_l | Impose minimum value of 100 m for PBL height (.true.) | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_vert\_diff\_s | Vertical diffusion option | undefined| | 'RPNPHY','RPNPHY_I','RPNPHY_U','NIL',| Set to 'RPNPHY' in RAQDPS023 |

#### Chemical Schemes Settings

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| chm\_model\_s | Name of chemistry model.<br>Options:<br>  - 'MACH' | undefined | character | Set to 'MACH' in RAQDPS022. |
| chm\_bkgd\_ch4 | Background CH4 concentration in ppmV | 1.8 ppm | real | Not present in RAQDPS022. |
| chm\_bkgd\_co2 | Background CO2 concentration in ppmV | 400.0 ppm | real | Not present in RAQDPS022. |
| chm\_pkg\_gas\_s | Name of gas-phase species package.<br>Options:<br>  - 'ADOM2' <br>  - 'ADOM2KPPB' <br>  - 'ADOM2KPP' <br>  - 'SAPRC07CS' <br>  - 'SAPRC07C' <br>  - 'NIL' | undefined | character | Set to 'ADOM2' in RAQDPS022. <br>When this key is set to 'NIL': <br>a. chm\_strato\_s and chm\_soa\_s are overwritten to 'NIL'; <br>b. chm\_messy\_jval\_l is overwritten to .false.; <br>c. diagnostics keys (i.e. chm\_aqhi\_l, chm\_diag\_accum\_l and chm\_diag\_column\_l) are overwritten to .false.; and <br>d. if chm\_pkg\_pm\_s is set to one of the two CAM versions, seven (7) gas species required by CAM (H2O2, HNO3, NH3, O3, ROOH, SO2 and SO4) are initiated <br>When this key is set to 'ADOM2KPP' and chm\_messy\_jval\_l is set to .false., chm\_messy\_jval\_l is overwritten to .true. |
| chm\_pkg\_pm\_s | Name of PM species package<br>Options:<br>  - 'CAM2BINS'  (i.e., PMfine + PMcf)<br>  - 'CAM12BINS' <br>  - 'GOCART\_SO2' <br>  - 'NIL' | undefined | character | Set to 'CAM2BINS' in RAQDPS022. <br>When this key is set to 'NIL' or 'GOCART\_SO2', chm\_soa\_s, chm\_seaflux\_s, chm\_winddust\_s, chm\_hetchem\_s, chm\_pm\_drydep\_s, chm\_aqueous\_s and chm\_met\_modulation\_s keys are all set to 'NIL' too, and chm\_aqhi\_l, chm\_blcld\_l, chm\_diag\_wetdep\_l, chm\_direct\_l, chm\_indirect\_l and chm\_diag\_aerosols\_l are overwritten to .false. |
| chm\_strato\_s | Name of stratospheric O3 scheme<br>Options:<br>  - 'LINOZ' <br>  - 'NIL' | 'NIL' | character | Set to 'LINOZ' in RAQDPS022. |
| chm\_gas\_drydep\_s | Name of gas dry deposition scheme<br>Options: <br>  - 'ROBICHAUD' (previously 'ROBICHAUD' with rc\_case=1) <br>  - 'ROBICHAUD2' (previously ROBICHAUD with rc\_case=2) <br>  - 'ROBICHAUD3' (updated ROBICHAUD)  <br>  - 'NIL' | undefined | character | Set to 'ROBICHAUD' in RAQDPS022. <br>When this key is set to 'NIL', chm\_diag\_drydep\_l gets overwritten to .false. |
| chm\_o3icedep | Specify different ozone deposition rate over ice/glaciers<br>Options: <br>  - 1 (specifies ozone dry deposition velocity over ice to be 1.0E-4 m/s) <br>  - 0 (calculates ozone dry deposition velocity | 0 | integer | Set to 1 in RAQDPS022. |
| chm\_mar\_halo\_l | Enable additional ozone removal over the oceans through parameterized reactions with marine halogens | .false. | logical | Not present in RAQDPS022 <br> Not implemented in mach\_gas\_drydep\_solver2.ftn90. |
| chm\_messy\_jval\_l | Invoke MESSY for photolysis rates | .false. | logical | Set to default in RAQDPS022. |
| chm\_active\_ch4\_l | Activate an oxidation of CH4 by OH independently of the gas-phase mechanism | .false. | logical | Not present in RAQDPS022. |
| aerosize (13) | Aerosol bin boundaries depend on number of bins (i.e., vector of length 13): <br>  - 2 bins: 0.005, 1.280, 5.12 <br>  - 12 bins: 0.005, 0.010, 0.020, 0.040, 0.080, 0.160, 0.320, 0.640, 1.280, 2.560, 5.120, 10.240, 20.480 | not defined | real | Set to values for 2 bins in RAQDPS022. |
| chm\_intrsec\_ndiv | Number of subdivisions per bin in mach\_cam\_intrsec | -1 | integer | Set to 6 in RAQDPS022. |
| chm\_intrsec\_ver | Version number of the code in mach\_cam\_intrsec<br> Values should be integer numbers between 1 and 6. | -1 | integer | Set to 1 in RAQDPS022. |
| chm\_pm\_coag\_step\_intvl | Interval of time steps to calculate aerosol coagulation coefficient | -1 | integer | Set to 1 in RAQDPS022. |
| chm\_soa\_s | Name of SOA scheme<br>Options:<br>  - 'ODUM' <br>  - 'JIANG' <br>  - 'NIL' | undefined | character | Set to 'JIANG' in RAQDPS022. <br> `ODUM` only available for `chm_pkg_gas_s[1:5]=['SAPRC']` |
| chm\_pm\_drydep\_s | Name of aerosol dry deposition scheme<br>Options:<br>  - 'ZHANG' <br>  - 'ZHANG\_MAKAR' <br>  - 'NIL' | undefined | character | Set to 'ZHANG' in RAQDPS022. <br>See Mantis [#4217] (https://ulysse.cmc.ec.gc.ca/mantis/view.php?id=4217) for more information |
| chm\_blcld\_l | Activate below cloud processes | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_split\_snow\_rain\_l | Use GEM rain-snow partitioning<br>Options:<br>  - .true. use GEM's native rain-snow partitioning<br>  - .false. based on temperature, compute rain-snow partitioning in MACH | .false. | logical | Not present in RAQDPS022. |
| chm\_wdep\_scav\_coef\_s | Below-cloud particle scavenging scheme<br>Options:<br>  - 'SLINN' is original scheme<br>  - 'WANG2014' is new semi-empirical scheme | 'SLINN' | character | Not present  in RAQDPS022. |
| chm\_hetchem\_s | Name of heterogeneous chemistry scheme<br>Options:<br>  - 'HETV' <br>  - 'ISO' <br>  - 'NIL'  | undefined | character | Set to HETV in RAQDPS022.<br>See [Makar et al. (2003)] (https://hpfx.science.gc.ca/~mim001/AURAMS/AURAMS\_Papers+Reports/2003/Makar\_et\_al\_2003a.pdf) for more information. |
| chm\_hetchem\_metstlb\_l | Impose heterogeneous metastable option | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_aqueous\_s | Name of aqueous-phase chemistry scheme<br>Options: <br>  - 'GONG' <br>  - 'NIL' | undefined | character | Set to 'GONG' in RAQDPS022. |
| chm\_direct\_l | Enable direct radiative feedback. Not fully functional yet. Missing changes to the physics code. | .false. | logical | Not present in RAQDPS022. <br>When this key is set to .true.:<br>  a. chm\_diag\_aero\_opt is overwritten to .true. too; and<br>  b. if less than four (4) mid-band wavelengths (elements of aero\_opt\_wavel) are provided, the model execution will abort. |
| chm\_indirect\_l | Enable indirect (cloud) feedback. Not fully functional yet. Missing changes to the physics code. | .false. | logical | Not present in RAQDPS022. <br>When this key is set to .true., a condensation scheme in physics (stcond key) has to be set to one of available microphysics schemes. Otherwise, it will be overwritten to .false. |
| chm\_canopy\_shading\_l | Enable forest canopy shading | .false. | logical | Not present in RAQDPS022. |
| chm\_ammonia\_bidi\_s | Ammonia bidirectional flux<br>Options:<br>  - 'OFF' <br> - 'GEP' <br>  - 'GEP2D' | 'OFF' | character | Ground emissions potential values are a function of land-use category only for option 'GEP'. For option 'GEP2D', the ground emission potential is set via 2D surface fields read from FST files. **Note: If emissions potentials are used to model a source type that is already included in the input emissions files, then the emissions from this source type must be removed from the input emissions files to avoid double counting.** |
| chm\_ammonia\_gep | Ground emissions potential for ammonia bidirectional flux scheme | [1000, 1000, 1000, 500, 1000, 2000, 800, 0, 20, 100, 20, 0, 0, 0, 0] | real, dimension(15) | Ground emissions potential values as a function of land-use category only. Only used when chm\_ammonia\_bidi\_s = 'GEP' |
| chm\_vit\_l | Enable additional vertical diffusion due to the vehicle-induced turbulence when on-road area emissions are provided separate from other area emissions (i.e. chm\_get\_on\_emis\_l=.true.). | .false. | logical | Not present in RAQDPS022.<br>For this process to be active in the model, additional hourly input fields are required (i.e. on-road Vehicle Kilometers Travelled (VKT) for cars, mid-sized vehicles and trucks).<br>If this key is set to .true., but chm\_get\_on\_emis to .false., model will abort. |
| chm\_dep\_lai2d\_l | Enable the use of the satellite-derived LAI field in gas dry deposition and correct LAI of needleleaf trees to be constant throughout the year at its maximum value | .false. | logical |  |
| chm\_sat\_seasons\_l | Enable satellite-LAI-based adjustment of seasons in biogenic emissions modules | .false. | logical |  |


#### Emissions Treatment Settings

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| chm\_emis\_master | Master key to activate reading of all types of emissions (mobile, biogenic, major).  No emissions are read if this key is set to .false. | .true. | logical | Set to .true. in RAQDPS022.<br>When this key is set to .false., all chm\\_get\\_\*\\_emis\\_l keys are overwritten to .false. too. <br> When this key is set to .true., but all chm\\_get\\_\*\\_emis\\_l keys are set to .false., this key gets overwritten to .false. too. |
| em\_nesdt | Emissions input timestep (s).<br>If this key is set to a value that is a multiple of the chemistry time step (cf. chm\_step\_factor), then the same emissions values will be applied for all of the chemistry time steps until the next set of emissions are input | 0.0 | real | Set to 3600. in RAQDPS022 (i.e., hourly emissions). |
| chm\_get\_ae\_emis\_l | Import "area" emissions (AE: effectively all anthropogenic surface emissions, including on- and off-road, area, and minor point sources) | .false. | logical | Set to .true. in RAQDPS022. <br>May exclude fugitive dust emissions depending upon setting for chm\_get\_fd\_emis\_l. <br>May exclude on-road emissions depending upon setting for chm\_get\_on\_emis\_l. |
| chm\_get\_be\_emis\_l | Import biogenic emissions (BE) | .false. | logical | Set to .true. in RAQDPS022. <br>When this key is set to .false., chm\_biog\_s gets overwritten to 'NIL'. |
| chm\_get\_mj\_emis\_l | Import major (MJ) point sources emissions. <br>These are emissions sources with physical stack characteristics (i.e., height, diameter, exit temperature, exit speed) for which plume rise values are calculated (cf. 'chm\_mj\_treatment\_s') | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_get\_wf\_emis\_l | Import CFFEPS wildfire (WF) emissions. <br>When set to .true., wf\_case has to be set to 4. <br>If FEPS wildfire emissions are being used, set to .false. | .false. | logical | Set to .false. in RAQDPS022 and to .true. in FireWork022.<br>When this key is set to .true. and wf\_case is set to 0, the model aborts. <br>When this key is set to .true. and wf\_case is set to 1, 2 or 3, the model only issues a warning. <br>When this key is set to .false. and wf\_case is set to 4, the model aborts.|
| chm\_get\_on\_emis\_l | Import on-road area emissions separately from the other area emissions | .false. | logical | Not present in RAQDPS022.<br>If the input files contain separate variables for on-road area emissions, this key needs to be set to .true. to treat them. |
| chm\_get\_fd\_emis\_l | Import fugitive dust emissions separately from the other area emissions | .false. | logical | Not present in RAQDPS022.<br>If the input files contain separate variables for fugitive dust area emissions, this key needs to be set to .true. to treat them.<br>When this key is set to .false. and chm\_met\_modulation\_s to 'ON', model aborts.<br>When this key is set to .true. and chm\_met\_modulation\_s to 'CM\_ONLY', model aborts. |
| chm\_htap\_emis\_l | Run with the HTAP set of emissions | .false. | logical | Set to default in RAQDPS022. |
| chm\_mj\_treatment\_s | Name of major point plume-rise algorithm.<br> Options: .<br> - 'PLUMERISE' (Briggs' plume rise)<br> - 'PLUMERISE2' (Briggs' layered atmospheric profile plume rise)<br>  - 'NIL'<br>When both chm\_get\_mj\_emis\_l and chm\_get\_wf\_emis\_l are set to .false., this key gets overwritten to 'NIL'. | undefined | character | Set to 'PLUMERISE' in RAQDPS022. |
| chm\_ae\_spread\_l | Impose the injection of area emissions into the lowest 2 layers instead of only the lowest layer | .false. | logical | Set to .true. in RAQDPS022. <br>This key was introduced as a short-term measure for use with GEM5.0 and GEM-MACH3.0 and the switch from L80 to L84.|
| wf\_case | Type of plumerise for wildfire emissions<br><br>Options:<br>  - **0** (Brigg's) <br>  - **1** (Landuse-based, with Gaussian distribution for flaming portion, and smoldering portion evenly distributed within PBL) <br>  - **2** (Landuse- and PBL-stability-based: In stable conditions same as 1, while in unstable conditions the flaming plume is also evenly distributed within PBL) <br>  - **3** (Evenly distributed within PBL) <br>  - **4** (CFFEPS distribution) | 0 | integer | Set to 0 in RAQDPS022 and to 4 in FireWork022. <br><br> Options 0, 1, 2 and 3 are valid only with FEPS input and chm\_get\_wf\_emis\_l has to be set to .false. <br><br> Option 4 is valid only with CFFEPS input; chm\_get\_wf\_emis\_l has to be set to .true. <br><br> For all options, anthropogenic emissions are treated as defined in chm\_mj\_treatment\_s |
| chm\_cffeps\_online\_l | Enable online CFFEPS calculations | .false. | logical | Not present in RAQDPS022 <br> Enabling 'chm\_cffeps\_online\_l' turns off chm\_get\_wf\_emis\_l. |
| chm\_biog\_s | Name of biogenic emissions package<br>Options: <br>  - 'BEIS3.09' <br>  - 'NIL' | undefined | character | Set to 'BEIS3.09' in RAQDPS022. |
| chm\_winddust\_s | Name of wind dust package<br>Options: <br>  - 'CAM\_WINDDUST' <br> - 'NIL' | undefined | character | Not currently implemented.  Should be removed. |
| chm\_seaflux\_s | Name of sea salt flux scheme<br>Options: <br>  - 'GONG\_MONAHAN' <br>  - 'GONG\_MONAHAN\_F' <br>  - 'SMITH' <br>  - 'NIL' | undefined | character | Set to 'GONG\_MONAHAN\_F' in RAQDPS022. |
| chm\_nosurfzoneflux\_l | Exclude surf-zone sea-salt flux | .false. | logical | Not yet utilized in the operational code. |
| chm\_met\_modulation\_s | Meteorological modulation of fugitive dust emissions <br>Options: <br>  - 'ON' (fugitive dust emissions are meteorologically modulated) <br>  - 'NIL' or 'OFF' (fugitive dust emissions are used without modulation) <br>  - 'CM\_ONLY' (when fugitive dust emissions are incorporated into area emissions, apply meteorological modulations only to the crustal-material anthropogenic area emissions) | 'OFF' | character | Not present in RAQDPS022. |


#### Diagnostic Output Settings

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| gas\_ppb\_out\_list\_s | List of gas species to have concentrations output in ppbV<br>Options: <br>  - a comma-separated list of any number of 4-character diagnostic output names of gas concentrations<br>  - 'LIST' (an abbreviated way of writing default 4 gases 'O3  ', 'N2  ', 'NO  ', 'S2  ')<br>  - 'TOUT' or 'ALL ' (an abbreviation to output all gas species concentrations in ppbV) | 'LIST' | character | Not present in RAQDPS022. |
| chm\_diag\_aerosols\_l | Diagnostic output for aerosol components in ug/m3 | .false. | logical | Set to default in RAQDPS022. |
| chm\_aqhi\_l | Diagnostic output for AQHI index: AQ25 and AQ10 | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_ss\_ao\_l | Enable including sea salt to the aerosol-concentration output (e.g. AF and AC) | .false. | logical | Not present in RAQDPS022. |
| chm\_diag\_wetdep\_l | Diagnostic output for wet deposition | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_diag\_drydep\_l | Diagnostic output for dry deposition | .false. | logical | Set to .true. in RAQDPS022. |
| chm\_diag\_accum\_l | Diagnostic output for accumulators | .false. | logical | Set to default in RAQDPS022. |
| chm\_diag\_colum\_l | Diagnostic output for column values | .false. | logical | Set to default in RAQDPS022. |
| chm\_moyhr | Accumulation interval in hours | 24 | integer | Set to default in RAQDPS022. |
| chm\_acchr | Running mean interval in hours | 24 | integer | Set to default in RAQDPS022. |
| aero\_opt\_wavel | List of up to 15 mid-band wavelengths (in meters) for AOD calculation | -1.0 | real | Not present in RAQDPS022. |
| chm\_diag\_aero\_opt\_l | Diagnostic output for aerosol optical properties | .false. | logical | Not present in RAQDPS022. <br>When this key is set to .true., at least one value has to be provided for aero\_opt\_wavel. Otherwise, the model execution will abort. |

#### Debugging Settings

| Name | Description | Default Value | Type | Comment |
| :------ | :----- | :------: | :-----: | :----- |
| chm\_timings\_l | Enable timing tracing in chemistry | .false. | logical | Set to default in RAQDPS022. |
| chm\_debug\_trace\_l | Enable debug tracing in chemistry | .false. | logical | Set to default in RAQDPS022. |
| chm\_debug\_2d\_i | Number of debug variables for 2D fields (max=18) | -1 | integer | Set to 0 in RAQDPS022. |
| chm\_debug\_3d\_i | Number of debug variables for 3D fields (max=18) | -1 | integer | Set to 0 in RAQDPS022. |
| dbg\_lat | Geographic latitude of debug location within the model domain | -99.0 | real | Not present in RAQDPS022. <br>(dbg\_lat,dbg\_lon) geographic point is mapped to the appropriate model grid point and saved internally as (dbg\_id,dbg\_jd) integers with default values of (-1,-1). <br>See note below. |
| dbg\_lon | Geographic longitude of debug location within the model domain | -999.0 | real | Not present in RAQDPS022. <br>See the comment for dbg\_lat. <br>See note below. |
| dbg\_lev | Model level number for debug listings | -1 | integer | Not present in RAQDPS022. <br>See note below. |
| dbg\_step | Time step for debug listings | -1 | integer | Not present in RAQDPS022. <br>See note below. |
| dbg\_itr | Tracer index of CAM species for debug listings | -1 | integer | Not present in RAQDPS022. |

Note: Although accessible from any CHEM code, dbg\_ variables are currently utilized only in CAM. Setting non-default values for them in gem\_settings.nml will result in debug listings for particular PM processes.


