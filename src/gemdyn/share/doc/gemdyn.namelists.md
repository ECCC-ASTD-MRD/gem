### grid Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| grd_dx | (LU only) Mesh length (resolution) in x-direction (degrees) | 0. | real |
| grd_dy | (LU only) Mesh length (resolution) in y-direction (degrees) | 0. | real |
| grd_iref | Reference Point I | -1 | integer |
| grd_jref | Reference Point J | -1 | integer |
| grd_latr | Latitude on rotated grid of reference point (degrees) | 0. | real |
| grd_lonr | Longitude on rotated grid of reference point (degrees) | 180. | real |
| grd_maxcfl | Max Supported Courrant number;<br>Pilot area = Grd_maxcfl_fact * Grd_maxcfl+Grd_bsc_base+Grd_bsc_ext1 | 1 | integer |
| grd_ni | Number of points along NI | 0 | integer |
| grd_nj | Number of points along NJ | 0 | integer |
| grd_overlap | (GY only) Overlap extent along latitude axis for GY grid (degrees) | 0. | real |
| grd_typ_s | Type of grid described using 2 characters:<br>- "GY" : Global Yin-Yang<br>- "LU" : LAM    Uniform | '' | character(len=2) |
| grd_xlat1 | Geographic longitude of the center of the computational domain (degrees) | 0. | real |
| grd_xlat2 | Geographic latitude of a point on the equator of the computational domain<br>east of  Grd_xlon1,Grd_xlat1  (degrees) | 0. | real |
| grd_xlon1 | Geographic latitude of the center of the computational domain (degrees) | 180. | real |
| grd_xlon2 | Geographic longitude of a point on the equator of the computational domain<br>east of Grd_xlon1,Grd_xlat1  (degrees) | 270. | real |


### vert_layers Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| g_nk |  |  |  |
| hyb | array of model levels (pressure),  0.0 < HYB < 1.0 |  |  |
| hyb_first_height |  | 40. | real |
| hyb_flat | a level (in the units of hyb or hyb_H above which<br>the vertical coordinate becomes flat | -1. | real |
| hyb_h | array of model levels (height  ),  hyb_H > 0.0 |  |  |
| hyb_lid_height |  | 60000. | real |
| hyb_lin_depth |  | -1. | real |
| hyb_nkequal |  | -1 | integer |
| hyb_rcoef | pair of coefficients (min,max) to control the flattenning of the<br>vertical coordinate |  |  |


### adz_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| adz_bc_lam_flux | Type of Flux when Bermejo-Conde LAM<br>- 0 -> NONE<br>- 1 -> Aranami et al. (2015)<br>- 2 -> Zerroukat and Shipway (2017) (ZLF) | 1 | integer |
| adz_bc_min_max_l | True-> Clipping Min Max + Proportional Mass-Fixer after Bermejo-Conde | .true. | logical |
| adz_bc_pexp_list |  | 1.0 | real |
| adz_ilmc_min_max_l | True-> Clipping Min Max after ILMC | .true. | logical |
| adz_ilmc_sweep_max | Number of neighborhood zones in ILMC | 2 | integer |
| adz_itraj | Number of iterations for trajectories computation | 3 | integer |
| adz_k0t_sub | Top boundary in GY for an embedded LAM | 1 | integer |
| adz_pil_sub_e_g |  | -1 | integer |
| adz_pil_sub_n_g |  | -1 | integer |
| adz_pil_sub_s_g |  | -1 | integer |
| adz_pil_sub_w_g |  | -1 | integer |
| adz_slt_winds | Use surface layer winds for advection of lowest thermodynamic level | .false. | logical |
| adz_verbose | Activate printing of conservation diagnostics if /=0 | 0 | integer |


### bubble_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| bubble_dx |  | 10. | real |
| bubble_dz |  | 10. | real |
| bubble_gaus_l |  | .false. | logical |
| bubble_ictr |  | -1 | integer |
| bubble_kctr |  | -1 | integer |
| bubble_ni |  | 101 | integer |
| bubble_nj |  | 1 | integer |
| bubble_nk |  | 100 | integer |
| bubble_rad |  | 25 | integer |
| bubble_theta |  | 303.16 | real |


### dyn_fisl Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| cstv_ba_8 | SL off-centering parameter for hydrostatic | 0.6 | real(kind=REAL64) |
| cstv_ba_m_8 | SL off-centering parameter for the momentum equations | 0.6 | real(kind=REAL64) |
| cstv_ba_nh_8 | SL off-centering parameter for nonhydrostatic | 0.5 | real(kind=REAL64) |
| cstv_psadj_8 | Fraction of adjustment to be given to the ocean | 1.d0 | real(kind=REAL64) |
| cstv_tstr_8 | T* basic state temperature (K) | 240.0 | real(kind=REAL64) |
| schm_advec | * 0   ->          NO advection<br>- 1   -> traditional advection<br>- 2   -> consistent advection with respect to off-centering<br>- 3   -> reversed advection with respect to off-centering | 1 | integer |
| schm_bitpattern_l | True-> Modify slightly code behaviour to ensure bitpattern<br>reproduction in restart mode using FST file | .false. | logical |
| schm_dry_mixing_ratio_l | True-> Tracers are mixing ratios with respect to dry air mass | .false. | logical |
| schm_hzdadw_l | True-> horizontal diffusion of momentum at each CN iteration | .false. | logical |
| schm_itcn | Number of iterations for Crank-Nicholson | 2 | integer |
| schm_itnlh | Number of iterations to solve non-linear Helmholtz problem | 2 | integer |
| schm_nblendyy | *  -1: no blending between Yin and Yang<br>-   0: blending at init only<br>- > 0: blending at every nblendyy timestep | -1 | integer |
| schm_orols_fromgeophy_l | True-> force the reading of MELS from geophysical file | .false. | logical |
| schm_orols_np | True-> Number of pass in the filter | 20 | integer |
| schm_orols_rc | True-> will retain 50% of Schm_orols_rc deltax waves | 10. | real |
| schm_phycpl_s | Physics coupling strategy | 'split' | character(len=16) |
| schm_psadj | * 0 -> No conservation of surface pressure<br>- 1 -> Conservation of Total air mass Pressure<br>- 2 -> Conservation of Dry air mass Pressure | 0 | integer |
| schm_psadj_print_l | True-> print dry/wet air masses | .false. | logical |
| schm_topo_l | True-> to use topography | .true. | logical |
| schm_wload_l | Apply water loading in the calculations | .false. | logical |


### dyn_kernel Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| dynamics_hydro_l | * True-> hydrostatic<br>- False-> non-hydrostatic | .false. | logical |
| dynamics_kernel_s | Main selector for dynamical kernel<br>- "DYNAMICS_FISL_P" : Fully implicit SL in pressure<br>- "DYNAMICS_FISL_H" : Fully implicit SL in height | 'DYNAMICS_FISL_P' | character(len=32) |
| schm_autobar_l | True-> auto barotropic option | .false. | logical |


### ensembles Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| ens_conf | Switch to activate generation of Markov chains, use of SKEB<br>and use of PTP | .false. | logical |
| ens_mc_seed | Seed of the random number generator usually we put DAY and member number<br>(3D MARKOV CHAINES) | -1 | integer |
| ens_ptp_cape | CAPE value in Kain-Fritsch scheme to stop perturbing the physical<br>tendencies<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_conf | switch to activate PTP (perturb tendencies of physics)<br>(2D MARKOV CHAINES) | .false. | logical |
| ens_ptp_crit_w | vertical velocity value (m/s) above which we stop perturbing the<br>physical tendencies<br>(2D MARKOV CHAINES) | 100.0 | real |
| ens_ptp_env_b | bottom value of transition zone of vertical envelope in sigma for PTP<br>(below that no perturbation)<br>(2D MARKOV CHAINES) | 1.0 | real |
| ens_ptp_env_u | upper value of transition zone of vertical envelope in sigma for PTP<br>(above that full perturbation)<br>(2D MARKOV CHAINES) | 1.0 | real |
| ens_ptp_fac_reduc | factor of reduction of the perturbation the physical tendencies (in PTP)<br>when convection occurs<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_l | (ignored) Ens_ptp_l = Ens_ptp_trnh-Ens_ptp_trnl+1<br>(2D MARKOV CHAINES) | 0 | integer |
| ens_ptp_latmax | (ignored) Ens_ptp_latmax = maxval(Ens_ptp_nlat)<br>(2D MARKOV CHAINES) | 0 | integer |
| ens_ptp_lmax | (ignored) Ens_ptp_lmax = maxval(Ens_ptp_l)<br>(2D MARKOV CHAINES) | 0 | integer |
| ens_ptp_m | (ignored) Ens_ptp_m = Ens_ptp_trnh+1<br>(2D MARKOV CHAINES) | 0 | integer |
| ens_ptp_max | maximum value of the 2D Markov chains<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_min | minimum value of the 2D Markov chain<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_mmax | (ignored) Ens_ptp_mmax = maxval(Ens_ptp_m)<br>(2D MARKOV CHAINES) | 0 | integer |
| ens_ptp_ncha | number of 2d Markov chains<br>(2D MARKOV CHAINES) | 0 | integer |
| ens_ptp_nlat | no. of latitudes for 2D Markov chains<br>(2D MARKOV CHAINES) | 8 | integer |
| ens_ptp_nlon | no. of longitudes for 2D Markov chains<br>(2D MARKOV CHAINES) | 16 | integer |
| ens_ptp_std | standard deviation value for 2D Markov chains<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_str | value of stretch for Markov chains<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_tau | decorrelation time (seconds) for 2D Markov chains<br>(2D MARKOV CHAINES) | 0.0 | real |
| ens_ptp_trnh | high wave number horizontal truncation limit for 2D Markov chains<br>(2D MARKOV CHAINES) | 8 | integer |
| ens_ptp_trnl | low wave number horizontal truncation limit for 2D Markov chains<br>(2D MARKOV CHAINES) | 1 | integer |
| ens_recycle_mc | Switch to activate recycling  of Markov chains parameters<br>in PTP & SKEB | .false. | logical |
| ens_skeb_alph | coefficient Alpha for momentum in SKEB<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_alpt | coefficient Alpha for temperature in SKEB<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_bfc | coefficient for Gaussian filter used in SKEB<br>(3D MARKOV CHAINES) | 1.0e-01 | real |
| ens_skeb_conf | Switch to activate SKEB<br>(3D MARKOV CHAINES) | .false. | logical |
| ens_skeb_dif | switch to do SKEB calculation based on diffusion<br>(3D MARKOV CHAINES) | .false. | logical |
| ens_skeb_div | switch to do the calculation of the divergence due to SKEB forcing<br>(3D MARKOV CHAINES) | .false. | logical |
| ens_skeb_gwd | switch to do SKEB calculation based on gravity wave drag<br>(3D MARKOV CHAINES) | .false. | logical |
| ens_skeb_lam | wavelength for Gaussian filter in SKEB<br>(3D MARKOV CHAINES) | 2.0e+05 | real |
| ens_skeb_max | maximum value of the 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_min | minimum value of the 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_ncha | <br>(3D MARKOV CHAINES) | 1 | integer |
| ens_skeb_nlat | number of latitudes of the gaussian grid used  for the 3D Markov chains<br>(used in the SKEB calculation)<br>(3D MARKOV CHAINES) | 8 | integer |
| ens_skeb_nlon | number of longitudes of the gaussian grid used for the 3D Markov chains<br>(in the SKEB calculation)<br>(3D MARKOV CHAINES) | 16 | integer |
| ens_skeb_std | std. dev. value for the 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_str | value of stretch for 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_tau | decorrelation time (seconds) for 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 0. | real |
| ens_skeb_trnh | high wave number truncation limit used in 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 8 | integer |
| ens_skeb_trnl | low wave number truncation limit used in 3D Markov chain (used by SKEB)<br>(3D MARKOV CHAINES) | 2 | integer |
| ens_spp_rhsint_lat | Set cutoff latitude (degrees from equator) for adv_rhsint SPP perturbations | 45. | real |
| ens_spplist | Activate SPP for specified parameters |  |  |
| ens_stat | switch to print global stat related to Markov chains, SKEB and PTP<br>(3D MARKOV CHAINES) | .false. | logical |


### gem_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| g_halox | number of points for the halo on X | 4 | integer |
| g_haloy | number of points for the halo on Y | 4 | integer |
| heap_nk | Heap memory will be painted to NaN using an array wrk01(G_ni,G_nj,Heap_nk) | -1 | integer |
| lctl_cktimeleft_l | True->to check for time left in job | .false. | logical |
| lctl_debug_l | True->to print more information to std output | .false. | logical |
| lctl_rxstat_s | precision in print glbstats<br>- 'LCL_4'<br>- 'GLB_8' | 'LCL_4' | character(len=6) |
| stat_liste | list of variables to do blocstat.<br>Any gmm variable name, or predefine lists :<br>- 'ALL'<br>- 'ALL_DYN_T0'<br>- 'ALL_DYN_T1'<br>- 'ALL_TR_T0'<br>- 'ALL_TR_T1' | ' ' | character(len=32) |
| tr3d_anydate_l | True-> tracers validity time does not have to match analysis | .false. | logical |
| tr3d_default_s | Override for default tracers attributes | ' ' | character(len=512) |
| tr3d_list_s | List of tracers and parameters to be read from analyse<br>- hzd = Horizontal Diffusion if 1<br>- wload = Water loading if 1<br>- min = Minimum admissible value<br>- max = Maximum admissible value<br>- mono = Type of monotonicity<br>-* 0 -> None<br>-* 1 -> Clipping QMSL<br>-* 2 -> ILMC<br>- mass = Type of mass conservation<br>-* 0 -> None<br>-* 1 -> Bermejo-Conde LEGACY=T<br>-* 100+weight*10+pexp_n -> Bermejo-Conde LEGACY=F<br>- intp  = Type of interpolation in advection (NOT case sensitive)<br>-* none -> No advection<br>-* tricub -> 3D Cubic Lagrange<br>-* bicubh_qv -> Horizontal Cubic + Vertical Quintic<br>List of sub-parameters<br>- weight = Weight in Bermejo-Conde<br>-* 1 -> Additive<br>-* 2 -> Multiplicative<br>-* 3 -> Additive+Factor pr_k/pr_s<br>- pexp_n = Rank in Adv_BC_pexp_list(PEXP_LIST_MAX) to choose P exponent REAL |  |  |
| vtopo_length_s | On which length of time to evolve topography | '' | character(len=16) |
| vtopo_start_s | Time at which to start evolving topography toward target | '' | character(len=16) |


### grdc Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| grdc_dx | x horizontal resolution of target cascade grid (degrees) | -1. | real |
| grdc_dy | y horizontal resolution of target cascade grid (degrees) | -1. | real |
| grdc_end_s | Time string (units D, H, M or S) from the start of the run to<br>stop producing the cascade files | ' ' | character(len=15) |
| grdc_fullgrid_l | TRUE to dump out permanent bus for cascade mode | .false. | logical |
| grdc_iref | Reference Point I on rotated cascade grid, 1 < Grdc_iref < Grdc_ni | -1 | integer |
| grdc_jref | Reference Point J on rotated cascade grid, 1 < Grdc_jref < Grdc_nj | -1 | integer |
| grdc_latr | Latitude on rotated grid of ref point, Grdc_iref,Grdc_jref (degrees) | 0. | real |
| grdc_lonr | Longitude on rotated grid of ref point, Grdc_iref,Grdc_jref (degrees) | 180. | real |
| grdc_maxcfl | Max Supported Courrant number;<br>Pilot area=Grdc_maxcfl +Grdc_bsc_base+Grdc_bsc_ext1 | 1 | integer |
| grdc_nbits | Number of bits for the packing factor |  |  |
| grdc_nfe | Nesting interval specified with digits ending with one character<br>for the units:<br>- S : seconds<br>- D : days<br>- M : minutes<br>- H : hours | ' ' | character(len=15) |
| grdc_ni | Number of points along X | 0 | integer |
| grdc_nj | Number of points along Y | 0 | integer |
| grdc_start_s | Time string (units D, H, M or S) from the start of the run to<br>start producing the cascade files | ' ' | character(len=15) |
| grdc_trnm_s | List of tracers to be written from piloting run | '@#$%' | character(len=4) |


### hvdif Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| eq_ramp_l | Latitudinal ramping of equatorial sponge | .false. | logical |
| eq_sponge | Coefficients that multiply KM to simulate sponge layer near the top<br>of the model. Warning! if this parameter is used, the EPONGE in the<br>physics namelist should be removed. |  |  |
| hzd_lnr | Background 2 delta-x removal ratio - range(0.0-1.0) | -1. | real |
| hzd_lnr_theta | Theta 2 delta-x removal ratio - range(0.0-1.0). | -1. | real |
| hzd_lnr_tr | Tracers 2 delta-x removal ratio - range(0.0-1.0) | -1. | real |
| hzd_pwr | Order of the background diffusion operator<br>2, 4, 6, 8 | -1 | integer |
| hzd_pwr_theta | Order of the background diffusion operator on theta<br>2, 4, 6, 8 | -1 | integer |
| hzd_pwr_tr | Order of the background diffusion operator on tracers | -1 | integer |
| hzd_smago_fric_heat | Frictional heating is considered when Hzd_smago_fric_heat>0. | 0. | real |
| hzd_smago_lev |  |  |  |
| hzd_smago_lnr | Coefficient of background diffusion added to the coefficient computed using the<br>Smagorinsky approach. Two options are available.<br>(i)  Fixed background diffusion: Just assign a single value, e.g., 0.2, and the resultant<br>diffusion will be similar to standard del-2 diffusion with Hzd_lnr=0.2.<br>(ii) Vertically variable diffusion: Takes three values. The first element of the array determines<br>the constant value of background diffusion coeff. below Hzd_smago_lev(1). The second<br>the value at Hzd_smago_lev(2). The third element determines the maximum coefficient<br>element represents at the model top. Two ramps of COS^2-type are used between<br>Hzd_smago_lev(1) and Hzd_smago_lev (2), and between Hzd_smago_lev(2) and the model lid. |  |  |
| hzd_smago_param | Main Smagorinsky control parameter (usual range 0.1-0.3) | -1. | real |
| hzd_smago_prandtl | Apply Smago diffusion on theta using Hzd_smago_param/Hzd_smago_prandtl parameter | -1. | real |
| hzd_smago_prandtl_hu | Apply Smago diffusion on HU using Hzd_smago_param/Hzd_smago_prandtl_hu parameter | -1. | real |
| hzd_smago_theta_base_l | If TRUE then background diffusion is applied to THETA and HU. | .false. | logical |
| p_lmvd_high_lat | Latitude at which the multiplication factor becomes P_lmvd_weigh_high_lat | 30.0 | real |
| p_lmvd_low_lat | latitude at which the multiplication factor becomes P_lmvd_weigh_low_lat | 5.0 | real |
| p_lmvd_weigh_high_lat | Multiplication factor of P_pbl_spng at latitude P_lmvd_high_lat | 1.0 | real |
| p_lmvd_weigh_low_lat | Multiplication factor of P_pbl_spng at latitude P_lmvd_low_lat | 1.0 | real |
| vspng_coeftop | Top coefficient for del-2 diffusion (m2/s) | -1. | real |
| vspng_nk | Number of levels from the top of the model | 0 | integer |
| vspng_riley_l | True-> Riley diffusion on vertical motion on Vspng_nk levels | .false. | logical |


### init Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| iau_cutoff | Filter cutoff period for Iau_weight_S='sin' in hours | 6. | real |
| iau_interval | The number of seconds between increment fields | -1. | real |
| iau_ninblocx | IAU Input PE blocking along npex | 1 | integer |
| iau_ninblocy | IAU Input PE blocking along npey | 1 | integer |
| iau_period | The number of seconds over which IAU will be  will be run<br>(typically the length of the assimilation window).<br>Default < 0 means that no IAUs are applied. | -1. | real |
| iau_stats_l | IAU Input Stats | .false. | logical |
| iau_tracers_s | An optional list of tracers to be incremented. |  |  |
| iau_weight_s | The type of weighting function to be applied to the analysis increments:<br>- 'constant' (default) uniform increments<br>- 'sin' DF-style weights (Fillion et al. 1995) | 'constant' | character(len=64) |
| init_balgm_l | true -> Digital filter initialization is performed | .false. | logical |
| init_dflength_s | number of points for digital filter (equals the number of timesteps +1) | '5p' | character(len=16) |
| init_dfpl_s | period limit of digital filter units D,H,M,S | '6h' | character(len=16) |
| init_dftr_l | * true -> passive tracers digitally filtered<br>- false-> passive tracers set to result obtained at mid-period during initialization period (no filtering) | .false. | logical |
| init_dfwin_l | true -> Windowing is applied | .true. | logical |
| zdot_divhlm_l | True-> divergence high level modulation in initial computation of Zdot | .false. | logical |


### inp Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| inp_blacklist_s | List of variables to NOT process during input |  |  |
| inp_npes | Number of PEs to use for input | 1 | integer |
| inp_vertintype_tracers_s | Type of vertical interpolation scheme | 'cubic' | character(len=8) |
| perturb_nbits | Number of bits to perturb on initial conditions | 0 | integer |
| perturb_npts | Stride for perturbation on initial conditions | 10 | integer |


### lam Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| lam_0ptend_l | True-> for blending to zero the physics tendency in blending area | .true. | logical |
| lam_blend_h | Number of points for horizontal blending | 10 | integer |
| lam_blend_t | Number of levels for top blending | 0 | integer |
| lam_blendoro_l | True-> to blend the model topography with the pilot topography | .true. | logical |
| lam_ctebcs_l | True-> to force constant (fixed) boundary conditions | .false. | logical |
| lam_gbpil_t | Number of levels for top piloting | -1 | integer |
| lam_hint_s | Type of horizontal interpolation to model grid<br>- 'CUB_LAG'<br>- 'LINEAR'<br>- 'NEAREST' | 'CUB_LAG' | character(len=16) |
| lam_toptt_l | True-> The plane of the top temperature layer is completely<br>overwritten from the 2D pilot data | .false. | logical |


### mtn_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| mtn_dx |  | 500. | real |
| mtn_dz |  | 300. | real |
| mtn_flo |  | 10. | real |
| mtn_hght |  | 250. | real |
| mtn_hwx |  | 10. | real |
| mtn_hwx1 |  | 8. | real |
| mtn_ni |  | 401 | integer |
| mtn_nj |  | 1 | integer |
| mtn_nk |  | 65 | integer |
| mtn_nstar |  | 0.01 | real |
| mtn_pos_seed |  |  |  |
| mtn_rcoef |  |  |  |
| mtn_tzero |  | 303.15 | real(kind=REAL64) |
| mtn_wind_seed |  | 0. | real |
| mtn_zblen_thk |  | 0. | real |


### ops_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| ops_configuration_s | Main selector for Operational configurations combo defaults | '' | character(len=32) |


### out Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| out3_cliph_l | True-> to clip humidity variables on output | .false. | logical |
| out3_close_interval_s | Interval of output file name change | '1h' | character(len=16) |
| out3_etik_s | 'etiket' used for output fields | 'GEMDM' | character(len=12) |
| out3_ip3 | Default value for IP3 is 0, -1 for IP3 to contain step number,<br>>0 for given IP3 | 0 | integer |
| out3_lieb_conv | Precision criteria for the Liebman procedure | 0.1 | real |
| out3_lieb_levels | List of levels for underground extrapolation |  |  |
| out3_lieb_maxite | Maximum number of iterations for the Liebman procedure | 100 | integer |
| out3_liebxch_iter | number of iterations to exchange halo for the Liebman procedure | 4 | integer |
| out3_linbot | Number of layers close to the bottom of the model within which a<br>linear interpolation of GZ will be performed | 0 | integer |
| out3_nbitg | Packing factor used for all variables except for those defined in<br>Out_xnbits_s | 16 | integer |
| out3_ndigits | Minimum of digits used to represent output units | 3 | integer |
| out3_npes | Total number of PEs for output using MFV collector | 1 | integer |
| out3_npex | Total number of PEs along npex for output using MID collector | -1 | integer |
| out3_npey | Total number of PEs along npey for output using MID collector | -1 | integer |
| out3_postproc_fact | Sortie jobs lauched every Out3_postproc_fact*Out3_close_interval_S | 6 | integer |
| out3_vinterp_type_s | Vertical interpolation scheme for output | 'linear' | character(len=12) |


### sol Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| sol_fgm_eps | Relative tolerance for convergence of the elliptic iterative solver,<br>norm(residual) < tol*norm(b) | 1.d-07 | real(kind=REAL64) |
| sol_fgm_maxits |  | 200 | integer |
| sol_gauss_niter | 3D preconditioner for iterative solver | 2 | integer |
| sol_im | size of Krylov subspace in iterative solver - should not exceed 100 | 15 | integer |
| sol_krylov3d_s | Krylov method for 3d iterative solver (FGMRES or FBICGSTAB) | 'FGMRES' | character(len=26) |
| sol_one_transpose_l | True => use the one-transpose solver | .false. | logical |
| sol_ovlpx | 3D preconditioner for iterative solver | 4 | integer |
| sol_ovlpy | 3D preconditioner for iterative solver | 4 | integer |
| sol_precond3d_s | 3D preconditioner for iterative solver | 'RAS' | character(len=26) |
| sol_type_s | Type of solver<br>- 'ITERATIF'<br>- 'DIRECT' | 'ITERATIVE_3D' | character(len=26) |
| sol_yyg_eps | Absolute tolerance for convergence of the classical Schwarz algorithm<br>in Yin Yang configurations, norm(residual) < tol | 1.d-04 | real(kind=REAL64) |
| sol_yyg_maxits | maximum number of iterations allowed for the Yin-Yang iterative solver | 40 | integer |


### spn Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| spn_cutoff_scale_large | The filter will be set to 0. for smaller scales (in km) | 300. | real |
| spn_cutoff_scale_small | The filter will be set to 1.0 for larger scales (in km)<br>Transition between Spn_cutoff_scale_small and Spn_cutoff_scale_large<br>will have a COS2 shape. | 100. | real |
| spn_freq | Nudging interval - in sec<br>Nudging is performed every every Spn_freq sec | -1 | integer |
| spn_relax_hours | Nudging relaxation timescale - in hours | 10. | real |
| spn_start_lev | Nudging profile lower end in hyb level (eg. 1.0 or 0.8)<br>profile will be set to 0. when hyb > 0.8 | 1.0 | real |
| spn_trans_shape_s | Nudging profile transition shape('COS2' or 'LINEAR')<br>Set the shape between Spn_start_lev and Spn_up_const_lev | 'LINEAR' | character(len=16) |
| spn_up_const_lev | Nudging profile upper end in hyb level (eg. 0.0 or 0.2)<br>profile wll be set to 1.0 when hyb < 0.2 | 0.0 | real |
| spn_weight_l | Nudging weight in temporal space (.true. or .false.).<br>If the driving fields are available every 6 hours and Spn_freq is<br>set to 30 minutes then nudging will have more weight every six hours<br>when the driving fields are available | .false. | logical |
| spn_wt_pwr | The weight factor when Spn_weight_L=.true.<br>(The weigh factor is COS2**(Spn_wt_pwr), Spn_wt_pwr could  be set as<br>0, 2, 4, 6. If Spn_wt_pwr = 2, weight factor is COS2) | 2 | integer |


### step Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| fcst_bkup_additional_s | Save a restart file + continue at that time | 'NIL' | character(len=16) |
| fcst_bkup_s | Save a restart file + continue every Fcst_bkup_S | 'NIL' | character(len=16) |
| fcst_end_s | End date for model run slice (yyyymmdd.hhmmss) | ' ' | character(len=16) |
| fcst_gstat_s | Output global stat (glbstat) every Fcst_gstat_S | ' ' | character(len=16) |
| fcst_nesdt_s | Read nesting data every Fcst_nesdt_S | ' ' | character(len=16) |
| fcst_rstrt_s | Save a restart file + stop every Fcst_rstrt_S | 'NIL' | character(len=16) |
| fcst_spinphy_s |  | ' ' | character(len=16) |
| fcst_start_s | Starting date for model run slice (yyyymmdd.hhmmss) | ' ' | character(len=16) |
| step_alarm | Setting for Fortran alarm time | 600 | integer |
| step_dt | Length of model timestep (sec) | -1. | real(kind=REAL64) |
| step_leapyears_l | Account for leap years | .true. | logical |
| step_runstrt_s | Starting date for model run  (yyyymmdd.hhmmss) | 'NIL' | character(len=16) |


### theo_cfgs Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| theo_case_s |  | 'NONE' | character(len=15) |


### dcmip Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| dcmip_case | Dcmip case selector<br>- Dcmip_case=0 (NONE)<br>- Dcmip_case=11 (3D deformational flow)<br>- Dcmip_case=12 (3D Hadley-like meridional circulation)<br>- Dcmip_case=13 (2D solid-body rotation of thin cloud-like tracer in the presence of orography)<br>- Dcmip_case=20 (Steady-state at rest in presence of oro.)<br>- Dcmip_case=21 (Mountain waves over a Schaer-type mountain)<br>- Dcmip_case=22 (As 21 but with wind shear)<br>- Dcmip_case=31 (Gravity wave along the equator)<br>- Dcmip_case=41X (Dry Baroclinic Instability Small Planet)<br>- Dcmip_case=43 (Moist Baroclinic Instability Simple physics)<br>- Dcmip_case=161 (Baroclinic wave with Toy Terminal Chemistry)<br>- Dcmip_case=162 (Tropical cyclone)<br>- Dcmip_case=163 (Supercell Small Planet) | 0 | integer |
| dcmip_lower_value | Set lower value of Tracer in Terminator<br>- Dcmip_lower_value=0 (free)<br>- Dcmip_lower_value=1 (0)<br>- Dcmip_lower_value=2 (1.0e-15) | 0 | integer |
| dcmip_moist | Account for moisture<br>- Dcmip_moist=0 (dry)<br>- Dcmip_moist=1 (moist) | 1 | integer |
| dcmip_nuz_th | Vertical Diffusion Theta (if <0,we remove REF) | 0. | real |
| dcmip_nuz_tr | Vertical Diffusion Tracers (if <0,we remove REF) | 0. | real |
| dcmip_nuz_wd | Vertical Diffusion Winds (if <0,we remove REF) | 0. | real |
| dcmip_pbl_type | Type of planetary boundary layer<br>- Dcmip_pbl_type=-1 (none)<br>- Dcmip_pbl_type=0 (Reed-Jablonowski Boundary layer)<br>- Dcmip_pbl_type=1 (Georges Bryan Boundary Layer) | -1 | integer |
| dcmip_prec_type | Type of precipitation/microphysics<br>- Dcmip_prec_type=-1 (none)<br>- Dcmip_prec_type=0 (Kessler Microphysics)<br>- Dcmip_prec_type=1 (Reed-Jablonowski Large-scale precipitation) | -1 | integer |
| dcmip_terminator_l | Do Terminator chemistry if T | .false. | logical |
| dcmip_x | Earth's radius reduction factor | 1.d0 | real(kind=REAL64) |


### williamson Namelist

| Name          | Description            |  Default Value | Type |
| ------------- | ---------------------- | -------------- | ---- |
| williamson_alpha | Rotation angle in DEGREE (W_NAIR=0) | 0. | real |
| williamson_case | Williamson case selector<br>- Williamson_case=0 (none)<br>- Williamson_case=1 (Advection 2D - Williamson_NAIR/Terminator_L)<br>- Williamson_case=2 (Steady-state nonlinear zonal geostrophic flow - Williamson et al.,1992,JCP,102,211-224)<br>- Williamson_case=5 (Zonal flow over an isolated mountain - Williamson et al.,1992,JCP,102,211-224)<br>- Williamson_case=6 (Rossby-Haurwitz wave - Williamson et al.,1992,JCP,102,211-224)<br>- Williamson_case=7 (The 21 December 1978 case - Williamson et al.,1992,JCP,102,211-224)<br>- Williamson_case=8 (Galewsky's barotropic wave - Galewsky et al.,2004,Tellus,56A,429-440) | 0 | integer |
| williamson_clat0 | LAT cosine Bell in DEGREE (W_NAIR=0) | 0. | real |
| williamson_clon0 | LON cosine Bell in DEGREE (W_NAIR=0) | 270. ! 3*pi/2 rad in the paper | real |
| williamson_k,williamson_n,williamson_amp_8,williamson_wave_type |  | 'ROSSBY' | character(len=6) |
| williamson_lower_value | Set lower value of Tracer in Terminator<br>- Williamson_lower_value=0 (free)<br>- Williamson_lower_value=1 (0)<br>- Williamson_lower_value=2 (1.0e-15) | 0 | integer |
| williamson_nair | Used when Williamson_case=1<br>- Williamson_NAIR=0 (Solid body rotation of a cosine bell - Williamson et al.,1992,JCP,102,211-224)<br>- Williamson_NAIR=1 (Deformational Non-divergent winds - Lauritzen et al.,2012,GMD,5,887-901)<br>- Williamson_NAIR=2 (Deformational divergent winds - Lauritzen et al.,2012,GMD,5,887-901)<br>- Williamson_NAIR=3 (Deformational Flow for Circular vortex - Nair and Machenhauer,2002,MWR,130,649-667) | 0 | integer |
| williamson_period | rotation period in SECS (W_NAIR=0) | 12.*24.*3600. | real |
| williamson_phi0 | Phi0 of cosine Bell (W_NAIR=0) | 1000. | real |
| williamson_radius | Scaling for radius of cosine Bell (W_NAIR=0) | 3. | real |
| williamson_rlat0 | LAT rot.POLE in DEGREE (W_NAIR=3) | 0. | real |
| williamson_rlon0 | LON rot.POLE in DEGREE (W_NAIR=0/3) | 0. | real |
| williamson_terminator_l | Do Terminator chemistry if T | .false. | logical |


