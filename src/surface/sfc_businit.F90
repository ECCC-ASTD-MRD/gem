!-------------------------------------- LICENCE BEGIN --------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ----------------------------

!/@*
subroutine sfc_businit(moyhr,ni,nk)
   use sfc_options
   use sfcbus_mod
   use svs_configs
   implicit none
!!!#include <arch_specific.hf>
   !@Object Establishes requirements in terms of variables in the 4 main buses
   !        (busent, busdyn, busper and busvol) for the TEB (Urban scheme).
   !@Arguments
   integer, intent(in) ::  moyhr,ni,nk !# horiz and vert dimensions
   !@Author M. Desgagne (Oct 1995)

   !@Revision
   ! 001      L. Spacek  (Aug 2010) - Complete rewrite
   ! 002      B. Dugas   (Oct 2010) - A few small corrections
   ! 003      L. Spacek  (Sep 2011) - Eliminate obsolete convection options
   ! 004      M. Abrahamowicz (May 2016) - Add SVS
   ! 005      M. MacKay       (Oct 2018/Sep 2022) - Add CSLM
   !*@/

#include "phymkptr.hf"

   integer :: nmos
   integer :: alb_road, alb_roaden, alb_roof, alb_roofen, alb_wall, &
        alb_wallen, azim, bld, blden, bld_height, bld_heighten, can_hw_ratio, &
        d_road, d_roaden, d_roof, d_roofen, d_wall, d_wallen, emis_road, &
        emis_roaden, emis_roof, emis_roofen, emis_wall, emis_wallen, fvapliq, &
        g_road, &
        g_roof, g_town, g_wall, h_industry, h_industryen, h_road, h_roof, &
        h_town, h_traffic, h_trafficen, h_wall, hc_road, hc_roaden, hc_roof, &
        hc_roofen, hc_wall, hc_wallen, le_industry, le_industryen, le_road, &
        le_roof, le_town, le_traffic, le_trafficen, le_wall, nat, naten, pav, &
        paven, q_canyon, q_canyonen, rn_road, rn_roof, rn_town, rn_wall, &
        sroad_alb, sroad_alben, sroad_emis, sroad_emisen, sroad_rho, &
        sroad_rhoen, sroad_scheme, sroad_t, sroad_ten, sroad_ts, sroad_tsen, &
        sroad_wsnow, sroad_wsnowen, sroof_alb, sroof_alben, sroof_emis, &
        sroof_emisen, sroof_rho, sroof_rhoen, sroof_scheme, sroof_t, &
        sroof_ten, sroof_ts, sroof_tsen, sroof_wsnow, sroof_wsnowen, &
        svf_road, svf_wall, t_canyon, t_canyonen, t_road, t_roaden, t_roof, &
        t_roofen, t_wall, t_wallen, tc_road, tc_roaden, tc_roof, tc_roofen, &
        tc_wall, tc_wallen, ti_bld, ti_blden, ti_road, ti_roaden, tsun, &
        u_canyon, wall_o_hor, wall_o_horen, ws_road, ws_roaden, ws_roof, &
        ws_roofen, &
        yradin, yradrfsun, yradrfshade, yutciin, yutcirfsun, &
        yutcirfshade, ytrfzt, ytrdzt, yurdzu, ywbgtrfsun, ywbgtrfshade, &
        yutcicin, yutcicsun, yutcicshade, yutcicrfsun, yutcicrfshade, &
        ytglbrfsun, ytglbrfshade, ytwetbrf, yq8, yq9, yq10, yq11, yq12, &
        yq13, &
        z0_road, z0_roaden, z0_roof, z0_roofen, z0_town, &
        z0_townen, zenith, emtw, alscatw, tsradtw
   integer :: acoef, alveg, bcoef, c1sat, c2ref, c3ref, clay, cveg, &
        eflux, emsvc, gamveg, husurf, hv, iceline, lai, melts,  &
        meltsr, pcoef, psn, psng, psnv, resa, rgl, rnet_s, rst, sand, &
        snoagen, snoalen, snoma, snoro, stomr, tsoil, vegdati, vegf, &
        vegfrac, vegf_evol, wfc, wsat, wsnow, wveg, wwilt
   integer :: cgsat, dsst, dtdiag, glacier, glsea0, &
        icedp, skin_depth, skin_inc, snoal, snoden, &
        snodp, tglacier, tmice, tnolim, &
        twater, urban, &
        yradsun, yradshade, yutcisun, yutcishade, &
        ywbgtsun, ywbgtshade, ytglbsun, ytglbshade, ytwetb, yQ1, yQ2, &
        yq3, yq4, yq5, yq6, yq7, &
        z0en, z0veg, z0tveg, qdiagtyp, tdiagtyp, udiagtyp, vdiagtyp, &
        qdiagtypv, tdiagtypv, udiagtypv, vdiagtypv, tddiagtyp, tddiagtypv
   character(len=2) :: nm, nagg, nrow
   !--------   FOR CSLM -----------------
   integer :: tke, hdpth, lkiceh, sniceh, expw, dtemp, delu, gred, rhomix, tsed, roficeh, &
              snol, rhosnol, tsnowl, albsnol, wsnowl, tlak, lst, hlaksil, gridarea, &
             ficl, lakd, lfxi, lfxo, lstd, lstf, lakearea, lakefr, riverfr, evlak
   character(len=2) :: nlklv
   !--------   FOR SVS -----------------
   character(len=2) :: ngl, nglp1, nstel, nstpl, iemib, iicel, izp, izvg2
   integer :: acroot, algr, alvl , alvh, avg_gwsol, clayen, co2i1, cvh, cvl, d50, d95, &
        deciduous, draindens, eg, emis, emisgr, emistg, emistgen, emisvh, emisvl, &
        er, etr, evergreen, &
        fbcof, frootd, gamvh, gamvl, grkef, grksat, hfluxsa, hfluxsv, &
        impervu, &
        khc, ksat, ksatc, laictem, laideci, laiva, laivf26, laivh, laivl, &
        lesv, psi, psisat, psngrvl, psnvh, psnvha, &
        rcctem, resagr, resavg, resasa, resasv, resaef, rglvh, rglvl, &
        rnetsa, rnetsv, rsnowsa, &
        rsnowsv, rveg, sanden, skyview, slop, snodpl, snval, &
        snvden,  snvdp, snvma,  snvro, stomrvh, stomrvl, svs_wta, &
        tground, tsa, tsnavg, tsnow, tsnowveg, &
        tsvavg, tvege,  vegh, vegl, vegtrans, vgctem, &
        wsoilm, wfcdp, wfcint, wsnv, &
        z0ha, z0hbg, z0hvg, z0mland, z0mlanden, z0mvg, z0mvh, z0mvhen, z0mvl, &
        conddry, condsld, quartz, rhosoil, soilhcapz, soilcondz, &
        tperm, tpsoil,wunfrz,watpond,maxpond

   !---------------------------------------------------------------------
   
   include "cslm.cdk"
   if (schmlake == 'CSLM') then
      !  nlklv is the maximum number of lake levels
      write(nlklv,'(i2)') nlakmax
   endif

   
   if (schmsol == 'SVS') then
      ! initialize levels for soil texture data
      call init_soil_text_levels()

      ! number of soil/"ground" layers in svs
      Write(ngl,'(i2)') nl_svs
      ! number of soil/"ground" layers PLUS 1
      Write(nglp1,'(i2)') nl_svs+1
      ! number of layer for ENTRY bus SVS clay and sand var.
      Write(nstel,'(i2)') nl_ste

      ! number of layer for PHYSICS bus SVS clay and sand var.
      Write(nstpl,'(i2)') nl_stp
   endif

   write(nagg,'(i2)') nsurf+1

   !# nm is the number of mosaic 'layers'
   nmos = 0
   write(nm,'(i2)') nmos !#TODO: delete

   !# nl is the number of levels in sea ice
   write(nrow,'(i2)') nl

   !# EMIB is only a required input if isba_soil_emiss=='CLIMATO'
   iemib = '0'
   if (isba_soil_emiss=='CLIMATO') iemib = '1'

   !# ICEL is only a required input if icelac=.true.
   iicel = '0'
   if (icelac) iicel = '1'

   !# ZP is only a required input if .not.z0veg_only
   izp = '1'
   if (z0veg_only) izp = '0'

   !# never read zvg2 for SVS, use Z0LD, with key svs_local_z0m instead
   izvg2 = '1'
   if (schmsol == 'SVS') izvg2 = '0'


   !#TODO: check if schmsol conditional
   PHYVAR2D1(lakefr,       'VN=lakefr       ;ON=LAKF;VD=Lake fraction in the grid cell                                       ;VB=p0        ;MIN=0')
   PHYVAR2D1(riverfr,      'VN=riverfr      ;ON=RIVF;VD=River fraction in the grid cell                                      ;VB=p0        ;MIN=0')
   PHYVAR2D1(cgsat,        'VN=cgsat        ;ON=6I  ;VD=thermal coef. at saturation                                          ;VB=p0')
   PHYVAR2D1(dsst,         'VN=dsst         ;ON=DSST;VD=warm layer diurnal SST increment                                     ;VB=p0')
   PHYVAR2D1(dtdiag,       'VN=dtdiag       ;ON=DLIM;VD=DeltaT at screen level of tdiaglim                                   ;VB=p0')
   PHYVAR2D1(glacier,      'VN=glacier      ;ON=2F  ;VD=continental ice fraction                                             ;VB=p1;IN=GA  ;MIN=0')
   PHYVAR2D1(glsea0,       'VN=glsea0       ;ON=GY  ;VD=sea ice fraction (unmodified)                                        ;VB=p1;IN=LG  ;MIN=0')
   PHYVAR2D1(icedp,        'VN=icedp        ;ON=I8  ;VD=sea ice thickness                                                    ;VB=p1        ;MIN=0')
   PHYVAR2D1(iceline,      'VN=iceline      ;ON=ICEL;VD=ice line                                                             ;VB=p'//iicel)
   PHYVAR3D1(qdiagtyp,     "VN=qdiagtyp     ;ON=DQST;VD=screen level specific humidity for each sfc type   ; MIN=0 ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(qdiagtypv,    "VN=qdiagtypv    ;ON=DQSZ;VD=qdiagtyp for z0 vegetation-only ; MIN=0 ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR2D1(skin_depth,   'VN=skin_depth   ;ON=SDEP;VD=sea surface cold skin depth                                          ;VB=p0')
   PHYVAR2D1(skin_inc,     'VN=skin_inc     ;ON=SINC;VD=sea surface cold skin SST increment                                  ;VB=p0')
   PHYVAR3D1(snodp,        'VN=snodp        ;ON=SD  ;VD=snow depth                                     ;VS=A*'//nagg//'    ;VB=p1        ;MIN=0')
   PHYVAR3D1(tddiagtyp,    "VN=tddiagtyp    ;ON=TDST;VD=screen level dew temperature for each sfc type    ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tddiagtypv,   "VN=tddiagtypv   ;ON=TDSZ;VD=tddiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tdiagtyp,     "VN=tdiagtyp     ;ON=TJST;VD=screen level temperature for each sfc type    ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tdiagtypv,    "VN=tdiagtypv    ;ON=TJSZ;VD=tdiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tglacier,     'VN=tglacier     ;ON=I9  ;VD=glaciers temperature                           ;VS=A*2             ;VB=p1')
   PHYVAR3D1(tmice,        'VN=tmice        ;ON=I7  ;VD=sea ice temperature                            ;VS=A*'//nrow//'    ;VB=p1')
   PHYVAR2D1(tnolim,       'VN=tnolim       ;ON=TNOL;VD=screen level temp without max on gradient                            ;VB=p0')
   PHYVAR2D1(twater,       'VN=twater       ;ON=TM  ;VD=sea surface temperature                                              ;VB=p1')
   PHYVAR3D1(udiagtyp,     "VN=udiagtyp     ;ON=UDST;VD=screen level U-wind for each sfc type         ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(udiagtypv,    "VN=udiagtypv    ;ON=UDSZ;VD=udiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(vdiagtyp,     "VN=vdiagtyp     ;ON=VDST;VD=screen level V-wind for each sfc type         ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(vdiagtypv,    "VN=vdiagtypv    ;ON=VDSZ;VD=vdiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(vegf,         'VN=vegf         ;ON=2V  ;VD=vegetation fractions                           ;VS=A*26            ;VB=p1        ;MIN=0; IN=VF')
   PHYVAR3D1(vegfrac,      'VN=vegfrac      ;ON=K1  ;VD=vegetation fraction                            ;VS=A*5             ;VB=p0        ;MIN=0')
   PHYVAR2D1(urban,        'VN=urban        ;ON=URBF;VD=urban fraction                                                     ;VB=p0')
 
   PHYVAR3DC(yradsun,      "VN=yradsun      ;ON=RTSU;VD=MRT in the exposed sunny street (K)                       ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yradshade,    "VN=yradshade    ;ON=RTHD;VD=MRT in the shaded street (K)                              ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yutcisun,     "VN=yutcisun     ;ON=DXSU;VD=UTCI in the exposed sunny street (C)                      ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yutcishade,   "VN=yutcishade   ;ON=DXHD;VD= UTCI in the shaded street (C)                            ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ywbgtsun,     "VN=ywbgtsun     ;ON=GXSU;VD= WBGT in the exposed sunny street (C)                     ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ywbgtshade,   "VN=ywbgtshade   ;ON=GXHD;VD= WBGT in the shaded street (C)                            ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ytglbsun,     "VN=ytglbsun     ;ON=GTSU;VD=TGlobe in the exposed sunny street (K)                    ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ytglbshade,   "VN=ytglbshade   ;ON=GTHD;VD=TGlobe in the shaded street (K)                           ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ytwetb,       "VN=ytwetb       ;ON=WBT ;VD=Wet-Bulb Temp at zt above the ground (K)                  ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ1,          "VN=yQ1          ;ON=QSSU;VD= Contribution of direct solar rad (W/m2)                  ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ2,          "VN=yQ2          ;ON=QSSK;VD= Contribution of sky SW rad (W/m2)                        ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ3,          "VN=yQ3          ;ON=QLSK;VD= Contribution of sky LW rad (W/m2)                        ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ4,          "VN=yQ4          ;ON=QSRD;VD= Contribution of ground SW rad (W/m2)                     ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ5,          "VN=yQ5          ;ON=QLRD;VD= Contribution of ground LW rad (W/m2)                     ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ6,          "VN=yQ6          ;ON=QSWL;VD= Contribution of facet SW rad (W/m2)                      ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ7,          "VN=yQ7          ;ON=QLWL;VD= Contribution of facet LW rad (W/m2)                      ; VS=A*"//nagg//" ; VB=v0", thermal_stress)

   PHYVAR2D1(z0veg,        'VN=z0veg        ;ON=ZVG2;VD=vegetation roughness length                                          ;VB=p'//izvg2//';MIN=0')
   PHYVAR2D1(z0tveg,       'VN=z0tveg       ;ON=ZVGT;VD=thermodynamic vegetation roughness length                            ;VB=p0;MIN=0')
   PHYVAR2D1(z0en,         'VN=z0en         ;ON=2B  ;VD=roughness length (E)                                                 ;VB=e'//izp//'; IN=ZP')

   if (moyhr > 0) &
        PHYVAR2D1(insmavg,   'VN=insmavg      ;ON=IMAV;VD=integrated soil moist avg over last moyhr hrs                      ;VB=p0        ;MIN=0')

   IF_ISBA: if (schmsol == 'ISBA') then
      PHYVAR2D1(acoef,        'VN=acoef        ;ON=1I  ;VD=a coef. in wgeq                                                   ;VB=p0')
      PHYVAR2D1(alveg,        'VN=alveg        ;ON=AX  ;VD=visible canopy albedo                                             ;VB=p0        ;MIN=0')
      PHYVAR2D1(bcoef,        'VN=bcoef        ;ON=1G  ;VD=slope of retention curve                                          ;VB=p0')
      PHYVAR2D1(c1sat,        'VN=c1sat        ;ON=3I  ;VD=c1 coef. at saturation                                            ;VB=p0')
      PHYVAR2D1(c2ref,        'VN=c2ref        ;ON=4I  ;VD=reference value of c2                                             ;VB=p0')
      PHYVAR2D1(c3ref,        'VN=c3ref        ;ON=5I  ;VD=drainage coef. to deeper soil                                     ;VB=p0')
      PHYVAR3D1(clay,         'VN=clay         ;ON=J2  ;VD=percentage of clay in soil                     ;VS=A*3            ;VB=p1        ;MIN=0')
      PHYVAR2D1(cveg,         'VN=cveg         ;ON=CV  ;VD=thermal coefficient for canopy                                    ;VB=p0')
      PHYVAR2D1(drain,        'VN=drain        ;ON=DR  ;VD=water drainage at bottom of soil layer                            ;VB=p0')
      PHYVAR2D1(drainaf,      'VN=drainaf      ;ON=O1  ;VD=accum. of base drainage                                           ;VB=p0')
      PHYVAR2D1(eflux,        'VN=eflux        ;ON=4F  ;VD=specific hum. flux (=-alfaq)                                      ;VB=v0')
      PHYVAR2D1(emsvc,        'VN=emsvc        ;ON=EMIB;VD=broadband climatological ground-veg emissivity isba               ;VB=p'//iemib//';MIN=0')
      PHYVAR2D1(fvapliq,      'VN=fvapliq      ;ON=HFLQ;VD=surf. evaporation (kg/m2 or mm)                                   ;VB=p0')
      PHYVAR2D1(fvapliqaf,    'VN=fvapliqaf    ;ON=AHFL;VD=accum. surf. evaporation (HFLQ) (kg/m2 or mm)                     ;VB=p0')
      PHYVAR2D1(gamveg,       'VN=gamveg       ;ON=GG  ;VD=stomatal resistance parameter                                     ;VB=p0')
      PHYVAR2D1(husurf,       'VN=husurf       ;ON=FH  ;VD=spec. humid. of the surface                                       ;VB=v0        ;MIN=0')
      PHYVAR2D1(hv,           'VN=hv           ;ON=HV  ;VD=relative humidity of veg. canopy                                  ;VB=v0        ;MIN=0')
      PHYVAR2D1(isoil,        'VN=isoil        ;ON=I2  ;VD=soil volumetric ice contents                                      ;VB=p1        ;MIN=0')
      PHYVAR2D1(lai,          'VN=lai          ;ON=J4  ;VD=leaf area index                                                   ;VB=p0        ;MIN=0')
      PHYVAR2D1(leg,          'VN=leg          ;ON=L2  ;VD=latent heat flux over bare grnd                                   ;VB=v0')
      PHYVAR2D1(legaf,        'VN=legaf        ;ON=O5  ;VD=accum. of bare ground LE flux                                     ;VB=p0')
      PHYVAR2D1(ler,          'VN=ler          ;ON=LR  ;VD=latent heat flux from leaves                                      ;VB=v0')
      PHYVAR2D1(leraf,        'VN=leraf        ;ON=O6  ;VD=accum. of direct veg LE flux                                      ;VB=p0')
      PHYVAR2D1(les,          'VN=les          ;ON=LS  ;VD=latent heat flux over snow                                        ;VB=v0')
      PHYVAR2D1(lesaf,        'VN=lesaf        ;ON=O7  ;VD=accum. of sublimation from snow                                   ;VB=p0')
      PHYVAR2D1(letr,         'VN=letr         ;ON=LT  ;VD=latent heat of evapotransp.                                       ;VB=v0')
      PHYVAR2D1(letraf,       'VN=letraf       ;ON=O8  ;VD=accum. of veg. transpiration                                      ;VB=p0')
      PHYVAR2D1(lev,          'VN=lev          ;ON=LV  ;VD=latent heat flux over vegetation                                  ;VB=v0')
      PHYVAR2D1(levaf,        'VN=levaf        ;ON=O9  ;VD=accum. of evaporation from veg.                                   ;VB=p0')
      PHYVAR2D1(melts,        'VN=melts        ;ON=MLTS;VD=accum. snow melting (kg/m2)                                       ;VB=p0')
      PHYVAR2D1(meltsr,       'VN=meltsr       ;ON=MLTR;VD=accum. snow melting due to rain (kg/m2)                           ;VB=p0')
      PHYVAR2D1(overfl,       'VN=overfl       ;ON=RO  ;VD=overland runoff                                                   ;VB=v0')
      PHYVAR2D1(overflaf,     'VN=overflaf     ;ON=N0  ;VD=accum. of surface runoff                                          ;VB=p0')
      PHYVAR2D1(pcoef,        'VN=pcoef        ;ON=7I  ;VD=p coef. in wgeq                                                   ;VB=p0')
      PHYVAR2D1(psn,          'VN=psn          ;ON=5P  ;VD=fraction of the grid covered by snow                              ;VB=v0        ;MIN=0')
      PHYVAR2D1(psng,         'VN=psng         ;ON=3P  ;VD=fraction of bare ground covered by snow                           ;VB=v0        ;MIN=0')
      PHYVAR2D1(psnv,         'VN=psnv         ;ON=4P  ;VD=fraction of vegetation covered by snow                            ;VB=v0        ;MIN=0')
      PHYVAR2D1(resa,         'VN=resa         ;ON=RD  ;VD=aerodynamic resistance                                            ;VB=p0')
      PHYVAR2D1(rgl,          'VN=rgl          ;ON=RG  ;VD=parameter stomatal resistance                                     ;VB=p0')
      PHYVAR2D1(rnet_s,       'VN=rnet_s       ;ON=NR  ;VD=net radiation (soil only)                                         ;VB=v0')
      PHYVAR2D1(rootdp,       'VN=rootdp       ;ON=D2  ;VD=rooting soil depth                                                ;VB=p0')
      PHYVAR2D1(rst,          'VN=rst          ;ON=R1  ;VD=stomatal resistance                                               ;VB=v0')
      PHYVAR3D1(runofftot,    'VN=runofftot    ;ON=TRUN;VD=total surface runoff                           ;VS=A*'//nagg//' ;VB=v0')
      PHYVAR3D1(runofftotaf,  'VN=runofftotaf  ;ON=TRAF;VD=accum. of total surface runoff                 ;VS=A*'//nagg//' ;VB=p0')
      PHYVAR3D1(sand,         'VN=sand         ;ON=J1  ;VD=percentage of sand in soil                     ;VS=A*3          ;VB=p1        ;MIN=0')
      if (snoalb_anl) then
         PHYVAR2D1(snoalen,      'VN=snoalen      ;ON=5H  ;VD=snow albedo (E)                                                ;VB=e1;IN=I6  ;MIN=0')
      else
         PHYVAR2D1(snoagen,      'VN=snoagen      ;ON=3H  ;VD=age of snow (E)                                                ;VB=e1;IN=XA  ;MIN=0')
      endif
      PHYVAR3D1(snoal,        'VN=snoal        ;ON=I6  ;VD=albedo of snow                                 ;VS=A@'//nm//'   ;VB=p0        ;MIN=0')
      PHYVAR3D1(snoden,       'VN=snoden       ;ON=DN  ;VD=snow density in kg/m3                          ;VS=A@'//nm//'   ;VB=p0        ;MIN=0; IN=DN0')
      PHYVAR2D1(snoma,        'VN=snoma        ;ON=I5  ;VD=snow mass                                                         ;VB=p0        ;MIN=0')
      PHYVAR2D1(snoro,        'VN=snoro        ;ON=7S  ;VD=relative snow density                                             ;VB=p1;IN=DN  ;MIN=0')
      PHYVAR2D1(stomr,        'VN=stomr        ;ON=RS  ;VD=minimum stomatal resistance                                       ;VB=p0')
      PHYVAR3D1(tsoil,        'VN=tsoil        ;ON=I0  ;VD=surface and soil temperatures                  ;VS=A*2          ;VB=p1')
      PHYVAR2D1(wfc,          'VN=wfc          ;ON=J5  ;VD=vol. water content at field cap.                                  ;VB=p0        ;MIN=0')
      PHYVAR2D1(wflux,        'VN=wflux        ;ON=M8  ;VD=water flux from surface to atm.                                   ;VB=v0')
      PHYVAR2D1(wfluxaf,      'VN=wfluxaf      ;ON=N7  ;VD=acc. of soil surface upward water flux                            ;VB=p0')
      PHYVAR2D1(wsat,         'VN=wsat         ;ON=J6  ;VD=vol. water content at saturation                                  ;VB=p0        ;MIN=0')
      PHYVAR3D1(wsnow,        'VN=wsnow        ;ON=I4  ;VD=water in the snow pack                         ;VS=A@'//nm//'   ;VB=p1        ;MIN=0')
      PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=I1  ;VD=soil volumetric water contents                 ;VS=A*2          ;VB=p1        ;MIN=0')
      PHYVAR2D1(wveg,         'VN=wveg         ;ON=I3  ;VD=water retained on the vegetation                                  ;VB=p1        ;MIN=0')
      PHYVAR2D1(wwilt,        'VN=wwilt        ;ON=J7  ;VD=vol. water cont. at wilting pt.                                   ;VB=p0        ;MIN=0')
   endif IF_ISBA

   IF_SVS: if (schmsol == 'SVS') then
! check/add min values !!!
      PHYVAR2D1(accevap,      'VN=accevap      ;ON=ACWF;VD=accum. of actual surf. evap. (kg/m2 or mm)                        ;VB=p0')
      PHYVAR3D1(acroot,       'VN=acroot       ;ON=ACRT;VD=active fraction of roots in soil layer         ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(algr,         'VN=algr         ;ON=ALGR;VD=visible albedo for bare ground                                    ;VB=p0        ;MIN=0')
      PHYVAR2D1(alvh,         'VN=alvh         ;ON=ALVH;VD=visible canopy albedo for high vegetation only                    ;VB=p0        ;MIN=0')
      PHYVAR2D1(alvl,         'VN=alvl         ;ON=ALVL;VD=visible canopy albedo for low vegetation only                     ;VB=p0        ;MIN=0')
      PHYVAR2D1(avg_gwsol,    'VN=avg_gwsol    ;ON=AGWS;VD=average soil moisture stress term                                 ;VB=p0        ;MIN=0')
      PHYVAR3D1(bcoef,        'VN=bcoef        ;ON=1G  ;VD=slope of retention curve                       ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(clay,         'VN=clay         ;ON=J2  ;VD=percentage of clay in soil                     ;VS=A*'//nstpl//';VB=p0        ;MIN=0')
      PHYVAR3D1(clayen,       'VN=clayen       ;ON=2H  ;VD=perc. of clay in soil (E)                      ;VS=A*'//nstel//';VB=e1;IN=J2  ;MIN=0')
      PHYVAR3D1(co2i1,        'VN=co2i1        ;ON=CO3 ;VD=CO2 CONCENTRATION   CTEM                       ;VS=A*9          ;VB=p0')
      PHYVAR3D1(conddry,      'VN=conddry      ;ON=CDRY;VD=dry thermal conductivity for soil              ;VS=A*'//ngl//'    ;VB=p0')
      PHYVAR3D1(condsld,      'VN=condsld      ;ON=CSLD;VD=thermal conductivity for soil solids           ;VS=A*'//ngl//'    ;VB=p0')   
      PHYVAR2D1(cveg,         'VN=cveg         ;ON=CV  ;VD=thermal coefficient for canopy                                    ;VB=p0')
      PHYVAR2D1(cvh,          'VN=cvh          ;ON=CVH ;VD=thermal coefficient for canopy of high veg                        ;VB=p0')
      PHYVAR2D1(cvl,          'VN=cvl          ;ON=CVL ;VD=thermal coefficient for canopy of low veg                         ;VB=p0')
      PHYVAR2D1(d50,          'VN=d50          ;ON=d50 ;VD=depth[m] above which 50% of roots are located                     ;VB=p0')
      PHYVAR2D1(d95,          'VN=d95          ;ON=d95 ;VD=depth[m] above which 95% of roots are located                     ;VB=p0')
      PHYVAR2D1(deciduous,    'VN=deciduous    ;ON=DECI;VD=frac. of high veg. that is deciduous.                             ;VB=p0')
      PHYVAR3D1(drain,        'VN=drain        ;ON=DR  ;VD=water drainage in deep soil layers             ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(drainaf,      'VN=drainaf      ;ON=O1  ;VD=accum. of base drainage                                           ;VB=p0')
      PHYVAR2D1(draindens,    'VN=draindens    ;ON=DRND;VD=drainage density (m/m2)                                           ;VB=p1')
      PHYVAR2D1(eflux,        'VN=eflux        ;ON=EFLX;VD=specific hum. flux (=-alfaq)                                      ;VB=v0')
      PHYVAR2D1(eg,           'VN=eg           ;ON=EG  ;VD=evapo. rate over bare grnd(no frac)                               ;VB=v0')
      PHYVAR2D1(emis,         'VN=emis         ;ON=EMI1;VD=emissivity of nat surface                                         ;VB=p0')
      PHYVAR2D1(emisgr,       'VN=emisgr       ;ON=EMGR;VD=emissivity of bare ground                                         ;VB=p0')
      PHYVAR2D1(emistg,       'VN=emistg       ;ON=EMTG;VD=emissivity land surface with no snow (read-in)                    ;VB=p0')
      if (read_emis) &
           PHYVAR2D1(emistgen,     'VN=emistgen     ;ON=ETG1;VD=avg. emissivity land surface with no snow (E)       ;VB=e1;IN=EMIB;MIN=0')
      PHYVAR2D1(emisvh,       'VN=emisvh       ;ON=EMVH;VD=emissivity of high vegetation                                     ;VB=p0')
      PHYVAR2D1(emisvl,       'VN=emisvl       ;ON=EMVL;VD=emissivity of low vegetation                                      ;VB=p0')
      PHYVAR2D1(er,           'VN=er           ;ON=ER  ;VD=evapo rate from leaves(no frac)                                   ;VB=v0')
      PHYVAR2D1(etr,          'VN=etr          ;ON=ETR ;VD=evapotranspiration rate (no frac)                                 ;VB=v0')
      PHYVAR2D1(evergreen,    'VN=evergreen    ;ON=EVER;VD=frac. of high veg. that is evergreen                              ;VB=p0')
      PHYVAR3D1(fbcof,        'VN=fbcof        ;ON=3G  ;VD=parameter derived from bcoef                   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(frootd,       'VN=frootd       ;ON=FRTD;VD=deep soil layer root density                   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(fvapliq,      'VN=fvapliq      ;ON=HFLQ;VD=surf. evaporation (kg/m2 or mm)                                   ;VB=p0')
      PHYVAR2D1(fvapliqaf,    'VN=fvapliqaf    ;ON=AHFL;VD=accum. surf. evaporation (HFLQ) (kg/m2 or mm)                     ;VB=p0')
      PHYVAR2D1(gamvh,        'VN=gamvh        ;ON=GGVH;VD=stomatal resistance parameter for high veg                        ;VB=p0')
      PHYVAR2D1(gamvl,        'VN=gamvl        ;ON=GGVL;VD=stomatal resistance parameter for low veg                         ;VB=p0')
      PHYVAR2D1(grkef,        'VN=grkef        ;ON=GKE; VD=WATDR parameter                                                   ;VB=p0')
      PHYVAR3D1(grksat,       'VN=grksat       ;ON=GKS  ;VD=sat. horiz. soil hydraulic conductivity       ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(hfluxsa,      'VN=hfluxsa      ;ON=HFSA;VD=sensible heat flux (snow only)                                    ;VB=p0')
      PHYVAR2D1(hfluxsv,      'VN=hfluxsv      ;ON=HFSV;VD=sensible heat flux (snow under veg. only)                         ;VB=p0')
      PHYVAR2D1(husurf,       'VN=husurf       ;ON=FH  ;VD=spec. humid. of the surface                                       ;VB=v0        ;MIN=0')
      PHYVAR2D1(hv,           'VN=hv           ;ON=HV  ;VD=relative humidity of veg. canopy                                  ;VB=v0        ;MIN=0')
      PHYVAR2D1(impervu,      'VN=impervu      ;ON=IMPU;VD=frac. of land sfc considered impervious (urban)                   ;VB=p0')
      PHYVAR3D1(isoil,        'VN=isoil        ;ON=ISOL;VD=soil volumetric ice contents per layer         ;VS=A*'//ngl//'  ;VB=p1        ;MIN=0')
      PHYVAR3D1(khc,          'VN=khc          ;ON=KHC ;VD=soil hydraulic conductivity                    ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(ksat,         'VN=ksat         ;ON=KSAT  ;VD=sat. soil hydraulic conductivity             ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(ksatc,        'VN=ksatc        ;ON=KSTC  ;VD=corrected sat. soil hydraulic conductivity   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(laictem,      'VN=laictem      ;ON=LC    ;VD=vegetation LAI for 9 CTEM plant classes      ;VS=A*9          ;VB=p0')
      PHYVAR2D1(laideci,      'VN=laideci      ;ON=LAID;VD=leaf area index for high deciduous veg. only                      ;VB=p0')
      PHYVAR2D1(laiva,        'VN=laiva        ;ON=LAIA;VD=avg. leaf area index seen from atm.                               ;VB=p0')
      PHYVAR3D1(laivf26,      'VN=laivf26      ;ON=LAVF;VD=lai for each vf class times fractipm           ;VS=A*26         ;VB=p0')
      PHYVAR2D1(laivh,        'VN=laivh        ;ON=LAIH;VD=leaf area index for high vegetation only                          ;VB=p0')
      PHYVAR2D1(laivl,        'VN=laivl        ;ON=LAIL;VD=leaf area index for low vegetation only                           ;VB=p0')
      PHYVAR2D1(latflaf,      'VN=latflaf      ;ON=ALAT;VD=Accum. of LATF at all levels (kg/m2 = mm)                         ;VB=p0')
      PHYVAR3D1(latflw,       'VN=latflw       ;ON=LATF;VD=Lateral flow                                   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(leg,          'VN=leg          ;ON=L2  ;VD=latent heat flux over bare grnd                                   ;VB=v0')
      PHYVAR2D1(ler,          'VN=ler          ;ON=LR  ;VD=latent heat flux from leaves                                      ;VB=v0')
      PHYVAR2D1(les,          'VN=les          ;ON=LS  ;VD=latent heat flux over snow                                        ;VB=v0')
      PHYVAR2D1(lesv,         'VN=lesv         ;ON=LSV ;VD=latent heat flux over snow-under-veg                              ;VB=v0')
      PHYVAR2D1(letr,         'VN=letr         ;ON=LT  ;VD=latent heat of evapotransp.                                       ;VB=v0')
      PHYVAR2D1(lev,          'VN=lev          ;ON=LV  ;VD=latent heat flux over vegetation                                  ;VB=v0')
      PHYVAR2D1(maxpond,      'VN=maxpond      ;ON=MAXP;VD=maximum depth[m] of ponded water at surface                       ;VB=p0')
      PHYVAR2D1(melts,        'VN=melts        ;ON=MLTS;VD=accum. snow melting (kg/m2)                                       ;VB=p0')
      PHYVAR2D1(meltsr,       'VN=meltsr       ;ON=MLTR;VD=accum. snow melting due to rain (kg/m2)                           ;VB=p0')
      PHYVAR3D1(psi,          'VN=psi          ;ON=PSI ;VD=soil water suction                             ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(psisat,       'VN=psisat       ;ON=D5  ;VD=sat. soil water suction                        ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(psngrvl,      'VN=psngrvl      ;ON=PSGL;VD=frac. of bare soil &/or low veg. cov. by snow                     ;VB=v0')
      PHYVAR2D1(psnvh,        'VN=psnvh        ;ON=PSVH;VD=fraction of ground -under high vegetation- covered by snow        ;VB=p0')
      PHYVAR2D1(psnvha,       'VN=psnvha       ;ON=PSVA;VD=snow fraction on ground under high veg., as seen through high veg.;VB=p0')
      PHYVAR3D1(quartz,       'VN=quartz       ;ON=QRTZ;VD=quartz contain of the soil layer               ;VS=A*'//ngl//'    ;VB=p0')
      PHYVAR2D1(rcctem,       'VN=rcctem       ;ON=RCC ;VD=stomatal resistance CTEM                                          ;VB=p0')
      PHYVAR2D1(resagr,       'VN=resagr       ;ON=RSGR;VD=aerodynamic resistance over bare ground                           ;VB=p0')
      PHYVAR2D1(resavg,       'VN=resavg       ;ON=RSVG;VD=aerodynamic resistance over veget.                                ;VB=p0')
      PHYVAR2D1(resasa,       'VN=resasa       ;ON=RSSA;VD=aerodynamic resistance for snow on bg/low veg                     ;VB=p0')
      PHYVAR2D1(resasv,       'VN=resasv       ;ON=RSSV;VD=aerodynamic resistance for snow under high veg                    ;VB=p0')
      PHYVAR2D1(resaef,       'VN=resaef       ;ON=RSEF;VD=effective aerodynamic resistance for SVS land tile                ;VB=p0')
      PHYVAR2D1(rglvh,        'VN=rglvh        ;ON=RGVH;VD=parameter stomatal resistance for high veg                        ;VB=p0')
      PHYVAR2D1(rglvl,        'VN=rglvl        ;ON=RGVL;VD=parameter stomatal resistance for low veg                         ;VB=p0')
      PHYVAR3D1(rhosoil,      'VN=rhosoil      ;ON=RHSL;VD=soil dry density                               ;VS=A*'//ngl//'    ;VB=p0')  
      PHYVAR2D1(rnet_s,       'VN=rnet_s       ;ON=NR  ;VD=net radiation (soil only)                                         ;VB=v0')
      PHYVAR2D1(rnetsa,       'VN=rnetsa       ;ON=RNSA;VD=net radiation (snow only)                                         ;VB=p0')
      PHYVAR2D1(rnetsv,       'VN=rnetsv       ;ON=RNSV;VD=net radiation (snow under veg. only)                              ;VB=p0')
      PHYVAR2D1(rootdp,       'VN=rootdp       ;ON=D2  ;VD=rooting soil depth                                                ;VB=p0')
      PHYVAR2D1(rsnowsa,      'VN=rsnowsa      ;ON=RSA ;VD=liquid water out of the snow pack                                 ;VB=p0')
      PHYVAR2D1(rsnowsv,      'VN=rsnowsv      ;ON=RSV ;VD=liquid water out of the snow-under-veg pack                       ;VB=p0')
      PHYVAR2D1(rst,          'VN=rst          ;ON=R1  ;VD=stomatal resistance                                               ;VB=v0')
      PHYVAR3D1(runofftot,    'VN=runofftot    ;ON=TRUN;VD=total surface runoff                           ;VS=A*'//nagg//'   ;VB=v0')
      PHYVAR3D1(runofftotaf,  'VN=runofftotaf  ;ON=TRAF;VD=accum. of total surface runoff                 ;VS=A*'//nagg//'   ;VB=p0')
      PHYVAR2D1(rveg,         'VN=rveg         ;ON=RVG ;VD=runoff from the vegetation (mm/s)                                 ;VB=p0')
      PHYVAR3D1(sand,         'VN=sand         ;ON=J1  ;VD=percentage of sand in soil                     ;VS=A*'//nstpl//'  ;VB=p0')
      PHYVAR3D1(sanden,       'VN=sanden       ;ON=2G  ;VD=perc. of sand in soil (E)                      ;VS=A*'//nstel//'  ;VB=e1;IN=J1  ;')
      PHYVAR2D1(skyview,      'VN=skyview      ;ON=SVF ;VD=sky view factor for tall vegetation                               ;VB=p0')
      PHYVAR2D1(slop,         'VN=slop         ;ON=SLOP;VD=average maximum subgrid-scale topo slope (nil)                    ;VB=p1')
      PHYVAR2D1(snoal,        'VN=snoal        ;ON=SNAL;VD=snow-over-low-veg/bare-ground albedo                              ;VB=p1        ;MIN=0')
      PHYVAR2D1(snoden,       'VN=snoden       ;ON=SNDN;VD=snow-over-low-veg/bare-ground density in kg/m3                    ;VB=p1        ;MIN=0')
      PHYVAR2D1(snodpl,       'VN=snodpl       ;ON=SNDP;VD=snow-over-low-veg/bare-ground depth                               ;VB=p1        ;MIN=0')
      PHYVAR2D1(snoma,        'VN=snoma        ;ON=SNM ;VD=snow-over-low-veg/bare-ground mass                                ;VB=p0        ;MIN=0')
      PHYVAR2D1(snoro,        'VN=snoro        ;ON=SNDR;VD=snow-over-low-veg/bare-ground relative density                    ;VB=p0        ;MIN=0')
      PHYVAR2D1(snval,        'VN=snval        ;ON=SVAL;VD=snow-under-high-veg albedo                                        ;VB=p1')
      PHYVAR2D1(snvden,       'VN=snvden       ;ON=SVDN;VD=snow-under-high-veg density in kg/m3                              ;VB=p1')
      PHYVAR2D1(snvdp,        'VN=snvdp        ;ON=SVDP;VD=snow-under-high-veg depth                                         ;VB=p1')
      PHYVAR2D1(snvma,        'VN=snvma        ;ON=SVM ;VD=snow-under-high-veg mass                                          ;VB=p0')
      PHYVAR2D1(snvro,        'VN=snvro        ;ON=SVDR;VD=snow-under-high-veg relative density                              ;VB=p0')
      PHYVAR3D1(soilhcapz,    'VN=soilhcapz    ;ON=SLCA;VD=soil heat capacity                             ;VS=A*'//ngl//'    ;VB=p0')         
      PHYVAR3D1(soilcondz,    'VN=soilcondz    ;ON=SLCO;VD=soil heat conductivity                         ;VS=A*'//ngl//'    ;VB=p0')  
      PHYVAR2D1(stomrvh,      'VN=stomrvh      ;ON=RSVH;VD=min. stomatal resistance for high vegetation                      ;VB=p0')
      PHYVAR2D1(stomrvl,      'VN=stomrvl      ;ON=RSVL;VD=min. stomatal resistance for low vegetation                       ;VB=p0')
      PHYVAR3D1(svs_wta,      'VN=svs_wta      ;ON=SVSW;VD=weight for svs used in aggregation **FROM SPACE;VS=A*5            ;VB=p0')
      PHYVAR3D1(tground,      'VN=tground      ;ON=TGR ;VD=skin and mean ground temp.                     ;VS=A*2            ;VB=p1')
      PHYVAR2D1(tsa,          'VN=tsa          ;ON=TSA ;VD=skin temp. of land surface as seen from atm                       ;VB=p0')
      PHYVAR2D1(tsnavg,       'VN=tsnavg       ;ON=ATSN;VD=snow-low-veg/bare-grnd avg temp. for melt/freez                   ;VB=p0')
      PHYVAR3D1(tsnow,        'VN=tsnow        ;ON=TSN ;VD=snow-low-veg/bare-grnd skin and mean temp.     ;VS=A*2            ;VB=p1')
      PHYVAR3D1(tsnowveg,     'VN=tsnowveg     ;ON=TSNV;VD=snow-under-high-veg skin and mean temp.        ;VS=A*2            ;VB=p1')
      PHYVAR2D1(tsvavg,       'VN=tsvavg       ;ON=ATSV;VD=snow-under-high-veg avg temp. for melt/freez                      ;VB=p0')
      PHYVAR2D1(tperm,        'VN=tperm        ;ON=TPRD;VD=constant deep soil temperature                                    ;VB=p0')
      PHYVAR3D1(tpsoil,       'VN=tpsoil       ;ON=TGRD;VD=soil temperature profil                        ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(tvege,        'VN=tvege        ;ON=TVG ;VD=skin and mean vegetation temp.                 ;VS=A*2            ;VB=p1')
      PHYVAR3D1(vegdati,      'VN=vegdati      ;ON=SPAR;VD=sparsness of each vegtation class                 ;VS=A*26        ;VB=p0')
      PHYVAR3D1(vegf_evol,    'VN=vegf_evol    ;ON=VFEV;VD=vegetation fraction * vegdat: actual veg. fraction;VS=A*26        ;VB=p0')
      PHYVAR2D1(vegh,         'VN=vegh         ;ON=VEGH;VD=fraction of grid covered by high vegetation                       ;VB=p0')
      PHYVAR2D1(vegl,         'VN=vegl         ;ON=VEGL;VD=fraction of grid covered by low vegetation                        ;VB=p0')
      PHYVAR2D1(vegtrans,     'VN=vegtrans     ;ON=VGTR;VD=transmissivity of tall vegetation                                 ;VB=p0')
      PHYVAR3D1(vgctem,       'VN=vgctem       ;ON=VGCT;VD=CTEM vegetation type fractions                 ;VS=A*9          ;VB=p0')
      PHYVAR3D1(watflow,      'VN=watflow      ;ON=WFL ;VD=waterflow between layers                       ;VS=A*'//nglp1//';VB=p0')
      PHYVAR2D1(watpond,      'VN=watpond      ;ON=WPON;VD=depth of water (m) ponding at the surface                           ;VB=p0')
      PHYVAR3D1(wfc,          'VN=wfc          ;ON=WFC ;VD=vol. water content at field cap.               ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(wfcdp,        'VN=wfcdp        ;ON=WFCD;VD=vol. water content at field cap. at lowst layer                   ;VB=p0')
      PHYVAR3D1(wfcint,       'VN=wfcint       ;ON=WFCI;VD=water content at field capacity along slope    ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(wflux,        'VN=wflux        ;ON=M8  ;VD=water flux from surface to atm.                                   ;VB=v0')
! wfluxaf --- seems to be replaced by accevap in SVS... 
!      PHYVAR2D1(wfluxaf,      'VN=wfluxaf      ;ON=N7  ;VD=acc. of soil surface upward water flux                            ;VB=p0')
      PHYVAR3D1(wsat,         'VN=wsat         ;ON=WSAT;VD=vol. water content at saturation               ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(wsnow,        'VN=wsnow        ;ON=WSN ;VD=water in low-veg/bare-grnd snowpack                               ;VB=p1        ;MIN=0')
      PHYVAR2D1(wsnv,         'VN=wsnv         ;ON=WSV ;VD=water in under-high-veg snowpack                                  ;VB=p1        ;MIN=0')
      PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=WSOL;VD=soil volm water content per layer              ;VS=A*'//ngl//'  ;VB=p1')
      PHYVAR2D1(wsoilm,       'VN=wsoilm       ;ON=WSLM;VD=mean soil volm watr cont for the whole column                     ;VB=p0')
      PHYVAR2D1(wveg,         'VN=wveg         ;ON=WVEG;VD=water retained on the vegetation                                  ;VB=p1        ;MIN=0')
      PHYVAR3D1(wunfrz,       'VN=wunfrz       ;ON=WUFR;VD=unfrozen residual liquid water content           ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(wwilt,        'VN=wwilt        ;ON=WWLT;VD=vol. water cont. at wilting pt.                ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(z0ha,         'VN=z0ha         ;ON=Z0HA;VD=thermal roughness for snowless veg.                               ;VB=p0')
      PHYVAR2D1(z0hbg,        'VN=z0hbg        ;ON=ZTBG;VD=thermal roughness for bare ground in SVS                          ;VB=p0')
      PHYVAR2D1(z0hvg,        'VN=z0hvg        ;ON=ZTVG;VD=thermal roughness for vegetation in SVS                           ;VB=p0')
      PHYVAR2D1(z0mland,      'VN=z0mland      ;ON=Z0LD;VD=local mom roughness length for SVS(snow-free,land-only,no oro.)   ;VB=p0')
      if (svs_dynamic_z0h .or. svs_local_z0m ) &
      PHYVAR2D1(z0mlanden,      'VN=z0mlanden  ;ON=SVS1;VD=local mom roughness length for SVS(snow-free,land-only,no oro.)(E);VB=e1; IN=Z0LD;MIN=0')
      PHYVAR2D1(z0mvg,        'VN=z0mvg        ;ON=Z0MV;VD=momentum roughness length for snow-free vegetation(SVS)           ;VB=p0')
      PHYVAR2D1(z0mvh,        'VN=z0mvh        ;ON=Z0VH;VD=local mom roughness length for high veg.                          ;VB=p0')
      if(read_z0vh) &
      PHYVAR2D1(z0mvhen,      'VN=z0mvhen    ;ON=ZVH1;VD=local mom roughness length for high veg. (E)                        ;VB=e1; IN=Z0VH;MIN=0')
      PHYVAR2D1(z0mvl,        'VN=z0mvl        ;ON=Z0VL;VD=local mom roughness length for low veg.                           ;VB=p0')
   endif IF_SVS

   IF_LAKES: if (schmlake == 'CSLM') then
      PHYVAR3D1(tlak,          'VN=tlak        ;ON=TLAK;   VD=lake temperature profile                 ;VS=A*'//nlklv//'    ;VB=p1')

      

      PHYVAR2D1(tke,          'VN=tke          ;ON=TKEL    ;VD=lake mixed layer tke                                          ;VB=p1')
      PHYVAR2D1(rofinlak,     'VN=rofinlak     ;ON=ROFL    ;VD=runoff input to lake                                          ;VB=p1')
      PHYVAR2D1(rofinlakaf,   'VN=rofinlakaf   ;ON=RLAF    ;VD=accumulated runoff input to lake                              ;VB=p0')
      PHYVAR2D1(lfxi,         'VN=lfxi         ;ON=LFXI    ;VD=accumulated lake water flux in (kg/m2)                        ;VB=p0')
      PHYVAR2D1(lfxo,         'VN=lfxo         ;ON=LFXO    ;VD=accumulated lake water flux out (kg/m2)                       ;VB=p0')
      PHYVAR2D1(evlak,        'VN=evlak        ;ON=EVLA    ;VD=accumulated lake evaporation removed from runoff (kg/m2)      ;VB=p0')
      PHYVAR2D1(lstd,         'VN=lstd         ;ON=LSTD    ;VD=initial lake water storage (kg/m2)                            ;VB=p0')
      PHYVAR2D1(lstf,         'VN=lstf         ;ON=LSTF    ;VD=final lake water storage (kg/m2)                              ;VB=p0')
      PHYVAR2D1(gridarea,     'VN=gridarea     ;ON=AREA    ;VD=surface area (m2) of grid cell                                ;VB=p0')
      PHYVAR2D1(lakearea,     'VN=lakearea     ;ON=LACS    ;VD=surface area (km2) of lake                                    ;VB=p1')
      PHYVAR2D1(lakd,         'VN=lakd         ;ON=LAKD;    VD=Mean lake depth (m)                                           ;VB=p1; IN=DEEP;')
      PHYVAR2D1(lst,          'VN=lst          ;ON=LST     ;VD=lake surface temperature                                      ;VB=p1')
      PHYVAR2D1(hlaksil,      'VN=hlaksil      ;ON=HSIL    ;VD=height of lake water above sill                               ;VB=p1')
      PHYVAR2D1(tsed,         'VN=tsed         ;ON=TSDL    ;VD=lake sediment temperature                                     ;VB=p1')
      PHYVAR2D1(snol,         'VN=snol         ;ON=SNOL    ;VD=snow on lake ice                                              ;VB=p1')
      PHYVAR2D1(ficl,         'VN=ficl         ;ON=FICL    ;VD=fractional lake ice cover                                     ;VB=p0')
      PHYVAR2D1(rhosnol,      'VN=rhosnol      ;ON=RSNL    ;VD=density ofsnow on lake ice                                    ;VB=p1')
      PHYVAR2D1(tsnowl,       'VN=tsnowl       ;ON=TSNL    ;VD=temp of snow on lake ice                                      ;VB=p1')
      PHYVAR2D1(albsnol,      'VN=albsnol      ;ON=ASNL    ;VD=alb. of snow on lake ice                                      ;VB=p1')
      PHYVAR2D1(wsnowl,       'VN=wsnowl       ;ON=WSNL    ;VD=liq.water of snow on lake ice                                 ;VB=p1')
      PHYVAR2D1(delu,         'VN=delu         ;ON=DELU    ;VD=current jump across lake m.l.                                 ;VB=p1')
      PHYVAR2D1(hdpth,        'VN=hdpth        ;ON=HDPL    ;VD=lake mixed layer depth                                        ;VB=p0')
      PHYVAR2D1(lkiceh,       'VN=lkiceh       ;ON=LICE    ;VD=total lake ice thickness                                      ;VB=p1')
      PHYVAR2D1(sniceh,       'VN=sniceh       ;ON=SICE    ;VD=lake snow-ice thickness                                       ;VB=p1')
      PHYVAR2D1(expw,         'VN=expw         ;ON=EXPWL   ;VD=expansivity of lake water                                     ;VB=p0')
      PHYVAR2D1(dtemp,        'VN=dtemp        ;ON=DTEMPL  ;VD=temp jump across lake m.l.                                    ;VB=p0')
      PHYVAR2D1(gred,         'VN=gred         ;ON=GREDL   ;VD=red. grav. across lake m.l.                                   ;VB=p0')
      PHYVAR2D1(rhomix,       'VN=rhomix       ;ON=RHOMIXL ;VD=lake mix. layer density                                       ;VB=p0')
      PHYVAR2D1(roficeh,      'VN=roficeh      ;ON=ROFICEHL ;VD=lake ice from snow melt                                      ;VB=p0')
   endif IF_LAKES

   IF_RIVERS: if (schmriver /= 'NIL') then
      !PHYVAR2D1(riverfr,     'VN=riverfr      ;ON=RIVF;VD=River fraction in the grid cell                    ;VB=p0        ;MIN=0')
   endif IF_RIVERS


   if (schmurb /= 'TEB') return

   PHYVAR2D1(alb_road,     'VN=alb_road     ;ON=ALRD;VD=road albedo                                               ;VB=p0         ;MIN=0')
   PHYVAR2D1(alb_roaden,   'VN=alb_roaden   ;ON=TB9 ;VD=road albedo (E)                                           ;VB=e1; IN=ALRD;MIN=0')
   PHYVAR2D1(alb_roof,     'VN=alb_roof     ;ON=ALRF;VD=roof albedo                                               ;VB=p0         ;MIN=0')
   PHYVAR2D1(alb_roofen,   'VN=alb_roofen   ;ON=TB10;VD=roof albedo (E)                                           ;VB=e1; IN=ALRF;MIN=0')
   PHYVAR2D1(alscatw,      'VN=alscatw      ;ON=ALSC;VD=Town albedo for scattered solar radiation                 ;VB=p0')
   PHYVAR2D1(alb_wall,     'VN=alb_wall     ;ON=ALWL;VD=wall albedo                                               ;VB=p0         ;MIN=0')
   PHYVAR2D1(alb_wallen,   'VN=alb_wallen   ;ON=TB11;VD=wall albedo (E)                                           ;VB=e1; IN=ALWL;MIN=0')
   PHYVAR2D1(azim,         'VN=azim         ;ON=AZIM;VD=solar azimuthal angle                                     ;VB=p0')
   PHYVAR2D1(bld,          'VN=bld          ;ON=BLDF;VD=building fraction                                         ;VB=p0         ;MIN=0')
   PHYVAR2D1(blden,        'VN=blden        ;ON=TB2 ;VD=building fraction (E)                                     ;VB=e1; IN=BLDF;MIN=0')
   PHYVAR2D1(bld_height,   'VN=bld_height   ;ON=BLDH;VD=building height                                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(bld_heighten, 'VN=bld_heighten ;ON=TB3 ;VD=building height (E)                                       ;VB=e1; IN=BLDH;MIN=0')
   PHYVAR2D1(can_hw_ratio, 'VN=can_hw_ratio ;ON=ASPC;VD=aspect ratio of the street                                ;VB=p0         ;MIN=0')
   PHYVAR3D1(d_road,       'VN=d_road       ;ON=DPRD;VD=depth of the road layers                       ;VS=A*3  ;VB=p0')
   PHYVAR3D1(d_roaden,     'VN=d_roaden     ;ON=TB12;VD=depth of the road layers (E)                   ;VS=A*3  ;VB=e1; IN=DPRD;')
   PHYVAR3D1(d_roof,       'VN=d_roof       ;ON=DPRF;VD=depth of the roof layers                       ;VS=A*3  ;VB=p0')
   PHYVAR3D1(d_roofen,     'VN=d_roofen     ;ON=TB13;VD=depth of the roof layers (E)                   ;VS=A*3  ;VB=e1; IN=DPRF;')
   PHYVAR3D1(d_wall,       'VN=d_wall       ;ON=DPWL;VD=depth of the wall layers                       ;VS=A*3  ;VB=p0')
   PHYVAR3D1(d_wallen,     'VN=d_wallen     ;ON=TB14;VD=depth of the wall layers (E)                   ;VS=A*3  ;VB=e1; IN=DPWL;')
   PHYVAR2D1(emis_road,    'VN=emis_road    ;ON=EMRD;VD=road emissivity                                           ;VB=p0')
   PHYVAR2D1(emis_roaden,  'VN=emis_roaden  ;ON=TB15;VD=road emissivity (E)                                       ;VB=e1; IN=EMRD;')
   PHYVAR2D1(emis_roof,    'VN=emis_roof    ;ON=EMRF;VD=roof emissivity                                           ;VB=p0')
   PHYVAR2D1(emis_roofen,  'VN=emis_roofen  ;ON=TB16;VD=roof emissivity (E)                                       ;VB=e1; IN=EMRF;')
   PHYVAR2D1(emis_wall,    'VN=emis_wall    ;ON=EMWL;VD=wall emissivity                                           ;VB=p0')
   PHYVAR2D1(emis_wallen,  'VN=emis_wallen  ;ON=TB17;VD=wall emissivity (E)                                       ;VB=e1; IN=EMWL;')
   PHYVAR2D1(emtw,         'VN=emtw         ;ON=EMTW;VD=Town emissivity                                           ;VB=p0')
   PHYVAR2D1(g_road,       'VN=g_road       ;ON=QGRD;VD=storage heat flux for road                                ;VB=p0')
   PHYVAR2D1(g_roof,       'VN=g_roof       ;ON=QGRF;VD=storage heat flux for roof                                ;VB=p0')
   PHYVAR2D1(g_town,       'VN=g_town       ;ON=QGTW;VD=storage heat flux for town                                ;VB=p0')
   PHYVAR2D1(g_wall,       'VN=g_wall       ;ON=QGWL;VD=storage heat flux for wall                                ;VB=p0')
   PHYVAR2D1(h_industry,   'VN=h_industry   ;ON=QHIN;VD=sensible heat flux from industry                          ;VB=p0')
   PHYVAR2D1(h_industryen, 'VN=h_industryen ;ON=TB18;VD=sensible heat flux from industry (E)                      ;VB=e1; IN=QHIN;')
   PHYVAR2D1(h_road,       'VN=h_road       ;ON=QHRD;VD=sensible heat flux over road                              ;VB=p0')
   PHYVAR2D1(h_roof,       'VN=h_roof       ;ON=QHRF;VD=sensible heat flux over roof                              ;VB=p0')
   PHYVAR2D1(h_town,       'VN=h_town       ;ON=QHTW;VD=sensible heat flux over town                              ;VB=p0')
   PHYVAR2D1(h_traffic,    'VN=h_traffic    ;ON=QHTR;VD=sensible heat flux from traffic                           ;VB=p0')
   PHYVAR2D1(h_trafficen,  'VN=h_trafficen  ;ON=TB28;VD=sensible heat flux from traffic (E)                       ;VB=e1; IN=QHTR;')
   PHYVAR2D1(h_wall,       'VN=h_wall       ;ON=QHWL;VD=sensible heat flux over wall                              ;VB=p0')
   PHYVAR3D1(hc_road,      'VN=hc_road      ;ON=HCRD;VD=road heat capacities                           ;VS=A*3  ;VB=p0')
   PHYVAR3D1(hc_roaden,    'VN=hc_roaden    ;ON=TB19;VD=road heat capacities (E)                       ;VS=A*3  ;VB=e1; IN=HCRD;')
   PHYVAR3D1(hc_roof,      'VN=hc_roof      ;ON=HCRF;VD=roof heat capacities                           ;VS=A*3  ;VB=p0')
   PHYVAR3D1(hc_roofen,    'VN=hc_roofen    ;ON=TB20;VD=roof heat capacities (E)                       ;VS=A*3  ;VB=e1; IN=HCRF;')
   PHYVAR3D1(hc_wall,      'VN=hc_wall      ;ON=HCWL;VD=wall heat capacities                           ;VS=A*3  ;VB=p0')
   PHYVAR3D1(hc_wallen,    'VN=hc_wallen    ;ON=TB21;VD=wall heat capacities (E)                       ;VS=A*3  ;VB=e1; IN=HCWL;')
   PHYVAR2D1(le_industry,  'VN=le_industry  ;ON=QEIN;VD=latent heat flux from industry                            ;VB=p0')
   PHYVAR2D1(le_industryen,'VN=le_industryen;ON=TB27;VD=latent heat flux from industry (E)                        ;VB=e1; IN=QEIN;')
   PHYVAR2D1(le_road,      'VN=le_road      ;ON=QERD;VD=latent heat flux over road                                ;VB=p0')
   PHYVAR2D1(le_roof,      'VN=le_roof      ;ON=QERF;VD=latent heat flux over roof                                ;VB=p0')
   PHYVAR2D1(le_town,      'VN=le_town      ;ON=QETW;VD=latent heat flux over town                                ;VB=p0')
   PHYVAR2D1(le_traffic,   'VN=le_traffic   ;ON=QETR;VD=latent heat flux from traffic                             ;VB=p0')
   PHYVAR2D1(le_trafficen, 'VN=le_trafficen ;ON=TB26;VD=latent heat flux from traffic (E)                         ;VB=e1; IN=QETR;')
   PHYVAR2D1(le_wall,      'VN=le_wall      ;ON=QEWL;VD=latent heat flux over wall                                ;VB=p0')
   PHYVAR2D1(nat,          'VN=nat          ;ON=NATF;VD=natural surface fraction in urban area                    ;VB=p0         ;MIN=0')
   PHYVAR2D1(naten,        'VN=naten        ;ON=TB1 ;VD=natural surface fraction in urban area (E)                ;VB=e1; IN=NATF;MIN=0')
   PHYVAR2D1(pav,          'VN=pav          ;ON=PAVF;VD=impervious fraction (road)                                ;VB=p0         ;MIN=0')
   PHYVAR2D1(paven,        'VN=paven        ;ON=TB4 ;VD=impervious fraction (road) (E)                            ;VB=e1; IN=PAVF;MIN=0')
   PHYVAR2D1(q_canyon,     'VN=q_canyon     ;ON=QCAN;VD=specific humi inside the canyon                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(q_canyonen,   'VN=q_canyonen   ;ON=2XEN; VD=canyon air humi (E)                                      ;VB=e1 ;IN=QCAN;MIN=0')
   PHYVAR2D1(rn_road,      'VN=rn_road      ;ON=RNRD;VD=Net radiation over road                                   ;VB=p0')
   PHYVAR2D1(rn_roof,      'VN=rn_roof      ;ON=RNRF;VD=Net radiation over roof                                   ;VB=p0')
   PHYVAR2D1(rn_town,      'VN=rn_town      ;ON=RNTW;VD=Net radiation over town                                   ;VB=p0')
   PHYVAR2D1(rn_wall,      'VN=rn_wall      ;ON=RNWL;VD=Net radiation over wall                                   ;VB=p0')
   PHYVAR2D1(sroad_alb,    'VN=sroad_alb    ;ON=SARD;VD=snow albedo for roads                                     ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroad_alben,  'VN=sroad_alben  ;ON=9ZEN; VD=snow albedo for roads (E)                                ;VB=e1 ;IN=SARD;MIN=0')
   PHYVAR2D1(sroad_emis,   'VN=sroad_emis   ;ON=SERD;VD=snow emissivity for roads                                 ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroad_emisen, 'VN=sroad_emisen ;ON=1MEN; VD=snow emmissivity for roads (E)                           ;VB=e1 ;IN=SERD;MIN=0')
   PHYVAR2D1(sroad_rho,    'VN=sroad_rho    ;ON=SDRD;VD=snow density for roads                                    ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroad_rhoen,  'VN=sroad_rhoen  ;ON=8ZEN; VD=snow density for roads (E)                               ;VB=e1 ;IN=SDRD;MIN=0')
   PHYVAR2D1(sroad_scheme, 'VN=sroad_scheme ;ON=SCRD;VD=snow scheme for roads                                     ;VB=p0')
   PHYVAR2D1(sroad_t,      'VN=sroad_t      ;ON=STRD;VD=snow temperature for roads                                ;VB=p0')
   PHYVAR2D1(sroad_ten,    'VN=sroad_ten    ;ON=7ZEN; VD=snow temp for roads (E)                                  ;VB=e1 ;IN=STRD;')
   PHYVAR2D1(sroad_ts,     'VN=sroad_ts     ;ON=SSRD;VD=snow surf temperature for roads                           ;VB=p0')
   PHYVAR2D1(sroad_tsen,   'VN=sroad_tsen   ;ON=2MEN; VD=Snow surface temp for roads (E)                          ;VB=e1 ;IN=SSRD;')
   PHYVAR2D1(sroad_wsnow,  'VN=sroad_wsnow  ;ON=SWRD;VD=water and snow content for roads                          ;VB=p0')
   PHYVAR2D1(sroad_wsnowen,'VN=sroad_wsnowen;ON=6ZEN; VD=water/snow content for roads (E)                         ;VB=e1 ;IN=SWRD;')
   PHYVAR2D1(sroof_alb,    'VN=sroof_alb    ;ON=SARF;VD=snow albedo for roofs                                     ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroof_alben,  'VN=sroof_alben  ;ON=1ZEN; VD=snow albedo for roofs (E)                                ;VB=e1 ;IN=SARF;MIN=0')
   PHYVAR2D1(sroof_emis,   'VN=sroof_emis   ;ON=SERF;VD=snow emissivity for roofs                                 ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroof_emisen, 'VN=sroof_emisen ;ON=2ZEN; VD=snow emmissivity for roofs (E)                           ;VB=e1 ;IN=SERF;MIN=0')
   PHYVAR2D1(sroof_rho,    'VN=sroof_rho    ;ON=SDRF;VD=snow density for roofs                                    ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroof_rhoen,  'VN=sroof_rhoen  ;ON=9YEN; VD=snow density for roofs (E)                               ;VB=e1 ;IN=SDRF;MIN=0')
   PHYVAR2D1(sroof_scheme, 'VN=sroof_scheme ;ON=SCRF;VD=snow scheme for roofs                                     ;VB=p0')
   PHYVAR2D1(sroof_t,      'VN=sroof_t      ;ON=STRF;VD=snow temperature for roofs                                ;VB=p0')
   PHYVAR2D1(sroof_ten,    'VN=sroof_ten    ;ON=8YEN; VD=snow temp for roofs (E)                                  ;VB=e1 ;IN=STRF;')
   PHYVAR2D1(sroof_ts,     'VN=sroof_ts     ;ON=SSRF;VD=snow surf temperature for roofs                           ;VB=p0')
   PHYVAR2D1(sroof_tsen,   'VN=sroof_tsen   ;ON=3ZEN; VD=Snow surface temp for roofs (E)                          ;VB=e1 ;IN=SSRF;')
   PHYVAR2D1(sroof_wsnow,  'VN=sroof_wsnow  ;ON=SWRF;VD=water and snow content for roofs                          ;VB=p0')
   PHYVAR2D1(sroof_wsnowen,'VN=sroof_wsnowen;ON=7YEN; VD=water/snow content for roofs (E)                         ;VB=e1 ;IN=SWRF;')
   PHYVAR2D1(svf_road,     'VN=svf_road     ;ON=SVRD;VD=road sky-view factor                                      ;VB=p0')
   PHYVAR2D1(svf_wall,     'VN=svf_wall     ;ON=SVWL;VD=wall sky-view factor                                      ;VB=p0')
   PHYVAR2D1(t_canyon,     'VN=t_canyon     ;ON=TCAN;VD=air temperature inside the canyon                         ;VB=p0')
   PHYVAR2D1(t_canyonen,   'VN=t_canyonen   ;ON=1XEN; VD=canyon air temp (E)                                      ;VB=e1 ;IN=TCAN;')
   PHYVAR3D1(t_road,       'VN=t_road       ;ON=TLRD;VD=temperatures of road layers                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(t_roaden,     'VN=t_roaden     ;ON=4XEN; VD=road temperatures (E)                         ;VS=A*3  ;VB=e1 ;IN=TLRD;')
   PHYVAR3D1(t_roof,       'VN=t_roof       ;ON=TLRF;VD=temperatures of roof layers                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(t_roofen,     'VN=t_roofen     ;ON=3XEN; VD=roof temperatures (E)                         ;VS=A*3  ;VB=e1 ;IN=TLRF;')
   PHYVAR3D1(t_wall,       'VN=t_wall       ;ON=TLWL;VD=temperatures of wall layers                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(t_wallen,     'VN=t_wallen     ;ON=5XEN; VD=wall temperatures (E)                         ;VS=A*3  ;VB=e1 ;IN=TLWL;')
   PHYVAR3D1(tc_road,      'VN=tc_road      ;ON=TCRD;VD=road thermal condcutivities                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(tc_roaden,    'VN=tc_roaden    ;ON=TB22;VD=road thermal condcutivities (E)                ;VS=A*3  ;VB=e1; IN=TCRD;')
   PHYVAR3D1(tc_roof,      'VN=tc_roof      ;ON=TCRF;VD=roof thermal condcutivities                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(tc_roofen,    'VN=tc_roofen    ;ON=TB23;VD=roof thermal condcutivities (E)                ;VS=A*3  ;VB=e1; IN=TCRF;')
   PHYVAR3D1(tc_wall,      'VN=tc_wall      ;ON=TCWL;VD=wall thermal condcutivities                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(tc_wallen,    'VN=tc_wallen    ;ON=TB24;VD=wall thermal condcutivities (E)                ;VS=A*3  ;VB=e1; IN=TCWL;')
   PHYVAR2D1(ti_bld,       'VN=ti_bld       ;ON=TBLD;VD=internal building temperature                             ;VB=p0')
   PHYVAR2D1(ti_blden,     'VN=ti_blden     ;ON=6QEN; VD=bld internal temp (E)                                    ;VB=e1 ;IN=TBLD;')
   PHYVAR2D1(ti_road,      'VN=ti_road      ;ON=TIRD;VD=internal road temperature                                 ;VB=p0')
   PHYVAR2D1(ti_roaden,    'VN=ti_roaden    ;ON=5QEN; VD=road internal temp (E)                                   ;VB=e1 ;IN=TIRD;')
   PHYVAR2D1(tsradtw,      'VN=tsradtw      ;ON=TSTW;VD=Town radiative surface temperature                        ;VB=p0')
   PHYVAR2D1(tsun,         'VN=tsun         ;ON=TSUN;VD=solar time (s)                                            ;VB=p0')
   PHYVAR2D1(u_canyon,     'VN=u_canyon     ;ON=UCAN;VD=wind in canyon                                            ;VB=p0')
   PHYVAR2D1(wall_o_hor,   'VN=wall_o_hor   ;ON=WHOR;VD=ratio vertical per horizontal surf                        ;VB=p0         ;MIN=0')
   PHYVAR2D1(wall_o_horen, 'VN=wall_o_horen ;ON=TB5 ;VD=ratio vertical per horizontal surf (E)                    ;VB=e1; IN=WHOR;MIN=0')
   PHYVAR2D1(ws_road,      'VN=ws_road      ;ON=WSRD;VD=water content of road reservoir                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(ws_roaden,    'VN=ws_roaden    ;ON=4QEN; VD=road water reservoir (E)                                 ;VB=e1 ;IN=WSRD;MIN=0')
   PHYVAR2D1(ws_roof,      'VN=ws_roof      ;ON=WSRF;VD=water content of roof reservoir                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(ws_roofen,    'VN=ws_roofen    ;ON=3QEN; VD=roof water reservoir (E)                                 ;VB=e1 ;IN=WSRF;MIN=0')

   PHYVAR2DC(yradin,       'VN=yradin       ;ON=RTIN;VD=MRT inside building  (K)                                  ;VB=v0', thermal_stress)
   PHYVAR2DC(yradrfsun,    'VN=yradrfsun    ;ON=RTFS;VD=MRT on the exposed sunny roof (K)                         ;VB=v0', thermal_stress)
   PHYVAR2DC(yradrfshade,  'VN=yradrfshade  ;ON=RTFD;VD=MRT on the shaded roof (K)                                ;VB=v0', thermal_stress)
   PHYVAR2DC(yutciin,      'VN=yutciin      ;ON=DXIN;VD=UTCI inside building  (C)                                 ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcirfsun,   'VN=yutcirfsun   ;ON=DXFS;VD=UTCI on the exposed sunny roof (C)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcirfshade, 'VN=yutcirfshade ;ON=DXFD;VD=UTCI on the shaded roof (C)                               ;VB=v0', thermal_stress)
   PHYVAR2DC(ytrfzt,       'VN=ytrfzt       ;ON=T2RF;VD=Temperature at zt above the roof                          ;VB=v0', thermal_stress)
   PHYVAR2DC(ytrdzt,       'VN=ytrdzt       ;ON=T2RD;VD=Temperature at zt above the ground                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yurdzu,       'VN=yurdzu       ;ON=UVRD;VD=wind speed at zu above the ground (m/s)                   ;VB=v0', thermal_stress)
   PHYVAR2DC(ywbgtrfsun,   'VN=ywbgtrfsun   ;ON=GXFS;VD=WBGT on the exposed sunny roof (C)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(ywbgtrfshade, 'VN=ywbgtrfshade ;ON=GXFD;VD=WBGT on the shaded roof (C)                               ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicin,     'VN=yutcicin     ;ON=DCIN;VD=cumulative UTCI inside building  (C)                      ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicsun,    'VN=yutcicsun    ;ON=DCSU;VD=cumulative UTCI in the exposed sunny street               ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicshade,  'VN=yutcicshade  ;ON=DCHD;VD=cumulative UTCI in the shaded street (C)                  ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicrfsun,  'VN=yutcicrfsun  ;ON=DCFS;VD=cumulative UTCI on the exposed sunny roof (C)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicrfshade,'VN=yutcicrfshade;ON=DCFD;VD=cumulative UTCI on the shaded roof (C)                    ;VB=v0', thermal_stress)
   PHYVAR2DC(ytglbrfsun,   'VN=ytglbrfsun   ;ON=GTFS;VD=TGlobe on the exposed sunny roof (K)                      ;VB=v0', thermal_stress)
   PHYVAR2DC(ytglbrfshade, 'VN=ytglbrfshade ;ON=GTFD;VD=TGlobe on the shaded roof (K)                             ;VB=v0', thermal_stress)
   PHYVAR2DC(ytwetbrf,     'VN=ytwetbrf     ;ON=WBRF;VD=Wet-Bulb Temperature at zt above the roof (K)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ8,          'VN=yQ8          ;ON=QSRF;VD=Contribution of roof SW rad (W/m2)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ9,          'VN=yQ9          ;ON=QLRF;VD=Contribution of roof LW rad (W/m2)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ10,         'VN=yQ10         ;ON=QSFK;VD=Contribution of sky on the roof SW rad (W/m2)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ11,         'VN=yQ11         ;ON=QLFK;VD=Contribution of sky on the roof LW rad (W/m2)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ12,         'VN=yQ12         ;ON=QSFW;VD=Contribution of wall on the roof SW rad (W/m2)            ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ13,         'VN=yQ13         ;ON=QLFW;VD=Contribution of wall on the roof LW rad (W/m2)            ;VB=v0', thermal_stress)

   PHYVAR2D1(z0_road,      'VN=z0_road      ;ON=Z0RD;VD=aerodyn roughness length for road                         ;VB=p0')
   PHYVAR2D1(z0_roaden,    'VN=z0_roaden    ;ON=TB8 ;VD=aerodyn roughness length for road (E)                     ;VB=e1; IN=Z0RD;')
   PHYVAR2D1(z0_roof,      'VN=z0_roof      ;ON=Z0RF;VD=aerodyn roughness length for roof                         ;VB=p0')
   PHYVAR2D1(z0_roofen,    'VN=z0_roofen    ;ON=TB7 ;VD=aerodyn roughness length for roof (E)                     ;VB=e1; IN=Z0RF;')
   PHYVAR2D1(z0_town,      'VN=z0_town      ;ON=Z0TW;VD=aerodyn roughness length for town                         ;VB=p0')
   PHYVAR2D1(z0_townen,    'VN=z0_townen    ;ON=TB6 ;VD=aerodyn roughness length for town (E)                     ;VB=e1; IN=Z0TW;')
   PHYVAR2D1(zenith,       'VN=zenith       ;ON=ZENI;VD=solar zenith angle                                        ;VB=p0')

   !---------------------------------------------------------------------
   return
end subroutine sfc_businit
