!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module sfcbus_mod
   implicit none
   public
   save

#define SFCVAR(MYVAR,MYNAME) type(SFCVAR_T) :: MYVAR = SFCVAR_T(-1,0,0,0,0,.false.,MYNAME)

   type :: sfcptr
      sequence
      real, pointer :: ptr(:,:)
   end type sfcptr

   type :: SFCVAR_T
      sequence
      integer :: i, agg, mul, niveaux, mosaik
      logical :: doagg_L
      character(len=32) :: n
   end type SFCVAR_T

   type :: SFCVARLIST_T
      sequence

      SFCVAR(umoins , 'PW_UU:M')
      SFCVAR(uplus  , 'PW_UU:P')
      SFCVAR(vmoins , 'PW_VV:M')
      SFCVAR(vplus  , 'PW_VV:P')
      SFCVAR(tmoins , 'PW_TT:M')
      SFCVAR(tplus  , 'PW_TT:P')
      SFCVAR(sigm   , 'PW_PM:P')
      SFCVAR(humoins, 'TR/HU:M')
      SFCVAR(huplus , 'TR/HU:P')

      SFCVAR(acoef, 'acoef')
      SFCVAR(acroot, 'acroot')
      SFCVAR(accevap, 'accevap')
      SFCVAR(alb_road, 'alb_road')
      SFCVAR(alb_roaden, 'alb_roaden')
      SFCVAR(alb_roof, 'alb_roof')
      SFCVAR(alb_roofen, 'alb_roofen')
      SFCVAR(alb_wall, 'alb_wall')
      SFCVAR(alb_wallen, 'alb_wallen')
      SFCVAR(alfaq, 'alfaq')
      SFCVAR(alfat, 'alfat')
      SFCVAR(algr, 'algr')
      SFCVAR(alscatw, 'alscatw')
      SFCVAR(alveg, 'alveg')
      SFCVAR(alvh, 'alvh')
      SFCVAR(alvis, 'alvis')
      SFCVAR(alvl, 'alvl')
      SFCVAR(alwater, 'alwater')
      SFCVAR(avg_gwsol, 'avg_gwsol')
      SFCVAR(azim, 'azim')
      SFCVAR(bcoef, 'bcoef')
      SFCVAR(bld, 'bld')
      SFCVAR(blden, 'blden')
      SFCVAR(bld_height, 'bld_height')
      SFCVAR(bld_heighten, 'bld_heighten')
      SFCVAR(bm, 'bm')
      SFCVAR(bt, 'bt')
      SFCVAR(c1sat, 'c1sat')
      SFCVAR(c2ref, 'c2ref')
      SFCVAR(c3ref, 'c3ref')
      SFCVAR(cang, 'cang')
      SFCVAR(can_hw_ratio, 'can_hw_ratio')
      SFCVAR(cgsat, 'cgsat')
      SFCVAR(clay, 'clay')
      SFCVAR(clayen, 'clayen')
      SFCVAR(co2i1, 'co2i1')
      SFCVAR(cosz, 'cosz')
      SFCVAR(cveg, 'cveg')
      SFCVAR(cvh, 'cvh')
      SFCVAR(cvl, 'cvl')
      SFCVAR(d50, 'd50')
      SFCVAR(d95, 'd95')
      SFCVAR(d_road, 'd_road')
      SFCVAR(d_roaden, 'd_roaden')
      SFCVAR(d_roof, 'd_roof')
      SFCVAR(d_roofen, 'd_roofen')
      SFCVAR(d_wall, 'd_wall')
      SFCVAR(d_wallen, 'd_wallen')
      SFCVAR(deciduous, 'deciduous')
      SFCVAR(dhdx, 'dhdx')
      SFCVAR(dhdxdy, 'dhdxdy')
      SFCVAR(dhdy, 'dhdy')
      SFCVAR(dlat, 'dlat')
      SFCVAR(dlon, 'dlon')
      SFCVAR(drain, 'drain')
      SFCVAR(draindens, 'draindens')
      SFCVAR(drainaf, 'drainaf')
      SFCVAR(dsst, 'dsst')
      SFCVAR(dtdiag, 'dtdiag')
      SFCVAR(eflux, 'eflux')
      SFCVAR(eg, 'eg')
      SFCVAR(emis, 'emis')
      SFCVAR(emisr, 'emisr')
      SFCVAR(emis_road, 'emis_road')
      SFCVAR(emis_roaden, 'emis_roaden')
      SFCVAR(emis_roof, 'emis_roof')
      SFCVAR(emis_roofen, 'emis_roofen')
      SFCVAR(emis_wall, 'emis_wall')
      SFCVAR(emis_wallen, 'emis_wallen')
      SFCVAR(emisgr, 'emisgr')
      SFCVAR(emistg, 'emistg')
      SFCVAR(emistgen, 'emistgen')
      SFCVAR(emsvc, 'emsvc')
      SFCVAR(emisvh, 'emisvh')
      SFCVAR(emisvl, 'emisvl')
      SFCVAR(emtw, 'emtw')
      SFCVAR(en, 'en')
      SFCVAR(er, 'er')
      SFCVAR(etr, 'etr')
      SFCVAR(evergreen, 'evergreen')
      SFCVAR(fbcof, 'fbcof')
      SFCVAR(fc, 'fc')
      SFCVAR(fcor, 'fcor')
      SFCVAR(fdsi, 'fdsi')
      SFCVAR(fdss, 'fdss')
      SFCVAR(fl, 'fl')
      SFCVAR(flusolis, 'flusolis')
      SFCVAR(fluslop, 'fluslop')
      SFCVAR(fq, 'fq')
      SFCVAR(frootd, 'frootd')
      SFCVAR(frv, 'frv')
      SFCVAR(fsd, 'fsd')
      SFCVAR(fsf, 'fsf')
      SFCVAR(ftemp, 'ftemp')
      SFCVAR(fv, 'fv')
      SFCVAR(fvap, 'fvap')
      SFCVAR(fvapliq, 'fvapliq')
      SFCVAR(fvapliqaf, 'fvapliqaf')
      SFCVAR(g_road, 'g_road')
      SFCVAR(g_roof, 'g_roof')
      SFCVAR(g_town, 'g_town')
      SFCVAR(g_wall, 'g_wall')
      SFCVAR(gamveg, 'gamveg')
      SFCVAR(gamvh, 'gamvh')
      SFCVAR(gamvl, 'gamvl')      
      SFCVAR(gc, 'gc')
      SFCVAR(glacier, 'glacier')
      SFCVAR(glsea, 'glsea')
      SFCVAR(glsea0, 'glsea0')
      SFCVAR(grkef, 'grkef')
      SFCVAR(grksat, 'grksat')
      SFCVAR(gztherm, 'gztherm')
      SFCVAR(h, 'h')
      SFCVAR(h_industry, 'h_industry')
      SFCVAR(h_industryen, 'h_industryen')
      SFCVAR(h_road, 'h_road')
      SFCVAR(h_roof, 'h_roof')
      SFCVAR(h_town, 'h_town')
      SFCVAR(h_traffic, 'h_traffic')
      SFCVAR(h_trafficen, 'h_trafficen')
      SFCVAR(h_wall, 'h_wall')
      SFCVAR(hc_road, 'hc_road')
      SFCVAR(hc_roaden, 'hc_roaden')
      SFCVAR(hc_roof, 'hc_roof')
      SFCVAR(hc_roofen, 'hc_roofen')
      SFCVAR(hc_wall, 'hc_wall')
      SFCVAR(hc_wallen, 'hc_wallen')
      SFCVAR(hfluxsa, 'hfluxsa')
      SFCVAR(hfluxsv, 'hfluxsv')
      SFCVAR(hst, 'hst')
      SFCVAR(husurf, 'husurf')
      SFCVAR(hv, 'hv')
      SFCVAR(icedp, 'icedp')
      SFCVAR(iceline, 'iceline')
      SFCVAR(ilmo, 'ilmo')
      SFCVAR(impervu, 'impervu')
      SFCVAR(isoil, 'isoil')
      SFCVAR(kcl, 'kcl')
      SFCVAR(khc, 'khc')
      SFCVAR(km, 'km')
      SFCVAR(ksat, 'ksat')
      SFCVAR(ksatc, 'ksatc')
      SFCVAR(kt, 'kt')
      SFCVAR(lai, 'lai')
      SFCVAR(laictem, 'laictem')
      SFCVAR(laideci, 'laideci')
      SFCVAR(laiva, 'laiva')
      SFCVAR(laivf26, 'laivf26')
      SFCVAR(laivh, 'laivh')
      SFCVAR(laivl, 'laivl')
      SFCVAR(latflaf, 'latflaf')
      SFCVAR(latflw, 'latflw')
      SFCVAR(le_industry, 'le_industry')
      SFCVAR(le_industryen, 'le_industryen')
      SFCVAR(le_road, 'le_road')
      SFCVAR(le_roof, 'le_roof')
      SFCVAR(le_town, 'le_town')
      SFCVAR(le_traffic, 'le_traffic')
      SFCVAR(le_trafficen, 'le_trafficen')
      SFCVAR(le_wall, 'le_wall')
      SFCVAR(leg, 'leg')
      SFCVAR(ler, 'ler')
      SFCVAR(les, 'les')
      SFCVAR(lesv, 'lesv')
      SFCVAR(letr, 'letr')
      SFCVAR(lev, 'lev')
      SFCVAR(lhtg, 'lhtg')
      SFCVAR(melts, 'melts')
      SFCVAR(meltsr, 'meltsr')
      SFCVAR(mf, 'mf')
      SFCVAR(mg, 'mg')
      SFCVAR(ml, 'ml')
      SFCVAR(mt, 'mt')
      SFCVAR(mtdir, 'mtdir')
      SFCVAR(nat, 'nat')
      SFCVAR(overfl, 'overfl')
      SFCVAR(pav, 'pav')
      SFCVAR(paven, 'paven')
      SFCVAR(pcoef, 'pcoef')
      SFCVAR(pmoins, 'pmoins')
      SFCVAR(pplus, 'pplus')
      SFCVAR(psi, 'psi')
      SFCVAR(psisat, 'psisat')
      SFCVAR(psn, 'psn')
      SFCVAR(psng, 'psng')
      SFCVAR(psngrvl, 'psngrvl')
      SFCVAR(psnv, 'psnv')
      SFCVAR(psnvh, 'psnvh')
      SFCVAR(psnvha, 'psnvha')
      SFCVAR(q_canyon, 'q_canyon')
      SFCVAR(q_canyonen, 'q_canyonen')
      SFCVAR(qdiag, 'qdiag')
      SFCVAR(qdiagstn, 'qdiagstn')
      SFCVAR(qdiagstnv, 'qdiagstnv')
      SFCVAR(qdiagtyp, 'qdiagtyp')
      SFCVAR(qdiagtypv, 'qdiagtypv')
      SFCVAR(qsurf, 'qsurf')
      SFCVAR(rainrate, 'rainrate')
      SFCVAR(rcctem, 'rcctem')
      SFCVAR(resa, 'resa')
      SFCVAR(resagr, 'resagr')
      SFCVAR(resavg, 'resavg')
      SFCVAR(resasa, 'resasa')
      SFCVAR(resasv, 'resasv')
      SFCVAR(resaef, 'resaef')
      SFCVAR(rgl, 'rgl')
      SFCVAR(rglvh, 'rglvh')
      SFCVAR(rglvl, 'rglvl')
      SFCVAR(rib, 'rib')
      SFCVAR(rn_road, 'rn_road')
      SFCVAR(rn_roof, 'rn_roof')
      SFCVAR(rn_town, 'rn_town')
      SFCVAR(rn_wall, 'rn_wall')
      SFCVAR(rnet_s, 'rnet_s')
      SFCVAR(rnetsa, 'rnetsa')
      SFCVAR(rnetsv, 'rnetsv')
      SFCVAR(rootdp, 'rootdp')
      SFCVAR(rsnowsa, 'rsnowsa')
      SFCVAR(rsnowsv, 'rsnowsv')
      SFCVAR(rst, 'rst')
      SFCVAR(rt, 'rt')
      SFCVAR(runofftot, 'runofftot')
      SFCVAR(runofftotaf, 'runofftotaf')
      SFCVAR(rveg, 'rveg')
      SFCVAR(sand, 'sand')
      SFCVAR(sanden, 'sanden')
      SFCVAR(sfcwgt, 'sfcwgt')
      SFCVAR(skin_depth, 'skin_depth')
      SFCVAR(skin_inc, 'skin_inc')
      SFCVAR(skyview, 'skyview')
      SFCVAR(slop, 'slop')
      SFCVAR(slope, 'slope')
      SFCVAR(snden, 'snden')
      SFCVAR(snoagen, 'snoagen')
      SFCVAR(snoal, 'snoal')
      SFCVAR(snoalen, 'snoalen')
      SFCVAR(snoden, 'snoden')
      SFCVAR(snodp, 'snodp')
      SFCVAR(snodpl, 'snodpl')
      SFCVAR(snoma, 'snoma')
      SFCVAR(snoro, 'snoro')
      SFCVAR(snowrate, 'snowrate')
      SFCVAR(snval, 'snval')
      SFCVAR(snvden, 'snvden')
      SFCVAR(snvdp, 'snvdp')
      SFCVAR(snvma, 'snvma')
      SFCVAR(snvro, 'snvro')  
      SFCVAR(sroad_alb, 'sroad_alb')
      SFCVAR(sroad_alben, 'sroad_alben')
      SFCVAR(sroad_emis, 'sroad_emis')
      SFCVAR(sroad_emisen, 'sroad_emisen')
      SFCVAR(sroad_rho, 'sroad_rho')
      SFCVAR(sroad_rhoen, 'sroad_rhoen')
      SFCVAR(sroad_t, 'sroad_t')
      SFCVAR(sroad_ten, 'sroad_ten')
      SFCVAR(sroad_ts, 'sroad_ts')
      SFCVAR(sroad_tsen, 'sroad_tsen')
      SFCVAR(sroad_wsnow, 'sroad_wsnow')
      SFCVAR(sroad_wsnowen, 'sroad_wsnowen')
      SFCVAR(sroof_alb, 'sroof_alb')
      SFCVAR(sroof_alben, 'sroof_alben')
      SFCVAR(sroof_emis, 'sroof_emis')
      SFCVAR(sroof_emisen, 'sroof_emisen')
      SFCVAR(sroof_rho, 'sroof_rho')
      SFCVAR(sroof_rhoen, 'sroof_rhoen')
      SFCVAR(sroof_t, 'sroof_t')
      SFCVAR(sroof_ten, 'sroof_ten')
      SFCVAR(sroof_ts, 'sroof_ts')
      SFCVAR(sroof_tsen, 'sroof_tsen')
      SFCVAR(sroof_wsnow, 'sroof_wsnow')
      SFCVAR(sroof_wsnowen, 'sroof_wsnowen')
      SFCVAR(stomr, 'stomr')
      SFCVAR(stomrvh, 'stomrvh')
      SFCVAR(stomrvl, 'stomrvl')
      SFCVAR(svf_road, 'svf_road')
      SFCVAR(svf_wall, 'svf_wall')
      SFCVAR(svs_wta, 'svs_wta')    
      SFCVAR(t_canyon, 't_canyon')
      SFCVAR(t_canyonen, 't_canyonen')
      SFCVAR(t_road, 't_road')
      SFCVAR(t_roaden, 't_roaden')
      SFCVAR(t_roof, 't_roof')
      SFCVAR(t_roofen, 't_roofen')
      SFCVAR(t_wall, 't_wall')
      SFCVAR(t_wallen, 't_wallen')
      SFCVAR(tc_road, 'tc_road')
      SFCVAR(tc_roaden, 'tc_roaden')
      SFCVAR(tc_roof, 'tc_roof')
      SFCVAR(tc_roofen, 'tc_roofen')
      SFCVAR(tc_wall, 'tc_wall')
      SFCVAR(tc_wallen, 'tc_wallen')
      SFCVAR(tddiagtyp, 'tddiagtyp')
      SFCVAR(tddiagtypv, 'tddiagtypv')
      SFCVAR(tdiag, 'tdiag')
      SFCVAR(tdiagstn, 'tdiagstn')
      SFCVAR(tdiagstnv, 'tdiagstnv')
      SFCVAR(tdiagtyp, 'tdiagtyp')
      SFCVAR(tdiagtypv, 'tdiagtypv')
      SFCVAR(tglacier, 'tglacier')
      SFCVAR(tground, 'tground')
      SFCVAR(thetaa, 'thetaa')
      SFCVAR(thetaap, 'thetaap')
      SFCVAR(ti_bld, 'ti_bld')
      SFCVAR(ti_blden, 'ti_blden')
      SFCVAR(ti_road, 'ti_road')
      SFCVAR(ti_roaden, 'ti_roaden')
      SFCVAR(tmice, 'tmice')
      SFCVAR(tnolim, 'tnolim')
      SFCVAR(tsa, 'tsa')
      SFCVAR(tsnavg, 'tsnavg')
      SFCVAR(tsnow, 'tsnow')
      SFCVAR(tsnowveg, 'tsnowveg')
      SFCVAR(tsoil, 'tsoil')
      SFCVAR(tsrad, 'tsrad')
      SFCVAR(tsradtw, 'tsradtw')
      SFCVAR(tss, 'tss')
      SFCVAR(tsun, 'tsun')
      SFCVAR(tsurf, 'tsurf')
      SFCVAR(tsvavg, 'tsvavg')
      SFCVAR(tve, 'tve')
      SFCVAR(tvege, 'tvege')
      SFCVAR(twater, 'twater')
      SFCVAR(u_canyon, 'u_canyon')
      SFCVAR(udiag, 'udiag')
      SFCVAR(udiagstn, 'udiagstn')
      SFCVAR(udiagstnv, 'udiagstnv')
      SFCVAR(udiagtyp, 'udiagtyp')
      SFCVAR(udiagtypv, 'udiagtypv')
      SFCVAR(urban, 'urban')
      SFCVAR(vdiag, 'vdiag')
      SFCVAR(vdiagstn, 'vdiagstn')
      SFCVAR(vdiagstnv, 'vdiagstnv')
      SFCVAR(vdiagtyp, 'vdiagtyp')
      SFCVAR(vdiagtypv, 'vdiagtypv')
      SFCVAR(vegf, 'vegf')
      SFCVAR(vegfrac, 'vegfrac')
      SFCVAR(vegh, 'vegh')
      SFCVAR(vegl, 'vegl')
      SFCVAR(vegtrans, 'vegtrans')
      SFCVAR(vgctem, 'vgctem')
      SFCVAR(wall_o_hor, 'wall_o_hor')
      SFCVAR(wall_o_horen, 'wall_o_horen')
      SFCVAR(watflow, 'watflow')
      SFCVAR(wfc, 'wfc')
      SFCVAR(wfcdp, 'wfcdp')
      SFCVAR(wfcint, 'wfcint')
      SFCVAR(wflux, 'wflux')
      SFCVAR(ws_road, 'ws_road')
      SFCVAR(ws_roaden, 'ws_roaden')
      SFCVAR(ws_roof, 'ws_roof')
      SFCVAR(ws_roofen, 'ws_roofen')
      SFCVAR(wsat, 'wsat')
      SFCVAR(wsnow, 'wsnow')
      SFCVAR(wsnv, 'wsnv')
      SFCVAR(wsoil, 'wsoil')
      SFCVAR(wsoilm, 'wsoilm')
      SFCVAR(wveg, 'wveg')
      SFCVAR(wwilt, 'wwilt')
      SFCVAR(xcent, 'xcent')
      SFCVAR(yradin, 'yradin')
      SFCVAR(yradsun, 'yradsun')
      SFCVAR(yradshade, 'yradshade')
      SFCVAR(yradrfsun, 'yradrfsun')
      SFCVAR(yradrfshade, 'yradrfshade')
      SFCVAR(yutciin, 'yutciin')
      SFCVAR(yutcisun, 'yutcisun')
      SFCVAR(yutcishade, 'yutcishade')
      SFCVAR(yutcirfsun, 'yutcirfsun')
      SFCVAR(yutcirfshade, 'yutcirfshade')
      SFCVAR(ywbgtsun, 'ywbgtsun')
      SFCVAR(ywbgtshade, 'ywbgtshade')
      SFCVAR(ywbgtrfsun, 'ywbgtrfsun')
      SFCVAR(ywbgtrfshade, 'ywbgtrfshade')
      SFCVAR(yutcicin, 'yutcicin')
      SFCVAR(yutcicsun, 'yutcicsun')
      SFCVAR(yutcicshade, 'yutcicshade')
      SFCVAR(yutcicrfsun, 'yutcicrfsun')
      SFCVAR(yutcicrfshade, 'yutcicrfshade')
      SFCVAR(ytglbsun, 'ytglbsun')
      SFCVAR(ytglbshade, 'ytglbshade')
      SFCVAR(ytglbrfsun, 'ytglbrfsun')
      SFCVAR(ytglbrfshade, 'ytglbrfshade')
      SFCVAR(ytwetb, 'ytwetb')
      SFCVAR(ytwetbrf, 'ytwetbrf')
      SFCVAR(ytrfzt, 'ytrfzt')
      SFCVAR(ytrdzt, 'ytrdzt')
      SFCVAR(yurdzu, 'yurdzu')
      SFCVAR(yQ1, 'yQ1')
      SFCVAR(yQ2, 'yQ2')
      SFCVAR(yQ3, 'yQ3')
      SFCVAR(yQ4, 'yQ4')
      SFCVAR(yQ5, 'yQ5')
      SFCVAR(yQ6, 'yQ6')
      SFCVAR(yQ7, 'yQ7')
      SFCVAR(yQ8, 'yQ8')
      SFCVAR(yQ9, 'yQ9')
      SFCVAR(yQ10, 'yQ10')
      SFCVAR(yQ11, 'yQ11')
      SFCVAR(yQ12, 'yQ12')
      SFCVAR(yQ13, 'yQ13')
      SFCVAR(z0, 'z0')
      SFCVAR(z0_road, 'z0_road')
      SFCVAR(z0_roaden, 'z0_roaden')
      SFCVAR(z0_roof, 'z0_roof')
      SFCVAR(z0_roofen, 'z0_roofen')
      SFCVAR(z0_town, 'z0_town')
      SFCVAR(z0_townen, 'z0_townen')
      SFCVAR(z0en, 'z0en')
      SFCVAR(z0ha, 'z0ha')
      SFCVAR(z0mvh, 'z0mvh')
      SFCVAR(z0mvl, 'z0mvl')
      SFCVAR(z0veg, 'z0veg')
      SFCVAR(z0t, 'z0t')
      SFCVAR(z0tveg, 'z0tveg')
      SFCVAR(za, 'za')
      SFCVAR(ze, 'ze')
      SFCVAR(zenith, 'zenith')
      SFCVAR(ztsl, 'ztsl')
      SFCVAR(zusl, 'zusl')

   end type SFCVARLIST_T

   integer, parameter :: INDX_SOIL    =  1
   integer, parameter :: INDX_GLACIER =  2
   integer, parameter :: INDX_WATER   =  3
   integer, parameter :: INDX_ICE     =  4
   integer, parameter :: INDX_AGREGE  =  5
   integer, parameter :: INDX_URB     =  6
   integer, parameter :: INDX_MAX     =  6

   type(SFCVARLIST_T), target  :: vd
   type(SFCVAR_T), allocatable :: vl(:)
   type(sfcptr), allocatable :: busptr(:)
   integer, allocatable :: statut(:,:)

   integer :: surfesptot = 0
   integer :: nvarsurf = 0  !# Number of surface bus var
   integer :: nsurf    = 0  !# Number of surface "types"
   integer :: tsrad_i=0, z0_i=0, z0t_i=0 !#TODO: remove, replace by vd%tsrad...

   integer :: drain=0
   integer :: drainaf=0
   integer :: insmavg=0
   integer :: isoil=0
   integer :: leg=0
   integer :: legaf=0
   integer :: ler=0
   integer :: leraf=0
   integer :: les=0
   integer :: lesaf=0
   integer :: letr=0
   integer :: letraf=0
   integer :: lev=0
   integer :: levaf=0
   integer :: overfl=0
   integer :: overflaf=0
   integer :: rootdp=0
   integer :: wflux=0
   integer :: wfluxaf=0
   integer :: wsoil=0

contains


   function sfcbus_init() result(F_istat)
      use phy_typedef, only: phymeta
      use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
      use sfc_options, only: schmurb
      implicit none
      integer :: F_istat

#include <msg.h>
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>

      integer :: i, istat, mulmax, idxmax
      type(SFCVAR_T) :: vl0(1)
      type(phymeta) :: mymeta
      type(phymetaplus) :: mymetaplus

      F_istat = RMN_ERR

      if (nsurf == 0) then
         idxmax = max(INDX_SOIL, INDX_GLACIER, INDX_WATER, INDX_ICE, INDX_AGREGE)
         if (schmurb /= 'NIL') idxmax = max(idxmax, INDX_URB)
         nsurf = idxmax - 1
      endif
 
      nvarsurf = size(transfer(vd, vl0))
      allocate(vl(nvarsurf))
      allocate(busptr(nvarsurf))
      vl = transfer(vd, vl)
      mulmax = 0
      do i = 1,nvarsurf
         vl(i)%i = i
         istat = clib_toupper(vl(i)%n)
         nullify(busptr(i)%ptr)
         istat = phygetmetaplus(mymetaplus, vl(i)%n, F_npath='V', &
              F_bpath='DPVE', F_quiet=.true., F_shortmatch=.false.)
         if (istat >= 0) then
            mymeta = mymetaplus%meta
            busptr(i)%ptr => mymetaplus%vptr
            vl(i)%doagg_L = (mymeta%bus(1:1) /= 'E')
            vl(i)%mul = mymeta%fmul
            vl(i)%niveaux = mymeta%nk
            vl(i)%mosaik = mymeta%mosaic + 1
            mulmax = max(mulmax, vl(i)%mul)
         endif

         select case(vl(i)%n)
         case('TSRAD')
            tsrad_i = i
         case('Z0')
            z0_i = i
         case('Z0T')
            z0t_i = i
         case('DRAIN')
            drain = mymetaplus%index
         case('DRAINAF')
            drainaf = mymetaplus%index
         case('INSMAVG')
            insmavg = mymetaplus%index
         case('ISOIL')
            isoil = mymetaplus%index
         case('LEG')
            leg = mymetaplus%index
         case('LEGAF')
            legaf = mymetaplus%index
         case('LER')
            ler = mymetaplus%index
         case('LERAF')
            leraf = mymetaplus%index
         case('LES')
            les = mymetaplus%index
         case('LESAF')
            lesaf = mymetaplus%index
         case('LETR')
            letr = mymetaplus%index
         case('LETRAF')
            letraf = mymetaplus%index
         case('LEV')
            lev = mymetaplus%index
         case('LEVAF')
            levaf = mymetaplus%index
         case('OVERFL')
            overfl = mymetaplus%index
         case('OVERFLAF')
            overflaf = mymetaplus%index
         case('ROOTDP')
            rootdp = mymetaplus%index
         case('WFLUX')
            wflux = mymetaplus%index
         case('WFLUXAF')
            wfluxaf = mymetaplus%index
         case('WSOIL')
            wsoil = mymetaplus%index
         end select

      enddo
      vd = transfer(vl, vd)
 
      allocate(statut(nvarsurf, mulmax))
      statut = 0

      F_istat = RMN_OK
      return
   end function sfcbus_init

end module sfcbus_mod
