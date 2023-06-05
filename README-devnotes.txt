TODO: remove bus access from lower level params
## egrep -Ri '\bv\(' $rpnphy/src | grep -v ':!' | cut -d: -f1| cut -d/ -f12 | sort -u

#### Look for unused subroutine
## cd rpnphy/src/ ; for item in $(grep -i subroutine [a-z]*/* | cut -d: -f2 | cut -d'(' -f1 | grep -vi 'end subroutine' | egrep -v '!' | grep -v '"'| grep -i subroutine | sed 's/SUBROUTINE//ig') ; do echo === $item : $(egrep -i "\b$item\b" [a-z]*/*.[fF]* [a-z]*/*.cdk90 2>/dev/null | grep -i call | egrep -v ':!' | egrep -v ':\s+!' | egrep -i "\b$item\b" | wc -l) ; done | grep ': 0'
=== phy_typedef_dummy : 0 ok
=== mp_p3_wrapper_wrf : 0  part of base/mp_p3
=== SEDI_main_1 : 0 part of base/my_sedi_mod
=== TTTT : 0 commented call in base/sun7
=== FESERI : 0 main
=== prphynml : 0 main
=== test_sfclayer : 0 main
=== testphy_main : 0 main
=== dumpini3 : 0 utils
=== statvps : 0 utils
=== statvps1 : 0 utils
=== COUPLING_TEB : 0 ??
=== SURFACE_CD : 0 ??
=== TEB : 0 ??
=== URBAN_DRAG : 0 ??
Done


#### Look for not included cdk
## for item in *.cdk ; do echo $item : $(grep $item *.[fF]* *.cdk90 2>/dev/null | grep -i include | egrep -v ':!' | egrep -v ':\s+!' | wc -l) ; done | grep ': 0'
Done

#### Look for surface.cdk var not used in base
for item in $(cat rpnphy/src/base/surface.cdk | grep -i common | cut -d/ -f3) ; do echo $item : $(egrep -i "\b${item}\b" rpnphy/src/base/*.F90 | egrep -v ':\s+!' | egrep -v ':!' | wc -l); done
Done

#### Look for gesdict surface var not used in base
for item in $(grep -i lisba rpnphy/src/base/phybusinit.F90  | grep -i gesdict | cut -d, -f3); do echo ==== $item ; egrep -i "\b${item}\b" rpnphy/src/base/* | egrep -v ':\s+!' | egrep -v ':!' | grep -v phybus.F90 | grep -v phybusinit.F90 ; done
Done

#### Look for gesdict surface var not used in base/convect
for item in $(grep -i gesdict rpnphy/src/base/phybusinit.F90  | egrep -v '^\s+!' | cut -d, -f3); do echo ==== $item : $(egrep -i "\b${item}\b" rpnphy/src/{api,base,convect,series,utils}/* | egrep -v ':\s+!' | egrep -v ':!' | grep -v phybus.F90 | grep -v phybusinit.F90 | wc -l) : $(egrep -i "\b${item}\b" rpnphy/src/surface/* | egrep -v ':\s+!' | egrep -v ':!' | grep -v phybus.F90 | grep -v phybusinit.F90 | wc -l) | grep ': 0 :'; done

mrkv : 0 : 0 #ens_marfield_cg.F90:      istat = phy_put
runoff #hydro1, teb, town, urban
Done

#### Look for surface.cdk var not used in base

for item in $(cat rpnphy/src/surface/sfcbus.cdk | grep DCL_PHYVAR | cut -d'(' -f2 | cut -d, -f1) ; do echo ==== $item : .... | grep -v ': 0 : 0 : 0 : 0' ; done | grep ': 0'
# tot : tot no xst : tot phybusini ; tot phybus 
==== umoins : 11 : 11 : 1 : 1
==== uplus : 20 : 19 : 1 : 1
==== vmoins : 11 : 11 : 1 : 1
==== vplus : 20 : 19 : 1 : 1
==== tmoins : 28 : 28 : 1 : 1
==== tplus : 36 : 34 : 1 : 1
==== sigm : 26 : 26 : 1 : 2
==== humoins : 25 : 25 : 1 : 1
==== huplus : 27 : 27 : 1 : 1
==== alfaq : 4 : 4 : 1 : 1
==== alfat : 6 : 6 : 1 : 1
==== alvis : 6 : 5 : 1 : 1
==== bm : 9 : 9 : 1 : 1
==== bt : 13 : 13 : 1 : 1
==== cosz : 4 : 4 : 1 : 1
==== dhdx : 3 : 3 : 1 : 1     #gwd, inicham/inisurf
==== dhdxdy : 3 : 3 : 1 : 1   #gwd, inicham/inisurf
==== dhdxdyen : 2 : 2 : 1 : 1 #gwd, inicham/inisurf
==== dhdxen : 2 : 2 : 1 : 1   #gwd, inicham/inisurf
==== dhdy : 3 : 3 : 1 : 1     #gwd, inicham/inisurf
==== dhdyen : 2 : 2 : 1 : 1   #gwd, inicham/inisurf
==== dlat : 39 : 38 : 1 : 1
==== dlon : 22 : 21 : 1 : 1
==== en : 55 : 47 : 3 : 1
==== epstfn : 3 : 3 : 1 : 1 #difver, pbl_difver, inisurf
==== fc : 34 : 14 : 3 : 1
==== fcor : 14 : 14 : 1 : 1
==== fdsi : 6 : 6 : 2 : 1
==== fdss : 11 : 5 : 3 : 1
==== fl : 13 : 9 : 1 : 1
==== flusolis : 10 : 4 : 1 : 1
==== fluslop : 4 : 2 : 1 : 1
==== fq : 8 : 4 : 3 : 1
==== frv : 10 : 10 : 1 : 1
==== ftemp : 8 : 8 : 1 : 1
==== fv : 26 : 6 : 3 : 1
==== fvap : 8 : 8 : 1 : 1
==== gc : 6 : 6 : 3 : 3
==== glsea : 5 : 4 : 1 : 1
==== gztherm : 5 : 5 : 1 : 1
==== h : 137 : 133 : 1 : 1
==== hst : 4 : 4 : 1 : 1
==== ilmo : 36 : 36 : 1 : 1
==== kcl : 14 : 14 : 1 : 1
==== km : 27 : 23 : 5 : 1
==== kt : 24 : 20 : 1 : 1
==== lhtg : 7 : 6 : 1 : 1   #gwd, inicham/inisurf
==== lhtgen : 2 : 2 : 1 : 1 #gwd, inicham/inisurf
==== mg : 59 : 58 : 2 : 1
==== ml : 34 : 34 : 1 : 1
==== mtdir : 12 : 12 : 1 : 1
==== pmoins : 80 : 50 : 1 : 1
==== qdiag : 6 : 6 : 1 : 1 #calcdiag, boundary_layer
==== qsurf : 5 : 5 : 1 : 1 #calcdiag, ccc*, ctmdiag, difver, pbl_difver
==== rainrate : 24 : 24 : 1 : 1
==== rib : 18 : 18 : 1 : 1
==== rt : 3 : 1 : 2 : 1     #calcdiag, extdiag... not in surface!
==== scl : 17 : 13 : 1 : 1
==== slope : 47 : 47 : 2 : 1
==== snowrate : 24 : 24 : 1 : 1
==== tdiag : 9 : 9 : 1 : 1
==== thetaa : 2 : 2 : 1 : 1 #cmtdiag
==== tsrad : 11 : 9 : 1 : 1
==== tss : 24 : 22 : 1 : 1
==== tsurf : 9 : 9 : 1 : 1
==== tve : 55 : 55 : 1 : 1
==== udiag : 6 : 6 : 1 : 1
==== vdiag : 6 : 6 : 1 : 1
==== xcent : 10 : 10 : 1 : 1
==== z0 : 50 : 49 : 1 : 1
==== z0t : 21 : 21 : 1 : 1
==== za : 6 : 6 : 1 : 1
==== ze : 74 : 67 : 1 : 1
==== ztsl : 12 : 12 : 1 : 1
==== zusl : 3 : 3 : 1 : 1 #calcdiag, ccc*, 


SERIES:
*TODO: serxst list base for all gesdict var
*TODO: mv serxst of surface from ccdiag to diagnosurf


DIAG:
*TODO: see paulV email about duplicate var
*TODO: calcdiag list base
*TODO: mv condensation accumul to calcdiag (?cnv_main modify them afterward?)


SURFACE API:
*TODO: itf_sfc_main, itf_cpl_init: replace phygrd.cdk needed for cpl ijdrv_phy
*TODO: rename files
		sfc_exch -> sfcexch
		sfc_debu -> sfc_init
*TODO: formalize sfc_ API 
		sfc_nml
      sfc_init     : was sfc_debu
		sfc_input?
		sfc_stepinit : was inichamp or inisurf
		sfc_step     : was itf_sfc_main
      wb_put/get   : to document list
      ?sfc_businit?
*TODO: move gesdict sfc call to surface dir s/r
TODO: remove mixing of base/surf init in inichamps 
		(remove non physics API calls from surface)
		radcons (from inichamps to phystepinit?)
		equivmount

OTHER:
*TODO: phy nml should have option no phy (with still no nml == no phy but if nml is there and empty use phy)

*TODO: update testphy for new API (COMPAT, TDMASK)

TODO: use vardict/gmmx instead of gesdict
      Bridge vardict/gmmx to gesdict during transition

TODO: rm tttt.F90 ? commented call in base/sun7

TODO: make sure phy_getmeta always return same ptr shape (folded) to sfc s/r (option?)

TODO: apply tendences consistency, make after param



phyexe.ftn90

+ inichamp.ftn90
. + inisurf.ftn90
. . + initown.ftn90

+ phystepinit1

+ radiation2
. + prep_cw_rad.ftn90
. + diagno_cw_rad.ftn90
. + cccmarad.ftn90 (cccmarad_ptr_as.cdk)
. + newrad.ftn90
. . + radir.ftn90
. + radslop.ftn90

+ itf_sfc_main
. + copybus.ftn90
. + isba/ebudget.ftn90
. + agrege.ftn90
. + diagnosurf.ftn90

+ metox2

+ gwd.ftn90
. + sgoflx.ftn90

+ apply_tendencies1

+ turbulence.ftn90
. + boundary_layer.ftn90

. . + pbl_turbul.ftn90
. . . + eturbl.ftn90
. . . . + absdvdz3.ftn90

. . + pbl_difver.ftn90

. . + apply_tendencies1

. . + turbul.ftn90
. . . + moistke.ftn90
. . . . + blcloud3.ftn90
. . . + eturbl.ftn90
. . . . + absdvdz3.ftn90
. . . + twind
. . . . + windgust.ftn90


. + boundary_layer_modlevs2
. . + turbul
. . . (...)
. . + difver.ftn90
. . . + ctmdiag.ftn90
. . + apply_tendencies1

+ shallconv3

+ precipitation
. + condensation.ftn90 (cnd_ptr_as.cdk)
. . + mp_my2_mod.f90
. . + my_dmom_mod.cdk90
. . ....
. . + ccdiagnostics
. + itf_cnv_main.ftn90 (cnv_ptr_as.cdk)

+ prep_cw.ftn90
+ tendency4
+ ens_ptp1
+ chm_exe

+ calcdiag (calcdiag_ptr_as.cdk)
. + refractivity.ftn90
. + lightning.ftn90

+ extdiag

testphy_phyexe.ftn90
dumpbus.ftn90       #dumpini3 not called

NOTE: 
* sfcbus  has index/order in the table
  pointeurs(indx) = indx
  ptdebut(indx)   = mymeta%index

phy_nml
		  phy_options_init
		  sfc_nml
					 sfc_options_init
					 cpl_nml
					 sl_put
		  cnv_nml
		  check_options
		  ser_nml
		  chm_nml
		  
phy_init
		  phydebu
					 litozon
					 litblrad
					 TABULATE_XCW
					 READ_ISCCPDATA
					 p3_init
					 phybusinit
					 #gmm_create buses
		  sfc_debu(sfc_getindx)
					 iniptsurf(sfc_getindx)  [sfcbus, DCLCHAR]
					 			external block data sfcbus_data [DATACHAR]
								phy_getmeta
					 sfc_exch_options
					 block data sfc_debu_data
		  mapping2drivergrid
		  ser_init5
		  itf_cpl_init
					 cpl_init


phy_input
		  phyfillbus
					 phyfold2
		  phyinputdiag #init
		  			 phy_getmeta
		  priv_ozone
					 intozon
		  chm_load_emissions
		  phy_getmeta
		  priv_fold
					 physimple_transforms3d
					 pre_fold_opr_clbk
					 statfld_dm
					 phyfold2
		  priv_checklist
					 phy_getmeta
		  phyent2per

phy_step
		  phydgflt
		  ser_ctrl
		  phyput_input_param
		  sfc_get_input_param
		  cpl_step
		  !$omp parallel
		  physlb1
					 testphy_phyexe
					 phyexe
								inichamp3  [sfcbus]
										  inisurf3  [sfcbus]
													 lacs3  [sfcbus]
													 inicover  [sfcbus]
													 coherence2  [sfcbus]
													 inisoili  [sfcbus]
													 initown  [sfcbus]
													 radcons
										  equivmount
								phystepinit1
										  sigmalev
										  serget
										  sersetm
										  mzoniv
										  tothermo
										  calcz0
										  phy_getmeta
										  phybusindx
										  mfotvt
										  surf_precip
										  surf_precip3
								radiation2
								itf_sfc_main  [sfcbus]
										  cpl_update
										  copybus2(bus_soil,.true.)  [sfcbus]
										  isba3  [sfcbus]
										  			fillagg   [sfcbus]
										  copybus2(bus_soil,.false.)
										  copybus2(bus_water,.true.)
										  water1  [sfcbus]
										  			fillagg   [sfcbus]
										  copybus2(bus_water,.false.)
										  copybus2(bus_ice,.true.)
										  seaice2  [sfcbus]
										  			fillagg   [sfcbus]
										  copybus2(bus_ice,.false.)
										  copybus2(bus_glacier,.true.)
										  glaciers1  [sfcbus]
										  			fillagg   [sfcbus]
										  copybus2(bus_glacier,.false.)
										  copybus2(bus_urb,.true.)
										  town  [sfcbus]
										  			fillagg   [sfcbus]
										  copybus2(bus_urb,.false.)
										  agrege1 [sfcbus]
										  diagnosurf  [sfcbus]
								metox2
								gwd8
								apply_tendencies1
								turbulence1
								shallconv3
								precipitation
								prep_cw2
								tendency4
								ens_ptp1
								chm_exe
								calcdiag
								extdiag
		  !$omp end parallel
		  ser_out
		  phydgflt
