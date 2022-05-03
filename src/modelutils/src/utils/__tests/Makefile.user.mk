ifneq (,$(DEBUGMAKE))
$(info ## ====================================================================)
$(info ## File: Makefile.user.mk)
$(info ## )
endif

## For details RDE Makefile variables see:
## https://wiki.cmc.ec.gc.ca/wiki/RDE/1.0/Ref#Makefile_Vars.2C_Rules_and_targets

VERBOSE = 1
# OPTIL   = 2
# OMP     = -openmp
# MPI     = -mpi
# LFLAGS  =     # User's flags passed to the linker
# ifneq (,$(filter intel%,$(COMP_ARCH))$(filter PrgEnv-intel%,$(COMP_ARCH)))
# FFLAGS  = -C -g -traceback -ftrapuv #-warn all
# FFLAGS  = -warn unused
# CFLAGS  = -C -g -traceback -ftrapuv -fp-model precise #-Wall
# else
# FFLAGS  = -C -g -traceback
# CFLAGS  = -C -g -traceback 
# endif
# LIBAPPL = 
# LIBPATH_USER = /users/dor/armn/sch/Data/ords/big_tmp/libs/phy_6.0.a13/$(EC_ARCH)
# COMP_RULES_FILE = 
# PROFIL  = -prof   # set to -prof to enable the profiler (default = "")


## Optionally for gem the following options can also be modified
##
## For details GEM additional Makefile variables see:
## https://wiki.cmc.ec.gc.ca/wiki/GEM/4.8/dev'
## 
## Note: MODELUTILS_COMP_RULES_DEBUG is the Compiler_rules file
##       used to produce the debug libs

# BUILDNAME          = my_exp_name
# RMN_VERSION_NUMBER = 015.2
# COMM_VERSION       = _4051606
# VGRID_VERSION      = _$(VGRIDDESCRIPTORS_VERSION)
# COMP_RULES_FILE    = $(MODELUTILS_COMP_RULES_DEBUG)


## For GEM developpers:
## code developpement should mandatory be done with
## the following options: (uncomment the lines below)

# ifneq (,$(filter intel%,$(COMP_ARCH))$(filter PrgEnv-intel%,$(COMP_ARCH)))
# FFLAGS = -C -g -traceback -ftrapuv #-warn all
# else
# FFLAGS = -C -g -traceback
# endif
# MODELUTILS_SFX = -d
# RPNPHY_SFX = -d
# GEMDYN_SFX = -d


## To build with a local version of all libraries
## 
## Remove files (include or w/t module/sub/function) with: 
## rderm filename.ext
## do NOT create a stopping stub (prefer to catch at load time than run time)
##
## make dep         #mandatory
## make -j objects  #mandatory
## make modelutils_libs rpnphy_libs gemdyn_libs  #mandatory
## make allbin_gem # allbin_modelutils allbin_rpnphy allbin_gemdyn 

# MODELUTILS_VERSION=$(USER)
# RPNPHY_VERSION=a5a13
# GEMDYN_VERSION=$(USER)

ifneq (,$(ATM_MODEL_USERLIBS))
ifeq (,$(COMP_RULES_FILE))
ifeq (,$(wildcard $(HOME)/userlibs/$(EC_ARCH)/Compiler_rules))
ifneq (,$(wildcard $(ATM_MODEL_USERLIBS)/$(EC_ARCH)/Compiler_rules))
COMP_RULES_FILE = $(ATM_MODEL_USERLIBS)/$(EC_ARCH)/Compiler_rules
endif
endif
endif
endif

# maintest_vgrid_wb=test_vgrid_wb.Abs
# test_vgrid_wb_abs: | test_vgrid_wb_rm $(BINDIR)/$(maintest_vgrid_wb)
# 	ls -l $(BINDIR)/$(maintest_vgrid_wb)
# test_vgrid_wb_rm:
# 	rm -f $(BINDIR)/$(maintest_vgrid_wb)
# $(BINDIR)/$(maintest_vgrid_wb): | $(MODELUTILS_VFILES)
# 	export MAINSUBNAME="test_vgrid_wb" ;\
# 	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
# 	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
# 	export RBUILD_LIBAPPL="$(MODELUTILS_LIBS_V) $(MODELUTILS_LIBS_DEP)" ;\
# 	export RBUILD_COMM_STUBS=$(LIBCOMM_STUBS) ;\
# 	$(RBUILD4objNOMPI)

maintest_vgrid_wb=test_vgrid_wb.Abs
test_vgrid_wb_abs: | test_vgrid_wb_rm $(BINDIR)/$(maintest_vgrid_wb)
	ls -l $(BINDIR)/$(maintest_vgrid_wb)
test_vgrid_wb_rm:
	rm -f $(BINDIR)/$(maintest_vgrid_wb)
$(BINDIR)/$(maintest_vgrid_wb): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_vgrid_wb" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

maintest_timestr=test_timestr.Abs
test_timestr: | test_timestr_rm $(BINDIR)/$(maintest_timestr)
	ls -l $(BINDIR)/$(maintest_timestr)
test_timestr_rm:
	rm -f $(BINDIR)/$(maintest_timestr)
$(BINDIR)/$(maintest_timestr): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_timestr" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(MODELUTILS_LIBS_V) $(MODELUTILS_LIBS_DEP)" ;\
	export RBUILD_COMM_STUBS=$(LIBCOMM_STUBS) ;\
	$(RBUILD4objNOMPI)

maintest_fst=test_fst.Abs
test_fst: | test_fst_rm $(BINDIR)/$(maintest_fst)
	ls -l $(BINDIR)/$(maintest_fst)
test_fst_rm:
	rm -f $(BINDIR)/$(maintest_fst)
$(BINDIR)/$(maintest_fst): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_fst" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

maintest_fstmpio=test_fstmpio.Abs
test_fstmpio: | test_fstmpio_rm $(BINDIR)/$(maintest_fstmpio)
	ls -l $(BINDIR)/$(maintest_fstmpio)
test_fstmpio_rm:
	rm -f $(BINDIR)/$(maintest_fstmpio)
$(BINDIR)/$(maintest_fstmpio): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_fstmpio" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

maintest_incfg2=test_incfg2.Abs
test_incfg2: | test_incfg2_rm $(BINDIR)/$(maintest_incfg2)
	ls -l $(BINDIR)/$(maintest_incfg2)
test_incfg2_rm:
	rm -f $(BINDIR)/$(maintest_incfg2)
$(BINDIR)/$(maintest_incfg2): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_incfg2" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

maintest_input7=test_input7.Abs
test_input7: | test_input7_rm $(BINDIR)/$(maintest_input7)
	ls -l $(BINDIR)/$(maintest_input7)
test_input7_rm:
	rm -f $(BINDIR)/$(maintest_input7)
$(BINDIR)/$(maintest_input7): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_input7" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

maintest_inputio7=test_inputio7.Abs
test_inputio7: | test_inputio7_rm $(BINDIR)/$(maintest_inputio7)
	ls -l $(BINDIR)/$(maintest_inputio7)
test_inputio7_rm:
	rm -f $(BINDIR)/$(maintest_inputio7)
$(BINDIR)/$(maintest_inputio7): | $(MODELUTILS_VFILES)
	export MAINSUBNAME="test_inputio7" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME}" ;\
	export ATM_MODEL_VERSION="$(MODELUTILS_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

## Sample MPI abs targets

mympiabsname: | mympiabsname_rm $(BINDIR)/mympiabsname
	ls -l $(BINDIR)/mympiabsname
mympiabsname_rm:
	rm -f $(BINDIR)/mympiabsname
$(BINDIR)/mympiabsname: my_main_sub_name.o
	export MAINSUBNAME="my_main_sub_name" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME} $(BUILDNAME)" ;\
	export ATM_MODEL_VERSION="$(GEM_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	$(RBUILD4objMPI)

## Sample no-MPI abs targets

myabsname: | myabsname_rm $(BINDIR)/myabsname
	ls -l $(BINDIR)/myabsname
myabsname_rm:
	rm -f $(BINDIR)/myabsname
$(BINDIR)/myabsname: my_main_sub_name.o
	export MAINSUBNAME="my_main_sub_name" ;\
	export ATM_MODEL_NAME="$${MAINSUBNAME} $(BUILDNAME)" ;\
	export ATM_MODEL_VERSION="$(GEM_VERSION)" ;\
	export RBUILD_LIBAPPL="$(GEMDYN_LIBS_V) $(GEMDYN_LIBS_DEP)" ;\
	export RBUILD_COMM_STUBS=$(LIBCOMM_STUBS) ;\
	$(RBUILD4NOMPI)


ifneq (,$(DEBUGMAKE))
$(info ## ==== Makefile.user.mk [END] ========================================)
endif

LIBHPCSPERF =
