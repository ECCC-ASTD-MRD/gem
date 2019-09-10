
ifeq ($(BASE_ARCH),$(EC_ARCH))
   $(info ERROR: please setup the compilation environment by typing the folowing:)
   $(info . setup.dot or . setup_gpsc.dot)
   $(error )
endif

ifeq ($(USE_EFENCE),1)
   LEFENCE=-lefence
endif

include version.dot

default: libdescrip.a

dependencies.mk: *.c *.F90 *.h *.hc
	ls *.hc *.c *.F90 | ./s.dependencies.pl > dependencies.mk ;\

include dependencies.mk

%.o: %.F90
	s.f90 -c -openmp $(MY_FFLAGS) $(FFLAGS) $(LEFENCE) $<

%.o: %.c
	s.cc -c -openmp $(MY_CFLAGS) $(CFLAGS) $(LEFENCE) $<

libdescrip.a: $(OBJECTS) dependencies.mk vgrid_version
	rm -f libdescrip.a && ar scru libdescrip.a *.o ;\

shared: libdescrip.a
	if [ $$ORDENV_PLAT != ubuntu-12.04-amd64-64 -a $$ORDENV_PLAT != aix-7.1-ppc7-64 ]; then \
	   s.f90 -shared -openmp -o libdescripshared.so *.o -lrmnshared -Wl,-rpath,$(RPN_SHARED_LIB);\
	else \
	   echo "No shared library on $$ORDENV_PLAT anymore due to new incompatible flag -static-intel";\
	fi

vgrid_version:
	echo "character(len=128) :: vgrid_descriptors_version = \"VGRID $(BH_PULL_SOURCE_GIT_BRANCH) $(COMP_ARCH) $(ORDENV_PLAT) $(shell date)\"" > vgrid_version.hf ;\
	echo                "char vgrid_descriptors_version[] = \"VGRID $(BH_PULL_SOURCE_GIT_BRANCH) $(COMP_ARCH) $(ORDENV_PLAT) $(shell date)\";" > vgrid_version.hc ;\

ssm:
	../bh_vgrid.py -p $(ORDENV_PLAT) -m $(BH_MODE) -v BH_PULL_SOURCE=$(PWD)/.. -v BH_PULL_SOURCE_GIT_BRANCH=$(BH_PULL_SOURCE_GIT_BRANCH) -v BH_PACKAGE_VERSION=$(BH_PULL_SOURCE_GIT_BRANCH)-$(COMP_ARCH) --local ;\

clean:
	rm -f *.o *.mod dependencies.mk;\

clean_ssm: clean_bh
	rm -f *.ssm ;\

clean_bh:
	rm -rf vgriddescriptors vgriddescriptors_*  log-*-*-*-*-* report-*-*-*-*-* work-*-*-*-*-* ;\

distclean: clean
	rm -f libdescrip.a libdescripshared.so ;\

clean_tilde:
	rm -f *~ ;\

clean_all: distclean clean_ssm clean_bh clean_tilde
