.SUFFIXES :

.SUFFIXES : .tmpl90 .inc .ftn .ftn90 .hf .ptn .f .f90 .c .s .o

SUPP_OPT =

RMNLIB_BASE_PATH=$(ARMNLIB)

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

TARRLS = beta

OBJECTS = gmm.o

FICHIERS =  gmm_internals.cdk90 gmm_checkpoint.tmpl90 gmm_create.tmpl90 gmm_directory.ftn90 \
            gmm_encodemeta.ftn90 gmm_get.tmpl90 gmm_keys.ftn90 \
            gmm_pointer_table.tmpl90 gmm_shuffle.ftn90 gmm_rename.ftn90 gmm_getmeta.ftn90 \
            gmm_updatemeta.ftn90 gmm_update_pti.tmpl90 gmm_verbosity.ftn90 gmm_dumpinfo.ftn90 \
            gmm_definitions.inc gmm_internal_constants.inc gmm_nulls.inc

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

gmm_code.o: gmm_code.f90

libgmm_$(EC_ARCH).a: gmm_checkpoint.o gmm_create.o gmm_directory.o gmm_internals.o  gmm_get.o gmm_keys.o gmm_pointer_table.o gmm_code.tmpl90
	ar rcv libgmm_$(EC_ARCH).a gmm_checkpoint.o gmm_create.o gmm_directory.o gmm_internals.o  gmm_get.o gmm_keys.o gmm_pointer_table.o

gmm_internals.cdk90: gmm_internal_constants.inc gmm_directory_interface.inc

gmm_internals.o: gmm_internals.cdk90
#	s.compile -openmp -src gmm_internals.cdk90
	s.ftn90 -c -openmp -src gmm_internals.cdk90

gmm.hf: 	gmm_constants.hf gmm_definitions.inc gmm_nulls.inc \
			gmm_checkpoint_interface.inc gmm_create_interface.inc \
         gmm_encodemeta_interface.inc gmm_get_interface.inc gmm_getmeta_interface.inc gmm_keys_interface.inc \
         gmm_pointer_table_interface.inc gmm_rename_interface.inc gmm_shuffle_interface.inc \
         gmm_updatemeta_interface.inc \
         gmm_update_pti_interface.inc gmm_verbosity_interface.inc gmm_dumpinfo_interface.inc
	cat gmm_constants.hf gmm_definitions.inc gmm_nulls.inc gmm_checkpoint_interface.inc \
		gmm_create_interface.inc \
		gmm_encodemeta_interface.inc gmm_get_interface.inc gmm_getmeta_interface.inc \
		gmm_keys_interface.inc gmm_rename_interface.inc gmm_shuffle_interface.inc \
		gmm_updatemeta_interface.inc \
		gmm_verbosity_interface.inc gmm_dumpinfo_interface.inc > gmm.hf

# gmm.inc: gmm_constants.inc gmm_definitions.inc gmm_keys_interface.inc gmm_shuffle_interface.inc \
#          gmm_getmeta_interface.inc gmm_rename_interface.inc
# 	cat gmm_constants.inc gmm_definitions.inc gmm_keys_interface.inc gmm_getmeta_interface.inc \
# 	gmm_rename_interface.inc gmm_shuffle_interface.inc  > gmm.inc

gmm.ftn90: gmm.hf gmm_checkpoint.o gmm_create.o gmm_directory.o  gmm_encodemeta.o gmm_keys.o gmm_get.o \
       gmm_rename.o gmm_getmeta.o gmm_pointer_table.o gmm_shuffle.o gmm_updatemeta.o gmm_update_pti.o gmm_verbosity.o gmm_dumpinfo.o
	rm -f gmm.f90
	cat gmm_internals.cdk90 gmm_pointer_table.f90 gmm_checkpoint.f90 gmm_create.f90 gmm_directory.f90  gmm_encodemeta.f90 gmm_get.f90 gmm_keys.f90 \
	gmm_shuffle.f90 gmm_rename.f90 gmm_getmeta.f90 gmm_updatemeta.f90 gmm_update_pti.f90 gmm_verbosity.f90 gmm_dumpinfo.f90 > gmm.ftn90


gmm.o: gmm.ftn90

interface: gmm_internals.o gmm_pointer_table_interface.inc gmm_checkpoint_interface.inc gmm_create_interface.inc gmm_directory_interface.inc gmm_encodemeta_interface.inc gmm_getmeta_interface.inc

test: a.out
	a.out

gmm_checkpoint.o: gmm_checkpoint_interface.inc gmm_pointer_table.o
#	s.compile -openmp -src gmm_checkpoint.ftn90
	s.ftn90 -c -openmp -src gmm_checkpoint.ftn90

gmm_checkpoint_interface.inc: gmm_checkpoint.tmpl90
	r.tmpl90.ftn90 gmm_checkpoint
#	s.compile -openmp -P -src gmm_checkpoint_interface.cdk90
	s.GPPF90 gmm_checkpoint_interface.cdk90 > gmm_checkpoint_interface.f90
	mv -f gmm_checkpoint_interface.f90 gmm_checkpoint_interface.inc

gmm_create.o: gmm_create_interface.inc gmm_pointer_table.o
#	s.compile -openmp -src gmm_create.ftn90
	s.ftn90 -c -openmp -src gmm_create.ftn90

gmm_create_interface.inc: gmm_create.tmpl90
	r.tmpl90.ftn90 gmm_create
#	s.compile -openmp -P -src gmm_create_interface.cdk90
	s.GPPF90 gmm_create_interface.cdk90 > gmm_create_interface.f90
	mv -f gmm_create_interface.f90 gmm_create_interface.inc

gmm_directory_interface.inc: gmm_directory.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_directory.ftn90 > gmm_directory_interface.inc

gmm_dumpinfo_interface.inc: gmm_dumpinfo.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_dumpinfo.ftn90 > gmm_dumpinfo_interface.inc

gmm_directory.o: gmm_directory.ftn90

gmm_get_interface.inc: gmm_get.tmpl90
	r.tmpl90.ftn90 gmm_get
#	s.compile -openmp -P -src gmm_get_interface.cdk90
	s.GPPF90 gmm_get_interface.cdk90 > gmm_get_interface.f90
	mv -f gmm_get_interface.f90 gmm_get_interface.inc

gmm_get.o: gmm_get_interface.inc gmm_pointer_table.o
#	s.compile -openmp -src gmm_get.ftn90
	s.ftn90 -c -openmp -src gmm_get.ftn90

gmm_getmeta_interface.inc: gmm_getmeta.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_getmeta.ftn90 > gmm_getmeta_interface.inc

gmm_getmeta.o: gmm_getmeta.ftn90

gmm_keys_interface.inc: gmm_keys.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_keys.ftn90 > gmm_keys_interface.inc

gmm_keys.o: gmm_keys.ftn90

gmm_encodemeta_interface.inc: gmm_encodemeta.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_encodemeta.ftn90 > gmm_encodemeta_interface.inc

gmm_encodemeta.o: gmm_encodemeta.ftn90

gmm_rename_interface.inc: gmm_rename.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY< gmm_rename.ftn90 > gmm_rename_interface.inc

gmm_rename.o: gmm_rename.ftn90

gmm_getmeta.o: gmm_getmeta.ftn90

gmm_shuffle_interface.inc: gmm_shuffle.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_shuffle.ftn90 > gmm_shuffle_interface.inc

gmm_shuffle.o: gmm_shuffle.ftn90

gmm_pointer_table_interface.inc: gmm_pointer_table.tmpl90
	r.tmpl90.ftn90 gmm_pointer_table
#	s.compile -openmp -P -src gmm_pointer_table_interface.cdk90
	s.GPPF90 gmm_pointer_table_interface.cdk90 > gmm_pointer_table_interface.f90
	mv -f gmm_pointer_table_interface.f90 gmm_pointer_table_interface.inc

gmm_pointer_table.o: gmm_internals.o gmm_pointer_table_interface.inc
#	s.compile -openmp -src gmm_pointer_table.ftn90
	s.ftn90 -c -openmp -src gmm_pointer_table.ftn90

gmm_updatemeta_interface.inc: gmm_updatemeta.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_updatemeta.ftn90 > gmm_updatemeta_interface.inc

gmm_updatemeta.o: gmm_updatemeta.ftn90

gmm_update_pti_interface.inc: gmm_update_pti.tmpl90
	r.tmpl90.ftn90 gmm_update_pti
#	s.compile -openmp -P -src gmm_update_pti_interface.cdk90
	s.GPPF90 gmm_update_pti_interface.cdk90 > gmm_update_pti_interface.f90
	mv -f gmm_update_pti_interface.f90 gmm_update_pti_interface.inc

gmm_update_pti.o: gmm_update_pti_interface.inc gmm_pointer_table.o
#	s.compile -openmp -src gmm_update_pti.ftn90
	s.ftn90 -c -openmp -src gmm_update_pti.ftn90

gmm_verbosity_interface.inc: gmm_verbosity.ftn90
	r.gppf -f90 -gpp -DINTERFACEONLY < gmm_verbosity.ftn90 > gmm_verbosity_interface.inc

gmm_verbosity.o:	gmm_verbosity.ftn90

gmm_dumpinfo.o: gmm_dumpinfo.ftn90

tarball: $(FICHIERS) Makefile
	ls -al $(FICHIERS) Makefile tests
	tar cfzv /data/armnraid1/www/ssm/sources/gmm_$(TARRLS)_all.tgz $(FICHIERS) Makefile

obj: $(OBJECTS)

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir \
	$(TEMPLIB)(gmm.o)
	rm -f $(OBJECTS)

updlib:
#mettre a jour la programmatheque
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) gmm.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) gmm.o ; \
	fi

genlib: $(OBJECTS)
	r.ar  rcv $(MYLIB) $(OBJECTS)
	rm -f $(OBJECTS)

clean:
	rm -f *.mod *.o a.out
	rm -rf .fo

veryclean:
	rm -f gmm.ftn90 gmm.fh *.mod *.o a.out gmm_definitions.ftn90 gmm_code.ftn90
	rm -f *.f90 *interface.inc
	rm -f gmm_checkpoint.ftn90 gmm_create.ftn90 gmm_get.ftn90 gmm_pointer_table.ftn90 \
	      gmm_update_pti.ftn90
	rm -f libgmm_*.a gmm_restart gmm.inc 
	rm -rf .fo
