include $(VPATH)/RPN_COMM_version.inc
$(info this is RPN_COMM version $(RPN_COMM_version_s))

# general building rules
include $(VPATH)/Makefile.common
# sources specific (mechanically generated) dependencies
include $(VPATH)/dependencies.mk

LIB      = rpn_comm$(SERIAL)
CLEAN    = rpn_comm_fortran_stubs.F90 rpn_comm_c_stubs.c mpi_stub.h \
           $(STUB_LIBRARY) $(LIBRARY) \
           $(VPATH)/rpn-comm_$(RPN_COMM_version_s)_multi.ssm $(VPATH)/RPN_COMM_interfaces.inc \
           $(VPATH)/RPN_COMM_interfaces_int.inc $(VPATH)/dependencies.mk $(VPATH)/dependencies+.mk

CLEANDIRS= $(VPATH)/rpn-comm_$(RPN_COMM_version_s)_multi $(LIBDIR)
TESTS    = TEST_000.Abs TEST_001.Abs TEST_002.Abs TEST_004.Abs \
           TEST_005.Abs TEST_006.Abs TEST_007.Abs TEST_008.Abs TEST_009.Abs \
           TEST_010.Abs TEST_011.Abs

FMODULES = RPN_COMM_mod.o
MPI_VERSION = $(shell $(VPATH)/../tools/mpi_version.sh $(VTAG) )
LIBNAME  = $(LIB)_$(RPN_COMM_version)$(MPI_VERSION)
LIBRARY  = $(LIBDIR)/lib$(LIBNAME).a
STUB_LIBRARY = $(LIBDIR)/lib$(LIB)stubs_$(RPN_COMM_version)$(MPI_VERSION).a
SOURCES  = $(INCDECKS) $(CDECKS) $(FDECKS) $(HDECKS) $(F90DECKS)

DISTINCLUDES = $(VPATH)/RPN_COMM_interfaces.inc $(VPATH)/RPN_COMM.inc $(VPATH)/rpn_comm.inc \
               $(VPATH)/RPN_COMM_types.inc $(VPATH)/RPN_COMM_constants.inc \
               $(VPATH)/RPN_COMM_ftoc.inc $(VPATH)/RPN_COMM_is_null.inc

lib: $(LIBRARY)

obj: $(OBJECTS)

all:  itf inc lib tests

full:  itf dep $(LIBRARY) $(TESTS) $(STUB_LIBRARY) $(LIBRARY).inc

tests:	$(TESTS)

stublib: $(STUB_LIBRARY)

#lib: $(LIBRARY)

inc: $(LIBRARY).inc

itf: $(VPATH)/RPN_COMM_interfaces.inc $(VPATH)/RPN_COMM_interfaces_int.inc

$(VPATH)/RPN_COMM_ptr.F90: $(VPATH)/../tools/gen_rpn_comm_ptr.sh
	(cd $(VPATH) ; ../tools/gen_rpn_comm_ptr.sh >$(VPATH)/RPN_COMM_ptr.F90)

$(VPATH)/RPN_COMM_interfaces.inc: $(wildcard $(VPATH)/RPN_COMM*.?90) $(wildcard $(VPATH)/*.c)
	(cd $(VPATH) ; cat RPN_COMM_*.?90 | grep '!InTfX!' | sed 's/^!![ ]*/      /' >RPN_COMM_interfaces.inc )
	(cd $(VPATH) ; cat RPN_COMM_*.?90 RPN_COMM_*.c | ../tools/extract_interface.sh >>RPN_COMM_interfaces.inc ; rm -f ../tools/wrap_code.exe)

$(VPATH)/RPN_COMM_interfaces_int.inc: $(wildcard $(VPATH)/RPN_COMM*.?90) $(wildcard $(VPATH)/*.c)
	(cd $(VPATH) ; rm -f RPN_COMM_interfaces_int.inc;)
	(cd $(VPATH) ; for target in RPN_COMM_*.?90; do grep -q '!InTfX!' $$target || continue ; ( echo "#if ! defined(IN_$${target%.*})" ; cat $$target | grep '!InTfX!' | sed 's/^!![ ]*/      /' ; echo "#endif" ) >>RPN_COMM_interfaces_int.inc ; done;)
	(cd $(VPATH) ; for target in RPN_COMM_*.?90 RPN_COMM_*.c ; do ../tools/extract_interface.sh $$target >>RPN_COMM_interfaces_int.inc ; done; rm -f ../tools/wrap_code.exe)

#dep: $(VPATH)/dependencies.mk

#$(VPATH)/dependencies.mk: $(wildcard $(VPATH)/*.f) $(wildcard $(VPATH)/*.f90) $(wildcard $(VPATH)/*.F90) $(wildcard $(VPATH)/*.c) \
#                          $(wildcard $(VPATH)/*.h) $(wildcard $(VPATH)/*.inc) \
#                          $(VPATH)/RPN_COMM_interfaces_int.inc $(VPATH)/RPN_COMM_interfaces.inc
$(VPATH)/dependencies.mk:
	rm -f $(VPATH)/dependencies.mk $(TMPDIR)/dependencies+.mk
#	-which gnu_find 2>/dev/null 1>/dev/null || (cd $(VPATH) ; find . -maxdepth 1 -type f | grep RPN_COMM | grep -v TEST_0 | ../tools/mk.dependencies.pl >$(TMPDIR)/dependencies+.mk ) || true
#	-(which gnu_find 2>/dev/null 1>/dev/null && (cd $(VPATH) ; gnu_find . -maxdepth 1 -type f | grep RPN_COMM | grep -v TEST_0 ../tools/mk.dependencies.pl >$(TMPDIR)/dependencies+.mk )) || true
	(cd $(VPATH) ; ../tools/rdedep.pl --flat_layout --out=$(TMPDIR)/dependencies+.mk $$(ls -1 RPN_COMM* | grep -v TEST_0)) || true
	mv $(TMPDIR)/dependencies+.mk $(VPATH)/dependencies.mk
	echo === touch  $(VPATH)/dependencies.mk ===
	touch  $(VPATH)/dependencies.mk
#
dep_rm:
	rm -f $(VPATH)/dependencies.mk
dep: dep_rm $(VPATH)/dependencies.mk 
#	rm -f $(VPATH)/dependencies.mk $(TMPDIR)/dependencies+.mk
#	-which gnu_find 2>/dev/null 1>/dev/null || (cd $(VPATH) ; find . -maxdepth 1 -type f | grep RPN_COMM | grep -v TEST_0 | ../tools/mk.dependencies.pl >$(TMPDIR)/dependencies+.mk ) || true
#	-(which gnu_find 2>/dev/null 1>/dev/null && (cd $(VPATH) ; gnu_find . -maxdepth 1 -type f | grep RPN_COMM | grep -v TEST_0 ../tools/mk.dependencies.pl >$(TMPDIR)/dependencies+.mk )) || true
#	mv $(TMPDIR)/dependencies+.mk $(VPATH)/dependencies.mk

ssm-package:
	rm -rf $(VPATH)/rpn-comm_${RPN_COMM_version_s}_multi
	(cd $(VPATH) ; tar zxf ssmtemplate_1.0_all.ssm ; mv ssmtemplate_1.0_all rpn-comm_$(RPN_COMM_version_s)_multi )
	(cd $(VPATH) ; cp rpn_comm_stubs.sh $(SOURCES) rpn-comm_$(RPN_COMM_version_s)_multi/src/.)
	(cd $(VPATH) ; cp Makefile Makefile.common Makefile.default Makefile.ECsetup *.mk \
	    rpn-comm_$(RPN_COMM_version_s)_multi/src/.)
	(cd $(VPATH) ; tar zcf rpn-comm_$(RPN_COMM_version_s)_multi.ssm rpn-comm_$(RPN_COMM_version_s)_multi)

rpn_comm_fortran_stubs.o: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh fortran

rpn_comm_c_stubs.o: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh c

$(STUB_LIBRARY): rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	mkdir -p $(LIBDIR)
	ar rcv $(STUB_LIBRARY) rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	(cd $(LIBDIR) ; ln -sf lib$(LIB)stubs_$(RPN_COMM_version)$(MPI_VERSION).a lib$(LIB)stubs$(MPI_VERSION).a)

$(LIBRARY): $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY)_ RPN_COMM*.o
#	ar rcv $(LIBRARY)_ $(OBJECTS)
#	ar dv $(LIBRARY)_ TEST_stubs.o rpn_comm_c_stubs.o rpn_comm_fortran_stubs.o
	mv $(LIBRARY)_ $(LIBRARY)
	ar t $(LIBRARY) | sort -u >$(VPATH)/objects.lst
	(cd $(LIBDIR) ; ln -sf lib$(LIBNAME).a  lib$(LIB)$(MPI_VERSION).a)
#	cp *.mod $(INCDIR)

checkref:
	sort -u $(VPATH)/REFERENCE.lst >$(VPATH)/SORTED.lst
	diff $(VPATH)/SORTED.lst $(VPATH)/objects.lst

#.PHONY:	$(VPATH)/includes

#$(VPATH)/includes: $(DISTINCLUDES)
$(LIBRARY).inc: $(DISTINCLUDES)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY).inc $(DISTINCLUDES)
	mkdir -p $(INCDIR)
	cp $(DISTINCLUDES) $(INCDIR)

TEST_000.Abs: $(LIBRARY) TEST_000.o
TEST_001.Abs: $(LIBRARY) TEST_001.o
TEST_002.Abs: $(LIBRARY) TEST_002.o
TEST_003.Abs: $(LIBRARY) TEST_003.o
TEST_004.Abs: $(LIBRARY) TEST_004.o
TEST_005.Abs: $(LIBRARY) TEST_005.o
TEST_006.Abs: $(LIBRARY) TEST_006.o
TEST_007.Abs: $(LIBRARY) TEST_007.o
TEST_008.Abs: $(LIBRARY) TEST_008.o
TEST_009.Abs: $(LIBRARY) TEST_009.o
TEST_010.Abs: $(LIBRARY) TEST_010.o
TEST_011.Abs: $(LIBRARY) TEST_011.o
