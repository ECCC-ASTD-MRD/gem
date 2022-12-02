#COMM_VERSION=301

TESTLIBS = $(LIBS)
TESTLIBSNOMPI = $(TESTLIBS) rpn_commstubs$(COMM_VERSION)
TESTRBUILD = $(RBUILD) -bidon -librmn $(RMN_VERSION) -libappl "$(TESTLIBSNOMPI)" -codebeta $(CODEBETA) -optf "=$(LFLAGS)" -obj *.o
TESTRBUILD2 =    cd $(LCLPO);$(RBUILD) $(OMP)      -bidon -libappl "$(TESTLIBSNOMPI)" -codebeta $(CODEBETA) -librmn $(RMN_VERSION) -optf "=$(LFLAGS)" -obj *.o
TESTRBUILD2MPI = cd $(LCLPO);$(RBUILD) $(OMP) -mpi -bidon -libappl "$(TESTLIBS)"      -codebeta $(CODEBETA) -librmn $(RMN_VERSION) -optf "=$(LFLAGS)" -obj *.o
CPL = cpl_stubs

test_str_split: testutils.o str_split.o test_str_split.o
	$(TESTRBUILD2) -main test_str_split -o test_str_split.Abs

test_fstinf_datev: test_fstinf_datev.o
	$(TESTRBUILD2) -main test_fstinf_datev -o test_fstinf_datev.Abs

test_fst: testutils.o fst.o test_fst.o
	$(TESTRBUILD2) -main test_fst -o test_fst.Abs

test_incfg: testutils.o str_split.o incfg.o test_incfg.o
	$(TESTRBUILD2) -main test_incfg -o test_incfg.Abs

test_input: testutils.o str_split.o fst.o incfg.o input.o test_input.o
	$(TESTRBUILD2) -main test_input -o test_input.Abs

test_input2:  test_input2.o
	$(TESTRBUILD2) -main test_input2 -o test_input2.Abs

test_input2_mpi:  test_input2_mpi.o
	$(TESTRBUILD2MPI) -main test_input2_mpi -o test_input2_mpi.Abs

test_time_interp: testutils.o time_interp.o test_time_interp.o
	$(TESTRBUILD2) -main test_time_interp -o test_time_interp.Abs

#test_ezgrid: testutils.o ezgrid.o test_ezgrid.o
test_ezgrid: test_ezgrid.o
	$(TESTRBUILD2) -main test_ezgrid -o test_ezgrid.Abs

#test_ezgrid_mpi: testutils.o ezgrid.o test_ezgrid_mpi.o
test_ezgrid_mpi: test_ezgrid_mpi.o
	$(TESTRBUILD2MPI) -main test_ezgrid_mpi -o test_ezgrid_mpi.Abs

test_inbloc_mpi: testutils.o test_inbloc_mpi.o
	$(TESTRBUILD2MPI) -main test_inbloc_mpi -o test_inbloc_mpi.Abs

tmp_print_list: tmp_create_lists.o tmp_print_list.o
	$(TESTRBUILD2) -main tmp_print_list -o tmp_print_list.Abs

#FFLAGS="'-C -g'"

test_file_utils.o: test_file_utils.ftn90 ../utils4tests.hf
test_handle_error.o: test_handle_error.ftn90 ../utils4tests.hf
test_handle_error2.o: test_handle_error2.ftn90 ../utils4tests.hf
test_stop_mpi.o: test_stop_mpi.ftn90 ../stop_mpi.h ../utils4tests.hf

copy_hgrid: copy_hgrid.o
	$(TESTRBUILD2) -main copy_hgrid -o copy_hgrid
