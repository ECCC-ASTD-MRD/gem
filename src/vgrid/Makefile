.SUFFIXES:

.SUFFIXES : .o .F90 .c .a

SHELL = /bin/bash

CFLAGS =
FFLAGS =
!!!FFLAGS="'-C -qopenmp -warn all'"

!!!MY_CFLAGS='-C'

SUPP_OPT = -openmp

OPTIMIZ = 2

MYLIB = libvgrid.a

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

DEFINE = -DFICHSTD98

.PRECIOUS:

# It is assumed that the environment makes available makefile_suffix_rules.inc; e.g.
#     . ssmuse-sh -p /ssm/net/rpn/libs/4.1b/environment-includes_4.3_all


OBJECTS = vgrid_utils.o vgrid.o vgrid_descriptors.o

default: genlib

obj:	$(OBJECTS)

vgrid.o:	vgrid.c \
	BODY_C_compute_heights_4001.hc  BODY_C_compute_pressure_2001.hc  armnlib.h  BODY_C_compute_pressure_5100.hc  \
	BODY_Cvgd_diag_withref_2ref.hc BODY_C_compute_pressure_5002_5003_5004_5005.hc  vgrid.h  \
	BODY_C_compute_heights_0001.hc  BODY_C_compute_pressure_1003_5001.hc  BODY_C_compute_pressure_1001_1002.hc   \
	BODY_C_compute_heights_21001.hc  

vgrid_descriptors.o:	vgrid_descriptors.F90 \
	vgrid_utils.o  BODY_F_dpidpis_withref.hf  BODY_F_dpidpis_withref_prof.hf  BODY_F_diag_withref.hf  \
	BODY_F_levels_withref_prof.hf   BODY_F_levels_withref.hf  BODY_F_diag_withref_prof.hf  vgrid_descriptors.hf  

vgrid_utils.o:	vgrid_utils.F90 vgrid_descriptors.hf  


genlib: $(OBJECTS)
	rm -f $(MYLIB)
	ar rcv $(MYLIB) $(OBJECTS)

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f libvgrid.a
