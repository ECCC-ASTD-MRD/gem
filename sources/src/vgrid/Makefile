SHELL = /bin/bash

include vgrid_version.inc


.SUFFIXES :

.SUFFIXES : .F90 .c .o .a

SHELL = /bin/bash

CPP = /lib/cpp

FFLAGS =

CFLAGS =

DEFINE = -DFICHSTD98

SUPP_OPT =

OPTIMIZ = 2 

#OPTIMIZ = -O 3 -fast
#C_OPTIMIZ = -O 2
#F_OPTIMIZ = -O 2

MYLIB = libvgrid_$(VGRID_VERSION).a

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

OBJET = vgrid_utils.o vgrid.o vgrid_descriptors.o

FICHIERS_F90 = \
vgrid_utils.F90 vgrid_descriptors.F90

FICHIERS_C = \
vgrid.c

FICHIERS_H = \
armnlib.h   vgrid.h

FICHIERS = $(FICHIERS_F90) $(FICHIERS_C)

obj: $(OBJET)
#Produire les fichiers objets (.o) pour tous les fichiers

genlib: $(OBJET)
	rm -f $(MYLIB)
	ar rv $(MYLIB) $(OBJET)

clean:
	rm -f *.o *.mod


