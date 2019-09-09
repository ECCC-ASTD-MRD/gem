.SUFFIXES :

.SUFFIXES : .ftn .f .c .o

SHELL = /bin/sh

CPP = /lib/cpp

##MASSVLIB = $(ARMNLIB)/lib/$(ARCH)$(ABI)/libmassvp4_safe.a

FFLAGS =

CFLAGS =

OPTIMIZ = 2

CPPFLAGS = 

.PRECIOUS: $(MASSVLIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

OBJECTS= \
vasin.o       vdiv.o        vscos.o       vsincos.o     vsqrt.o       vssqrt.o \
vatan2.o      vdnint.o      vscosisin.o   vslog.o       vsrec.o       vstan.o \
vcos.o        vexp.o        vrec.o        vsdiv.o       vsrsqrt.o     vtan.o \
vcosisin.o    vlog.o        vrsqrt.o      vsexp.o       vssin.o       vpow_non_ibm.o \
vdint.o       vsatan2.o     vsin.o        vssincos.o

tarball:
	tar cfzv /data/armnraid1/www/ssm/sources/libmassvp4_1.0_all.tgz v*.f Makefile

obj: $(OBJECTS)
#Produire les fichiers objets (.o) pour tous les fichiers

##genlib: $(OBJECTS)
###Creer ou mettre a jour la programmatheque 
##	r.ar rcv $(MASSVLIB) $(OBJECTS)

locallib: obj
	mkdir -p $(EC_ARCH)
	ar rcv $(EC_ARCH)/lib_local.a $(OBJECTS)

##updlib: 
###mettre a jour la programmatheque 
##	r.ar rcv $(MASSVLIB) *.o

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	rm -f *.o
