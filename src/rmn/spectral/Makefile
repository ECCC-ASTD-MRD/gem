.SUFFIXES :

.SUFFIXES : .F .c .o .a

SHELL = /bin/sh

CPP = /lib/cpp

FFLAGS =

CFLAGS =

SUPP_OPT =

OPTIMIZ = 2

RMNLIB_BASE_PATH=$(ARMNLIB)

CPPFLAGS = -I$(RPN_TEMPLATE_LIBS)/include

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

LIB_X = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

UPDATEX =

TARRLS = beta

.PRECIOUS: $(MYLIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

FDECKS= \
         alpas2.F       alpas8.F       alpdl2.F       alpdl8.F \
         alpdr2.F       alpdr8.F       alpnm2.F       alpnm8.F \
         dimcal.F       epsil2.F       epsil8.F       fasp8.F \
         gauss.F        gauss8.F       ggasp0.F \
         ggasp8.F       gwaqd0.F       gwaqd8.F       ordleg8.F \
         perm.F         qqqfasp.F      scof2.F        spagg0.F \
         spagg8.F       fasp2.F        ewfdr2.F       spaf2.F \
         spaf8.F

OBJECTS= \
         alpas2.o       alpas8.o       alpdl2.o       alpdl8.o \
         alpdr2.o       alpdr8.o       alpnm2.o       alpnm8.o \
         dimcal.o       epsil2.o       epsil8.o       fasp8.o \
         gauss.o        gauss8.o       ggasp0.o \
         ggasp8.o       gwaqd0.o       gwaqd8.o       ordleg8.o \
         perm.o         qqqfasp.o      scof2.o        spagg0.o \
         spagg8.o       fasp2.o        ewfdr2.o       spaf2.o \
         spaf8.o


FICHIERS = $(FDECKS)

tarball:  $(FICHIERS) Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/spectral_$(TARRLS)_all.tgz $(FICHIERS) Makefile 

obj: $(OBJECTS) 
#Produire les fichiers objets (.o) pour tous les fichiers

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir \
         $(TEMPLIB)(alpas2.o)       $(TEMPLIB)(alpas8.o)       $(TEMPLIB)(alpdl2.o)       $(TEMPLIB)(alpdl8.o) \
         $(TEMPLIB)(alpdr2.o)       $(TEMPLIB)(alpdr8.o)       $(TEMPLIB)(alpnm2.o)       $(TEMPLIB)(alpnm8.o) \
         $(TEMPLIB)(dimcal.o)       $(TEMPLIB)(epsil2.o)       $(TEMPLIB)(epsil8.o)       $(TEMPLIB)(fasp8.o) \
         $(TEMPLIB)(gauss.o)        $(TEMPLIB)(gauss8.o)       $(TEMPLIB)(ggasp0.o) \
         $(TEMPLIB)(ggasp8.o)       $(TEMPLIB)(gwaqd0.o)       $(TEMPLIB)(gwaqd8.o)       $(TEMPLIB)(ordleg8.o) \
         $(TEMPLIB)(perm.o)         $(TEMPLIB)(qqqfasp.o)      $(TEMPLIB)(scof2.o)        $(TEMPLIB)(spagg0.o) \
         $(TEMPLIB)(spagg8.o)       $(TEMPLIB)(fasp2.o)        $(TEMPLIB)(ewfdr2.o)       $(TEMPLIB)(spaf2.o) \
         $(TEMPLIB)(spaf8.o)

genlib: $(OBJECTS)
#Creer ou mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) $(OBJECTS)

updlib: 
#mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) *.o ; \
	fi

copie:
	for i in $(FDECKS) ;\
	do \
	cp /users/dor/armn/lib/rmnlib/spmo/$$i $$i;\
	done

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	rm -f *.o
