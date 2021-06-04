.SUFFIXES :

.SUFFIXES : .f90 .ftn .f .c .o .a

SHELL = /bin/sh

CPP = /lib/cpp

FFLAGS =

CFLAGS = 

SUPP_OPT =

OPTIMIZ = 2

CPPFLAGS = -I$(RPN_TEMPLATE_LIBS)/include

TEMPLIB = ./$(EC_ARCH)/lib_local.a

MYLIB = $(ARMNLIB)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

LIB_X = $(ARMNLIB)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

UPDATEX =

TARRLS = beta

.PRECIOUS: $(MYLIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

OBJECTS= armn_compress_32.o zfstlib.o c_zfstlib.o

tarball:  armn_compress_32.c c_zfstlib.c  zfstlib.f90  zfstlib.h  Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/compresseur_$(TARRLS)_all.tgz armn_compress_32.c c_zfstlib.c  zfstlib.f90  zfstlib.h Makefile 

obj: $(OBJECTS) 
#Produire les fichiers objets (.o) pour tous les fichiers

obj2: $(OBJECTS) 
#Produire les fichiers objets (.o) pour tous les fichiers

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir $(TEMPLIB)(armn_compress_32.o) $(TEMPLIB)(zfstlib.o) $(TEMPLIB)(c_zfstlib.o)

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

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	-if [ "*.ftn" != "`echo *.ftn`" ] ; \
	then \
	for i in *.ftn ; \
	do \
	fn=`r.basename $$i '.ftn'`; \
	rm -f $$fn.f; \
	done \
	fi
	rm -f *.o
