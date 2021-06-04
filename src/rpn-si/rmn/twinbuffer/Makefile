.SUFFIXES :

.SUFFIXES : .ftn .f .c .o .a

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

.PRECIOUS: (MYLIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

OBJECT = I_O_2Buffer.o

OBJECT2 = test.o

FICHIERS = I_O_2Buffer.c 

tarball:  $(FICHIERS) Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/twinbuffer_$(TARRLS)_all.tgz $(FICHIERS) Makefile 

obj: $(OBJECT) 
#Produire les fichiers objets (.o) pour tous les fichiers

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir $(TEMPLIB)(I_O_2Buffer.o)

genlib: $(OBJECT)
#Creer ou mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) $(OBJECT)

updlib: 
#mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) *.o ; \
	fi

test: $(OBJECT2)
	r.build -o $@ -obj $(OBJECT2) $(OBJECT) -arch $(EC_ARCH) -abi $(ABI) -librmn

clean:
#Faire le grand menage. On enleve tous les fichiers sources inutiles et les .o 
	-if [ "*.ftn" != "`echo *.ftn`" ] ; \
	then \
	for i in *.ftn ; \
	do \
	fn=`r.basename $$i '.ftn'`; \
	rm -f $$fn.f; \
	done \
	fi
	rm -f *.o
