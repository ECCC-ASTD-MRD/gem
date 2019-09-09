.SUFFIXES :

.SUFFIXES : .c .F .f .o .a .f90

SHELL = /bin/sh

CPP = /lib/cpp

RMNLIB_BASE_PATH=$(ARMNLIB)

RMNLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

FFLAGS =

CFLAGS =

SUPP_OPT =

OPTIMIZ = 2

CPPFLAGS = -I$(RPN_TEMPLATE_LIBS)/include

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

LIB_X = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

UPDATEX =

TARRLS = beta

.PRECIOUS: $(RMNLIB) $(MALIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

CDECKS= \
         burp98.c   fstd98.c   xdf98.c   c_burp.c excdes_new.c \
         memcpy_16_32.c  memcpy_8_16.c if_fstd98.c if_xdf98.c if_burp98.c FC_string.c c_fstgrib_helper.c \
         strgr4a.c r4astrg.c fst_missing.c

FDECKS= \
         convip.f   fstcvt2.f

FTNDECKS= \
          hybref.F

COMDECKS =  proto.h  qstdir.h ip_kind.h

OBJECTS= \
	 burp98.o   fstd98.o   xdf98.o    convip.o \
         fstcvt2.o  c_burp.o   strgr4a.o  r4astrg.o \
         hybref.o   FC_string.o memcpy_8_16.o \
         memcpy_16_32.o c_fstgrib_helper.o fst_missing.o \
         convert_ip123.o convip_plus.o excdes_new.o


FICHIERS = $(FDECKS) $(CDECKS) $(FTNDECKS) $(COMDECKS)

strgr4a.o:	strgr4a.c
r4astrg.o:	r4astrg.c
burp98.o: 	burp98.c if_xdf98.c qstdir.h proto.h
fstd98.o:	fstd98.c if_fstd98.c qstdir.h proto.h
xdf98.o:	xdf98.c if_xdf98.c qstdir.h proto.h
fstdtst2.o:	fstdtst2.f
tstxdf.o:	tstxdf.f
hybref.o:	hybref.F
convip_plus.o:  convip_plus.f90 convert_ip123.o

tarball:  $(FICHIERS) Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/fstd98_$(TARRLS)_all.tgz $(FICHIERS) Makefile 

obj: $(OBJECTS) 
#Produire les fichiers objets (.o) pour tous les fichiers

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir \
         $(TEMPLIB)(burp98.o)   $(TEMPLIB)(fstd98.o)   $(TEMPLIB)(xdf98.o)    $(TEMPLIB)(convip.o) \
         $(TEMPLIB)(fstcvt2.o)  $(TEMPLIB)(c_burp.o)   $(TEMPLIB)(strgr4a.o)  $(TEMPLIB)(r4astrg.o) \
         $(TEMPLIB)(hybref.o)   $(TEMPLIB)(FC_string.o) $(TEMPLIB)(memcpy_8_16.o) \
         $(TEMPLIB)(memcpy_16_32.o) $(TEMPLIB)(c_fstgrib_helper.o) $(TEMPLIB)(fst_missing.o) \
         $(TEMPLIB)(convert_ip123.o) $(TEMPLIB)(convip_plus.o) $(TEMPLIB)(excdes_new.o)

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

debug:
	echo FC=${FC} define=${define} AR=${AR}

fstdtst: $(OBJECTS) fstdtst2.o
	r.build -o $@ -obj fstdtst2.o -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta

test_fstd98: $(OBJECTS) test_fstd98.o
	r.build -o $@ -obj test_fstd98.o $(OBJECTS) -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta -libappl efence
#test_fstd98: $(OBJECTS) test_fstd98.o
#	r.build -o $@ -obj test_fstd98.o -arch $(EC_ARCH) -abi $(ABI) -librmn rmn_005

infol: $(OBJECTS) infol.o
	r.build -o $@ -obj infol.o $(OBJECTS) -arch $(EC_ARCH) -abi $(ABI) -librmn rmn_x -libappl efence

tstxdf: $(OBJECTS) tstxdf.o
	r.build -o $@ -obj tstxdf.o $(OBJECTS) -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	rm -f *.o *.mod
