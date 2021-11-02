SHELL = /bin/sh

CPP = /lib/cpp

RMNLIB_BASE_PATH=$(ARMNLIB)

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

LIB_X = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

UPDATEX =

COMPILE = compile

FFLAGS =

CFLAGS =

SUPP_OPT =

OPTIMIZ = 2
#OPTIMIZ = -debug

CPPFLAGS = -I$(RPN_TEMPLATE_LIBS)/include

TARRLS = beta

.PRECIOUS: $(MYLIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc


OBJECTS= \
         pak.o compact.o   compact_integer.o     compact_rle.o\
         unpackWrapper.o   compact_IEEEblock.o   f_pow.o\
	       float_packer.o


OBJ_Test= \
          cmainInteger.o   cmainFloat.o          cmainDouble.o\
          cmainRLE.o       cmainIEEE.o           cmainIEEEDouble.o\
          tpak.o           mins.o

OBJ_Test2= \
          cmainInteger.o   cmainFloat.o          cmainDouble.o\
          cmainRLE.o       cmainIEEE.o           cmainIEEEDouble.o\
          test_I.o           mins.o

FICHIERS = pak.c compact.c compact_h.c compact_integer.c compact_rle.c \
           unpackWrapper.c compact_IEEEblock.c compact_IEEEblock_h.c \
           f_pow.f float_packer.c
          
tarball:  $(FICHIERS) Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/packers_$(TARRLS)_all.tgz $(FICHIERS) Makefile 
          
obj: $(OBJECTS) $(OBJ_NEC)
#Produire les fichiers objets (.o) pour tous les fichiers

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir \
         $(TEMPLIB)(pak.o) $(TEMPLIB)(compact.o)   $(TEMPLIB)(compact_integer.o)     $(TEMPLIB)(compact_rle.o)\
         $(TEMPLIB)(unpackWrapper.o)   $(TEMPLIB)(compact_IEEEblock.o)   $(TEMPLIB)(f_pow.o)\
         $(TEMPLIB)(float_packer.o)

genlib: $(OBJECTS) $(OBJ_NEC)
#Creer ou mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) $(OBJECTS) $(OBJ_NEC)

updlib: 
#mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) *.o ; \
	fi

mytest: $(OBJECTS) $(OBJ_Test) $(LIB)
	r.build -o $@ -obj $(OBJECTS) $(OBJ_Test) -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta -debug

libtest: $(OBJ_Test)$(LIB)
	r.build -o $@ -obj $(OBJ_Test) -arch $(EC_ARCH) -abi $(ABI) -librmn rmn_rc008 -debug

test_I: $(OBJECTS) $(OBJ_Test2) $(LIB)
	r.build -o $@ -obj $(OBJECTS) $(OBJ_Test2) -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta -debug

clean:
	rm -f *.o mytest
