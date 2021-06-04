.SUFFIXES :

.SUFFIXES : .f90 .F .c .o .a

SHELL = /bin/sh

CPP = /lib/cpp

RMNLIB_BASE_PATH=$(ARMNLIB)

RMNLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

FFLAGS =

CFLAGS =

SUPP_OPT =

OPTIMIZ = 2

# SUPP_OPT = -mpi -openmp

CPPFLAGS = -I$(RPN_TEMPLATE_LIBS)/include

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

LIB_X = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

UPDATEX =

TARRLS = beta

.PRECIOUS: $(RMNLIB) $(MALIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

FDECKS= \
         bmf_statfld.F

F90DECKS= \
                              bmf_catalog.f90     bmf_gobe.f90       bmf_splitend.f90 \
         bmf_splitname.f90    bmf_splitwrite.f90  bmf_splitwritey.f90 \
                              bmf_clear.f90       bmf_get.f90         bmf_init.f90       bmf_splithalo.f90 \
         bmf_splitoptions.f90 bmf_splitwriteh.f90 bmf_write.f90                       0  bmf_connect.f90 \
         bmf_get2.f90         bmf_list.f90        bmf_splithole.f90   bmf_splitstart.f90 bmf_splitwritex.f90 \
         bmf_write2.f90       bmf_errorlvl.f90    bmf_geth.f90        bmf_perturb.f90 \
         bmf_splitinit.f90    bmf_splitwrall.f90  bmf_splitwritexy.f90 \
         bmf_splitwrite2.f90  bmf_splitwriteh2.f90  bmf_splitwritex2.f90  bmf_splitwritexy2.f90 \
         bmf_splitwritey2.f90 bmf_mod.f90         bmf_modsplit.f90

CDECKS= \
         bmf_char2i.c   

DOTMOD= \
         bmf_mod.o           bmf_modsplit.o

OBJECTS= \
         bmf_statfld.o \
                            bmf_catalog.o     bmf_gobe.o       bmf_splitend.o \
         bmf_splitname.o    bmf_splitwrite.o  bmf_splitwritey.o \
                            bmf_clear.o       bmf_get.o         bmf_init.o       bmf_splithalo.o \
         bmf_splitoptions.o bmf_splitwriteh.o bmf_write.o                        bmf_connect.o \
         bmf_get2.o         bmf_list.o        bmf_splithole.o   bmf_splitstart.o bmf_splitwritex.o \
         bmf_write2.o       bmf_errorlvl.o    bmf_geth.o        bmf_perturb.o \
         bmf_splitinit.o    bmf_splitwrall.o  bmf_splitwritexy.o \
         bmf_char2i.o \
         bmf_mod.o          bmf_modsplit.o \
         bmf_splitwrite2.o  bmf_splitwriteh2.o  bmf_splitwritex2.o  bmf_splitwritexy2.o \
         bmf_splitwritey2.o

FICHIERS = $(FDECKS) $(F90DECKS) $(CDECKS)
tarball: $(FICHIERS)
	tar cfzv /data/armnraid1/www/ssm/sources/bmf_$(TARRLS)_all.tgz $(FICHIERS) Makefile 

obj: $(DOTMOD) $(OBJECTS)
#Produire les fichiers objets (.o)

genlib: $(OBJECTS)
#Creer ou mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) $(OBJECTS)

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir $(DOTMOD) \
         $(TEMPLIB)(bmf_statfld.o) \
                                        $(TEMPLIB)(bmf_catalog.o) \
         $(TEMPLIB)(bmf_gobe.o)         $(TEMPLIB)(bmf_splitend.o) \
         $(TEMPLIB)(bmf_splitname.o)    $(TEMPLIB)(bmf_splitwrite.o)  $(TEMPLIB)(bmf_splitwritey.o) \
                                        $(TEMPLIB)(bmf_clear.o)       $(TEMPLIB)(bmf_get.o) \
         $(TEMPLIB)(bmf_init.o)         $(TEMPLIB)(bmf_splithalo.o) \
         $(TEMPLIB)(bmf_splitoptions.o) $(TEMPLIB)(bmf_splitwriteh.o) $(TEMPLIB)(bmf_write.o) \
                                        $(TEMPLIB)(bmf_connect.o) \
         $(TEMPLIB)(bmf_get2.o)         $(TEMPLIB)(bmf_list.o)        $(TEMPLIB)(bmf_splithole.o) \
         $(TEMPLIB)(bmf_splitstart.o)   $(TEMPLIB)(bmf_splitwritex.o) \
         $(TEMPLIB)(bmf_write2.o)       $(TEMPLIB)(bmf_errorlvl.o) \
         $(TEMPLIB)(bmf_geth.o)         $(TEMPLIB)(bmf_perturb.o) \
         $(TEMPLIB)(bmf_splitinit.o)    $(TEMPLIB)(bmf_splitwrall.o)  $(TEMPLIB)(bmf_splitwritexy.o) \
         $(TEMPLIB)(bmf_char2i.o) \
         $(TEMPLIB)(bmf_mod.o) \
         $(TEMPLIB)(bmf_modsplit.o) \
         $(TEMPLIB)(bmf_splitwrite2.o)  $(TEMPLIB)(bmf_splitwriteh2.o) $(TEMPLIB)(bmf_splitwritex2.o) \
         $(TEMPLIB)(bmf_splitwritexy2.o)                              $(TEMPLIB)(bmf_splitwritey2.o)
       
updlib: 
#mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) *.o ; \
	fi

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	rm -f *.o *.mod
