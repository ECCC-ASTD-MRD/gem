.SUFFIXES :

.SUFFIXES : .F90 .f .o .a

SHELL = /bin/sh

CPP = /lib/cpp

RMNLIB_BASE_PATH=$(ARMNLIB)

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

COMPILE = compile

FFLAGS = 

CFLAGS = 

SUPP_OPT =

OPTIMIZ = 2

CPPFLAGS = -I$(RPN_TEMPLATE_LIBS)/include

LIB_X = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

UPDATEX =

TARRLS = beta

.PRECIOUS: $(MYLIB)

default: burp98_f.o

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

F90DECKS= \
         bufrchr.F90       mrbcol.F90       mrbcov.F90       mrbcvt.F90 \
         mrbdcv.F90       mrbini.F90       mrblocx.F90       mrbprm.F90 \
         mrbprml.F90       mrbrpt.F90       mrbtyp.F90       mrbupd.F90 \
         mrfcls.F90       mrfdel.F90       mrfgoc.F90       mrfgor.F90 \
         mrfloc.F90       mrfmxl.F90       mrfnbr.F90       mrfopc.F90 \
         mrfopn.F90       mrfopr.F90       mrfprm.F90       mrfvoi.F90 \
         qbrptri.F90       qrbnbdt.F90       qrbsct.F90     mrfbfl.F90

COMDECKS= \
         bpl.cdk       burpopt.cdk       codes.cdk       defi.cdk \
         enforc8.cdk   masques.cdk

FICHIERS = $(F90DECKS)

OBJET = burp98_f.o

burp98_f.F90: $(F90DECKS) $(COMDECKS)
#Produire le source du progiciel burp et le regrouper\ndans le fichier burp98_f.f
	cat $(F90DECKS) > burp98_f.F90

tarball: $(F90DECKS) $(COMDECKS)
	tar cfzv /data/armnraid1/www/ssm/sources/burp98_f_$(TARRLS)_all.tgz $(F90DECKS) $(COMDECKS) Makefile 

obj: $(OBJET) 
#Produire les fichiers objets (.o) pour tous les fichiers

burptest.o: burptest.f
#Produire un objet du programme test fortran

burptest: burptest.o $(OBJET) 
#creer un executable burptest
	r.build -o $@ -obj burptest.o -abi $(ABI) -librmn

testc.o: testc.c

ctest: testc.o $(OBJET)
#Produire un executable pour tester l'interface C de burp
	r.build -o $@ -obj testc.o -abi $(ABI) -librmn

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir $(TEMPLIB)(burp98_f.o)
        
genlib: $(OBJET)
#Creer ou mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) $(OBJET)

updlib: 
#mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) *.o ; \
	fi

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	rm -f *.o burp98_f.F90
