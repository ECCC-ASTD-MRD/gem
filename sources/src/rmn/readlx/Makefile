.SUFFIXES :

.SUFFIXES : .f .o .a

SHELL = /bin/sh

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

OBJET = readlx.o

OBJET2 = tstrdlx.o 

TARRLS = beta

#CFLAGS = -I$(ARMNLIB)/include/macros -DC910=c910 -g
#FFLAGS = -g

# FFLAGS = -O  (pour cidsv08)

# FFLAGS = -O1 (pour cidsv02)

# FFLAGS = -O1

# DEFINE = -DC910=c910 -DRATFOR=ratfor

FICHIERS = \
 bufftok.cdk   parmadr.cdk   qlxbuff.cdk   qlxfmt.cdk\
 aaaamac.dck   argdims.dck   argdope.dck   lexins.dck    qlxadi.dck\
 qlxadr.dck    qlxasg.dck    qlxbak.dck    qlxcall.dck\
 qlxchr.dck    qlxdbg.dck    qlxdtyp.dck   qlxerr.dck\
 qlxflsh.dck   qlxfnd.dck    qlxind.dck    qlxins.dck\
 qlxnum.dck    qlxnvar.dck   qlxopr.dck    qlxopt.dck\
 qlxpri.dck    qlxprnt.dck   qlxrpn.dck    qlxskp.dck\
 qlxtok.dck    qlxundf.dck   qlxval.dck    qlxxpr.dck\
 readlx.dck

FILES = readlx.f

default: obj

.dck.rat:
	cat *.dck | r.ftnpp - -P $(CPPFLAGS) -I$(RPN_TEMPLATE_LIBS)/include $(DEFINE) > $*.rat

.rat.f:
	rat77 -f $< -l $*.of -o $*.f

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

#readlx.rat: $(FICHIERS)
#readlx.f: readlx.rat
readlx.o: readlx.f

tstrdlx.o: tstrdlx.f

#genrat: $(FICHIERS)
#	cat *.dck | /lib/cpp -P $(CPPFLAGS) -I$(ARMNLIB)/include $(DEFINE) > readlx.rat

tarball:  $(FILES) Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/readlx_$(TARRLS)_all.tgz $(FILES) Makefile 

obj: $(OBJET) 
#Produire les fichiers objets (.o) pour tous les fichiers

tstrdlx: $(OBJET2) $(OBJET)
	r.build -o $@ -obj $(OBJET2) $(OBJET) -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta

tstrdlx12: $(OBJET2) $(OBJET)
	r.build -o $@ -obj $(OBJET2) $(OBJET) /users/dor/armn/lib/obj/node12/*.o -arch $(EC_ARCH) -abi $(ABI) -librmn rmnbeta -mpi

tstc: tstc.o c_readlx.o readln.o
	r.build -o $@ -obj tstc.o c_readlx.o readln.o -arch $(EC_ARCH) -abi $(ABI) -librmn

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir $(TEMPLIB)(readlx.o)        
        
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
	-rm -f *.o *.of tstrdlx

