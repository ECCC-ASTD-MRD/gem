#!!mak Interpv - makefile for the interp1D package
#!
#!AUTHOR
#!     J.W. Blezius JULY 2002
#!
#!REVISION
#! v1_0    Blezius J.W.          - new code
#! v1_2    Blezius J.W. OCT 2002 - add missing dependencies
#! v1_3    Blezius J.W. OCT 2003 - add Extrap1D_Surface
#!         Blezius J.W. MAY 2014 - add Extrap1D_Abort
#!
#!OBJECT
#!        To compile the components of the 1-D interpolation package, for
#!        inclusion in a library
#!
#!!

.SUFFIXES:

.SUFFIXES : .o .F90 .f90 .inc .a

SHELL = /bin/sh

RMNLIB_BASE_PATH=$(ARMNLIB)

RMNLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn.a

TEMPLIB = ./$(EC_ARCH)/lib_local.a

FFLAGS =

CFLAGS =

SUPP_OPT =

OPTIMIZ = 2
#OPTIMIZ = 0 -debug

CPPFLAGS = -I$(ARMNLIB)/include

MYLIB = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmnbeta.a

LIB_X = $(RMNLIB_BASE_PATH)/lib/$(EC_ARCH)$(ABI)/librmn_x.a

UPDATEX =

TARRLS = beta

.PRECIOUS: $(RMNLIB) $(MALIB)

default: obj

include $(RPN_TEMPLATE_LIBS)/include/makefile_suffix_rules.inc

Interp1D_Arguments_X.inc : Interp1D_Arguments.inc
Interp1D_Constants.o : Interp1D_Constants.F90 Interp1D_Arguments_X.inc
Extrap1D_LapseRate.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Extrap1D_LapseRate.F90 \
	Extrap1D_LapseRate_Body.inc Interp1D_PreX_Shell.inc
Extrap1D_Fixed.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Extrap1D_Fixed.F90 \
	Extrap1D_Fixed_Body.inc Interp1D_PreX_Shell.inc
Extrap1D_Abort.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Extrap1D_Abort.F90 \
	Extrap1D_Abort_Body.inc Interp1D_PreX_Shell.inc
Extrap1D_Surface.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Extrap1D_Surface.F90 \
	Extrap1D_Surface_Body.inc
Extrap1D_SurfaceWind.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Extrap1D_SurfaceWind.F90 \
	Extrap1D_SurfaceWind_Body.inc
Interp1D_CubicLagrange.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Interp1D_CubicLagrange.F90 \
	Interp1D_CubicLagrange_Body.inc Interp1D_PreX_Shell.inc
Interp1D_CubicWithDerivs.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Interp1D_CubicWithDerivs.F90 \
	Interp1D_CubicWithDerivs_Body.inc Interp1D_PreX_Shell.inc
Interp1D_FindPos.o : Interp1D_Constants.o Interp1D_FindPos.F90 Interp1D_FindPos_Body.inc
Interp1D_Linear.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Interp1D_Linear.F90 Interp1D_Linear_Body.inc \
	Interp1D_PreX_Shell.inc
Interp1D_NearestNeighbour.o : Interp1D_Constants.o Interp1D_Arguments_X.inc Interp1D_NearestNeighbour.F90 \
	Interp1D_NearestNeighbour_Body.inc Interp1D_PreX_Shell.inc

OBJECTS= \
        Interp1D_Constants.o Extrap1D_LapseRate.o Extrap1D_Fixed.o Extrap1D_Abort.o Extrap1D_Surface.o \
        Extrap1D_SurfaceWind.o Interp1D_CubicLagrange.o Interp1D_CubicWithDerivs.o \
        Interp1D_FindPos.o Interp1D_Linear.o Interp1D_NearestNeighbour.o

FICHIERS = *.F90 *.f90 *.inc

Test_Interp1D.F90:
	cp test/Test_Interp1D.F90 .

SpeedExtrap.F90:
	cp test/SpeedExtrap.F90 .

SpeedFindPos.F90:
	cp test/SpeedFindPos.F90 .

SpeedInterp.F90:
	cp test/SpeedInterp.F90 .

SpeedFindPos:	SpeedFindPos.o $(OBJECTS)
	s.compile -o SpeedFindPos_$(COMP_ARCH).Abs -obj SpeedFindPos.o $(OBJECTS) \
	-abi$(ABI) -librmn rmn_016.2

SpeedInterp:	SpeedInterp.o $(OBJECTS)
	s.compile -o SpeedInterp_$(COMP_ARCH).Abs -obj SpeedInterp.o $(OBJECTS) \
	-abi$(ABI) -librmn rmn_016.2

SpeedInterp_TestLib:	SpeedInterp.o
	. s.ssmuse.dot /users/dor/arma/jbl/home1/armnlib_sandbox/sandbox_armnlib/sandbox_1/test_interpv; \
	s.compile -O $(OPTIMIZ) $(SUPP_OPT) -o SpeedInterp_TestLib_$(COMP_ARCH).Abs -obj SpeedInterp.o \
	-abi$(ABI) -librmn rmn_015.2 \
	-libpath /users/dor/arma/jbl/home1/armnlib_sandbox/sandbox_armnlib/sandbox_1 \
	-libappl rmn_Alpha_015.3

SpeedExtrap:	SpeedExtrap.o $(OBJECTS)
	s.compile -o SpeedExtrap_$(COMP_ARCH).Abs -obj SpeedExtrap.o $(OBJECTS) \
	-abi$(ABI) -librmn rmn_016.2

Test_Interp1D:	Test_Interp1D.o $(OBJECTS)
	s.compile -o Test_Interp1D_$(COMP_ARCH).Abs -obj Test_Interp1D.o $(OBJECTS) \
	-abi$(ABI) -librmn rmn_016.2

Test_Interp1D_TestLib:	Test_Interp1D.o
	. s.ssmuse.dot /users/dor/arma/jbl/home1/armnlib_sandbox/sandbox_armnlib/sandbox_1/test_interpv; \
	s.compile -O $(OPTIMIZ) $(SUPP_OPT) -o Test_Interp1D_TestLib_$(COMP_ARCH).Abs -obj Test_Interp1D.o \
	-abi$(ABI) -librmn rmn_015.2 \
	-libpath /users/dor/arma/jbl/home1/armnlib_sandbox/sandbox_armnlib/sandbox_1 \
	-libappl rmn_Alpha_015.3

tarball:  $(FICHIERS) Makefile
	tar cfzv /data/armnraid1/www/ssm/sources/interpv_$(TARRLS)_all.tgz $(FICHIERS) Makefile 

obj: $(OBJECTS)
#Produire les fichiers objets (.o) pour tous les fichiers

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 
	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir \
        $(TEMPLIB)(Interp1D_Constants.o) $(TEMPLIB)(Extrap1D_LapseRate.o) \
        $(TEMPLIB)(Extrap1D_Fixed.o) $(TEMPLIB)(Extrap1D_Abort.o) \
        $(TEMPLIB)(Extrap1D_Surface.o) $(TEMPLIB)(Extrap1D_SurfaceWind.o) \
        $(TEMPLIB)(Interp1D_CubicLagrange.o) $(TEMPLIB)(Interp1D_CubicWithDerivs.o) $(TEMPLIB)(Interp1D_FindPos.o) \
        $(TEMPLIB)(Interp1D_Linear.o) $(TEMPLIB)(Interp1D_NearestNeighbour.o)

genlib: $(OBJECTS)
#Creer ou mettre a jour la programmatheque 
	r.ar rcv $(MYLIB) $(OBJECTS)

updlib: 
#mettre a jour la programmatheque 
	r.ar rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar rcv $(LIB_X) *.o ; \
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
	rm -f *.o *.f90 *.mod *.MOD *.stb  *.Abs
	rm -f Test_Interp1D.F90 SpeedExtrap.F90 SpeedFindPos.F90 SpeedInterp.F90
