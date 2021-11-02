.SUFFIXES :

.SUFFIXES : .f90 .F .f .c .o .a .F90

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

FDECKS= \
         dgauss.F       grgg.F       idxes.F       llfxy.F \
         mzonxst.F      npspin.F     ordleg.F      sglo.F \
         sglogg.F       sll.F        sps.F         vglo.F \
         vglogg.F       vll.F        vps.F         xyflcc.F \
         xyfll.F        igaxg95.F    legvar.F      mxma.F \
         grid_to_grid.F

CDECKS= \
         c_time.c         ipsort.c       movr8r4.c  

fDECKS= \
         afix.f           afix8.f        amax.f          amean.f \
         amin.f           bilin.f        corcof.f        coriol.f \
         coupe.f          cvt_base40.f   d1.f            d1int1.f \
         d1int2.f         d1intr.f       datec.f         datmgp.f \
         defvec.f         fconw.f \
         fd1.f            fdm.f          fft_m8.f        gdadcn.f \
         gdadgd.f         gdmpgd.f       gdsqrt.f        \
         grlcc.f          grll.f         grps.f          idatmg.f \
         igaxg.f          incdat.f        int1d1.f \
         inzono.f         isort.f        jdatec.f        jjulien.f \
         movlev.f         movlev8.f \
         mscale.f         mwerms.f       \
         permut.f       polagr.f        rossr3.f \
         s1scor.f         \
         sort.f           spd.f          splat.f \
         sqadsq.f         valide.f \
         vpolagr.f \
         windv.f        xgaig.f \
         system_time.F90  mxm.f\
	 cvmgm.f          cvmgn.f        cvmgp.f         cvmgt.f\
         cvmgz.f          nr_rand.f      fft_m4.f

OBJECTS= \
         afix.o           afix8.o        amax.o          amean.o \
         amin.o           bilin.o        c_time.o        corcof.o \
         coriol.o         coupe.o        cvt_base40.o    d1.o \
         d1int1.o         d1int2.o       d1intr.o        datec.o \
         datmgp.o         defvec.o       dgauss.o        \
         fconw.o          fd1.o          fdm.o           fft_m8.o\
         gdadcn.o         gdadgd.o       gdmpgd.o        gdsqrt.o \
         grgg.o           grlcc.o        grll.o          grps.o \
         idatmg.o         idxes.o        igaxg.o         incdat.o \
         int1d1.o         inzono.o       ipsort.o        isort.o \
         jdatec.o         jjulien.o      llfxy.o         moduledate.o \
         movlev.o         movlev8.o      mscale.o        mwerms.o \
         mzonxst.o        npspin.o       ordleg.o        permut.o \
         polagr.o         rossr3.o       s1scor.o        sglo.o \
         sglogg.o         sll.o          sort.o          spd.o \
         splat.o          sps.o          sqadsq.o        valide.o \
         vglo.o           vglogg.o       vll.o           \
                          vpolagr.o      vps.o           windv.o \
         xgaig.o          xyflcc.o       xyfll.o         system_time.o\
         mxma.o           mxm.o          igaxg95.o       legvar.o\
         movr8r4.o        nr_rand.o      grid_to_grid.o  fft_m4.o\
	 set_user_lock.o

COMDECKS= \
         cubic.cdk       kases.cdk       lagrange.cdk       ordr.cdk \
         pi.cdk       piprm.cdk  qqqcgb1.cdk  qqqgb2.cdk

FICHIERS = $(fDECKS) $(FDECKS) $(CDECKS) $(COMDECKS)

dgauss.o:	dgauss.F 	kases.cdk 	pi.cdk 
grgg.o:		grgg.F 	kases.cdk 	pi.cdk 
idxes.o:	idxes.F 	kases.cdk 
llfxy.o:	llfxy.F 	pi.cdk 
mxma.o:         mxma.F
	s.compile -O 1 -defines "=$(DEFINE)" -optf "=$(FFLAGS)" $(SUPP_OPT) $(FTN_SUPP_OPT) -src $<

mzonxst.o:	mzonxst.F 
npspin.o:       npspin.F 	kases.cdk 	pi.cdk 
ordleg.o:       ordleg.F 	pi.cdk 
sglo.o:		sglo.F 	cubic.cdk 	kases.cdk 	ordr.cdk 
sglogg.o:	sglogg.F 	cubic.cdk 	kases.cdk 	lagrange.cdk \
		ordr.cdk 	pi.cdk 
sll.o:		sll.F 	cubic.cdk 	lagrange.cdk 	ordr.cdk 
sps.o: 		sps.F 	cubic.cdk 	lagrange.cdk 	ordr.cdk 
vglo.o:		vglo.F 	cubic.cdk 	kases.cdk 	pi.cdk 
vglogg.o:	vglogg.F 	cubic.cdk 	kases.cdk 	lagrange.cdk \
		pi.cdk 
vll.o:		vll.F 	pi.cdk 
vps.o:		vps.F 	pi.cdk 
xyflcc.o: 	xyflcc.F 	piprm.cdk 
xyfll.o: 	xyfll.F 	pi.cdk 
legvar.o:       legvar.F
grid_to_grid.o: grid_to_grid.F

tarball: $(FICHIERS)
	tar cfzv /data/armnraid1/www/ssm/sources/base_$(TARRLS)_all.tgz $(FICHIERS) Makefile 

obj: $(OBJECTS) $(OBJ_FRONT)
#Produire les fichiers objets (.o) pour tous les fichiers

genlib: $(OBJECTS) $(OBJ_FRONT)
#Creer ou mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) $(OBJECTS) $(OBJ_FRONT)

gen_ec_arch_dir:
#Creer le repertoire $EC_ARCH 

	mkdir -p ./$(EC_ARCH)

locallib: gen_ec_arch_dir \
         $(TEMPLIB)(afix.o)           $(TEMPLIB)(afix8.o)        $(TEMPLIB)(amax.o)          $(TEMPLIB)(amean.o) \
         $(TEMPLIB)(amin.o)           $(TEMPLIB)(bilin.o)        $(TEMPLIB)(c_time.o)        $(TEMPLIB)(corcof.o) \
         $(TEMPLIB)(coriol.o)         $(TEMPLIB)(coupe.o)        $(TEMPLIB)(cvt_base40.o)    $(TEMPLIB)(d1.o) \
         $(TEMPLIB)(d1int1.o)         $(TEMPLIB)(d1int2.o)       $(TEMPLIB)(d1intr.o)        $(TEMPLIB)(datec.o) \
         $(TEMPLIB)(datmgp.o)         $(TEMPLIB)(defvec.o)       $(TEMPLIB)(dgauss.o)        \
         $(TEMPLIB)(fconw.o)          $(TEMPLIB)(fd1.o)          $(TEMPLIB)(fdm.o)           $(TEMPLIB)(fft_m8.o)\
         $(TEMPLIB)(gdadcn.o)         $(TEMPLIB)(gdadgd.o)       $(TEMPLIB)(gdmpgd.o)        $(TEMPLIB)(gdsqrt.o) \
         $(TEMPLIB)(grgg.o)           $(TEMPLIB)(grlcc.o)        $(TEMPLIB)(grll.o)          $(TEMPLIB)(grps.o) \
         $(TEMPLIB)(idatmg.o)         $(TEMPLIB)(idxes.o)        $(TEMPLIB)(igaxg.o)         $(TEMPLIB)(incdat.o) \
         $(TEMPLIB)(int1d1.o)         $(TEMPLIB)(inzono.o)       $(TEMPLIB)(ipsort.o)        $(TEMPLIB)(isort.o) \
         $(TEMPLIB)(jdatec.o)         $(TEMPLIB)(jjulien.o)      $(TEMPLIB)(llfxy.o)         $(TEMPLIB)(moduledate.o) \
         $(TEMPLIB)(movlev.o)         $(TEMPLIB)(movlev8.o)      $(TEMPLIB)(mscale.o)        $(TEMPLIB)(mwerms.o) \
         $(TEMPLIB)(mzonxst.o)        $(TEMPLIB)(npspin.o)       $(TEMPLIB)(ordleg.o)        $(TEMPLIB)(permut.o) \
         $(TEMPLIB)(polagr.o)         $(TEMPLIB)(rossr3.o)       $(TEMPLIB)(s1scor.o)        $(TEMPLIB)(sglo.o) \
         $(TEMPLIB)(sglogg.o)         $(TEMPLIB)(sll.o)          $(TEMPLIB)(sort.o)          $(TEMPLIB)(spd.o) \
         $(TEMPLIB)(splat.o)          $(TEMPLIB)(sps.o)          $(TEMPLIB)(sqadsq.o)        $(TEMPLIB)(valide.o) \
         $(TEMPLIB)(vglo.o)           $(TEMPLIB)(vglogg.o)       $(TEMPLIB)(vll.o)           \
         $(TEMPLIB)(vpolagr.o)        $(TEMPLIB)(vps.o)          $(TEMPLIB)(windv.o) \
         $(TEMPLIB)(xgaig.o)          $(TEMPLIB)(xyflcc.o)       $(TEMPLIB)(xyfll.o)         $(TEMPLIB)(system_time.o)\
         $(TEMPLIB)(mxma.o)           $(TEMPLIB)(mxm.o)          $(TEMPLIB)(igaxg95.o)       $(TEMPLIB)(legvar.o)\
         $(TEMPLIB)(movr8r4.o)        $(TEMPLIB)(nr_rand.o)      $(TEMPLIB)(grid_to_grid.o)  $(TEMPLIB)(fft_m4.o)\
         $(TEMPLIB)(set_user_lock.o)
       
updlib: 
#mettre a jour la programmatheque 
	r.ar -arch $(EC_ARCH) rcv $(MYLIB) *.o
	if [ "$(UPDATEX)" = "1" ] ; \
	then \
	r.ar -arch $(EC_ARCH) rcv $(LIB_X) *.o ; \
	fi

clean:
#Faire le grand menage. On enleve tous les fichiers sources\ninutiles et les .o 
	rm -f *.o
