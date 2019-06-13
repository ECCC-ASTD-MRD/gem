#ifdef __DOC__
/*
!   definition des fonctions thermodynamiques de base
!   pour les constantes, utiliser le common /consphy/
!     note: toutes les fonctions travaillent avec les unites s.i.
!           i.e. ttt en deg k, prs en pa, qqq en kg/kg
!          *** n. brunet - mai 90 ***
!          * revision 01 - mai 94 - n. brunet
!                          nouvelle version pour faibles pressions
!          * revision 02 - aout 2000 - j-p toviessi
!                          calcul en real(REAL64) ::
!          * revision 03 - sept 2000 - n. brunet
!                          ajout de nouvelles fonctions
!          * revision 04 - janv 2000 - j. mailhot
!                          fonctions en phase mixte
!          * revision 05 - dec 2001 - g. lemay
!                          double precision pour phase mixte
!          * revision 06 - avr 2002 - a. plante
!                          ajout des nouvelles fonctions fottvh et fotvht
!
!* FOEW(ttt):
!    calcule la tension de vapeur saturante (formule de Tetens);
!    on obtient ew (tension de vapeur par rapport à l’eau) ou ei
!    (tension de vapeur par rapport à la glace) selon la température.
!* FODLE(ttt):
!    calcule la dérivée selon la température du ln ew ou ln ei
!* FOQST(ttt, prs):
!    calcule l’humidité spécifique saturante (qst) à partir de la
!    température et de la pression.
!* FODQS(qst, ttt):
!    calcule la dérivée de qst selon la température. qst est la sortie de FOQST
!* FOEFQ(qqq, prs):
!    calcule la tension de vapeur à partir de l’humidité spécifique et
!    de la pression.
!* FOQFE(eee, prs):
!    calcule l’humidité spécifique à partir de la tension de vapeur et
!    de la pression.
!* FOTVT(ttt, qqq):
!    calcule la température virtuelle (tvi) à partir de la température et
!    de l’humidité spécifique.
!* FOTVHT(ttt,qqq,qqh):
!    calcule la temp virt. (tvi) de temp (ttt), hum sp (qqq) et masse sp
!    des hydrometeores.
!* FOTTV(tvi, qqq):
!    calcule la température à partir de la température virtuelle et
!    de l’humidité spécifique.
!* FOTTVH(tvi,qqq,qqh):
!    calcule ttt de temp virt. (tvi), hum sp (qqq) et masse sp des
!    hydrometeores (qqh)
!* FOHR(qqq, ttt, prs):
!    calcule l’humidité relative à partir de l’humidité spécifique,
!    de la température et de la pression. Le résultat est en fraction
!    (et non en pourcentage).
!* FOLV(ttt):
!    calcule la chaleur latente de condensation, fonction de la
!    température (J/kg).
!* FOLS(ttt):
!    calcule la chaleur latente de sublimation, fonction de la température (J/kg)
!* FOPOIT(t0, p0, pf):
!    résout l’équation de Poisson pour la température; 
!    si pf=100000 pa, on obtient le theta standard (K).
!* FOPOIP(t0, tf, p0):
!    résout l’équation de Poisson pour la pression (pa).
!
!Les cinq fonctions internes suivantes sont valides dans le contexte où
!on ne désire pas tenir compte de la phase glace dans les calculs de saturation,
!c’est-à-dire que la phase eau seulement est considérée quelque soit la température.
!
!* FOEWA(ttt):
!    tension de vapeur saturante, ew.
!* FODLA(ttt):
!    dérivée du ln ew selon la température.
!* FOQSA(ttt, prs):
!    humidité spécifique saturante. 11 janvier, 2001 3
!* FODQA(qst, prs):
!    la dérivée de FOQSA selon la température.
!* FOHRA(qqq, ttt, prs):
!    l’humidité relative (résultat en fraction et non en pourcentage)
!
!Definition of basic thermodynamic functions in mixed-phase mode
!fff is the fraction of ice and ddff its derivative w/r to t
!* FESI(ttt):
!    saturation calculations in presence of liquid phase only function for saturation vapor pressure (tetens)
!* ...
!
!Dans toutes les fonctions internes ci-dessus, on a le symbolisme suivant:
!* eee : la tension de vapeur en pa
!* prs ou p0 (p zéro) ou pf : la pression en pa.
!* qqq : l’humidité spécifique en kg/kg
!* qst : l’humidité spécifique saturante en kg/kg
!* tvi : la température virtuelle en deg K
!* ttt ou t0 (t zér0) ou tf: la température en deg K
!
!Reference:
!https://wiki.cmc.ec.gc.ca/images/9/9e/RPNPhy_Thermodynamic_functions_brunet.pdf
!https://wiki.cmc.ec.gc.ca/images/4/4c/RPNPhy_Thermodynamic_functions.pdf
*/
#endif

#ifndef __TDPACK_FUNC__
#define __TDPACK_FUNC__

#ifdef __FORTRAN__

#define _DSIGN(x,y)  dsign(x,y)
#define _DABS(x)     dabs(x)
#define _DMIN1(x,y)  dmin1(x,y)
#define _DMAX1(x,y)  dmax1(x,y)
#define _DLOG(x)     dlog(x)
#define _DEXP(x)     dexp(x)
#define _POW(x,y)    ((x)**(y))

#else

#define _DSIGN(x,y)  ((y >= 0) ? abs(x) : (-abs(x)))
#define _DABS(x)     abs(x)
#define _DMIN1(x,y)  (((x) < (y)) ? (x) : (y))
#define _DMAX1(x,y)  (((x) > (y)) ? (x) : (y))
#define _DLOG(x)     log(x)
#define _DEXP(x)     exp(x)
#define _POW(x,y)    pow(x,y)

#endif

#ifndef _DBLE
#ifdef __FORTRAN__
#define _DBLE(x)     dble(x)
#else
#define _DBLE(x)     ((double)x)
#endif
#endif

#define _ZERO  _DBLE(0.)
#define _ONE   _DBLE(1.)

#define _CAPPA _DBLE(CAPPA)
#define _CHLC  _DBLE(CHLC)
#define _CHLF  _DBLE(CHLF)
#define _CPV   _DBLE(CPV)
#define _DELTA _DBLE(DELTA)
#define _EPS1  _DBLE(EPS1)
#define _EPS2  _DBLE(EPS2)
#define _TRPL  _DBLE(TRPL)

#define _CTE1 (AERK2W*(_TRPL-AERK3W))
#define _CTE2 (AERK2I*(_TRPL-AERK3I)-_CTE1)
#define _CTE3 _DBLE(2317.)
#define _CTE4 _DBLE(7.24)
#define _CTE5 _DBLE(128.4)

#define _DIFTRPL(ttt)          (_DBLE(ttt)-_TRPL)

#define _FN1(ttt)              (_DBLE(ttt)-AERK3W+_DMAX1(_ZERO,_DSIGN(AERK3W-AERK3I,-_DIFTRPL(ttt))))

#define MASKT(ttt)             _DMAX1(_ZERO,_DSIGN(_ONE,_DIFTRPL(ttt)))

#define FOMULTS(ddd,ttt)       ((AERK1W*MASKT(ttt)+(_ONE-MASKT(ttt))*AERK1I)*ddd)

#define FOEWF(ttt)             (_DMIN1(_DSIGN(AERK2W,_DIFTRPL(ttt)),_DSIGN(AERK2I,_DIFTRPL(ttt)))*_DABS(_DIFTRPL(ttt))/_FN1(ttt))
#define FOEW(ttt)              (FOMULTS(_DEXP(FOEWF(ttt)),ttt))
#define FOEWAF(ttt)            (AERK2W*(_DIFTRPL(ttt))/(_DBLE(ttt)-AERK3W))
#define FOEWA(ttt)             (AERK1W*_DEXP(FOEWAF(ttt)))

#define FODLE(ttt)             ((_CTE1+_DMAX1(_ZERO,_DSIGN(_CTE2,-_DIFTRPL(ttt))))/_POW(_FN1(ttt),2.))

#define FOQST(ttt,prs)         (_EPS1/(_DMAX1(_ONE,_DBLE(prs)/FOEW(ttt))-_EPS2))
#define FOQSTX(prs,ddd)        (_EPS1/(_DMAX1(_ONE,_DBLE(prs)/ddd)-_EPS2))
#define FOQSA(ttt,prs)         (_EPS1/(_DMAX1(_ONE,_DBLE(prs)/FOEWA(ttt))-_EPS2))
#define FQSMX(ttt,prs,fff)     (_EPS1/(_DMAX1(_ONE,_DBLE(prs)/FESMX(ttt,fff))-_EPS2))
#define FQSMXX(fesmx8,prs)     (_EPS1/(_DMAX1(_ONE,_DBLE(prs)/fesmx8)-_EPS2))

#define FODQS(qst,ttt)         (_DBLE(qst)*(_ONE+_DELTA*_DBLE(qst))*FODLE(ttt))

#define FOEFQ(qqq,prs)         (_DMIN1(_DBLE(prs),(_DBLE(qqq)*_DBLE(prs))/(_EPS1+_EPS2*_DBLE(qqq))))

#define FOQFE(eee,prs)         (_DMIN1(_ONE,_EPS1*_DBLE(eee)/(_DBLE(prs)-_EPS2*_DBLE(eee))))

#define FOTVT(ttt,qqq)         (_DBLE(ttt)*(_ONE+_DELTA*_DBLE(qqq)))
#define FOTTV(tvi,qqq)         (_DBLE(tvi)/(_ONE+_DELTA*_DBLE(qqq)))

#define FOTVHT(ttt,qqq,qqh)    (_DBLE(ttt)*(_ONE+_DELTA*_DBLE(qqq)-_DBLE(qqh)))
#define FOTTVH(tvi,qqq,qqh)    (_DBLE(tvi)/(_ONE+_DELTA*_DBLE(qqq)-_DBLE(qqh)))

#define FODQA(qst,ttt)         (_DBLE(qst)*(_ONE+_DELTA*_DBLE(qst))*FODLA(ttt))
#define FDQSMX(qsm,dlemx)      (_DBLE(qsm)*(_ONE+_DELTA*_DBLE(qsm))*_DBLE(dlemx))

#define FOHR(qqq,ttt,prs)      (_DMIN1(_DBLE(prs),FOEFQ(qqq,prs))/FOEW(ttt))
#define FOHRX(qqq,prs,ddd)     (_DMIN1(_DBLE(prs),FOEFQ(qqq,prs))/ddd)

#define FOLV(ttt)              (_CHLC-(_CTE3*_DIFTRPL(ttt)))

#define FOLS(ttt)              (_CHLC+_CHLF+(_CPV-(_CTE4*_DBLE(ttt)+_CTE5))*_DIFTRPL(ttt))

#define FOPOIT(t00,pr0,pf)     (_DBLE(t00)*_POW(_DBLE(pr0)/_DBLE(pf),-_CAPPA))

#define FOPOIP(t00,tf,pr0)     (_DBLE(pr0)*_DEXP(-(_DLOG(_DBLE(t00)/_DBLE(tf))/_CAPPA)))

#define FODLA(ttt)             (_CTE1/_POW(_DBLE(ttt)-AERK3W,2.))


#define FOHRA(qqq,ttt,prs)     (_DMIN1(_DBLE(prs),FOEFQ(qqq,prs))/FOEWA(ttt))

#define FESIF(ttt)             (AERK2I*_DIFTRPL(ttt)/(_DBLE(ttt)-AERK3I))
#define FESI(ttt)              (AERK1I*_DEXP(FESIF(ttt)))

#define FDLESI(ttt)            (AERK2I*(_TRPL-AERK3I)/_POW(_DBLE(ttt)-AERK3I,2.))

#define FESMX(ttt,fff)             ((_ONE-_DBLE(fff))*FOEWA(ttt)+_DBLE(fff)*FESI(ttt))

#define FESMXX(fff,fesi8,foewa8)   ((_ONE-_DBLE(fff))*foewa8+_DBLE(fff)*fesi8)

#define FDLESMX(ttt,fff,ddff)      (((_ONE-_DBLE(fff))*FOEWA(ttt)*FODLA(ttt)+_DBLE(fff)*FESI(ttt)*FDLESI(ttt)+_DBLE(ddff)*(FESI(ttt)-FOEWA(ttt)))/FESMX(ttt,fff))
#define FDLESMXX(ttt,fff,ddff,foewa8,fesi8,fesmx8)   (((_ONE-_DBLE(fff))*foewa8*FODLA(ttt)+_DBLE(fff)*fesi8*FDLESI(ttt)+_DBLE(ddff)*(fesi8-foewa8))/fesmx8)

#endif
