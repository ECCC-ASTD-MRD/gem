
#ifdef __DOC__
/*
!# fpp tdpack_const.h provides the pre-processing result
!# fpp -D__FORTRAN__ tdpack_const.h provides the pre-processing result

!  AI     : pour fn htvocp
!  AW     : pour fn htvocp
!  BI     : pour fn htvocp
!  BW     : pour fn htvocp
!  CAPPA  : rgasd/cpd
!  CHLF   : ch. lat. fusion       [J kg-1]
!  CHLC   : ch. lat. condens.(0C) [J kg-1]
!  CONSOL : constante solaire     [W m-2]
!  CONSOL2: constante solaire (corrigee)    [W m-2]
!  CPD    : chal. spec. air sec (pression constante) [J kg-1 K-1]
!  CPV    : chal. spec. vap eau   [J kg-1 K-1]
!  CPI    : chal. spec. glace     [J kg-1 K-1]
!  CPW    : chal. spec. eau liq.  [J kg-1 K-1]
!  CVD    : chal. spec. air sec (volume constant) [J kg-1 K-1]
!  DELTA  : 1/eps1 - 1
!  EPS1   : rgasd/rgasv
!  EPS2   : 1 - eps1
!  GRAV   : acc. de gravite       [m s-2]
!  KARMAN : cte de von karman
!  KNAMS  : passage kt a m/s
!  OMEGA  : rotation terre        [s-1]
!  PI     : cte pi:acos(-1)
!  RAUW   : densite eau liq       [kg m-3]
!  RAYT   : rayon moy. terre      [m]
!  RGASD  : cte gaz - air sec     [J kg-1 K-1]
!  RGASV  : cte gaz - vap eau     [J kg-1 K-1]
!  RIC    : cte richardson crit.
!  SLP    : pour fn htvocp
!  STEFAN : cte stefan-boltzmann  [J m-2 s-1 K-4]
!  STLO   : schuman-newell l.r.   [K s2 m-2]
!  T1S    : pour fn htvocp        [K]
!  T2S    : pour fn htvocp        [K]
!  TCDK   : passage k a c         [C]
!  TGL    : temp glace dans atm   [K]
!  TRPL   : point triple - eau    [K]
!
!  !# Tetens coefficients in saturation vapor pressure formulas
!  TTNS1  :
!  TTNS3W :
!  TTNS3I :
!  TTNS4W :
!  TTNS4I :
!
!  !# Alduchov and Eskridge (1995) coefficients in saturation vapor pressure formulas
!  AERK1W : voir Alduchov and Eskridge (1995) -AERK
!  AERK2W : voir Alduchov and Eskridge (1995) -AERK
!  AERK3W : voir Alduchov and Eskridge (1995) -AERK
!  AERK1I : voir Alduchov and Eskridge (1995) -AERKi
!  AERK2I : voir Alduchov and Eskridge (1995) -AERKi
!  AERK3I : voir Alduchov and Eskridge (1995) -AERKi
 */
#endif

#ifndef __TDPACK_CONST__
#define __TDPACK_CONST__

#define AI       0.2864887713087e+04
#define AW       0.3135012829948e+04
#define BI       0.1660931315020e+00
#define BW       0.2367075766316e+01
#define CHLF     0.3340000000000e+06
#define CHLC     0.2501000000000e+07
#define CONSOL   0.1367000000000e+04
#define CONSOL2  0.1361000000000e+04
#define CPD      0.1005460000000e+04
#define CPV      0.1869460000000e+04
#define CPI      0.2115300000000e+04
#define CPW      0.4218000000000e+04
#define DELTA    0.6077686814144e+00
#define EPS1     0.6219800221014e+00
#define EPS2     0.3780199778986e+00
#define GRAV     0.9806160000000e+01
#define KARMAN   0.4000000000000e+00
#define KNAMS    0.5147910000000e+00
#define OMEGA    0.7292000000000e-04
#define PI       3.14159265358979323846
#define RAUW     0.1000000000000e+04
#define RAYT     0.6371220000000e+07
#define RGASD    0.2870500000000e+03
#define RGASV    0.4615100000000e+03
#define RIC      0.2000000000000e+00
#define SLP      0.6666666666667e-01
#define STEFAN   0.5669800000000e-07
#define STLO     0.6628486583943e-03
#define T1S      0.2731600000000e+03
#define T2S      0.2581600000000e+03
#define TCDK     0.2731500000000e+03
#define TGL      0.2731600000000e+03
#define TRPL     0.2731600000000e+03

#define CAPPA    (RGASD/CPD)
#define CVD      (CPD - RGASD)

#ifdef __FORTRAN__

#define TTNS1    610.78d0
#define TTNS3W    17.269d0
#define TTNS3I    21.875d0
#define TTNS4W    35.86d0
#define TTNS4I     7.66d0
#define AERK1W   610.94d0
#define AERK2W    17.625d0
#define AERK3W    30.11d0
#define AERK1I   611.21d0
#define AERK2I    22.587d0
#define AERK3I    -0.71d0

#else

#define TTNS1    610.78
#define TTNS3W    17.269
#define TTNS3I    21.875
#define TTNS4W    35.86
#define TTNS4I     7.66
#define AERK1W   610.94
#define AERK2W    17.625
#define AERK3W    30.11
#define AERK1I   611.21
#define AERK2I    22.587
#define AERK3I    -0.71

#endif

#endif

#ifdef __FORTRAN__
#define _TDVAR4(N, V) real,   parameter :: N = V
#define _TDVAR8(N, V) real(REAL64), parameter :: N = V

   _TDVAR4(ai, AI)
   _TDVAR4(aw, AW)
   _TDVAR4(bi, BI)
   _TDVAR4(bw, BW)
   _TDVAR4(cappa, CAPPA)
   _TDVAR4(chlf, CHLF)
   _TDVAR4(chlc, CHLC)
   _TDVAR4(consol, CONSOL)
   _TDVAR4(consol2, CONSOL2)
   _TDVAR4(cpd, CPD)
   _TDVAR4(cpv, CPV)
   _TDVAR4(cpi, CPI)
   _TDVAR4(cpw, CPW)
   _TDVAR4(cvd, CVD)
   _TDVAR4(delta, DELTA)
   _TDVAR4(eps1, EPS1)
   _TDVAR4(eps2, EPS2)
   _TDVAR4(grav, GRAV)
   _TDVAR4(karman, KARMAN)
   _TDVAR4(knams, KNAMS)
   _TDVAR4(omega, OMEGA)
   _TDVAR4(pi, PI)
   _TDVAR4(rauw, RAUW)
   _TDVAR4(rayt, RAYT)
   _TDVAR4(rgasd, RGASD)
   _TDVAR4(rgasv, RGASV)
   _TDVAR4(ric, RIC)
   _TDVAR4(slp, SLP)
   _TDVAR4(stefan, STEFAN)
   _TDVAR4(stlo, STLO)
   _TDVAR4(t1s, T1S)
   _TDVAR4(t2s, T2S)
   _TDVAR4(tcdk, TCDK)
   _TDVAR4(tgl, TGL)
   _TDVAR4(trpl, TRPL)

   _TDVAR8(ai_8, AI)
   _TDVAR8(aw_8, AW)
   _TDVAR8(bi_8, BI)
   _TDVAR8(bw_8, BW)
   _TDVAR8(cappa_8, CAPPA)
   _TDVAR8(chlf_8, CHLF)
   _TDVAR8(chlc_8, CHLC)
   _TDVAR8(consol_8, CONSOL)
   _TDVAR8(consol2_8, CONSOL2)
   _TDVAR8(cpd_8, CPD)
   _TDVAR8(cpv_8, CPV)
   _TDVAR8(cpi_8, CPI)
   _TDVAR8(cpw_8, CPW)
   _TDVAR8(cvd_8, CVD)
   _TDVAR8(delta_8, DELTA)
   _TDVAR8(eps1_8, EPS1)
   _TDVAR8(eps2_8, EPS2)
   _TDVAR8(grav_8, GRAV)
   _TDVAR8(karman_8, KARMAN)
   _TDVAR8(knams_8, KNAMS)
   _TDVAR8(omega_8, OMEGA)
   _TDVAR8(pi_8, PI)
   _TDVAR8(rauw_8, RAUW)
   _TDVAR8(rayt_8, RAYT)
   _TDVAR8(rgasd_8, RGASD)
   _TDVAR8(rgasv_8, RGASV)
   _TDVAR8(ric_8, RIC)
   _TDVAR8(slp_8, SLP)
   _TDVAR8(stefan_8, STEFAN)
   _TDVAR8(stlo_8, STLO)
   _TDVAR8(t1s_8, T1S)
   _TDVAR8(t2s_8, T2S)
   _TDVAR8(tcdk_8, TCDK)
   _TDVAR8(tgl_8, TGL)
   _TDVAR8(trpl_8, TRPL)

   _TDVAR8(ttns1, TTNS1)
   _TDVAR8(ttns3w, TTNS3W)
   _TDVAR8(ttns3i, TTNS3I)
   _TDVAR8(ttns4w, TTNS4W)
   _TDVAR8(ttns4i, TTNS4I)
   _TDVAR8(aerk1w, AERK1W)
   _TDVAR8(aerk2w, AERK2W)
   _TDVAR8(aerk3w, AERK3W)
   _TDVAR8(aerk1i, AERK1I)
   _TDVAR8(aerk2i, AERK2I)
   _TDVAR8(aerk3i, AERK3I)

#endif
