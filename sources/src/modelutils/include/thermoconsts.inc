!module module_thermoconsts
!implicit none

  real, parameter ::  AI     = 0.2864887713087e+04  ! pour fn htvocp
  real, parameter ::  AW     = 0.3135012829948e+04  ! pour fn htvocp
  real, parameter ::  BI     = 0.1660931315020e+00  ! pour fn htvocp
  real, parameter ::  BW     = 0.2367075766316e+01  ! pour fn htvocp
  real, parameter ::  CAPPA  = 0.2854912179500e+00  ! rgasd/cpd
  real, parameter ::  CHLF   = 0.3340000000000e+06  ! ch. lat. fusion       J kg-1
  real, parameter ::  CHLC   = 0.2501000000000e+07  ! ch. lat. condens.(0C) J kg-1
  real, parameter ::  CONSOL = 0.1367000000000e+04  ! constante solaire     W m-2
  real, parameter ::  CONSOL2= 0.1361000000000e+04  ! constante solaire (corrigee)    W m-2
  real, parameter ::  CPD    = 0.1005460000000e+04  ! chal. spec. air sec   J kg-1 K-1
  real, parameter ::  CPV    = 0.1869460000000e+04  ! chal. spec. vap eau   J kg-1 K-1
  real, parameter ::  CPI    = 0.2115300000000e+04  ! chal. spec. glace     J kg-1 K-1
  real, parameter ::  DELTA  = 0.6077686814144e+00  ! 1/eps1 - 1
  real, parameter ::  EPS1   = 0.6219800221014e+00  ! rgasd/rgasv
  real, parameter ::  EPS2   = 0.3780199778986e+00  ! 1 - eps1
  real, parameter ::  GRAV   = 0.9806160000000e+01  ! acc. de gravite       m s-2
  real, parameter ::  KARMAN = 0.4000000000000e+00  ! cte de von karman
  real, parameter ::  KNAMS  = 0.5147910000000e+00  ! passage kt a m/s
  real, parameter ::  OMEGA  = 0.7292000000000e-04  ! rotation terre        s-1
  real, parameter ::  PI     = 0.3141592653590e+01  ! cte pi=acos(-1)
  real, parameter ::  RAUW   = 0.1000000000000e+04  ! densite eau liq       kg m-3
  real, parameter ::  RAYT   = 0.6371220000000e+07  ! rayon moy. terre      m
  real, parameter ::  RGASD  = 0.2870500000000e+03  ! cte gaz - air sec     J kg-1 K-1
  real, parameter ::  RGASV  = 0.4615100000000e+03  ! cte gaz - vap eau     J kg-1 K-1
  real, parameter ::  RIC    = 0.2000000000000e+00  ! cte richardson crit.
  real, parameter ::  SLP    = 0.6666666666667e-01  ! pour fn htvocp
  real, parameter ::  STEFAN = 0.5669800000000e-07  ! cte stefan-boltzmann  J m-2 s-1 K-4
  real, parameter ::  STLO   = 0.6628486583943e-03  ! schuman-newell l.r.   K s2 m-2
  real, parameter ::  T1S    = 0.2731600000000e+03  ! pour fn htvocp        K
  real, parameter ::  T2S    = 0.2581600000000e+03  ! pour fn htvocp        K
  real, parameter ::  TCDK   = 0.2731500000000e+03  ! passage k a c         C
  real, parameter ::  TGL    = 0.2731600000000e+03  ! temp glace dans atm   K
  real, parameter ::  TRPL   = 0.2731600000000e+03  ! point triple - eau    K

  real*8, parameter :: TTNS1  = 610.78D0
  real*8, parameter :: TTNS3W = 17.269D0
  real*8, parameter :: TTNS3I = 21.875D0
  real*8, parameter :: TTNS4W = 35.86D0
  real*8, parameter :: TTNS4I =  7.66D0

!  logical :: bit = .true.
!end module module_thermoconsts
