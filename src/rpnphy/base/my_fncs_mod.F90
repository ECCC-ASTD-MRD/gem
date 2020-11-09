module my_fncs_mod
   use, intrinsic :: iso_fortran_env, only: REAL64

!==============================================================================!
!  The following functions are used by the schemes in the multimoment package. !
!                                                                              !
!  Package version:  2.18.0      (internal bookkeeping)                        !
!  Last modified  :  2009-04-27                                                !
!==============================================================================!

   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   private

   public  :: NccnFNC,SxFNC,gamma,NccnFNC_v33,SxFNC_v33,gammaDP,diagAlpha_v33, &
              solveAlpha_v33,gser,gammln,gammp,cfg,gamminc

   contains

!==============================================================================!

 real function NccnFNC(Win,Tin,Pin,CCNtype)

!---------------------------------------------------------------------------!
! This function returns number concentration (activated aerosols) as a
! function of w,T,p, based on polynomial approximations of detailed
! approach using a hypergeometric function, following Cohard and Pinty (2000a).
!---------------------------------------------------------------------------!

  implicit none

! PASSING PARAMETERS:
  real,    intent(in) :: Win, Tin, Pin
  integer, intent(in) :: CCNtype

! LOCAL PARAMETERS:
  real :: T,p,x,y,a,b,c,d,e,f,g,h,T2,T3,T4,x2,x3,x4,p2

  x= log10(Win*100.);   x2= x*x;  x3= x2*x;  x4= x2*x2
  T= Tin - 273.15;      T2= T*T;  T3= T2*T;  T4= T2*T2
  p= Pin*0.01;          p2= p*p

  if (CCNtype==1) then  !Maritime

     a= 1.47e-9*T4 -6.944e-8*T3 -9.933e-7*T2 +2.7278e-4*T -6.6853e-4
     b=-1.41e-8*T4 +6.662e-7*T3 +4.483e-6*T2 -2.0479e-3*T +4.0823e-2
     c= 5.12e-8*T4 -2.375e-6*T3 +4.268e-6*T2 +3.9681e-3*T -3.2356e-1
     d=-8.25e-8*T4 +3.629e-6*T3 -4.044e-5*T2 +2.1846e-3*T +9.1227e-1
     e= 5.02e-8*T4 -1.973e-6*T3 +3.944e-5*T2 -9.0734e-3*T +1.1256e0
     f= -1.424e-6*p2 +3.631e-3*p -1.986
     g= -0.0212*x4 +0.1765*x3 -0.3770*x2 -0.2200*x +1.0081
     h= 2.47e-6*T3 -3.654e-5*T2 +2.3327e-3*T +0.1938
     y= a*x4 + b*x3 + c*x2 + d*x + e + f*g*h
     NccnFNC= 10.**min(2.,max(0.,y)) *1.e6                ![m-3]

  else if (CCNtype==2) then  !Continental

     a= 0.
     b= 0.
     c=-2.112e-9*T4 +3.9836e-8*T3 +2.3703e-6*T2 -1.4542e-4*T -0.0698
     d=-4.210e-8*T4 +5.5745e-7*T3 +1.8460e-5*T2 +9.6078e-4*T +0.7120
     e= 1.434e-7*T4 -1.6455e-6*T3 -4.3334e-5*T2 -7.6720e-3*T +1.0056
     f= 1.340e-6*p2 -3.5114e-3*p  +1.9453
     g= 4.226e-3*x4 -5.6012e-3*x3 -8.7846e-2*x2 +2.7435e-2*x +0.9932
     h= 5.811e-9*T4 +1.5589e-7*T3 -3.8623e-5*T2 +1.4471e-3*T +0.1496
     y= a*x4 +b*x3 +c*x2 + d*x + e + (f*g*h)
     NccnFNC= 10.**max(0.,y) *1.e6

  else

    print*, '*** STOPPED in MODULE ### NccnFNC  *** '
    print*, '    Parameter CCNtype incorrectly specified'
    stop

  endif

 end function NccnFNC
!======================================================================!

   real function SxFNC(Win,Tin,Pin,Qsw,Qsi,CCNtype,WRT)

!---------------------------------------------------------------------------!
! This function returns the peak supersaturation achieved during ascent with
! activation of CCN aerosols as a function of w,T,p, based on polynomial
! approximations of detailed approach using a hypergeometric function,
! following Cohard and Pinty (2000a).
!---------------------------------------------------------------------------!

 implicit none

! PASSING PARAMETERS:
  integer, intent(IN) :: WRT
  integer, intent(IN) :: CCNtype
  real,    intent(IN) :: Win, Tin, Pin, Qsw, Qsi

! LOCAL PARAMETERS:
  real   ::  Si,Sw,Qv,T,p,x,a,b,c,d,f,g,h,Pcorr,T2corr,T2,T3,T4,x2,x3,x4,p2
  real, parameter :: TRPL= 273.15

  x= log10(max(Win,1.e-20)*100.);   x2= x*x;  x3= x2*x;  x4= x2*x2
  T= Tin;                           T2= T*T;  T3= T2*T;  T4= T2*T2
  p= Pin*0.01;                      p2= p*p

  if (CCNtype==1) then  !Maritime

     a= -5.109e-7*T4 -3.996e-5*T3 -1.066e-3*T2 -1.273e-2*T +0.0659
     b=  2.014e-6*T4 +1.583e-4*T3 +4.356e-3*T2 +4.943e-2*T -0.1538
     c= -2.037e-6*T4 -1.625e-4*T3 -4.541e-3*T2 -5.118e-2*T +0.1428
     d=  3.812e-7*T4 +3.065e-5*T3 +8.795e-4*T2 +9.440e-3*T +6.14e-3
     f= -2.012e-6*p2 + 4.1913e-3*p    - 1.785e0
     g=  2.832e-1*x3 -5.6990e-1*x2 +5.1105e-1*x -4.1747e-4
     h=  1.173e-6*T3 +3.2174e-5*T2 -6.8832e-4*T +6.7888e-2
     Pcorr= f*g*h
     T2corr= 0.9581-4.449e-3*T-2.016e-4*T2-3.307e-6*T3-1.725e-8*T4

  else if (CCNtype==2) then  !Continental [computed for -35<T<-5C]

     a=  3.80e-5*T2 +1.65e-4*T +9.88e-2
     b= -7.38e-5*T2 -2.53e-3*T -3.23e-1
     c=  8.39e-5*T2 +3.96e-3*T +3.50e-1
     d= -1.88e-6*T2 -1.33e-3*T -3.73e-2
     f= -1.9761e-6*p2 + 4.1473e-3*p - 1.771e0
     g=  0.1539*x4 -0.5575*x3 +0.9262*x2 -0.3498*x -0.1293
     h=-8.035e-9*T4+3.162e-7*T3+1.029e-5*T2-5.931e-4*T+5.62e-2
     Pcorr= f*g*h
     T2corr= 0.98888-5.0525e-4*T-1.7598e-5*T2-8.3308e-8*T3

  else

    print*, '*** STOPPED in MODULE ### SxFNC  *** '
    print*, '    Parameter CCNtype incorrectly specified'
    stop

  endif

  Sw= (a*x3 + b*x2 +c*x + d) + Pcorr
  Sw= 1. + 0.01*Sw
  Qv= Qsw*Sw
  Si= Qv/Qsi
  Si= Si*T2corr
  if (WRT.eq.1) then
     SxFNC= Sw
  else
     SxFNC= Si
  endif
  if (Win.le.0.) SxFNC= 1.

 end function SxFNC
!======================================================================!

 real function gamma(xx)

!  Modified from "Numerical Recipes"

  implicit none

! PASSING PARAMETERS:
  real, intent(IN) :: xx

! LOCAL PARAMETERS:
  integer  :: j
  real(REAL64) :: ser,stp,tmp,x,y,cof(6),gammadp


  save cof,stp
  data cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=dble(xx)
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammadp=tmp+log(stp*ser/x)
  gammadp= exp(gammadp)

  gamma  = sngl(gammadp)

 end function gamma
!======================================================================!
! ! !
! ! ! -- USED BY DIAGNOSTIC-ALPHA DOUBLE-MOMENT (SINGLE-PRECISION) VERSION --
! ! !      FOR FUTURE VERSIONS OF M-Y PACKAGE WITH, THIS S/R CAN BE USED
! ! !
! ! !  real FUNCTION diagAlpha(Dm,x)
! ! !
! ! !   implicit none
! ! !
! ! !   integer :: x
! ! !   real    :: Dm
! ! !   real, dimension(5) :: c1,c2,c3,c4
! ! !   real, parameter    :: pi = 3.14159265
! ! !   real, parameter    :: alphaMAX= 80.e0
! ! !   data c1 /19.0, 12.0, 4.5, 5.5, 3.7/
! ! !   data c2 / 0.6,  0.7, 0.5, 0.7, 0.3/
! ! !   data c3 / 1.8,  1.7, 5.0, 4.5, 9.0/
! ! !   data c4 /17.0, 11.0, 5.5, 8.5, 6.5/
! ! !   diagAlpha= c1(x)*tanh(c2(x)*(1.e3*Dm-c3(x)))+c4(x)
! ! !   if (x==5.and.Dm>0.008) diagAlpha= 1.e3*Dm-2.6
! ! !   diagAlpha= min(diagAlpha, alphaMAX)
! ! !
! ! !  END function diagAlpha
! ! !
! ! ! !======================================================================!
! ! !
! ! ! -- USED BY DIAGNOSTIC-ALPHA DOUBLE-MOMENT (SINGLE-PRECISION) VERSION --
! ! !      FOR FUTURE VERSIONS OF M-Y PACKAGE WITH, THIS S/R CAN BE USED
! ! !
! ! !  real FUNCTION solveAlpha(Q,N,Z,Cx,rho)
! ! !
! ! !  implicit none
! ! !
! ! ! ! PASSING PARAMETERS:
! ! !   real, intent(IN) :: Q, N, Z, Cx, rho
! ! !
! ! ! ! LOCAL PARAMETERS:
! ! !   real             :: a,g,a1,g1,g2,tmp1
! ! !   integer          :: i
! ! !   real, parameter  :: alphaMax= 40.
! ! !   real, parameter  :: epsQ    = 1.e-14
! ! !   real, parameter  :: epsN    = 1.e-3
! ! !   real, parameter  :: epsZ    = 1.e-32
! ! !
! ! ! !  Q         mass mixing ratio
! ! ! !  N         total concentration
! ! ! !  Z         reflectivity
! ! ! !  Cx        (pi/6)*RHOx
! ! ! !  rho       air density
! ! ! !  a         alpha (returned as solveAlpha)
! ! ! !  g         function g(a)= [(6+a)(5+a)(4+a)]/[(3+a)(2+a)(1+a)],
! ! ! !              where g = (Cx/(rho*Q))**2.*(Z*N)
! ! !
! ! !
! ! !   if (Q==0. .or. N==0. .or. Z==0. .or. Cx==0. .or. rho==0.) then
! ! !   ! For testing/debugging only; this module should never be called
! ! !   ! if the above condition is true.
! ! !     print*,'*** STOPPED in MODULE ### solveAlpha *** '
! ! !     print*,'*** : ',Q,N,Z,Cx*1.9099,rho
! ! !     stop
! ! !   endif
! ! !
! ! !   IF (Q>epsQ .and. N>epsN .and. Z>epsZ ) THEN
! ! !
! ! !      tmp1= Cx/(rho*Q)
! ! !      g   = tmp1*Z*tmp1*N    ! g = (Z*N)*[Cx / (rho*Q)]^2
! ! !
! ! !  !Note: The above order avoids OVERFLOW, since tmp1*tmp1 is very large
! ! !
! ! ! !----------------------------------------------------------!
! ! ! ! !Solve alpha numerically: (brute-force; for testing only)
! ! ! !      a= 0.
! ! ! !      g2= 999.
! ! ! !      do i=0,4000
! ! ! !         a1= i*0.01
! ! ! !         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
! ! ! !         if(abs(g-g1)<abs(g-g2)) then
! ! ! !            a = a1
! ! ! !            g2= g1
! ! ! !         endif
! ! ! !      enddo
! ! ! !----------------------------------------------------------!
! ! !
! ! ! !Piecewise-polynomial approximation of g(a) to solve for a:  [2004-11-29]
! ! !      if (g>=20.) then
! ! !        a= 0.
! ! !      else
! ! !        g2= g*g
! ! !        if (g<20.  .and.g>=13.31) a= 3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0
! ! !        if (g<13.31.and.g>=7.123) a= 1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0
! ! !        if (g<7.123.and.g>=4.200) a= 1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0
! ! !        if (g<4.200.and.g>=2.946) a= 5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1
! ! !        if (g<2.946.and.g>=1.793) a= 4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1
! ! !        if (g<1.793.and.g>=1.405) a= 4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2
! ! !        if (g<1.405.and.g>=1.230) a= 3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2
! ! !        if (g<1.230) a= alphaMax
! ! !      endif
! ! !
! ! !      solveAlpha= max(0.,min(a,alphaMax))
! ! !
! ! !   ELSE
! ! !
! ! !      solveAlpha= 0.
! ! !
! ! !   ENDIF
! ! !
! ! !  END FUNCTION solveAlpha
!======================================================================!

! The following functions are used only by 'my_main_full.ftn90'.  They are somewhat
! redundant from above routines, though there are small differences.  Eventually,
! 'my_main_full.ftn90' should be modified to use the same functions as the other
! versions of the scheme.
! 2008-04-15

!======================================================================!
 real function NccnFNC_v33(Win,Tin,Pin,AIRTYPE)

!---------------------------------------------------------------------------!
! This function returns number concentration (activated aerosols) as a
! function of w,T,p, based on polynomial approximations of detailed
! approach using a hypergeometric function, following Cohard and Pinty (2000a).
!---------------------------------------------------------------------------!

  implicit none

! PASSING PARAMETERS:
  real,    intent(IN) :: Win, Tin, Pin
  integer, intent(IN) :: AIRTYPE

! LOCAL PARAMETERS:
  real :: T,p,x,y,a,b,c,d,e,f,g,h,T2,T3,T4,x2,x3,x4,p2

  x= log10(Win*100.);   x2= x*x;  x3= x2*x;  x4= x2*x2
  T= Tin - 273.15;      T2= T*T;  T3= T2*T;  T4= T2*T2
  p= Pin*0.01;          p2= p*p

  if (AIRTYPE==1) then  !Maritime

     a= 1.47e-9*T4 -6.944e-8*T3 -9.933e-7*T2 +2.7278e-4*T -6.6853e-4
     b=-1.41e-8*T4 +6.662e-7*T3 +4.483e-6*T2 -2.0479e-3*T +4.0823e-2
     c= 5.12e-8*T4 -2.375e-6*T3 +4.268e-6*T2 +3.9681e-3*T -3.2356e-1
     d=-8.25e-8*T4 +3.629e-6*T3 -4.044e-5*T2 +2.1846e-3*T +9.1227e-1
     e= 5.02e-8*T4 -1.973e-6*T3 +3.944e-5*T2 -9.0734e-3*T +1.1256e0
     f= -1.424e-6*p2 +3.631e-3*p -1.986
     g= -0.0212*x4 +0.1765*x3 -0.3770*x2 -0.2200*x +1.0081
     h= 2.47e-6*T3 -3.654e-5*T2 +2.3327e-3*T +0.1938
     y= a*x4 + b*x3 + c*x2 + d*x + e + f*g*h
     NccnFNC_v33= 10.**min(2.,max(0.,y)) *1.e6                ![m-3]

  else if (AIRTYPE==2) then  !Continental

     a= 0.
     b= 0.
     c=-2.112e-9*T4 +3.9836e-8*T3 +2.3703e-6*T2 -1.4542e-4*T -0.0698
     d=-4.210e-8*T4 +5.5745e-7*T3 +1.8460e-5*T2 +9.6078e-4*T +0.7120
     e= 1.434e-7*T4 -1.6455e-6*T3 -4.3334e-5*T2 -7.6720e-3*T +1.0056
     f= 1.340e-6*p2 -3.5114e-3*p  +1.9453
     g= 4.226e-3*x4 -5.6012e-3*x3 -8.7846e-2*x2 +2.7435e-2*x +0.9932
     h= 5.811e-9*T4 +1.5589e-7*T3 -3.8623e-5*T2 +1.4471e-3*T +0.1496
     y= a*x4 +b*x3 +c*x2 + d*x + e + (f*g*h)
     NccnFNC_v33= 10.**max(0.,y) *1.e6

  else

    print*, '*** STOPPED in MODULE ### NccnFNC  *** '
    print*, '    Parameter AIRTYPE incorrectly specified'
    stop

  endif

 end function NccnFNC_v33
!======================================================================!

   real(REAL64) function SxFNC_v33(Win,Tin,Pin,Qsw,Qsi,AIRTYPE,WRT)

!---------------------------------------------------------------------------!
! This function returns the peak supersaturation achieved during ascent with
! activation of CCN aerosols as a function of w,T,p, based on polynomial
! approximations of detailed approach using a hypergeometric function,
! following Cohard and Pinty (2000a).
!---------------------------------------------------------------------------!

 implicit none

! PASSING PARAMETERS:
  integer, intent(IN) :: WRT
  integer, intent(IN) :: AIRTYPE
  real,    intent(IN) :: Win, Tin, Pin, Qsw, Qsi

! LOCAL PARAMETERS:
  real   ::  FOQSA,FOQST,Si,Sw,Qv,T,p,x,a,b,c,d,f,g,h,Pcorr,T2corr,   &
             T2,T3,T4,x2,x3,x4,p2
  real, parameter :: TRPL= 273.15

  x= log10(max(Win,1.e-20)*100.);   x2= x*x;  x3= x2*x;  x4= x2*x2
  T= Tin;                           T2= T*T;  T3= T2*T;  T4= T2*T2
  p= Pin*0.01;                      p2= p*p

  if (AIRTYPE==1) then  !Maritime

     a= -5.109e-7*T4 -3.996e-5*T3 -1.066e-3*T2 -1.273e-2*T +0.0659
     b=  2.014e-6*T4 +1.583e-4*T3 +4.356e-3*T2 +4.943e-2*T -0.1538
     c= -2.037e-6*T4 -1.625e-4*T3 -4.541e-3*T2 -5.118e-2*T +0.1428
     d=  3.812e-7*T4 +3.065e-5*T3 +8.795e-4*T2 +9.440e-3*T +6.14e-3
     f= -2.012e-6*p2 + 4.1913e-3*p    - 1.785e0
     g=  2.832e-1*x3 -5.6990e-1*x2 +5.1105e-1*x -4.1747e-4
     h=  1.173e-6*T3 +3.2174e-5*T2 -6.8832e-4*T +6.7888e-2
     Pcorr= f*g*h
     T2corr= 0.9581-4.449e-3*T-2.016e-4*T2-3.307e-6*T3-1.725e-8*T4

  else if (AIRTYPE==2) then  !Continental [computed for -35<T<-5C]

     a=  3.80e-5*T2 +1.65e-4*T +9.88e-2
     b= -7.38e-5*T2 -2.53e-3*T -3.23e-1
     c=  8.39e-5*T2 +3.96e-3*T +3.50e-1
     d= -1.88e-6*T2 -1.33e-3*T -3.73e-2
     f= -1.9761e-6*p2 + 4.1473e-3*p - 1.771e0
     g=  0.1539*x4 -0.5575*x3 +0.9262*x2 -0.3498*x -0.1293
     h=-8.035e-9*T4+3.162e-7*T3+1.029e-5*T2-5.931e-4*T+5.62e-2
     Pcorr= f*g*h
     T2corr= 0.98888-5.0525e-4*T-1.7598e-5*T2-8.3308e-8*T3

  else

    print*, '*** STOPPED in MODULE ### SxFNC  *** '
    print*, '    Parameter AIRTYPE incorrectly specified'
    stop

  endif

  Sw= (a*x3 + b*x2 +c*x + d) + Pcorr
  Sw= 1. + 0.01*Sw
  Qv= Qsw*Sw
  Si= Qv/Qsi
  Si= Si*T2corr
  if (WRT.eq.1) then
     SxFNC_v33= Sw
  else
     SxFNC_v33= Si
  endif
  if (Win.le.0.) SxFNC_v33= 1.

 end function SxFNC_v33
!======================================================================!

 function gammaDP(xx)

!  Modified from "Numerical Recipes"

  implicit none

! PASSING PARAMETERS:
  double precision, intent(IN) :: xx

! LOCAL PARAMETERS:
  double precision  :: gammaDP
  integer  :: j
  double precision  :: ser,stp,tmp,x,y,cof(6)


  save cof,stp
  data cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammaDP=tmp+log(stp*ser/x)
  gammaDP= exp(gammaDP)

 end function gammaDP
!======================================================================!

 function diagAlpha_v33(Dm,x)

  implicit none

  integer :: x
  real(REAL64) :: diagAlpha_v33,Dm
  real(REAL64), dimension(5) :: c1,c2,c3,c4
  real(REAL64), parameter    :: pi = 3.14159265d0
  real(REAL64), parameter    :: alphaMAX= 80.d0
  data c1 /19.0d0, 12.0d0, 4.5d0, 5.5d0, 3.7d0/
  data c2 / 0.6d0,  0.7d0, 0.5d0, 0.7d0, 0.3d0/
  data c3 / 1.8d0,  1.7d0, 5.0d0, 4.5d0, 9.0d0/
  data c4 /17.0d0, 11.0d0, 5.5d0, 8.5d0, 6.5d0/
  diagAlpha_v33= c1(x)*tanh(c2(x)*(1.d3*Dm-c3(x)))+c4(x)
  if (x==5.and.Dm>0.008d0) diagAlpha_v33= 1.d3*Dm-2.6d0
  diagAlpha_v33= min(diagAlpha_v33, alphaMAX)

 end function diagAlpha_v33

!======================================================================!

 function solveAlpha_v33(Q,N,Z,Cx,rho)

 implicit none

! PASSING PARAMETERS:
  real, intent(IN) :: Q, N, Z, Cx, rho

! LOCAL PARAMETERS:
  real(REAL64) :: solveAlpha_v33
  real   :: a,g,a1,g1,g2,tmp1
  integer :: i
  real, parameter :: alphaMax= 40.
  real, parameter :: epsQ    = 1.e-14
  real, parameter :: epsN    = 1.e-3
  real, parameter :: epsZ    = 1.e-32

!  Q         mass mixing ratio
!  N         total concentration
!  Z         reflectivity
!  Cx        (pi/6)*RHOx
!  rho       air density
!  a         alpha (returned as solveAlpha)
!  g         function g(a)= [(6+a)(5+a)(4+a)]/[(3+a)(2+a)(1+a)],
!              where g = (Cx/(rho*Q))**2.*(Z*N)


  if (Q==0. .or. N==0. .or. Z==0. .or. Cx==0. .or. rho==0.) then
  ! For testing/debugging only; this module should never be called
  ! if the above condition is true.
    print*,'*** STOPPED in MODULE ### solveAlpha *** '
    print*,'*** : ',Q,N,Z,Cx*1.9099,rho
    stop
  endif

  if (Q>epsQ .and. N>epsN .and. Z>epsZ ) then

     tmp1= Cx/(rho*Q)
     g   = tmp1*Z*tmp1*N    ! g = (Z*N)*[Cx / (rho*Q)]^2

 !Note: The above order avoids OVERFLOW, since tmp1*tmp1 is very large

!----------------------------------------------------------!
! !Solve alpha numerically: (brute-force; for testing only)
!      a= 0.
!      g2= 999.
!      do i=0,4000
!         a1= i*0.01
!         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
!         if(abs(g-g1)<abs(g-g2)) then
!            a = a1
!            g2= g1
!         endif
!      enddo
!----------------------------------------------------------!

!Piecewise-polynomial approximation of g(a) to solve for a:  [2004-11-29]
     if (g>=20.) then
       a= 0.
     else
       g2= g*g
       if (g<20.  .and.g>=13.31) a= 3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0
       if (g<13.31.and.g>=7.123) a= 1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0
       if (g<7.123.and.g>=4.200) a= 1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0
       if (g<4.200.and.g>=2.946) a= 5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1
       if (g<2.946.and.g>=1.793) a= 4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1
       if (g<1.793.and.g>=1.405) a= 4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2
       if (g<1.405.and.g>=1.230) a= 3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2
       if (g<1.230) a= alphaMax
     endif

     solveAlpha_v33= max(0.,min(a,alphaMax))

  else

     solveAlpha_v33= 0.

  endif

 end function solveAlpha_v33

!======================================================================!

 subroutine gser(gamser,a,x,gln)

! USES gammln

!   Returns the incomplete gamma function P(a,x) evaluated by its series
!   representation as gamser.  Also returns GAMMA(a) as gln.

 implicit none

 integer :: itmax
 real    :: a,gamser,gln,x,eps
 parameter (itmax=100, eps=3.e-7)
 integer :: n
 real :: ap,de1,summ
!!$ real, external :: gammln

 gln=gammln(a)
 if(x.le.0.)then
    if(x.lt.0.)pause 'x <0 in gser'
    gamser=0.
    return
 endif
 ap=a
 summ=1./a
 de1=summ
 do n=1,itmax
    ap=ap+1.
    de1=de1*x/ap
    summ=summ+de1
    if(abs(de1).lt.abs(summ)*eps) goto 1
 enddo
 pause 'a too large, itmax too small in gser'
1 gamser=summ*exp(-x+a*log(x)-gln)
 return

end subroutine gser
!======================================================================!

 real function gammln(xx)

!  Returns value of ln(GAMMA(xx)) for xx>0
!   (modified from "Numerical Recipes")

  implicit none

! PASSING PARAMETERS:
  real, intent(IN) :: xx

! LOCAL PARAMETERS:
  integer  :: j
  real(REAL64) :: ser,stp,tmp,x,y,cof(6),gammadp

  save cof,stp
  data cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=dble(xx)
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6   !original
!  do j=1,4
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammln= sngl( tmp+log(stp*ser/x)  )

 end function gammln
!======================================================================!

 real function gammp(a,x)

! USES gcf,gser

! Returns the incomplete gamma function P(a,x)

 implicit none

 real :: a,x,gammcf,gamser,gln

 if(x.lt.0..or.a.le.0.) pause 'bad arguments in gammq'
 if(x.lt.a+1.)then
    call gser(gamser,a,x,gln)
    gammp=gamser
 else
    call cfg(gammcf,a,x,gln)
    gammp=1.-gammcf
 endif
 return

 end function gammp
!======================================================================!

 subroutine cfg(gammcf,a,x,gln)

! USES gammln

! Returns the incomplete gamma function (Q(a,x) evaluated by tis continued fraction
! representation as gammcf.  Also returns ln(GAMMA(a)) as gln.  ITMAX is the maximum
! allowed number of iterations; EPS is the relative accuracy; FPMIN is a number near
! the smallest representable floating-point number.

 implicit none

 integer :: i
 real    :: a,gammcf,gln,x,fpmin
 real    :: an,b,c,d,de1,h
 integer, parameter :: itmax=100
 real, parameter :: eps=3.e-7
!!$ real, external :: gammln

 gln=gammln(a)
 b=x+1.-a
 c=1./fpmin
 d=1./b
 h=d
 do i= 1,itmax
   an=-i*(i-a)
   b=b+2.
   d=an*d+b
   if(abs(d).lt.fpmin)d=fpmin
   c=b+an/c
 if(abs(c).lt.fpmin) c=fpmin
   d=1./d
   de1=d*c
   h=h*de1
   if(abs(de1-1.).lt.eps) goto 1
 enddo
 pause 'a too large, itmax too small in gcf'
1 gammcf=exp(-x+a*log(x)-gln)*h
 return

end subroutine cfg
!======================================================================!

 real function gamminc(p,xmax)

! USES gammp, gammln
! Returns incomplete gamma function, gamma(p,xmax)= P(p,xmax)*GAMMA(p)
 real :: p,xmax
 gamminc= gammp(p,xmax)*exp(gammln(p))

 end function gamminc

!======================================================================!
!  real function x_tothe_y(x,y)
!
!     implicit none
!     real, intent(in) :: x,y
!     x_tothe_y= exp(y*log(x))
!
!  end function x_tothe_y
!======================================================================!

end module my_fncs_mod
