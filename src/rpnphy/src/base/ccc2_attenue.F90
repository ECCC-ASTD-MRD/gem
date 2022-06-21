!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
!**S/P ATTENUE - CALCULATES THE DOWNWARD FLUX ATTENUATION
!

      subroutine CCC2_attenue4 (atten, coef1, s1, dp, dip, dt, dt0, &
                          rmu, inpt, lc, isl, il1, il2, ilg)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lc, isl, il1, il2
      real atten(ilg), coef1(5,lc), dp(ilg), dip(ilg), &
           dt(ilg), dt0(ilg), rmu(ilg), s1(ilg)
      integer inpt(ilg)

!Author
!
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001   M. Lazare (May 05,2006) - previous version attenue2 for gcm15e/f:
!                               - implement rpn fix for inpt.
! 002   M.Lazare,J.Li (Dec 05,2007) - new version for gcm15g:
!                               - work array "s1" now local instead
!                                 of passed in.
!                               - methane contribution added through
!                                 ng=4.
! 003  M.Lazare(Apr 18,2008)    - previous version attenue3 for gcm15g:
!                               - cosmetic change of n->lc passed
!                                 in and used in dimension of coef
!                                 array(s), so that it is not redefined
!                                 and passed back out to calling
!                                 routine, changing value in data
!                                 statement!
!                               - cosmetic change to add threadprivate
!                                 for common block "trace", in support
!                                 of "radforce" model option.
! 004  M.Lazare,K.Winger(Apr 08) - correct bug - lc now used for
!                              dimension of array and n as a variable integer
! 005  J.Li (Feb 09,2009)       - new version for gcm15h:
!                               - 3d ghg implemented, thus no need
!                                 for "trace" common block or
!                                 temporary work arrays to hold
!                                 mixing ratios of ghg depending on
!                                 a passed, specified option.
! 006  P. Vaillancourt, A-M. Leduc. (Apr 2010). Integrating the Li/Lazare 2009 version. (attenu4).
!                                - Replace qq and O3 by S1. remove ng.

!Object
!
!        This subroutine calculates the downward flux attenuation above
!        the model top level
!        isl = 1 for solar, isl = 2 for infrared.
!        ng = 1, h2o; ng = 2, o3; ng = 3, co2; ng = 6, o2                     PLUS NECESSAIRE    AML
!        assuming the temperature at 0.0005 mb is 210 k
!
!Arguments
!
! atten  for solar: the attenuation factor for downward flux from
!        toa to the model top level; for longwave: the optical
!        / diffuse factor
! dp     air mass path for a model layer (explained in raddriv).
! s1     gas mass mixing ratio
! dip    interpolation factor for pressure between two
!        neighboring standard input data pressure levels
! dt     layer temperature - 250 k
! dt0    temperature in moon layer - 250 k
! rmu    cos of solar zenith angle
!
!Implicites
!
#include "ccc_tracegases.cdk"
!
!Modules
!
!*
!
      integer i,n,nm1
      real ru, x1, x2, tau
!
      data ru / 1.6487213 /
!***********************************************************************
!
      if (isl .eq. 1)                                               then
       if (inpt(1).lt.950)                                          then
        do 1000 i = il1, il2
          n =  inpt(i)
          nm1 =  max (n - 1, 1)
          x1       =   coef1(1,nm1) + dt0(i) * (coef1(2,nm1) + dt0(i) * &
                      (coef1(3,nm1) + dt0(i) * (coef1(4,nm1) + dt0(i) * &
                       coef1(5,1))))
          x2       =   coef1(1,n) + dt(i) * (coef1(2,n) + dt(i) * &
                      (coef1(3,n) + dt(i) * (coef1(4,n) + dt(i) * &
                       coef1(5,n))))
!
          tau      =  (x1 + (x2 - x1) * dip(i)) * 1.02 * s1(i) * dp(i)
          atten(i) =  exp(- tau / rmu(i))
 1000   continue
       else
        n =  inpt(1) - 1000
        nm1 =  max (n - 1, 1)
        do 1002 i = il1, il2
          x1       =   coef1(1,nm1) + dt0(i) * (coef1(2,nm1) + dt0(i) * &
                      (coef1(3,nm1) + dt0(i) * (coef1(4,nm1) + dt0(i) * &
                       coef1(5,1))))
          x2       =   coef1(1,n) + dt(i) * (coef1(2,n) + dt(i) * &
                      (coef1(3,n) + dt(i) * (coef1(4,n) + dt(i) * &
                       coef1(5,n))))
!
          tau      =  (x1 + (x2 - x1) * dip(i)) * 1.02 * s1(i) * dp(i)
          atten(i) =  exp(- tau / rmu(i))
 1002   continue
       endif

      else
       if (inpt(1).lt.950)                                          then
        do 2000 i = il1, il2
          n =  inpt(i)
          nm1 =  max (n - 1, 1)
          x1       =   coef1(1,nm1) + dt0(i) * (coef1(2,nm1) + dt0(i) * &
                      (coef1(3,nm1) + dt0(i) * (coef1(4,nm1) + dt0(i) * &
                       coef1(5,1))))
          x2       =   coef1(1,n) + dt(i) * (coef1(2,n) + dt(i) * &
                      (coef1(3,n) + dt(i) * (coef1(4,n) + dt(i) * &
                       coef1(5,n))))
!
          tau      =  (x1 + (x2 - x1) * dip(i)) * 1.02 * s1(i) * dp(i)
!
          atten(i) =   ru * tau
 2000   continue
       else
        n =  inpt(1) - 1000
        nm1 =  max (n - 1, 1)
        do 2002 i = il1, il2
          x1       =   coef1(1,nm1) + dt0(i) * (coef1(2,nm1) + dt0(i) * &
                      (coef1(3,nm1) + dt0(i) * (coef1(4,nm1) + dt0(i) * &
                       coef1(5,1))))
          x2       =   coef1(1,n) + dt(i) * (coef1(2,n) + dt(i) * &
                      (coef1(3,n) + dt(i) * (coef1(4,n) + dt(i) * &
                       coef1(5,n))))
!
          tau      =  (x1 + (x2 - x1) * dip(i)) * 1.02 * s1(i) * dp(i)
!
          atten(i) =   ru * tau
 2002   continue
       endif
      endif
!
      return
      end
