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
!**S/P TCONTL - INFRARED WATER VAPOR CONTINUUM
!
      subroutine ccc2_tcontl1 (taug, coef1, coef2, s1, dp, dip, dt, &
                         lc, inpt, mcont, gh, il1, il2, ilg, lay)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lc, il1, il2, mcont(ilg)
      real taug(ilg,lay), coef1(5,lc), coef2(5,lc)
      real s1(ilg,lay), dp(ilg,lay), dip(ilg,lay), dt(ilg,lay)
      integer inpt(ilg,lay)
      logical gh
!
!Authors
!
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001   M. Lazare (May 05,2006) - new version for gcm15e:
!                               - implement rpn fix for inpt.
! 002  P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the Li/Lazare 2009 version. (tcontl1)
!
!
!Object
!        Infrared water vapor continuum, coef1 is the coefficient for
!        self, coef2 is the coefficient for foreign. The continuum only
!        applies to the layers below 138.9440 mb or even lower region
!        depending on each band. lc is number of level for standard
!        pressure considered in calculating the continuum.
!        1.608 = 28.97 / 18.016, a fctr for water vapor partial pressure
!
!Arguments
!
! taug   gaseous optical depth
! s1     input h2o mixing ratio for each layer
! dp     air mass path for a model layer (explained in raddriv)
! dip    interpolation factor for pressure between two
!        neighboring standard input data pressure levels
! dt     layer temperature - 250 k
! inpt   number of the level for the standard input data pressures
! mcont  the highest level for water vapor continuum calculation
!
!*
      integer  k, i, j, m, n, nc
      real x1, y1, x2, y2

      if (gh)                                                       then
        nc =  29 - lc
      else
        nc =  19 - lc
      endif
!
      do 300 i = il1, il2
      do 200 k = mcont(i), lay
        if (inpt(1,k) .lt. 950)                                     then
            j =  inpt(i,k)
            if (j .ge. nc)                                          then
              m  =  j - nc + 1
              n  =  m + 1
              x1        =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                           dt(i,k) * (coef1(3,m) + dt(i,k) * &
                          (coef1(4,m) + dt(i,k) * coef1(5,m))))
!
              x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                           dt(i,k) * (coef1(3,n) + dt(i,k) * &
                          (coef1(4,n) + dt(i,k) * coef1(5,n))))
!
              y1        =  coef2(1,m) + dt(i,k) * (coef2(2,m) + &
                           dt(i,k) * (coef2(3,m) + dt(i,k) * &
                          (coef2(4,m) + dt(i,k) * coef2(5,m))))
!
              y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + &
                           dt(i,k) * (coef2(3,n) + dt(i,k) * &
                          (coef2(4,n) + dt(i,k) * coef2(5,n))))
!
              taug(i,k) =  taug(i,k) + &
                          ( (x1 - y1 + (x2 - x1 - y2 + y1) * &
                           dip(i,k)) * 1.608 * s1(i,k) + &
                             y1 + (y2 - y1) * dip(i,k) ) * &
                            s1(i,k) * dp(i,k)
            endif
        else
          j =  inpt(1,k) - 1000
          m  =  j - nc + 1
          n  =  m + 1
          if (j .ge. nc)                                          then
              x1        =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                           dt(i,k) * (coef1(3,m) + dt(i,k) * &
                          (coef1(4,m) + dt(i,k) * coef1(5,m))))
!
              x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                           dt(i,k) * (coef1(3,n) + dt(i,k) * &
                          (coef1(4,n) + dt(i,k) * coef1(5,n))))
!
              y1        =  coef2(1,m) + dt(i,k) * (coef2(2,m) + &
                           dt(i,k) * (coef2(3,m) + dt(i,k) * &
                          (coef2(4,m) + dt(i,k) * coef2(5,m))))
!
              y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + &
                           dt(i,k) * (coef2(3,n) + dt(i,k) * &
                          (coef2(4,n) + dt(i,k) * coef2(5,n))))
!
              taug(i,k) =  taug(i,k) + &
                          ( (x1 - y1 + (x2 - x1 - y2 + y1) * &
                           dip(i,k)) * 1.608 * s1(i,k) + &
                             y1 + (y2 - y1) * dip(i,k) ) * &
                            s1(i,k) * dp(i,k)
          endif
        endif
  200   continue
  300 continue
!
      return
      end
