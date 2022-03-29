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
!**S/P TLINEHC - OPTICAL DEPTH FOR H2O AND CO2
!
      subroutine ccc1_tlinehc (taug, coef1u, coef1d, &
                          coef2u, coef2d, qq, dp, dip, dir, &
                          dt, inptr, inpt, lev1, il1, il2, ilg, lay)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev1, il1, il2, k, i, m, n, l, lp1, r
      real x1, y1, x2, y2, x11, x21, y21, y11, x12, x22, y12, y22
      real taug(ilg,lay), coef1u(5,11), coef1d(5,5,7), &
                          coef2u(5,11), coef2d(5,5,7)
      real qq(ilg,lay), dp(ilg,lay), dip(ilg,lay), dir(ilg,lay), &
           dt(ilg,lay)
      integer inptr(ilg,lay), inpt(ilg,lay)
!
!Authors
!
!        J. Li, M. Lazarre, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!         JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001
!
!Object
!
!        This subroutine determines the optical depth for h2o and co2 in
!        the region of 540-800 cm^-1
!
!Arguments
!
! taug   gaseous optical depth
! qq     input h2o mixing ratio for each layer
! dp     air mass path for a model layer (explained in raddriv)
! dip    interpolation factor for pressure between two neighboring
!        standard input data pressure levels
! dir    interpolation factor for mass ratio of h2o / co2 between
!        two neighboring standard input ratios
! dt     layer temperature - 250 k
! inpr   number of the ratio level for the standard 5 ratios
! inpt   number of the level for the standard input data pressures
!
!Implicites
!
#include "ccc_tracegases.cdk"
!
!*
      do 200 k = lev1, lay
      if (inpt(1,k) .lt. 950)                                       then
      do 101 i = il1, il2
        m       =  inpt(i,k)
        if (m .le. 11)                                              then
          n     =  m + 1
          if (m .gt. 0)                                             then
            x1  =  coef1u(1,m) + dt(i,k) * (coef1u(2,m) + dt(i,k) * &
                  (coef1u(3,m) + dt(i,k) * (coef1u(4,m) + &
                   dt(i,k) * coef1u(5,m))))
            y1  =  coef2u(1,m) + dt(i,k) * (coef2u(2,m) + dt(i,k) * &
                  (coef2u(3,m) + dt(i,k) * (coef2u(4,m) + &
                   dt(i,k) * coef2u(5,m))))
          else
            x1  =  0.0
            y1  =  0.0
          endif
          if (m .lt. 11)                                            then
            x2  =  coef1u(1,n) + dt(i,k) * (coef1u(2,n) + dt(i,k) * &
                  (coef1u(3,n) + dt(i,k) * (coef1u(4,n) + &
                   dt(i,k) * coef1u(5,n))))
            y2  =  coef2u(1,n) + dt(i,k) * (coef2u(2,n) + dt(i,k) * &
                  (coef2u(3,n) + dt(i,k) * (coef2u(4,n) + &
                   dt(i,k) * coef2u(5,n))))
          else
            x2  =  coef1d(1,1,1) + dt(i,k) * (coef1d(2,1,1) + dt(i,k) * &
                  (coef1d(3,1,1) + dt(i,k) * (coef1d(4,1,1) + &
                   dt(i,k) * coef1d(5,1,1))))
            y2  =  coef2d(1,1,1) + dt(i,k) * (coef2d(2,1,1) + dt(i,k) * &
                  (coef2d(3,1,1) + dt(i,k) * (coef2d(4,1,1) + &
                   dt(i,k) * coef2d(5,1,1))))
          endif
        else
          m     =  m - 11
          n     =  m + 1
          l     =  inptr(i,k)
          if (l .lt. 1)                                             then
            x1  =  coef1d(1,1,m) + dt(i,k) * (coef1d(2,1,m) + dt(i,k) * &
                  (coef1d(3,1,m) + dt(i,k) * (coef1d(4,1,m) + &
                   dt(i,k) * coef1d(5,1,m))))
            x2  =  coef1d(1,1,n) + dt(i,k) * (coef1d(2,1,n) + dt(i,k) * &
                  (coef1d(3,1,n) + dt(i,k) * (coef1d(4,1,n) + &
                   dt(i,k) * coef1d(5,1,n))))
!
            y1  =  coef2d(1,1,m) + dt(i,k) * (coef2d(2,1,m) + dt(i,k) * &
                  (coef2d(3,1,m) + dt(i,k) * (coef2d(4,1,m) + &
                   dt(i,k) * coef2d(5,1,m))))
            y2  =  coef2d(1,1,n) + dt(i,k) * (coef2d(2,1,n) + dt(i,k) * &
                  (coef2d(3,1,n) + dt(i,k) * (coef2d(4,1,n) + &
                   dt(i,k) * coef2d(5,1,n))))
!
          else if (l .lt. 5)                                        then
            lp1 =  l + 1
            x11 =  coef1d(1,l,m) + dt(i,k) * (coef1d(2,l,m) + dt(i,k) * &
                  (coef1d(3,l,m) + dt(i,k) * (coef1d(4,l,m) + &
                   dt(i,k) * coef1d(5,l,m))))
            x21 =  coef1d(1,l,n) + dt(i,k) * (coef1d(2,l,n) + dt(i,k) * &
                  (coef1d(3,l,n) + dt(i,k) * (coef1d(4,l,n) + &
                   dt(i,k) * coef1d(5,l,n))))
!
            y11 =  coef2d(1,l,m) + dt(i,k) * (coef2d(2,l,m) + dt(i,k) * &
                  (coef2d(3,l,m) + dt(i,k) * (coef2d(4,l,m) + &
                   dt(i,k) * coef2d(5,l,m))))
            y21 =  coef2d(1,l,n) + dt(i,k) * (coef2d(2,l,n) + dt(i,k) * &
                  (coef2d(3,l,n) + dt(i,k) * (coef2d(4,l,n) + &
                   dt(i,k) * coef2d(5,l,n))))
!
            x12 =  coef1d(1,lp1,m) + dt(i,k) * (coef1d(2,lp1,m) + &
                                     dt(i,k) * (coef1d(3,lp1,m) + &
                                     dt(i,k) * (coef1d(4,lp1,m) + &
                                     dt(i,k) * coef1d(5,lp1,m))))
            x22 =  coef1d(1,lp1,n) + dt(i,k) * (coef1d(2,lp1,n) + &
                                     dt(i,k) * (coef1d(3,lp1,n) + &
                                     dt(i,k) * (coef1d(4,lp1,n) + &
                                     dt(i,k) * coef1d(5,lp1,n))))
!
            y12 =  coef2d(1,lp1,m) + dt(i,k) * (coef2d(2,lp1,m) + &
                                     dt(i,k) * (coef2d(3,lp1,m) + &
                                     dt(i,k) * (coef2d(4,lp1,m) + &
                                     dt(i,k) * coef2d(5,lp1,m))))
            y22 =  coef2d(1,lp1,n) + dt(i,k) * (coef2d(2,lp1,n) + &
                                     dt(i,k) * (coef2d(3,lp1,n) + &
                                     dt(i,k) * (coef2d(4,lp1,n) + &
                                     dt(i,k) * coef2d(5,lp1,n))))
!
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1d(1,5,m) + dt(i,k) * (coef1d(2,5,m) + dt(i,k) * &
                  (coef1d(3,5,m) + dt(i,k) * (coef1d(4,5,m) + &
                   dt(i,k) * coef1d(5,5,m))))
            x2  =  coef1d(1,5,n) + dt(i,k) * (coef1d(2,5,n) + dt(i,k) * &
                  (coef1d(3,5,n) + dt(i,k) * (coef1d(4,5,n) + &
                   dt(i,k) * coef1d(5,5,n))))
            y1  =  coef2d(1,5,m) + dt(i,k) * (coef2d(2,5,m) + dt(i,k) * &
                  (coef2d(3,5,m) + dt(i,k) * (coef2d(4,5,m) + &
                   dt(i,k) * coef2d(5,5,m))))
            y2  =  coef2d(1,5,n) + dt(i,k) * (coef2d(2,5,n) + dt(i,k) * &
                  (coef2d(3,5,n) + dt(i,k) * (coef2d(4,5,n) + &
                   dt(i,k) * coef2d(5,5,n))))
          endif
        endif
!
        taug(i,k) = ( (x1 + (x2 - x1) * dip(i,k)) * qq(i,k) + &
                      (y1 + (y2 - y1) * dip(i,k)) * rmco2 ) * dp(i,k)
 101  continue
      else
      m       =  inpt(1,k) - 1000
      do 102 i = il1, il2
        if (m .le. 11)                                              then
          n     =  m + 1
          if (m .gt. 0)                                             then
            x1  =  coef1u(1,m) + dt(i,k) * (coef1u(2,m) + dt(i,k) * &
                  (coef1u(3,m) + dt(i,k) * (coef1u(4,m) + &
                   dt(i,k) * coef1u(5,m))))
            y1  =  coef2u(1,m) + dt(i,k) * (coef2u(2,m) + dt(i,k) * &
                  (coef2u(3,m) + dt(i,k) * (coef2u(4,m) + &
                   dt(i,k) * coef2u(5,m))))
          else
            x1  =  0.0
            y1  =  0.0
          endif
          if (m .lt. 11)                                            then
            x2  =  coef1u(1,n) + dt(i,k) * (coef1u(2,n) + dt(i,k) * &
                  (coef1u(3,n) + dt(i,k) * (coef1u(4,n) + &
                   dt(i,k) * coef1u(5,n))))
            y2  =  coef2u(1,n) + dt(i,k) * (coef2u(2,n) + dt(i,k) * &
                  (coef2u(3,n) + dt(i,k) * (coef2u(4,n) + &
                   dt(i,k) * coef2u(5,n))))
          else
            x2  =  coef1d(1,1,1) + dt(i,k) * (coef1d(2,1,1) + dt(i,k) * &
                  (coef1d(3,1,1) + dt(i,k) * (coef1d(4,1,1) + &
                   dt(i,k) * coef1d(5,1,1))))
            y2  =  coef2d(1,1,1) + dt(i,k) * (coef2d(2,1,1) + dt(i,k) * &
                  (coef2d(3,1,1) + dt(i,k) * (coef2d(4,1,1) + &
                   dt(i,k) * coef2d(5,1,1))))
          endif
        else
          r     =  m - 11
          n     =  r + 1
          l     =  inptr(i,k)
          if (l .lt. 1)                                             then
            x1  =  coef1d(1,1,r) + dt(i,k) * (coef1d(2,1,r) + dt(i,k) * &
                  (coef1d(3,1,r) + dt(i,k) * (coef1d(4,1,r) + &
                   dt(i,k) * coef1d(5,1,r))))
            x2  =  coef1d(1,1,n) + dt(i,k) * (coef1d(2,1,n) + dt(i,k) * &
                  (coef1d(3,1,n) + dt(i,k) * (coef1d(4,1,n) + &
                   dt(i,k) * coef1d(5,1,n))))
!
            y1  =  coef2d(1,1,r) + dt(i,k) * (coef2d(2,1,r) + dt(i,k) * &
                  (coef2d(3,1,r) + dt(i,k) * (coef2d(4,1,r) + &
                   dt(i,k) * coef2d(5,1,r))))
            y2  =  coef2d(1,1,n) + dt(i,k) * (coef2d(2,1,n) + dt(i,k) * &
                  (coef2d(3,1,n) + dt(i,k) * (coef2d(4,1,n) + &
                   dt(i,k) * coef2d(5,1,n))))
!
          else if (l .lt. 5)                                        then
            lp1 =  l + 1
            x11 =  coef1d(1,l,r) + dt(i,k) * (coef1d(2,l,r) + dt(i,k) * &
                  (coef1d(3,l,r) + dt(i,k) * (coef1d(4,l,r) + &
                   dt(i,k) * coef1d(5,l,r))))
            x21 =  coef1d(1,l,n) + dt(i,k) * (coef1d(2,l,n) + dt(i,k) * &
                  (coef1d(3,l,n) + dt(i,k) * (coef1d(4,l,n) + &
                   dt(i,k) * coef1d(5,l,n))))
!
            y11 =  coef2d(1,l,r) + dt(i,k) * (coef2d(2,l,r) + dt(i,k) * &
                  (coef2d(3,l,r) + dt(i,k) * (coef2d(4,l,r) + &
                   dt(i,k) * coef2d(5,l,r))))
            y21 =  coef2d(1,l,n) + dt(i,k) * (coef2d(2,l,n) + dt(i,k) * &
                  (coef2d(3,l,n) + dt(i,k) * (coef2d(4,l,n) + &
                   dt(i,k) * coef2d(5,l,n))))
!
            x12 =  coef1d(1,lp1,r) + dt(i,k) * (coef1d(2,lp1,r) + &
                                     dt(i,k) * (coef1d(3,lp1,r) + &
                                     dt(i,k) * (coef1d(4,lp1,r) + &
                                     dt(i,k) * coef1d(5,lp1,r))))
            x22 =  coef1d(1,lp1,n) + dt(i,k) * (coef1d(2,lp1,n) + &
                                     dt(i,k) * (coef1d(3,lp1,n) + &
                                     dt(i,k) * (coef1d(4,lp1,n) + &
                                     dt(i,k) * coef1d(5,lp1,n))))
!
            y12 =  coef2d(1,lp1,r) + dt(i,k) * (coef2d(2,lp1,r) + &
                                     dt(i,k) * (coef2d(3,lp1,r) + &
                                     dt(i,k) * (coef2d(4,lp1,r) + &
                                     dt(i,k) * coef2d(5,lp1,r))))
            y22 =  coef2d(1,lp1,n) + dt(i,k) * (coef2d(2,lp1,n) + &
                                     dt(i,k) * (coef2d(3,lp1,n) + &
                                     dt(i,k) * (coef2d(4,lp1,n) + &
                                     dt(i,k) * coef2d(5,lp1,n))))
!
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1d(1,5,r) + dt(i,k) * (coef1d(2,5,r) + dt(i,k) * &
                  (coef1d(3,5,r) + dt(i,k) * (coef1d(4,5,r) + &
                   dt(i,k) * coef1d(5,5,r))))
            x2  =  coef1d(1,5,n) + dt(i,k) * (coef1d(2,5,n) + dt(i,k) * &
                  (coef1d(3,5,n) + dt(i,k) * (coef1d(4,5,n) + &
                   dt(i,k) * coef1d(5,5,n))))
            y1  =  coef2d(1,5,r) + dt(i,k) * (coef2d(2,5,r) + dt(i,k) * &
                  (coef2d(3,5,r) + dt(i,k) * (coef2d(4,5,r) + &
                   dt(i,k) * coef2d(5,5,r))))
            y2  =  coef2d(1,5,n) + dt(i,k) * (coef2d(2,5,n) + dt(i,k) * &
                  (coef2d(3,5,n) + dt(i,k) * (coef2d(4,5,n) + &
                   dt(i,k) * coef2d(5,5,n))))
          endif
        endif
!
        taug(i,k) = ( (x1 + (x2 - x1) * dip(i,k)) * qq(i,k) + &
                      (y1 + (y2 - y1) * dip(i,k)) * rmco2 ) * dp(i,k)
 102  continue
      endif
 200  continue
!
      return
      end
