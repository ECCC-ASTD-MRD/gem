!**S/P TLINEHC - OPTICAL DEPTH FOR H2O AND CO2

module ccc2_tlinehc_m
   implicit none
   private
   public :: ccc2_tlinehc2
   
contains
   
      subroutine ccc2_tlinehc2 (taug, coef1u, coef1d, &
                          coef2u, coef2d, qq, co2, dp, dip, dir, &
                          dt, inptr, inpt, lev1, ni, lay)
      implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: ni, lay, lev1
      real, intent(inout) :: taug(ni,lay)
      real, intent(in) :: coef1u(5,11), coef1d(5,5,7), &
                          coef2u(5,11), coef2d(5,5,7)
      real, intent(in) :: qq(ni,lay), co2(ni,lay), dp(ni,lay), dip(ni,lay), &
           dir(ni,lay), dt(ni,lay)
      integer, intent(in) :: inptr(ni,lay), inpt(ni,lay)

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
! 001  M. Lazare  (May 05,2006) - previoius version tlinehc2 for gcm15e:
!                         - implement rpn fix for inpt.
! 002  M. Lazare  (Apr 18,2008)-  previous version tlinehc3 for gcm15g:
!                               - cosmetic change to add threadprivate
!                                 for common block "trace", in support
!                                 of "radforce" model option.
! 003  J. Li  (Feb 09,2009) -   new version for gcm15h:
!                               - 3d ghg implemented, thus no need
!                                 for "trace" common block or
!                                 temporary work arrays to hold
!                                 mixing ratios of ghg depending on
!                                 a passed, specified option.
! 004  P. vaillancourt, A-M. Leduc (Apr 2010) - Integrating the Li/Lazare 2009 version.(tlinehc4)
!                                  add co2
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
      integer  k, i, m, n, l, lp1, r
      real x1, y1, x2, y2, x11, x21, y21, y11, x12, x22, y12, y22


!*
      DO_K: do k = lev1, lay

      IF_INPT950: if (inpt(1,k) .lt. 950)                                       then

      DO_I1: do i = 1, ni
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
                      (y1 + (y2 - y1) * dip(i,k)) * co2(i,k) ) * dp(i,k)

      enddo DO_I1
     
      else  !IF_INPT950
      
      m       =  inpt(1,k) - 1000
      DO_I2: do i = 1, ni
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
                      (y1 + (y2 - y1) * dip(i,k)) * co2(i,k) ) * dp(i,k)
      enddo DO_I2
      endif IF_INPT950
      enddo DO_K

      return
   end subroutine ccc2_tlinehc2

end module ccc2_tlinehc_m
