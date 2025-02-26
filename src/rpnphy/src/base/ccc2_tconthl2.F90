!**S/P TCONTHL - WATER VAPOR CONTINUUM

module ccc2_tconthl
   implicit none
   private
   public :: ccc2_tconthl2
   
contains

      subroutine ccc2_tconthl2(taug, coef1, coef2, s1, dp, dip, dir, dt, &
                          inptr, inpt, mcont, ni, lay)
      implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: ni, lay
      integer, intent(in) :: mcont(ni)
      real, intent(inout) :: taug(ni,lay)
      real, intent(in) :: coef1(5,5,4), coef2(5,5,4)
      real, intent(in) :: s1(ni,lay), dp(ni,lay), dip(ni,lay), dir(ni,lay), &
           dt(ni,lay)
      integer, intent(in) :: inptr(ni,lay), inpt(ni,lay)

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
! 001   M. Lazare (Jun 19/2006) - previous version tconthl1 for gcm15e:
!                               - implement rpn fix for inpt.
!                               - use variable instead of constant
!                                 in intrinsics such as "max",
!                                 so that can compile in 32-bit mode
!                                 with real  .
! 002   J. Li (Jun 21/2006)     - previous version tconthl2 for gcm15f:
!                               - new scheme to properly account for the
!                                 h2o/c)2 ratio.
! 003  L. Solheim (Apr 21/2008) - previous version tconthl3 for gcm15g:
!                               - cosmetic change to add threadprivate
!                                 for common block "trace", in support
!                                 of "radforce" model option.
! 004  J. Li(Feb 9,2009)       -  new version for gcm15h:
!                              - remove unused trace common block.
! 005  P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the Li/Lazare 2009 version. (tconthl4)

!Object
!
!        Water vapor continuum for 540-800 cm-1. different from tcontl,
!        variation of mass mixing ratio for h2o and co2 is consider.
!        lc = 4, but the inptut data are with 6 group, 2 of them are used
!        for mass mixing ratio changes
!
!Arguments
!
! taug   gaseous optical depth
! s1     input h2o mixing ratio for each layer
! dp     air mass path for a model layer (explained in raddriv)
! dip    interpolation factor for pressure between two
!        neighboring standard input data pressure levels
! dir    interpolation factor for mass ratio of h2o / co2 between
!        two neighboring standard input ratios
! dt     layer temperature - 250 k
! inptr  number of the ratio level for the standard 5 ratios
! inpt   number of the level for the standard input data pressures
! mcont  the highest level for water vapor continuum calculation
!
!Implicites
#include "ccc_tracegases.cdk"
!
!*
      integer  k, i, j, m, n, l, lp1
      real x1, y1, x2, y2, x11, x21, y21, y11, x12, x22, y12, y22

      DO_I: do i = 1, ni
      do k = mcont(i), lay
       if (inpt(1,k) .lt. 950)                                      then
         m =  inpt(i,k) - 14
         n =  m + 1
         if (m .ge. 1)                                              then
          l   =  inptr(i,k)
          if (l .lt. 1)                                             then
            x1  =  coef1(1,1,m) + dt(i,k) * (coef1(2,1,m) + dt(i,k) * &
                  (coef1(3,1,m) + dt(i,k) * (coef1(4,1,m) + &
                   dt(i,k) * coef1(5,1,m))))
            x2  =  coef1(1,1,n) + dt(i,k) * (coef1(2,1,n) + dt(i,k) * &
                  (coef1(3,1,n) + dt(i,k) * (coef1(4,1,n) + &
                   dt(i,k) * coef1(5,1,n))))
!
            y1  =  coef2(1,1,m) + dt(i,k) * (coef2(2,1,m) + dt(i,k) * &
                  (coef2(3,1,m) + dt(i,k) * (coef2(4,1,m) + &
                   dt(i,k) * coef2(5,1,m))))
            y2  =  coef2(1,1,n) + dt(i,k) * (coef2(2,1,n) + dt(i,k) * &
                  (coef2(3,1,n) + dt(i,k) * (coef2(4,1,n) + &
                   dt(i,k) * coef2(5,1,n))))
!
          else if (l .lt. 5)                                        then
            lp1 =  l + 1
            x11 =  coef1(1,l,m) + dt(i,k) * (coef1(2,l,m) + dt(i,k) * &
                  (coef1(3,l,m) + dt(i,k) * (coef1(4,l,m) + &
                   dt(i,k) * coef1(5,l,m))))
            x21 =  coef1(1,l,n) + dt(i,k) * (coef1(2,l,n) + dt(i,k) * &
                  (coef1(3,l,n) + dt(i,k) * (coef1(4,l,n) + &
                   dt(i,k) * coef1(5,l,n))))
!
            y11 =  coef2(1,l,m) + dt(i,k) * (coef2(2,l,m) + dt(i,k) * &
                  (coef2(3,l,m) + dt(i,k) * (coef2(4,l,m) + &
                   dt(i,k) * coef2(5,l,m))))
            y21 =  coef2(1,l,n) + dt(i,k) * (coef2(2,l,n) + dt(i,k) * &
                  (coef2(3,l,n) + dt(i,k) * (coef2(4,l,n) + &
                   dt(i,k) * coef2(5,l,n))))
!
            x12 =  coef1(1,lp1,m) + dt(i,k) * (coef1(2,lp1,m) + &
                                    dt(i,k) * (coef1(3,lp1,m) + &
                                    dt(i,k) * (coef1(4,lp1,m) + &
                                    dt(i,k) * coef1(5,lp1,m))))
            x22 =  coef1(1,lp1,n) + dt(i,k) * (coef1(2,lp1,n) + &
                                    dt(i,k) * (coef1(3,lp1,n) + &
                                    dt(i,k) * (coef1(4,lp1,n) + &
                                    dt(i,k) * coef1(5,lp1,n))))
!
            y12 =  coef2(1,lp1,m) + dt(i,k) * (coef2(2,lp1,m) + &
                                    dt(i,k) * (coef2(3,lp1,m) + &
                                    dt(i,k) * (coef2(4,lp1,m) + &
                                    dt(i,k) * coef2(5,lp1,m))))
            y22 =  coef2(1,lp1,n) + dt(i,k) * (coef2(2,lp1,n) + &
                                    dt(i,k) * (coef2(3,lp1,n) + &
                                    dt(i,k) * (coef2(4,lp1,n) + &
                                    dt(i,k) * coef2(5,lp1,n))))
!
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1(1,5,m) + dt(i,k) * (coef1(2,5,m) + dt(i,k) * &
                  (coef1(3,5,m) + dt(i,k) * (coef1(4,5,m) + &
                   dt(i,k) * coef1(5,5,m))))
            x2  =  coef1(1,5,n) + dt(i,k) * (coef1(2,5,n) + dt(i,k) * &
                  (coef1(3,5,n) + dt(i,k) * (coef1(4,5,n) + &
                   dt(i,k) * coef1(5,5,n))))
            y1  =  coef2(1,5,m) + dt(i,k) * (coef2(2,5,m) + dt(i,k) * &
                  (coef2(3,5,m) + dt(i,k) * (coef2(4,5,m) + &
                   dt(i,k) * coef2(5,5,m))))
            y2  =  coef2(1,5,n) + dt(i,k) * (coef2(2,5,n) + dt(i,k) * &
                  (coef2(3,5,n) + dt(i,k) * (coef2(4,5,n) + &
                   dt(i,k) * coef2(5,5,n))))
          endif
!
          taug(i,k) =  taug(i,k) + ((x1 - y1 + (x2 - x1 - y2 + y1) * &
                       dip(i,k)) * 1.608 * s1(i,k) + y1 + (y2 - y1) * &
                       dip(i,k)) * s1(i,k) * dp(i,k)
         endif
       else
        j =  inpt(1,k) - 1000
        m =  j - 14
        n =  m + 1
         if (m .ge. 1)                                              then
          l   =  inptr(i,k)
          if (l .lt. 1)                                             then
            x1  =  coef1(1,1,m) + dt(i,k) * (coef1(2,1,m) + dt(i,k) * &
                  (coef1(3,1,m) + dt(i,k) * (coef1(4,1,m) + &
                   dt(i,k) * coef1(5,1,m))))
            x2  =  coef1(1,1,n) + dt(i,k) * (coef1(2,1,n) + dt(i,k) * &
                  (coef1(3,1,n) + dt(i,k) * (coef1(4,1,n) + &
                   dt(i,k) * coef1(5,1,n))))
!
            y1  =  coef2(1,1,m) + dt(i,k) * (coef2(2,1,m) + dt(i,k) * &
                  (coef2(3,1,m) + dt(i,k) * (coef2(4,1,m) + &
                   dt(i,k) * coef2(5,1,m))))
            y2  =  coef2(1,1,n) + dt(i,k) * (coef2(2,1,n) + dt(i,k) * &
                  (coef2(3,1,n) + dt(i,k) * (coef2(4,1,n) + &
                   dt(i,k) * coef2(5,1,n))))
!
          else if (l .lt. 5)                                        then
            lp1 =  l + 1
            x11 =  coef1(1,l,m) + dt(i,k) * (coef1(2,l,m) + dt(i,k) * &
                  (coef1(3,l,m) + dt(i,k) * (coef1(4,l,m) + &
                   dt(i,k) * coef1(5,l,m))))
            x21 =  coef1(1,l,n) + dt(i,k) * (coef1(2,l,n) + dt(i,k) * &
                  (coef1(3,l,n) + dt(i,k) * (coef1(4,l,n) + &
                   dt(i,k) * coef1(5,l,n))))
!
            y11 =  coef2(1,l,m) + dt(i,k) * (coef2(2,l,m) + dt(i,k) * &
                  (coef2(3,l,m) + dt(i,k) * (coef2(4,l,m) + &
                   dt(i,k) * coef2(5,l,m))))
            y21 =  coef2(1,l,n) + dt(i,k) * (coef2(2,l,n) + dt(i,k) * &
                  (coef2(3,l,n) + dt(i,k) * (coef2(4,l,n) + &
                   dt(i,k) * coef2(5,l,n))))
!
            x12 =  coef1(1,lp1,m) + dt(i,k) * (coef1(2,lp1,m) + &
                                    dt(i,k) * (coef1(3,lp1,m) + &
                                    dt(i,k) * (coef1(4,lp1,m) + &
                                    dt(i,k) * coef1(5,lp1,m))))
            x22 =  coef1(1,lp1,n) + dt(i,k) * (coef1(2,lp1,n) + &
                                    dt(i,k) * (coef1(3,lp1,n) + &
                                    dt(i,k) * (coef1(4,lp1,n) + &
                                    dt(i,k) * coef1(5,lp1,n))))
!
            y12 =  coef2(1,lp1,m) + dt(i,k) * (coef2(2,lp1,m) + &
                                    dt(i,k) * (coef2(3,lp1,m) + &
                                    dt(i,k) * (coef2(4,lp1,m) + &
                                    dt(i,k) * coef2(5,lp1,m))))
            y22 =  coef2(1,lp1,n) + dt(i,k) * (coef2(2,lp1,n) + &
                                    dt(i,k) * (coef2(3,lp1,n) + &
                                    dt(i,k) * (coef2(4,lp1,n) + &
                                    dt(i,k) * coef2(5,lp1,n))))
!
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1(1,5,m) + dt(i,k) * (coef1(2,5,m) + dt(i,k) * &
                  (coef1(3,5,m) + dt(i,k) * (coef1(4,5,m) + &
                   dt(i,k) * coef1(5,5,m))))
            x2  =  coef1(1,5,n) + dt(i,k) * (coef1(2,5,n) + dt(i,k) * &
                  (coef1(3,5,n) + dt(i,k) * (coef1(4,5,n) + &
                   dt(i,k) * coef1(5,5,n))))
            y1  =  coef2(1,5,m) + dt(i,k) * (coef2(2,5,m) + dt(i,k) * &
                  (coef2(3,5,m) + dt(i,k) * (coef2(4,5,m) + &
                   dt(i,k) * coef2(5,5,m))))
            y2  =  coef2(1,5,n) + dt(i,k) * (coef2(2,5,n) + dt(i,k) * &
                  (coef2(3,5,n) + dt(i,k) * (coef2(4,5,n) + &
                   dt(i,k) * coef2(5,5,n))))
          endif
!
          taug(i,k) =  taug(i,k) + ((x1 - y1 + (x2 - x1 - y2 + y1) * &
                       dip(i,k)) * 1.608 * s1(i,k) + y1 + (y2 - y1) * &
                       dip(i,k)) * s1(i,k) * dp(i,k)
         endif
        endif
      enddo
      enddo DO_I
      return
   end subroutine ccc2_tconthl2

end module ccc2_tconthl
