!**S/P TCONTL - INFRARED WATER VAPOR CONTINUUM

module ccc2_tcontl
   implicit none
   private
   public :: ccc2_tcontl2
   
contains

      subroutine ccc2_tcontl2(taug, coef1, coef2, s1, dp, dip, dt, &
                         lc, inpt, mcont, gh, nig, ni, lay)
      implicit none
!!!#include <arch_specific.hf>
!
      integer, intent(in) :: ni, lay, lc, nig
      integer, intent(in) :: mcont(ni)
      real, intent(inout) :: taug(ni,lay)
      real, intent(in) :: coef1(5,lc), coef2(5,lc)
      real, intent(in) :: s1(ni,lay), dp(ni,lay), dip(ni,lay), dt(ni,lay)
      integer, intent(in) :: inpt(ni,lay)
      logical, intent(in) :: gh

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
      integer  :: k, i, j, m, n, nc, w1, j1
      real :: x1, y1, x2, y2

      if (gh) then
        nc =  29 - lc
      else
        nc =  19 - lc
      endif

      !#opt V6b - a bit faster than 5b... more ugliness
      do k = minval(mcont(1:nig)), lay
         IF950: if (inpt(1,k) < 950) then
            
            do i = 1, nig
               j = inpt(i,k)
               if (j >= nc .and. k >= mcont(i)) then
                  m  =  j - nc + 1
                  n  =  m + 1
                  x1 =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                       dt(i,k) * (coef1(3,m) + dt(i,k) * &
                       (coef1(4,m) + dt(i,k) * coef1(5,m))))
                  x2 =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                       dt(i,k) * (coef1(3,n) + dt(i,k) * &
                       (coef1(4,n) + dt(i,k) * coef1(5,n))))
                  y1 =  coef2(1,m) + dt(i,k) * (coef2(2,m) + &
                       dt(i,k) * (coef2(3,m) + dt(i,k) * &
                       (coef2(4,m) + dt(i,k) * coef2(5,m))))
                  y2 =  coef2(1,n) + dt(i,k) * (coef2(2,n) + &
                       dt(i,k) * (coef2(3,n) + dt(i,k) * &
                       (coef2(4,n) + dt(i,k) * coef2(5,n))))
                  taug(i,k) =  taug(i,k) + &
                       ( (x1 - y1 + (x2 - x1 - y2 + y1) * &
                       dip(i,k)) * 1.608 * s1(i,k) + &
                       y1 + (y2 - y1) * dip(i,k) ) * &
                       s1(i,k) * dp(i,k)
               endif
            enddo
            

         else  !IF950
            
            j =  inpt(1,k) - 1000
            if (j >= nc)  then
               m  =  j - nc + 1
               n  =  m + 1
               do i = 1, nig
                  if (k >= mcont(i)) then
                     x1 =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                          dt(i,k) * (coef1(3,m) + dt(i,k) * &
                          (coef1(4,m) + dt(i,k) * coef1(5,m))))
                     x2 =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                          dt(i,k) * (coef1(3,n) + dt(i,k) * &
                          (coef1(4,n) + dt(i,k) * coef1(5,n))))
                     y1 =  coef2(1,m) + dt(i,k) * (coef2(2,m) + &
                          dt(i,k) * (coef2(3,m) + dt(i,k) * &
                          (coef2(4,m) + dt(i,k) * coef2(5,m))))
                     y2 =  coef2(1,n) + dt(i,k) * (coef2(2,n) + &
                          dt(i,k) * (coef2(3,n) + dt(i,k) * &
                          (coef2(4,n) + dt(i,k) * coef2(5,n))))
                     taug(i,k) =  taug(i,k) + &
                          ( (x1 - y1 + (x2 - x1 - y2 + y1) * &
                          dip(i,k)) * 1.608 * s1(i,k) + &
                          y1 + (y2 - y1) * dip(i,k) ) * &
                          s1(i,k) * dp(i,k)
                  endif
               enddo
            endif
            
         endif IF950
      enddo
      
   return
   end subroutine ccc2_tcontl2
   
end module ccc2_tcontl
