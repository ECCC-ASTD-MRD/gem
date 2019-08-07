!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

subroutine ccc2_tline1z(taug, coef1, s1,  dp, dip, &
     dt, inpt, lev1, gh, lc, iplus, &
     il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, lc, lev1, iplus, il1, il2
   real taug(ilg,lay), coef1(5,lc), s1(ilg,lay), dp(ilg,lay), &
        dip(ilg,lay), dt(ilg,lay)
   integer inpt(ilg,lay)
   logical gh

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001    M. Lazare   (May 05,2006)  -   previous version tline1x for gcm15e/f:
   !                                   - cosmetic cleanup to eliminate duplicate
   !                                        coding between iplus=1 and  iplus=2, by
   !                                        initializing taug if iplus=2 (this was
   !                                        the only difference between the two sections).
   !                                     - implement rpn fix for inpt.
   ! 002   M. Lazare (Apr 18,2008)       -  previous version tline1y for gcm15g:
   !                               - cosmetic change of n->lc passed
   !                                 in and used in dimension of coef
   !                                 array(s), so that it is not redefined
   !                                 and passed back out to calling
   !                                 routine, changing value in data
   !                                 statement!
   !                               - cosmetic change to add threadprivate
   !                                 for common block "trace", in support
   !                                 of "radforce" model option.
   !                               - work array "s1" now local instead of
   !                                 passed-in workspace.
   ! 003    M.Lazarre,K.Winger   (Apr 08) - correct bug - lc now used for dimension of array and n as a variable integer
   !
   ! 004    J.Li     (Feb 09,2009) -  new version for gcm15h:
   !                               - 3d ghg implemented, thus no need
   !                                 for "trace" common block or
   !                                 temporary work arrays to hold
   !                                 mixing ratios of ghg depending on
   !                                 a passed, specified option.
   ! 005  P. vaillancourt, A-M. Leduc (Apr 2010) - Integrating the Li/Lazare 2009 version.(tline1z)
   !                                  remove arguments ng , s1

   !@Object
   !        Calculation of optical depth for one gas (line contribution only)
   !        the gaseous absorption coefficients in units of cm^2 / gram.
   !        s in mass mixing ratio. absorption coefficient are calculated
   !        at the temperature t for the 18 or 26 (28) pressure levels.
   !        isl = 1 for solar, isl = 2 for infrared.

   !@Arguments
   ! taug  gaseous optical depth
   ! s1    input gas mixing ratio for each layer
   ! dp    air mass path for a model layer (explained in raddriv)
   ! dip   interpolation factor for pressure between two neighboring
   !       standard input data pressure levels
   ! dt    layer temperature - 250 k
   ! inpt  number of the level for the standard input data pressures

#include "ccc_tracegases.cdk"

   integer k, i, m, n
   integer lay1, lay2
   real x1, x2

   if (gh) then
      lay1 =  1
   else
      lay1 =  lev1
   endif
   lay2   =  lay

   !     initialize taug if iplus=2.

   if (iplus .eq. 2) then
      do k = lay1, lay2
         do i = il1, il2
            taug(i,k)     =  0.
         enddo
      enddo
   endif

   do k = lay1, lay2
      if (inpt(1,k) .lt. 950) then
         do i = il1, il2
            m  =  inpt(i,k)
            n  =  m + 1
            x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                 dt(i,k) * (coef1(3,n) + dt(i,k) * &
                 (coef1(4,n) + dt(i,k) * coef1(5,n))))
            if (m .gt. 0) then
               x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                    dt(i,k) * (coef1(3,m) + dt(i,k) * &
                    (coef1(4,m) + dt(i,k) * coef1(5,m))))
            else
               x1      =  0.0
            endif

            taug(i,k) =  taug(i,k) + (x1 + (x2 - x1) * dip(i,k)) * &
                 s1(i,k) * dp(i,k)
         enddo
      else
         m  =  inpt(1,k) - 1000
         n  =  m + 1
         do i = il1, il2
            x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                 dt(i,k) * (coef1(3,n) + dt(i,k) * &
                 (coef1(4,n) + dt(i,k) * coef1(5,n))))
            if (m .gt. 0) then
               x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                    dt(i,k) * (coef1(3,m) + dt(i,k) * &
                    (coef1(4,m) + dt(i,k) * coef1(5,m))))
            else
               x1      =  0.0
            endif

            taug(i,k) =  taug(i,k) + (x1 + (x2 - x1) * dip(i,k)) * &
                 s1(i,k) * dp(i,k)
         enddo
      endif
   enddo

   return
end subroutine ccc2_tline1z
