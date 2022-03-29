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

subroutine ccc1_tline2(taug, coef1, coef2, qq, o3, &
     ng2, dp, dip, dt, inpt, &
     lev1, gh, lc, il1, il2, ilg, lay, s2)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, lc, lev1, il1, il2, lay1, lay2
   integer k, i, m, n, ng2
   real x1, y1, x2, y2
   real taug(ilg,lay), coef1(5,lc), coef2(5,lc), qq(ilg,lay), &
        o3(ilg,lay), dp(ilg,lay), dip(ilg,lay), dt(ilg,lay), &
        s2(ilg,lay)
   integer inpt(ilg,lay)
   logical gh

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !@Revisions
   ! 001    M.Lazarre,K.Winger   (Apr 08) - correct bug - lc now used for dimension of array and n as a variable integer

   !@Object The same as tline1, but with two mixed gases, one must be h2o
   ! taug   gaseous optical depth
   ! qq     input h2o mixing ratio for each layer
   ! o3     input another gas mixing ratio for each layer
   ! dp     air mass path for a model layer (explained in raddriv)
   ! dip    interpolation factor for pressure between two neighboring
   !        standard input data pressure levels
   ! dt     layer temperature - 250 k
   ! inpt   number of the level for the standard input data pressures

#include "ccc_tracegases.cdk"

   if (gh) then
      lay1 =  1
   else
      lay1 =  lev1
   endif
   lay2   =  lay

   if (ng2 .eq. 2) then

      do k = lay1, lay2
         do i = il1, il2
            s2(i,k)     =  o3(i,k)
         enddo
      enddo

   else if (ng2 .eq. 3) then

      do k = lay1, lay2
         do i = il1, il2
            s2(i,k)     =  rmco2
         enddo
      enddo

   else if (ng2 .eq. 5) then

      do k = lay1, lay2
         do i = il1, il2
            s2(i,k)     =  rmn2o
         enddo
      enddo

   else if (ng2 .eq. 6) then

      do k = lay1, lay2
         do i = il1, il2
            s2(i,k)     =  rmo2
         enddo
      enddo

   endif

   DO2000: do k = lay1, lay2

      if (inpt(1,k) .lt. 950) then

         do i = il1, il2
            m  =  inpt(i,k)
            n  =  m + 1
            x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + dt(i,k) * &
                 (coef1(3,n) + dt(i,k) * (coef1(4,n) + &
                 dt(i,k) * coef1(5,n))))

            y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + dt(i,k) * &
                 (coef2(3,n) + dt(i,k) * (coef2(4,n) + &
                 dt(i,k) * coef2(5,n))))
            if (m .gt. 0) then
               x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + dt(i,k) * &
                    (coef1(3,m) + dt(i,k) * (coef1(4,m) + &
                    dt(i,k) * coef1(5,m))))

               y1      =  coef2(1,m) + dt(i,k) * (coef2(2,m) + dt(i,k) * &
                    (coef2(3,m) + dt(i,k) * (coef2(4,m) + &
                    dt(i,k) * coef2(5,m))))
            else
               x1      =  0.0
               y1      =  0.0
            endif

            taug(i,k) = ( (x1 + (x2 - x1) * dip(i,k)) * qq(i,k) + &
                 (y1 + (y2 - y1) * dip(i,k)) * s2(i,k) ) * &
                 dp(i,k)
         enddo

      else

         m  =  inpt(1,k) - 1000
         n  =  m + 1
         do i = il1, il2
            x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + dt(i,k) * &
                 (coef1(3,n) + dt(i,k) * (coef1(4,n) + &
                 dt(i,k) * coef1(5,n))))

            y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + dt(i,k) * &
                 (coef2(3,n) + dt(i,k) * (coef2(4,n) + &
                 dt(i,k) * coef2(5,n))))
            if (m .gt. 0) then
               x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + dt(i,k) * &
                    (coef1(3,m) + dt(i,k) * (coef1(4,m) + &
                    dt(i,k) * coef1(5,m))))

               y1      =  coef2(1,m) + dt(i,k) * (coef2(2,m) + dt(i,k) * &
                    (coef2(3,m) + dt(i,k) * (coef2(4,m) + &
                    dt(i,k) * coef2(5,m))))
            else
               x1      =  0.0
               y1      =  0.0
            endif

            taug(i,k) = ( (x1 + (x2 - x1) * dip(i,k)) * qq(i,k) + &
                 (y1 + (y2 - y1) * dip(i,k)) * s2(i,k) ) * &
                 dp(i,k)
         enddo

      endif

   enddo DO2000

   return
end subroutine ccc1_tline2
