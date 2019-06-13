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
!-------------------------------------- LICENCE END --------------------------

subroutine equivmount1(gc,sxx,syy,sxy, &
     ilg,il1,il2, &
     slope,xcent,mtdir)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ilg,il1,il2
   real gc(ilg),sxx(ilg),syy(ilg),sxy(ilg)
   real slope(ilg),xcent(ilg),mtdir(ilg)
   !@Author
   !         A. Zadra - Feb 2002
   !                    based on formulae by Lott & Miller (1997)
   !@Object compute elliptical mountain parameters
   !         To calculate the slope, eccentricity and direction
   !         of the equivalent/representative elliptical mountain
   !         using the launching height and the gradient correlation
   !         tensor of the unresolved topography.
   !@Revisions
   ! 001     Add protection to calculation of sppr, to avoid negative
   !         values that might occur in case the geophysical fields
   !         are interpolated
   !@Arguments
   !         - Input -
   !GC       land-sea mask (between 0(sea) and -1(land))
   !SXX      gradient correlation, x-direction
   !SYY      gradient correlation, y-direction
   !SXY      gradient cross-correlation, xy-directions
   !ILG      total horizontal dimension
   !IL1      index of 1st point to start calculation
   !IL2      index of last point to stop calculation
   !         - Output -
   !SLOPE    mountain slope
   !XCENT    mountain eccentricity
   !MTDIR    mountain direction

   integer i
   real(REAL64) :: zero,unit,piotwo,hmin
   real(REAL64) :: spls,smns,scrs,sppr,smpr

   zero = 0.0
   unit = 1.0
   piotwo = .5*acos(-1.)
   hmin = 3.0

   !     Initialize mountain parameters:
   do i=1,ilg
      slope(i) = zero
      xcent(i) = unit
      mtdir(i) = zero
   enddo

   !     Mountain slope, eccentricity and direction:
   do i=il1,il2
      if (gc(i).eq.-1.) then
         spls = .5*( sxx(i) + syy(i) )
         smns = .5*( sxx(i) - syy(i) )
         scrs = sxy(i)
         sppr =     spls + sqrt(smns**2 + scrs**2)
         smpr = abs(spls - sqrt(smns**2 + scrs**2))

         sppr = max(zero,sppr)
         slope(i) = sqrt(sppr)

         if (sppr .lt. 1.e-10) then
            xcent(i) = unit
         else
            xcent(i) = sqrt(smpr/sppr)
         endif

         if ( abs(scrs) .lt. 1.e-10 .and. abs(smns) .lt. 1.e-10) then
            mtdir(i) = zero
         else
            mtdir(i) = .5*atan2(scrs,smns)
            if (mtdir(i).gt.  piotwo ) mtdir(i) = mtdir(i) - 2.*piotwo
            if (mtdir(i).lt.(-piotwo)) mtdir(i) = mtdir(i) + 2.*piotwo
         endif
      endif
   enddo

   return
end subroutine equivmount1
