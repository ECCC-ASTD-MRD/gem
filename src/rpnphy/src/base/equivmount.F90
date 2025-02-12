
module equivmount
   implicit none
   private
   public :: equivmount2

contains

subroutine equivmount2(gc, sxx, syy, sxy, slope, xcent, mtdir, ni)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: ni
   real, intent(in) :: gc(ni), sxx(ni), syy(ni), sxy(ni)
   real, intent(out) :: slope(ni), xcent(ni), mtdir(ni)
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
   !GC       land-sea mask (between 0(sea) and 1(land))
   !SXX      gradient correlation, x-direction
   !SYY      gradient correlation, y-direction
   !SXY      gradient cross-correlation, xy-directions
   !NI      total horizontal dimension
   !         - Output -
   !SLOPE    mountain slope
   !XCENT    mountain eccentricity
   !MTDIR    mountain direction

   real(REAL64), parameter :: zero = 0.0
   real(REAL64), parameter :: unit = 1.0
   real(REAL64), parameter :: piotwo = .5*acos(-1.)
   real(REAL64), parameter :: hmin = 3.0
   
   integer :: i
   real(REAL64) :: spls,smns,scrs,sppr,smpr

   !     Initialize mountain parameters:
   do i=1,ni
      slope(i) = zero
      xcent(i) = unit
      mtdir(i) = zero
   enddo

   !     Mountain slope, eccentricity and direction:
   do i=1,ni
      if (abs(nint(gc(i))) == 1.) then
         spls = .5*( sxx(i) + syy(i) )
         smns = .5*( sxx(i) - syy(i) )
         scrs = sxy(i)
         sppr =     spls + sqrt(smns**2 + scrs**2)
         smpr = abs(spls - sqrt(smns**2 + scrs**2))

         sppr = max(zero,sppr)
         slope(i) = sqrt(sppr)

         if (sppr < 1.e-10) then
            xcent(i) = unit
         else
            xcent(i) = sqrt(smpr/sppr)
         endif

         if (abs(scrs) < 1.e-10 .and. abs(smns) < 1.e-10) then
            mtdir(i) = zero
         else
            mtdir(i) = .5*atan2(scrs,smns)
            if (mtdir(i) >  piotwo) mtdir(i) = mtdir(i) - 2.*piotwo
            if (mtdir(i) < -piotwo) mtdir(i) = mtdir(i) + 2.*piotwo
         endif
      endif
   enddo

   return
end subroutine equivmount2

end module equivmount



subroutine equivmount1(rgc,sxx,syy,sxy, &
     ni,il1,il2, &
     slope,xcent,mtdir)
   use equivmount, only: equivmount2
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: ni,il1,il2
   real, intent(in) :: rgc(ni), sxx(ni), syy(ni), sxy(ni)
   real, intent(out) :: slope(ni), xcent(ni), mtdir(ni)

   real :: gc(ni)

   gc(:) = - rgc(:)
   call equivmount2(gc, sxx, syy, sxy, slope,xcent, mtdir, ni)
   
   return
end subroutine equivmount1

