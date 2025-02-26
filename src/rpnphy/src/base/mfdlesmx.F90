
subroutine MFDLESMX(RES,TT,FF,DF,NI,NK)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack, only: FOEWAF, FESIF, AERK1W, AERK1I, FESMXX, FDLESMXX
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer NI,NK
   real RES(NI,NK),TT(NI,NK),FF(NI,NK),DF(NI,NK)
   !
   !@Author A. Plante (May 2003), based on FDLESMX from J. Mailhot
   !@Object MIXED PHASE SATURATION VAPOR PRESSURE CALCULATION
   !        calculate mixed phase saturation vapor pressure
   !@Arguments
   !          - Output -
   ! RES      mixed phase saturation vapor pressure
   !          - Input -
   ! TT       temperature
   ! FF       ice fraction
   ! DF       value of derivative w/r to T for saturation w/r to ice or liquid
   ! NI       horizontal dimension
   ! NK       vertical dimension
   !*
   integer i,k
   real(REAL64) :: foewad, fesid, fesmxd

   !***********************************************************************
   do k=1,nk
      do i=1,ni
         foewad = aerk1w * exp(dble(FOEWAF(tt(i,k))))
         fesid  = aerk1i * exp(dble(FESIF(tt(i,k))))
         fesmxd = FESMXX(ff(i,k),fesid,foewad)
         res(i,k) = FDLESMXX(tt(i,k),ff(i,k),df(i,k),foewad,fesid,fesmxd)
      enddo
   enddo

   return
end subroutine MFDLESMX
