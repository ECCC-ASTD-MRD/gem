

subroutine SERPAT2(ST , VT , T , NT , NK , &
     NSURF2, NPROF2, NSTAT2)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   integer T,NT,NK,NSURF2,NPROF2,NSTAT2
   real ST(NT,NSURF2,NSTAT2),VT(NK,NT,NPROF2,NSTAT2)

   !@Author R. Benoit
   !@Revision
   ! 001      V.Alex.(Feb 87)
   !                   Documentation
   ! 002      B. Bilodeau (Jul 2002) - Optimization
   !@Object TRANSFERER DE POINTS A TEMPS
   !        reorder the time-series with respect to time
   !        instead of station number
   !@Arguments
   !          - Output -
   ! ST       table of surface variables ordered in time
   ! VT       table of profile variables in time
   !          - Input -
   ! T        time-pointer referring to timestep
   ! NT       number of timesteps
   ! NK       vertical dimension
   ! NSURF2   number of surface variables
   ! NPROF2   number of profile variables
   ! NSTAT2   number of stations

#include "series.cdk"

   integer K,L,M

   if (NSTAT .le. 0) return
   if (.not. INITOK) return

   do L=1,NSTAT

      do M=1,NSURF
         ST(T,M,L) = SERS(L,M)
      enddo

      do M=1,NPROF
         do  K=1,NK
            VT(K,T,M,L) = SERP(K,L,M)
         enddo
      enddo
   enddo

   return
end subroutine SERPAT2
