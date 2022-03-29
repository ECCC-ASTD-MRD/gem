!-------------------------------------- LICENCE BEGIN ------------------------
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
