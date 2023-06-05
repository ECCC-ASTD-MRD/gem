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
!-------------------------------------- LICENCE END --------------------------

subroutine SERACC(VS,HH,NSURF,NT,SURFACE,TSMOYHR,SRWRI)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer   NSURF
   integer   NT, TSMOYHR,SRWRI
   character(len=*) :: SURFACE(*)
   real      VS(NT,*)
   real(REAL64) :: HH(NT)

   !@Author K. Winger (UQAM 2006)
   !@Object  CORRECTS ACCUMULATORS FOR ONE STATION
   !          Accumulator get reinitialized in the model every 'P_out_moyhr' hours
   !          if the model runs in climate mode. When the time series are written
   !          at a higher frequency than that, the value of an accumulator does not
   !          correspond to the time series output interval but to the model output
   !          interval. This routine recalculates the accumulators for the time
   !          series output interval.
   !@Arguments
   !          - Input/Output -
   ! VS       time-serie values of surface variables requested
   !          - Input -
   ! HH       hour relative to start of the forecast (central time)
   ! NSURF    number of surface variables requested
   ! NT       timestep number
   ! SURFACE  names of time serie surface variables requested
   ! TSMOYHR  output interval for model output (in hours)
   ! SRWRI    output interval for time series output (in seconds)

   integer KSURF,A,T
   real    ACCUM

#include "acclist.cdk"

   !     If the model output interval is shorter than the
   !     time series output interval **OR** if TSMOYHR is
   !     not a multiple of SRWRI, nothing can be done.
   if (TSMOYHR*3600.lt.SRWRI .or. mod(TSMOYHR*3600,SRWRI).ne.0) return

   !     If the model and time series output interval is the same
   !     only the second output step must be taken care of
   if (TSMOYHR*3600.eq.SRWRI .and. NT.ge.2) then

      !       If the time series output interval is larger than delta t
      if ( NT .ge. 3 ) then
         if ( HH(2)-HH(1) .ne. HH(3)-HH(2) ) then

            !           Loop over all surface variables
            do KSURF=1,NSURF

               !             Loop over all accumulators
               do A=1,NBR

                  !               If variable is an accumulator substract first time step
                  !               from second time series output step
                  !               (The model always outputs the first time step in each job
                  !                regardless of the output interval.)
                  if (SURFACE(KSURF).eq.PERMIS(2,A)) then
                     VS(2,KSURF) = VS(2,KSURF) - VS(1,KSURF)
                  end if

               enddo
            enddo
         end if
      end if

      !     If the model output interval is longer than the time series output interval
      !     recalculate all accumulators
   else

      !       Loop over all surface variables
      do KSURF=1,NSURF

         !         Loop over all accumulators
         do A=1,NBR

            !           If variable is an accumulator
            if (SURFACE(KSURF).eq.PERMIS(2,A)) then
               ACCUM = 0.0

               !             Loop over all output time steps
               do T=1,NT
                  VS(T,KSURF) = VS(T,KSURF) - ACCUM
                  if (modulo(HH(T),dble(TSMOYHR)).lt.-0.001 .or. &
                       modulo(HH(T),dble(TSMOYHR)).gt.0.001) then
                     ACCUM = ACCUM + VS(T,KSURF)
                  else
                     ACCUM = 0.0
                  endif
               enddo
            endif
         enddo
      enddo
   end if


   return
end subroutine SERACC
