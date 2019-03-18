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
!**S/R SERACC  CORRECTS ACCUMULATORS FOR ONE STATION
!
      SUBROUTINE SERACC(VS,HH,NSURF,NT,SURFACE,TSMOYHR,SRWRI)
!
      implicit none
#include <arch_specific.hf>
      INTEGER   NSURF
      INTEGER   NT, TSMOYHR,SRWRI
      CHARACTER *(*) SURFACE(*)
      REAL      VS(NT,*)
      REAL*8    HH(NT)
!
!Author
!          K. Winger (UQAM 2006)
!
!Revision
! 001      None
!
!Object
!          Accumulator get reinitialized in the model every 'P_out_moyhr' hours
!          if the model runs in climate mode. When the time series are written
!          at a higher frequency than that, the value of an accumulator does not
!          correspond to the time series output interval but to the model output
!          interval. This routine recalculates the accumulators for the time
!          series output interval.
!
!Arguments
!
!          - Input/Output -
! VS       time-serie values of surface variables requested
!
!          - Input -
! HH       hour relative to start of the forecast (central time)
! NSURF    number of surface variables requested
! NT       timestep number
! SURFACE  names of time serie surface variables requested
! TSMOYHR  output interval for model output (in hours)
! SRWRI    output interval for time series output (in seconds)
!
!
!*
      INTEGER KSURF,A,T,I
      REAL    ACCUM
!
#include "acclist.cdk"
!
!     If the model output interval is shorter than the
!     time series output interval **OR** if TSMOYHR is
!     not a multiple of SRWRI, nothing can be done.
      IF (TSMOYHR*3600.LT.SRWRI .OR. MOD(TSMOYHR*3600,SRWRI).NE.0) RETURN

!     If the model and time series output interval is the same
!     only the second output step must be taken care of
      IF (TSMOYHR*3600.EQ.SRWRI .AND. NT.GE.2) THEN
!
!       If the time series output interval is larger than delta t
        IF ( NT .ge. 3 ) THEN
          IF ( HH(2)-HH(1) .NE. HH(3)-HH(2) ) THEN
!
!           Loop over all surface variables
            DO 10 KSURF=1,NSURF
!
!             Loop over all accumulators
              DO 10 A=1,NBR
!
!               If variable is an accumulator substract first time step
!               from second time series output step
!               (The model always outputs the first time step in each job
!                regardless of the output interval.)
                IF (SURFACE(KSURF).EQ.PERMIS(2,A)) THEN
                  VS(2,KSURF) = VS(2,KSURF) - VS(1,KSURF)
                END IF
!
10          CONTINUE
          END IF
        END IF
!
!     If the model output interval is longer than the time series output interval
!     recalculate all accumulators
      ELSE
!
!       Loop over all surface variables
        DO 20 KSURF=1,NSURF
!
!         Loop over all accumulators
          DO 20 A=1,NBR
!
!           If variable is an accumulator
            IF (SURFACE(KSURF).EQ.PERMIS(2,A)) THEN
              ACCUM = 0.0
!
!             Loop over all output time steps
              DO 30 T=1,NT
                VS(T,KSURF) = VS(T,KSURF) - ACCUM
                IF (modulo(HH(T),dble(TSMOYHR)).LT.-0.001 .OR. &
                    modulo(HH(T),dble(TSMOYHR)).GT.0.001) THEN
                  ACCUM = ACCUM + VS(T,KSURF)
                ELSE
                  ACCUM = 0.0
                ENDIF
30            CONTINUE
            ENDIF
20      CONTINUE
      END IF
!
!
      RETURN
      END
