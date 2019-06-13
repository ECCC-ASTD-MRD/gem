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
!**S/R modhui3  -  MODIFICATION DE HU EN HAUT DE 300 MB
!

SUBROUTINE modhui3(HU,TX,PS,SWPH,NI,NK,N)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>
   INTEGER NI, NK, N
   REAL HU(NI,NK), TX(NI,NK)
   REAL PS(NI,*)
   LOGICAL SWPH

   if (ni < 0) print *,ps(:,1)  !#Note: to keep old itf with PS w/o the warning
   call modhui4(HU,TX,SWPH,NI,NK,N)

   return
end SUBROUTINE modhui3


SUBROUTINE modhui4(HU,TX,SWPH,NI,NK,N)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>

   INTEGER :: NI, NK, N
   REAL :: HU(NI,NK), TX(NI,NK)
   REAL :: PN1,PN2
   LOGICAL :: SWPH
 
   !Author
   !     C. GIRARD - MAR 93
   !
   !Revision
   ! 001      B. Dugas (Apr 1996) - Modifications for hybrid
   !                                coordinate in sef model
   ! 002      G. Pellerin (Mar 1996)   - Revised interpolation above 300mb
   !                                     level with limit above 150mb
   ! 003      B. Bilodeau (Jan 2001) - Automatic arrays
   !
   !Object
   !     To recalculate the specific humidity over 300MB.
   !     by using the zeroth order interpolation
   !     and keeping one of the following constant:
   !     A) The relative humidity when the temperature decreases
   !        to simulate what happens below the tropopause.
   !     or
   !     B) The specific humidity when the temperature increases
   !        to simulate what happens above the tropopause.
   !
   !Arguments
   !
   !          -Output-
   ! HU       Specific humidity in kg/kg
   !
   !          -Input-
   ! TX       temperature in K
   ! PS       pressure in pa
   ! SWPH     if .TRUE., water and ice phase considered
   !          if .FALSE.,water phase only for all temperatures.
   ! NI       Horizontal dimension
   ! NK       Vertical dimension
   ! N        Number of points to process
   !*
   !--------------------------------------------------------------------

   REAL :: Qsat,Qnot
   REAL :: ALPHA, HREL, HRPRIM, QPRIM
   Real, Dimension(ni,nk) :: PN
   Real, Dimension(ni,nk) :: PN0

   INTEGER K,K0, I
!--------------------------------------------------------------------

      DO 10 K0=NK-1,1,-1
         K=K0
         DO I=1,N
            PN0(I,K)=PN(I,K)
         ENDDO
         K=K0+1
         K=K0
         Qnot = 2.5E-6
         DO 300 I=1,N
!--------------------------------------------------------------------
           IF(PN0(I,K).LT.3.E4) THEN
!--------------------------------------------------------------------
            IF(SWPH) THEN
              Qsat=FOQST(TX(I,K),PN0(I,K))
            ELSE
              Qsat=FOQSA(TX(I,K),PN0(I,K))
            ENDIF
!--------------------------------------------------------------------
           IF(PN0(I,K).LT.1.50E4) THEN
!--------------------------------------------------------------------
            HU(I,K) = min(Qnot,0.8*qsat)
           ELSE

            PN1 = PN0(I,K)
            PN2 = PN(I,K)
            ALPHA = (PN1 - 1.5E4) /1.5E4
            QPRIM = ALPHA * HU(I,K+1) + (1.-ALPHA) * Qnot

            IF(SWPH) THEN
              HREL = FOHR ( HU(I,K+1) , TX(I,K+1), PN2)
            ELSE
              HREL = FOHRA( HU(I,K+1) , TX(I,K+1), PN2)
            ENDIF

            HRPRIM = ALPHA * HREL + (1.-ALPHA) * 0.8
            HU(I,K) = min( QPRIM , HRPRIM * Qsat)
           ENDIF
!--------------------------------------------------------------------
           ENDIF
!--------------------------------------------------------------------
300      CONTINUE

10    CONTINUE

   RETURN
END SUBROUTINE modhui4
