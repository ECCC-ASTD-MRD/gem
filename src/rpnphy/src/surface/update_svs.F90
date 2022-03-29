!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%
!!! S/P UPDATE_SVS
!
      SUBROUTINE UPDATE_SVS( WDT, WFT, WRT, LATFLW, WATFLOW, &
           WD, WF, WR, WDM, LATFLAF, DRAINAF, N )


        use svs_configs
        use sfc_options
      implicit none
!!!#include <arch_specific.hf>

!
!
      INTEGER N

      real, dimension(n,nl_svs) :: wdt, wft, wd, wf, latflw
      real, dimension(n) :: wrt, wr, wdm, latflaf, drainaf 
      real, dimension(n,nl_svs+1) :: watflow
!
!Author
!          S. Belair,M.Abrahamowicz,S.Z.Husain (2015)
!Revisions
! 001      Name (date) - Comment
!
!Object
!          Update the prognostic variables
!
!Arguments
!
!           - Input -
! WDT(NL)    soil volumetric water content in soil layer (NL layers) at time +
! WFT(NL)    frozen soil water  in soil layer (NL layers) at time +   
! WRT        water content retained by the vegetation canopy at time +
! LATFLW(NL) lateral flow (interflow) 
! WATFLOW(NL+1) water flux between soil layers and through bottom 
!
!          - Output -
! WD(NL)     updated prognostic var.: soil volumetric water content per layer
! WF(NL)     updated prognostic var.: frozen soil volum. water content per layer
! WR         updated prognostic var.: water content retained by the veg. canopy
! WDM        Mean soil moisture for the soil layers (NL soil layers)
! LATFLAF    Accum. of LATF at all levels (kg/m2 = mm)
! DRAINAF    Accum. of base drainage
      INTEGER I,K
!
!
!
      DO I=1,N
         DO K=1,NL_SVS         
            WD(I,K) = MAX ( WDT(I,K) , 0.001 )
            WF(I,K) = WFT(I,K)
         ENDDO   
         WR(I) = WRT(I)
         !     
!           Calculate mean soil moisture 
!           Soil moisture weighted by depth of each layer...        
         WDM(I) = WD(I,1) * DL_SVS(1)
         DO K=2,NL_SVS
            WDM(I) = WDM(I) + WD(I,K) * ( DL_SVS(K) - DL_SVS(K-1) )
         ENDDO
         WDM(I) = WDM(I) / DL_SVS(NL_SVS)
         !
         !       ACCUMULATE LATERAL FLOW FOR EACH SOIL LAYER
         do K=1,KHYD
            LATFLAF(i)  = LATFLAF(i) + LATFLW(I,K)       
         enddo

         !       ACCUMULATION OF DRAINAGE (BASE FLOW)
         ! Drainage is the vertical flow across the bottom of the lowermost 
         ! active layer: level = # khyd + 1
         DRAINAF(I) = DRAINAF(I) + WATFLOW(I,KHYD+1)
        
         
      END DO
!
!
      RETURN
      END
