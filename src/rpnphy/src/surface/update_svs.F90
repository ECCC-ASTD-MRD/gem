!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%
!!! S/P UPDATE_SVS
!
      SUBROUTINE UPDATE_SVS( WDT, WFT, WRT, &
           WD, WF, WR, WDM, N )


        use svs_configs
        use sfc_options
      implicit none
!!!#include <arch_specific.hf>

!
!
      INTEGER N

      real, dimension(n,nl_svs) :: wdt, wft, wd, wf
      real, dimension(n) :: wrt, wr, wdm
     
!
!Author
!          S. Belair,M.Abrahamowicz,S.Z.Husain (2015)
!Revisions
! 001      E. Gaborit, 2022 - Move accumulators drainaf, latlaf to avoid double
!             counting (to sfc_calcdiag.F90 - M.A.)    
!     
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
!
!          - Output -
! WD(NL)     updated prognostic var.: soil volumetric water content per layer
! WF(NL)     updated prognostic var.: frozen soil volum. water content per layer
! WR         updated prognostic var.: water content retained by the veg. canopy
! WDM        Mean soil moisture for the soil layers (NL soil layers)
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
         
      END DO
!
!
      RETURN
      END
