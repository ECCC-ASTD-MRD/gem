!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%
!!! S/P UPDATE_SVS
!
module update_svs2_mod
  implicit none
  public
contains
      SUBROUTINE UPDATE_SVS2( WDT, WFT, WR_VLT, WR_VHT, WFLT, &
           WD, WF, WR_VL, WR_VH, WFL, WDM, N )


        use svs_configs
        use sfc_options
      implicit none
!!!#include <arch_specific.hf>

!
!
      INTEGER N

      real, dimension(n,nl_svs) :: wdt, wft, wd, wf
      real, dimension(n) :: wr_vlt, wr_vht, wr_vl, wr_vh, wdm, wflt, wfl
     
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
! WR_VLT     water content retained by the low vegetation canopy at time +
! WR_VHT     water content retained by the high vegetation canopy at time +
! WFLT       water content retained by the forest litter at time +
!
!          - Output -
! WD(NL)     updated prognostic var.: soil volumetric water content per layer
! WF(NL)     updated prognostic var.: frozen soil volum. water content per layer
! WR_VL      updated prognostic var.: water content retained by the low veg. canopy
! WR_VH      updated prognostic var.: water content retained by the high veg. canopy
! WFL        updated prognostic var.: water content retained by the forest litter
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
         WR_VL(I) = WR_VLT(I)
         WR_VH(I) = WR_VHT(I)
         IF (LFORLIT) THEN
            WFL(I) = WFLT(I)
          ENDIF
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
    END SUBROUTINE UPDATE_SVS2
  end module update_svs2_mod
