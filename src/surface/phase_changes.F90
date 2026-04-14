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
!**S/P PHASE_CHANGES
!
module phase_changes_mod
  implicit none
  public
contains
      SUBROUTINE PHASE_CHANGES ( DT, WTG, LAI, &
           CAP, DZ_FL, WSAT, PSISAT, BCOEF, &
           TG, WF, WFL_ICE, WFLT,WRMAX_FL, TFL, PFLHCAP, WDT,WFT, WDTT, DELWATER, DELICE, &
           FCD, N, PHASEF, PHASEM, DELTAT, APPHEATCAP, TMAX  )
!    -----------------------------------------------------------------------
!
      use sfc_options
      use svs_configs
      use tdpack
      use svs2_tile_configs

      implicit none
! #include <arch_specific.hf>
!
!    declarations of arguments
      INTEGER N,TRNCH
      REAL DT, LAI(N), FCD(N,NL_SVS), BCOEF(N,NL_SVS), PSISAT(N,NL_SVS), WSAT(N,NL_SVS)
      REAL TG(N,NL_SVS), WDT(N,NL_SVS), WDTT(N,NL_SVS), WF(N,NL_SVS), WFT(N,NL_SVS), CAP(N,NL_SVS)
      REAL DELWATER(N,NL_SVS), DELICE(N,NL_SVS), TFL(N), WFL_ICE(N), WFLT(N), PFLHCAP(N), WRMAX_FL(N)
      REAL PHASEF(N,NL_SVS), PHASEM(N,NL_SVS),DELTAT(N,NL_SVS), APPHEATCAP(N,NL_SVS), TMAX(N,NL_SVS)
      REAL WTG(N,svs2_tilesp1), DZ_FL(N)
!
!    ------------------------------------------------------------------------
!Author
!          A. Boone :   Meteo-France 
!                       ice_soildif.F90 subroutine of SURFEX
!
!          Adapted by N. Gauthier for SVS code
!Revisions
! 001      Include refreezing/melting in forest litter layer (N.Leroux)
!
!Object
!     This subroutine calculates soil water phase changes using the
!     available/excess energy approach. Soil temperature and volumetric
!     ice content are adjusted due to phase changes. See the references
!     Boone et al., 39, JAM, 2000 and Giard and Bazile, 128, MWR, 2000.
!     NOTE that more recently a modification was made: freeze/thaw follows
!     a relationship between liquid water and temperature derriving from
!     the Clausius Clapeyron Eq. This results in little to no freezing for
!     sufficiently dry but cold (below freezing) soils. Scatter about this
!     curve results due to 'phase change efficiencies' and the surface 
!     insolation coefficient.
!     
!Arguments
!
!          - Input 
! DT           timestep
! WTG       Weights for SVS2 surface types as seen from GROUND
! FCD(NL_SVS)  root fraction within soil layer (NL_SVS soil layers) bus(x(frootd ,1,1)
! LAI          leaf area index
! NL_SVS       number of soil layers
! WDT(NL_SVS)  soil volumetric water content in soil layers (m3/m3) 
! D(NL_SVS)    depth of soil layer (NL_SVS soil layers) in METERS
!              soil heat capacity [J/(m3 K)]
! WSAT         soil porosity (m3/m3) or volumetric water content at saturation
! CAP(NL_SVS)  soil heat capacity [J/(m3 K)]
! PSISAT       matric potential at saturation (m)
! TG(N)        soil temperature (K)
! WF(NL_SVS)   soil volumetric ice content (m3/m3)
! WFL_ICE      Ice retained in the forest litter layer (kg/m2)
! WFLT         Liquid water retained in the forest litter layer (kg/m2)
! WRMAX_FL         Max water holding capacity in the forest litter layer (kg/m2)
! TFL          Temperature of the forest litter layer (K)
! PFLHCAP   Heat mass for the forest litter layer (J K-1 m-2)
! DZ_FL        Thickness of the forest litter layer (m)
! 
!           - Output -
! TG          soil temperature (K)
! WFT         soil volumetric ice content (m3/m3)
! WDTT        soil volumetric liquid water content (m3/m3)

!           - Local -
! KVEG(NL_SVS) effect of surface layer insolation on phase changes
!              Giard and Bazile (2000): non-dimensional
! TAUICE       soil phase change characteristic time scale (s)
!
!   -------------------------------------- ----------------------------------------


!     declarations of local variables 
      INTEGER I, K
!

      real WORK, WORKLOG, WLMAX, PSIMAX, PSI, ZK, WGM, WGIM, TGM
!      real PHASEM, PHASEF, PHASE, DELTAT, APPHEATCAP, PSISATZ
      real PSISATZ
      real, DIMENSION(N) :: PHI_FL, & ! The phase change rate for forest litter (kg m−2 s−1)
                            EXCESS_FL, & ! Variable to adjust mass balance due to refreezing/melting of forest litter
                            TFLT  ! Temperature of the forest litter layer at the end of the time step (K)
!
      REAL, DIMENSION(N,NL_SVS) :: EXCES
!
!
! Initialization:
! ---------------
      EXCES(:,:)=0.0
      EXCESS_FL(:)=0.0

!     1.  SURFACE LAYER VEGETATION INSULATION COEFFICIENT
!         -----------------------------------------------

!        a activer plus tard au besoin
!!$      DO I=1,N
!!$         DO K=1,NL_SVS
!!$            KVEG(I,K) = (1.0-INSOLFRZ_VEG) * (1.0-((LAI(I)* (1.0 - FCD(I,K)) )/INSOLFRZ_LAI))
!!$            KVEGBARE(I,K) = 1.0  !no vegetation
!!$         ENDDO
!!$      ENDDO
      ZK = 1.0

!     2.  SOIL ICE EVOLUTION COMPUTATION
!         ------------------------------
      DO I=1,N
         DO K=1,NL_SVS  
!
            WGIM = WF(I,K)
            WGM  = WDT(I,K)
            TGM  = TG(I,K)

            PSISATZ = -PSISAT(I,K)

!                      The maximum liquid water content ( wlmax ) as
!                      as function of temperature (sub-freezing)
!                      based on Gibbs free energy (Fuchs et al., 1978):
!
            PSIMAX  = MIN(PSISATZ, CHLF * (TG(I,K) - TRPL) / (GRAV * TG(I,K)))

            WORK    = PSIMAX/PSISATZ
            WORKLOG = LOG(WORK)/BCOEF(I,K)
            WLMAX   = WSAT(I,K)*EXP(-WORKLOG)
! 
!                      Calculate maximum temperature for ice based on Gibbs free energy:
!                      first compute soil water potential using Brook and Corey (1966) model:
!                      psi=mpotsat*(w/wsat)**(-bcoef)
!
            WORK = WDT(I,K)/WSAT(I,K)
            WORKLOG  = BCOEF(I,K)*LOG(WORK)
            PSI   = PSISATZ*EXP(-WORKLOG)
!
!                      Maximum temperature at which soil ice is present
            TMAX(I,K) = CHLF*TRPL/(CHLF-GRAV*PSI)

            DELTAT(I,K) = TG(I,K) - CHLF*TRPL/(CHLF-GRAV*PSI)

!  for tests
!!$           DELTAT(I,K) = TG(I,K) - TRPL
!
!                      Compute apparent heat capacity. This is considered
!                      only when there is available energy (cold) and liquid water
!                      available...freezing front.
!                      This also has the secondary effect of increasing numerical stability
!                      during freezing (as there is a strong temperature dependence) by
!                      i) potentially significantly increasing the "apparent" heat capacity and
!                      ii) this part is also treated implicitly herein.
!
            WORK = (CICE*RHOI/(CHLF*RHOW))*ZK*MAX(0.0,-DELTAT(I,K))

!        
            APPHEATCAP(I,K)=0.0
            PHASEM(I,K) = 0.0
            PHASEF(I,K) = 0.0

            IF (DELTAT(I,K)<0.0.AND.WGM>=WLMAX.AND.WORK>=MAX(0.0, WGM-WLMAX)) THEN
               APPHEATCAP(I,K) = -(TRPL*RHOW*CHLF*CHLF/GRAV) * WLMAX/(PSIMAX*BCOEF(I,K)*TGM*TGM)
            ENDIF

!
!                     * MELT *   ice if energy and ice available:  WGIM = WF
            PHASEM(I,K)  = (DT/TAUICE)*MIN(ZK*CICE*RHOI*MAX(0.0,DELTAT(I,K)), WGIM*CHLF*RHOW)
!
!                     * FREEZE * liquid water if energy and water available:
            PHASEF(I,K)  = (DT/TAUICE)*MIN(ZK*CICE*RHOI*MAX(0.0, -DELTAT(I,K)), MAX(0.0, WGM-WLMAX)*CHLF*RHOW)
!
!                     Update heat content if melting or freezing
            TG(I,K) = TGM + (PHASEF(I,K) - PHASEM(I,K))/(CAP(I,K)+APPHEATCAP(I,K))

!                     Adjust ice and liquid water contents (m3/m3) accordingly :
            WFT(I,K) = WGIM + (PHASEF(I,K) - PHASEM(I,K))/(CHLF*RHOW)   
            WDTT(I,K) = WGM  - (PHASEF(I,K) - PHASEM(I,K))/(CHLF*RHOW)

            DELICE(I,K)   =  (PHASEF(I,K) - PHASEM(I,K))/(CHLF*RHOW) 
            DELWATER(I,K) = - DELICE(I,K)

         ENDDO
      ENDDO

! Prevent keeping track of very small numbers for ice content: (melt it)
! and conserve energy:
!
      DO I=1,N
         DO K=1,NL_SVS 
            IF(WFT(I,K)>0.0.AND.WFT(I,K)<1.0E-6)THEN
               WDTT  (I,K) = WDTT(I,K) + WFT(I,K)
               EXCES(I,K) = EXCES(I,K) + WFT(I,K)
               WFT  (I,K) = 0.0
               DELWATER(I,K) = DELWATER(I,K) + EXCES(I,K)
               DELICE(I,K)   = -DELWATER(I,K)
            ENDIF
            TG (I,K) = TG(I,K) - EXCES(I,K)*CHLF*RHOW/CAP(I,K) 
         ENDDO
      ENDDO
!

! Refreezing and melting in the forest litter layer
!
      IF (LFORLIT ) THEN
            DO I=1,N
                  IF (WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN
                        IF (TFL(I) .gt. TRPL) THEN  ! Melting

                              PHI_FL(I) = MIN(WFL_ICE(I), (RHOI*CICE * DZ_FL(I) / CHLF * (TFL(I) - TRPL))) / (TAUICE)

                              WFLT(I) = WFLT(I) + DT * PHI_FL(I)
                              WFL_ICE(I) = WFL_ICE(I) - DT * PHI_FL(I)
                              TFLT(I) = TFL(I) - DT / PFLHCAP(I) * CHLF * PHI_FL(I)

                              ! Liquid water in forest litter should not exceed maximum holding capacity
                              IF(WFLT(I) .gt. WRMAX_FL(I))THEN
                                    EXCESS_FL(I) = (WRMAX_FL(I) - WFLT(I)) / DT
                                    WFLT(I) = WFLT(I) - DT * EXCESS_FL(I)
                                    WFL_ICE(I) = WFL_ICE(I) + DT * EXCESS_FL(I)
                                    TFLT(I) = TFLT(I) + DT / PFLHCAP(I) * CHLF * EXCESS_FL(I)
                              ENDIF

                              ! Prevent keeping track of very small numbers for ice content: (melt it)
                              ! and conserve energy:
                              !
                              IF(WFLT(I) > 0. .AND. WFL_ICE(I) < 1.0E-4)THEN
                                    WFLT(I) = WFLT(I) + WFL_ICE(I) 
                                    EXCESS_FL(I) = WFL_ICE(I) / DT
                                    WFL_ICE(I) = 0.
                                    TFL(I) = TFL(I) - DT / PFLHCAP(I) * CHLF * EXCESS_FL(I)
                              ENDIF

                        ELSE ! Freezing

                              PHI_FL(I) =  MIN(WFLT(I), (RHOI*CICE * DZ_FL(I) / CHLF * (TRPL - TFL(I)))) / (TAUICE)

                              WFLT(I) = WFLT(I) - DT * PHI_FL(I)
                              WFL_ICE(I) = WFL_ICE(I) + DT * PHI_FL(I)
                              TFLT(I) = TFL(I) + DT / PFLHCAP(I) * CHLF * PHI_FL(I)

                              ! Liquid water in forest litter should not be below a chosen threshol of 0.01 mm
                              IF(WFLT(I) < 0.01)THEN
                                    EXCESS_FL(I) = (0.01 - WFLT(I)) / DT
                                    WFLT(I) = 0.01
                                    WFL_ICE(I) = WFL_ICE(I) - DT * EXCESS_FL(I)
                                    TFLT(I) = TFLT(I) - DT / PFLHCAP(I) * CHLF * EXCESS_FL(I)
                              ENDIF
                        ENDIF
                  ENDIF
            ENDDO


            DO I=1,N
                  TFL(I) = TFLT(I)
            ENDDO
      ENDIF
!

      RETURN
    END SUBROUTINE PHASE_CHANGES
end module phase_changes_mod
