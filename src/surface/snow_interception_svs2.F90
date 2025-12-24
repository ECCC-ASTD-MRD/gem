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

      SUBROUTINE SNOW_INTERCEPTION_SVS2 (DT, TVEG, T, HU, PS, WIND_TOP, ISWR,  RHOA,    &
                     RR, SR, SNCMA, WRMAX_VH, SKYVIEW, ESUBSNC,SUBSNC_CUM, LAIVH, WTG, HM_CAN,   &
                     VGH_DENS, SCAP, WR_VH, RR_VEG, SR_VEG,UNLOAD_RATE, FCANS, N)


      use tdpack
      use sfc_options
      use svs_configs
      USE MODE_THERMOS
      USE MODD_CSTS
      USE CANOPY_CSTS, only : ALPHA, ZVENT, RADIUS_ICESPH, CLUMPING, &
            ALBEDO_ICESPH, KS, FRACT, MWAT, RGAZ, CICE, TCNC, TCNM, ZBETA
      use svs2_tile_configs
      implicit none

      INTEGER N


      REAL SNCMA(N), RR(N), SR(N),  RR_VEG(N), SR_VEG(N), UNLOAD_RATE(N), ISWR(N)
      REAL PS(N), RHOA(N), HU(N), WTG(N, svs2_tilesp1), SKYVIEW(N)
      REAL LAIVH(N), TVEG(N), WIND_TOP(N),VGH_DENS(N), T(N)
      REAL ESUBSNC(N), SUBSNC_CUM(N), SCAP(N), FCANS(N), HM_CAN(N)
      REAL DT, SNCMA_INI(N), WRMAX_VH(N), WR_VH(N) 

!
!
!Author
!             V. Vionnet (April 2023)
!
!Object
!             Module to simulate interception of snow by high vegetation and
!             sublimation and unloading of intercepted snow using the
!             approach proposed in the FSM2 model
!
!Reference    Hedstrom and Pomeroy (1998) for the snow interception module (referred below as HP98)
!
!Arguments

!             - Input/Output (Prognostic variables)
!
! SNCMA       mass of intercepted snow in the canopy [kg/m2]

!             - Input (Forcing) -

! RR          liquid precipitation rate at the top of the high-vegetation [kg/m2/s]
! SR          solid  precipitation rate at the top of the high-vegetation [kg/m2/s]
! T           air temperature within the high vegetation [K]
! HU          Specific humidity of air within the high vegetation [kg kg-1]
! TVEG           Canopy temperature within the high vegetation [K]
! WIND_TOP    Wind speed at canopy top [m/s]
! PS          Surface pressure [Pa]
! ISWR        Solar radiation incident at the top of the high-vegetation [W m-2]
! ESUBSNC      sublimation rate from incercepted snow [kg m-2 s-1]
!
!             - Input (other parameters)
! DT          time step [s]
! LAIVH       Vegetation leaf area index for HIGH vegetation only [m2/m2]
! VGH_DENS    tree cover density in  HIGH vegetation [0-1]
! SCAP       Vegetation layer snow capacities (kg m-2)
! HM_CAN    Heat mass for the high vegetation layer (J K-1 m-2)
! WRMAX_VH          max water content retained on hihj vegetation [kg/m2]
!
!             - Input-Output
! WR_VH          water content retained by high vegetation canopy [kg/m2]
! SUBSNC_CUM   Cumulated mass loss due to sublimation of incercepted snow [kg m-2]
!             - Output
!
! RR_VEG       liquid precipitation rate below high-vegetation [kg/m2/s]
! SR_VEG       solid  precipitation rate below high-vegetation from throufall [kg/m2/s]
! UNLOAD_RATE  rate of solid unloading below high-vegetation [kg/m2/s]       
! FCANS        Canopy layer snowcover fractions from FSM2

      INTEGER I

      REAL TUNL  ! Unloading time constant


      REAL, DIMENSION(N) :: DRIP_CPY ! Mass of liquid water dripping below the canopy (kg/m^2)
      REAL, DIMENSION(N) :: NET_SNOW ! Snow mass sent to the snowpack below the canopy  (kg/m^2)
      REAL, DIMENSION(N) :: UNLOAD_MASS ! Solid mass unloading below the canopy  (kg/m^2)
      REAL, DIMENSION(N) :: PQSAT ! Specific humidity at saturation [kg kg-1]
      REAL, DIMENSION(N) :: PPSAT ! Vapor presure at saturation [Pa]
      REAL, DIMENSION(N) :: HM_CAN_INI !  Heat mass of the canopy at the beginning of the subroutine

      REAL  SNW_UNLOAD !  Unloaded mass in solid form (kg/m^2)
      REAL UNLOAD !  Unloaded mass (kg/m^2)
      REAL DIRECT_SNOW ! Snow mass directly transferred through the canopy (kg/m^2)
      REAL MELT, REFREEZE
      REAL INTCPT !  Canopy interception (kg/m^2)
      REAL   SUB_CPY ! Intercepted snow mass lost by sublimation (kg/m^2)
      REAL VSUBL ! Ventilation wind speed used in the calcultation of sublimation

      REAL XI2,EXT2, WINDEXT2, EVAP, LAMBDA, SSTAR, A1, B1,  LS, J
      REAL C1, SIGMA, SVDENS, SUB_RATE,SUB_POT, MPM, CE
      REAL NU, NR, MU,DVAP

      !
      ! 0. Initialise variables
      !
      ! Heat mass of the canopy without the snow
      HM_CAN_INI(:) = HM_CAN(:)

      !
      ! 1. Evolution of snow mass intercepted by the canopy
      !

      ! Initialize snowfall mass directly transfered through the canopy
      NET_SNOW(:) = 0.

      ! Initialize mass of liquid water dripping below the canopy 
      DRIP_CPY(:) = 0.

      ! Initialize solid mass unloading below the canopy
      UNLOAD_MASS(:) = 0.
      
      SNCMA_INI(:) = SNCMA(:)

      ! Compute specific humidity and vapor pressure at saturation
      PQSAT(:) = QSATI( T(:), PS(:))
      PPSAT(:) = PSAT(T(:))

      DO I=1,N
         IF (CANO_REF_FORCING .EQ.'O2F') THEN
	    ! Initialize sublimation of intercepted snow
            ESUBSNC(I) = 0.

            ! For O2F, do the phase change here before intercepted snow sublimation is calculated
            IF(TVEG(I) .LE. 273.15) THEN ! refreezing of liquid intercepted waterto intercepted snow 
               SNCMA(I) = SNCMA(I) + WR_VH(I)
               WR_VH(I) = 0.
            ENDIF
         ENDIF

         IF (WTG(I,indx_svs2_vh).GE.EPSILON_SVS) THEN
            IF(SNCMA(I) .GT. 0. .OR. SR(I) .GT. 0.) THEN  ! Snow is present on the canopy or occurrence of snowfall


               !!!!!!!!
               ! Interception Code (Simple version of HP98)
               !!!!!!!!

               ! Compute Canopy interception, check Pomeroy VGH_DENS position
               INTCPT = (SCAP(I) - SNCMA(I))*(1. - exp(-(1.-SKYVIEW(I))* SR(I)*DT/SCAP(I))) 

               ! Update intercepted snow mass
               SNCMA(I) = SNCMA(I) + INTCPT


               ! Update the amount of snow that falls below high vegetation
               DIRECT_SNOW  = SR(I)*DT - INTCPT

               ! Remove mass of intercepted lost by sublimation
               ! Sublimation calculated in the previous time step so should be
               ! applied to snow already on canopy

               IF (CANO_REF_FORCING .EQ.'ABV') THEN

                  ! Sublimation rate of intercepted snow on canopy
                  ! is calculated in ebudget
                  ! Sublimation of intercepted snow is computed in watsurf_budget when ABV is selected

               ELSE ! Open 2 Forest

                  ! Fraction of the entire forest height [-]
                  XI2 = 1.-ZVENT

                  ! Canopy wind speed extinction coefficient [-]
                  ! From Marke et al., (2016); Liston and Elder (2006), consistent with canopy_met_svs2 and drag_svs2
                  ! VEG_DENS is included here as the wind speed would be higher in more sparse forest
                  EXT2 = ZBETA * CLUMPING * LAIVH(I) * VGH_DENS(I)

                  ! Computation of ventilation wind speed of intercepted snow derived from above-caopny wind speed  [Eq 8 in EL10]
                  ! Estimated within canopy wind speed at fraction XI2 of the entire tree height [Eq 8 in EL10]  [m s-1]
                  VSUBL = WIND_TOP(I) * EXP(-1. * EXT2 * XI2)

                  ! Vaport density at saturation
                  SVDENS = PPSAT(I) * MWAT / (RGAZ * T(I))

                  ! Derive partial vapor pressure from speficic humidity and pressure [Pa]
                  EVAP = HU(I) * PS(I) / (0.622+0.378 * HU(I) )

                  ! Sutherland's equation for kinematic viscosity
                  MU=1.8325e-5*416.16/( T(I)+120)*(T(I)/296.16)*SQRT(T(I)/296.16)/RHOA(I)

                  ! Compute thermal conductivity of air [J m-1 s-1 K-1]
                  LAMBDA = 0.000063 * T(I) + 0.00673

                  ! Compute latent heat of sublimation [J kg-1]
                  ! We use here the formulation of Roges and Yau (1989) as in Harder and Pomeroy (2013)
                  ! Note that the canopy module in CRHM used a constant value
                  ! TODO: Crocus approach should be used for consistency
                  IF(TVEG(I) .LT. 273.15) THEN
                     LS =  1000.0 * (2834.1 - 0.29 *(T(I)-273.15) - 0.004*(T(I)-273.15)**2.)
                  ELSE
                     LS = 1000.0 * (2501.0 - (2.361 * (T(I)-273.15)))
                  ENDIF

                  ! Compute diffusivity of water vapour in air [m2 s-1]m 'HP13'
                  DVAP = 2.063e-5 * (T(I)/273.15)**1.75

                  ! Compute Reynolds Number     [-]
                  NR  = 2.0 * RADIUS_ICESPH * VSUBL / MU

                  ! Compute the Nusselt Number  [-]
                  NU = 1.79 + 0.606 * SQRT(NR)

                  ! Incoming shortwave radiation to the ideal snow particle [W m-2]
                  SSTAR = PI* RADIUS_ICESPH**2. * (1.0 - ALBEDO_ICESPH) * ISWR(I)

                  ! Compute the term used in equation for the sublimation rate of an ice sphere
                  ! from Thorpe and Masson (1966)
                  A1  = LAMBDA * T(I) *NU
                  B1 = ( LS * MWAT /(RGAZ * T(I)) ) - 1.0
                  J = B1/A1
                  C1  = 1.0 / (DVAP * NU * SVDENS )

                  ! Compute water vapour deficit with respect to  ice  [-]
                  SIGMA = EVAP/PPSAT(I) - 1.0

                  ! Compute the mean mass of a snow particle assuming that the size distribution of
                  ! intercepted snow particle follow a gamma distribution with alpha  = 5
                  ! and a mean radius equals to RADIUS_ICESPH
                  ! This approach is used in PBSM3D (see EQ 23 in Pomeroy et al (1993))
                  ! However, it is never described in the papers about the canopy module in CRHM (J. Pomeroy, personal
                  ! communication, 2022)
                  MPM = 4.0 / 3.0 * PI * 917. * RADIUS_ICESPH**3. * (1.0 + 3.0 / ALPHA + 2.0 / ALPHA**2.)

                  ! Sublimation rate coefficient of single 'ideal' ice sphere [s-1]
                  SUB_RATE = (2.0 *PI * RADIUS_ICESPH *SIGMA - SSTAR *J)/ ( LS *J + C1) / MPM

                  ! Compute the intercepted snow exposure coefficient
                  !
                  IF ((SNCMA(I)/SCAP(I)) .LE.  0.0) THEN
                     CE = 0.07
                  ELSE
                    CE = KS * (SNCMA(I)/SCAP(I))**(-1.0*FRACT)
                  ENDIF

                  ! Calculate 'potential' canopy sublimation [s-1]
                  SUB_POT = VGH_DENS(I) * SUB_RATE * CE

                  ! Limit sublimation to canopy snow available and take sublimated snow away from canopy snow at timestep start
                  SUB_CPY = MAX(0.,-SNCMA(I)*SUB_POT*DT)  ! Ensures that only sublimation is computed (neglect solid condensation)

                  !!!! Remove mass of intercepted snow lost by sublimation
                  IF(SUB_CPY .GT. SNCMA(I)) THEN
                     SUB_CPY = SNCMA(I)
                     SNCMA(I) = 0.
                  ELSE
                     SNCMA(I) = SNCMA(I) - SUB_CPY
                  ENDIF

                  ESUBSNC(I) = SUB_CPY/DT
                  SUBSNC_CUM(I) = SUBSNC_CUM(I) + ESUBSNC(I) * DT

               ENDIF

               SNCMA(I) = MAX(0.,SNCMA(I))

               !!!!!!!!
               !! Unloading Code
               !!!!!!!!

               ! Unloading time constant
               IF(TVEG(I) .GE. 273.15) THEN
                  TUNL = TCNM
               ELSE
                  TUNL = TCNC
               ENDIF

               ! Compute unloaded snow mass
               UNLOAD = SNCMA(I) * DT / TUNL

               ! Update intercepted snow mass
               SNCMA(I)  = SNCMA(I) - UNLOAD

               !!!!!!!!
               ! Unload and Melting of snow canopy
               !!!!!!!!
               IF (CANO_REF_FORCING .EQ.'ABV') THEN ! From Boone et al. (2017)

                  HM_CAN(I) = HM_CAN_INI(I) + CPI * SNCMA(I) + CPW *WR_VH(I)

                  IF(TVEG(I) .GT. 273.15) THEN ! Melting of intercepted snow

                     MELT = 5.56E-6 * DT * (SNCMA(I)/SCAP(I)) * (TVEG(I) - 273.15)

                     IF (MELT .GT. SNCMA(I)) THEN
                        MELT = SNCMA(I)
                     ENDIF

                     WR_VH(I) = WR_VH(I) + MELT

                     DRIP_CPY(I) = MAX(0., WR_VH(I) - WRMAX_VH(I) )

                     WR_VH(I) = MIN(WR_VH(I), WRMAX_VH(I))

                     TVEG(I) = TVEG(I) -   CHLF * MELT / HM_CAN(I)

                     SNCMA(I) = SNCMA(I) - MELT

                  ELSE !Refreezing of intercepted liquid water
                     ! Boone et al. (2017) Eq. 83 still uses SNCMA for refreezing of liquid water, so not working if there is no intercepted snow

                     REFREEZE = 5.56E-6 * DT * (WR_VH(I)/WRMAX_VH(I)) * (TVEG(I) - 273.15)  ! < 0
                     
                     IF (ABS(REFREEZE) .GT. WR_VH(I)) THEN
                        REFREEZE = -WR_VH(I)
                     ENDIF
                     WR_VH(I) = WR_VH(I) + REFREEZE 
                     WR_VH(I) = MAX(WR_VH(I), 0.)

                     SNCMA(I) = SNCMA(I) - REFREEZE

                     TVEG(I) = TVEG(I) - CHLF * REFREEZE / HM_CAN(I)
                     DRIP_CPY(I) = 0.    ! Unload considered as solid

                  ENDIF

                  ! 
                  NET_SNOW(I) = DIRECT_SNOW  ! Throughfall only
                  UNLOAD_MASS(I) =  UNLOAD

               ELSE ! Open 2 Forest
                  !!!!!!!
                  !Update snowfall and rainfall rate sent to snow below high-veh
                  !!!!!!!
                  IF(T(I) .GT. 273.15) THEN ! All unload is liquid, first intercepted as liquid intercepted water
                     WR_VH(I) = WR_VH(I) + UNLOAD
                     DRIP_CPY(I) = MAX(0., WR_VH(I) - WRMAX_VH(I) )
                     WR_VH(I) = MIN(WR_VH(I), WRMAX_VH(I))
                     SNW_UNLOAD = 0.
                  ELSE ! All unload is solid + refreezing of liquid intercepted waterto intercepted snow
                     DRIP_CPY(I) = 0.    ! Unload considered as solid
                     SNW_UNLOAD = UNLOAD
                     SNCMA(I) = SNCMA(I) + WR_VH(I)
                     WR_VH(I) = 0.
                  ENDIF

                  !!!! Update the total snow mass sent to the snowpack below
                  !!!! high veg via throughfall and unloading
                  NET_SNOW(I) = DIRECT_SNOW
                  UNLOAD_MASS(I) =  SNW_UNLOAD

               ENDIF

               ! Canopy layer snowcover fractions Boone
               IF (CANO_REF_FORCING .EQ.'ABV') THEN
                  FCANS(I) = (SNCMA(I)/SCAP(I))**0.67
                  ! Update canopy heat mass with new intercepted snow
                  HM_CAN(I) = HM_CAN_INI(I) + CPI * SNCMA(I) + CPW *WR_VH(I)
               ELSE
                  FCANS(I) = 0. ! To avoid computation of sublimation from canopy energy balance in ebudget_svs2_oneprofile_skin
               ENDIF

            ENDIF
         ENDIF
      END DO

      DO I=1,N

         IF (WTG(I,indx_svs2_vh).GE.EPSILON_SVS) THEN

            ! Rain is first intercepted by high veg (rain * (1.-SKYVIEW))
            ! Excess of WRMAX_VH drips + drip from melting intercepted snow + rain*SKYVIEW reach ground

            WR_VH(I) = WR_VH(I) + DT * (1.-SKYVIEW(I)) * RR(I)
            RR_VEG(I) = MAX(0., (WR_VH(I) - WRMAX_VH(I))/DT) + DRIP_CPY(I)/DT + SKYVIEW(I) * RR(I)
                WR_VH(I) = MIN(WR_VH(I), WRMAX_VH(I))

            ! Snowfall rate below high vegetation via throughfall 
            SR_VEG(I) = NET_SNOW(I)/DT

            ! Rate of solid unloading below high-vegetation 
            UNLOAD_RATE(I) = UNLOAD_MASS(I)/DT

         ELSE
            ! no high veg, no modification of rainfall and snowfall

            RR_VEG(I) = RR(I)
            SR_VEG(I) = SR(I)
            UNLOAD_RATE(I) = 0.            

         ENDIF

      ENDDO

      END SUBROUTINE SNOW_INTERCEPTION_SVS2
