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

      SUBROUTINE CANOPY_MET_SVS2 (T, HU,VMOD, ISW, ILW, TVEG, ZU,ZT,          &
                                 SUNCOS,VGH_HEIGHT,VGH_DENS,    &
                                 Z0MVH,Z0SNOW, WTG,LAIVH,SKYVIEW,EMISVH, &
                                 ISW_CAN,ILW_CAN,VMOD_CAN,T_CAN,HU_CAN,   &
                                 VMOD_TOP, PUREF_VEG, PTREF_VEG,  &
                                 PWIND_DRIFT, N )

        use tdpack
        use sfc_options, only: cano_ref_forcing
        use CANOPY_CSTS
        USE MODE_THERMOS 
        use svs_configs
        use svs2_tile_configs
        implicit none

      INTEGER N

      REAL T(N), HU(N), VMOD(N)
      REAL ISW(N), ILW(N),SUNCOS(N) 
      REAL TVEG(N)
      REAL ZU(N), ZT(N),Z0MVH(N),WTG(N,svs2_tilesp1)
      REAL SKYVIEW(N), LAIVH(N),VGH_HEIGHT(N)
      REAL VGH_DENS(N), Z0SNOW(N)
      REAL ISW_CAN(N), ILW_CAN(N), EMISVH(N)
      REAL VMOD_CAN(N),T_CAN(N),HU_CAN(N),PUREF_VEG(N),PTREF_VEG(N),VMOD_TOP(N), PWIND_DRIFT(N)

!
!
!Author
!             V. Vionnet (Mai 2023)
!
!Object
!             Module to simulate the impact of forest of near surface
!             meteorological variables in SVS2:
!                 - radiative fluxes
!                 - wind speed
!                 - air temperature and relative humidity
!
!Arguments
!
!
!          - Input/Output -
!
!          - Input -
! EMISVH    emissivity of high vegetation


      
      INTEGER I,J


      REAL KT ! Sky clearness index
      REAL DFRAC ! Fraction of diffuse radiation 
      REAL DIR_SW,SCA_SW ! Incoming direct and diffuse radiation
      REAL TDIF, TDIR ! Canopy transmittivity for direct and diffuse radiation
      LOGICAL LRAD_SVS ! Use the formulation of SVS for the impact of forest on radiative fluxes 
                       ! Else use the FSM approach

      REAL VMOD_SUB ! Wind speed below canopy at sub-canopy reference height
      REAL VMOD_CLEARING ! Wind speed in clearing at sub-canopy reference height
      REAL DH       ! Displacement height
      REAL ZUREF     ! Height of the wind speed forcing above the ground
      !LOGICAL LFORCING_ABOVE ! True if wind forcing height is given above the canopy  (forcing from NWP system or reanalysis)
                                 ! False is wind forcing height is the absolute height above the ground (mdoel driven by observed forcing) 

      REAL WCAN          ! Canopy wind decay coefficient
      REAL USTAR,  &
           Z0H, & ! Canopy roughness length for heat
           FSURF ! Function used in the calculation of RSURF

     !
     ! 0. Initialize parameters
     !


      DO I=1,N
         

        IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell

           ! Displacement height 
           DH = VGH_HEIGHT(I)*RCHD 

           ZUREF =  VGH_HEIGHT(I) + ZU(I)

           ! Temperature and HU do not change, only the heights do
           T_CAN(I) = T(I)
           HU_CAN(I) = HU(I)

           IF (CANO_REF_FORCING == 'ABV') THEN

              !
              !  xxxx. All forcings are above the canopy and using the surface resistance in the computation of the turbulent fluxes
              !
              !
               PUREF_VEG(I) =  ZUREF
               PTREF_VEG(I) = ZT(I) + VGH_HEIGHT(I)

               VMOD_CAN(I) =  VMOD(I)

               ! Calculation of the wind speed for snow drift under the canopy
               ! Wind speed for snow drift calculation under the canopy at PUREF_VEG

               ! Wind speed at canopy top, assuming logarithmic profile above canopy, dense canopy 
               VMOD_TOP(I) = VMOD(I) * LOG((VGH_HEIGHT(I)-DH)/Z0MVH(I))/LOG((ZUREF-DH)/Z0MVH(I))

               IF (LWIND_FOREST .EQ. 'VDENS_WCAN') THEN
 
                   ! Wind speed at canopy base height, dense canopy, assuming exp profile between canopy top and canopy base height
                   WCAN = WIND_CANO_COEF(LAIVH(I),VGH_DENS(I))  ! From Marke et al., (2016); Liston and Elder (2006)
                   ! Impose a minimum value of 0.2 m/s
                   VMOD_SUB =MAX(0.2,VMOD_TOP(I)*EXP(WCAN*(HSUBCANO/VGH_HEIGHT(I)-1.)))

                   PWIND_DRIFT(I) = VMOD_SUB

                   ! PWIND_DRIFT re-calculated at PUREF (above canopy here) from HSUBCANO (reverse of equations in SNOWFALL_UPGRID() and SNOWDRIFT()
                   PWIND_DRIFT(I) = PWIND_DRIFT(I) *LOG(PUREF_VEG(I)/Z0SNOW(I))/LOG(HSUBCANO/Z0SNOW(I))

               ELSE IF  (LWIND_FOREST  .EQ. 'WEIGHT_AVG') THEN

                   ! Wind speed at canopy base height, dense canopy, assuming exp profile between canopy top and canopy base height
                   WCAN = ZBETA * LAIVH(I)   ! From Marke et al., (2016); Liston and Elder (2006)
                   ! Impose a minimum value of 0.2 m/s                   
                   VMOD_SUB =MAX(0.2,VMOD_TOP(I)*EXP(WCAN*(HSUBCANO/VGH_HEIGHT(I)-1.)))

                   ! PWIND_DRIFT re-calculated at PUREF as weighthed average of wind below the canopy (re-calculated at PUREF) 
                   !and wind in the open as in Mazzotti et al. (2024)
                   PWIND_DRIFT(I) = VGH_DENS(I)**(0.5) *VMOD_SUB*LOG(PUREF_VEG(I)/Z0SNOW(I))/LOG(HSUBCANO/Z0SNOW(I)) + & 
                                    (1.-VGH_DENS(I)**(0.5)) *VMOD(I)
               ENDIF 

           ELSE IF (CANO_REF_FORCING == 'O2F') THEN! Reference height below the canopy

               !
               !  xxxx. Compute impact of forest on wind speed
               !
               !

               PUREF_VEG(I) =  HSUBCANO
               PTREF_VEG(I) = ZT(I)

               ! Wind speed at canopy top, assuming logarithmic profile above canopy, dense canopy
               VMOD_TOP(I) = VMOD(I) * LOG((VGH_HEIGHT(I)-DH)/Z0MVH(I))/LOG((ZUREF-DH)/Z0MVH(I))

             ! Wind speed at canopy base height, dense canopy, assuming exp profile between canopy top and canopy base height
               IF (LWIND_FOREST .EQ. 'VDENS_WCAN') THEN

                    WCAN = WIND_CANO_COEF(LAIVH(I),VGH_DENS(I))
                    VMOD_CAN(I) = MAX(0.2,VMOD_TOP(I)*EXP(WCAN*(HSUBCANO/VGH_HEIGHT(I)-1.)))

               ELSE IF  (LWIND_FOREST .EQ. 'WEIGHT_AVG') THEN

                    WCAN = ZBETA * LAIVH(I) 
                    VMOD_SUB = VMOD_TOP(I)*EXP(WCAN*(HSUBCANO/VGH_HEIGHT(I)-1.))
                    VMOD_CLEARING =VMOD(I)*LOG(PUREF_VEG(I)/Z0SNOW(I))/LOG(ZUREF/Z0SNOW(I)) ! Transfer from ZUREF to PUREF_VEG(I)
                    VMOD_CAN(I) = MAX(0.2,VMOD_SUB * VGH_DENS(I)**(0.5)+(1.-VGH_DENS(I)**(0.5))* VMOD_CLEARING)

               ENDIF


               ! Wind speed for snow drift calculation under the canopy at PUREF_VEG below the canopy
               PWIND_DRIFT(I) = VMOD_CAN(I)

           ELSE  ! Forcing already in forest below canopy
               PUREF_VEG(I) =  ZU(I)
               PTREF_VEG(I) =  ZT(I)
               VMOD_CAN(I) = VMOD(I)
               VMOD_TOP(I) = VMOD(I) ! In that case, should not be used as interception should be off
               T_CAN(I) = T(I)
               HU_CAN(I) = HU(I)

           ENDIF

         ELSE  
             ! No high vegetation present 
             PUREF_VEG(I) =  ZU(I)
             PTREF_VEG(I) =  ZT(I)
             VMOD_CAN(I) = VMOD(I)
             VMOD_TOP(I) = VMOD(I)
             T_CAN(I) = T(I)
             HU_CAN(I) = HU(I)

         ENDIF    

      ENDDO

      ! Radiation under canopy
      DO I=1,N
        IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell

          IF (CANO_REF_FORCING == 'FOR') THEN  ! Forcing already in forest below canopy
               ILW_CAN(I) = ILW(I)
               ISW_CAN(I) = ISW(I)
          ELSE ! 'O2F' or 'ABV', radiation should be determined under canopy

     !
     !  xxxx. Compute impact of forest on incoming radiative fluxes 
     !


              ! Compute sky clearness
              ! SUNCOS: cos of zenithal angle (max when sun is at its
              ! highest above the horizon).

              IF(SUNCOS(I)>0) THEN
                  KT = ISW(I)/ (1367.*SUNCOS(I))
              ELSE
                  KT = 0.
              ENDIF

              ! Compute fraction of diffuse radiation 
              IF(KT>0.8) THEN
                   DFRAC  = 0.165
              ELSE IF(KT>0.22) THEN
                   DFRAC = 0.95 - 0.16*KT + 4.39*KT**2. - 16.64*KT**3. +12.34*KT**4. 
              ELSE
                    DFRAC = 1 - 0.09*KT
              ENDIF

              ! Compute solar diffuse and direct radiation
              SCA_SW = DFRAC*ISW(I)  ! Diffuse incoming SW
              DIR_SW = (1.-DFRAC)*ISW(I) ! Direct incoming SW

              ! Compute transmissivity of diffuse radiation 
              ! Need to double check if VEGDENS should be employed here
              ! to derive an effective LAI. 
              TDIF = SKYVIEW(I)

              ! Compute transmissivity of direct radiation 
              IF(SUNCOS(I)>0) THEN
                 TDIR = EXP(-KEXT*VGH_DENS(I)*LAIVH(I)*CLUMPING/SUNCOS(I))
              ELSE
                 TDIR = TDIF
              ENDIF

              ! Subcanopy total shortwave radiation
              ISW_CAN(I) = TDIF * SCA_SW + TDIR *DIR_SW

              IF (CANO_REF_FORCING == 'ABV') THEN
                  ! Subcanopy longwave radiation with skin canopy temperature
                  ILW_CAN(I) = TDIF * ILW(I) + (1.-TDIF)*STEFAN *EMISVH(I)*TVEG(I)**4
              ELSE IF (CANO_REF_FORCING == 'O2F') THEN
               ! Subcanopy longwave radiation with canopy temperature
               ! taken as forcing air temperature as a proxy
                  ILW_CAN(I) = TDIF * ILW(I) + (1.-TDIF)*STEFAN *T_CAN(I)**4
              ENDIF

          ENDIF

        ELSE
          ILW_CAN(I) = ILW(I)
          ISW_CAN(I) = ISW(I)
        ENDIF

      ENDDO


      END SUBROUTINE CANOPY_MET_SVS2
