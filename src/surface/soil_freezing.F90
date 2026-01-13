!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%

      SUBROUTINE SOIL_FREEZING(DT, TSOIL, VEGL, VEGH, PSN, PSNVH,  &
                                SOILCONDZ, SOILHCAPZ , TGRS, TVEGS,   &
                                WSOIL, ISOIL,  &
                                SNORO, SNODP, TSNO, TSKIN_BG, &
                                SNVRO, SNVDP, TSNV, TSKIN_VEG, &
                                TDEEP, WUNFRZ, WSAT, &
                                DWATERDT_SURF,DWATERDT_DEEP,  N)


      USE TDPACK
      USE SFC_OPTIONS
      USE SVS_CONFIGS

      IMPLICIT NONE

      ! Input
      INTEGER N ! Number of grid points

      REAL DT

      REAL, DIMENSION(N)        :: VEGL, VEGH, PSN, PSNVH, TGRS, TDEEP,TVEGS
      REAL, DIMENSION(N)        :: SNORO, SNODP, TSNO, TSKIN_BG 
      REAL, DIMENSION(N)        :: SNVRO, SNVDP, TSNV, TSKIN_VEG 
      REAL, DIMENSION(N)        :: DWATERDT_SURF,DWATERDT_DEEP
      REAL, DIMENSION(N,NL_SVS) :: TSOIL,SOILCONDZ, SOILHCAPZ,WSOIL,ISOIL, WUNFRZ, WSAT

      !
      !Author
      !          V. Vionnet, V. Fortin, K. Rasouli (April 2020)
      !Revisions
      !
      !Object
      ! Simulate the evolution of soil freezing and thawing 
      ! using the simple heat conduction method proposed by
      ! Hayashi et al. (2007) and Mohammed et al. (2012)
      ! The upper boundary conditions are taken from the FR schemes and
      ! the coupling is done as the IFS code (see IFS techincal
      ! documentation)
      
      !
      !Arguments
      !
      !          - INPUT -
      !
      ! DT       timestep
      
      !          --- (Surface) Cover Fraction  ---
      !
      ! VEGL           fraction of LOW vegetation [0-1]
      ! VEGH           fraction of HIGH vegetation [0-1]
      ! PSN            fraction of bare ground or low veg. covered by snow [0-1]
      ! PSNVH          fraction of ground covered by snow below high veg.  [0-1]
      ! 
      !          ---  Soil thermal properties   ---
      !
      ! SOILCONDZ      soil thermal conductivity (per layer) [W K-1 m-1]
      ! SOILHCAPZ      soil heat capacity (per layer) [J m-3 K-1]
      ! TDEEP          constant deep soil temperature [K]
      ! WUNFRZ         unfrozen residual water content [m3/m3]
      ! WSAT           Saturated hydraulic conductivity or porosity [m3/m3]

      !          --- Prognostic variables of SVS not modified by SOIL_FREEZING ---
      !
      ! TGRS          bare ground surface temperature from Force Restore
      ! TVEGS         surface vegetation temperature from Force Restore
      ! SNODP         snow depth for snow over bare ground/low veg
      ! SNORO         snow density for snow over bare ground/low veg
      ! TSKIN_BG      skin snow temperature over bare ground 
      ! TSNO          deep snow temperature for snow over bare ground/low veg
      ! SNVDP         snow depth for snow over under high veg
      ! SNVRO         snow density for snow under high veg
      ! TSNV          deep snow temperature for snow under high veg
      ! TSKIN_VEG     skin snow temperature under high vegetation 
      !
      !          - INPUT/OUTPUT  -
      !
      !          --- Prognostic variables of SVS modified by SOIL_FREEZING ---
      !
      ! TSOIL          Soil temperature (per layer) [K]
      ! WSOL (NL_SVS)    soil volumetric water content (per layer) [m3/m3]
      ! ISOL (NL_SVS)    frozen soil volumetric water (per layer) [m3/m3]
      !
      !          - OUTPUT  -
      ! DWATERDT_SURF  net tendency of melting-freezing of soil water for the
      !                surface layer (for ebudget_svs) [kg/m2/s]
      ! DWATERDT_DEEP  net tendency of melting-freezing of soil water for the
      !                deep layer of the FR scheme (for ebudget_svs) [kg/m2/s]
      !
      !          -  DIMENSIONS  -
      !
      ! N              number of grid cells

      ! Local Variable and arrays
      INTEGER I, K


      INTEGER OPT_FRAC    ! Option to compute the snow cover fraction        
      INTEGER OPT_LIQWAT  ! Option to compute the unfrozen redisudal water content  
      INTEGER OPT_VEGCOND ! Option to compute the skin conductivity from the snow-free vegetation
      
      REAL LAM_VEGL_STAB, LAM_VEGH_STAB,LAM_VEGL_UNSTAB, LAM_VEGH_UNSTAB

      REAL HNET,HNETR,TTEST, TTEST2, UFWC,DFWC, FWCTEST, QLAT
      REAL RTH_GRND, RTH_SNO,RTH_SNV,FAC_SNW
      REAL LAM_GRND
      REAL CHI
      REAL, DIMENSION(N) :: HFLUX_GRND, HFLUX_SNO,HFLUX_SNV, HFLUX_VEG
      REAL, DIMENSION(N, NL_SVS+1) :: RTH, HFLUX
      REAL, DIMENSION(N, NL_SVS) :: WC, RFS, ISOILT, TSOILT
      REAL, DIMENSION(NL_SVS)   :: ZLAYER, WSURF,WDEEP!, LHEAT_RELEASE
      LOGICAL LHEAT_RELEASE

      REAL, DIMENSION(N)           :: TBTM, DBTM, LAMS, LAMSV, FRAC_SNWL, FRAC_SNWH
      REAL KDIFFU,KDIFFUV 
      REAL, DIMENSION(N)           :: DAMPD,DAMPDV,LAM_VEG
      REAL, DIMENSION(N,SVS_TILESP1)           :: WTG
      REAL, DIMENSION(N)           :: DHEAT


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   0. Initialize bottom temperature and depth of the bottom layer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !

      OPT_FRAC = 2 ! Option to compute the snow cover fraction
                   ! 1: use a fraction = SWE/1 mm
                   ! 2: use the formulation of Niu and Yang (2007) (Recommended) 

      OPT_LIQWAT = 2 ! Option to handle unfrozen liquid water content
                     ! 1: use a constant value of 0.06 (as in VSMB)
                     ! 2: use a value that depends on soil texture (Recommended for large scale simulations) 

      OPT_VEGCOND = 1 ! Option to compute the skin conductivity for the vegetation 
                      ! 1: same value of 10 for low and high veg. in both stable and unstable conditions (Most recent version of EC Land, Boussetta et al., 2021; Recommended)
                      ! 2: different values for low and high veg. in stable and unstable conditions (See Trigo et al., JGR, 2015)

                   
      IF(soilsnowhf_svs1 .EQ. 'DST_HD') THEN
              FAC_SNW = 0.5
      ELSE IF(soilsnowhf_svs1 .EQ. 'DST_FD') THEN
              FAC_SNW = 1.0
      ENDIF                   

      ! Define skin conductivity for low and high vegetation (W m-2 K-1)
      IF(OPT_VEGCOND==1) THEN 
          ! Effect of stable and unstable stratification are not taken into account 
          ! Values taken from Tab 1.2 in sup. material of Boussetta et al (2021)
          LAM_VEGH_STAB = 10 
          LAM_VEGL_STAB = 10  
          LAM_VEGH_UNSTAB = 10 
          LAM_VEGL_UNSTAB = 10  
      ELSE IF(OPT_VEGCOND ==2) THEN
          ! Effect of stable and unstable stratification are  taken into account
          ! Values taken from Tab 3 in Trigo et al. (2015)
          LAM_VEGH_STAB = 15
          LAM_VEGL_STAB = 10  
          LAM_VEGH_UNSTAB = 20 
          LAM_VEGL_UNSTAB = 10  
      ENDIF

      IF(soilgrndhf_svs1 .EQ. 'LAM_BOU2021') THEN
          LAM_GRND = 15
      ENDIF

      ! Compute layer depth
      ZLAYER(1) =  DELZ(1)
      DO K =2, NL_SVS
        ZLAYER(K) = ZLAYER(K-1) + DELZ(K)
      ENDDO

      ! Compute weight for the calculation of the net tendency of melting-freezing 
      ! For the suface and the deep layer
      DO K =1, NL_SVS
        IF(ZLAYER(K) .LE. HSURF) THEN
           WSURF(K) = 1.0
        ELSE IF( ZLAYER(K)> HSURF ) THEN
           IF(K==1) THEN
              WSURF(K) = HSURF/ZLAYER(K)
           ELSE IF(ZLAYER(K-1)<=HSURF) THEN
              WSURF(K) = (HSURF-ZLAYER(K-1))/(ZLAYER(K)-ZLAYER(K-1))
           ELSE
              WSURF(K) = 0.
           ENDIF
        ENDIF

        IF(ZLAYER(K) .LE. HDEEP) THEN
           WDEEP(K) = 1.0
        ELSE IF( ZLAYER(K)> HDEEP ) THEN
           IF(K==1) THEN
              WDEEP(K) = HDEEP/ZLAYER(K)
           ELSE IF(ZLAYER(K-1)<=HDEEP) THEN
              WDEEP(K) = (HDEEP-ZLAYER(K-1))/(ZLAYER(K)-ZLAYER(K-1))
           ELSE
              WDEEP(K) = 0.
           ENDIF
        ENDIF
        
      ENDDO

      
      DO  I=1,N
        TBTM(I) = TDEEP(I) ! K
        IF (soildbtm_svs1 .EQ. 'MID' .OR. soildbtm_svs1 .EQ. 'NOFL') THEN
            DBTM(I) = ZLAYER(NL_SVS) + 0.5* DELZ(NL_SVS) ! m
        ELSE IF (soildbtm_svs1 .EQ. 'DEEP') THEN
            IF (ZLAYER(NL_SVS) < 3.0) THEN                                        ! If the soil column is thinner than 3 m, set DBTM to 8.5 m (min condition)
                DBTM(I) = 7.5 ! m
            ELSE IF (ZLAYER(NL_SVS) >= 5.0 .AND. ZLAYER(NL_SVS) < 12.5) THEN      ! if the soil column thickness is between 5 and 12.5 m, set DBTM to 20 m (max condition)
                DBTM(I) = 12.5 ! m
            ELSE IF (ZLAYER(NL_SVS) >= 12.5) THEN                                 ! if the soil column is thicker than 12.5 m, use the default parameterization to compute DBTM (to avoid computational errors)
                DBTM(I) = ZLAYER(NL_SVS) + 0.5* DELZ(NL_SVS) ! m
            ELSE                                                                  ! if soil layer thickness is between 3 and 5 m, DBTM corresponds to 2.5 times the thickness of the soil column.
                DBTM(I) = ZLAYER(NL_SVS)*2.5 ! m
            ENDIF
        ENDIF

        ! Initialize values for net tendency of thawing-freezing of soil water
        DWATERDT_SURF(I) = 0
        DWATERDT_DEEP(I) = 0.

         ! Snow thermal conductitivy
        LAMS(I) = LAMI * SNORO(I)**1.88
        LAMSV(I) = LAMI * SNVRO(I)**1.88

        ! Vegetation average skin conductivity
        IF(VEGL(I)+VEGH(I)> 0.) THEN
            IF(TVEGS(I) > TSOIL(I,1)) THEN ! Stable case
                 LAM_VEG(I) =(VEGL(I)*LAM_VEGL_STAB+VEGH(I)*LAM_VEGH_STAB)/(VEGL(I)+VEGH(I))
            ELSE  ! Unstable case
                 LAM_VEG(I) =(VEGL(I)*LAM_VEGL_UNSTAB+VEGH(I)*LAM_VEGH_UNSTAB)/(VEGL(I)+VEGH(I))
            ENDIF
        ELSE
            LAM_VEG(I) = 1. ! Set default value to make sure code is running           
        ENDIF

        ! Snow cover fraction used for the exchanges with the surface
        IF(OPT_FRAC==1) THEN
              FRAC_SNWL(I) =  MIN(RAUW*SNORO(I)*SNODP(I)/1.0,1.0)
              FRAC_SNWH(I) =  MIN(RAUW*SNVRO(I)*SNVDP(I)/1.0,1.0)
        ELSE
             !
             ! Use the approach of Niu and Yang (2007) by default in the soil freezing scheme
             FRAC_SNWL(I) = 0.
             FRAC_SNWH(I) = 0.
             IF(  SNODP(I)>0.) THEN
                FRAC_SNWL(I) = TANH(SNODP(I)/(2.5*Z0_NIU*(RAUW*SNORO(I)/RHONEW)**MFAC))
             ENDIF
             IF( SNVDP(I)>0.) THEN
                FRAC_SNWH(I) = TANH(SNVDP(I)/(2.5*Z0_NIU*(RAUW*SNVRO(I)/RHONEW)**MFAC))
             ENDIF
             !
        ENDIF 

        DO K =1, NL_SVS
            IF(OPT_LIQWAT==1) THEN
                  RFS(I,K) = 0.06 ! Residual unfrozen content
            ELSE
                  RFS(I,K) = WUNFRZ(I,K) ! Residual unfrozen content
            ENDIF
            WC(I,K) = WSOIL(I,K) + ISOIL(I,K) ! Total water content
            ISOILT(I,K) = ISOIL(I,K)  ! Initialize frozen soil volumetric water (to be updated in the routine) 
            TSOILT(I,K) = TSOIL(I,K)  ! Initialize soil temperature 
        ENDDO

        IF(soilsnowhf_svs1 .EQ. 'DST_MAXD' .OR. soilsnowhf_svs1 .EQ. 'ST_D_DD') THEN !Calculation of damping depth
! 
             IF( SNODP(I)>0.) THEN

                KDIFFU =  LAMS(I) / ( CICE * RAUW*SNORO(I)) ! Thermal diffusivity (You et al., 2014)
                DAMPD(I) = SQRT(  2.0 * KDIFFU / MYOMEGA  ) !  Damping depth in m , assuming diurnal forcing dominates
             ELSE
                DAMPD(I) = 0.
             ENDIF

             IF( SNVDP(I)>0.) THEN

                KDIFFUV =  LAMSV(I) / ( CICE * RAUW*SNVRO(I)) ! Thermal diffusivity (You et al., 2014)
                DAMPDV(I) = SQRT( 2.0 * KDIFFUV / MYOMEGA  ) !  Damping depth in m , assuming diurnal forcing dominates
             ELSE
                DAMPDV(I) = 0.
             ENDIF

        ENDIF        
      ENDDO

      ! Compute weights of surface type in SVS:
      !               - snow-free bare ground WTG_2
      !               - snow-free vegetation WTG_3
      !               - snow-covered bare ground and low vegetation WTG_4
      !               - snow below high-vegetation WTG_5

      CALL WEIGHTS_SVS(VEGH, VEGL,FRAC_SNWL,FRAC_SNWH,N,WTG)

      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     1. Compute the thermal resistances and the heat flux between
      !        adjacent layers
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      DO  I=1,N
        !
        ! Treatment of surface layer!
        !

        ! Upper boundary condition for snow-free bare ground
        IF (soilgrndhf_svs1 .EQ. 'RTH_GRND') THEN   
            RTH_GRND = 0.5*DELZ(1)/SOILCONDZ(I,1)
            HFLUX_GRND(I) = (TGRS(I) - TSOIL(I,1)) / RTH_GRND
        ELSE IF (soilgrndhf_svs1 .EQ. 'LAM_BOU2021') THEN
            HFLUX_GRND(I) = (TGRS(I) - TSOIL(I,1)) * LAM_GRND
        ENDIF

        ! Upper boundary condition for snow-free vegetation
        HFLUX_VEG(I) = (TVEGS(I) - TSOIL(I,1)) * LAM_VEG(I)              

        ! Upper boundary condition for snow over bare ground and low veg.
        IF(SNODP(I) > 0.) THEN ! Snow is present
            IF (soilsnowhf_svs1 .EQ. 'DST_HD' .OR. soilsnowhf_svs1 .EQ. 'DST_FD') THEN
                RTH_SNO = FAC_SNW*SNODP(I)/LAMS(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
            ELSE IF (soilsnowhf_svs1 .EQ. 'DST_MAXD') THEN
                RTH_SNO = MAX(SNODP(I)/2., SNODP(I)-DAMPD(I))/LAMS(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
            ELSE
                IF (SNODP(I) > DAMPD(I)) THEN
                    RTH_SNO = MAX(SNODP(I)/2., SNODP(I)-DAMPD(I))/LAMS(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
                ELSE
                    RTH_SNO = SNODP(I)/LAMS(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
                ENDIF
            ENDIF
        
            ! Heat flux at the snow/soil interface
            IF (soilsnowhf_svs1 .EQ. 'ST_D_DD' .AND. SNODP(I) <= DAMPD(I)) THEN
                HFLUX_SNO(I) = (TSKIN_BG(I) - TSOIL(I,1)) / RTH_SNO
            ELSE
                HFLUX_SNO(I) = (TSNO(I) - TSOIL(I,1)) / RTH_SNO
            ENDIF
        ELSE ! No snow 
            HFLUX_SNO(I) = 0.
        ENDIF

                 
        ! Upper boundary condition for snow below high vegetation. 
        IF(SNVDP(I) > 0.) THEN ! Snow is present
                
            IF (soilsnowhf_svs1 .EQ. 'DST_HD' .OR. soilsnowhf_svs1 .EQ. 'DST_FD') THEN
                RTH_SNV = FAC_SNW*SNVDP(I)/LAMSV(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
            ELSE IF (soilsnowhf_svs1 .EQ. 'DST_MAXD') THEN 
                RTH_SNV = MAX(SNVDP(I)/2., SNVDP(I)-DAMPDV(I))/LAMSV(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
            ELSE
                IF (SNVDP(I) > DAMPDV(I)) THEN
                    RTH_SNV = MAX(SNVDP(I)/2., SNVDP(I)-DAMPDV(I))/LAMSV(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
                ELSE
                    RTH_SNV = SNVDP(I)/LAMSV(I) + 0.5*DELZ(1)/SOILCONDZ(I,1)
                ENDIF
            ENDIF
!
!            ! Heat flux at the snow/soil interface
            IF (soilsnowhf_svs1 .EQ. 'ST_D_DD' .AND. SNVDP(I) <= DAMPDV(I)) THEN
                HFLUX_SNV(I) = (TSKIN_VEG(I) - TSOIL(I,1)) / RTH_SNV
            ELSE
                HFLUX_SNV(I) = (TSNV(I) - TSOIL(I,1)) / RTH_SNV
            ENDIF
!
        ELSE ! No snow 
           HFLUX_SNV(I) = 0.
        ENDIF
                

        ! Compute average surface heat flux using weights for each surface tile
        HFLUX(I,1) = WTG(I,2) * HFLUX_GRND(I) + WTG(I,3) * HFLUX_VEG(I) + &
                               WTG(I,4) * HFLUX_SNO(I) + WTG(I,5) * HFLUX_SNV(I)
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
        !
        ! Treatment of the following layers
        !
        DO K =2, NL_SVS
             RTH(I,K) = 0.5*DELZ(K-1)/SOILCONDZ(I,K-1) + 0.5*DELZ(K)/SOILCONDZ(I,K)
             HFLUX(I,K) = (TSOIL(I,K-1) - TSOIL(I,K))/ RTH(I,K)
        ENDDO
        !
        ! Treatment of the bottom layer
        ! Use thermal conductivity of the deepest SVS layer
        !
        RTH(I,NL_SVS+1) = 0.5* DELZ(NL_SVS)/ SOILCONDZ(I,NL_SVS) + (DBTM(I) - ZLAYER(NL_SVS)) / SOILCONDZ(I,NL_SVS)

        IF (soildbtm_svs1 .EQ. 'MID' .OR. soildbtm_svs1 .EQ. 'DEEP') THEN                                       !Flux below the soil column based on DBTM.
            HFLUX(I,NL_SVS+1) = ( TSOIL(I,NL_SVS) - TBTM(I)) / RTH(I,NL_SVS+1)
        ELSE IF (soildbtm_svs1 .EQ. 'NOFL') THEN                                                     !Zero flux condition
            HFLUX(I,NL_SVS+1) = 0
        ENDIF

      ENDDO
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !     2. Compute the evolution of soil temperature
      !        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !
      ! 

      DO I=1, N
        DO K =1, NL_SVS

          ! Set the flag for latent heat relase due to freezing/thawing
          ! to false (default value) 
           LHEAT_RELEASE = .FALSE.

           HNET = (HFLUX(I,K)- HFLUX(I,K+1))*DT ! Heat flux received by layer K

           !Compute the efficiency factor for freezing and thawing
           IF (lphase_change_eff_svs1) THEN ! From surfex (Adjusted)
               IF (HNET .LT. 0.0) THEN
                   CHI = MIN(MAX((WSOIL(I,K)-RFS(I,K))/WSAT(I,K),CHI_MIN),1.0) 
               ELSE
                   CHI = MIN(MAX(ISOIL(I,K)/(WSAT(I,K)-RFS(I,K)),CHI_MIN),1.0)
               ENDIF
           ELSE
               CHI = 1.0
           ENDIF
           
           IF(TSOILT(I,K) - TRPL .GT. EPSILON_SVS_TK) THEN
              !TSOIL POSITIVE
                TTEST = TSOILT(I,K) + HNET/(SOILHCAPZ(I,K)*DELZ(K))
                IF(TTEST .LT. TRPL) THEN
                     UFWC = MAX(WSOIL(I,K) - RFS(I,K), 0.) !Maximum liquid water available for freezing
                     IF(UFWC>0.) THEN 
                        ! if have unfrozen water available for freezing
                        HNETR = HNET + (TSOILT(I,K)-TRPL) * SOILHCAPZ(I,K)*DELZ(K)

                        ! Maximum ice content that could be potentially formed
                        DFWC = -1.0*CHI*HNETR/(RAUW*CHLF*DELZ(K)) ! Maximum ice content that could be potentially formed
                        IF(UFWC>DFWC) THEN  !  Enough liquid water for freezing, temperature stay constant
                           ! All energy will be used to freeze water
                           ! because max created ice < max liquid water that can be frozen
                           TSOILT(I,K) = TRPL
                           ISOILT(I,K) = DFWC + ISOIL(I,K)
                        ELSE ! All available liquid water is frozen and temperature keep decreasing
                           ! Freeze all available water, and remaining energy flux will decrease temperature
                           !Remaining energy for temperature change (cooling)
                           HNETR  = HNETR +UFWC* RAUW*CHLF*DELZ(K)/CHI
                           TSOILT(I,K) =  TRPL + HNETR/(SOILHCAPZ(I,K)*DELZ(K))
                           ISOILT(I,K) = UFWC + ISOIL(I,K)
                        ENDIF

                        ! Latent heat is released due to freezing
                        LHEAT_RELEASE = .TRUE.

                     ELSE
                        TSOILT(I,K) = TTEST ! No enough liquid water for freezing, temperature keep decreasing. 
                     ENDIF                     
               ELSE  
                     TSOILT(I,K) = TTEST
               ENDIF

            ELSE IF( abs(TSOILT(I,K)-TRPL) .LE. EPSILON_SVS_TK) THEN

               ! TSOIL within "epsilon" of TRPL
               DFWC = -1.0 * CHI * HNET/(RAUW*CHLF*DELZ(K))
               UFWC = MAX(WSOIL(I,K) - RFS(I,K) , 0.)
               
               FWCTEST = ISOIL(I,K) + DFWC
               IF(FWCTEST.LE. 0.0) THEN 
                  ! Total melting of frozen content and ground heating 
                  ! with the remaining energy
                  HNETR = HNET -ISOIL(I,K) * RAUW*CHLF*DELZ(K)/CHI
                  ISOILT(I,K) = 0.0    
                  TSOILT(I,K) = TSOILT(I,K) + HNETR/(SOILHCAPZ(I,K)*DELZ(K))
               ELSE 
                  IF(DFWC.GT.UFWC) THEN
                      !Total freezing of soil layer and ground cooling
                      ! with the remaining energy
                      HNETR = HNET + UFWC * RAUW*CHLF*DELZ(K)/CHI
                      ISOILT(I,K) = ISOIL(I,K) + UFWC
                      TSOILT(I,K) = TSOILT(I,K) + HNETR/(SOILHCAPZ(I,K)*DELZ(K))
                   ELSE
                      ! layer is still partially frozen and T = 0 deg
                      ISOILT(I,K) = FWCTEST
                      TSOILT(I,K) = TRPL
                   ENDIF
                ENDIF

                ! Latent heat is released due to freezing
                LHEAT_RELEASE = .TRUE.
                

             ELSE  ! Soil at negative temperature

                ! Temperature that would be reached without phase change
                TTEST = TSOILT(I,K) + HNET/(SOILHCAPZ(I,K)*DELZ(K)) 

                IF(TTEST .GT. TRPL) THEN
                     !
                     ! Enough energy is brought to heat the soil to 0 deg and 
                     ! melt part of the ice if ice is present
                     !
                     IF(ISOIL(I,K)>0.) THEN 
                        ! If ice is present, compute the energy left after warming the soil
                        ! temp. to 0 degC
                        HNETR = HNET + (TSOILT(I,K)-TRPL) * SOILHCAPZ(I,K)*DELZ(K)

                        ! Maximum ice content that could be potentially melted with such amount of energy
                        DFWC = CHI * HNETR/(RAUW*CHLF*DELZ(K))

                        IF(DFWC<ISOIL(I,K)) THEN 
                            ! All energy is used to melt ice and some ice remains 
                            TSOILT(I,K) = TRPL
                            ISOILT(I,K) = ISOIL(I,K)-DFWC
                        ELSE         
                            ! All ice is melted and remaining energy is used to warm the layer above 0 degC
                            ! Remove the energy required to melt the ice 
                            HNETR  = HNETR -ISOIL(I,K)* RAUW*CHLF*DELZ(K)/CHI 
                            ! Update the temperature and the ice content
                            TSOILT(I,K) =  TRPL + HNETR/(SOILHCAPZ(I,K)*DELZ(K))
                            ISOILT(I,K) = 0. 
                        ENDIF

                        ! Latent heat is released due to melting
                        LHEAT_RELEASE = .TRUE.

                     ELSE
                        TSOILT(I,K) = TTEST ! No melting of ice and no contribution from phase change  
                     ENDIF

                ELSE
                     !
                     ! The temperature remains negative 
                     ! 
                     ! Check if liquid water is present above the residual unfrozen water content
                     ! and freeze this water if it is the case. Such situation should not be encoutered at model runtime
                     ! but may be present if data assimilation has changed in an unconsistent way the soil temperature and/or water content. 
                     UFWC = MAX(WSOIL(I,K) - RFS(I,K), 0.)

                     IF(UFWC .GT. 0.) THEN
                         ! There is liquid water that can be frozen.
                         ! Compute the energy that would be released by the freezing of this amount of water 
                         QLAT  = UFWC* RAUW*CHLF*DELZ(K)/CHI
                         ! Compute the temperature that would be reached
                         TTEST2 = TTEST + QLAT/(SOILHCAPZ(I,K)*DELZ(K))

                         IF(TTEST2 .GT. TRPL) THEN 
                              ! Too much energy would be released 
                              QLAT =  (TRPL-TTEST)*SOILHCAPZ(I,K)*DELZ(K) ! Compute the energy which is actually relasead
                              ISOILT(I,K)  = ISOIL(I,K) + CHI * QLAT/(RAUW*CHLF*DELZ(K)) ! Update ice content
                              TSOILT(I,K) =TRPL
                         ELSE
                              ! All the available liquid water is melting and the temperature reamins below 0 deg 
                              ISOILT(I,K) = ISOIL(I,K) + UFWC 
                              TSOILT(I,K) = TTEST2
                         ENDIF

                        ! Latent heat is released due to freezing
                        LHEAT_RELEASE = .TRUE.

                     ELSE
                          ! No liquid water is available for freezing
                          ! The ice content does not change and the temperature remains negative
                          TSOILT(I,K) = TTEST  
                     ENDIF

                  ENDIF

             ENDIF

             IF(LHEAT_RELEASE) THEN
                 DHEAT(I) = HNET - (TSOILT(I,K)- TSOIL(I,K) ) *SOILHCAPZ(I,K)*DELZ(K)
                 DWATERDT_SURF(I) = DWATERDT_SURF(I)-1.0*CHI*WSURF(K)*DHEAT(I)/(CHLF*DT)
                 DWATERDT_DEEP(I) = DWATERDT_DEEP(I) -1.0*CHI*WDEEP(K)*DHEAT(I)/(CHLF*DT)
             ENDIF

            ! Update soil temperature
            TSOIL(I,K) = TSOILT(I,K)
             
            ! Update frozen water content
            ISOIL(I,K) = ISOILT(I,K)
    
            ! Update liquid water content
            WSOIL(I,K) = WC(I,K) - ISOIL(I,K)
       
        ENDDO
      ENDDO

      ! Final adjustement to make sure that DWATERDT_DEEP does not include the contribution from DWATERDT_SURF
      DWATERDT_DEEP(:) = DWATERDT_DEEP(:) -DWATERDT_SURF(:) 


!

      END SUBROUTINE SOIL_FREEZING
