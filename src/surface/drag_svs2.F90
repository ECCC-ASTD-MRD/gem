!-------------------------------------- LICENCE BEGIN -------------------------
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

      SUBROUTINE DRAG_SVS2 ( TGRS, TGRVS, TVGLS, TVGHS, TSV,TS,  WD1, &
                              WR_VL, WR_VH,  THETAA, VMOD, VDIR, HU, RHOA, &
                              PS, RS, Z0, Z0LOC, Z0VG, WFC, WSAT, CLAY1,  &
                              SAND1, LAI_VL, LAI_VH, WRMAX_VL,WRMAX_VH,&
                              ZUSL, ZTSL, LAT, PSNVH, &
                              FCOR, Z0HA, WTG, &
                              VGH_DENS, Z0MVH,Z0MVL, Z0SNOW, Z0HSN,VGH_HEIGHT,  &
                              LAIVH, ZVCAN, FCANS,SNCMA, &
                              RESAGR,RESAGRV, RESA_VL, RESA_VH, RES_SNCA, RESA_SV, &
                              HUSURF,HUSURFGV, &
                              HRSURF,HRSURFGV, &
                              HV_VL, HV_VH, HVSN_VH, DEL_VL, DEL_VH,  &
                              Z0HBG, Z0HVL, Z0HVH,Z0HGV,  N)
      use tdpack
      use sfc_options
      use sfclayer, only: sl_sfclayer,SL_OK
      use svs_configs
      use MODE_THERMOS
      use MODD_CSTS
      use svs2_tile_configs
      use canopy_csts
!
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N
      REAL TGRS(N), TVGLS(N), TVGHS(N), WR_VL(N), WR_VH(N), THETAA(N), VMOD(N), VDIR(N), HU(N), CLAY1(N)
      REAL SAND1(N), PS(N), RS(N), Z0(N),Z0LOC(N), Z0VG(N), WFC(N,NL_SVS), WSAT(N,NL_SVS)
      REAL LAI_VL(N), LAI_VH(N),WRMAX_VH(N),WRMAX_VL(N), ZUSL(N), ZTSL(N), LAT(N)
      REAL FCOR(N), Z0HA(N), Z0HBG(N), Z0HVL(N), Z0HVH(N), Z0HGV(N), FCANS(N)
      REAL VGH_DENS(N), Z0MVH(N), Z0MVL(N), VGH_HEIGHT(N), LAIVH(N), ZVCAN(N)
      REAL RESAGR(N),RESAGRV(N), RESA_VL(N), RESA_VH(N), RES_SNCA(N), RESA_SV(N)
      REAL HUSURF(N),HUSURFGV(N), HV_VL(N), HV_VH(N), HVSN_VH(N), DEL_VL(N), DEL_VH(N)
      REAL HRSURF(N),HRSURFGV(N), WD1(N), Z0SNOW(N), Z0HSN(N), PSNVH(N)
      REAL WTG(N,svs2_tilesp1), TGRVS(N), SNCMA(N), RHOA(N)
      REAL TSV(N,NSL),TS(N,NSL)
!
!Author
!          S. Belair, M.Abrahamowicz, S.Z.Husain, N.Alavi, S.Zhang (June 2015) 
!Revisions
! 001      Name (date) - Comment
!
!
!Object
!
!     Calculates the drag coefficients for heat and momentum transfers
!     over ground and vegetation(i.e., Ch and Cd).
!
!
!Method
!
!
!     1) computes hu, hv, and DEL
!
!     2) use this to find qsoil, the grid-averaged relative humidity
!        of the soil
!
!     3) find the transfer and resistance coefficients Ch, Cd, and Ra
!        Calculate the surface fluxes of heat, moisture,
!        and momentum over water surfaces.
!
!Arguments
!
!          - Input/Output -
! RESAGR    aerodynamical surface resistance for bare ground
! RESAGRV    aerodynamical surface resistance for bare ground
! RESA_VL    aerodynamical surface resistance for low vegetation
! RESA_VH    aerodynamical surface resistance for high vegetation
!
!          - Input -
! TGRS      skin (surface) temperature of bare ground
! TGRVS      skin (surface) temperature of ground below high veg.
! TVGLS      skin (surface) temperature of low vegetation
! TVGHS      skin (surface) temperature of high vegetation
! TSV       Temperature of the snow layers in the high vegetation
! WD1       Soil volumetric water content (first level)
! WR_VL     water content retained by the low vegetation canopy
! WR_VH     water content retained by the high vegetation canopy
! THETAA    potential temperature at the lowest level
! VMOD      wind speed at the lowest level
! VDIR      wind direction at the lowest level
! HU        specific humidity of air at the lowest level
! PS        surface pressure
! RS        surface or stomatal resistance
! Z0        momentum roughness length (no snow)
! Z0LOC     local (no orography) land momentum roughness
! Z0VG      averaged vegetation-only momentum roughness
! WFC       volumetric water content at the field capacity
! WSAT      volumetric water content at saturation
! CLAY1     percentage of clay in first soil layer
! SAND1     percentage of sand in first soil layer
! LAI_VL    leaf Area Index of low vegetation
! LAI_VH    leaf Area Index of high vegetation
! WRMAX_VL  max volumetric water content retained on low vegetation
! WRMAX_VH  max volumetric water content retained on high vegetation
! ZTSL      reference height for temperature and humidity input
! ZUSL      reference height for wind input
! LAT       latitude
! FCOR      Coriolis factor
! Z0HA      AVERAGED Local roughness associated with exposed (no snow)
!           vegetation only (also no orography), for heat transfer
! WTG       Weights for SVS2 surface types as seen from GROUND
! VGH_DENS   Density of trees in areas of high vegeation (-)
! VGH_HEIGHT Height of trees in areas of high vegeation (m)
! ZVCAN   Wind speed within or above the canopy depending on CANO_REF_FORCING (m/s)
! FCANS        Canopy layer snowcover fractions
! SNCMA   Snow intercepted in high vegetation (kg m-2)
! RHOA   Air density (kg m-3)
! Z0SNOW  Roughness length for momentum of snow (m)
! Z0HSN Roughness length for heat of snow (m)
! PSNVH    fraction of HIGH vegetation covered by snow
!
!           - Output -
! HRSURF   relative humidity of the bare ground surface (1st layer)
! HUSURF    specific humidity of the bare ground surface
! HRSURFGV   relative humidity of the the snow-free ground below high veg (1st layer)
! HUSURFGV   specific humidity of the snow-free ground below high veg
! HV_VL        Halstead coefficient of low vegetation canopy
! HV_VH        Halstead coefficient of the high vegetation canopy
! HVSN_VH        Halstead coefficient of the high vegetation canopy accounting for intercepted snow
! DEL_VL    fraction of low veg. canopy covered by intercepted water
! DEL_VH    fraction of high veg canopy covered by intercepted water
! Z0HBG     Bare ground thermal roughness
! Z0HGV     Ground below high veg thermal roughness
! Z0HVL     LOW Vegetation thermal roughness
! Z0HVH     HIGH Vegetation thermal roughness
! Z0MVH      local mom roughness length for high veg.
! Z0MVL      local mom roughness length for high veg.
! RES_SNCA  Resistance for sublimation of the intercepted snow in high canopy (
  
!
      INTEGER I, zopt

      real, dimension(n) :: temp, coef_vh, qsatgr, qsat_vl, qsat_vh, qsat_sv, &
           zqs_vl, zqs_vh, ctugr, ctugrv, ctuvh, ctuvl, wcrit_hrsurf, z0bg_n,ra,&
           z0gv_n, qsatgrv,wcrit_hrsurfgv, z0hg, zz0hgv, ZZ0HVH, ZZ0HVL, TSV_CORR
     real, dimension(n) :: ZUGV, ZTGV, ZZ0MGV, ZDH, QSATI_VH, VSUBL, Z0MVH_EFF, &
                           ZDISPLCAN, ZBELOW, ZUREF_VH,ZTREF_VH,ZDIAG_TOP,ZZ0_TOP, &
                           STABM_VH, ZUSTAR_VH,LZZ0M_VH, & ! Related to VGH
                           STABT_VH,LZZ0T_VH, & ! Related to VGH
                           ZRSURF, VMOD_BELOW, & ! Related to surface below VGH
                           ZU_TOP,ZV_TOP,VMOD_TOP, &
                           ZKH, ZWCAN, ZRUPPER_CAN, ZRLOWER_CAN, &
                           VMOD_DISPL,ZRABV_CAN
     REAL :: NU, MU, NR, DVAP
     REAL :: XI2,EXT2
!
!***********************************************************************
!
!
!
!
!
!------------------------------------------------------------------------
!
!
!
!         BARE GROUND LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.
   
      DO I=1,N
         Z0BG_N(I) = Z0MBG
         Z0HG(I) = Z0M_TO_Z0H * Z0BG_N(I)
      END DO         
!
!         GROUND BELOW HIGH VEG LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.

      DO I=1,N
         Z0GV_N(I) = Z0MBG ! Value is the same as BG
         ZZ0HGV(I) = Z0M_TO_Z0H * Z0GV_N(I)
      END DO
!
!         HIGH AND LOW VEG LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.

      DO I=1,N
!         Calculate an effective roughness length for the sparse canopy that accounts for VEG_DENS
!         as well as does a weighted average for the open canopy of the Z0 of snow and BGV (from Essery)    
         Z0MVH_EFF(I) =  VGH_DENS(I) * Z0MVH(I) + (1.-VGH_DENS(I))*(Z0SNOW(I)**PSNVH(I) * Z0GV_N(I)**(1.-PSNVH(I))) 
         ZZ0HVH(I) = Z0M_TO_Z0H * Z0MVH_EFF(I)
         ZZ0HVL(I) = Z0M_TO_Z0H * Z0MVL(I)
      END DO

!          Make sure that wind speed for the canopy is not equal to 0

      DO I=1,N
         ZVCAN(I) = MAX(ZVCAN(I), 0.1)
      END DO
!
!
!
!*       1.A     RELATIVE HYMIDITY OF THE EXPOSED BARE GROUND AND OF THE
!                GROUND BELOW HIGH VEG
!               -------------------------------------------------
!
!                        This relative humidity is related to
!                        the superficial soil moisture and the
!                        field capacity of the ground
!                        ** If the 1st soil layer is very shallow (under 5cm)
!                        might need to change the calc. to use a deeper layer
!
!                        Same value since only one soil profile is used
!

      ! formulation based on field capacity
      DO I=1,N
         TEMP(I)   = PI*WD1(I)/WFC(I,1)
         HRSURF(I) = 0.5 * ( 1.-COS(TEMP(I)) )  
         wcrit_hrsurf(i) = wfc(i,1)
         HRSURFGV(I) = HRSURF(I)
         WCRIT_HRSURFGV(I) = WCRIT_HRSURF(I)
      END DO
      
!
!
!*       1.B    SPECIFIC HUMIDITY OF THE EXPOSED BARE GROUND (HU)
!               -------------------------------------------------
!                         there is a specific treatment for dew
!                         (see Mahfouf and Noilhan, jam, 1991)
!
!                         first calculate the saturation vapor
!                         pressure and specific humidity
!
      DO I=1,N
        QSATGR(I) = FOQST( TGRS(I), PS(I) )
      END DO
!
!
      DO I=1,N
!
!
!                         when hu*qsat < qa, there are two
!                         possibilities
!
!                         low-level air is dry, i.e.,
!                         qa < qsat
!
!

        IF ( HRSURF(I)*QSATGR(I).LT.HU(I).AND.QSATGR(I).GT.HU(I) )& 
                HRSURF(I) = HU(I) / QSATGR(I)

!
!
!                          b) low-level air is humid, i.e.,
!                          qa >= qsat
!
        IF ( HRSURF(I)*QSATGR(I).LT.HU(I).AND.QSATGR(I).LE.HU(I) )& 
                  HRSURF(I) = 1.0

!
!                          for very humid soil (i.e., wg > wfc ),
!                          we take hu=1
!
        IF ( WD1(I).GT.WCRIT_HRSURF(I) ) HRSURF(I) = 1.0

!
      END DO
!
!                           Calculate specific humidity over ground 
      DO I=1,N
        HUSURF(I) = HRSURF(I) * QSATGR(I)
      END DO
!
!
!
!*       1.C   SPECIFIC HUMIDITY OF THE GROUND BELOW HIGH VEGETATION (HU)
!               -------------------------------------------------
!
!                         there is a specific treatment for dew
!                         (see Mahfouf and Noilhan, jam, 1991)
!
          DO I=1,N
!
              IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN ! High vegetation is present

!                         first calculate the saturation vapor
!                         pressure and specific humidity
                          QSATGRV(I) = FOQST( TGRVS(I), PS(I) )

!                         when hu*qsat < qa, there are two
!                         possibilities
!
!                         low-level air is dry, i.e.,
!                         qa < qsat
!
                  IF ( HRSURFGV(I)*QSATGRV(I).LT.HU(I).AND.QSATGRV(I).GT.HU(I) )&
                       HRSURFGV(I) = HU(I) / QSATGRV(I)
!
!                          b) low-level air is humid, i.e.,
!                          qa >= qsat
!
                  IF ( HRSURFGV(I)*QSATGRV(I).LT.HU(I).AND.QSATGRV(I).LE.HU(I) )&
                      HRSURFGV(I) = 1.0

!
!                          for very humid soil (i.e., wg > wfc ),
!                          we take hu=1
!
                  IF ( WD1(I).GT.WCRIT_HRSURFGV(I) ) HRSURFGV(I) = 1.0

!
!
!                           Calculate specific humidity over ground
                  HUSURFGV(I) = HRSURFGV(I) * QSATGRV(I)
!
            ELSE      ! No high vegetation
!
!               ! Use value from bare ground to avoid empty arrays
!
             HUSURFGV(I) = HUSURF(I)
!
            ENDIF
!
          END DO

!
!
!
!
!
!**     2.A     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR LOW VEGETATION
!*             ------------------------------------------------------------
!
!                         first calculate the saturation vapor
!                         pressure over low vegetation
!
!
      DO I=1,N
        QSAT_VL(I) = FOQST( TVGLS(I), PS(I) )
      END DO
!
!
!*                         then calculate the fraction of the foliage
!                          covered by intercepted water (DEL)
!
!
      DO I=1,N
!
!                          Calculate the maximum value of
!                          equivalent water content in the
!                          vegetation canopy
!
!                          calculate DEL
!
         IF( WTG(I,indx_svs2_vl) .GE. EPSILON_SVS) THEN   ! Low vegetation present in the grid cell
!
             DEL_VL(I) =   (WR_VL(I)/WRMAX_VL(I))**(2./3.)
!
             DEL_VL(I) = MIN(DEL_VL(I),0.25) ! from Boone et al. (2017), same as DEL_VH

         ELSE
             DEL_VL(I) = 0.
         ENDIF
!
      END DO
!
!
!                         calculate Hv based on previous time
!                         step resavg to use in flxsurf4
!
      DO I=1,N
!
!                         calculate Hv based on previous time
!                         step resavg to calculate specific
!                         humidity of low vegetation
!
         IF(WTG(I,indx_svs2_vl) .GE. EPSILON_SVS) THEN   ! Low vegetation present in the grid cell
             HV_VL(I) = 1. - MAX(0.,SIGN(1.,QSAT_VL(I)-HU(I)))&
                   *RS(I)*(1.-DEL_VL(I)) / (RESA_VL(I)+RS(I))

!      Atmospheric resistence for exchange between the the foliage
!      and the air within the canopy space (Dearrorff, 1978)
!        RA(I)  = 100.0 / (0.3*VMOD(I) + 0.3)  !
!       Equivalent of Hv defined with respect to the mean flow inside the canopy space
!        RPP(I) = 1. - MAX(0.,SIGN(1.,QSATVG(I)-QAF(I)))&
!                 *RS(I)*(1.-DEL(I)) / (RA(I)+RS(I))

!
!                         calculate specific humidity of vegetation
!
             ZQS_VL(I) = HV_VL(I) * QSAT_VL(I) + ( 1. - HV_VL(I) ) * HU(I)
!
         ELSE
             ZQS_VL(I) = 0.
             HV_VL(I) = 0.
         ENDIF
      END DO
!
!
!
!
      if ( svs_dynamic_z0h ) then
         zopt=9

         ! TO_DO NL: which values of Z0LOC, Z0HA to use here

         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGLS, ZQS_VL, Z0LOC, Z0HA, LAT, FCOR, z0mloc=z0loc, &
              optz0=zopt, L_min=sl_Lmin_soil, &
              coeft=CTUVL, z0t_optz0=Z0HVL )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else
        ! NL: Updated with roughn. lengths for VH instead of averaged roughn. length for veg
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGLS, ZQS_VL, Z0MVL, ZZ0HVL, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUVL )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif

         do i=1,n ! TO_DO NL: delete or update that with roughn. length for heat for VL?
            z0hvl(i)=zz0hvl(i)
         enddo
      endif

      DO I=1,N
         IF(WTG(I,indx_svs2_vl) .GE. EPSILON_SVS) THEN   ! Low vegetation present in the grid cell
             RESA_VL(I) = 1. / CTUVL(I)
         ELSE
             RESA_VL(I) = 1.
         ENDIF
      END DO
!
!
!
!
!**     2.B     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR HIGH VEGETATION
!*             ------------------------------------------------------------
!
!                         first calculate the saturation vapor
!                         pressure over high vegetation
!
!
      DO I=1,N
        QSAT_VH(I) = FOQST( TVGHS(I), PS(I) )
      END DO
      QSATI_VH(:) = QSATI( TVGHS(:), PS(:) )
!
!
!*                         then calculate the fraction of the foliage
!                          covered by intercepted water (DEL_VH)    
!
!
      DO I=1,N
!
!                          Calculate the maximum value of
!                          equivalent water content in the
!                          vegetation canopy
!
!
!                          calculate DEL_VH
!
         IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
             COEF_VH(I) = 1. + 2.*LAI_VH(I)
!
             IF (HU(I) .GT. QSAT_VH(I)) THEN ! Condensation
                DEL_VH(I) = 1.
             ELSE
             	!DEL_VH(I) =   MIN(WR_VH(I),WRMAX_VH(I)) / &
                !   ( (1.-COEF_VH(I))*MIN(WR_VH(I),WRMAX_VH(I)) + COEF_VH(I)*WRMAX_VH(I) )

!               From Deardoff (1978), also used in Gouttevin et al. 2015 and same equation as for FCANS
!               Avoids using COEFF_VH as we do not know where it comes from
                DEL_VH(I) = (WR_VH(I) / WRMAX_VH(I))**(2./3.)
                DEL_VH(I) = MIN(DEL_VH(I),0.25) ! from Boone et al. (2017)
             ENDIF
          
         ELSE
             COEF_VH(I) = 1.
             DEL_VH(I) = 0.
         ENDIF
!
      END DO
!
!
!
!                         calculate HV_VH based on previous time
!                         step resa_vh to use in flxsurf4
!
      DO I=1,N
!
!                         calculate HV_VH based on previous time
!                         step RESA_VH to calculate specific 
!                         humidity of vegetation
!
         IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
               HV_VH(I) = 1. - MAX(0.,SIGN(1.,QSAT_VH(I)-HU(I)))&
                  *RS(I)*(1.-DEL_VH(I)) / (RESA_VH(I)+RS(I))

      !      Atmospheric resistance for exchange between the the foliage
!      and the air within the canopy space (Dearrorff, 1978)
      !        RA(I)  = 100.0 / (0.3*VMOD(I) + 0.3)  !
!       Equivalent of Hv defined with respect to the mean flow inside the canopy space
      !        RPP(I) = 1. - MAX(0.,SIGN(1.,QSATVG(I)-QAF(I)))&
      !                 *RS(I)*(1.-DEL(I)) / (RA(I)+RS(I))

!                        
      !                         calculate specific humidity of high vegetation
!

             ZQS_VH(I) = (1.-FCANS(I)) * (HV_VH(I) * QSAT_VH(I) + ( 1. - HV_VH(I) ) * HU(I)) &
                         + FCANS(I) * QSATI_VH(I)

!
         ELSE
               ZQS_VH(I) = 0.
               HV_VH(I) = 0.
         ENDIF
      END DO   
!
!
!          Calculate aerodynamic resistance of high vegetation
!
!

      do i=1,n 

        IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
           
           ! Compute displacement height
           ZDISPLCAN(I) = VGH_HEIGHT(I)*RCHD

           ! Compute height of the level below the canopy 
           ZBELOW(I) = HSUBCANO

           ! Compute reference forcing level for momentum for high veg.
           ! including the displacement height (as in urban_drag, TEB)
           ZUREF_VH(I) = ZUSL(I) +  VGH_HEIGHT(I) - ZDISPLCAN(I)

           ! Compute reference forcing level for heat for high veg.
           ! including the displacement height  (as in urban_drag, TEB)
           ZTREF_VH(I) = ZTSL(I) +  VGH_HEIGHT(I) - ZDISPLCAN(I)

           ! Compute height use to compute the wind speed at the top of
           ! the canopy (taking into account the displacement height)
           ZDIAG_TOP(I) = VGH_HEIGHT(I)-ZDISPLCAN(I)

           ! Compute height use derive the stability term used
           ! in the resistance above the canopy 
           ZZ0_TOP(I) = VGH_HEIGHT(I)-ZDISPLCAN(I)

        ELSE
           ZDISPLCAN(I) = 0.
           ZUREF_VH(I)  = ZUSL(I) 
           ZTREF_VH(I)  = ZTSL(I)
           ZDIAG_TOP(I) = 10.
           ZZ0_TOP(I)   =  ZZ0HVH(I)   
        ENDIF 
             
      end do
      
      if ( svs_dynamic_z0h ) then
         zopt=9
         ! TO_DO NL: Which values of Z0LOC and Z0HA to use here
         ! Which height to use? Should be above the canopy
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGHS, ZQS_VH, Z0LOC, Z0HA, LAT, FCOR, z0mloc=z0loc, &
              optz0=zopt, L_min=sl_Lmin_soil, &
              coeft=CTUVH, z0t_optz0=Z0HVH )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else


          ! First call to sfc_layer to derive the friction velocity and
          ! the wind speed at the top of the Canopy (VGH_HEIGHT). The
          ! term VGH_HEIGHT-ZDISPLCAN is used to account for the
          ! displacement height  
          i = sl_sfclayer( THETAA, HU, VMOD, VDIR,ZUREF_VH, ZTREF_VH, &
              TVGHS, ZQS_VH, Z0MVH_EFF, ZZ0HVH, LAT, FCOR, &
              hghtm_diag_row =ZDIAG_TOP,L_min=sl_Lmin_soil, &
              coeft=CTUVH, stabm=STABM_VH,lzz0m=LZZ0M_VH, &
              ue = ZUSTAR_VH, u_diag = ZU_TOP, v_diag = ZV_TOP )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif

         ! Second call to sfc_layer to determine the stability term used
         ! in the resistance above the canopy (see left term in Eq 60 of
         ! Esery et al. (2024)
          i = sl_sfclayer( THETAA, HU, VMOD, VDIR,ZUREF_VH, ZTREF_VH, &
              TVGHS, ZQS_VH, Z0MVH_EFF, ZZ0_TOP , LAT, FCOR, &
              L_min=sl_Lmin_soil,stabt=STABT_VH,lzz0t=LZZ0T_VH)

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif

         
         do i=1,n 
            z0hvh(i)=zz0hvh(i)

            IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN

             ! Compute wind speed at the canopy top from wind components
             ! at canopy top provided by sfclayer
             VMOD_TOP(I) = (ZU_TOP(I)**2. + ZV_TOP(I)**2.)**0.5

             ! Aerodynamic resistance above the canopy
             ! See left term in Eq 60 of Essery et al., 2024
             ! Stability term is derived from the call to sfc_layer
             ZRABV_CAN(I) = 1./(KARMAN*ZUSTAR_VH(I))*STABT_VH(I)

             ! Compute eddy diffusion coefficient for the canopy at
             ! height VGH_HEIGHT (eq 24 in Mahat et al, WRR, 2013)
             ZKH(I)  = KARMAN * ZUSTAR_VH(I) * (VGH_HEIGHT(I)-ZDISPLCAN(I))    

             ! Compute wind speed attenuation coefficient in the canopy
             ZWCAN(I) = WIND_CANO_COEF(LAIVH(I),VGH_DENS(I))

             ! Compute additional resistance term from the upper part of
             ! the canopy where the wind follws an exponential profile
             ! (2nd term in Eq 25 in Mahat et al, WRR, 2013)
             ZRUPPER_CAN(I) =  ( VGH_HEIGHT(I) )/(ZWCAN(I)*ZKH(I)) *&
                    ( EXP(ZWCAN(I) -ZWCAN(I)*(Z0MVH_EFF(I)+ZDISPLCAN(I))/VGH_HEIGHT(I)) - 1.)

             ! Wind speed at the base of high vegetation (height = HSUBCANO)
             ! Assuming an exponential profile in the canopy
              ! A minimum value of 0.2 m/s is used to ensure certain turbulent exchanges below the canopy
              VMOD_BELOW(I) =MAX(0.2,VMOD_TOP(I)*EXP(ZWCAN(I)*(HSUBCANO/VGH_HEIGHT(I)-1.)))

             ! Wind speed in the middle of the canopy (height = ZDISPLCAN)
             ! Assuming an exponential profile in the canopy
              VMOD_DISPL(I) =MAX(0.2,VMOD_TOP(I)*EXP(ZWCAN(I)*(ZDISPLCAN(I)/VGH_HEIGHT(I)-1.)))

              ! Compute resistance from the lower part of the canopy
              ! between  HSUBCANO and ZDISPLCAN + Z0V
              ZRLOWER_CAN(I) = ( VGH_HEIGHT(I) *EXP(ZWCAN(I)))/(ZWCAN(I)*ZKH(I)) *&
                        ( EXP(-ZWCAN(I)*HSUBCANO/VGH_HEIGHT(I)) -  &
                          EXP(-ZWCAN(I)*(Z0MVH_EFF(I)+ZDISPLCAN(I))/VGH_HEIGHT(I))) 
         
           ENDIF

         enddo
      endif

      DO I=1,N
         IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
           RESA_VH(I) = ZRABV_CAN(I) + ZRUPPER_CAN(I) 

!          TODO: should VEG_DENS impact the aero resistance for a sparse canopy?
!          IF (VGH_DENS(I) .GT. EPSILON_SVS) THEN
!             ! From Mazzotti et al. (2020). Decreases the resistance (higher flux) as VGH_DENS decreases (more open)
!             RESA_VH(I) = RESA_VH(I) / VGH_DENS(I)**0.5
!          ENDIF
         ELSE
             RESA_VH(I) = 1.
         ENDIF
      END DO
!
!
!
!
!**     2.C     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR BARE GROUND
!*             ---------------------------------------------------------------
!
!                      *************************************
!                      DANS LE CALL DE FLXSURF, LE Z0H DEVRAIT ETRE
!                      REMPLACER PAR UN Z0H_LOCAL JUSTE POUR LE SOL NU
!                      *************************************
!
      if ( svs_dynamic_z0h ) then
         zopt=9
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TGRS, HUSURF, Z0LOC, Z0HG, LAT, FCOR, optz0=zopt ,&
              z0mloc=z0loc, L_min=sl_Lmin_soil, &
              coeft=CTUGR, z0t_optz0=Z0HBG )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif

      else
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TGRS, HUSURF, Z0, Z0HG, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUGR )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif
         
         do i=1,N
            z0hbg(i)=z0hg(i)
         enddo

      endif


      DO I=1,N
        RESAGR(I) = 1. / CTUGR(I)
      END DO
!
!
!
!
!
!
!
!
!**     2.D     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR GROUND
!                 BELOW HIGH VEG.
!
!*             ---------------------------------------------------------------
!
!                      *************************************
!                      DANS LE CALL DE FLXSURF, LE Z0H DEVRAIT ETRE
!                      REMPLACER PAR UN Z0H_LOCAL JUSTE POUR LE SOL NU
!                      *************************************
!

      if ( svs_dynamic_z0h ) then
         zopt=9
         ! TO_DO NL: check the roughness lengths used
         i = sl_sfclayer( THETAA, HU, ZVCAN, VDIR, ZUSL, ZTSL, &
              TGRVS, HUSURFGV, Z0LOC, ZZ0HGV, LAT, FCOR, optz0=zopt ,&
              z0mloc=z0loc, L_min=sl_Lmin_soil, &
              coeft=CTUGRV, z0t_optz0=Z0HGV )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif

      else

         IF (CANO_REF_FORCING == 'ABV') THEN ! Reference height above the canopy. In this case, z0 should be the canopy roughness lengths and the heights above canopy

             DO I=1,N
                IF (WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN

                   ! Compute the below-canopy aerodynamic resistance with two components: 
                   !   -  exponential wind profile in the canopy (from HSUBCANO to ZDISPLCAN) ZRLOWER_CAN (no stability included )
                   !   - logarithmic wind profile below the canopy (below HSUBCANO) (no stability included) 
                   ZRSURF(I) = ZRLOWER_CAN(I)     &
                             + 1./(KARMAN**2. * VMOD_BELOW(I))   &
                             * LOG(HSUBCANO/Z0GV_N(I))*LOG(HSUBCANO/ZZ0HGV(I))

                ELSE !  This is placed to fill the values if there is no vegetation to run sl_sfclayer
                    ZRSURF(I) = 1.
                    VMOD_DISPL(I) = VMOD(I)
                    ZDISPLCAN(I) = ZUSL(I)
                ENDIF


             ENDDO

             DO I=1,N
                 ! Compute final resistance between the ground below the canopy and the atmospheric forcing level
                 CTUGRV(I) = 1. / (RESA_VH(I) + ZRSURF(I))
                 Z0HGV(I) = ZZ0HGV(I)
             ENDDO

         ELSE ! O2F or FOREST

             i = sl_sfclayer( THETAA, HU, ZVCAN, VDIR, ZUSL, ZTSL, &
                  TGRVS, HUSURFGV, Z0GV_N, ZZ0HGV, LAT, FCOR, &
                  L_min=sl_Lmin_soil, &
                  coeft=CTUGRV )

             if (i /= SL_OK) then
                call physeterror('drag_svs', 'error returned by sl_sfclayer()')
                return
             endif


             do i=1,N
                Z0HGV(i) = ZZ0HGV(i)
             enddo
         ENDIF

      endif


      DO I=1,N
        IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN
            RESAGRV(I) = 1. / CTUGRV(I)
         ELSE
            RESAGRV(I) = 1.
         ENDIF
      END DO
!
!
!
!
!
!**     2.E     SURFACE AERO RESISTANCE FOR TURBULENT FLUXES FOR SNOW
!                 BELOW HIGH VEG.
!
!*             ---------------------------------------------------------------
!
!                      *************************************
!                      DANS LE CALL DE FLXSURF, LE Z0H DEVRAIT ETRE
!                      REMPLACER PAR UN Z0H_LOCAL JUSTE POUR LE SOL NU
!                      *************************************
!


     DO I=1,N
         IF (TSV(I,1) .LT. EPSILON_SVS) THEN ! If there is no snow, TSV == 0 and it creates issues
            TSV_CORR(I) = 273.15
         ELSE
            TSV_CORR(I) = TSV(I,1)
         ENDIF   
     ENDDO

     QSAT_SV(:) = QSATI( TSV_CORR(:), PS(:) )

     IF (CANO_REF_FORCING == 'ABV') THEN ! Reference height above the canopy. 

         DO I=1,N

             IF (WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN

                 ! Compute the below-canopy aerodynamic resistance with two components: 
                 !   -  exponential wind profile in the canopy (from HSUBCANO to ZDISPLCAN) (no stability included)
                 !   - logarithmic wind profile below the canopy (below HSUBCANO) (no stability included) 
                ZRSURF(I) =  ZRLOWER_CAN(I) +                   &
                             1./(KARMAN**2. * VMOD_BELOW(I))  * &
                           LOG(HSUBCANO/Z0SNOW(I))*LOG(HSUBCANO/Z0HSN(I))

            ELSE !  This is placed to fill the values if there is no vegetation to run sl_sfclayer
                ZRSURF(I) = 1.
                VMOD_DISPL(I) = VMOD(I)
                ZDISPLCAN(I) = ZUSL(I)
            ENDIF

         ENDDO

         do i=1,N
            IF(WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN
                 ! Total resistance between the surface and the atmospheric forcing level 
                 RESA_SV(I) = RESA_VH(I) + ZRSURF(I)
            ELSE
                 ! Set to 0 in that case. If it's 0, aero resistance for snow is calculated in surface_aero_cond
                 RESA_SV(I) = 0. 
            ENDIF
         ENDDO


     ELSE ! O2F or FOREST 

         DO i=1,N
            ! Set to 0 in that case. If it's 0, aero resistance for snow is calculated in surface_aero_cond
            RESA_SV(I) = 0.
         ENDDO

     ENDIF



!*       3.     Resistance of the snow intercepted in the high vegetation (RES_SNCA)
!               ---------------------------------------------------------------
!

      DO I=1,N

         IF (FCANS(I) .GT. 0. .AND. lres_snca) THEN

            ! Fraction of the entire forest height [-]
            XI2 = 1.-ZVENT

            ! Canopy wind speed extinction coefficient [-]
            ! Ellis et al (2010) (EL10) refers to Eagleson (2002) to justify the formulation of this coefficient
            EXT2 = ZBETA * LAI_VH(I)

            ! Computation of ventilation wind speed of intercepted snow derived from above-caopny wind speed  [Eq 8 in EL10]
            ! Estimated within canopy wind speed at fraction XI2 of the entire tree height [Eq 8 in EL10]  [m s-1]
            VSUBL(I) = ZVCAN(I) * EXP(-1. * EXT2 * XI2)

            ! Sutherland's equation for kinematic viscosity
            MU = 1.8325e-5*416.16/( THETAA(I)+120.)*(THETAA(I)/296.16)*SQRT(THETAA(I)/296.16)/RHOA(I)

            ! Compute Reynolds Number     [-]
            NR = 2.0 * RADIUS_ICESPH * ZVCAN(I) / MU

            ! Compute the Nusselt Number  [-]
            NU = 1.79 + 0.606 * SQRT(NR)

            ! Compute diffusivity of water vapour in air [m2 s-1]
            DVAP = 2.063e-5 * (TVGHS(I)/273.15)**(1.75)

            ! Resistance for snow within the canopy
            RES_SNCA(I) = 2.*917.*RADIUS_ICESPH**2 /(3.*0.02*FCANS(I)**(-0.4/0.67)*SNCMA(I)*NU*DVAP)

         ELSE
            RES_SNCA(I) = 0.
         ENDIF


      END DO
!


!*       4.     HALSTEAD COEFFICIENT (RELATIVE HUMIDITY OF THE VEGETATION) (HV)
!               ---------------------------------------------------------------
!
!
!                          Here we calculate "true" new HV based on updated
!                          resavg
!
      DO I=1,N


        HV_VL(I) = 1. - MAX(0.,SIGN(1.,QSAT_VL(I)-HU(I)))&
                 *RS(I)*(1.-DEL_VL(I)) / (RESA_VL(I)+RS(I))

        HV_VH(I) = 1. - MAX(0.,SIGN(1.,QSAT_VH(I)-HU(I)))&
                 *RS(I)*(1.-DEL_VH(I)) / (RESA_VH(I)+RS(I))

        ! Account for intercepted snow in the high vegetation
        HVSN_VH(I) = FCANS(I) * (CHLF + CHLC) *  RESA_VH(I) / (RESA_VH(I) + RES_SNCA(I)) &
                    + (1.-FCANS(I)) * HV_VH(I) * CHLC


      END DO
!


      RETURN
      END
