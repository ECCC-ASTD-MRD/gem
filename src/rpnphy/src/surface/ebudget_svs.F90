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

      SUBROUTINE EBUDGET_SVS(TSA, WD1, WF , &
                   TGRS,TGRD,TVGS,TVGD, & 
                   DT, VMOD, VDIR, LAT, & 
                   RG, ALVG, LAI, GAMVEG, & 
                   ALGR,EMGR, & 
                   RAT, THETAA, FCOR, ZUSL, ZTSL, HU, PS, &  
                   RHOA, WTA, Z0, Z0LOC, Z0H, & 
                   HRSURF, HV, DEL, RS, & 
                   CG,CVP, EMIS, PSNG, &  
                   RESAGR, RESAVG, RESASA, RESASV, &
                   RNETSN, HFLUXSN, LESNOFRAC, ESNOFRAC, & 
                   ALPHAS, &  
                   TSNS, & 
                   RNETSV, HFLUXSV, LESVNOFRAC, ESVNOFRAC,  &
                   ALPHASV, & 
                   TSVS, & 
                   VEGH, VEGL, PSNVH,PSNVHA, & 
                   RR,WR,SNM,SVM, &
                   ALBT, & 
                   RNET, HFLUX, LE, LEG, LEV, LES,LESV, & 
                   LER, LETR, EG, ER, ETR, GFLUX, EFLUX, & 
                   BM, FQ, BT, RESAEF, & 
                   LEFF, DWATERDT, & 
                   FTEMP, FVAP, ZQS, FRV, & 
                   ALFAT, ALFAQ, ILMO, HST, TRAD, N)
      use tdpack
      use sfclayer, only: sl_sfclayer,SL_OK
      use sfc_options
      use svs_configs
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N
      
      REAL TSA(N),  DT, VMOD(N)
      REAL VDIR(N), LAT(N)
      REAL TGRS(N), TGRD(N), TVGS(N), TVGD(N)
      REAL TSNS(N)
      REAL WD1(N), WF(N)
      REAL RG(N), ALVG(N), RAT(N), THETAA(N), FCOR(N)
      REAL ZUSL(N), ZTSL(N)
      REAL HU(N), PS(N), RHOA(N), WTA(N,svs_tilesp1), Z0(N)
      REAL Z0LOC(N), Z0H(N)
      REAL HV(N), DEL(N), RS(N)
      REAL CG(N), CVP(N),  PSNG(N), EMIS(N)
      REAL LAI(N), GAMVEG(N), ALGR(N), EMGR(N)
      REAL RNET(N), HFLUX(N), LE(N), ALPHAS(N)
      REAL ALBT(N)
      REAL LEG(N), LEV(N), LER(N), LETR(N), GFLUX(N)
      REAL EFLUX(N), BM(N), FQ(N), BT(N), LES(N)
      REAL FTEMP(N), FVAP(N), ER(N), ETR(N)
      REAL LEFF(N), DWATERDT(N), ZQS(N), FRV(N)
      REAL EG(N), HRSURF(N)
      REAL RESAGR(N), RESAVG(N), RESASA(N), RESASV(N), RESAEF(N)
      REAL RNETSN(N), HFLUXSN(N), LESNOFRAC(N), ESNOFRAC(N)
      REAL RNETSV(N), HFLUXSV(N), LESVNOFRAC(N), ESVNOFRAC(N)
      REAL TSVS(N), ALPHASV(N)
      REAL ALFAT(N), ALFAQ(N), LESV(N)
      REAL VEGH(N), VEGL(N), PSNVH(N), PSNVHA(N)
      REAL ILMO(N), HST(N), TRAD(N)
      REAL WR(N),RR(N),SNM(N),SVM(N)
!
!Author
!          S. Belair et al. (January 2016)
!Revisions
! 001      Bug fixes: M. Abrahamowicz, S. Z. Husain, N. Gauthier, E. Gaborit, 
!          V. Vionnet, D. Deacu      
!
!Object
!
!     Calculates the evolution of the surface and deep-soil temperature
!     (i.e., Ts and T2), as well as all the surface fluxes.
!
!
!!!!  METHOD
!!    ------
!
!     1- find the grid-averaged albedo, emissivity, and roughness length
!     2- compute the za, zb, and zc terms involved in the numerical
!        resolution of the equations for Ts and T2.
!     3- find Ts(t) and T2(t).
!     4- derive the surface fluxes.!
!
!Arguments
!
!
!          - Input/Output -
! TGRS      (bare) ground temperature -- S for "skin"
! TGRD      mean ground temperature -- D for "deep"
! TVGS      vegetation temperature -- S for "skin"
! TVGD      mean vegetation temperature -- D for "deep"
! TS        surface  temperature (new) as seen from ground
! Z0H       agg. thermal roughness length for land surface
!           Output only when svs_dynamic_z0h=.true.
!          - Input -
! VMOD      module of the low-level wind
! VDIR      wind direction at the lowest level
! LAT       latitude
! WD1       soil water content of the first deep layer 
! WF        frozen soil water
! DT        timestep
! RG        global radiation (downward solar)
! ALVG      AVERAGED surface albedo associated with vegetation type
! LAI       AVERAGED vegetation leaf area index
! GAMVEG    AVERAGED parameter related to the vegetation height
! ALGR      albedo of bare ground (soil)
! EMGR      emissivity of bare ground (soil)   
! RAT       atmospheric radiation incident on the ground (NIR)
! THETAA    air potential temperature at the lowest level
! FCOR      Coriolis factor
! ZTSL      reference height for temperature and humidity input
! ZUSL      reference height for wind input
! HU        specific humidity of air at the lowest level
! PS        surface pressure
! RHOA      air density near the surface
! WTA       Weights for SVS surface types as seen from SPACE
! Z0        momentum roughness length (no snow)
! Z0LOC     local land momentum roughness length (no orography) 
! HRSURF    relative humidity of the bare ground surface (1st soil layer)
! HV        Halstead coefficient (relative humidity of veg. canopy)
! DEL       portion of the leaves covered by water
! RS        stomatal resistance
! CG        soil thermal coefficient
! CVP       AVERAGED vegetation thermal coefficient (with LAI effect)
! EMIS      AVERAGED surface emissivity when vegetation fully grown
! PSNG      fraction of bare ground or low veg. covered by snow
! RESAGR    aerodynamical surface resistance for bare ground
! RESAVG    aerodynamical surface resistance for vegetation
! RESASA    aerodynamical surface resistance for snow on bare ground/low veg
! RESASV    aerodynamical surface resistance for snow under high veg.
! RNETSN    net radiation over snow 
! HFLUXSN   sensible heat flux over snow 
! ESNOFRAC  water vapor flux from the snow surface
! LESNOFRAC latent heat flux from the snow surface 
! ALPHAS    albedo of snow
! TSNS      surface (skin) snow temperature at time t+dt (update in snow_alone.ftn)
!
! RNETSV    net radiation over snow-under-vegetation
! HFLUXSV   sensible heat flux over snow-under-vegetation
! ESVNOFRAC water vapor flux from the snow-under-vegetation
! LESVNOFRAC latent heat flux from the snow-under-vegetation 
! ALPHASV   albedo of snow-under-veg
! TSVS      surface (skin) snow-under-veg temperature at time t+dt (update in snow_alone.ftn)
!
! VEGH     fraction of HIGH vegetation
! VEGL     fraction of LOW vegetation
! PSNVH    fraction of HIGH vegetation covered by snow
! PSNVHA   fraction of HIGH vegetation covered by snow as seen from
!          the ATMOSPHERE 
! RR       Liquid precipitation rate at the surface in [mm/s]
! WR       Water retained by vegetation
! SVM      snow water equivalent (SWE) for snow under high veg [kg/m2]
!
!           - Output -
! ALBT      total surface albedo (snow + vegetation + bare ground)
! RNET      net radiation
! HFLUX     sensible heat flux
! LE        latent heat flux
! LEG       latent heat flux over bare ground
! LEV       latent heat flux over vegetation
! LES       latent heat flux over snow
! LESV      latent heat flux over snow-under-veg
! LER       direct latent heat flux from vegetation leaves
! LETR      evapotranspiration latent heat flux
! EG         evaporation rate (no fraction) over bare ground (1st soil layer)
! ER        direct evaporation rate (no fraction) from veg. leaves
! ETR       evapotranspiration rate (no faction) 
!
! GFLUX     ground flux
! EFLUX     water vapor flux
! DWATERDT  net tendency of melting-freezing of soil water
! TS        surface  temperature (new) as seen from ground
! TD        mean soil temperature
! TSA       surface  temperature (new) as seen from space
! LEFF      effective latent heat
! FTEMP      land sfc avg turbulent surface flux of temperature
! FVAP       land sfc avg turbulent surface flux of vapor
! ILMO       land sfc avg (1/length of Monin-Obukov) 
! HST        land sfc avg height of the boundary layer 
! FRV        land sfc average friction velocity 
! BM         homogeneous boundary condition term in the
!            diffusion equation for U and V
! BT         homogeneous boundary condition term in the
!            diffusion equation for T and Q
! FQ         land sfc average momentum flux
! RESAEF     effective aerodynamic resistance for land sfc
! ALFAT      inhomogeneous boundary term in the diffusion equation for Theta
! ALFAQ      inhomogeneous boundary term in the diffusion equation for Q
! ZQS       area-averaged specific  humidity of a model tile
! TRAD      averaged radiative temperature

!
!
      INTEGER I,zopt
!
!
      REAL EMISSN, EMSOIL, KCOEF, RHOW
      REAL BFREEZ, RAIN1, RAIN2
      REAL ABARK
!
!     MULTIBUDGET VARIABLES 
!     GR:ground, SN:snow, VG:vegetation, AG: aggregated 
       real, dimension(n) :: a2, b2, c2, a3, b3, c3, zhv, freezfrac, emvg, &
            alvglai, zqsatgr, zdqsatgr, zqsatvg, zdqsatvg, zqsatgrt, zqsatvgt, &
            rnetgr, rnetvg, hfluxgr, hfluxvg, roragr, roravg,  &
            zqsatsno, tgrst, tgrdt, tvgst, tvgdt, esf, esvf, evf, &
            egf, ev, zqsatsnv, levnofrac, legnofrac, frach,  &
            cmu, cm, ctu, vmod_lmin

!************************************************************************
!
!
!
!                                THE FOLLOWING SHOULD BE PUT IN 
!                                A COMMON COMDECK
!
      EMISSN = 0.97
      EMSOIL = 0.94
      RHOW   = 1000.  
      KCOEF  = 1.E-6
      BFREEZ = 4.
!                                Albedo of Bark (S. Wang, Ecological Modelling, 2005)
      ABARK  = 0.15
!
!
      RAIN1  = 2.8e-8
      RAIN2  = 2.8e-7
! 
!
!
!!       1.     GRID-AVERAGED ALBEDO, EMISSIVITY, AND ROUGHNESS LENGTH
!       ------------------------------------------------------
!                          (considering snow surfaces)
!
      DO I=1,N
!
!
!                               Calculate grid-averaged albedo 
!
!                               Using PSNVHA
!        
        ALBT(I) = AG( WTA(I,indx_svs_bg), WTA(I,indx_svs_vg), &
                      WTA(I,indx_svs_sn), WTA(I,indx_svs_sv), &
                      ALGR(I),ALVG(I),ALPHAS(I),ALPHASV(I))

! 
!                               Recalculate vegetation-only albedo to take LAI
!                               effect into account, can only consider it 
!                               if have high Vegetation (i.e. Gamveg greater
!                               or equal to 0.01). In the high vegetation case,
!                               weight albedo of leaves vs. bark according to LAI
!
            IF(GAMVEG(I).GE.0.01)THEN 
               ALVGLAI(I) = MIN( LAI(I)      , LAI0 )   * ALVG(I) / LAI0 &  
                          + MAX( LAI0-LAI(I) , 0.0    ) * ABARK  / LAI0
            ELSE
               ALVGLAI(I) = ALVG(I)
            ENDIF
!
!
!                            
!                               Calculate vegetation-only emissivity,
!                               the only bark emissivity i could find
!                               was for eucalyptus trees which is not
!                               very helpful as the bark is usually very
!                               bright
        EMVG(I)  = EMIS(I) 
!
      END DO
!        2.     LATENT HEAT COEFFICIENTS - CORRECTION DUE TO FREEZING
!               AND MELTING OF SOIL WATER
!               -----------------------------------------------------
!
!                               Using the fraction of frozen water
!                               in the soil, calculate the "effective"
!                               latent heat of evaporation/sublimation
!
        DO I=1,N
           FREEZFRAC(I) = WF(I) / (WD1(I)+WF(I)+EPSILON_SVS)
           LEFF(I)      = FREEZFRAC(I)      * (CHLC+CHLF)  &
                         + (1.-FREEZFRAC(I)) *  CHLC
        END DO
!
!
!
!
!
!!       3A.     COEFFICIENTS FOR THE TIME INTEGRATION OF  TGRS
!                     (i.e. BARE SOIL/GROUND SKIN TEMPERATURE)
!               --------------------------------------------
!
!                            Thermodynamic functions
!
!
       DO I=1,N
          ZQSATGR(I)  = FOQST( TGRS(I),PS(I) )
          ZDQSATGR(I) = FODQS( ZQSATGR(I),TGRS(I) )
       END DO
!
!
!                              function zrsra
!
       DO I=1,N
          RORAGR(I) = RHOA(I) / RESAGR(I)
       END DO
!
!                          
!
!                                        terms za, zb, and zc for the
!                                              calculation of tgrs(t)
!
       DO I=1,N

          A2(I) = 1. / DT + CG(I) * & 
                 (4. * EMGR(I) * STEFAN * (TGRS(I)**3) &   
                 +  RORAGR(I) * ZDQSATGR(I) * LEFF(I)* HRSURF(I) &  
                 +  RORAGR(I) * CPD) &  
                 + 2. * PI / 86400.

          B2(I) = 1. / DT + CG(I) *  &  
                 (3. * EMGR(I) * STEFAN * (TGRS(I)**3) &   
                 + RORAGR(I) * ZDQSATGR(I) * LEFF(I) * HRSURF(I) )


          C2(I) = 2. * PI * TGRD(I) / 86400. &   
                 + CG(I) *  & 
                 (RORAGR(I) * CPD * THETAA(I) &   
                 + RG(I) * (1. - ALGR(I)) + EMGR(I)*RAT(I)  &  
                 - RORAGR(I) *  & 
                 LEFF(I)* (HRSURF(I)*ZQSATGR(I)-HU(I)))
       END DO
!
!
!
       DO I=1,N
          TGRST(I) = ( TGRS(I)*B2(I) + C2(I) ) / A2(I) 
       ENDDO
!
!!       3B.     COEFFICIENTS FOR THE TIME INTEGRATION OF  TVGS
!                     (i.e. VEGETATION SKIN TEMPERATURE)
!               --------------------------------------------
!         CALCULATE ONLY IF VEGETATION NON-ZERO PRESENT, OTHERWISE USE BARE GROUND TO IMPOSE DEFAULT PHYSICAL VALUE
!
       DO I=1,N      
          IF ( (VEGH(I)+VEGL(I)*(1-PSNG(I))).GE.EPSILON_SVS ) THEN
             ! VEGETATION PRESENT
!
!
!                            Thermodynamic functions
!
             
             ZQSATVG(I)  = FOQST( TVGS(I),PS(I) )
             ZDQSATVG(I) = FODQS( ZQSATVG(I),TVGS(I) )
!
!                              function zrsra      
!
             RORAVG(I) = RHOA(I) / RESAVG(I)
!
!                              Fraction of high vegetation to total vegetation
!    
             FRACH(I) =  MIN(   (VEGH(I)*PSNVH(I))/(VEGH(I)+VEGL(I)*(1-PSNG(I))) , 1.0) 
!
!                                        terms za, zb, and zc for the
!                                              calculation of tvgs(t)
             A3(I) = 1. / DT + CVP(I) *  & 
                    (4. * EMVG(I) * STEFAN * (TVGS(I)**3)  &  
                    +  RORAVG(I) * ZDQSATVG(I) * CHLC * HV(I) &  
                    +  RORAVG(I) * CPD )  & 
                    + 2. * PI / 86400.
!
             B3(I) = 1. / DT + CVP(I) *   &
                    (3. * EMVG(I) * STEFAN * (TVGS(I)** 3)  &   
                    + RORAVG(I) * ZDQSATVG(I) * CHLC* HV(I) )
           
!
             C3(I) = 2. * PI * TVGD(I) / 86400. &   
                     + CVP(I) *  & 
                     ( RORAVG(I) * CPD * THETAA(I)  &  
                     + RG(I) * (1. - ALVGLAI(I)) + EMVG(I)*RAT(I)  & 
                     - RORAVG(I)  & 
                     * CHLC * HV(I) * (ZQSATVG(I)-HU(I)) )



             TVGST(I) =  ( TVGS(I)*B3(I) + C3(I) ) / A3(I)

          ELSE
             ! NO VEGETATION -- USE BARE GROUND VALUES or ZERO to fill arrays to avoid numerical errors
             ZQSATVG(I)  =  ZQSATGR(I)
             ZDQSATVG(I) = ZDQSATGR(I)
             RORAVG(I) = RORAGR(I)
             FRACH(I) = 0.0
             TVGST(I) = TGRST(I)
          ENDIF


       ENDDO

!
!
!!$!        4.     FREEZING AND MELTING TENDENCIES FOR SOIL WATER
!!$!               ----------------------------------------------


!     SET DWATERDT to zero, just to initialize array
      DWATERDT=0.0

!
!
!!       5.A     TGRD AT TIME 'T+DT'
!               -----------------
!
      DO I=1,N
        TGRDT(I) = (TGRD(I) + DT*TGRST(I)/86400.) /   &  
                      (1.+DT/86400.)
!
!
!
!                            Include the effect of soil freeze/thaw
! note dwaterdt already multiplied by DT above, so have DT^2 ...below
!chekc if this is ok
! ** DO NOT INCLUDE THE EFFECT OF DWATERDT BECAUSE SET TO ZERO ABOVE !
!

        TGRDT(I) = TGRDT(I) 

      END DO
!
!!       5.B     TVGD AT TIME 'T+DT'
!               -----------------
!
      DO I=1,N
!                   Note that as an added precaution,
!                   we set the vegetation temperature to
!                   that of the ground, when no vegetation is present
!
         IF(VEGH(I)+VEGL(I)*(1.-PSNG(I)).ge.EPSILON_SVS)THEN
            TVGDT(I) = (TVGD(I) + DT*TVGST(I)/86400.) / (1.+DT/86400.)
         ELSE
            TVGDT(I) = TGRDT(I)
         ENDIF
            
      END DO
!
!
!
!!       7.     FLUX CALCULATIONS
!               -----------------
!
!

      DO I=1,N
!                                            recalculate the qsat functions
!
        ZQSATGRT(I)  = FOQST(  TGRST(I)  ,  PS(I)   )
        ZQSATVGT(I) =  FOQST(  TVGST(I)  ,  PS(I)   )
!
      ENDDO
!                     
!
!
      DO I=1,N
!
!                                            NET RADIATION 
!                                            ---------------
!
!                                            Net radiation over bare ground
!
        RNETGR(I) = (1. - ALGR(I)) * RG(I) + EMGR(I) *&  
                 (RAT(I) - STEFAN * (TGRST(I)** 4))
!
!                                            Net radiation over vegetation
!
        RNETVG(I) = (1. - ALVGLAI(I)) * RG(I) + EMVG(I) *&  
                 (RAT(I) - STEFAN * (TVGST(I)** 4))
!
!                                            AGGREGATED net radiation (including snow)
!                                    
        RNET(I) = AG( WTA(I,indx_svs_bg), WTA(I,indx_svs_vg), &
                      WTA(I,indx_svs_sn), WTA(I,indx_svs_sv), &
                      RNETGR(I),RNETVG(I),RNETSN(I),RNETSV(I) )
!        
!
!                                            SENSIBLE HEAT FLUX 
!                                            ---------------
!
!
!                                            Sensible heat flux from the ground
!
        HFLUXGR(I) = RHOA(I) * CPD * (TGRST(I) - THETAA(I)) / RESAGR(I)
!
!                                            Sensible heat flux from the vegetation
!
        HFLUXVG(I) = RHOA(I) * CPD * (TVGST(I) - THETAA(I)) / RESAVG(I)
!
!                                             AGGREGATED sensible heat flux (including snow)
        HFLUX(I)   =  AG( WTA(I,indx_svs_bg), WTA(I,indx_svs_vg), &
                          WTA(I,indx_svs_sn), WTA(I,indx_svs_sv), &
                          HFLUXGR(I),HFLUXVG(I),HFLUXSN(I),HFLUXSV(I))

!
!
!
!                                             AGGREGATED turbulent surface flux of temperature
!
!        FTEMP(I) = HFLUX(I) / ( RHOA(I) * CPD )
!
!
!
!                                            LATENT HEAT FLUXES 
!                                            ---------------
!
!
!                                            Latent heat of evaporation from
!                                            the ground
!
        LEGNOFRAC(I) = RHOA(I) * LEFF(I) * (HRSURF(I)* ZQSATGRT(I) - HU(I)) / RESAGR(I)
        LEG(I) = WTA(I,indx_svs_bg) * LEGNOFRAC(I)          
!
!
!                                            Water vapor flux from ground
        EGF(I) = WTA(I,indx_svs_bg) * (HRSURF(I)* ZQSATGRT(I) - HU(I)) / RESAGR(I)
!
!                                            Evaporation rate from ground (for hydro_svs.ftn)
!
        EG(I) = RHOA(I)*(HRSURF(I)* ZQSATGRT(I) - HU(I)) / RESAGR(I)

!
! 
!                                            FOR VEGETATION --- SET FLUXES TO ZERO if NO VEGETATION
          IF(VEGH(I)+VEGL(I)*(1.-PSNG(I)).ge.EPSILON_SVS) THEN
!

!                                            Latent heat of evaporation from
!                                            vegetation
!
             LEVNOFRAC(I) = RHOA(I) * CHLC * HV(I) * (ZQSATVGT(I) - HU(I)) / RESAVG(I)



!

!                                            Latent heat of Transpiration
!
!
             ZHV(I) = MAX(0.0 , SIGN(1.,ZQSATVGT(I) - HU(I)))
             LETR(I) = ZHV(I) * (1. - DEL(I)) * RHOA(I)  &
                       * CHLC * (VEGH(I) * (1. - PSNVHA(I)) + VEGL(I) * (1. - PSNG(I))) & 
                       * (ZQSATVGT(I)  - HU(I)) / (RESAVG(I) + RS(I))
!
!                                           Transpiration rate (for hydro_svs.ftn)
!
             ETR(I) = RHOA(I)*ZHV(I)*(1. - DEL(I))*(ZQSATVGT(I) - HU(I))/(RESAVG(I) + RS(I))


!                                            Evapotranspiration rate from vege. (for hydro_svs.ftn)

             EV(I) =  RHOA(I)*HV(I) * (ZQSATVGT(I) - HU(I)) / RESAVG(I)
!
!                                            Latent heat of evapotranspiration from
!                                            vegetation
!
             LEVNOFRAC(I) = RHOA(I) * CHLC * HV(I) * (ZQSATVGT(I) - HU(I)) / RESAVG(I)

             !  EV is limited to WR/DT+RR+ETR to avoid negative WR in hydro_svs when direct evaporation exceeds rainrate
             !  When snow is present, rain falls through vegetation to snow bank... so is not considered in evaporation... This is to conserve water budget.


             IF( SNM(I).GE.CRITSNOWMASS .AND. SVM(I).GE. CRITSNOWMASS ) THEN
                !both snow packs exists, rain falls directly to snow
                EV(I) = MIN (EV(I),(WR(I)/DT+ETR(I)))

             ELSE IF ( SNM(I).GE.CRITSNOWMASS .AND. SVM(I).LT.CRITSNOWMASS ) THEN
                ! only low vegetation snow pack present, rain for low vegetation portion
                ! falls through to snow, rain for high vegetation portion is considered in evaporation
                EV(I) = MIN (EV(I),(WR(I)/DT+ETR(I)+ & 
                     RR(I)*VEGH(I)/(VEGH(I)+VEGL(I)*(1.-PSNG(I)))))
   
             ELSE IF ( SNM(I).LT.CRITSNOWMASS .AND. SVM(I).GE.CRITSNOWMASS ) THEN
                ! only high vegetation snow pack present, rain for high vegetation portion
                ! falls through to snow, rain for low vegetation portion is considered in evaporation
                EV(I) = MIN (EV(I),(WR(I)/DT+ETR(I)+ &
                     RR(I)*(VEGL(I)*(1-PSNG(I))/(VEGH(I)+VEGL(I)*(1.-PSNG(I))))))

             ELSE 
                ! no snow present, all rain is considered evaporation
                EV(I) = MIN (EV(I),(WR(I)/DT+ETR(I)+RR(I)))
             ENDIF
        

        ELSE
           ! NO VEGETATION --- SET FLUXES TO ZERO
           LEVNOFRAC(I) = 0.0
           ZHV(I) = 0.0
           LETR(I) = 0.0
           ETR(I) = 0.0
           EV(I) = 0.0
           LEVNOFRAC(I) = 0.0
          
        ENDIF

!                                            Water vapor flux from vegetation
        EVF(I) =  WTA(I,indx_svs_vg)  * EV(I)/RHOA(I)
!                                            Latent heat of evaporation from vegetation
        LEV(I) = RHOA(I) * CHLC * EVF(I)
!
!
!                                            Direct evapo. rate from veg.(for hydro_svs.ftn)
!
        ER(I) = EV(I) - ETR(I)
!
!                                            Latent heat of direct evaporation
!
        LER(I)  = LEV(I) - LETR(I)
!
!                                            Calculate latent heat snow weighted
!                                            by grid-cell snow-coverage fraction
!       
        LES(I)  =  WTA(I,indx_svs_sn) *  LESNOFRAC(I)
        ESF(I)  =  WTA(I,indx_svs_sn) *  ESNOFRAC(I)
!
!                                            Same for snow-under-vegetation
!
        LESV(I) =  WTA(I,indx_svs_sv)   *  LESVNOFRAC(I)
        ESVF(I) =  WTA(I,indx_svs_sv)   *  ESVNOFRAC(I)
!
!                                            Total latent heat of evaporation
!                                            (Including snow contribution)
!
        LE(I) = LEG(I) + LEV(I) + LES(I) + LESV(I)
!
!                                            Total water vapor flux
!                                            (Including snow contribution)
        EFLUX(I) = EGF(I) + EVF(I) + ESF(I) + ESVF(I)
!        FVAP(I)  = EFLUX(I)
!
!                                            Heat flux into the ground
!
        GFLUX(I) = RNET(I) - HFLUX(I) - LE(I)
!
      ENDDO
!
!
!
!
!*       8.     NEW "land-tile"-AVERAGED QUANTITIES  (FRV, ZQS, BM, FQ)
!               -----------------------------------------------  
!
!
!
      DO I=1,N
!                                             Re-calculate snow saturation humidity
!                                             instead of passing it from snow_alone.ftn
!                                             (N.B. snow always saturated 
!                                             i.e., specific=saturation humidity)
!
        ZQSATSNO(I)=FOQST( TSNS(I), PS(I) )
        ZQSATSNV(I)=FOQST( TSVS(I), PS(I) )
        !
     ENDDO

     IF ( .NOT. use_eff_surf_tq  ) THEN
        !  Area-average weighted mean calculation for land sfc temperature and humidity

        !
        DO i=1,n
           ! Calculate land-tile-averaged specific humidity
!
           ZQS(I) =     WTA(I,indx_svs_bg)      *    HRSURF(I)        * ZQSATGRT(I) & 
                      + WTA(I,indx_svs_sn)                            * ZQSATSNO(I) &  
                      + WTA(I,indx_svs_sv)                            * ZQSATSNV(I) &  
                      + WTA(I,indx_svs_vg)      *        HV(I)        * ZQSATVGT(I) &
                      + WTA(I,indx_svs_vg)      *    (1.-HV(I))       * HU(I)
!
!              Calculate land-tile-averaged surface temperature
!              i.e., aggregate skin temperatures of diff. surfaces 
!
           TSA(I) = AG( WTA(I,indx_svs_bg), WTA(I,indx_svs_vg), &
                        WTA(I,indx_svs_sn), WTA(I,indx_svs_sv), &
                        TGRST(I),TVGST(I),TSNS(I),TSVS(I) )

        ENDDO
              
     else
        ! Consider aerodynamic resistance in calc. of effective ("land-tile-averaged") land sfc 
        ! temperature and humidity (see svs_configs.ftn90)

!
!        Calculate effective areodynamic resistance
!         
        DO i=1,n
           RESAEF(I) = 1. / ( WTA(I,indx_svs_bg)/RESAGR(I) + WTA(I,indx_svs_vg)/RESAVG(I) + &
                              WTA(I,indx_svs_sn)/RESASA(I) + WTA(I,indx_svs_sv)/RESASV(I) )
!
!          Calculate effective land sfc specific humdity
!         
           ZQS(I) = RESAEF(I) *                                                             &
                 (  WTA(I,indx_svs_bg)      *    HRSURF(I)        * ZQSATGRT(I)/RESAGR(I) & 
                 +  WTA(I,indx_svs_sn)                            * ZQSATSNO(I)/RESASA(I) &  
                 +  WTA(I,indx_svs_sv)                            * ZQSATSNV(I)/RESASV(I) &  
                 + (WTA(I,indx_svs_vg)      *        HV(I)        * ZQSATVGT(I)           & 
                   +  WTA(I,indx_svs_vg)      *    (1.-HV(I))       * HU(I) )/RESAVG(I) )
!
!          Calculate effective land sfc temperature
!         
           TSA(I) = RESAEF(I) *  AG( WTA(I,indx_svs_bg), WTA(I,indx_svs_vg), &
                                     WTA(I,indx_svs_sn), WTA(I,indx_svs_sv), &
                                     TGRST(I)/RESAGR(I), TVGST(I)/RESAVG(I), &
                                     TSNS(I)/RESASA(I),  TSVS(I)/RESASV(I) )

        ENDDO
     endif


!    Calculated averaged radiative Temperature

     DO i=1,n

        TRAD(I)= AG( WTA(I,indx_svs_bg), WTA(I,indx_svs_vg), &
                     WTA(I,indx_svs_sn), WTA(I,indx_svs_sv), &
                     TGRST(I)**4,TVGST(I)**4,TSNS(I)**4,TSVS(I)**4 )
        
        TRAD(I)=TRAD(I)**(1./4.)
     ENDDO
!
!
!                                             Calculate surface layer transfer coefficients
!                                             and fluxes using ZQS and TSA. 
!                                             Here want FRV (average friction velocity for 
!                                             grid-cell, used in diasurf2.ftn), CMU and 
!                                             CTU(Turbulent transfer coeff. for thermodynamics) 
!
!

! DO we still need FRV (UE), ILMO, HST ???? !!!! as output
!
!
     if( svs_dynamic_z0h ) then
        zopt=9
        i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
                         TSA, ZQS, Z0, Z0H, LAT, FCOR, &
                         optz0=zopt, z0mloc=Z0LOC, &
                         L_min=sl_Lmin_soil,spdlim=vmod_lmin, &
                         coefm=CMU, coeft=CTU, &
                         flux_t=FTEMP, flux_q=FVAP, &
                         ilmo=ILMO, ue=FRV, h=HST, &
                         z0t_optz0=Z0H )


        if (i /= SL_OK) then
           call physeterror('ebudget_svs', 'error returned by sl_sfclayer()')
           return
        endif
     else
      
        i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
                         TSA, ZQS, Z0, Z0H, LAT, FCOR, &
                         L_min=sl_Lmin_soil,spdlim=vmod_lmin, &
                         coefm=CMU, coeft=CTU, &
                         flux_t=FTEMP, flux_q=FVAP, &
                         ilmo=ILMO, ue=FRV, h=HST )


        if (i /= SL_OK) then
           call physeterror('ebudget_svs', 'error returned by sl_sfclayer()')
           return
        endif
        
   endif


      DO I=1,N
!
!
!        TERMS SAME FOR IMPLICIT AND EXPLICIT FORMULATION
        
         CM(I)  = CMU(I) / FRV(I)

         if ( sl_Lmin_soil > 0.) then
            ! use wind module consistent with imposed Monin-Obukhov length
            bm(i) = vmod_lmin(i) * (cm(i)**2)  

            fq(i) = rhoa(i) * bm(i) * vmod_lmin(i)

         else
            bm(i) = vmod(i) * (cm(i)**2)  

            fq(i) = rhoa(i) * bm(i) * vmod(i)
         endif



      ENDDO

      
      ! -------  EXPLICIT FORMULATION
      IF (.NOT.IMPFLX) THEN
         DO I=1,N   
 
         !   inhomogeneous boundary term in the diffusion equation for Theta
         !   alfat = ctu * (ta -t_surf) = -ftemp
            ALFAT(I)   =  -FTEMP(I)
!            inhomogeneous boundary term in the diffusion equation for Q
         !   alfaq = ctu * (qa -q_surf) = -fvap
            ALFAQ(I)   =  -FVAP(I)

         ! homogeneous boundary condition term in the
         ! diffusion equation for T and Q          
            BT(I) = 0.0

         ENDDO
      ENDIF


      ! -------  IMPLICIT FORMULATION
      IF (IMPFLX) THEN
         DO I=1,N   
      
         !  inhomogeneous boundary term in the diffusion equation for Theta
            ALFAT(I)   =  -CTU(I)  * TSA(I)

         ! inhomogeneous boundary term in the diffusion equation for Q
            ALFAQ(I)   =  -CTU(I)  * ZQS(I)

         ! homogeneous boundary condition term in the
         ! diffusion equation for T and Q          
            BT(I) = CTU(I)

         ENDDO
      ENDIF

!
!*       9.     UPDATE TEMPERATURE VARIABLES
!              -----------------------------
!
      DO I=1,N
        TGRS(I)   = TGRST(I)
        TGRD(I)   = TGRDT(I)
        TVGS(I)   = TVGST(I)
        TVGD(I)   = TVGDT(I)
      ENDDO
!
!
      RETURN
      END
