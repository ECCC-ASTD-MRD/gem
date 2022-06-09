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

      SUBROUTINE DRAG_SVS ( TGRS, TVGS, WD1, &
                              WR, THETAA, VMOD, VDIR, HU,  &  
                              PS, RS, Z0, Z0LOC, Z0VG, WFC, WSAT, CLAY1,  &
                              SAND1, LAI, WRMAX, ZUSL, ZTSL, LAT, & 
                              FCOR, Z0HA, & 
                              RESAGR, RESAVG, &  
                              HUSURF, & 
                              HRSURF, &       
                              HV, DEL,  &
                              Z0HBG, Z0HVG, &
                              N)
      use tdpack
      use sfc_options
      use sfclayer, only: sl_sfclayer,SL_OK
      use svs_configs
!
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N
      REAL TGRS(N), TVGS(N), WR(N), THETAA(N), VMOD(N), VDIR(N), HU(N), CLAY1(N)
      REAL SAND1(N), PS(N), RS(N), Z0(N),Z0LOC(N), Z0VG(N), WFC(N,NL_SVS), WSAT(N,NL_SVS)
      REAL LAI(N), WRMAX(N), ZUSL(N), ZTSL(N), LAT(N)
      REAL FCOR(N), Z0HA(N), Z0HBG(N), Z0HVG(N)
      REAL RESAGR(N), RESAVG(N)
      REAL HUSURF(N), HV(N), DEL(N)
      REAL HRSURF(N), WD1(N)
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
! RESAVG    aerodynamical surface resistance for vegetation
!
!          - Input -
! TGRS      skin (surface) temperature of bare ground
! TVGS      skin (surface) temperature of vegetation
! WD1       Soil volumetric water content (first level)
! WR        water content retained by the vegetation canopy
! THETAA    potential temperature at the lowest level
! VMOD      wind speed at the lowest level
! VDIR      wind direction at the lowest level
! HU        specific humidity of air at the lowest level
! PS        surface pressure
! RS        surface or stomatal resistance
! Z0        momentum roughness length (no snow)
! Z0LOC     local (no orography) land momentum roughness
! Z0VG      vegetation-only momentum roughness      
! WFC       volumetric water content at the field capacity
! WSAT      volumetric water content at saturation
! CLAY1     percentage of clay in first soil layer
! SAND1     percentage of sand in first soil layer
! LAI       AVERAGED leaf area index
! WRMAX     Max volumetric water content retained on vegetation
! ZTSL      reference height for temperature and humidity input
! ZUSL      reference height for wind input
! LAT       latitude
! FCOR      Coriolis factor
! Z0HA      AVERAGED Local roughness associated with exposed (no snow)
!           vegetation only (also no orography), for heat transfer
!
!           - Output -
! HRSURF   relative humidity of the bare ground surface (1st layer)
! HUSURF    specific humidity of the bare ground surface
! HV        Halstead coefficient (i.e., relative humidity of the
!           vegetation canopy)
! DEL       fraction of canopy covered by intercepted water
! Z0HBG     Bare ground thermal roughness
! Z0HVG     Vegetation thermal roughness     
!
!
!
!
      INTEGER I, zopt

      real, dimension(n) :: temp, coef, qsatgr, qsatvg, &
           zqsvg, ctugr, ctuvg, z0hg, wcrit_hrsurf, z0bg_n
           
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
         Z0HG(I) = Z0M_TO_Z0H * Z0MBG
      END DO         
!
!
!
!*       1.     RELATIVE AND SPECIFIC HUMIDITY OF THE GROUND (HU)
!               -------------------------------------------------
!
!                        This relative humidity is related to
!                        the superficial soil moisture and the
!                        field capacity of the ground
!                        ** If the 1st soil layer is very shallow (under 5cm)
!                        might need to change the calc. to use a deeper layer

      if( svs_hrsurf_sltext ) then
         !use hrsurf formulation based on soil texture

         DO I=1,N
             ! Set soil water max for hrsurf calc. based on clay percentage)
            if ( clay1(i) .lt. 1.0 ) then
               wcrit_hrsurf(i) = wfc(i,1)
            else if ( clay1(i)  .lt. 40.0 ) then
               wcrit_hrsurf(i) =  (sand1(i)/((sand1(i)+clay1(i)))) * wfc(i,1) & 
                    +             (clay1(i)/((sand1(i)+clay1(i)))) * wsat(i,1) 
            else
               wcrit_hrsurf(i)= wsat(i,1)
            endif

            TEMP(I)   = PI*WD1(I)/WCRIT_HRSURF(I)
            HRSURF(I) = 0.5 * ( 1.-COS(TEMP(I)) )  
            
         END DO
      else
         ! formulation based on field capacity
         DO I=1,N
    
            TEMP(I)   = PI*WD1(I)/WFC(I,1)
            HRSURF(I) = 0.5 * ( 1.-COS(TEMP(I)) )  
            wcrit_hrsurf(i) = wfc(i,1)
         END DO
      endif
!
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
!**     2.A     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR BARE GROUND 
!*             ---------------------------------------------------------------
!
!
!
!
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
!**     2.B     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR VEGETATION 
!*             ------------------------------------------------------------
!
!
!
!                         first calculate the saturation vapor
!                         pressure over vegetation
!
!
      DO I=1,N
        QSATVG(I) = FOQST( TVGS(I), PS(I) )
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
!
!                          calculate DEL
!
         COEF(I) = 1. + 2.*LAI(I)
!
 
         DEL(I) =   MIN(WR(I),WRMAX(I)) / &
               ( (1.-COEF(I))*MIN(WR(I),WRMAX(I)) + COEF(I)*WRMAX(I) )
!
         DEL(I) = MIN(DEL(I),0.1) 
          

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
!                         humidity of vegetation
!
        HV(I) = 1. - MAX(0.,SIGN(1.,QSATVG(I)-HU(I)))&  
                 *RS(I)*(1.-DEL(I)) / (RESAVG(I)+RS(I))
!                        
!                         calculate specific humidity of vegetation
!
        ZQSVG(I) = HV(I) * QSATVG(I) + ( 1. - HV(I) ) * HU(I)
!
      END DO   
!
!
!
!
      if ( svs_dynamic_z0h ) then
         zopt=9
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGS, ZQSVG, Z0LOC, Z0HA, LAT, FCOR, z0mloc=z0loc, &
              optz0=zopt, L_min=sl_Lmin_soil, &
              coeft=CTUVG, z0t_optz0=Z0HVG )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGS, ZQSVG, Z0, Z0HA, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUVG )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
         
         do i=1,n
            z0hvg(i)=z0ha(i)
         enddo
      endif
   
      DO I=1,N
         RESAVG(I) = 1. / CTUVG(I)
      END DO

!*       3.     HALSTEAD COEFFICIENT (RELATIVE HUMIDITY OF THE VEGETATION) (HV)
!               ---------------------------------------------------------------
!
!
!                          Here we calculate "true" new HV based on updated
!                          resavg
!
      DO I=1,N

        HV(I) = 1. - MAX(0.,SIGN(1.,QSATVG(I)-HU(I)))& 
                 *RS(I)*(1.-DEL(I)) / (RESAVG(I)+RS(I))

      END DO
!
      RETURN
      END
