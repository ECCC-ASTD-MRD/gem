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


subroutine drag7(TS, WG, WR, THETAA, VMOD, VDIR, HU, &
     PS, RS, VEG, Z0H, Z0TOT, WFC, &
     PSNG, PSNV, LAI, ZUSL, ZTSL, LAT, FCOR, &
     RESA, ILMO, HST, FRV, FTEMP, FVAP, &
     CH, CD, HRSURF, HUSURF, HV, DEL, ZQS, &
     CTU, N)
   use tdpack
   use sfclayer, only: sl_sfclayer,SL_OK
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>
   !@Object
   !     Calculates the drag coefficients for heat and momentum transfers
   !     over ground (i.e., Ch and Cd).
   !@Arguments
   !          - Input/Output -
   ! RESA      aerodynamical surface resistance
   ! ILMO
   ! HST
   ! FRV
   ! FTEMP
   ! FVAP
   !          - Input -
   ! TS        surface temperature
   ! WG        superficial volumetric water content
   ! WR        water content retained by the vegetation canopy
   ! THETAA    potential temperature at the lowest level
   ! HU        specific humidity of air at the lowest level
   ! VMOD      wind speed at the lowest level
   ! VDIR      wind direction at the lowest level
   ! PS        surface pressure
   ! RS        surface or stomatal resistance
   ! VEG       fraction of a model grid area covered by vegetation
   ! Z0H       roughness length for heat transfers
   ! Z0TOT     roughness length including the effect of snow
   ! WFC       volumetric water content at the field capacity
   ! PSNG      fraction of bare ground covered by snow
   ! PSNV      fraction of vegetation covered by snow
   ! LAI       leaf area index
   ! ZTSL      reference height for temperature and humidity input
   ! ZUSL      reference height for wind input
   ! LAT       latitude
   ! FCOR      Coriolis factor
   !           - Output -
   ! CH        drag coefficient for heat
   ! CD        drag coefficient for momentum
   ! HRSURF    relative humidity of the surface
   ! HUSURF    specific humidity of the surface
   ! HV        Halstead coefficient (i.e., relative humidity of the
   !           vegetation canopy)
   ! DEL       fraction of canopy covered by intercepted water
   ! ZQS       area-averaged relative humidity of a model tile
   ! CTU       homogeneous boundary condition term in the
   !           diffusion equation for Theta and Q

   integer N
   real TS(N), WG(N), WR(N), THETAA(N), VMOD(N), VDIR(N), HU(N)
   real PS(N), RS(N), VEG(N), Z0TOT(N), WFC(N)
   real Z0H(N)
   real PSNG(N), PSNV(N), LAI(N), ZUSL(N), ZTSL(N)
   real LAT(N), FCOR(N)
   real RESA(N), ILMO(N), HST(N), FRV(N), FTEMP(N), FVAP(N)
   real CH(N), CD(N), HUSURF(N), HV(N), DEL(N), ZQS(N)
   real HRSURF(N), CTU(N)

   !@Author S. Belair (January 1997)
   !@Revisions
   ! 001      S. Belair (November 1998)
   !             Use FLXSURF1 to calculate the surface transfer
   !             coefficients instead of the MOMCOEF and HEATCOEF
   !             subroutines (Mascart method).
   !
   ! 002      S. Belair (November 1998)
   !             Remove Z0 from the arguments (because now we
   !             use FLXSURF1.
   !
   ! 003      B. Bilodeau (January 2001)
   !             Automatic arrays
   !
   ! 004      S. Belair (September 2001)
   !             Add protection to calculation of DEL
   !
   ! 005      Y. Delage (September 2004)
   !             Replace ZA by ZUSL and ZTSL, FLXSURF3 by FLXSURF3,
   !             UE2 by FRV and rename subroutine DRAG2
   !
   ! 006      M. Abrahamowicz (May 2013)
   !             Replace call to flxsurf3 by call to flxsurf4, with optz0=0 (neutral change)
   !
   !@Notes
   ! Method
   !     1) computes hu, hv, and DEL
   !
   !     2) use this to find qsoil, the grid-averaged relative humidity
   !        of the soil
   !
   !     3) find the transfer and resistance coefficients Ch, Cd, and Ra
   !        Calculate the surface fluxes of heat, moisture,
   !        and momentum over water surfaces.

   integer :: i, zopt
   real :: ue
   real, dimension(n) :: TEMP, WRMAX, QSAT, COEF, CMU

   !------------------------------------------------------------------------

   ! 1. RELATIVE AND SPECIFIC HUMIDITY OF THE GROUND (HU)
   !    -------------------------------------------------

   ! this relative humidity is related to
   ! the superficial soil moisture and the
   ! field capacity of the ground

   do I=1,N
      TEMP(I)   = PI*WG(I)/WFC(I)
      HRSURF(I) = 0.5 * ( 1.-cos(TEMP(I)) )
   end do

   ! there is a specific treatment for dew
   ! (see Mahfouf and Noilhan, jam, 1991)

   ! first calculate the saturation vapor
   ! pressure and specific humidity

   do I=1,N
      QSAT(I) = FOQST( TS(I), PS(I) )
   end do

   do I=1,N

      !  when hu*qsat < qa, there are two possibilities

      !  low-level air is dry, i.e., qa < qsat

      if ( HRSURF(I)*QSAT(I).lt.HU(I).and.QSAT(I).gt.HU(I) ) &
           HRSURF(I) = HU(I) / QSAT(I)

      ! b) low-level air is humid, i.e., qa >= qsat

      if ( HRSURF(I)*QSAT(I).lt.HU(I).and.QSAT(I).le.HU(I) ) then
         HRSURF(I) = 1.0
      endif

      !  for very humid soil (i.e., wg > wfc ), we take hu=1

      if ( WG(I).gt.WFC(I) ) &
           HRSURF(I) = 1.0

   end do

   do I=1,N
      HUSURF(I) = HRSURF(I) * QSAT(I)
   end do

   ! 2.     FRACTION OF THE FOLIAGE COVERED BY INTERCEPTED WATER (DEL)
   !   ------------------------------------------------------------

   ! first calculate the maximum value of equivalent water content in the
   ! vegetation canopy

   do I=1,N

      WRMAX(I) = 0.2 * VEG(I) * LAI(I)

      !   calculate DEL

      COEF(I) = 1. + 2.*LAI(I)

      if ( VEG(I).gt.0.0.and.WRMAX(I).gt.0.0 ) then
         DEL(I) =   min(WR(I),WRMAX(I))                             / &
              ( (1.-COEF(I))*min(WR(I),WRMAX(I)) + COEF(I)*WRMAX(I) )
      else
         DEL(I) = 0.0
      end if

   end do

   ! 3. HALSTEAD COEFFICIENT (RELATIVE HUMIDITY OF THE VEGETATION) (HV)
   !    ---------------------------------------------------------------

   do I=1,N
      HV(I) = 1. - max(0.,sign(1.,QSAT(I)-HU(I))) &
           *RS(I)*(1.-DEL(I)) / (RESA(I)+RS(I))
   end do

   ! 4. GRID-AVERAGED HUMIDITY OF THE SOIL (ZQS)
   !    ----------------------------------------

   do I=1,N
      ZQS(I) = ( (1.-VEG(I))*(1.-PSNG(I))*HRSURF(I) &
           +     (1.-VEG(I))*    PSNG(I) &
           +     VEG(I) *(1.-PSNV(I))*HV(I) &
           +     VEG(I) *    PSNV(I)  )*QSAT(I) &
           +     VEG(I) *(1.-PSNV(I))*(1.-HV(I))*HU(I)
   end do

   ! 5. SURFACE TRANSFER COEFFICIENTS FOR HEAT AND MOMENTUM (CH and CD)
   !    ---------------------------------------------------------------
   if (z0tevol == 'ZILI95') then
      zopt = 9
   elseif (z0tevol == 'FIXED') then
      zopt = 0
   else
      call physeterror('drag', 'unknown option for z0tevol='//trim(z0tevol))
      return
   endif

   i = sl_sfclayer(THETAA,HU,VMOD,VDIR,ZUSL,ZTSL,TS,ZQS, &
        Z0TOT,Z0H,LAT,FCOR,optz0=zopt,L_min=sl_Lmin_soil,spdlim=VMOD, &
        coefm=CMU,coeft=CTU,flux_t=FTEMP,flux_q=FVAP,ilmo=ILMO,ue=FRV,h=HST)
   if (i /= SL_OK) then
      call physeterror('drag', 'error returned by sl_sfclayer()')
      return
   endif

   do I=1,N
      UE       =  FRV(I)
      CMU(I)   = CMU(I) / UE

      CD(I) = CMU(I) * CMU(I)
      CH(I) = CMU(I) * CTU(I)/UE

      RESA(I) = 1. / CH(I) / VMOD(I)
   end do

   ! 7. HALSTEAD COEFFICIENT (WITH THE NEW VALUES OF RESA)
   !    --------------------------------------------------

   do I=1,N
      HV(I) = 1. - max(0.,sign(1.,QSAT(I)-HU(I))) &
           *RS(I)*(1.-DEL(I)) / (RESA(I)+RS(I))
   end do

   return
end subroutine drag7
