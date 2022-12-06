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

subroutine vegi(RG, T, TS, HU, PS, W2, RGL, LAI, RSMIN, GAMMA, WWILT, &
     WFC, RS, N)
   use tdpack
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>
   !@Object Calculates the surface stomatal resistance Rs
   !@Arguments
   !          - Input -
   ! RG       solar radiation
   ! T        low-level temperature of air
   ! TS       surface temperature
   ! HU       low-level specific humidity of air
   ! PS       surface pressure
   ! W2       soil volumetric water content
   ! RGL      constant for the calculation of the stomatal resistance
   ! LAI      Leaf area index
   ! RSMIN    minimum stomatal resistance
   ! GAMMA    other constant for RS
   ! WWILT    volumetric water content at the wilting point
   ! WFC      volumetric water content at the field capacity
   !          - Output -
   ! RS       Surface or stomatal resistance

   integer :: N
   real :: RG(N), T(N), HU(N), PS(N), W2(N), TS(N)
   real :: RGL(N), LAI(N), RSMIN(N), GAMMA(N), WWILT(N)
   real :: WFC(N), RS(N)

   !@Author S. Belair (January 1997)
   !@Revisions
   ! 001      B. Bilodeau (January 2001)

   integer :: i
   real, dimension(n) :: f, f1,f2,f3,f4, qsat

   !***********************************************************************
   
   ! 1.     THE 'ZF1' FACTOR
   !         ---------------
   !   This factor measures the influence
   !   of the photosynthetically active radiation
   do I=1,N
      F(I)  = 0.55*2.*RG(I) / (RGL(I)+1.E-6) / &
           ( LAI(I)+1.E-6 )
      F1(I) = ( F(I) + RSMIN(I)/5000. ) / ( 1. + F(I) )
   end do

   ! 2.     THE 'ZF2' FACTOR
   !        ----------------
   !  This factor takes into account the effect
   !  of the water stress on the surface resistance
   !
   !  NOTE that W2 (liquid portion of soil water) is
   !  used here instead of W2+WF.  Thus, when soil water
   !  freezes (W2 --> 0), ZF2 becomes small and the
   !  surface increases increases (no transpiration when soils are frozen).
   do I=1,N

      !            For humid soils, this factor does not
      !            increase the stomatal resistance

      if (W2(I).ge.WFC(I)) F2(I) = 1.0

      !            The stomatal resistance should be large
      !            when the soil is very dry

      if (W2(I).lt.WFC(I).and.W2(I).le.WWILT(I)) &
           F2(I) = 1.E-5


      !            For intermediate soils:

      if (W2(I).lt.WFC(I).and.W2(I).gt.WWILT(I)) &
           F2(I) = ( W2(I)-WWILT(I) ) / &
           ( WFC(I)-WWILT(I) + 1.E-6 )

   end do

   ! 3.     THE 'ZF3' FACTOR
   !         ----------------
   !  This factor represents the effect of
   !  vapor pressure deficit of the atmosphere.
   !  For very humid air, the stomatal resistance
   !  is small, whereas it increases as the
   !  air becomes drier.
   do I=1,N
      QSAT(I) = FOQST( TS(I), PS(I) )
      F3(I) = max( 1. - GAMMA(I)*( QSAT(I) &
           - HU(I) )*1000. , 1.E-3 )
   end do

   ! 4.     THE 'ZF4' FACTOR
   !        ----------------
   !  This factor introduces an air temperature
   !  dependance on the surface stomatal resistance
   do I=1,N
      F4(I) = max( 1.0 - 0.0016*(298.15-T(I))**2, 1.E-3 )
   end do


   ! 5.     THE SURFACE STOMATAL RESISTANCE
   !        -------------------------------
   do I=1,N
      RS(I) = veg_rs_mult * RSMIN(I) / ( LAI(I)+1.E-6 ) &
           / F1(I) / F2(I) / F3(I) / F4(I)

      RS(I) = min( RS(I),5000.  )
      RS(I) = max( RS(I), 1.E-4 )
   end do

   return
end subroutine vegi
