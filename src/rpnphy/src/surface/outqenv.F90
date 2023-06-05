
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

!/@*
subroutine OUTQENV2(PSCA_SW, PREF_SW_FAC, PREF_SW_GRND,    &
     PEMIT_LW_FAC, PEMIT_LW_GRND, PLW_RAD, &
     PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,    &
     PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,N,OPT,OPT_BODY,          &
     ZEB,ZAB,ZHB)
   use tdpack_const, only: PI
   implicit none
!!!#include <arch_specific.hf>

   !    PURPOSE       : Computes the energy budget received by a body over a slab surface
   !    AUTHOR        : S. Leroyer   (Original  10/2016)
   !    REFERENCE     : Leroyer et al. (2018)
   !    MODIFICATIONS :
   !    METHOD        :

   !*       0.     DECLARATIONS

   !*      0.1    declarations of arguments
   integer, intent(IN)  :: N
   integer, intent(IN)  :: OPT
   integer, intent(IN)  :: OPT_BODY
   real, dimension(N), intent(IN)  :: PREF_SW_GRND !Solar radiation reflected by ground [road + garden] (W/m²)
   real, dimension(N), intent(IN)  :: PREF_SW_FAC !Solar radiation reflected by facade [wall + glazing] (W/m²)
   real, dimension(N), intent(IN)  :: PSCA_SW !Diffuse solar radiation (W/m²)
   real, dimension(N), intent(IN)  :: PDIR_SW !Direct solar radiation (W/m²)
   real, dimension(N), intent(IN)  :: PZENITH !solar zenithal angle (rad from vert.)
   real, dimension(N), intent(IN)  :: PEMIT_LW_FAC !Longwave radiation emitted by the facade [wall + glazing] (W/m²)
   real, dimension(N), intent(IN)  :: PEMIT_LW_GRND !Longwave radiation emitted by the ground [road + garden] (W/m²)
   real, dimension(N), intent(IN)  :: PLW_RAD !Atmospheric longwave radiation (W/m²)
   real, dimension(N), intent(IN)  :: PBLD !Building surface fraction
   real, dimension(N), intent(IN)  :: PBLD_HEIGHT !Building surface fraction
   real, dimension(N), intent(IN)  :: PWALL_O_HOR !Building surface fraction
   real, dimension(N), intent(OUT)  :: PQ1
   real, dimension(N), intent(OUT)  :: PQ2
   real, dimension(N), intent(OUT)  :: PQ3
   real, dimension(N), intent(OUT)  :: PQ4
   real, dimension(N), intent(OUT)  :: PQ5
   real, dimension(N), intent(OUT)  :: PQ6
   real, dimension(N), intent(OUT)  :: PQ7

   real, intent(IN) :: ZHB  ! average height of the body
   real, intent(IN) :: ZAB  !absorption coef of solar radiation by the  body
   real, intent(IN) :: ZEB  !emissivity of  body
   !   ##########################################################################

   real, dimension(N) :: ZWROAD !width of the road (m)
   real, dimension(N) :: ZWROOF !width of the roof (m)                                !!! INPUT ?????
   real, dimension(N) :: ZL1, ZL2, ZL4 !lengths for view factor calculation
   real, dimension(N) :: ZHD, ZL3, zz, yy, kk !lengths for view factor calculation
   real, dimension(N) :: ZFFAC !facade view factor of human body
   real, dimension(N) :: ZFGRND !ground view factor of human body
   real, dimension(N) :: ZFROOF! roof view factor of human body
   real, dimension(N) :: ZFSKY !sky view factor of human body
   real, dimension(N) :: ZDIRSWBODY !solar radiation received by human body
   real, dimension(N) :: ZELEV !solar elevation angle

   integer :: JJ

   do JJ = 1, N

      !*  1 - calculation of view factors
      ZFFAC(JJ) =  0.
      ZFGRND(JJ) = 0.
      ZFROOF(JJ) = 0.
      ZFSKY(JJ) = 0.
      PQ1(JJ)=0.0
      PQ2(JJ)=0.0
      PQ3(JJ)=0.0
      PQ4(JJ)=0.0
      PQ5(JJ)=0.0
      PQ6(JJ)=0.0
      PQ7(JJ)=0.0

      !=====================
      if (OPT .eq. 1) then    ! flat surface
         !======================
         ZFGRND(JJ) = 0.5
         ZFSKY(JJ) = 1.0 - ZFGRND(JJ)

         !======================
      elseif (OPT .eq. 2) then  ! street
         !=======================

         ZWROAD(JJ) = PBLD_HEIGHT(JJ) * 2. * (1. - PBLD(JJ)) / PWALL_O_HOR(JJ)

         ZL1(JJ) = sqrt(ZHB**2                   + (ZWROAD(JJ)/2.)**2)
         ZL2(JJ) = sqrt( PBLD_HEIGHT(JJ)**2      + (ZWROAD(JJ)/2.)**2)
         ZL4(JJ) = sqrt((PBLD_HEIGHT(JJ)-ZHB)**2 + (ZWROAD(JJ)/2.)**2)

         ! wall view factor
         ZFFAC (JJ) = (ZL1(JJ) + ZL2(JJ) - ZWROAD(JJ)/2. - ZL4(JJ)) / (2. * ZHB)
         ! ground view factor
         ZFGRND(JJ) = 0.5*ZWROAD(JJ)/ZHB
         ZFGRND(JJ) = 0.5 * (ZFGRND(JJ) + 1. - sqrt(ZFGRND(JJ)**2 + 1.) )
         ! sky view factor
         ZFSKY (JJ) = 1. - ZFFAC(JJ) - ZFGRND(JJ)

         !========================
      elseif (OPT .eq. 3) then   ! rooftop
         !========================

         ZWROAD(JJ) = PBLD_HEIGHT(JJ) * 2. * (1. - PBLD(JJ)) / PWALL_O_HOR(JJ)
         ZWROOF(JJ) = PBLD_HEIGHT(JJ) * 2. * PBLD(JJ) / PWALL_O_HOR(JJ)

         ! roof view factor
         ZFROOF(JJ) = 0.5 * ZWROOF(JJ)/ZHB
         ZFROOF(JJ) = 0.5 * (ZFROOF(JJ) + 1. - sqrt(ZFROOF(JJ)**2 + 1.) )

         ! ground (road) view factor --> simplification assume 0 although could be >0 if ZHD>BLDH
         ZFGRND(JJ) = 0.0

         ! height d of the facade seen by the body on the roof
         ZHD(JJ) = min(PBLD_HEIGHT(JJ), 2.0* ZHB* ZWROAD(JJ)/ZWROOF(JJ) )

         ! coefficients to compute facade view factor of one portion of the body zz
         yy(JJ)=0.5*ZHD(JJ)* ZWROOF(JJ)/ZWROAD(JJ)
         if (yy(JJ) .gt. ZHB) then
            yy(JJ)=ZHB
         endif
         zz(JJ)=ZHB - yy(JJ)

         kk(JJ) = ZWROAD(JJ) + (ZWROOF(JJ)/2.)

         ZL1(JJ) = sqrt( yy(JJ)**2            + kk(JJ)**2 )
         ZL2(JJ) = sqrt( (ZHD(JJ)+ZHB)**2     + kk(JJ)**2 )
         ZL3(JJ) = sqrt( ZHB**2              + kk(JJ)**2 )
         ZL4(JJ) = sqrt( (ZHD(JJ)+yy(JJ))**2   + kk(JJ)**2 )

         ! wall view factor (ZHD seen)
         if (zz(JJ) .lt.  0.1) then
            ZFFAC(JJ)  = 0.0
         else
            ZFFAC(JJ)  = ( ZL1(JJ) + ZL2(JJ) - ZL3(JJ) - ZL4(JJ) ) / (2. * ZZ(JJ))
         endif

         ! sky view factor  (can be a slighty different from 0.5)
         ZFSKY(JJ)  = 1.0 - ZFROOF(JJ) - ZFFAC(JJ) - ZFGRND(JJ)
         !========================
      endif

      ! Energy Fluxes (short-wave / longwave)

      ! sky
      PQ2(JJ)=ZAB*PSCA_SW(JJ)*ZFSKY(JJ)
      PQ3(JJ)=ZEB*PLW_RAD(JJ)*ZFSKY(JJ)
      ! ground
      PQ4(JJ)=ZAB*PREF_SW_GRND (JJ)*ZFGRND(JJ)
      PQ5(JJ)=ZEB*PEMIT_LW_GRND(JJ)*ZFGRND(JJ)
      ! facade
      PQ6(JJ)=ZAB*PREF_SW_FAC (JJ)*ZFFAC(JJ)
      PQ7(JJ)=ZEB*PEMIT_LW_FAC(JJ)*ZFFAC(JJ)
      ! roof => see urban_outqenv

      !=====================================
      ! ADD CONTRIBUTION FROM THE SUN
      ! OPT_BODY = 1  human body, for UTCI
      ! OPT_BODY = 2  globe sensor, for WBGT
      !======================================

      ZELEV(JJ) = ( PI/2. - PZENITH(JJ)) *180./PI
      if (ZELEV(JJ) .lt. 1E-6) then
         ZELEV(JJ) = 0.
      endif

      if (OPT_BODY .eq. 1) then
         ! the direct solar radiation is weighted by a projected area factor which can be expressed by this equation
         ! for a rotationally symmetric human being (Fanger, 1970)
         ZDIRSWBODY(JJ) = PDIR_SW(JJ) / max( cos(PZENITH(JJ)) ,0.1)              &
              * 0.308 * cos( PI/180. *ZELEV(JJ)* (1.-ZELEV(JJ)**2/48402.) )
      elseif (OPT_BODY .eq. 2) then
         ! sphere Gaspar and Quintela 2009
         ZDIRSWBODY(JJ) = PDIR_SW(JJ) * 0.25 /  max( cos(PZENITH(JJ)) ,0.1)
      endif

      !==
      PQ1(JJ)=ZAB* ZDIRSWBODY(JJ)

   enddo

   return
end subroutine OUTQENV2
