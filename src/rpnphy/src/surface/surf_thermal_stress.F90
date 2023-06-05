
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
subroutine SURF_THERMAL_STRESS(PTA, PQA,                &
     PU10, PUSR, PPS,                                   &
     PDIR_SW, PSCA_SW, PLW_RAD, PZENITH,                &
     PREF_SW_SURF, PEMIT_LW_SURF,                       &
     PUTCI_OUTSUN, PUTCI_OUTSHADE,                      &
     WBGT_SUN, WBGT_SHADE,                              &
     PTRAD_HSUN, PTRAD_HSHADE,                          &
     PTGLOBE_SUN, PTGLOBE_SHADE, PTWETB,                &
     PQ1_H,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H,PQ7_H ,N)
   use tdpack_const, only: TCDK
   implicit none

   !    PURPOSE       : COMPUTES THERMAL STRESS INDICATORS over a slab surface
   !    AUTHOR        : S. Leroyer   (Original  10/2016)
   !    REFERENCE     : Leroyer et al. (2018) , urban climate
   !    MODIFICATIONS : S Leroyer (2020), wetbulb in C
   !    METHOD        :
   !--------------------------------------------------------------------------

   !*       0.     DECLARATIONS
   !               ------------
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !*      0.1    declarations of arguments
   integer, intent(IN)  :: N

   real, dimension(N), intent(IN)  :: PTA        ! Air temperature (K)
   real, dimension(N), intent(IN)  :: PQA        !  Air specific humidity  (kg/kg)
   real, dimension(N), intent(IN)  :: PPS        !  air pressure at the surface

   real, dimension(N), intent(IN)  :: PU10       !  wind speed at 10m (m/s)
   real, dimension(N), intent(IN)  :: PUSR       !  Air wind speed 10-m over the roof at the sensor level (m/s)

   real, dimension(N), intent(IN)  :: PZENITH       ! solar zenithal angle (rad from vert.)
   real, dimension(N), intent(IN)  :: PDIR_SW       ! Direct solar radiation (W/m²)
   real, dimension(N), intent(IN)  :: PSCA_SW       ! Diffuse solar radiation (W/m²)
   real, dimension(N), intent(IN)  :: PLW_RAD       ! Longwave radiation (W/m²)

   real, dimension(N), intent(IN)  :: PREF_SW_SURF  ! Solar radiation reflected by ground [road + garden] (W/m²)
   real, dimension(N), intent(IN)  :: PEMIT_LW_SURF  ! Longwave radiation emitted by the ground [road + garden] (W/m²)

   real, dimension(N), intent(OUT) :: PUTCI_OUTSUN   ! UTCI for outdoor person at sun (°C)
   real, dimension(N), intent(OUT) :: PUTCI_OUTSHADE ! UTCI for outdoor person in shade (°C)

   real, dimension(N), intent(OUT) :: WBGT_SUN         ! WBGT  wet bulb globe temperature in the street (C)
   real, dimension(N), intent(OUT) :: WBGT_SHADE       ! WBGT  wet bulb globe temperature in the street (C)

   ! optional to get out
   real, dimension(N), intent(OUT), optional :: PTRAD_HSUN        ! for a human body
   real, dimension(N), intent(OUT), optional :: PTRAD_HSHADE
   real, dimension(N), intent(OUT), optional :: PTGLOBE_SUN      ! Globe Temperature in the exposed street (K)
   real, dimension(N), intent(OUT), optional :: PTGLOBE_SHADE    ! Globe Temperature in the shaded street (K)

   real, dimension(N), intent(OUT) :: PTWETB          ! wet-bulb temperature in the street (C)

   real, dimension(N), intent(OUT)   :: PQ1_H  ! energy components for the standing standard clothed human
   real, dimension(N), intent(OUT)   :: PQ2_H
   real, dimension(N), intent(OUT)   :: PQ3_H
   real, dimension(N), intent(OUT)   :: PQ4_H
   real, dimension(N), intent(OUT)   :: PQ5_H
   real, dimension(N), intent(OUT)   :: PQ6_H
   real, dimension(N), intent(OUT)   :: PQ7_H
   ! REAL, DIMENSION(N), INTENT(OUT), OPTIONAL   :: PQ8_H
   ! REAL, DIMENSION(N), INTENT(OUT), OPTIONAL   :: PQ9_H
   ! REAL, DIMENSION(N), INTENT(OUT), OPTIONAL   :: PQ10_H
   ! REAL, DIMENSION(N), INTENT(OUT), OPTIONAL   :: PQ11_H
   ! REAL, DIMENSION(N), INTENT(OUT), OPTIONAL   :: PQ12_H
   ! REAL, DIMENSION(N), INTENT(OUT), OPTIONAL   :: PQ13_H

   !  declarations of local variables
   real, dimension(N) :: ZEHPA !water vapour pressure (hPa)
   real, dimension(N) :: ZUNDEF
   ! energy components for the globe
   real, dimension(N) :: PQ1_G
   real, dimension(N) :: PQ2_G
   real, dimension(N) :: PQ3_G
   real, dimension(N) :: PQ4_G
   real, dimension(N) :: PQ5_G
   real, dimension(N) :: PQ6_G
   real, dimension(N) :: PQ7_G
   ! MRT for the globe
   real, dimension(N) :: PTRAD_GSUN
   real, dimension(N) :: PTRAD_GSHADE

   !#WARNING: by initializing these var, they are authomatically "saved", make it explicit save (or parameter)
   real, save :: ZEB_G = 0.957 !emissivity of a globe sensor          WBGT
   real, save :: ZEB_H = 0.97  !emissivity of clothed human body      UTCI
   real, save :: ZAB_H = 0.7   !absorption coef of solar radiation by human body        UTCI
   real, save :: ZAB_G = 0.957 !absorption coef of solar radiation by globe sensor      WBGT
   real, save :: ZHB_H = 1.7  !average height of human person (m)                       UTCI & WBGT
   real, save :: ZHB_G = 1.7  ! 2.5 !panam  !average height of the globe sensor
   real, save :: ZGD = 0.148  ! black globe sensor diameter in m (value given by Matt Wright for PanAm2015)
   integer :: ZOPT
   integer :: ZOPT_BODY
   integer :: JJ

   real, external :: TGLOBE_BODY_SURF
   real, external :: WETBULBT_SURF
   real, external :: MRT_BODY_SURF
   real, external :: UTCI_APPROX_SURF

   do JJ = 1, N
      ZUNDEF(JJ)=0.0
   enddo

   !========================================================
   ! COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY (standard clothed standing human)
   !========================================================
   ZOPT_BODY=1

   ZOPT=1   ! CANOPY
   call OUTQENV2(PSCA_SW, ZUNDEF, PREF_SW_SURF,  &
        ZUNDEF, PEMIT_LW_SURF, PLW_RAD,     &
        ZUNDEF, ZUNDEF, ZUNDEF, PDIR_SW, PZENITH,   &
        PQ1_H,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H,PQ7_H,  &
        N,ZOPT,ZOPT_BODY,ZEB_H,ZAB_H, ZHB_H )

   !========================================================
   ! COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY (globe sensor)
   !========================================================
   ZOPT_BODY=2

   ZOPT=1   ! canopy
   call OUTQENV2(PSCA_SW, ZUNDEF, PREF_SW_SURF,  &
        ZUNDEF, PEMIT_LW_SURF, PLW_RAD,     &
        ZUNDEF, ZUNDEF, ZUNDEF, PDIR_SW, PZENITH,   &
        PQ1_G,PQ2_G,PQ3_G,PQ4_G,PQ5_G,PQ6_G,PQ7_G,  &
        N,ZOPT,ZOPT_BODY,ZEB_G,ZAB_G, ZHB_G )

   do JJ = 1, N
      !========================================================
      ! COMPUTE THE MEAN RADIANT TEMPERATURES (K)
      !========================================================

      ! 1-calculation of mean radiant temperature values for a standard clothed standig human (eg, for UTCI)

      PTRAD_HSUN(JJ)     = MRT_BODY_SURF(ZEB_H,PQ2_H(JJ),PQ3_H(JJ),PQ4_H(JJ), &
           PQ5_H(JJ),PQ6_H(JJ),PQ7_H(JJ),PQ1_H(JJ) )
      PTRAD_HSHADE(JJ)   = MRT_BODY_SURF(ZEB_H,PQ2_H(JJ),PQ3_H(JJ),PQ4_H(JJ), &
           PQ5_H(JJ),PQ6_H(JJ),PQ7_H(JJ),ZUNDEF(JJ))

      ! 2-calculation of mean radiant temperature values for a black globe sensor (eg, for WBGT)

      PTRAD_GSUN(JJ)     = MRT_BODY_SURF(ZEB_G,PQ2_G(JJ),PQ3_G(JJ),PQ4_G(JJ),   &
           PQ5_G(JJ),PQ6_G(JJ),PQ7_G(JJ),PQ1_G(JJ))
      PTRAD_GSHADE(JJ)   = MRT_BODY_SURF(ZEB_G,PQ2_G(JJ),PQ3_G(JJ),PQ4_G(JJ),   &
           PQ5_G(JJ),PQ6_G(JJ),PQ7_G(JJ),ZUNDEF(JJ) )

      !========================================================
      ! COMPUTE THE GLOBE TEMPERATURE (K)
      !========================================================

      PTGLOBE_SUN(JJ)  = TGLOBE_BODY_SURF(PTRAD_GSUN(JJ), PTA(JJ), PUSR(JJ),    &
           ZGD, ZEB_G)
      PTGLOBE_SHADE(JJ)= TGLOBE_BODY_SURF(PTRAD_GSHADE(JJ), PTA(JJ), PUSR(JJ),   &
           ZGD, ZEB_G)

      !========================================================
      ! compute the (psychometric) wet-bulb temperatures (C)
      !========================================================

      PTWETB(JJ)       = WETBULBT_SURF(PPS(JJ), PTA(JJ)-TCDK, PQA(JJ))

      !========================================================
      ! compute the Universal Thermal and Climate Index UTCI (C)
      !========================================================

      ZEHPA(JJ) = PQA(JJ)* PPS(JJ)/ (0.622 + 0.378 * PQA(JJ)) /100.

      PUTCI_OUTSUN(JJ) = UTCI_APPROX_SURF(PTA(JJ) -TCDK, ZEHPA(JJ),       &
           PTRAD_HSUN(JJ) -TCDK, PU10(JJ))
      PUTCI_OUTSHADE(JJ) = UTCI_APPROX_SURF(PTA(JJ) -TCDK, ZEHPA(JJ),     &
           PTRAD_HSHADE(JJ)-TCDK, PU10(JJ) )

      !========================================================
      ! compute the wet bulb globe temperature indices  (WBGT) (C)
      !========================================================

      WBGT_SUN(JJ)   = 0.2 * (PTGLOBE_SUN(JJ) -TCDK)     +   &
           0.7 *  PTWETB(JJ)                +   &
           0.1 * (PTA(JJ)-TCDK)

      WBGT_SHADE(JJ) = 0.3 * (PTGLOBE_SHADE(JJ) -TCDK)   +   &
           0.7 * PTWETB(JJ)

      !========================================================
   enddo

end subroutine SURF_THERMAL_STRESS


!===================================================================
real function TGLOBE_BODY_SURF(ZTRAD,ZTA,ZUMOD,ZGD,ZGE)
use, intrinsic :: iso_fortran_env, only : REAL64
implicit none
!!!#include <arch_specific.hf>

   !  @Author : S. Leroyer (sept. 2014)
   !            MODIFICATIONS : S. Leroyer and V. Lee (sept 2020): optimization 
   
   !  Computes the black globe temperature that would be measured by a black globe sensor
   !  from the radiant temperature equivalent to the total radiation received by the human body
   !  Analytical solution of the equation:
   !  x**4 + a * x - b =0
   !  with  x=G_T
   !  b=Tmrt**4+a*A_T
   !  a= 1.335* 1E8 * va**0.71   / (em  *  D**0.4)
   !  4 solutions but 1 only in the desired range

   real :: ZTRAD,ZTA,ZUA,ZGD,ZGE
   real :: ZUMOD    ! bounded wind speed
   real(REAL64) :: ZWORKA ! Term for the resolution of the equation
   real(REAL64) :: ZWORKB ! Term for the resolution of the equation
   real(REAL64) :: ZWORKK ! Term for the resolution of the equation
   real(REAL64) :: ZWORKE ! Term for the resolution of the equation
   real :: ZWORKJ ! Term for the resolution of the equation
   real :: ZWORKI ! Term for the resolution of the equation
! optional for verification 4.
!!$   real :: ZTRAD2 ! body MRT for verification
!!$   real :: ZDIFF ! body MRT for verification
! real, save :: zcf = 1.73205     ! f = 3.0**0.5
! real, save :: zcg = 0.381571    ! g = ( 2.0**(1./3.) * 3.0**(2./3.) )**(-1.) 
! real, save :: zcq = 3.4943      ! cq = 4.0 * (2./3.)**(1./3.)

!*       1.    set limits for the  wind value
!               ---------------------------
!  note: some tests suggest a lower limit of 0.001 
   ZUA = min(max(0.2,ZUMOD),15.0) 

   !*       2.    convection coefficient
   !               ---------------------------
! convection coeff set to =  1.100 * 1.0E8 * U**0.6    ISO-ASHRAE - Kuehn et al. 1970
   ZWORKA =  1.1* 1.0E8 * ZUA **0.6   / (ZGE  *  ZGD**0.4)

   !*       3.    analytic resolution
   !               -------------------
   ZWORKB = ZTRAD **4. + ZWORKA * ZTA

 ZWORKE= ( (9.d0*ZWORKA**2.) +  1.73205d0* &
             ( (27.d0*ZWORKA**4.)+(256.d0*ZWORKB**3.) )**0.5 )**(1./3.)
 ZWORKK=ZWORKE*0.381571d0 - ((3.4943d0 * ZWORKB)/ZWORKE)
 ZWORKI= SNGL(0.5d0 *  ( 2.0d0 * ZWORKA /  ZWORKK**0.5 - ZWORKK)**0.5)
 ZWORKJ= SNGL(0.5d0 * ZWORKK**0.5)

 TGLOBE_BODY_SURF = -1.0 * ZWORKJ + ZWORKI

   !*       4.    optional : verification (reciprocity)
   !               ---------------------------

!   ZTRAD2 = ( TGLOBE_BODY_SURF**4.0 +   &
!        ZWORKa*(TGLOBE_BODY_SURF -ZTA) )**0.25

!   ZDIFF= abs(ZTRAD2 -ZTRAD)

!   if (ZDIFF .gt. 0.1) then
      !#TODO: use of print in the physics must be used only for debugging purpose at it renders listings files unusable (too big) in operational configs
!      print*,'** Tglobe diff is  /',ZDIFF ,'i tr t u ',ZTRAD,ZTA,ZUA
!   endif

   return
end function TGLOBE_BODY_SURF


!===================================================================
real function WETBULBT_SURF(PPA,PTA,PQA)
   use tdpack_const, only: TCDK
   implicit none
!!!#include <arch_specific.hf>

   !  Author : S. Leroyer, EC (sept. 2014)
   !  Computes the psychometric Wet-Bulb Temperature in K (diff with natural ?)


   ! !*      0.1    declarations of arguments
   real  :: PPA          ! Air pressure (Pa)
   real  :: PTA          ! Air temperature (C)
   real  :: PQA          ! Air specific humidity (kg/kg)
   ! !*      0.2    declarations of local variables
   real :: ZRH      ! Relative humidity (%)
   real :: ZTD      ! Dew-point temperature (C) at T,P
   real :: ZTI      ! temperature (C) on dry adiabat at P
   real :: ZP       ! pressure (hPa)
   real :: ZPI      ! pressure (hPa) intersec. dry adiabat/mixing ratio curves
   real :: ZAW      ! Mixing ratio (g/kg) at saturation
   real :: ZAO      ! Dry adiabat (C) at T,P (Theta)
   real :: ZAOS     ! Saturation dry adiabat (Theta Equivalent)
   real :: ZDIFF    ! absolute temperature difference
   real :: qsat
   integer :: JITER

   ! EXTERNAL FUNCTIONS (see wet_bulb.cdk90)
   real, external :: TPDD
   real, external :: N_QSAT
   real, external :: DWPT
   real, external :: W
   real, external :: O
   real, external :: TMR
   real, external :: TDA
   real, external :: OS
   real, external :: TSA

   !    1. Compute the dew-point temperature ZTD in C (INPUT: ZT,ZRH)

   ZP=0.01*PPA         ! convert in hPa
   qsat=N_QSAT(PTA,PPA)
   ZRH=100.*max(min(PQA,qsat),0.)/qsat

   ZTD=DWPT(PTA,ZRH)
   !    2. Compute the saturation mixing ratio in g/kg
   ZAW = W(ZTD,ZP)
   !    3. Compute the Dry Adiabat ZAO in C (INPUT:ZT,ZP)
   ZAO = O(PTA,ZP)
   ZPI = ZP
   !    4. Iterate to determine local pressure ZPI (hPa)
   !       at the intersection of the two curves of saturation mixing ratio
   !       and dry adiabat (first guess is ZP)

   do JITER= 1,10
      ZDIFF= 0.02*(TMR(ZAW,ZPI)-TDA(ZAO,ZPI))
      if (abs(ZDIFF).lt.0.01) GO TO 5
      ZPI= ZPI*(2.**ZDIFF)
   enddo
   !    5. Compute temperature on the dry adiabat at ZPI  (Intesrsection point)
5  ZTI= TDA(ZAO,ZPI)

   !    6. Compute a saturation adiabat from the intersection
   !       ==> equivalent potential temperature of a parcel saturated
   !           at temperature ZTI and pressure ZPI
   ZAOS= OS(ZTI,ZPI)

   !    7. Compute the Wet-bulb temperature (K) of a parcel at ZP given its
   !       equivalent potential temperature
   WETBULBT_SURF= TSA(ZAOS,ZP)

   return
end function WETBULBT_SURF


!===================================================================
real function UTCI_APPROX_SURF(PTA,PEHPA,PTMRT,PVA)
   implicit none
!!!#include <arch_specific.hf>
   real :: PTA, PEHPA, PTMRT, PVA

   !local variables
   real :: ZPA, ZD_TMRT

   real, dimension(7,7,7,7) :: ZZ

   real :: Z0, Z1, Z2, Z3, Z4, Z5, Z6, ZF, ZS
   ZS(Z0,Z1,Z2,Z3,Z4,Z5,Z6,ZF) = Z0 + Z1*ZF + Z2*ZF**2 + Z3*ZF**3 + Z4*ZF**4 + Z5*ZF**5 + Z6*ZF**6
   real, dimension(7) :: ZC_TA, ZC_VA, ZC_TMRT, ZC_PA
   integer :: J1, J2, J3, J4

   ZZ(:,:,:,:)=0.
   ! va
   !                  Ta**0,         Ta**1,           Ta**2,         Ta**3,          Ta**4,           Ta**5,          Ta**6
   ZZ(1,1,1,1:7) = (/ 6.07562052D-01,-2.27712343D-02, 8.06470249D-04,-1.54271372D-04,-3.24651735D-06,&
        7.32602852D-08, 1.35959073D-09/) ! va**0
   ZZ(1,1,2,1:6) = (/-2.25836520D+00, 8.80326035D-02, 2.16844454D-03,-1.53347087D-05,-5.72983704D-07,&
        -2.55090145D-09/) ! va**1
   ZZ(1,1,3,1:5) = (/-7.51269505D-01,-4.08350271D-03,-5.21670675D-05, 1.94544667D-06, 1.14099531D-08/) ! va**2
   ZZ(1,1,4,1:4) = (/ 1.58137256D-01,-6.57263143D-05, 2.22697524D-07,-4.16117031D-08/)                 ! va**3
   ZZ(1,1,5,1:3) = (/-1.27762753D-02, 9.66891875D-06, 2.52785852D-09/)                                 ! va**4
   ZZ(1,1,6,1:2) = (/ 4.56306672D-04,-1.74202546D-07/)                                                 ! va**5
   ZZ(1,1,7,1)   = -5.91491269D-06                                                                     ! va**6
   ! D_Tmrt / va
   ZZ(1,2,1,1:6) = (/ 3.98374029D-01, 1.83945314D-04,-1.73754510D-04,-7.60781159D-07, 3.77830287D-08, 5.43079673D-10/)
   ZZ(1,2,2,1:5) = (/-2.00518269D-02, 8.92859837D-04, 3.45433048D-06,-3.77925774D-07,-1.69699377D-09/)
   ZZ(1,2,3,1:4) = (/ 1.69992415D-04,-4.99204314D-05, 2.47417178D-07, 1.07596466D-08/)
   ZZ(1,2,4,1:3) = (/ 8.49242932D-05, 1.35191328D-06,-6.21531254D-09/)
   ZZ(1,2,5,1:2) = (/-4.99410301D-06,-1.89489258D-08/)
   ZZ(1,2,6,1)   = 8.15300114D-08
   ! D_Tmrt**2 / va
   ZZ(1,3,1,1:5) = (/ 7.55043090D-04,-5.65095215D-05,-4.52166564D-07, 2.46688878D-08, 2.42674348D-10/)
   ZZ(1,3,2,1:4) = (/ 1.54547250D-04, 5.24110970D-06,-8.75874982D-08,-1.50743064D-09/)
   ZZ(1,3,3,1:3) = (/-1.56236307D-05,-1.33895614D-07, 2.49709824D-09/)
   ZZ(1,3,4,1:2) = (/ 6.51711721D-07, 1.94960053D-09/)
   ZZ(1,3,5,1)   = -1.00361113D-08
   !D_Tmrt**3 / va
   ZZ(1,4,1,1:4) = (/-1.21206673D-05,-2.18203660D-07, 7.51269482D-09, 9.79063848D-11/)
   ZZ(1,4,2,1:3) = (/ 1.25006734D-06,-1.81584736D-09,-3.52197671D-10/)
   ZZ(1,4,3,1:2) = (/-3.36514630D-08, 1.35908359D-10/)
   ZZ(1,4,4,1)   = 4.17032620D-10
   !D_Tmrt**4 / va
   ZZ(1,5,1,1:3) = (/-1.30369025D-09, 4.13908461D-10, 9.22652254D-12/)
   ZZ(1,5,2,1:2) = (/-5.08220384D-09,-2.24730961D-11/)
   ZZ(1,5,3,1)   = 1.17139133D-10
   !D_Tmrt**5 / va
   ZZ(1,6,1,1:2) = (/6.62154879D-10, 4.03863260D-13/)
   ZZ(1,6,2,1)   = 1.95087203D-12
   !D_Tmrt**6
   ZZ(1,7,1,1)   = -4.73602469D-12
   ! Pa / va
   ZZ(2,1,1,1:6) = (/ 5.12733497D+00,-3.12788561D-01,-1.96701861D-02, 9.99690870D-04, 9.51738512D-06,-4.66426341D-07/)
   ZZ(2,1,2,1:5) = (/ 5.48050612D-01,-3.30552823D-03,-1.64119440D-03,-5.16670694D-06, 9.52692432D-07/)
   ZZ(2,1,3,1:4) = (/-4.29223622D-02, 5.00845667D-03, 1.00601257D-06,-1.81748644D-06/)
   ZZ(2,1,4,1:3) = (/-1.25813502D-03,-1.79330391D-04, 2.34994441D-06/)
   ZZ(2,1,5,1:2) = (/ 1.29735808D-04, 1.29064870D-06/)
   ZZ(2,1,6,1)   = -2.28558686D-06
   ! Pa / D_Tmrt / va
   ZZ(2,2,1,1:5) = (/-3.69476348D-02, 1.62325322D-03,-3.14279680D-05, 2.59835559D-06,-4.77136523D-08/)
   ZZ(2,2,2,1:4) = (/ 8.64203390D-03,-6.87405181D-04,-9.13863872D-06, 5.15916806D-07/)
   ZZ(2,2,3,1:3) = (/-3.59217476D-05, 3.28696511D-05,-7.10542454D-07/)
   ZZ(2,2,4,1:2) = (/-1.24382300D-05,-7.38584400D-09/)
   ZZ(2,2,5,1)   = 2.20609296D-07
   ! Pa / D_Tmrt**2 / va
   ZZ(2,3,1,1:4) = (/-7.32469180D-04,-1.87381964D-05, 4.80925239D-06,-8.75492040D-08/)
   ZZ(2,3,2,1:3) = (/ 2.77862930D-05,-5.06004592D-06, 1.14325367D-07/)
   ZZ(2,3,3,1:2) = (/ 2.53016723D-06,-1.72857035D-08/)
   ZZ(2,3,4,1)   = -3.95079398D-08
   ! Pa / D_Tmrt**3 / va
   ZZ(2,4,1,1:3) = (/-3.59413173D-07, 7.04388046D-07,-1.89309167D-08/)
   ZZ(2,4,2,1:2) = (/-4.79768731D-07, 7.96079978D-09/)
   ZZ(2,4,3,1)   = 1.62897058D-09
   ! Pa / D_Tmrt**4 / va
   ZZ(2,5,1,1:2) = (/ 3.94367674D-08,-1.18566247D-09/)
   ZZ(2,5,2,1)   = 3.34678041D-10
   ! Pa / D_Tmrt**5
   ZZ(2,6,1,1)   = -1.15606447D-10
   ! Pa**2 / va
   ZZ(3,1,1,1:5) = (/-2.80626406D+00, 5.48712484D-01,-3.99428410D-03,-9.54009191D-04, 1.93090978D-05/)
   ZZ(3,1,2,1:4) = (/-3.08806365D-01, 1.16952364D-02, 4.95271903D-04,-1.90710882D-05/)
   ZZ(3,1,3,1:3) = (/ 2.10787756D-03,-6.98445738D-04, 2.30109073D-05/)
   ZZ(3,1,4,1:2) = (/ 4.17856590D-04,-1.27043871D-05/)
   ZZ(3,1,5,1)   = -3.04620472D-06
   ! Pa**2 / D_Tmrt / va
   ZZ(3,2,1,1:4) = (/ 5.14507424D-02,-4.32510997D-03, 8.99281156D-05,-7.14663943D-07/)
   ZZ(3,2,2,1:3) = (/-2.66016305D-04, 2.63789586D-04,-7.01199003D-06/)
   ZZ(3,2,3,1:2) = (/-1.06823306D-04, 3.61341136D-06/)
   ZZ(3,2,4,1)   = 2.29748967D-07
   ! Pa**2 / D_Tmrt**2 / va
   ZZ(3,3,1,1:3) = (/3.04788893D-04,-6.42070836D-05, 1.16257971D-06/)
   ZZ(3,3,2,1:2) = (/7.68023384D-06,-5.47446896D-07/)
   ZZ(3,3,3,1)   = -3.59937910D-08
   ! Pa**2 / D_Tmrt**3 / va
   ZZ(3,4,1,1:2) = (/-4.36497725D-06, 1.68737969D-07/)
   ZZ(3,4,2,1)   = 2.67489271D-08
   ! Pa**2 / D_Tmrt**4
   ZZ(3,5,1,1) =  3.23926897D-09
   ! Pa**3 / va
   ZZ(4,1,1,1:4) = (/-3.53874123D-02, -2.21201190D-01, 1.55126038D-02, -2.63917279D-04/)
   ZZ(4,1,2,1:3) = (/4.53433455D-02, -4.32943862D-03, 1.45389826D-04/)
   ZZ(4,1,3,1:2) = (/2.17508610D-04, -6.66724702D-05/)
   ZZ(4,1,4,1)   = 3.33217140D-05
   ! Pa**3 / D_Tmrt / va
   ZZ(4,2,1,1:3) = (/-2.26921615D-03, 3.80261982D-04, -5.45314314D-09/)
   ZZ(4,2,2,1:2) = (/-7.96355448D-04, 2.53458034D-05/)
   ZZ(4,2,3,1)   = -6.31223658D-06
   ! Pa**3 / D_Tmrt**2 / va
   ZZ(4,3,1,1:2) = (/3.02122035D-04, -4.77403547D-06/)
   ZZ(4,3,2,1)   = 1.73825715D-06
   ! Pa**3 / D_Tmrt**3
   ZZ(4,4,1,1) = -4.09087898D-07
   !  Pa**4 / va
   ZZ(5,1,1,1:3) = (/6.14155345D-01, -6.16755931D-02, 1.33374846D-03/)
   ZZ(5,1,2,1:2) = (/3.55375387D-03, -5.13027851D-04/)
   ZZ(5,1,3,1)   = 1.02449757D-04
   !  Pa**4 / D_Tmrt / va
   ZZ(5,2,1,1:2) = (/-1.48526421D-03, -4.11469183D-05/)
   ZZ(5,2,2,1)   = -6.80434415D-06
   !  Pa**4 / D_Tmrt**2 / va
   ZZ(5,3,1,1) = -9.77675906D-06
   !  Pa**5 / va
   ZZ(6,1,1,1:2) = (/8.82773108D-02, -3.01859306D-03/)
   ZZ(6,1,2,1)   = 1.04452989D-03
   !  Pa**5 / D_Tmrt
   ZZ(6,2,1,1) = 2.47090539D-04
   ! Pa**6
   ZZ(7,1,1,1) = 1.48348065D-03

   ZD_TMRT = PTMRT - PTA
   ZPA = PEHPA / 10.0; !~ use vapour pressure in kPa

   ZC_TA(:) = 0.
   ZC_VA(:) = 0.
   ZC_TMRT(:) = 0.
   ZC_PA(:) = 0.

   do J4 = 1,7
      do J3 = 1,7
         do J2 = 1,7
            do J1 = 1,7
               ZC_TA(J1) = ZZ(J4,J3,J2,J1)
            enddo
            !ZC_VA(J2) = ZC_TA(1)+ZC_TA(2)*PTA+ZC_TA(3)*PTA**2+ZC_TA(4)*PTA**3+ZC_TA(5)*PTA**4 &
            !               +ZC_TA(6)*PTA**5+ZC_TA(7)*PTA**6
            ZC_VA(J2) = ZS(ZC_TA(1),ZC_TA(2),ZC_TA(3),ZC_TA(4),ZC_TA(5),ZC_TA(6),ZC_TA(7),PTA)
         enddo
         !ZC_TMRT(J3) = ZC_VA(1)+ZC_VA(2)*PVA+ZC_VA(3)*PVA**2+ZC_VA(4)*PVA**3+ZC_VA(5)*PVA**4 &
         !               +ZC_VA(6)*PVA**5+ZC_VA(7)*PVA**6
         ZC_TMRT(J3) = ZS(ZC_VA(1),ZC_VA(2),ZC_VA(3),ZC_VA(4),ZC_VA(5),ZC_VA(6),ZC_VA(7),PVA)
      enddo
      !ZC_PA(J4) = ZC_TMRT(1)+ZC_TMRT(2)*ZD_TMRT+ZC_TMRT(3)*ZD_TMRT**2+ZC_TMRT(4)*ZD_TMRT**3 &
      !               +ZC_TMRT(5)*ZD_TMRT**4 +ZC_TMRT(6)*ZD_TMRT**5+ZC_TMRT(7)*ZD_TMRT**6
      ZC_PA(J4) = ZS(ZC_TMRT(1),ZC_TMRT(2),ZC_TMRT(3),ZC_TMRT(4),ZC_TMRT(5),ZC_TMRT(6),ZC_TMRT(7),ZD_TMRT)
   enddo

   !PUTCI_APPROX = PTA + ZC_PA(1)+ZC_PA(2)*ZPA+ZC_PA(3)*ZPA**2+ZC_PA(4)*ZPA**3 &
   !                  +ZC_PA(5)*ZPA**4 + ZC_PA(6)*ZPA**5+ZC_PA(7)*ZPA**6
   UTCI_APPROX_SURF = PTA + ZS(ZC_PA(1),ZC_PA(2),ZC_PA(3),ZC_PA(4),ZC_PA(5),ZC_PA(6),ZC_PA(7),ZPA)


   return
end function UTCI_APPROX_SURF


!===================================================================
real function MRT_BODY_SURF(PEMISS_BODY,PQ2,PQ3,PQ4,  &
     PQ5,PQ6,PQ7,PQ1)
   use tdpack_const, only: STEFAN
   implicit none
!!!#include <arch_specific.hf>

   !    PURPOSE
   !    Computes the mean radiant temperature from the energy budget received by a body

   !* --  declarations of arguments
   real :: PQ1
   real :: PQ2
   real :: PQ3
   real :: PQ4
   real :: PQ5
   real :: PQ6
   real :: PQ7
   real :: ZQSHADE
   real :: PEMISS_BODY

   ! -- energy budget in the shade
   ZQSHADE = PQ2+PQ3+PQ4+PQ5+PQ6+PQ7

   ! -- MRT
   MRT_BODY_SURF = ((PQ1+ZQSHADE)/PEMISS_BODY/STEFAN)**0.25

   return
end function MRT_BODY_SURF
