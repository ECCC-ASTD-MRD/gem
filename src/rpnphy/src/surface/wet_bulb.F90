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

module MODI_WETBULBT

   !#TODO: no need for an interface here, put the funtion WETBULBT in the
   !       module to avoid non consistency between module and function interface

   interface
      function WETBULBT(PPA, PTA, PQA) result(PTWBT)
         real, dimension(:), intent(IN)    :: PPA     ! pressure (Pa)
         real, dimension(:), intent(IN)    :: PTA     ! Air  temperature (C)
         real, dimension(:), intent(IN)    :: PQA     ! Air spedcific humidity (kg/kg)
         real, dimension(size(PPA)) :: PTWBT          ! Wet-Bulb  temperature in deg C
      end function WETBULBT
   end interface
end module MODI_WETBULBT


!-----------------------------------------------------------------------------
function WETBULBT(PPA, PTA, PQA) result(PTWBT)
   !---------------------------------------------------------------------------
   !    PURPOSE       : Computes the Natural Wet-Bulb Temperature
   !    AUTHOR        : S. Leroyer (03/2014)
   !    REFERENCE     : Stipanuk (1973): ALGORITHMS FOR GENERATING A SKEW-T, LOG P
   !                    DIAGRAM AND COMPUTING SELECTED METEOROLOGICAL QUANTITIES."
   !                    ATMOSPHERIC SCIENCES LABORATORY,U.S. ARMY ELECTRONICS COMMAND
   !                    WHITE SANDS MISSILE RANGE, NEW MEXICO 88002, 33 PAGES
   !    MODIFICATIONS : L. Spacek (06/2014)
   !    MODIFICATIONS : S. Leroyer et V. Vionnet (09/2020) : Temperature in Celcius
   !---------------------------------------------------------------------------
   !
   !*       0.     DECLARATIONS
   !               ------------
   !
   use MODD_CSTS
   ! USE MODE_THERMOS
   ! !
   implicit none
!!!#include <arch_specific.hf>
   ! !
   ! !*      0.1    declarations of arguments
   real, dimension(:), intent(IN)  :: PPA         ! Air pressure (Pa)
   real, dimension(:), intent(IN)  :: PTA         ! Air temperature (C)
   real, dimension(:), intent(IN)  :: PQA         ! Air specific humidity (kg/kg)
   real, dimension(size(PPA))      :: PTWBT       ! Wet-bulb temperature in deg C
   ! !*      0.2    declarations of local variables
   real, dimension(size(PPA)) :: ZQA      ! modified Air specific humidity (kg/kg)
   real, dimension(size(PPA)) :: ZRH      ! Relative humidity (%)
   real, dimension(size(PPA)) :: ZTD      ! Dew-point temperature (C) at T,P
   real, dimension(size(PPA)) :: ZTI      ! temperature (C) on dry adiabat at P
   real, dimension(size(PPA)) :: ZP       ! pressure (hPa)
   real, dimension(size(PPA)) :: ZPI      ! pressure (hPa) intersec. dry adiabat/mixing ratio curves
   real, dimension(size(PPA)) :: ZAW      ! Mixing ratio (g/kg) at saturation
   real, dimension(size(PPA)) :: ZAO      ! Dry adiabat (C) at T,P (Theta)
   real, dimension(size(PPA)) :: ZAOS     ! Saturation dry adiabat (Theta Equivalent)
   real, dimension(size(PPA)) :: ZDIFF    ! absolute temperature difference
   real                       :: qsat
   integer :: JITER
   integer :: JJ,i,nn

   ! EXTERNAL FUNCTIONS
   real, external :: TPDD
   real, external :: N_QSAT
   real, external :: DWPT
   real, external :: W
   real, external :: O
   real, external :: TMR
   real, external :: TDA
   real, external :: OS
   real, external :: TSA

   nn=size(ppa)

   !    1. Compute the dew-point temperature ZTD in C (INPUT: PTA,ZRH)

   ZP=0.01*PPA         ! convert in hPa
   ZQA=PQA             ! local variable
   do i=1,nn
      qsat=n_qsat(PTA(i),PPA(i))
      if(ZQA(i)<0)ZQA(i)=0.0
      if(qsat<ZQA(i))ZQA(i)=qsat
      ZRH(i)=100.*ZQA(i)/qsat
      !?! ZRH(i)=100.*max(min(ZQA(i),qsat),0.)/qsat

      ZTD(i)=DWPT(PTA(i),ZRH(i))
      !    2. Compute the saturation mixing ratio in g/kg
      ZAW(i) = W(ZTD(i),ZP(i))
      !    3. Compute the Dry Adiabat ZAO in C (INPUT:PTA,ZP)
      ZAO(i) = O(PTA(i),ZP(i))
      ZPI(i) = ZP(i)
   enddo
   !    4. Iterate to determine local pressure ZPI (hPa)
   !       at the intersection of the two curves of saturation mixing ratio
   !       and dry adiabat (first guess is ZP)

   do JJ=1,size(PPA)
      do JITER= 1,10
         ZDIFF(JJ)= 0.02*(TMR(ZAW(JJ),ZPI(JJ))-TDA(ZAO(JJ),ZPI(JJ)))
         if (abs(ZDIFF(JJ)).lt.0.01) GO TO 5
         ZPI(JJ)= ZPI(JJ)*(2.**(ZDIFF(JJ)))
      enddo
      !    5. Compute temperature on the dry adiabat at ZPI  (Intesrsection point)
5     ZTI(JJ)= TDA(ZAO(JJ),ZPI(JJ))
   enddo

   !    6. Compute a saturation adiabat from the intersection
   !       ==> equivalent potential temperature of a parcel saturated
   !           at temperature ZTI and pressure ZPI
   do i=1,nn
      ZAOS(i)= OS(ZTI(i),ZPI(i))
      ! print*,'SATURATION DRY ADIABAT thetaEq OS(TI,PI) => ZAOS= ',ZAOS

      !    7. Compute the Wet-bulb temperature (C) of a parcel at ZP given its
      !       equivalent potential temperature
      PTWBT(i)= TSA(ZAOS(i),ZP(i))
   enddo
   return
end function WETBULBT


!!==============================================
function n_qsat(T,P)
   ! THIS FUNCTION RETURNS HUMIDITY AT SATURATION (KG/KG) GIVEN
   ! THE TEMPERATURE (C) AND SPECIFIC HUMIDITY (KG/KG)
   ! AND PRESSURE (Pa)
   implicit none
!!!#include <arch_specific.hf>
   real T,P
   real zeps,zfoes,n_qsat
   real TA

   real, save :: XRD    = 287.05967   ! gaz constant for dry air
   real, save :: XRV    = 461.524993  ! gaz constant for water vapor
   real, save :: XALPW=60.22416    ! constants in saturation pressure over liquid water
   real, save :: XBETAW=6822.459384
   real, save :: XGAMW=5.13948

   TA = T + 273.15 ! Convert T in K
   
   ZEPS      = XRD / XRV
   !*       1.    COMPUTE SATURATION VAPOR PRESSURE
   ZFOES = exp( XALPW - XBETAW/TA - XGAMW*log(TA)  )
   !*       2.    COMPUTE SATURATION HUMIDITY
   N_QSAT = ZEPS*ZFOES/P   &
        / (1.+(ZEPS-1.)*ZFOES/P)
   return
end function N_qsat


!==============================================
real function DWPT(T,RH)
   implicit none
!!!#include <arch_specific.hf>
   real :: t,rh,x,dpd
   !    THIS FUNCTION RETURNS THE DEW POINT (CELSIUS) GIVEN THE TEMPERATURE
   !    (CELSIUS) AND RELATIVE HUMIDITY (%). T
   X = 1.-0.01*RH
   !    COMPUTE DEW POINT DEPRESSION.
   DPD =(14.55+0.114*T)*X+((2.5+0.007*T)*X)**3+(15.9+0.117*T)*X**14
   DWPT = T-DPD

   return
end function DWPT


!===================================================================
real function TSA(OS,P)
   implicit none
!!!#include <arch_specific.hf>
   !    THIS FUNCTION RETURNS THE TEMPERATURE TSA (CELSIUS) ON A SATURATION
   !    ADIABAT AT PRESSURE P (MILLIBARS). OS IS THE EQUIVALENT POTENTIAL
   !    TEMPERATURE OF THE PARCEL (CELSIUS). SIGN(A,B) REPLACES THE
   !    ALGEBRAIC SIGN OF A WITH THAT OF B.
   !    B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO 0.001 OF THE LATENT
   !    HEAT OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
   !    PRESSURE FOR DRY AIR.
   integer :: i
   real :: os,p,a,b,tq,d,tqk,x
   real, external :: w
   data B/2.6518986/
   A= OS+273.15

   !    TQ IS THE FIRST GUESS FOR TSA.

   TQ= 253.15

   !    D IS AN INITIAL VALUE USED IN THE ITERATION BELOW.

   D= 120.

   !    ITERATE TO OBTAIN SUFFICIENT ACCURACY....SEE TABLE 1, P.8
   !    OF STIPANUK (1973) FOR EQUATION USED IN ITERATION.

   do I= 1,12
      TQK= TQ-273.15
      D= D/2.
      X= A*exp(-B*W(TQK,P)/TQ)-TQ*((1000./P)**.286)
      if (abs(X).lt.1E-7) GO TO 2
      TQ= TQ+sign(D,X)
   enddo
2  TSA= TQ-273.15
   return
end function TSA


!===================================================================
real function W(T,P)
   implicit none
!!!#include <arch_specific.hf>
   real :: t,p,x
   real, external :: esat
   !   THIS FUNCTION RETURNS THE MIXING RATIO (GRAMS OF WATER VAPOR PER
   !   KILOGRAM OF DRY AIR) GIVEN THE DEW POINT (CELSIUS) AND PRESSURE
   !   (MILLIBARS). IF THE TEMPERTURE  IS INPUT INSTEAD OF THE
   !   DEW POINT, THEN SATURATION MIXING RATIO (SAME UNITS) IS RETURNED.
   !   THE FORMULA IS FOUND IN MOST METEOROLOGICAL TEXTS.
   X= ESAT(T)
   X= min( X , 0.5*P )
   W= 622.*X/(P-X)
   W=min( W , 50.)
   W=max( W , 1.E-3)
   return
end function W


!===================================================================
real function O(T,P)
   implicit none
!!!#include <arch_specific.hf>
   !    THIS FUNCTION RETURNS POTENTIAL TEMPERATURE (CELSIUS) GIVEN
   !    TEMPERATURE T (CELSIUS) AND PRESSURE P (MB) BY SOLVING THE POISSON
   !    EQUATION.
   real :: t,p,tk,ok
   TK= T+273.15
   OK= TK*((1000./P)**.286)
   O= OK-273.15
   return
end function O


!===================================================================
real function TDA(O,P)
   implicit none
!!!#include <arch_specific.hf>
   real :: o,p,ok,tdak
   !!    THIS FUNCTION RETURNS THE TEMPERATURE TDA (CELSIUS) ON A DRY ADIABAT
   !    AT PRESSURE P (MILLIBARS). THE DRY ADIABAT IS GIVEN BY
   !    POTENTIAL TEMPERATURE O (CELSIUS). THE COMPUTATION IS BASED ON
   !    POISSON'S EQUATION.

   OK= O+273.15
   TDAK= OK*((P*.001)**.286)
   TDA= TDAK-273.15
   return
end function TDA


!===================================================================
real function OS(T,P)
   implicit none
!!!#include <arch_specific.hf>
   real :: t,p,b,tk,osk
   real, external :: w
   !    THIS FUNCTION RETURNS THE EQUIVALENT POTENTIAL TEMPERATURE OS
   !    (CELSIUS) FOR A PARCEL OF AIR SATURATED AT TEMPERATURE T (CELSIUS)
   !    AND PRESSURE P (MILLIBARS).
   data B/2.6518986/
   !    B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO THE LATENT HEAT
   !    OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
   !    PRESSURE FOR DRY AIR.

   TK = T+273.15
   OSK= TK*((1000./P)**.286)*(exp(B*W(T,P)/TK))
   OS= OSK-273.15
   return
end function OS


!===================================================================
real function TMR(W,P)
   implicit none
!!!#include <arch_specific.hf>
   real :: w,p,c1,c2,c3,c4,c5,c6,x,tmrk
   !   THIS FUNCTION RETURNS THE TEMPERATURE (CELSIUS) ON A MIXING
   !   RATIO LINE W (G/KG) AT PRESSURE P (MB). THE FORMULA IS GIVEN IN
   !   TABLE 1 ON PAGE 7 OF STIPANUK (1973).
   !
   !   INITIALIZE CONSTANTS

   data C1/.0498646455/,C2/2.4082965/,C3/7.07475/
   data C4/38.9114/,C5/.0915/,C6/1.2035/

   X= ALOG10(W*P/(622.+W))
   TMRK= 10.**(C1*X+C2)-C3+C4*((10.**(C5*X)-C6)**2.)
   TMR= TMRK-273.15
   return
end function TMR


!===================================================================
real function ESAT(T)
   implicit none
!!!#include <arch_specific.hf>
   real :: t,tk,p1,p2,c1
   !   THIS FUNCTION RETURNS THE SATURATION VAPOR PRESSURE OVER
   !   WATER (MB) GIVEN THE TEMPERATURE (CELSIUS).
   !   THE ALGORITHM IS DUE TO NORDQUIST, W.S.,1973:

   TK = T+273.15
   P1 = 11.344-0.0303998*TK
   P2 = 3.49149-1302.8844/TK
   C1 = 23.832241-5.02808*ALOG10(TK)
   ESAT = 10.**(C1-1.3816E-7*10.**P1+8.1328E-3*10.**P2-2949.076/TK)
   return
end function ESAT
