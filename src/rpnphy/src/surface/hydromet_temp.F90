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

module MODI_HYDROMET_TEMP

   !#TODO: no need for an interface here, put the funtion WETBULBT in the module to avoid non consistency between module and function interface

   interface
      function HYDROMET_TEMP(PPA, PTA, PQA) result(PTI)
         real, dimension(:), intent(IN)    :: PPA     ! pressure (Pa)
         real, dimension(:), intent(IN)    :: PTA     ! Air  temperature (deg C) 
         real, dimension(:), intent(IN)    :: PQA     ! Air spedcific humidity (kg/kg)
         real, dimension(size(PPA)) :: PTI          ! Hydrometeo temperature in deg C
      end function HYDROMET_TEMP
   end interface
end module MODI_HYDROMET_TEMP


!-----------------------------------------------------------------------------
function HYDROMET_TEMP(PPA, PTA, PQA) result(PTI)
   !-----------------------------------------------------------------------------
   !    PURPOSE       : Computes the temperature of a falling hydrometeor
   !    AUTHOR        : V. Vionnet (11/2020)
   !    REFERENCE     :  Harder, P., Pomeroy, J. (2013).
   !                     Estimating precipitation phase using a psychrometric energy balance method
   !                    Hydrological Processes  27(13), 1901-1914. https://dx.doi.org/10.1002/hyp.9799
   !-------------------------------------------------------------------------------
   !
   !*       0.     DECLARATIONS
   !               ------------
   !
   !use MODD_CSTS
   ! USE MODE_THERMOS
   ! !
   implicit none
#include <arch_specific.hf>
   ! !
   ! !*      0.1    declarations of arguments
   real, dimension(:), intent(IN)  :: PPA         ! Air pressure (Pa)
   real, dimension(:), intent(IN)  :: PTA          ! Air temperature (deg C) 
   real, dimension(:), intent(IN)  :: PQA          ! Air specific humidity (kg/kg) 
   real, dimension(size(PPA))      :: PTI         ! Hydrometeo temprtature in deg C
   ! !*      0.2    declarations of local variables
   real, dimension(size(PPA)) :: ZD    !diffusivity of water vapour in air [m^2 s-1]
   real, dimension(size(PPA)) :: ZLAMBDAT !thermal conductivity of air [J m^-1 s^-1 K^-1]
   real, dimension(size(PPA)) :: ZL    !latent heat of sublimation of vaporisation[J  kg^-1]
   real, dimension(size(PPA)) :: ZRHODA   !density of dry air [kg m-3]
   real, dimension(size(PPA)) :: ZRH   !relative humidity [-]

   real ZT,ZTINI,ZF,ZFDIFF
   integer :: JITER
   integer :: JJ,I,NN

   ! EXTERNAL FUNCTIONS
   real, external :: DVAP
   real, external :: EVSAT
   real, external :: LAMBDA
   real, external :: LS
   real, external :: RHODA
   real, external :: RHOVSAT
   real, external :: RHOVSAT_DIFF

   NN=SIZE(PPA)


 !  write(*,*) 'On entre'
   DO I=1,NN
      ! Compute diffusivity of water vapour in air 
      ZD(I) = DVAP(PTA(I))

      ! Compute thermal conductivity of air 
      ZLAMBDAT(I) = LAMBDA(PTA(I))

     ! Compute latent heat of sublimation or vaporisation (depending on air temperature)
      ZL(I) = LS(PTA(I))

     ! Compute density of dry air
      ZRHODA(I) = RHODA(PTA(I),PPA(I))

      ! Compute relative humidity
      ZRH(I) = 100. * PQA(I)*PPA(I)/0.622 /EVSAT(PTA(I)) 
   ENDDO


 !  write(*,*) 'OK calc thermo'

   ! 4.  Solve iteratively to get Ti
   ! using a Newton-Raphston approach
   
   DO I=1,NN
     ZT = PTA(I)  ! Value used the 1st guess in the integration
     DO JITER = 1,10

        ZTINI = ZT   ! 
        ZF = ZT - PTA(I) - ZD(I)*ZL(I)/ZLAMBDAT(I) * ( PQA(I)*ZRHODA(I) - RHOVSAT(ZT))
        ZFDIFF = 1 +  ZD(I)*ZL(I)/ZLAMBDAT(I) * RHOVSAT_DIFF(ZT)
        ZT = ZTINI - ZF/ZFDIFF
       ! write(*,*) JITER,ZT
        IF(ABS(ZT- ZTINI) .LT. 0.01) EXIT        
     ENDDO
     PTI(I) = ZT
    ! IF(ZT>-4 .AND. ZT< -3.) THEN
    !     write(*,*) I,PTA(I), PQA(I), PPA(I)
    !     write(*,*) ZRH(I),PTI(I)
    !     write(*,*) ZD(I),ZL(I),ZLAMBDAT(I),ZRHODA(I)
    ! ENDIF
   ENDDO
    
END FUNCTION HYDROMET_TEMP


!!============================================== 
FUNCTION DVAP(T) 
   ! Compute diffusivity of water vapour in air [m2 s-1]
   ! Following Thorpe and Mason (1966)
   implicit none
   REAL T ! Air temperature  (deg C)
   REAL DVAP
   DVAP = 2.063e-5 * ((T+273.15)/273.15)**1.75
   RETURN
END FUNCTION DVAP

!!============================================== 
FUNCTION LAMBDA(T)
   !  Compute thermal conductivity of air [J m-1 s-1 K-1]
   implicit none
   REAL T ! Air temperature  (deg C)
   REAL LAMBDA
   LAMBDA = 0.000063 * (T+273.15) + 0.00673
   RETURN
END FUNCTION LAMBDA

!!============================================== 
FUNCTION LS(T) 
   ! Compute latent heat of sublimation or vaporisation (depending on air temperature)
   implicit none
   REAL T ! Air temperature I (deg C)
   REAL LS
   IF(T <0.) THEN
      LS = 1000.0 * (2834.1 - 0.29 *T - 0.004*T**2.)
   ELSE
      LS = 1000.0 * (2501.0 - (2.361 * T))
   ENDIF
   RETURN
END FUNCTION LS

!!============================================== 
FUNCTION RHODA(T,P)
   ! Compute dry air density [kg m-3]
   implicit none
   REAL T ! Air temperature  (deg C)
   REAL P ! Pressure  (Pa)
   REAL RHODA

   REAL, SAVE :: XRDA = 287.04 ! Gas constant for dry air [J K-1 kg-1]

   RHODA = P/(XRDA*(T+273.15))
  
   RETURN
END FUNCTION RHODA

!!============================================== 
FUNCTION EVSAT(T)
   ! Compute saturated water vapour pressure [Pa]
   ! with respect to ice if  T<=0 degC
   ! with respect to liqui water if T> 0degC
   implicit none
   REAL T ! Air temperature  (deg C)

   REAL  EVSAT

   REAL, SAVE :: XAW = 17.27 ! Constant for esat with respect to liquid water
   REAL, SAVE :: XBW = 237.3 ! Constant for esat with respect to liquid water

   REAL, SAVE :: XAI = 21.87 ! Constant for esat with respect to ice
   REAL, SAVE :: XBI = 265.5 ! Constant for esat with respect to ice

   IF(T>0) THEN
       EVSAT = 611.0*EXP(XAW*T/(T+XBW))
   ELSE
       EVSAT = 611.0*EXP(XAI*T/(T+XBI))
   ENDIF

   RETURN
END FUNCTION EVSAT


!!============================================== 
FUNCTION RHOVSAT(T)
   ! Compute saturated water vapour density [kg m-3]
   implicit none
   REAL T ! Air temperature  (deg C)

   REAL RHOVSAT, ESAT 
   REAL, external :: EVSAT

   REAL, SAVE :: XRDV = 461.5 ! Gas constant for water vapor [J K-1 kg-1]

   ESAT = EVSAT(T) ! Saturated water vapour pressure

   RHOVSAT  = ESAT/(XRDV*(T+273.15)) ! Saturated water vapour density

   RETURN
END FUNCTION RHOVSAT

!!============================================== 
FUNCTION RHOVSAT_DIFF(T)
   ! Compute derivative of the function given saturated water vapour density [kg m-3]
   implicit none
   REAL T ! Air temperature  (deg C)

   REAL RHOVSAT_DIFF

   REAL, SAVE ::  XAW = 17.27 ! Constant for esat with respect to liquid water
   REAL, SAVE ::  XBW = 237.3 ! Constant for esat with respect to liquid water

   REAL, SAVE ::  XAI = 21.87 ! Constant for esat with respect to ice
   REAL, SAVE ::  XBI = 265.5 ! Constant for esat with respect to ice

   REAL, SAVE ::  XRDV = 461.5 ! Gas constant for water vapor [J K-1 kg-1]

   IF(T>0) THEN
      RHOVSAT_DIFF  = 611.0/(XRDV*(T+273.15)) * EXP(XAW*T/(T+XBW)) *  &
                      (-1/(T+273.15) + XAW*XBW/((T+XBW))**2.)
   ELSE
      RHOVSAT_DIFF  = 611.0/(XRDV*(T+273.15)) * EXP(XAI*T/(T+XBI)) *  &
                      (-1/(T+273.15) + XAI*XBI/((T+XBI))**2.)
   ENDIF 

   RETURN
END FUNCTION RHOVSAT_DIFF
