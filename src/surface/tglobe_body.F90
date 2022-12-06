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
module modi_tglobe_body 

!#TODO: no need for an interface here, put the funtion TGLOBE_BODY in the module to avoid non consistency between module and function interface

interface
function tglobe_body(ptrad, pta, pua, zgd, zge) result(ptglobe_body)
real, dimension(:), intent(in)  :: PTRAD       ! body MRT  (K)
real, dimension(:), intent(in)  :: PTA          !Air  temperature (K) 
real, dimension(:), intent(in)  :: PUA          ! wind speed at 10m (m/s)
real, INTENT(IN)  :: ZGD          ! black globe sensor diameter
real, INTENT(IN)  :: ZGE          ! black globe sensor emissivity 
real, dimension(size(ptrad)) :: PTGLOBE_BODY
end function tglobe_body
end interface
end module modi_tglobe_body

!-----------------------------------------------------------------------------
function tglobe_body(ptrad, pta, pua, zgd, zge) result(ptglobe_body)
!-----------------------------------------------------------------------------
!    PURPOSE       : Computes the black globe temperature that would be measured by a black globe sensor 
!                    from the radiant temperature equivalent to the total radiation received by the human body
!    AUTHOR        :  S. Leroyer   (Original  03/2014)
!    REFERENCE     :  Leroyer et al. (2018)
!    MODIFICATIONS :  
!    METHOD        :
!  Analytical solution of the equation: 
!  x**4 + a * x - b =0
!
!  with  x=G_T
!  b=Tmrt**4+a*A_T
!  a= 1.335* 1E8 * va**0.71   / (em  *  D**0.4)
!
!  4 solutions but 1 only in the desired range
!-----------------------------------------------------------------------------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
! INPUT:
! PTRAD      Mean Radiant Temperature   (K)
! A_T      Air temperature     (K)
! A_U      Wind speed          (m/s)
! fixed here
! ZGD      Globe diameter      (mm)
! ZGE      Globe emissivity    (-)
! OUTPUT
! G_T      Globe temperature   (K)
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!!$USE MODD_CSTS, ONLY : XTT
!
use, intrinsic :: iso_fortran_env, only : REAL64
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
real, dimension(:), intent(in)  :: PTRAD       ! body MRT  (K)
real, dimension(:), intent(in)  :: PTA          ! Air temperature (K) 
real, dimension(:), intent(in)  :: PUA          ! Air wind speed (m/s)
real, intent(in)  :: ZGD          ! black globe sensor diameter       (PanAm: 148 mm)
real, intent(in)  :: ZGE          ! black globe sensor emissivity     (PanAm: 0.902 )
 
real, dimension(size(ptrad))    :: PTGLOBE_BODY

!*      0.2    declarations of local variables
! real, save :: zcf = 1.73205     ! f = 3.0**0.5
! real, save :: zcg = 0.381571    ! g = ( 2.0**(1./3.) * 3.0**(2./3.) )**(-1.) 
! real, save :: zcq = 3.4943      ! cq = 4.0 * (2./3.)**(1./3.)

real, dimension(size(ptrad)) :: ZUA    ! bounded wind speed
real(REAL64), dimension(size(ptrad)) :: ZWORKA ! Term for the resolution of the equation
real(REAL64), dimension(size(ptrad)) :: ZWORKB ! Term for the resolution of the equation
!real, dimension(size(ptrad)) :: ZWORKM ! Term for the resolution of the equation
!real, dimension(size(ptrad)) :: ZWORKN ! Term for the resolution of the equation
!real, dimension(size(ptrad)) :: ZWORKP ! Term for the resolution of the equation
!real, dimension(size(ptrad)) :: ZWORKQ ! Term for the resolution of the equation
real(REAL64), dimension(size(ptrad)) :: ZWORKK ! Term for the resolution of the equation
real(REAL64), dimension(size(ptrad)) :: ZWORKE ! Term for the resolution of the equation
real, dimension(size(ptrad)) :: ZWORKJ ! Term for the resolution of the equation
real, dimension(size(ptrad)) :: ZWORKI ! Term for the resolution of the equation
! optional for verification 4. 
!!$real, dimension(size(ptrad)) :: PTRAD2 ! body MRT for verification
!!$real, dimension(size(ptrad)) :: ZDIFF ! body MRT for verification
!!$integer :: jj
! 
!*       1.    set limits for the  wind value
!               ---------------------------
!  note: some tests suggest a lower limit of 0.001 
ZUA(:) =  min(MAX(0.2,PUA(:)) ,15.0) 

!*       2.    convection coefficient
!               ---------------------------
! convection coeff set to =  1.100 * 1.0E8 * U**0.6    ISO-ASHRAE - Kuehn et al. 1970
 ZWORKA(:) =  1.1* 1.0E8 * ZUA(:) **0.6   / (ZGE  *  ZGD**0.4)

!*       3.    analytic resolution
!               ---------------------------
 ZWORKB(:) = PTRAD(:) **4. + ZWORKA(:) * PTA(:) 

!!ZWORKm(:)=9.*ZWORKA(:)**2.
!!ZWORKN(:)=27.*ZWORKA(:)**4.
!!ZWORKp(:)=256.d0*ZWORKB(:)**3.
!!ZWORKq(:)= 3.4943 * ZWORKB(:)         

 ZWORKE(:)= ( (9.d0*ZWORKA(:)**2.) +  1.73205d0* &
             ( (27.d0*ZWORKA(:)**4.)+(256.d0*ZWORKB(:)**3.) )**0.5 )**(1./3.)
!! ZWORKE(:)= (ZWORKM(:)+1.73205*(ZWORKN(:)+ZWORKP(:))**0.5)**(1./3.)

 ZWORKK(:)=ZWORKE(:)*0.381571d0 - ((3.4943d0 * ZWORKB(:))/ZWORKE(:))
!!ZWORKK(:)=ZWORKE(:)*0.381571 -ZWORKQ(:)/ZWORKE(:)

ZWORKI(:)= SNGL(0.5d0 *  ( 2.0d0 * ZWORKA(:) /  ZWORKK(:)**0.5 - ZWORKK(:))**0.5)
ZWORKJ(:)= SNGL(0.5d0 * ZWORKK(:)**0.5)

PTGLOBE_BODY(:) = -1.0 * ZWORKJ(:) + ZWORKI(:)

!*       4.    optional : verification (reciprocity)
!               ---------------------------
!
!PTRAD2(:) = ( PTGLOBE_BODY(:)**4.0 +   &
!             ZWORKa(:)*(PTGLOBE_BODY(:) -PTA(:)) )**0.25

!ZDIFF(:)= abs(PTRAD2(:) -PTRAD(:))

!      DO JJ=1,SIZE(PTRAD)
! IF (ZDIFF(JJ) .gt. 0.1) THEN
!  print*,'** Tglobe diff is  /',ZDIFF(JJ) ,'i tr t u ',JJ, PTRAD(JJ), PTA(JJ),ZUA(JJ)
! ENDIF
!     ENDDO 

  end

