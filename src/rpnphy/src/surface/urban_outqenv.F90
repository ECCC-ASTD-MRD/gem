
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
!   ##########################################################################
SUBROUTINE URBAN_OUTQENV(PSCA_SW, PREF_SW_FAC, PREF_SW_GRND, PREF_SW_ROOF,  &
                       PEMIT_LW_FAC, PEMIT_LW_GRND, PEMIT_LW_ROOF, PLW_RAD, &
                       PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,    &
                       PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,                         &
                       PQ8,PQ9,PQ10,PQ11,PQ12,PQ13,OPT,OPT_BODY,            &
                       ZEB,ZAB,ZHB                                         )
!
!    PURPOSE       : Computes the energy budget received by a body over a slab surface in a street and over rooftop 
!    AUTHOR        :  S. Leroyer   (Original  10/2016)
!    REFERENCE     :  Leroyer et al. (2018)
!    MODIFICATIONS :  
!    METHOD        :  
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_CSTS, ONLY : XPI
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
INTEGER, INTENT(IN)  :: OPT
INTEGER, INTENT(IN)  :: OPT_BODY
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_GRND !Solar radiation reflected by ground [road + garden] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_FAC !Solar radiation reflected by facade [wall + glazing] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PSCA_SW !Diffuse solar radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PDIR_SW !Direct solar radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH !solar zenithal angle (rad from vert.)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_FAC !Longwave radiation emitted by the facade [wall + glazing] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_GRND !Longwave radiation emitted by the ground [road + garden] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD !Atmospheric longwave radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_ROOF !Solar radiation reflected by the roof (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF !Longwave radiation emitted by the roof (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PBLD !Building surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PBLD_HEIGHT !Building surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PWALL_O_HOR !Building surface fraction
REAL, DIMENSION(:), INTENT(OUT) :: PQ1
REAL, DIMENSION(:), INTENT(OUT) :: PQ2
REAL, DIMENSION(:), INTENT(OUT) :: PQ3
REAL, DIMENSION(:), INTENT(OUT) :: PQ4
REAL, DIMENSION(:), INTENT(OUT) :: PQ5
REAL, DIMENSION(:), INTENT(OUT) :: PQ6
REAL, DIMENSION(:), INTENT(OUT) :: PQ7
REAL, DIMENSION(:), INTENT(OUT) :: PQ8
REAL, DIMENSION(:), INTENT(OUT) :: PQ9
REAL, DIMENSION(:), INTENT(OUT) :: PQ10
REAL, DIMENSION(:), INTENT(OUT) :: PQ11
REAL, DIMENSION(:), INTENT(OUT) :: PQ12
REAL, DIMENSION(:), INTENT(OUT) :: PQ13
!
REAL :: ZHB  ! average height of the body
REAL :: ZAB  !absorption coef of solar radiation by the  body       
REAL :: ZEB  !emissivity of  body                              
!   ##########################################################################

!*      0.2    declarations of local variables

REAL, DIMENSION(SIZE(PDIR_SW)) :: ZWROAD !width of the road (m)
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZWROOF !width of the roof (m)                                !!! INPUT ?????
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZL1, ZL2, ZL4 !lengths for view factor calculation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZHD, ZL3, zz, yy, kk !lengths for view factor calculation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZFFAC !facade view factor of human body
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZFGRND !ground view factor of human body
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZFROOF! roof view factor of human body
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZFSKY !sky view factor of human body
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIRSWBODY !solar radiation received by human body
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZELEV !solar elevation angle
  !
!*  1 - calculation of view factors
ZFFAC(:) =  0.
ZFGRND(:) = 0.
ZFROOF(:) = 0.
ZFSKY(:) = 0.
!
!=====================
 IF (OPT .eq. 1) THEN    ! flat surface
!======================
ZFGRND(:) = 0.5 
ZFSKY(:) = 1.0 - ZFGRND(:)
! sky
PQ2(:)=ZAB*PSCA_SW(:)*ZFSKY(:)
PQ3(:)=ZEB*PLW_RAD(:)*ZFSKY(:)
! ground
PQ4(:)=ZAB*PREF_SW_GRND (:)*ZFGRND(:)
PQ5(:)=ZEB*PEMIT_LW_GRND(:)*ZFGRND(:)

!======================
  ELSEIF (OPT .eq. 2) THEN  ! street
!=======================
  ZWROAD(:) = PBLD_HEIGHT(:) * 2. * (1. - PBLD(:)) / PWALL_O_HOR(:)
!
  ZL1(:) = SQRT(ZHB**2                  + (ZWROAD(:)/2.)**2)
  ZL2(:) = SQRT( PBLD_HEIGHT(:)**2      + (ZWROAD(:)/2.)**2)
  ZL4(:) = SQRT((PBLD_HEIGHT(:)-ZHB)**2 + (ZWROAD(:)/2.)**2)
!
! wall view factor
  ZFFAC (:) = (ZL1(:) + ZL2(:) - ZWROAD(:)/2. - ZL4(:)) / (2. * ZHB)
! ground view factor
  ZFGRND(:) = 0.5*ZWROAD(:)/ZHB
  ZFGRND(:) = 0.5 * (ZFGRND(:) + 1. - SQRT(ZFGRND(:)**2 + 1.) ) 
! sky view factor 
  ZFSKY (:) = 1. - ZFFAC(:) - ZFGRND(:)

! sky
PQ2(:)=ZAB*PSCA_SW(:)*ZFSKY(:)
PQ3(:)=ZEB*PLW_RAD(:)*ZFSKY(:)
! ground
PQ4(:)=ZAB*PREF_SW_GRND (:)*ZFGRND(:)
PQ5(:)=ZEB*PEMIT_LW_GRND(:)*ZFGRND(:)
! facade
PQ6(:)=ZAB*PREF_SW_FAC (:)*ZFFAC(:)
PQ7(:)=ZEB*PEMIT_LW_FAC(:)*ZFFAC(:)

!========================
  ELSEIF (OPT .eq. 3) THEN   ! rooftop
!========================
  ZWROAD(:) = PBLD_HEIGHT(:) * 2. * (1. - PBLD(:)) / PWALL_O_HOR(:)
  ZWROOF(:) = PBLD_HEIGHT(:) * 2. * PBLD(:) / PWALL_O_HOR(:)
!
! roof view factor 
  ZFROOF(:) = 0.5 * ZWROOF(:)/ZHB
  ZFROOF(:) = 0.5 * (ZFROOF(:) + 1. - SQRT(ZFROOF(:)**2 + 1.) )
!
! ground (road) view factor --> simplification assume 0 although could be >0 if ZHD>BLDH
  ZFGRND(:) = 0.0   
!
! height d of the facade seen by the body on the roof
  ZHD(:) = MIN(PBLD_HEIGHT(:), 2.0* ZHB* ZWROAD(:)/ZWROOF(:) )
!
! coefficients to compute facade view factor of one portion of the body zz
  yy(:)=0.5*ZHD(:)* ZWROOF(:)/ZWROAD(:)
  WHERE (yy(:)>ZHB) yy(:)=ZHB
  zz(:)=ZHB - yy(:)

  kk(:) = ZWROAD(:) + (ZWROOF(:)/2.)

  ZL1(:) = SQRT( yy(:)**2            + kk(:)**2 )
  ZL2(:) = SQRT( (ZHD(:)+ZHB)**2     + kk(:)**2 )  
  ZL3(:) = SQRT( ZHB**2              + kk(:)**2 )
  ZL4(:) = SQRT( (ZHD(:)+yy(:))**2   + kk(:)**2 )

! wall view factor (ZHD seen)
 WHERE (zz(:) < 0.1) ZFFAC(:)  = 0.0
 WHERE (zz(:) > 0.0999) ZFFAC(:)  = ( ZL1(:) + ZL2(:) - ZL3(:) - ZL4(:) ) / (2. * ZZ(:))
!                                                                                      
! sky view factor  (can be a slighty different from 0.5)
ZFSKY(:)  = 1.0 - ZFROOF(:) - ZFFAC(:) - ZFGRND(:)

! roof
PQ8(:)=ZAB*PREF_SW_ROOF(:)*ZFROOF(:) 
PQ9(:)=ZEB*PEMIT_LW_ROOF(:)*ZFROOF(:)
! sky
PQ10(:)=ZAB*PSCA_SW(:)*ZFSKY(:)
PQ11(:)=ZEB*PLW_RAD(:)*ZFSKY(:)
! wall
PQ12(:)=ZAB*PREF_SW_FAC (:)*ZFFAC(:) 
PQ13(:)=ZEB*PEMIT_LW_FAC(:)*ZFFAC(:)

 !========================
 ENDIF

!=====================================
! ADD CONTRIBUTION FROM THE SUN
! OPT_BODY = 1  human body, for UTCI
! OPT_BODY = 2  globse sensor, for WBGT
!======================================
 ZELEV(:) = ( XPI/2. - PZENITH(:)) *180./XPI
 WHERE (ZELEV(:) < 1E-6) ZELEV(:) = 0.

 IF (OPT_BODY .eq. 1) THEN 
! the direct solar radiation is weighted by a projected area factor which can be expressed by this equation
! for a rotationally symmetric human being (Fanger, 1970)               
 ZDIRSWBODY(:) = PDIR_SW(:) / MAX( COS(PZENITH(:)) ,0.1)              &
                 * 0.308 * COS( xPI/180. *ZELEV(:)* (1.-ZELEV(:)**2/48402.) )
 ELSEIF (OPT_BODY .eq. 2) THEN
! sphere Gaspar and Quintela 2009                      
      ZDIRSWBODY(:) = PDIR_SW(:) * 0.25 /  MAX( COS(PZENITH(:)) ,0.1) 
 ENDIF
!
  PQ1(:)=ZAB* ZDIRSWBODY(:) 
!
!
END SUBROUTINE URBAN_OUTQENV
