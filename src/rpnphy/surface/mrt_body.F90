
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
MODULE MODI_MRT_BODY

!#TODO: no need for an interface here, put the funtion MRT_BODY in the module to avoid non consistency between module and function interface

INTERFACE
FUNCTION MRT_BODY(PEMISS_BODY, PQ1, PQ2, PQ3, PQ4, PQ5, PQ6,PQ7, PQ8   &
                  ,PQ9, PQ10, PQ11, PQ12) RESULT(PTRAD_BODY)
REAL, DIMENSION(:), INTENT(IN) :: PQ1
REAL, DIMENSION(:), INTENT(IN) :: PQ2
REAL, DIMENSION(:), INTENT(IN) :: PQ3
REAL, DIMENSION(:), INTENT(IN) :: PQ4
REAL, DIMENSION(:), INTENT(IN) :: PQ5
REAL, DIMENSION(:), INTENT(IN) :: PQ6
REAL, DIMENSION(:), INTENT(IN) :: PQ7
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ8
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ9
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ10
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ11
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ12
REAL, DIMENSION(SIZE(PQ1))    :: PTRAD_BODY
REAL, DIMENSION(SIZE(PQ1))    :: ZQSHADE
REAL :: PEMISS_BODY 
END FUNCTION MRT_BODY
END INTERFACE
END MODULE MODI_MRT_BODY


!   ##########################################################################
FUNCTION MRT_BODY(PEMISS_BODY, PQ1, PQ2, PQ3, PQ4, PQ5, PQ6,PQ7, PQ8   &
                  ,PQ9, PQ10,PQ11, PQ12) RESULT(PTRAD_BODY)
!   ##########################################################################
!
!    PURPOSE       : computes the mean radiant temperature from the energy budget received by a body 
!    AUTHOR        :  S. Leroyer   (Original  07/2018)
!    REFERENCE     :  Leroyer et al. (2018)
!    MODIFICATIONS :  
!    METHOD        :  
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS, ONLY : XSTEFAN
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
REAL, DIMENSION(:), INTENT(IN) :: PQ1
REAL, DIMENSION(:), INTENT(IN) :: PQ2
REAL, DIMENSION(:), INTENT(IN) :: PQ3
REAL, DIMENSION(:), INTENT(IN) :: PQ4
REAL, DIMENSION(:), INTENT(IN) :: PQ5
REAL, DIMENSION(:), INTENT(IN) :: PQ6
REAL, DIMENSION(:), INTENT(IN) :: PQ7
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ8
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ9
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ10
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ11
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PQ12
REAL, DIMENSION(SIZE(PQ1))    :: PTRAD_BODY
REAL, DIMENSION(SIZE(PQ1))    :: ZQSHADE
!
REAL :: PEMISS_BODY
INTEGER :: JJ
!
DO JJ = 1, SIZE(PQ1) 

 !*  1 - calculation of MRT street

   IF (PRESENT(PQ8) .AND. PRESENT(PQ9)  & 
   .AND. PRESENT(PQ10) .AND. PRESENT(PQ11) .AND. PRESENT(PQ12) ) THEN
   ZQSHADE(JJ) = PQ8(JJ)+PQ9(JJ)+PQ10(JJ)+PQ11(JJ)+PQ12(JJ) 
 else
   ZQSHADE(JJ) = PQ2(JJ)+PQ3(JJ)+PQ4(JJ)+PQ5(JJ)+PQ6(JJ)+PQ7(JJ)
 ENDIF

   PTRAD_BODY(JJ) = ((PQ1(JJ)+ZQSHADE(JJ))/PEMISS_BODY/XSTEFAN)**0.25

ENDDO
!
END FUNCTION MRT_BODY
