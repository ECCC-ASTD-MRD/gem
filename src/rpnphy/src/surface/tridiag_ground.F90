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
!      #########################################
       SUBROUTINE TRIDIAG_GROUND(PA,PB,PC,PY,PX)
!      #########################################
!
!
!!****   *TRIDIAG_GROUND* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to resolve the linear system:
!
!       A.X = Y
!
!      where A is a tridiagonal matrix, and X and Y two vertical vectors.
!     However, the computations are performed at the same time for all
!     the verticals where an inversion of the system is necessary.
!     This explain the dimansion of the input variables.
!
!!**   METHOD
!!     ------
!!                      
!!        Then, the classical tridiagonal algorithm is used to invert the 
!!     implicit operator. Its matrix is given by:
!!
!!     (  b(1)      c(1)      0        0        0         0        0        0  )
!!     (  a(2)      b(2)     c(2)      0  ...    0        0        0        0  ) 
!!     (   0        a(3)     b(3)     c(3)       0        0        0        0  ) 
!!      .......................................................................
!!     (   0   ...   0      a(k)      b(k)     c(k)       0   ...  0        0  ) 
!!      .......................................................................
!!     (   0         0        0        0        0 ...  a(n-1)   b(n-1)   c(n-1))
!!     (   0         0        0        0        0 ...     0      a(n)     b(n) )
!!
!!
!!       All these computations are purely vertical and vectorizations are 
!!     easely achieved by processing all the verticals in parallel.
!!
!!     EXTERNAL
!!     --------
!!
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!
!!     REFERENCE
!!     ---------
!!
!!     AUTHOR
!!     ------
!!       V. Masson
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original        May 13, 1998
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
!
implicit none
!!!#include <arch_specific.hf>
!
!
!*       0.1 declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA  ! lower diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB  ! main  diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PC  ! upper diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PY  ! r.h.s. term   
!
REAL,    DIMENSION(:,:), INTENT(OUT) :: PX  ! solution of A.X = Y 
!
!*       0.2 declarations of local variables
!
INTEGER           :: JK             ! vertical loop control
INTEGER           :: IN             ! number of vertical levels
!
REAL, DIMENSION(SIZE(PA,1)           ) :: ZDET ! work array
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: ZW   ! work array
! ---------------------------------------------------------------------------
!
IN=SIZE(PX,2)
!
!*       1.  levels going up
!            ---------------
!
!*       1.1 first level
!            -----------
!
ZDET(:)   = PB(:,1)

PX  (:,1) = PY(:,1) / ZDET(:)
!
!*       1.2 other levels
!            ------------
!
DO JK=2,IN
  ZW(:,JK)    = PC(:,JK-1)/ZDET(:)
  ZDET(:)       = PB(:,JK  ) - PA(:,JK)*ZW(:,JK)

  PX  (:,JK)    = ( PY(:,JK) - PA(:,JK)*PX(:,JK-1) ) / ZDET(:)
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.  levels going down
!            -----------------
!
DO JK=IN-1,1,-1
  PX  (:,JK) = PX(:,JK) - ZW(:,JK+1)*PX(:,JK+1)
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_GROUND
