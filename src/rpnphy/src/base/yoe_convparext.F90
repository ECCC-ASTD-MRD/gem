!     ######################
      MODULE YOE_CONVPAREXT
!     ######################

IMPLICIT NONE

integer, SAVE :: JCVEXB ! start vertical computations at
                        ! 1 + JCVEXB = 1 + ( KBDIA - 1 )
integer, SAVE :: JCVEXT ! limit vertical computations to
                        ! KLEV - JCVEXT = KLEV - ( KTDIA - 1 )

END MODULE YOE_CONVPAREXT

