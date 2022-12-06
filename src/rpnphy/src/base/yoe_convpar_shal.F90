!     ########################
      MODULE YOE_CONVPAR_SHAL
!     ########################

!!****  *YOE_CONVPAR_SHAL* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare  the
!!      constants in the deep convection parameterization.
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (YOE_CONVPAR_SHAL)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  04/10/98
!-------------------------------------------------------------------------------

!*       0.   DECLARATIONS
!             ------------

IMPLICIT NONE

real, SAVE :: XCDEPTH     ! minimum necessary cloud depth
real, SAVE :: XCDEPTH_D   ! maximum allowed cloud thickness
real, SAVE :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)

real, SAVE :: XZLCL       ! maximum allowed allowed height
                          ! difference between departure level and surface
real, SAVE :: XZPBL       ! minimum mixed layer depth to sustain convection
real, SAVE :: XWSTARMIN   ! minimum convective velocity scale value (m/s)


real, SAVE :: XNHGAM      ! accounts for non-hydrost. pressure
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
real, SAVE :: XTFRZ1      ! begin of freezing interval
real, SAVE :: XTFRZ2      ! end of freezing interval


real, SAVE :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
real, SAVE :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE

END MODULE YOE_CONVPAR_SHAL

