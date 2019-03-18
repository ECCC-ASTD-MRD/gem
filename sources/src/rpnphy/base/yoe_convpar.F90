!-------------------------------------------------------------------------------

!     ###################
      MODULE YOE_CONVPAR
!     ###################

!!****  *YOE_CONVPAR* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!      The purpose of this declarative module is to declare  the
!      constants in the deep convection parameterization.

!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (YOE_CONVPAR)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------

!*       0.   DECLARATIONS
!             ------------

IMPLICIT NONE

real, SAVE :: XA25        ! 25 km x 25 km reference grid area

real, SAVE :: XCRAD       ! cloud radius
real, SAVE :: XCDEPTH     ! minimum necessary cloud depth
real, SAVE :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)

real, SAVE :: XZLCL       ! maximum allowed allowed height
                          ! difference between departure level and surface
real, SAVE :: XZPBL       ! minimum mixed layer depth to sustain convection
real, SAVE :: XWTRIG      ! constant in vertical velocity trigger
real, SAVE :: XDTHPBL     ! temperature perturbation in PBL for trigger
real, SAVE :: XDRVPBL     ! moisture perturbation in PBL for trigger


real, SAVE :: XNHGAM      ! accounts for non-hydrost. pressure
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
real, SAVE :: XTFRZ1      ! begin of freezing interval
real, SAVE :: XTFRZ2      ! end of freezing interval

real, SAVE :: XRHDBC      ! relative humidity below cloud in downdraft

real, SAVE :: XRCONV      ! constant in precipitation conversion
real, SAVE :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
real, SAVE :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
real, SAVE :: XUSRDPTH    ! pressure thickness used to compute updraft
                          ! moisture supply rate for downdraft
real, SAVE :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                          ! allowed below  melting level
real, SAVE :: XUVDP       ! constant for pressure perturb in momentum transport

END MODULE YOE_CONVPAR

!-------------------------------------------------------------------------------

