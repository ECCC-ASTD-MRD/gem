!     #######################
      SUBROUTINE SU_CONVPAR1
!     #######################

!!****  *SU_CONVPAR * - routine to initialize the constants modules
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules YOE_CONVPAR


!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------


!*       0.    DECLARATIONS
!              ------------

USE YOE_CONVPAR

IMPLICIT NONE
!!!#include <arch_specific.hf>

!-------------------------------------------------------------------------------

!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------


XA25     = 625.E6    ! 25 km x 25 km reference grid area

XCRAD    =  500.     ! cloud radius (m)
!XCRAD    = 1500.     ! cloud radius
XCDEPTH  = 3.E3      ! minimum necessary cloud depth
XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)

XZLCL    = 3.5E3     ! maximum allowed allowed height
                          ! difference between the surface and the DPL (m)
XZPBL    = 60.E2     ! minimum mixed layer depth to sustain convection
XWTRIG   = 6.00      ! constant in vertical velocity trigger
XDTHPBL  = .1        ! Temp. perturbation in PBL for trigger (K)
XDRVPBL  = 1.e-4     ! moisture  perturbation in PBL for trigger (kg/kg)


XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
XTFRZ1   = 273.16    ! begin of freezing interval (K)
XTFRZ2   = 250.16    ! end of freezing interval (K)

XRHDBC   = 0.9       ! relative humidity below cloud in downdraft

XRCONV   = 0.015     ! constant in precipitation conversion
XRCONV   = 0.01     ! constant in precipitation conversion
XSTABT   = 0.75      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
XUSRDPTH = 165.E2    ! pressure thickness used to compute updraft
                          ! moisture supply rate for downdraft
XMELDPTH = 200.E2    ! layer (Pa) through which precipitation melt is
                          ! allowed below downdraft
XUVDP    = 0.7       ! constant for pressure perturb in momentum transport
!
!
END SUBROUTINE SU_CONVPAR1

