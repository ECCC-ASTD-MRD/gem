!     ###########################
      SUBROUTINE SU_CONVPAR_SHAL
!     ###########################

!!****  *SU_CONVPAR * - routine to initialize the constants modules
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to initialize  the constants
!!     stored in  modules YOE_CONVPAR_SHAL
!!
!!
!!**  METHOD
!!    ------
!!      The shallow convection constants are set to their numerical values
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR_SHAL   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module YOE_CONVPAR_SHAL, routine SU_CONVPAR)
!!
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
use cnv_options
USE YOE_CONVPAR_SHAL

IMPLICIT NONE
!!!#include <arch_specific.hf>

!-------------------------------------------------------------------------------

!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------

XCDEPTH     = 0.5E3  ! minimum necessary shallow cloud depth
XCDEPTH_D   = 3.0E3  ! maximum allowed shallow cloud depth (m)
XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)

XZLCL    = 0.5E3     ! maximum allowed height (m)
                          ! difference between the DPL and the surface
XZPBL    = 40.E2     ! minimum mixed layer depth to sustain convection
XWSTARMIN= 0.1       ! minimum convective velocity scale value (m/s)


XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
XTFRZ1   = 273.16    ! begin of freezing interval (K)
XTFRZ2   = 250.16    ! end of freezing interval (K)


XSTABT   = 0.95      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE

!PV - max for shallow = min for deep
XCDEPTH_D   = kfcdepth  ! maximum allowed shallow cloud depth (m)

END SUBROUTINE SU_CONVPAR_SHAL

