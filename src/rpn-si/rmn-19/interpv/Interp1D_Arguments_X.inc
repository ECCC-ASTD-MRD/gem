!!!file Interp1D_Arguments_X - The arguments are common to all interpolaton and
!                              extrapolation functions in the 1-D package, except
!                              for Interp1D_FindPos
!
!AUTHOR
!     J.W. Blezius OCT 2003 extend the interface
!
!REVISION
! v1_3    Blezius J.W. OCT 2003 - initial code
!
!OBJECT
!        To minimize errors, the declaration of the arguments are written down
!        only once:  here.  It is by design that the arguments to all functions
!        are identical so that the function that is executed can be changed
!        simply by changing the name, or a pointer to the function if it is
!        called from C.
!
!ARGUMENTS
#include "Interp1D_Arguments.inc"

  external flux                         ! s/r to calculate flux at one vert level
  integer, intent(in) :: numExtArraysIn ! size of the array, ExtArraysIn
  integer, intent(in) :: numExtArraysOut! size of the array, ExtArraysOut
  real(real48), dimension(src_ijDim, numExtArraysIn), intent(in)  :: ExtArraysIn
  real(real48), dimension(dst_ijDim, numExtArraysOut),intent(out) :: ExtArraysOut

!WHERE:
                                        ! These are indices into the ExtArraysIn
  integer, parameter :: iZ0=1, &
                        iILMO=2, &
                        iH_BOUND=3, &
                        iFT=4, &
                        iANGLE_TOP=5, &
                        iLATITUDE=6
!
!  if 2nd index of ExtArraysIn
!     [The first three values are simply passed on to the call-back routines]
!      = iZ0:              roughness length
!        iILMO:            inverse of Monin-Obukhov length
!        iH_BOUND:         height of the boundary layer
!        iFT:              vertical flux at the surface
!                          [in (units of the state) / (units returned by flux())]
!                          In the case of wind:  normalized surface stress
!        iANGLE_TOP:       ft direction at surface layer top [used only for wind]
!        iLATITUDE:        [used only for wind]

                                        ! These are indices into the ExtArraysOut
  integer, parameter :: iY_DESTINATION = 1
!
!  if 2nd index of ExtArraysOut:
!      = 1  to destnumlevels:  y-component of stateDestn [used only for wind]
!      = to 2* destnumlevels:  y-component of stateDerivDestn[used only for wind]

                                        ! These values are used to check the
                                        ! supplied array sizes
#define MAX_EXT_ARRAYS_IN  (6)
#define MAX_EXT_ARRAYS_OUT (2*destnumlevels)
!
!NOTES
!        This file is to be included by the line,
!                    '#include "Interp1D_Arguments_X.inc" '.
!
!        It is assumed that wherever this file is included, the inclusion is
!        preceded by 'use Interp1D_Constants'.
!
!        This file defines the parameters that are unique to the 'extended'
!        interface.  The original interface has its parameters defined in
!        Interp1D_Arguments.inc.
!
!!

