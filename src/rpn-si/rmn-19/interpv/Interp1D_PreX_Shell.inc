!!!s/r Shell - a general shell to support the pre-extension interface
!subroutine RoutineWithOldInterface  &
                      (numInterpSets, srcNumLevels, destNumLevels, &
                       src_ijDim, dst_ijDim, &

                       vLevelSource, stateSource, stateDerivSource, &

                       posnDestInSrc, vLevelDestn, stateDestn, stateDerivDestn, &

                       extrapEnableDown, extrapEnableUp, &
                       extrapGuideDown, extrapGuideUp &
                      )
!
!AUTHOR
!     J.W. Blezius OCT 2003 support the old interface
!
!REVISION
! v1_3    Blezius J.W. OCT 2003 - new code
!
!OBJECT
!        Revision 1.3 has an extended interface.  The object of this module is to
!        transparently support the interface that existed previous to that.  This
!        file is included in other files so as to create a shell around a routine
!        with the extended interface.  The including file declares the name of
!        the old interface, and passes the equivalent name of the new interface
!        to this file as a defined name, DestnFunc.  Thus, this file is included
!        as many times as there are required shell routines.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments.inc"
!
!NOTES
!        See the called routine for more details.
!
!!

  external DestnFunc                    ! any routine would do

                                        ! supply null arrays, to satisfy compiler
  integer :: numExtArraysIn
  integer :: numExtArraysOut
  real(real48), dimension(1, 1)  :: ExtArraysIn
  real(real48), dimension(1, 1) :: ExtArraysOut

  numExtArraysIn = 0
  numExtArraysOut = 0


  call DestnFunc(numInterpSets, srcNumLevels, destNumLevels, &
                 src_ijDim, dst_ijDim, &

                 vLevelSource, stateSource, stateDerivSource, &

                 posnDestInSrc, vLevelDestn, stateDestn, stateDerivDestn, &

                 extrapEnableDown, extrapEnableUp, &
                 extrapGuideDown, extrapGuideUp, &

                 DestnFunc, numExtArraysIn, numExtArraysOut, &
                 ExtArraysIn, ExtArraysOut &
                )

end subroutine
