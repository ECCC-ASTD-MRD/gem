!!!file Interp1D_Arguments - The arguments are common to all interpolaton and
!                            extrapolation functions in the 1-D package, except
!                            for Interp1D_FindPos
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - initial code
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output
!
!OBJECT
!        To minimize errors, the declaration of the arguments are written down
!        only once:  here.  It is by design that the arguments to all functions
!        are identical so that the function that is executed can be changed
!        simply by changing the name, or a pointer to the function if it is
!        called from C.
!
!ARGUMENTS
  ! elements used (containing data):
  integer, intent(in) :: numInterpSets  ! no. sets (horiz. pts) for interpolation
  integer, intent(in) :: srcNumLevels   ! no. source vLevels containing data
  integer, intent(in) :: destNumLevels  ! no. destination vLevels containing data

  ! elements dimensioned:
  integer, intent(in) :: src_ijDim      ! horizontal dimension of source arrays
  integer, intent(in) :: dst_ijDim      ! horizontal dimension of dest'n arrays

                                        ! source co-ord's (vertical levels)
  real(real48), dimension(src_ijDim, srcNumLevels), intent(in) :: vLevelSource

                                        ! source state values (e.g. temp) at the
                                        !   source co-ordinates
  real(real48), dimension(src_ijDim, srcNumLevels), intent(in) :: stateSource
      
                                        ! source derivative state values at the
                                        !   source co-ordinates
  real(real48), dimension(src_ijDim, srcNumLevels), intent(in) ::stateDerivSource


                                        ! position (index) of vLevelSource that
                                        !   is just before the vLevelDestn
  integer, dimension(dst_ijDim, destNumLevels), intent(in) :: posnDestInSrc

                                        ! destination co-ord's (vert. levels)
  real(real48), dimension(dst_ijDim, destNumLevels), intent(in) :: vLevelDestn

                                        ! destination state values (e.g. temp))
                                        !   at the source co-ordinates
  real(real48), dimension(dst_ijDim, destNumLevels), intent(out) :: stateDestn

                                        ! destination derivative state values at
                                        !   the source co-ordinates
  real(real48), dimension(dst_ijDim, destNumLevels), intent(out)::stateDerivDestn


                                        ! data controlling extrapolation
  logical, intent(in) :: extrapEnableDown, extrapEnableUp
  real(real48), intent(in) :: extrapGuideDown, extrapGuideUp
!
!NOTES
!        This file is to be included by the line,
!                    '#include "Interp1D_Arguments.inc" '.
!
!        It is assumed that wherever this file is included, the inclusion is
!        preceded by 'use Interp1D_Constants'.
!
!!
