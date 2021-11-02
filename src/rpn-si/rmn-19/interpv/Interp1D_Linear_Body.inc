!!!s/r Interp1D_Linear - interpolation:  linear method
! subroutine Interp1D_Linear_X  &
! subroutine Interp1D_Linear_X8 &
                      (numInterpSets, srcNumLevels, destNumLevels, &
                       src_ijDim, dst_ijDim, &

                       vLevelSource, stateSource, stateDerivSource, &

                       posnDestInSrc, vLevelDestn, stateDestn, stateDerivDestn, &

                       extrapEnableDown, extrapEnableUp, &
                       extrapGuideDown, extrapGuideUp, &

                       flux, numExtArraysIn, numExtArraysOut, &
                       ExtArraysIn, ExtArraysOut &
                      )
  ! These are dummy arguments for this function:  stateDerivSource,
  ! stateDerivDestn, extrapGuideDown, extrapGuideUp, flux, numExtArraysIn,
  ! numExtArraysOut, ExtArraysIn,  ExtArraysOut
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - new code
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output! v1_3    Blezius J.W. OCT 2003 - extend the interface, but ignore the extra data
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Given the state at a set of srcNumLevels unevenly spaced levels,
!        vLevelSource, this routine interpolates the state at a different set of
!        destNumLevels unevenly spaced levels, vLevelDestn, based on linear
!        interpolation between two known values that straddle the sought value.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"
!
!NOTES
!        Although this routine does not use the six dummy parameters, they are
!        included in the parameter list for the convenience of the user.  Each
!        routine of the set of interpolation routines takes exactly the same
!        parameters, allowing the user to use the routines interchangeably, as
!        with a pointer.
!
!        In a general context, the (vertical) levels could be interpreted to be
!        a domain with any physical meaning, and the state (here thought of as
!        any state variable of the atmosphere) could be interpreted to be a range
!        with any physical meaning.
!
!!

  real(real48) :: lev1, s1, &          ! the two points used for an interpolation
                  lev2, s2
  real(real48) :: slope

  integer i                            ! horizontal loop index
  integer vt                           ! vertical loop index into the target grid

  integer smallestLevel                 ! index of smallest source vertical level
  integer largestLevel                  ! index of largest source vertical level

  ! Verify that enough source levels are supplied; the algorithm requires 2
  if(srcNumlevels < 2) then
    write(6,*)'ERROR:  the linear algorithm requires 2 source points.'
    return
  end if

  ! Determine whether the source levels are increasing or decreasing
  if( vlevelSource(1,1) < vLevelSource(1,2) ) then
    ! source vertical levels are in increasing order
    smallestLevel = 1
    largestLevel  = srcNumlevels
  else
    ! source vertical levels are in decreasing order
    smallestLevel = srcNumlevels
    largestLevel = 1
  endif

  ! Use the interpolation algorithm to estimate values for all points, even if
  ! it is actually EXTRApolation.  Doing so is faster than avoiding the
  ! duplication of calculation.  Such extrapolated estimates can be overwritten
  ! by calling an extrapolation routine.
  !
!$OMP parallel do private(i,s1,s2,lev1,lev2,slope)
  do vt=1,destNumLevels                 ! for each target vertical point
    do i=1,numInterpSets                ! for each horizontal point

      ! Identify the 2 known points that straddle the sought point
      lev1 = vLevelSource(i, posnDestInSrc(i,vt) )
      s1 = stateSource(i, posnDestInSrc(i,vt) )

      lev2 = vLevelSource(i, posnDestInSrc(i,vt)+1 )
      s2 = stateSource(i, posnDestInSrc(i,vt)+1 )

      ! 'Draw' a straight line through those two points to the sought level
      slope = (s2 - s1) / (lev2 - lev1)
      stateDestn(i,vt) = slope * (vLevelDestn(i,vt) - lev1) + s1

      ! clamp any extrapolated values to the nearest source value.
      if( .not. extrapEnableDown .and. &
                         vLevelDestn(i,vt) < vLevelSource(i,smallestLevel) ) then
        stateDestn(i,vt) = stateSource(i,smallestLevel)
      else if( .not. extrapEnableUp .and. &
                         vLevelDestn(i,vt) > vLevelSource(i,largestLevel) ) then
        stateDestn(i,vt) = stateSource(i,largestLevel)
      end if

    end do ! i
  end do ! vt
!$OMP END parallel do

end subroutine ! Interp1D_Linear_X
