!!!s/r Interp1D_NearestNeighbour - interpolation:  nearest-neighbour method
!subroutine Interp1D_NearestNeighbour_X  &
!subroutine Interp1D_NearestNeighbour_X8 &
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
  ! stateDerivDestn, extrapEnableDown, extrapGuideUp, extrapGuideDown,
  ! extrapGuideUp, flux, numExtArraysIn, numExtArraysOut, ExtArraysIn, 
  ! ExtArraysOut
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - new code
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output
! v1_3    Blezius J.W. OCT 2003 - extend the interface, but ignore the extra data
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Given the state at a set of srcNumLevels unevenly spaced levels,
!        vLevelSource, this routine interpolates the state at a different set of
!        destNumLevels unevenly spaced levels, vLevelDestn.  The algorithm
!        consists merely of assigning the state value of the source level that
!        is closest to the sought level.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"
!
!NOTES
!        Although this routine does not use the eight dummy parameters, they are
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
!        The use of ABS() is necessary in order to accept ascending or descending
!        source levels.
!
!!

  integer i                            ! horizontal loop index
  integer vt                           ! vertical loop index into the target grid
  real(real48) :: stateSourceBelow, stateSourceAbove, &
                  levelBelow, levelAbove

  ! Use the interpolation algorithm to estimate values for all points, even if
  ! it is actually EXTRApolation.  Doing so is faster than avoiding the
  ! duplication of calculation.  Such extrapolated estimates can be overwritten
  ! by calling an extrapolation routine.
  !
!$OMP parallel do private(i,levelBelow,levelAbove,stateSourceAbove,stateSourceBelow)
  do vt=1,destNumLevels                 ! for each target vertical point
    do i=1,numInterpSets                ! for each horizontal point

      ! Identify the 2 known points that straddle the sought point
      levelBelow = vLevelSource(i, posnDestInSrc(i,vt) )
      levelAbove = vLevelSource(i, posnDestInSrc(i,vt)+1 )

      ! Obtain the 2 corresponding known state values.
      stateSourceAbove = stateSource( i, posnDestInSrc(i,vt)+1 )
      stateSourceBelow = stateSource( i, posnDestInSrc(i,vt) )

      ! Select the closest source level.
      if(  abs(vLevelDestn(i,vt) - levelBelow) &
         < abs(vLevelDestn(i,vt) - levelAbove) ) then
        stateDestn(i,vt) = stateSourceBelow
      else
        stateDestn(i,vt) = stateSourceAbove
      endif

      ! In the case of this particular interpolation algorithm, clamping
      ! extrapolated values to the nearest source value is identical with
      ! extrapolation.  Therefore, no special attention is paid to the choice
      ! of the client.

    end do ! i
  end do ! vt
!$OMP END parallel do

#undef levelBelow
#undef levelAbove
end subroutine ! Interp1D_NearestNeighbour_X
