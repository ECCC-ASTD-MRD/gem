!!!s/r Interp1D_CubicWithDerivs - interpolation: cubic with derivatives supplied
!subroutine Interp1D_CubicWithDerivs_X  &
!subroutine Interp1D_CubicWithDerivs_X8 &
                      (numInterpSets, srcNumLevels, destNumLevels, &
                       src_ijDim, dst_ijDim, &

                       vLevelSource, stateSource, stateDerivSource, &

                       posnDestInSrc, vLevelDestn, stateDestn, stateDerivDestn, &

                       extrapEnableDown, extrapEnableUp, &
                       extrapGuideDown, extrapGuideUp, &

                       flux, numExtArraysIn, numExtArraysOut, &
                       ExtArraysIn, ExtArraysOut &
                      )
  ! These are dummy arguments for this function:  extrapGuideDown, extrapGuideUp,
  ! flux, numExtArraysIn, numExtArraysOut, ExtArraysIn,  ExtArraysOut
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - this routine was copied from vterp2, and
!                                 the comments were added
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output! v1_3    Blezius J.W. OCT 2003 - extend the interface, but ignore the extra data
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Given the state, S, and its derivative, S', at a set of
!        srcNumLevels unevenly spaced levels, vLevelSource, this
!        routine interpolates S and S' at a different set of destNumLevels
!        unevenly spaced levels, vLevelDestn, based on estimates of the
!        derivatives of S up to the third degree; hence the description, cubic.
!
!        Each target level is considered separately to estimate a state value at
!        that location.  The known state values (and derivatives) at levels just
!        above and just below the target level are used to calculate various
!        derivatives at the centre level (midway between the two known levels)
!        and at the target level.  The final result is a 'best' estimate of S and
!        S' at each target level.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"

!NOTES
!        Each routine of the set of interpolation routines takes exactly the same
!        parameters, allowing the user to use the routines interchangeably, as
!        with a pointer.
!
!        The values that are calculated at the centre level are effectively
!        averages.  However, any reference to average is omitted from the
!        notation for the sake of clarity.
!
!        In a general context, the (vertical) levels could be interpreted to be
!        a domain with any physical meaning, and the state (here thought of as
!        any state variable of the atmosphere) could be interpreted to be a range
!        with any physical meaning.
!
!        For now, extrapolation is not performed.  It is presumed that it is not
!        requested.
!
!!


  ! VARIABLES:
  !        The known state values (and derivatives) at levels just
  !        above and just below the target level are used to calculate various
  !        values at the centre level (midway between the two known levels) and
  !        at the target level.  The variable names that represent these values
  !        all begin with s, for 'state'.  The second character is either m for
  !        'mean' or t for 'target'.  The third is a numeral that indicates the
  !        degree of the derivative that that variable represents.  Each
  !        variable in this group represents an estimate of the value
  !        described:
  real(double) :: sc0 ! state, at the centre, 0th order derivative
  real(double) :: sc1 ! state, at the centre, 1st order derivative - half of it
  real(double) :: sc2 ! state, at the centre, 2nd order derivative - half of it
  real(double) :: sc3 ! state, at the centre, 3rd order derivative - 1/4 of it

  real(double) :: st0 ! state, at the target, 0th order derivative
  real(double) :: st1 ! state, at the target, 1st order derivative
  real(double) :: st2 ! state, at the target, 2nd order derivative



  ! OTHER VARIABLES:

  ! In what follows, BELOW refers to the level index which is smaller, although
  ! the level that is indexed may actually be greater.  This does not change 
  ! the mathematics.
  !
  ! posnDestInSrc(h,vt)   !index to the input level just BELOW the target level
  ! posnDestInSrc(h,vt)+1 !index to the input level just ABOVE the target level

  ! These three levels are in pressure units (mb)
  real(double) :: levelBelow  ! input vertical level just BELOW target level
  real(double) :: levelAbove  ! input vertical level just ABOVE target level
  real(double) :: levelCentre ! level midway between levelAbove & levelBelow

  real(double) :: sb0    ! state at levelBelow (units known only to the caller)
  real(double) :: sb1    ! state derivative at levelBelow
  real(double) :: sa0    ! state at levelAbove (units known only to the caller)
  real(double) :: sa1    ! state derivative at levelAbove

  ! the difference in level from the target to the input-level centre
  real(double) :: deltaLTar

  ! This is twice the difference of input level from the bounds to the centre,
  ! for CPU convenience.  Half of this value is analogous with deltaLTar.
  real(double) :: deltaLIn2

  real(double) :: over2DeltaLIn    ! 1 / deltaLIn2
  real(double) :: diffDeltaLSquares ! = ( (deltaLIn2/2)**2 - deltaLTar**2 ) / 2

  integer i              ! horizontal loop index
  integer vt             ! vertical loop index into the target grid

  integer smallestLevel                 ! index of smallest source vertical level
  integer largestLevel                  ! index of largest source vertical level


  ! Verify that enough source levels are supplied; the algorithm requires 2
  if(srcNumlevels < 2) then
    write(6,*)'ERROR:  the cubic algorithm with derivatives requires 2 source points.'
    return
  end if

  ! Determine whether the source levels are increasing or decreasing
  if( vLevelSource(1,1) < vLevelSource(1,2) ) then
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
!$OMP parallel do private(i,levelBelow,levelAbove,levelCentre,sb0,sa0,sb1,sa1,deltaLIn2,over2DeltaLIn,deltaLTar,diffDeltaLSquares,sc0,sc1,sc2,sc3,st2,st1,st0)
  do vt=1,destNumLevels                 ! for each target vertical point
    do i=1,numInterpSets                ! for each horizontal point

      ! identify the levels that bound the interpolation domain
      levelBelow = vLevelSource(i, posnDestInSrc(i,vt) )
      levelAbove = vLevelSource(i, posnDestInSrc(i,vt)+1 )
      levelCentre = 0.5 * (levelBelow + levelAbove)

      ! identify the state and its derivative at those levels
      sb0 = stateSource( i, posnDestInSrc(i,vt) )
      sa0 = stateSource( i, posnDestInSrc(i,vt)+1 )
      sb1 = stateDerivSource( i, posnDestInSrc(i,vt) )
      sa1 = stateDerivSource( i, posnDestInSrc(i,vt)+1 )

      ! calculate the differences from the centre point
      deltaLIn2 = levelAbove - levelBelow
      over2DeltaLIn = 1 / deltaLIn2
      deltaLTar = vLevelDestn(i,vt) - levelCentre

      ! The deltaLTar term of diffDeltaLSquares is used to complete the Taylor
      ! expansion; the deltaLIn2 is used to add a correction for the estimated
      ! state at the centre point that results from the knowledge of the
      ! derivatives at the level bounds.
      !
      diffDeltaLSquares = 0.125*deltaLIn2*deltaLIn2 - 0.5*deltaLTar*deltaLTar

      ! estimate the derivatives at the centre level
      sc0 = 0.5 * (sa0 + sb0)
      sc1 = (sa0 - sb0) * over2DeltaLIn
      sc2 = (sa1 - sb1) * over2DeltaLIn
      sc3 = ( (sa1 - sc1) - (sc1 - sb1) ) * over2DeltaLIn * over2DeltaLIn

      ! estimate the derivatives at the target level
      st2 = sc2 + 2. * deltaLTar * sc3
      st1 = sc1 + deltaLTar * st2
      st0 = sc0 + deltaLTar * sc1

      ! make the final estimates
      stateDestn(i, vt) = st0 - diffDeltaLSquares * st2
      stateDerivDestn(i, vt) = st1 - 2. * diffDeltaLSquares * sc3

      ! clamp any extrapolated values to the nearest source value.
      if( .not. extrapEnableDown .and. &
                         vLevelDestn(i,vt) < vLevelSource(i,smallestLevel) ) then
        stateDestn(i,vt) = stateSource(i,smallestLevel)
        stateDerivDestn(i,vt) = stateDerivSource(i,smallestLevel)
      else if( .not. extrapEnableUp .and. &
                         vLevelDestn(i,vt) > vLevelSource(i,largestLevel) ) then
        stateDestn(i,vt) = stateSource(i,largestLevel)
        stateDerivDestn(i,vt) = stateDerivSource(i,largestLevel)
      end if

    end do ! i
  end do ! vt
!$OMP END parallel do

end subroutine ! Interp1D_CubicWithDerivs_X
