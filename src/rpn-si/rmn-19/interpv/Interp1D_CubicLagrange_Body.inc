!!!s/r Interp1D_CubicLagrange - interpolation:  cubic Lagrange method
!subroutine Interp1D_CubicLagrange_X  &
!subroutine Interp1D_CubicLagrange_X8 &
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
! v1_0    Blezius J.W.          - this routine was copied from ez_irgdint_3_wnnc
!                                 and adapted for a single dimension
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output
! v1_3    Blezius J.W. OCT 2003 - extend the interface, but ignore the extra data
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Given the state at a set of srcNumLevels unevenly spaced levels,
!        vLevelSource, this routine interpolates the state at a different set of
!        destNumLevels unevenly spaced levels, vLevelDestn, based on a Lagrange
!        polynomial that is fitted to four known values that straddle the sought
!        value.
!
!        Each target level is considered separately to estimate a state level at
!        that location. Four known state values straddling the value in question,
!        two above it and two below it if possible, are used to define a Lagrange
!        polynomial.  The value of that polynomial at the desired location is
!        then reported.
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
!        The algorithm was simply copied from the routine, ez_irgdint_3_wnnc, of
!        the ezscint library package, with no particular attempt to illuminate
!        it.  The interested reader is referred to the literature for background
!        information.
!
!        In a general context, the (vertical) levels could be interpreted to be
!        a domain with any physical meaning, and the state (here thought of as
!        any state variable of the atmosphere) could be interpreted to be a range
!        with any physical meaning.
!
!!

  real(double) :: x1,x2,x3,x4 ! the 4 points used to interpolate current point
  real(double) :: z1,z2,z3,z4 ! the 4 known values of state at x1,x2,x3,x4
  real(double) :: xt          ! position of the value currently sought
  real(double) :: cx1, cx2, cx3, cx4, cx5, cx6 ! difference table of xn's
  real(double) :: a1,a2,a3,a4 ! Lagrange interpolation co-efficients

  integer i                           ! horizontal loop index
  integer vt                          ! vertical loop index into the target grid

  integer firstLevel                  ! the index of the first level of four
                                      ! that are used for the current
                                      ! interpolation

  integer smallestLevel                 ! index of smallest source vertical level
  integer largestLevel                  ! index of largest source vertical level

  ! Verify that enough source levels are supplied; the algorithm requires 4
  if(srcNumlevels < 4) then
    write(6,*)'ERROR:  the cubic Lagrange algorithm requires 4 source points.'
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
!$OMP parallel do private(i,x1,x2,x3,x4,z1,z2,z3,z4,cx1,cx2,cx3,cx4,cx5,cx6,a1,a2,a3,a4,xt)
  do vt=1,destNumLevels                 ! for each target vertical point
    do i=1,numInterpSets                ! for each horizontal point

      ! Identify the 4 known points that straddle the sought point
      firstLevel = max( min(posnDestInSrc(i,vt)-1, srcNumLevels-3), 1 )

      ! These are the four levels that will be used to interpolate this point
      x1=vLevelSource(i, firstLevel)
      x2=vLevelSource(i, firstLevel+1)
      x3=vLevelSource(i, firstLevel+2)
      x4=vLevelSource(i, firstLevel+3)

      ! These are the values of the state at those levels
      z1=stateSource(i, firstLevel)
      z2=stateSource(i, firstLevel+1)
      z3=stateSource(i, firstLevel+2)
      z4=stateSource(i, firstLevel+3)

      ! Calculate some differences between those levels
      cx1 = 1.0/(x2-x1)
      cx2 = 1.0/(x3-x1)
      cx3 = 1.0/(x3-x2)
      cx4 = 1.0/(x4-x1)
      cx5 = 1.0/(x4-x2)
      cx6 = 1.0/(x4-x3)

      ! Calculate some co-efficients for the interpolation
      a1 = z1
      a2 = cx1*(z2-z1)
      a3 = cx2*(cx3*(z3-z2)-cx1*(z2-z1))
      a4 = cx4*(  cx5*(cx6*(z4-z3) - cx3*(z3-z2)) &
                - cx2*(cx3*(z3-z2) - cx1*(z2-z1)) )

      ! Finally, plug in the position of the sought value
      xt = vLevelDestn(i,vt)
      stateDestn(i, vt) = a1+(xt-x1)*(a2+(xt-x2)*(a3+a4*(xt-x3)))

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

end subroutine ! Interp1D_CubicLagrange_X
