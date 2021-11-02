!!!s/r Extrap1D_Fixed - Extrapolate, to a fixed (selectable) value
!subroutine Extrap1D_Fixed_X  &
!subroutine Extrap1D_Fixed_X8 &
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
  ! stateDerivDestn, flux, numExtArraysIn, numExtArraysOut, ExtArraysIn, 
  ! ExtArraysOut
!
!AUTHOR
!     J.W. Blezius NOV 2008 First revision for this file
!
!REVISION
! v1_4    Blezius J.W. NOV 2008 - add fixed extrapolation (this file)
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Extrapolate values at vertical levels above (below) the highest (lowest)
!        vertical level.  This algorithm is very simple.  The extrapolated value
!        is set to the user-selected fixed value.  The extrapolated derivative of
!        the state observable is, by virtue of the fixed state value, equal to
!        zero.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"
!
!NOTES
!        In the context of this routine the parameters, extrapGuideDown and
!        extrapGuideUp, are the fixed values below and above (respectively) the
!        vLevelSource values.  The units of the fixed values are:  state units.
!
!        It is assumed that stateDerivSource and particularly stateDerivDestn are
!        real arrays.  The values need not make any sense since they are not used
!        beyond copying them.  However, stateDerivDestn is intended for output,
!        and as such it had better be legal to write to it.
!
!!

  integer i              ! horizontal loop index
  integer vt             ! vertical loop index into the target grid

  integer :: indexLowestLevel, &        ! Index of the lowest vertical level
             indexHighestLevel, &       ! Index of the highest vertical level
             highestPossiblePosn, &     ! Highest value for posnDestInSrc
             lowestPossiblePosn         ! Lowest value for posnDestInSrc


  ! Determine which index represents the highest and lowest vertical levels
  !
  if( vLevelSource(1,2) > vLevelSource(1,1) ) then
    ! The vertical level values INcrease with the index.
    indexLowestLevel    = 1
    indexHighestLevel   = srcNumLevels
    highestPossiblePosn = srcNumLevels - 1
    lowestPossiblePosn  = 1
  else ! not increasing levels
    ! The vertical level values DEcrease with the index.
    indexLowestLevel = srcNumLevels
    indexHighestLevel = 1
    highestPossiblePosn = 1
    lowestPossiblePosn  = srcNumLevels - 1
  endif ! increasing levels

  ! Each destination point must be checked because no particular order has been
  ! assumed for these points.
  !
!$OMP parallel do private(i)
  do vt=1,destNumLevels                 ! for each target vertical point
    do i=1,numInterpSets                ! for each horizontal point

                                        ! This test gives faster code
      if( extrapEnableDown .and. posnDestInSrc(i,vt) == lowestPossiblePosn) then
        if( vLevelDestn(i,vt) < vLevelSource(i,indexLowestLevel) ) then
          ! extrapolate below the lowest source vertical level
          stateDestn(i,vt) =   extrapGuideDown
        endif

      else if(      extrapEnableUp &
              .and. posnDestInSrc(i,vt) == highestPossiblePosn) then
        if( vLevelDestn(i,vt) > vLevelSource(i,indexHighestLevel) ) then
          ! extrapolate above the highest source vertical lavel
          stateDestn(i,vt) = extrapGuideUp
        end if
      end if ! extrapEnableDown
    end do ! i
  end do ! vt
!$OMP END parallel do


end subroutine ! Extrap1D_Fixed_X
