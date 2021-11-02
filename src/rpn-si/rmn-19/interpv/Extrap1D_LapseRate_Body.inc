!!!s/r Extrap1D_LapseRate - Extrapolate, based on lapse rates
!subroutine Extrap1D_LapseRate_X  &
!subroutine Extrap1D_LapseRate_X8 &
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
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - new code
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output
! v1_2    Blezius J.W. OCT 2002 - extrapolation failed miserably in the case
!                                 where the vertical levels are in DEcreasing
!                                 order
! v1_3    Blezius J.W. OCT 2003 - extend the interface, but ignore the extra data
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Extrapolate values at vertical levels above (below) the highest (lowest)
!        vertical level.  This algorithm is simple.  The extrapolated value
!        of the state observable is calculated from that at the vertical level
!        that is closest to the one sought, plus the difference that is
!        determined by the product of the lapse rate and the difference in
!        vertical levels. The extrapolated derivative of the state observable is
!        equal to the supplied lapse rate (by definition of lapse rate).
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"
!
!NOTES
!        In the context of this routine the parameters, extrapGuideDown and
!        extrapGuideUp, are the lapse rates below and above (respectively) the
!        vLevelSource values.  The units of the lapse rates are:  state units /
!        vLevel units.
!
!        It is assumed that stateDerivSource and particularly stateDerivDestn are
!        real arrays.  The values need not make any sense since they are not used
!        beyond copying them.  However, stateDerivDestn is intended for output,
!        and as such it had better be legal to write to it.
!
!        Derivatives are simply not extrapolated.  The estimates of their values
!        would be too poor to be of use anyway.
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
          stateDestn(i,vt) =   stateSource(i,indexLowestLevel) &
                         + extrapGuideDown * (  vLevelDestn (i,vt) &
                                              - vLevelSource(i,indexLowestLevel)&
                                             )
        endif

      else if(      extrapEnableUp &
              .and. posnDestInSrc(i,vt) == highestPossiblePosn) then
        if( vLevelDestn(i,vt) > vLevelSource(i,indexHighestLevel) ) then
          ! extrapolate above the highest source vertical lavel
          stateDestn(i,vt) = stateSource(i,indexHighestLevel) &
                          + extrapGuideUp * (  vLevelDestn(i,vt) &
                                             - vLevelSource(i,indexHighestLevel)&
                                            )
        end if
      end if ! extrapEnableDown
    end do ! i
  end do ! vt
!$OMP END parallel do


end subroutine ! Extrap1D_LapseRate_X
