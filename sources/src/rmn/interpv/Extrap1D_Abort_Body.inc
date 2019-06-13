!!!s/r Extrap1D_Abort - Abort if an extrapolation is required
!subroutine Extrap1D_Abort_X  &
!subroutine Extrap1D_Abort_X8 &
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
!     J.W. Blezius MAY 2014 First revision for this file
!
!REVISION
!         Blezius J.W. MAY 2014 - add abortion on extrapolation (this file)
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Test for the need to extrapolate values at vertical levels above (below)
!        the highest (lowest) vertical level.  If extrapolation is required, do
!        not extrapolate.  Abort execution instead.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"
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
          ! extrapolation is required:  abort execution
          write(*,*)
          write(*,*)"ERROR:  Attempting extrapolation to level ", &
                    vLevelDestn(i,vt)
          write(*,*)"        This is below the lowest level ", &
                    vLevelSource(i,indexLowestLevel)
          write(*,*)"        Abortion on extrapolation has been selected ... ", &
                             "ABORTING"
          call exit(2)
        endif

      else if(      extrapEnableUp &
              .and. posnDestInSrc(i,vt) == highestPossiblePosn) then
        if( vLevelDestn(i,vt) > vLevelSource(i,indexHighestLevel) ) then
          ! extrapolation is required:  abort execution
          write(*,*)
          write(*,*)"ERROR:  Attempting extrapolation to level ", &
                    vLevelDestn(i,vt)
          write(*,*)"        This is above the highest level ", &
                    vLevelSource(i,indexHighestLevel)
          write(*,*)"        Abortion on extrapolation has been selected ... ", &
                             "ABORTING"
          call exit(2)
        end if
      end if ! extrapEnableDown
    end do ! i
  end do ! vt
!$OMP END parallel do


end subroutine ! Extrap1D_Abort_X
