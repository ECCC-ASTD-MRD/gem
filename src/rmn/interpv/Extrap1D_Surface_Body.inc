!!!s/r Extrap1D_Surface - extrapolation:  accounts for the physics of the Earth's
!                                         surface
! subroutine Extrap1D_Surface_X  &
! subroutine Extrap1D_Surface_X8 &
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
  ! stateDerivDestn, extrapGuideDown, extrapGuideUp
!
!AUTHOR
!     J.W. Blezius OCT 2003 new kind of extrapolation
!
!REVISION
! v1_3    Blezius J.W. OCT 2003 - new code
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        Extrapolate values at vertical levels below the lowest vertical level,
!        not counting the Earth's surface itself.  Technically, that makes this
!        an interpolation, or perhaps a guided extrapolation of the results from
!        the forecast model.
!
!        This extrapolation takes into account the physics that governs diffusion
!        near the Earth's surface.  Thus, it is not applied to extrapolation
!        above the highest vertical level.  Moreover, it should not be used for a
!        general mathematical extrapolation; the algorithm makes sense only when
!        considering the surface of the Earth.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none
#include "Interp1D_Arguments_X.inc"
!
!NOTES
!        Each routine of the set of interpolation routines takes exactly the same
!        parameters, allowing the user to use the routines interchangeably, as
!        with a pointer.
!
!        Derivatives are simply not extrapolated.  The estimates of their values
!        would be too poor to be of use anyway.
!
!        This routine contains no physics knowledge itself.  The physics is
!        obtained from the callback function that is passed as an argument.
!
!        N.B.:  The physics of the callback function assumes that the source and
!               destination vertical levels are distances (in m) above the
!               Earth's surface.
!
!        N.B.:  stateSource(indexSurface) is assumed to contain true values, not
!               those that are obtained from GEM.
!
!!

  integer i                             ! horizontal loop index
  integer vt                            ! vertical loop index into target grid

  integer :: indexLowestLevel, &        ! Index of the lowest vertical level
             lowestPossiblePosn, &      ! Lowest value for posnDestInSrc
             indexSurface               ! Index of level that is Earth's surface

  real(real48), dimension(src_ijDim) :: f
  real(real48), dimension(ubound(ExtArraysIn,1)) :: z0, ilmo, hBound, ft

  ! Verify the number of extension arguments that were supplied
  if(     numExtArraysIn  > MAX_EXT_ARRAYS_IN &
     .or. numExtArraysOut > MAX_EXT_ARRAYS_OUT) then
    write(6,*)'ERROR:  too many extension arrays supplied for this revision.'
    return
  end if

  if(numExtArraysIn < 4) then
    write(6,*)'ERROR:  too few extension arrays supplied:  ', numExtArraysIn
    return
  end if

  ! Verify that enough source levels are supplied; the algorithm requires 2
  if(srcNumlevels < 2) then
    write(6,*)'ERROR:  the surface algorithm requires 2 source points.'
    return
  end if

  ! Create array sections for the required extension arrays
  z0    = ExtArraysIn(:, iZ0)           ! roughness length
  ilmo  = ExtArraysIn(:, iILMO)         ! inverse of Monin-Obukhov length
  hBound= ExtArraysIn(:, iH_BOUND)      ! height of the boundary layer

                                        ! vertical flux at the surface 
                                        ! [in (units of the state)
                                        !           / (units returned by flux())]
  ft    = ExtArraysIn(:, iFT)

  ! Determine which index represents the highest and lowest vertical levels
  !
  if( vLevelSource(1,2) > vLevelSource(1,1) ) then
    ! The vertical level values INcrease with the index.
    ! N.B.:  The lowest level is assumed to be the surface; take the next level
    indexLowestLevel    = 2
    lowestPossiblePosn  = 1
    indexSurface        = 1
  else ! not increasing levels
    ! The vertical level values DEcrease with the index.
    ! N.B.:  The lowest level is assumed to be the surface; take the next level
    indexLowestLevel    = srcNumLevels - 1
    lowestPossiblePosn  = srcNumLevels - 1
    indexSurface        = srcNumLevels
  endif ! increasing levels

  ! Each destination point must be checked because no particular order has been
  ! assumed for these points.
  !
!$OMP parallel do private(i,f)
  do vt=1,destNumLevels                 ! for each target vertical point
    call flux(f, vLevelDestn(:,vt), z0, ilmo,hBound,numInterpSets)!calculate f(i)
    do i=1,numInterpSets                ! for each horizontal point

                                        ! This test gives faster code
      if( extrapEnableDown .and. posnDestInSrc(i,vt) == lowestPossiblePosn) then
        if( vLevelDestn(i,vt) < vLevelSource(i,indexLowestLevel) ) then
          ! extrapolate below the second-lowest source vertical level
          stateDestn(i,vt) = stateSource(i,indexSurface) + ft(i)*f(i)
        endif
      end if ! extrapEnableDown
    end do ! i
  end do ! vt
!$OMP END parallel do

end subroutine ! Extrap1D_Surface_X
