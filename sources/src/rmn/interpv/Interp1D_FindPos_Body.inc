!!!s/r Interp1D_FindPos - Find array indices between which interpolation will
!                         occur
!subroutine Interp1D_FindPos  &
!subroutine Interp1D_FindPos8 &
                      (numInterpSets, srcNumLevels, destNumLevels, &
                       src_ijDim, dst_ijDim, &

                       vLevelSource, posnDestInSrc, vLevelDestn &
                      )
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - this routine was copied from find_pos1
! v1_1    Blezius J.W. SEP 2002 - make all input dim's the same; ditto for output
!         Blezius J.W. DEC 2015 - take advantage of OpenMP
!
!OBJECT
!        A part of the 1-D interpolation package, this function compares an array
!        of source co-ordinte values (usually vertical levels), at which state
!        values (usually temperature, wind, etc.) are known, with an array of
!        destination co-ordinate values, at which state values are sought.  For
!        each element of the destination array, the result of this comparison is
!        an index in the source array whose co-ordinate value (vertical level) is
!        just before that of the destination array.  The binary search repeated
!        on a series of independent interpolation sets.
!
!ARGUMENTS
  use Interp1D_Constants
  implicit none

  ! elements used (containing data):
  integer, intent(in) :: numInterpSets  ! no. sets (horiz. pts) for interpolation
  integer, intent(in) :: srcNumLevels   ! no. pts in vLevelSource containing data
  integer, intent(in) :: destNumLevels  ! no. pts in vLevelDestn containing data

  ! elements dimensioned:
  integer, intent(in) :: src_ijDim      ! horizontal dimension of source arrays
  integer, intent(in) :: dst_ijDim      ! horizontal dimension of dest'n arrays

                                        ! source co-ord's (vertical levels)
  real(real48), dimension(src_ijDim, srcNumLevels), intent(in) :: vLevelSource

                                        ! position (index) of vLevelSource that
                                        !   is just before the vLevelDestn
  integer, dimension(dst_ijDim, destNumLevels), intent(out) :: posnDestInSrc

                                        ! destination co-ord's (vert. levels)
  real(real48), dimension(dst_ijDim, destNumLevels), intent(in) :: vLevelDestn
!
!NOTES
!        It is the responsibility of the user to ensure that the arrays that he
!        supplies are large enough to contain the number of points that he
!        indicates; specifically, numInterpSets is the minimum first dimension
!        on all arrays.
!
!        It is assumed that the source vertical levels are in either ascending or
!        descending order.
!
!        The destination vertical levels may be in any order.  The determination
!        of the location of each destination level is completely independent of
!        all the others.
!
!!

  real, dimension(numInterpSets) :: indexReal! vertical index estimate, as a real
  integer indexInt                      ! index estimate (indexReal truncated)
  real uncertainty                      ! uncertainty of the value of indexReal
  integer t                             ! index of the target levels
  integer s                             ! index of the interpolation set


!$OMP parallel do private(s,indexReal,uncertainty,indexInt)
  do t=1,destNumLevels                  ! for each (vertical) target level

    ! initial guess is in the middle of the field, with an uncertainty of half
    ! the field; i.e. the sought index is somewhere between the lowest and
    ! highest indices
    !
    do s = 1, numInterpSets
      indexReal(s) = 0.5 * (srcNumLevels + 1)
    end do
    uncertainty = 0.5 * (srcNumLevels - 1)


    if( vLevelSource(1,2) > vLevelSource(1,1) ) then
      ! The vertical level values INcrease with the index.

      ! Because the index is an integer, an uncertainty less than 0.5 is as good
      ! as an uncertainty of 0.  However, due to a subtlety, an uncertainty less
      ! than 1 is also sufficient.  It should be noted that the two adjustments
      ! of indexReal have a special effect in the case where indexInt is just
      ! below the sought value:  indexReal is both incremented and decremented.
      ! However, it is not decremented if indexInt+1 is also below the sought
      ! value.  This means that on the final pass (with an uncertainty of 0.95 or
      ! so), indexInt+1 has been sampled to see whether it is the right answer. 
      ! The result is that it is not necessary to reduce the uncertainty to less
      ! than 0.5; less than 1 is sufficient.
      !
      do while(uncertainty > 0.95)

        ! It is known at this point that the sought index is somewhere between
        ! indexReal+uncertainty and indexReal-uncertainty. To state the obvious,
        ! the sought index must be in either the upper half or the lower half of
        ! this range.  Whichever half it is in, the next guess is in the middle
        ! of that half (i.e. the guess is either increased or decreased by half
        ! of the previous uncertainty), give or take a half of the previous
        ! uncertainty.
        !
        uncertainty = 0.5 * uncertainty
        do s = 1, numInterpSets
          indexInt = indexReal(s)
          if(vLevelDestn(s,t) >= vLevelSource(s, indexInt)) &
                                          indexReal(s)=indexReal(s) + uncertainty
          if(vLevelDestn(s,t) <= vLevelSource(s, indexInt+1)) &
                                          indexReal(s)=indexReal(s) - uncertainty
        end do ! s
      end do ! uncertainty .gt. 0.95

    else ! not increasing levels
      ! The vertical level values DEcrease with the index.

      ! This code block is the mirror image of the previous increasing-level
      ! block.

      do while(uncertainty > 0.95)
        uncertainty = 0.5 * uncertainty
        do s = 1, numInterpSets
          indexInt = indexReal(s)
          if(vLevelDestn(s,t) >= vLevelSource(s, indexInt)) &
                                          indexReal(s)=indexReal(s) - uncertainty
          if(vLevelDestn(s,t) <= vLevelSource(s, indexInt+1)) &
                                          indexReal(s)=indexReal(s) + uncertainty
        end do ! s
      end do ! uncertainty .gt. 0.95
    endif ! increasing levels

    ! The indices just below the current target_p%vLevel have been found for each
    ! (horizontal) point in the set.  Record them.
    posnDestInSrc(1:numInterpSets, t)=indexReal
  end do ! t
!$OMP END parallel do

end subroutine ! Interp1D_FindPos
