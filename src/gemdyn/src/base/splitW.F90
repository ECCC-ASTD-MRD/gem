      subroutine splitW (F_me,F_rank,F_glb,F_g0,F_lcl,F_start,F_end)
      implicit none
      
      integer, intent(IN ) :: F_me, F_rank,F_glb,F_g0
      integer, intent(OUT) :: F_lcl,F_start,F_end

      integer :: mpx,irest
!     
!--------------------------------------------------------------------
      mpx   = mod( F_me, F_rank )
      F_lcl = F_glb / F_rank
      irest = F_glb - F_lcl * F_rank
      F_start   = mpx * F_lcl + F_g0
      if ( mpx < irest ) then
         F_lcl   = F_lcl + 1
         F_start = F_start + mpx
      else
         F_start = F_start + irest
      endif
      F_end= F_start + F_lcl - 1
!     
!     
!--------------------------------------------------------------------
!
      return
      end subroutine splitW
