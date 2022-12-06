      integer*4 function get_file_unit (lu_max)
!
!   get_file_unit returns a unit number that is not in use
      integer*4 lu_max,  lu, m, iostat
      logical   opened
!
      m = lu_max  ;  if (m < 1) m = 97
      do lu = m,1,-1
         inquire (unit=lu, opened=opened, iostat=iostat)
         if (iostat.ne.0) cycle
         if (.not.opened) exit
      end do
!
      get_file_unit = lu
      return
      end function get_file_unit
