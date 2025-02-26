
!/@*
function bvstrip(F_name) result(F_norm)
  use phymem, only: PHY_NAMELEN
  implicit none
  !@objective Strip bus variable names of extraneous information
  !@arguments
  character(len=*), intent(in) :: F_name        !Bus variable name
  !@return
  character(len=PHY_NAMELEN) :: F_norm          !Normalized name
  !@author Ron McTaggart-Cowan, 2014-03
  !@description
  integer :: i
  character(len=1), dimension(2), parameter :: ext = (/'w','W'/)
  !*@/
  F_norm = F_name
  if (len_trim(F_norm) < 2) return  !minimum length for separator
  do i=1,size(ext)
     if (F_name(len_trim(F_name)-1:) == ','//ext(i)) &
          F_norm = F_name(:len_trim(F_name)-2)
  enddo
  return
end function bvstrip
