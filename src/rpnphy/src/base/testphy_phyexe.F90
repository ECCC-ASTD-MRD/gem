subroutine testphy_phyexe(d,f,v,dsiz,fsiz,vsiz,trnch,kount,ni,nk)
  ! Test double for RPN Physics parameterizations: increment time-plus 
  ! temperatures throughout the column by 1K and put 100 into the
  ! boundary layer height.
  use phybus
  implicit none

  ! Architecture-dependent content
!!!#include <arch_specific.hf>

  ! Input arguments
  integer, intent(in) :: dsiz,fsiz,vsiz,trnch,kount,ni,nk
  real, intent(inout) :: d(dsiz),f(fsiz),v(vsiz)

  ! External declarations

  ! Local variables
  integer :: i,k

  ! Increment temperatures
  do k=1,nk
     do i=1,ni
        d(tplus+(k-1)*ni+(i-1)) = d(tplus+(k-1)*ni+(i-1)) + 1.
     enddo
  enddo

  ! Specify boundary layer height
  do i=1,ni
     f(h+(i-1)) = 100.
  enddo

  if (ni<0) print *,'called testphy_phyexe with:',d(1),f(1),v(1),dsiz,fsiz,vsiz,trnch,kount,ni,nk
  ! End of test double
  return
end subroutine testphy_phyexe
