
module testtiming_mod
   public

#include <rmn/msg.h>

contains

   subroutine testtiming_test1(sum1, val1, val2, ni, nj, nk)
      integer, intent(in) :: ni, nj, nk
      real :: sum1(ni,nj,nk), val1(ni,nj,nk), val2(ni,nj,nk)
      integer :: j1,j2, j3
      do j3 = 1,nk
         do j2 = 1,nj
            do j1 = 1,ni
               sum1(j1,j2,j3) = 0.
            enddo
         enddo
      enddo
      do j3 = 1,nk
         do j2 = 1,nj
            do j1 = 1,ni
               sum1(j1,j2,j3) = sum1(j1,j2,j3) + val1(j1,j2,j3) + val2(j1,j2,j3)
            enddo
         enddo
      enddo
      return
   end subroutine testtiming_test1

   subroutine testtiming_test2(sum1, val1, val2, ni, nj, nk)
      integer, intent(in) :: ni, nj, nk
      real :: sum1(ni,nj,nk), val1(ni,nj,nk), val2(ni,nj,nk)
      sum1 = val1 + val2
      return
   end subroutine testtiming_test2


   subroutine testtiming_test1p(sum1, val1, val2, ni, nj, nk)
      integer, intent(in) :: ni, nj, nk
      real,pointer :: sum1(:,:,:), val1(:,:,:), val2(:,:,:)
      integer :: j1,j2, j3
      do j3 = 1,nk
         do j2 = 1,nj
            do j1 = 1,ni
               sum1(j1,j2,j3) = 0.
            enddo
         enddo
      enddo
      do j3 = 1,nk
         do j2 = 1,nj
            do j1 = 1,ni
               sum1(j1,j2,j3) = sum1(j1,j2,j3) + val1(j1,j2,j3) + val2(j1,j2,j3)
            enddo
         enddo
      enddo
      return
   end subroutine testtiming_test1p

   subroutine testtiming_test2p(sum1, val1, val2, ni, nj, nk)
      integer, intent(in) :: ni, nj, nk
      real,pointer :: sum1(:,:,:), val1(:,:,:), val2(:,:,:)
      if (ni == 0) print *, ni, nj, nk  !# Avoid compiler unused warning/remark
      sum1 = val1 + val2
      return
   end subroutine testtiming_test2p


   subroutine testtiming_test1po(sum1, val1, val2, ni, nj, nk)
      integer, intent(in) :: ni, nj, nk
      real,pointer :: sum1(:,:,:), val1(:,:,:), val2(:,:,:)
      integer :: j1,j2, j3
!$omp parallel private(j1,j2,j3)
!$omp do
     do j3 = 1,nk
         do j2 = 1,nj
            do j1 = 1,ni
               sum1(j1,j2,j3) = 0.
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel
!$omp parallel private(j1,j2,j3)
!$omp do
      do j3 = 1,nk
         do j2 = 1,nj
            do j1 = 1,ni
               sum1(j1,j2,j3) = sum1(j1,j2,j3) + val1(j1,j2,j3) + val2(j1,j2,j3)
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel
      return
   end subroutine testtiming_test1po

   !#NOTE: Apparently there is no wide compiler support for OMP with array notation, it is up to the vendor.
!!$   subroutine testtiming_test2po(sum1, val1, val2, ni, nj, nk)
!!$      integer, intent(in) :: ni, nj, nk
!!$      real,pointer :: sum1(:,:,:), val1(:,:,:), val2(:,:,:)
!!$!$omp parallel
!!$      sum1 = val1 + val2
!!$!$omp end parallel
!!$      return
!!$   end subroutine testtiming_test2po

end module testtiming_mod



subroutine testtiming()
   use testtiming_mod
   implicit none
#include <rmnlib_basics.hf>
 
  integer, parameter :: niter=30
  integer, parameter :: ni=2000, nj=1000, nk=100
  integer :: myproc, mynum, j1, j2, j3
  character(len=32) :: msg, name
  real,target :: sum1(ni,nj,nk), val1(ni,nj,nk), val2(ni,nj,nk)
  real, pointer :: sum1p(:,:,:), val1p(:,:,:), val2p(:,:,:)

  sum1p => sum1
  val1p => val1
  val2p => val2

  myproc=0

  do j3 = 1,nk
     do j2 = 1,nj
        do j1 = 1,ni
           val1(j1,j2,j3) = float(j2*2+j1+j3)
           val2(j1,j2,j3) = float(j1*2+j2+j3)
        enddo
     enddo
  enddo

  msg='testtiming'
  call timing_init2(myproc, msg)

  call timing_start2(1, 'All', 1)
  do j1 = 1,niter
     print *,'iter:',j1,'/',niter ; call flush(6)
     mynum=2; name='do'
     call timing_start2(mynum, name, 1)
     call testtiming_test1(sum1p, val1p, val2p, ni, nj, nk)
     call timing_stop(mynum)
     print *,'sum2=',maxval(sum1) ; call flush(6)

     mynum=3; name='arr'
     call timing_start2(mynum, name, 1)
     call testtiming_test2(sum1p, val1p, val2p, ni, nj, nk)
     call timing_stop(mynum)
     print *, 'sum3=', maxval(sum1) ; call flush(6)

     mynum=4; name='do ptr'
     call timing_start2(mynum, name, 1)
     call testtiming_test1p(sum1p, val1p, val2p, ni, nj, nk)
     call timing_stop(mynum)
     print *,'sum4=',maxval(sum1) ; call flush(6)

     mynum=5; name='arr ptr'
     call timing_start2(mynum, name, 1)
     call testtiming_test2p(sum1p, val1p, val2p, ni, nj, nk)
     call timing_stop(mynum)
     print *,'sum5=',maxval(sum1) ; call flush(6)

     mynum=6; name='do ptr omp'
     call timing_start2(mynum, name, 1)
     call testtiming_test1po(sum1p, val1p, val2p, ni, nj, nk)
     call timing_stop(mynum)
     print *,'sum6=',maxval(sum1) ; call flush(6)

!!$     mynum=7; name='arr ptr omp'
!!$     call timing_start2(mynum, name, 1)
!!$     call testtiming_test2po(sum1p, val1p, val2p, ni, nj, nk)
!!$     call timing_stop(mynum)
!!$     print *,'sum7=',maxval(sum1) ; call flush(6)

  enddo
  call timing_stop(1)

  call timing_terminate2(myproc, msg)
  return
end subroutine testtiming
