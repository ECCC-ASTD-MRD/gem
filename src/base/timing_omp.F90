   
module timing_omp
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_getenv, clib_toupper, CLIB_IS_OK
   implicit none
   private
      
#include <rmnlib_basics.hf>

   integer, external :: omp_get_thread_num
   DOUBLE PRECISION, external :: omp_get_wtime

   public :: timing_init_omp, timing_start1, timing_start_omp, &
        timing_stop1, timing_stop_omp, timing_terminate1

   integer, parameter :: MAX_instrumented=512
   integer, parameter :: MAX_threads=128

   integer(INT64),save :: timer_cnt(MAX_instrumented, MAX_threads)
   integer,save :: timer_level(MAX_instrumented)

   real(REAL64),save :: tb(MAX_instrumented, MAX_threads)
   real(REAL64),save :: sum_tb(MAX_instrumented, MAX_threads)

   character(len=16) :: Timing_S, nam_subr_S(MAX_instrumented)
   
contains

   
   subroutine timing_init_omp(myproc, msg)
      implicit none
      !@arguments
      character(len=*), intent(in) :: msg
      integer, intent(in) :: myproc
      !@author M. Desgagne   -- Winter 2012 --
      !@revision
      ! v4_40 - Desgagne - initial version
      ! v4_80 - Desgagne - introduce timer_level and timer_cnt

      character(len=16) :: dumc_S
      !-------------------------------------------------------------------
      !$omp single
      if (.not.CLIB_IS_OK(clib_getenv('TMG_ON', Timing_S))) Timing_S = 'NO'
      call low2up(Timing_S, dumc_S)
      Timing_S = dumc_S

      if (Timing_S == 'YES') call tmg_init(myproc, msg)

      sum_tb = 0.D0; timer_cnt = 0 ; timer_level = 0 ; nam_subr_S = ''
      !$omp end single nowait
      !-------------------------------------------------------------------
      return
   end subroutine timing_init_omp


   subroutine timing_start1(mynum, myname_S, mylevel)
      implicit none
      integer, intent(in) :: mynum, mylevel
      character(len=*), intent(in) :: myname_S

      integer :: t, mynum2
      !--------------------------------------------------------
      call timing_start_omp(mynum, myname_S, mylevel, 1)
      !--------------------------------------------------------
      return
   end subroutine timing_start1


   subroutine timing_start_omp(mynum, myname_S, mylevel,mythread)
      implicit none
      integer, intent(in) :: mynum, mylevel
      character(len=*), intent(in) :: myname_S
      integer, intent(in),optional :: mythread

      integer :: t, mynum2
      !--------------------------------------------------------
      if (present(mythread)) then
         t = mythread
      else
         t = min(max(0, omp_get_thread_num()) + 2, MAX_threads)
      endif
      mynum2 = max(1, min(mynum, MAX_instrumented))
      if (mynum /= mynum2) then
         print *, 'WARNING: (timing_start_omp) called with mnum=', &
              mynum, ' > MAX_instrumented=', MAX_instrumented
      endif

!$omp single
      if (Timing_S == 'YES') call tmg_start(mynum2, myname_S)
      nam_subr_S(mynum2)  = myname_S
      timer_level(mynum2) = mylevel
!$omp end single nowait

      tb(mynum2,t) = omp_get_wtime()
      timer_cnt(mynum2,t) = timer_cnt(mynum2,t) + 1
      !--------------------------------------------------------
      return
   end subroutine timing_start_omp


   subroutine timing_stop1(mynum)
      implicit none
      integer, intent(in) :: mynum
      !--------------------------------------------------------
      call timing_stop_omp(mynum,1)
      !--------------------------------------------------------
      return
   end subroutine timing_stop1


   subroutine timing_stop_omp(mynum,mythread)
      implicit none
      integer, intent(in) :: mynum
      integer, intent(in),optional :: mythread

      integer :: t, mynum2
      !--------------------------------------------------------
      mynum2 = max(1, min(mynum, MAX_instrumented))
      if (present(mythread)) then
         t = mythread
      else
         t = min(max(0, omp_get_thread_num()) + 2, MAX_threads)
      endif

!$omp single
      if (Timing_S == 'YES') call tmg_stop(mynum2)
!$omp end single nowait

      sum_tb(mynum2,t) = sum_tb(mynum2,t) + (omp_get_wtime() - tb(mynum2,t))
      !--------------------------------------------------------
      return
   end subroutine timing_stop_omp

   
   subroutine timing_terminate1(myproc, msg)
      implicit none
      !@arguments
      character(len=*), intent(in) :: msg
      integer, intent(in) :: myproc
      !@author M. Desgagne   -- Winter 2012 --

      character(len=64) :: fmt,nspace,nspace2,nspace3,tmp1_s,tmp2_s
      logical flag(MAX_instrumented)
      integer i,j,elem,lvl,lvlel(0:100),err,maxlen
      integer(INT64),dimension(MAX_instrumented) :: timer_cnt2
      real(REAL64),dimension(MAX_instrumented) :: sum_tb_mn, sum_tb_mx, sum_tb_mn2, sum_tb_mx2
      real(REAL64) :: mymax, mymin
      real :: sum_tb_per_mn, sum_tb_per_mx
      !--------------------------------------------------------
      if (Timing_S=='YES') call tmg_terminate ( myproc, msg )

      maxlen = 0
      do i=1,MAX_instrumented
         mymax = maxval(sum_tb(i,2:MAX_threads))
         sum_tb_mx(i) = sum_tb(i,1) + mymax
         mymin = mymax
         do j=2,MAX_threads
            if (sum_tb(i,j) > 0.) mymin = min(mymin, sum_tb(i,j))
         enddo
         sum_tb_mn(i)   = sum_tb(i,1)    + mymin
         timer_cnt(i,1) = timer_cnt(i,1) + maxval(timer_cnt(i,2:MAX_threads))
         maxlen = max(maxlen, len_trim(nam_subr_S(i)))
      enddo

      sum_tb_mn2 = 0.D0
      sum_tb_mx2 = 0.D0
      timer_cnt2 = 0
      call rpn_comm_reduce(sum_tb_mn,    sum_tb_mn2,    MAX_instrumented, &
           RPN_COMM_REAL8,    RPN_COMM_MIN, RPN_COMM_MASTER, RPN_COMM_GRID, err)
      call rpn_comm_reduce(sum_tb_mx,    sum_tb_mx2,    MAX_instrumented, &
           RPN_COMM_REAL8,    RPN_COMM_MAX, RPN_COMM_MASTER, RPN_COMM_GRID, err)
      call rpn_comm_reduce(timer_cnt(:,1), timer_cnt2, MAX_instrumented, &
           RPN_COMM_INTEGER8, RPN_COMM_MAX, RPN_COMM_MASTER, RPN_COMM_GRID, err)

      if (myproc.ne.0) return

      !#TODO: Add Mean and Var/Std to timings

      write(6,'(a)') '________________________________________________________________________________________'
      write(6,'(a)') '|____ TIMINGS __________________________________________________________________________|'
      write(6,'(a)') '|   |                               |   Wallclock [%]| Wallclock [Sec]         |        |'
      write(6,'(a)') '| ID| NAME                          |   Min%:    Max%| MinSec     : MaxSec     |  Count |'
      write(6,'(a)') '|---|-------------------------------|----------------|-------------------------|--------|'
      flag=.false.
      mymax = maxval(sum_tb_mx2)
      write (nspace3,'(i3)') max(5, maxlen)
      do i = 1,MAX_instrumented
         lvl= 0 ; elem= i
55       if ( (trim(nam_subr_S(elem)).ne.'') .and. (.not.flag(elem)) ) then

            err = clib_toupper(nam_subr_S(elem))
            sum_tb_per_mn = real(sum_tb_mn2(elem)) / max(0.00001,real(mymax))
            sum_tb_per_mx = real(sum_tb_mx2(elem)) / max(0.00001,real(mymax))
            write (nspace,'(i3)') 3*lvl+1
            write (nspace2,'(i3)') max(1,18-(3*lvl+1))
            fmt = '(a,i3,a,'//trim(nspace)//'x,a'//trim(nspace3)//','// &
                 trim(nspace2)//'x,a,f6.2,a,f6.2,a)'
            write(tmp1_S, fmt) &
                 "|", elem, "|", nam_subr_S(elem), &
                 "|", 100.*sum_tb_per_mn, "%: ", 100.*sum_tb_per_mx, "%|"
            fmt = '(1pe10.4,a,1pe10.4,a,i8,a)'
            write(tmp2_S, fmt) &
                 sum_tb_mn2(elem), " : ", sum_tb_mx2(elem), &
                 ' |', timer_cnt2(elem), '|'
            write(6, '(a,1x,a)') trim(tmp1_S),trim(tmp2_S)
            flag(elem) = .true. ; lvlel(lvl) = elem
65          do j = 1,MAX_instrumented
               if ((timer_level(j) .eq. elem) .and. (.not.flag(j)) )then
                  lvl= lvl+1
                  elem= j
                  goto 55
               endif
            end do
            lvl= lvl - 1
            if (lvl .ge. 0) then
               elem= lvlel(lvl)
               goto 65
            endif
         endif
      enddo

      write(6,'(a)') '________________________________________________________________________________________'
      !--------------------------------------------------------
      return
   end subroutine timing_terminate1
   
end module timing_omp


subroutine timing_init2(myproc, msg)
   use timing_omp
   implicit none
   character(len=*), intent(in) :: msg
   integer, intent(in) :: myproc
   call timing_init_omp(myproc, msg)
   return
end subroutine timing_init2


subroutine timing_start2(mynum, myname_S, mylevel)
   use timing_omp
   implicit none
   integer, intent(in) :: mynum,mylevel
   character(len=*), intent(in) :: myname_S
   call timing_start1(mynum, myname_S, mylevel)
   return
end subroutine timing_start2


subroutine timing_stop(mynum)
   use timing_omp
   implicit none
   integer, intent(in) :: mynum
   call timing_stop1(mynum)
   return
end subroutine timing_stop


subroutine timing_terminate2(myproc, msg)
   use timing_omp
   implicit none
   character(len=*), intent(in) :: msg
   integer, intent(in) :: myproc
   call timing_terminate1(myproc, msg)
   return
end subroutine timing_terminate2

