subroutine test_phyutils
   use str_mod, only: str_normalize, str_toint
   use str_split_mod
   use neark, only: neark_dp_opti, neark_dp_orig
   implicit none
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Local parameters
   integer, parameter :: MYPROC=0
   character(len=64), parameter :: MYNAME='test_integrals'
   integer, parameter :: NITER0=20000
   integer, parameter :: LONG_CHAR=2048
   integer, parameter :: STDOUT=6
   integer, parameter :: STDERR=0
   integer, parameter :: NI0=100
   integer, parameter :: NK0=80
   real, parameter :: PI=3.14159
   
   ! Local variable declaration
   integer, parameter :: NCLOCK = 64
   real, save :: tmg_elapsed(NCLOCK) = 0.
   character(len=10) :: tmg_name(NCLOCK)
   integer :: tmg_time0(NCLOCK)

   integer :: ni,nk,niter

   !-------------------------------------------------------------------
   tmg_name = ''

   ni = NI0
   nk = NK0
   niter = NITER0
   call get_input_args2()

   call test_neark(MYNAME,ni,nk,niter)
   
   call ttmg_print()

   !-------------------------------------------------------------------
   return

contains
   
   subroutine test_neark(name,ni,nk,niter)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(in) :: ni,nk,niter
      
      real :: sig(ni,nk), ps(ni), p0,p1
      integer :: i,k,kidx1(ni), kidx0(ni), iter, k0,k1
      real, parameter :: DPREF = 201.
      !-------------------------------------------------------------------

      print *,'[BEGIN] test_neark'

      do i=1,ni
         ps(i) = 1000. - (float(i-(ni/2))/float(ni)) * 50.
      enddo
      do k=1,nk
!!$         sig(:,k) =  (1000. * (float(nk-k+1)/float(nk))) / ps(:)
         sig(:,k) =  float(k)/float(nk) ! (1000. * (float(k)/float(nk))) / ps(:)
      enddo

      do iter=1,niter
         call ttmg_start(2,"nrkorig")
         kidx0 = neark_dp_orig(sig,ps,DPREF,ni,nk)
         call ttmg_stop(2)

         call ttmg_start(1,"nrkopti")
         kidx1 = neark_dp_opti(sig,ps,DPREF,ni,nk)
         call ttmg_stop(1)
      enddo

      do i=1,ni
         if (kidx0(i) /= kidx1(i)) then
            p0 = ps(i)*(1.-sig(i,kidx0(i)))
            p1 = ps(i)*(1.-sig(i,kidx1(i)))
            write(STDOUT,'(i5,a," orig[",i5,"]",f12.6," (",f12.6,") : opt[",i5,"]",f12.6," (",f12.6,")")') &
                 i, ' ERROR diff:', &
                 kidx0(i), p0, p0 - DPREF, &
                 
                 kidx1(i), p1, p1 - DPREF

         endif
      enddo
      k0 = max(1,min(kidx0(1),kidx1(1))-1)
      k1 = min(max(kidx0(1),kidx1(1))+1,nk)
      i=1
      do k=k0,k1
         write(STDOUT,'("k=",i5," s=",f9.6," p=",f12.6," dp=",f12.6," (",f12.6,")")') &
              k, sig(i,k), ps(i)*sig(i,k), ps(i)*(1-sig(i,k)), DPREF-ps(i)*(1-sig(i,k))
      enddo

      print *,'[END] test_neark'
      !-------------------------------------------------------------------
      return
   end subroutine test_neark

   
   subroutine ttmg_start(i, name)
      integer, intent(in) :: i
      character(len=*), intent(in)  :: name
      call system_clock(count=tmg_time0(i))
      tmg_name(i) = name
   end subroutine ttmg_start

   
   subroutine ttmg_stop(i)
      integer, intent(in) :: i
      integer :: count_rate, count_end
      call system_clock(count_rate=count_rate)
      call system_clock(count=count_end)
      tmg_elapsed(i) = tmg_elapsed(i) + &
           real(count_end - tmg_time0(i)) / real(count_rate)
   end subroutine ttmg_stop

   
   subroutine ttmg_print()
      integer :: i
      print *,'==============================================='
      do i=1,size(tmg_elapsed)
         if (tmg_name(i) /= '') then
            print *,'Timing: ', tmg_name(i), tmg_elapsed(i)
         endif
      enddo
   end subroutine ttmg_print

   
   subroutine get_input_args2()
      integer :: iarg, nslabs, nsubsteps2, istat
      character(len=64) :: str, key, val
      logical :: ok_L
      real :: dtin
   !-------------------------------------------------------------------
      do iarg = 1, command_argument_count()
         call get_command_argument(iarg, str)
         if (str(1:2) == '-h') then
            print *,'Help: Can specify params with key=value pairs'
            print *,'  Known keys:'
            print *,'  '
            call flush(6)
            call exit(1)
         endif
         call str_normalize(str)
         call str_split(key, val, str, '=')
!!$         istat = clib_tolower(key)
         ok_L = .true.
         select case(key(1:5))
         case ('ni   ')
            istat = str_toint(ni,val)
         case ('nk   ')
            istat = str_toint(nk,val)
         case ('niter')
            istat = str_toint(niter,val)
         case default
            ok_L = .false.
         end select
         if (.not.ok_L) then
            print *, 'IGNORING UNKNOWN ARG: ',key(1:16),' : ',trim(str)
         else
            print *, 'INPUT ARG: ',key(1:16),' = ',trim(val)
         endif
      enddo

      print *, '### Input Args -------------------------------------------'
      print *, 'ni    = ', ni
      print *, 'nk    = ', nk
      print *, 'niter = ', niter
      print *,'----------------------------------------------------------'
      !-------------------------------------------------------------------
      return
   end subroutine get_input_args2
   
end subroutine test_phyutils
