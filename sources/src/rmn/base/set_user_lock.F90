      subroutine set_user_lock(owner_thread,lock)
! manage user locks
!
! lock = .true. attempt to acquire lock, .false. release lock
! owner_thread : static integer variable used as lock, MUST be initialized to a negative value initially
!
! call set_user_lock(lock_variable,.true. ) to acquire lock
! call set_user_lock(lock_variable,.false.) to release lock
!
      IMPLICIT NONE
      logical, intent(IN) :: lock
      integer, intent(INOUT) :: owner_thread
#if defined(_OPENMP)
      integer, external :: omp_get_thread_num
      integer :: current_pid, owner_pid, owner_count
      logical :: ok

!     if (owner_thread < 0) owner_thread = 0
      current_pid = omp_get_thread_num()
      current_pid = iand(current_pid,x"FFFFFF")   ! keep lower 24 bits (8-31)

!     print *,'thread=',current_pid,' lock=',lock,' owner_pid=',owner_pid,' count=',owner_count
      ok = .false.
      do while(.not. ok)
!$OMP CRITICAL
      owner_pid   = iand(owner_thread,x"FFFFFF")  ! lower 24 bits (8-31)
      owner_count = ishft(owner_thread,-24)       ! bits 1-7
      owner_count = iand(owner_count,127)
      if(lock) then
!       print *,'attempt to lock by thread',current_pid
        if(owner_count == 0 .or. owner_pid == current_pid) then      ! lock is not owned or already owned by this thread, acquire it
          owner_count = owner_count + 1                               ! add one to depth
          owner_thread = current_pid + ishft(owner_count,24)
          ok = .true.
        endif
!       if(.not. ok) print *,'unsuccessful attempt to lock by thread',current_pid
      else
!       print *,'attempt to unlock by thread',current_pid
        if(current_pid .ne. owner_pid .and. owner_thread > 0) then  ! lock is owned by another thread
          print *,'ERROR: ',current_pid,' attempting to release a lock owned by',owner_pid
          call qqexit(1)
        endif
        owner_count = owner_count - 1                              ! subtract one from depth
        owner_thread = current_pid + ishft(owner_count,24)
        if(owner_count <= 0) owner_thread = 0   ! this thread's depth count for the lock is zero, release it totally
        ok = .true.
      endif
!$OMP END CRITICAL
      enddo
!     if(lock)      print *,'locked by thread',current_pid,' count=',owner_count
!     if(.not. lock)print *,'unlocked by thread',current_pid,' count=',owner_count
#endif
      return
      end subroutine set_user_lock
