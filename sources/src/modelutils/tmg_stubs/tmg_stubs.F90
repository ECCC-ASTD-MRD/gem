
subroutine tmg_init(myproc, msg)
   integer :: myproc
   character(len=*) :: msg
   if (myproc == 0) print *, 'tmg_init (stub):',msg
   return
end subroutine tmg_init

subroutine tmg_start(mynum, myname_S)
   integer :: mynum
   character(len=*) :: myname_S
   if (mynum == 0) print *, 'tmg_start (stub):',myname_S, mynum
   return
end subroutine tmg_start

subroutine tmg_stop(mynum)
   integer :: mynum
   if (mynum == 0) print *, 'tmg_stop (stub)', mynum
   return
end subroutine tmg_stop

subroutine tmg_terminate(myproc, msg)
   integer :: myproc
   character(len=*) :: msg
   if (myproc == 0) print *, 'tmg_terminate (stub):',msg
   return
end subroutine tmg_terminate
