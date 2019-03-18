
subroutine tmg_init(myproc, msg)
   integer :: myproc
   character(len=*) :: msg
   return
end subroutine tmg_init

subroutine tmg_start(mynum, myname_S)
   integer :: mynum
   character(len=*) :: myname_S
   return
end subroutine tmg_start

subroutine tmg_stop(mynum)
   integer :: mynum
   return
end subroutine tmg_stop

subroutine tmg_terminate(myproc, msg)
   integer :: myproc
   character(len=*) :: msg
   return
end subroutine tmg_terminate
