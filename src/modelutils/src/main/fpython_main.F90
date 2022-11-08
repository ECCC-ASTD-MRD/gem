subroutine fpython_main()
   use ISO_C_BINDING
   implicit none

   interface
      subroutine main(argc, argv) BIND(C,name='cpython_interp')
         use ISO_C_BINDING
         integer(C_INT), intent(IN), value :: argc
         type(C_PTR), intent(IN), dimension(0:argc) :: argv
      end subroutine main
   end interface

   external :: ftn_iargc, ftn_getarg
   type(C_PTR), pointer, dimension(:) :: argv
   integer :: argc, i, len, stat, ip, j
   character(len=4096) :: arg
   character(kind=C_CHAR), dimension(32768), target :: args

   argc = command_argument_count()
   allocate(argv(0:argc))
   ip = 0
   do i=0,argc
      call get_command_argument(i, arg, len, stat)
      argv(i) = c_loc(args(ip+1))
      do j=1,len
         args(ip+j) = arg(j:j)
      enddo
      ip = ip + len + 1
      args(ip) = char(0)
   enddo
   call main(argc+1,argv)
   return
end subroutine fpython_main
