subroutine test_integrals
   use integrals
   ! Test integral equation solver
   implicit none
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Local parameters
   integer, parameter :: NI=1
   integer, parameter :: LONG_CHAR=2048
   integer, parameter :: STDOUT=6
   real, parameter :: PI=3.14159
   
   ! Local variable declaration
   integer :: nk,k,fd_sol,fd_prof,fd_const,narg,dir,istat
   real :: alpha,beta,l,dz,czdep,ca,zexp
   real, dimension(NI) :: zdep,a,zval
   real, dimension(:,:), allocatable :: z,df
   character(len=LONG_CHAR) :: const_file,prof_file,sol_file,cdir

   ! Retrieve command line inputs
   narg = command_argument_count()
   if (narg < 3) then
      call msg(MSG_ERROR,'(checkint) Expected three file names as arguments')
      stop
   endif
   call get_command_argument(1,const_file)
   call get_command_argument(2,prof_file)
   call get_command_argument(3,sol_file)

   ! Read in constants
   fd_const = 0
   if (fnom(fd_const,const_file,'SEQ+OLD',0) /= RMN_OK) then
      call msg(MSG_ERROR,'(checkint) Cannot open constants file')
      stop
   endif
   read(fd_const,*)
   read(fd_const,*) alpha,beta,l,nk,dz,czdep,ca,dir,zexp
   if (fclos(fd_const) /= RMN_OK) &
        call msg(MSG_ERROR,'(checkint) Cannot close constants file')
   zdep(:) = czdep
   a(:) = ca
   cdir = 'UP'
   if (dir < 0) cdir = 'DOWN'

   ! Set up vertical profile
   allocate(z(NI,nk),df(NI,nk))
   do k=1,NK
      z(:,k) = (NK-k+1)*dz*exp(zexp*(NK-k))
      df(:,k) = alpha*sin(2*PI/l*z(:,k)) + beta*z(:,k)
   enddo

   ! Solve integral equation
   fd_prof = 0
   if (fnom(fd_prof,prof_file,'SEQ',0) /= RMN_OK) &
        call msg(MSG_WARNING,'(checkint) Cannot open profile output file')
   istat = int_solve(zval,dir*df,z,zdep,a,cdir,fd_unittest=fd_prof)   
   if (istat == INT_ERR) call msg(MSG_ERROR,'(checkint) Error returned by int_solve')
   write(STDOUT, '(i0)') istat
   if (fclos(fd_prof) /= RMN_OK) &
        call msg(MSG_WARNING,'(checkint) Cannot close profile output file')

   ! Generate result for plotting
   fd_sol = 0
   if (fnom(fd_sol,sol_file,'SEQ',0) /= RMN_OK) &
        call msg(MSG_WARNING,'(checkint) Cannot open solution output file')
   write(fd_sol, '(f0.3)') zdep + dir*zval
   if (fclos(fd_sol) /= RMN_OK) &
        call msg(MSG_WARNING,'(checkint) Cannot close solution output file')

   ! Shutdown
   deallocate(z,df)

end subroutine test_integrals
