subroutine test_integrals
   use integrals
   ! Test integral equation solver
   implicit none
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Local parameters
   integer, parameter :: MYPROC=0
   character(len=64), parameter :: MYNAME='test_integrals'
   integer, parameter :: NITER=20000
   integer, parameter :: LONG_CHAR=2048
   integer, parameter :: STDOUT=6
   integer, parameter :: STDERR=0
   integer, parameter :: NI0=250
   integer, parameter :: NK0=80
   real, parameter :: PI=3.14159
   
   ! Local variable declaration
   integer :: i,nk,k,fd_sol,fd_prof,fd_const,narg,dir,istat,n,cdir
   real :: alpha,beta,l,dz,czdep,ca,zexp,zmid,fac
   real, dimension(:,:), allocatable :: z,df,xzval,xzval2
   character(len=LONG_CHAR) :: const_file,prof_file,sol_file

!!$   ! Retrieve command line inputs
!!$   narg = command_argument_count()
!!$   if (narg < 3) then
!!$      call msg(MSG_ERROR,'(checkint) Expected three file names as arguments')
!!$      stop
!!$   endif
!!$   call get_command_argument(1,const_file)
!!$   call get_command_argument(2,prof_file)
!!$   call get_command_argument(3,sol_file)
!!$
!!$   !call timing_init2(MYPROC, MYNAME)
!!$   !call timing_start2(1," test_integrals", 0 )
!!$
!!$   ! Read in constants
!!$   fd_const = 0
!!$   if (fnom(fd_const,const_file,'SEQ+OLD',0) /= RMN_OK) then
!!$      call msg(MSG_ERROR,'(checkint) Cannot open constants file')
!!$      stop
!!$   endif
!!$   read(fd_const,*)
!!$   read(fd_const,*) alpha,beta,l,nk,dz,czdep,ca,dir,zexp
!!$   if (fclos(fd_const) /= RMN_OK) &
!!$        call msg(MSG_ERROR,'(checkint) Cannot close constants file')

!!$   call run_loop(const_file, alpha, beta, l, ni, nk, dz, czdep, ca, dir, zexp)

   !#                     alpha,  beta,     l,  ni,  nk,    dz,   zdep,     a, dir, zexp
   call run_loop('default', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  100.0,  250.0,  1, 0.00)
   call run_loop('dmx    ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  520.0,  100.0, -1, 0.00)
   call run_loop('umxf   ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  100.0,   50.0,  1, 0.00)
   call run_loop('umnf   ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  210.0,  -50.0,  1, 0.00)
   call run_loop('uincf  ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  250.0,  100.0,  1, 0.00)
   call run_loop('udecf  ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  150.0, -120.0,  1, 0.00)
   call run_loop('umx    ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  180.0,  400.0,  1, 0.00)
   call run_loop('umn    ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,   80.0, -125.0,  1, 0.00)
   call run_loop('uinc   ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  100.0,  300.0,  1, 0.00)
   call run_loop('udec   ', 5.0, 0.002, 450.0, NI0,  NK0, 100.0,  280.0, -300.0,  1, 0.00)
   call run_loop('dmxf   ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  395.0,   20.0, -1, 0.00)
   call run_loop('dmnf   ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  780.0, -100.0, -1, 0.00)
   call run_loop('dincf  ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  480.0,   75.0, -1, 0.00)
   call run_loop('ddecf  ', 5.0, 0.002, 400.0, NI0,  NK0, 100.0,  580.0, -200.0, -1, 0.00)
   call run_loop('dmx    ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  520.0,  100.0, -1, 0.00)
   call run_loop('dmn    ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  460.0, -200.0, -1, 0.00)
   call run_loop('dinc   ', 5.0, 0.002, 500.0, NI0,  NK0, 100.0,  460.0,  400.0, -1, 0.00)
   call run_loop('ddec   ', 5.0, 0.002, 250.0, NI0,  NK0, 100.0,  740.0, -200.0, -1, 0.00)
   call run_loop('vuinc  ', 5.0, 0.002, 600.0, NI0,  NK0, 100.0,  400.0,  250.0,  1, 0.05)
   call run_loop('vdmnf  ', 5.0, 0.002, 550.0, NI0,  NK0, 100.0, 1120.0, -150.0, -1, 0.05)

  
!!$   !call timing_stop(1)
!!$   !call timing_terminate2(MYPROC, MYNAME)

   return

contains

   subroutine run_loop(name,alpha,beta,l,ni,nk,dz,czdep,ca,dir,zexp)
      character(len=*) :: name
      real :: alpha,beta,l,dz,czdep,ca,zexp
      integer :: ni,nk,dir

      real :: zmid,fac,dz1,dz2,dzw,w
      real, dimension(ni) :: zdep,a,zval,z1,z2

      zdep(:) = czdep
      a(:) = ca
      cdir = INT_DIR_UP
      if (dir < 0) cdir = INT_DIR_DOWN
      
      dz1 = dz*1.05
      dz2 = dz*0.95
      
      ! Set up vertical profile
      allocate(z(ni,nk),df(ni,nk),xzval(ni,nk),xzval2(ni,nk))
      do k=1,NK
         do i=1,ni
            w = float(i)/float(ni)
            dzw = w*dz1 + (1-w)*dz2
            z(i,k) = (NK-k+1)*dzw*exp(zexp*(NK-k))
            df(i,k) = dir*alpha*sin(2*PI/l*z(i,k)) + beta*z(i,k)
            xzval(i,k) = 0.
            xzval2(i,k) = 0.
         enddo
      enddo
      do i=1,ni
         zmid = (z(i,1) + z(i,nk))/2.
         fac = float((i-1))/float(ni)
         z1(i) = z(i,1) + fac*zmid
         z2(i) = z(i,nk) - fac*zmid
         zval(i) = 0.
      enddo

      ! Solve integral equation
      fd_prof = 0
!!$   if (fnom(fd_prof,prof_file,'SEQ',0) /= RMN_OK) &
!!$        call msg(MSG_WARNING,'(checkint) Cannot open profile output file')

      call loop_integrals(zval,xzval,xzval2,zdep,a,z1,z2,ni,nk)

      if (istat == INT_ERR) call msg(MSG_ERROR,'(checkint) Error returned by int_solve')
      write(STDOUT, '(i0)') istat
!!$   if (fclos(fd_prof) /= RMN_OK) &
!!$        call msg(MSG_WARNING,'(checkint) Cannot close profile output file')

      write(STDERR,*) ''
      write(STDERR,*) trim(name)//' int_solve:',zval(1),zval((ni+1)/2)
      write(STDERR,*) trim(name)//' int_profile pc:',xzval(1,:),xzval((ni+1)/2,:)
      write(STDERR,*) trim(name)//' int_profile li:',xzval2(1,:),xzval2((ni+1)/2,:)

!!$   ! Generate result for plotting
!!$   fd_sol = 0
!!$   if (fnom(fd_sol,sol_file,'SEQ',0) /= RMN_OK) &
!!$        call msg(MSG_WARNING,'(checkint) Cannot open solution output file')
!!$   write(fd_sol, '(f0.3)') zdep(1) + dir*zval(1)
!!$   if (fclos(fd_sol) /= RMN_OK) &
!!$        call msg(MSG_WARNING,'(checkint) Cannot close solution output file')

      ! Shutdown
      deallocate(z,df,xzval,xzval2)
      return
   end subroutine run_loop
   
   subroutine loop_integrals(zval,xzval,xzval2,zdep,a,z1,z2,ni,nk)
      integer, intent(in) :: ni,nk
      real, dimension(ni), intent(in) :: zdep,a,z1,z2
      real, dimension(ni), intent(out) :: zval
      real, dimension(ni,nk), intent(out) :: xzval,xzval2

      integer :: nt
      
      do nt=1,NITER
!!$      call timing_start2(2, "int_solve", 1)
         zval = 0.
!!$         istat = int_solve(zval,df,z,zdep,a,cdir,fd_unittest=fd_prof)
         istat = int_solve(zval,df,z,zdep,a,cdir,ni,nk)
!!$      call timing_stop(2)
!!$      call timing_start2(3, "int_profile", 1)
         istat = int_profile(xzval,df,z,z1,z2)
         istat = int_profile(xzval2,df,z,z1,z2,INT_TYPE_LINEAR)
!!$      call timing_stop(3)
      enddo
      return
   end subroutine loop_integrals
   
end subroutine test_integrals
