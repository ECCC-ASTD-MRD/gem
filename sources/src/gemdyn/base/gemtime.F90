!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

module gem_timing

      use, intrinsic :: iso_fortran_env
      implicit none

      private
      public :: gemtime_init,gemtime_start,gemtime_stop,&
                gemtime_terminate,gemtime,Gem_timing_dyn_L

      integer, parameter :: MAX_instrumented=400

      character(len=16) nam_subr_S(MAX_instrumented)

      logical, save :: Gem_timing_dyn_L=.true.

      integer timer_cnt(MAX_instrumented), timer_level(MAX_instrumented)

      real(kind=REAL64) tb(MAX_instrumented), sum_tb(MAX_instrumented), total_time

contains

      subroutine gemtime_init ( myproc, msg )

      character(len=*), intent(in) :: msg
      integer       , intent(in) :: myproc
#include <clib_interface_mu.hf>

      character(len=16) dumc_S
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (clib_getenv ('GEMDYN_TIMING',dumc_S) < 0) Gem_timing_dyn_L=.false.
      if (Gem_timing_dyn_L) then
         sum_tb= 0; timer_cnt= 0 ; timer_level= 0 ; nam_subr_S= ''
         total_time= omp_get_wtime()
      else
         call timing_init2 ( myproc, msg )
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_init

      subroutine gemtime_start ( mynum, myname_S, mylevel )

      character(len=*), intent(in) :: myname_S
      integer       , intent(in) :: mynum,mylevel

      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (Gem_timing_dyn_L) then
         nam_subr_S (mynum) = myname_S
         timer_level(mynum) = mylevel
         tb         (mynum) = omp_get_wtime()
         timer_cnt  (mynum) = timer_cnt(mynum) + 1
      else
         call timing_start2 ( mynum, myname_S, mylevel )
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_start

      subroutine gemtime_stop (mynum)

      integer, intent(in) :: mynum

      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (Gem_timing_dyn_L) then
         sum_tb(mynum)= sum_tb (mynum) + (omp_get_wtime()    - tb (mynum))
      else
         call timing_stop (mynum)
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_stop

      subroutine gemtime_terminate ( myproc, msg )

      character(len=*), intent(in) :: msg
      integer       , intent(in) :: myproc

      character(len=16) name
      character(len=64) fmt,nspace
      logical flag(MAX_instrumented)
      integer i,j,k,elem,lvl,lvlel(0:100)

      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (Gem_timing_dyn_L) then

      if (myproc /= 0) return

      print *,'___________________________________________________________'
      print *,'__________________TIMINGS ON PE #0_________________________'

      flag=.false.
      total_time= omp_get_wtime() - total_time

      do i = 1,MAX_instrumented
         lvl= 0 ; elem= i
 55      if ( (trim(nam_subr_S(elem)) /= '') .and. (.not.flag(elem)) ) then

            write (nspace,'(i3)') 5*lvl+1
            fmt='(f6.2,"%",'//trim(nspace)//'x,a,1x,a,i3,a,3x,a,1pe13.6,2x,a,i8)'
            do k=1,len(name)
               name(k:k) = '.'
            end do
            name (len(name)-len(trim(nam_subr_S(elem)))+1:len(name))= &
            trim(nam_subr_S(elem))
            write (output_unit,trim(fmt)) sum_tb(elem)/total_time*100.,name,'(',elem,')','Wall clock= ',&
                   sum_tb(elem),'count= ',timer_cnt(elem)
            flag(elem) = .true. ; lvlel(lvl) = elem
 65         do j = 1,MAX_instrumented
               if ((timer_level(j) == elem) .and. (.not.flag(j)) )then
                  lvl= lvl+1
                  elem= j
                  goto 55
               endif
            end do
            lvl= lvl - 1
            if (lvl >= 0) then
               elem= lvlel(lvl)
               goto 65
            endif
         endif
      end do

      print *,'___________________________________________________________'

      else
         call timing_terminate2 ( myproc, msg )
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_terminate

      subroutine gemtime ( unf, from, last )
      use step_options

      character(len=*), intent(in) :: from
      logical, intent(in) :: last
      integer, intent(in) :: unf

      integer, external :: get_max_rss

      character date_S*8 ,time_S*10, jour_S*11, heure_S*10
      logical, save :: timini_L = .false.
      real          users,systs
      real, save :: user0=0.0, syst0=0.0
      real(kind=REAL64), save ::  START=-1.d0, END, avgtime(10)=0.d0, &
                                 ACCUM_w=0.d0, ACCUM_u=0.d0, ACCUM_s=0.d0
      real(kind=REAL64) omp_get_wtime
!
!----------------------------------------------------------------
!
      if (unf <= 0) return

      if (START < 0.d0) START = omp_get_wtime()

      call date_and_time( date_S,time_S )
      jour_S  = date_S(1:4)//'/'//date_S(5:6)//'/'//date_S(7:8)//' '
      heure_S = time_S(1:2)//'h'//time_S(3:4)//'m'//time_S(5:6)//'s,'

#if defined (AIX)
      call setrteopts('cpu_time_type=total_usertime')
#endif
      call cpu_time( users )

#if defined (AIX)
      call setrteopts('cpu_time_type=total_systime')
#endif
      call cpu_time( systs )

      END = omp_get_wtime()

      if (timini_L) then
         if (trim(from) =='CURRENT TIMESTEP') then
            avgtime(1:9) = avgtime(2:10)
            avgtime(10)  = END-START
            Step_maxwall = sum(avgtime) / 10.d0
         endif
         write(unf,1000) 'TIME: '//jour_S//heure_S, END-START, users-user0, &
                          systs-syst0,get_max_rss(),from
         ACCUM_w= ACCUM_w + END-START
         ACCUM_u= ACCUM_u + users-user0
         ACCUM_s= ACCUM_s + systs-syst0
         if (last) write(unf,1001) ACCUM_w,ACCUM_u,ACCUM_s
      else
         Step_maxwall = -1.d0
      endif

      user0= users  ;  syst0= systs  ;  START= END

      timini_L = .true.
!
!----------------------------------------------------------------
!
 1000    format(/A,' W: ',1pe13.6, &
                   ' U: ',1pe13.6, &
                   ' S: ',1pe13.6, &
                  ', Mem: ',i7,' (Kbytes/PE) ',a)
 1001    format(/'ACCUMULATED TIME: W: ',1pe13.6, &
                   ' U: ',1pe13.6, &
                   ' S: ',1pe13.6)

      return
      end subroutine gemtime

end module gem_timing
