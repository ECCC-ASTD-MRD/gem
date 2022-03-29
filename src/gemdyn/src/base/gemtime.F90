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

      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      implicit none
#include <time_trace.hf>

      public :: gemtime_init,gemtime_start,gemtime_stop,&
                gemtime_terminate,gemtime,Gem_timing_dyn_L,&
                gemlcltime_start,gemlcltime_stop

      external :: MPI_barrier

      type(time_context) :: gem_time_trace
      logical :: Gem_trace_barr

      integer, parameter :: MAX_instrumented=400
      integer, parameter :: MAX_event=500

      character(len=16) nam_subr_S(MAX_instrumented)
      character(len= 4) lcltime_name_S(MAX_instrumented)

      logical, save :: Gem_timing_dyn_L, New_timing_dyn_L, &
                       timing_barrier_L, lcltime_flag_L

      integer timer_cnt (MAX_instrumented), timer_level (MAX_instrumented)
      integer lcltime_id(MAX_instrumented), lcltime_rset(MAX_instrumented),&
              lcltime_step  (MAX_instrumented,MAX_event),&
              lcltime_nevent(MAX_instrumented),lcltime_last(MAX_instrumented)

      real(kind=REAL64) tb(MAX_instrumented), sum_tb(MAX_instrumented), total_time
      real(kind=REAL64) lcltime_tb(MAX_instrumented),lcltime_acc(MAX_instrumented)
      real(kind=REAL64) lcltime_event((MAX_instrumented),MAX_event)

contains

      subroutine gemtime_init ( myproc, msg )

      character(len=*), intent(in) :: msg
      integer       , intent(in) :: myproc
#include <clib_interface_mu.hf>

      character(len=16) timer_type_S,dumc_S
      integer err
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (clib_getenv ('GEMDYN_TIMING',timer_type_S) < 0) then
         Gem_timing_dyn_L = .false.
         New_timing_dyn_L = .true.
         timing_barrier_L = .false.
      else
         call low2up (timer_type_S,dumc_S)
         timer_type_S= dumc_S
         Gem_timing_dyn_L = timer_type_S(1:3) == 'DYN'
         New_timing_dyn_L = timer_type_S(1:3) == 'NEW'
         timing_barrier_L = timer_type_S(4:6) == '_WB'
      endif

      if (Gem_timing_dyn_L) then
         sum_tb= 0.; timer_cnt= 0 ; timer_level= 0 ; nam_subr_S= ''
         lcltime_id=0. ; lcltime_acc= 0. ; lcltime_name_S= ''
         lcltime_rset= -1 ; lcltime_last=0 ; lcltime_nevent=0
         lcltime_flag_L= .false.
         total_time= omp_get_wtime()
      else
         call timing_init2 ( myproc, msg )
      endif

      call rpn_comm_barrier ("GRID", err)
!MD      call time_trace_init (gem_time_trace)
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_init

      subroutine gemtime_start ( mynum, myname_S, mylevel )

      character(len=*), intent(in) :: myname_S
      integer         , intent(in) :: mynum,mylevel

      integer err
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (timing_barrier_L) call rpn_comm_barrier ("GRID", err)
      if (Gem_timing_dyn_L) then
         nam_subr_S (mynum) = myname_S
         timer_level(mynum) = mylevel
         tb         (mynum) = omp_get_wtime()
         timer_cnt  (mynum) = timer_cnt(mynum) + 1
      else if(New_timing_dyn_L) then
         call timing_start2 ( mynum, myname_S, mylevel )
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_start

      subroutine gemtime_stop (mynum)

      integer, intent(in) :: mynum

      integer err
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (timing_barrier_L) call rpn_comm_barrier ("GRID", err)
      if (Gem_timing_dyn_L) then
         sum_tb(mynum)= sum_tb (mynum) + (omp_get_wtime() - tb (mynum))
      else if(New_timing_dyn_L) then
         call timing_stop (mynum)
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemtime_stop

      subroutine gemtime_terminate ( myproc, msg )
      use ptopo
      use glb_ld
      use path
      use out_mod
      use out_options
      use gem_options
      use step_options

      character(len=*), intent(in) :: msg
      integer       , intent(in) :: myproc

#include <rmnlib_basics.hf>

      character(len=16) name
      character(len=64) fmt,nspace
      logical flag(MAX_instrumented)
      integer i,j,k,elem,lvl,lvlel(0:100),err,unf

      real(kind=REAL64) omp_get_wtime
      real(kind=REAL64) mytime(Ptopo_npex,Ptopo_npey),gmytime(Ptopo_npex,Ptopo_npey)
      real lcl_timing(Ptopo_npex,Ptopo_npey),w1(l_minx:l_maxx,l_miny:l_maxy),hyb0(1)
!
!-------------------------------------------------------------------
!
!MD      if ((Gem_trace_ctrl>0).and.(mod(Step_kount,Gem_trace_freq)>0)) &
!MD                                           call gemtime_trace_dump ()

      if (.not.New_timing_dyn_L .and. .not.Gem_timing_dyn_L) return

      if (Gem_timing_dyn_L) then

         if (lcltime_flag_L) then
            call out_open_file ('ti')
            call out_href ( 'Mass_point',1,G_ni,1,1,G_nj,1)
            if (Ptopo_myproc==0) then
               unf= 0
               err= fnom  ( unf, trim(Path_work_S)//'/lcl_timing.fst', 'STD+RND', 0 )
               err= fstouv( unf, 'RND' )
            endif
            do i = 1,MAX_instrumented
               if ( lcltime_nevent(i) > 0 ) then
                  do j=1,lcltime_nevent(i)
                     mytime=0. ; mytime(Ptopo_mycol+1,Ptopo_myrow+1)= lcltime_event(i,j)
                     call rpn_comm_REDUCE ( mytime, gmytime, Ptopo_npex*Ptopo_npey, &
                     "MPI_DOUBLE_PRECISION","MPI_SUM",0,"grid",err )
                     if (Ptopo_myproc==0) then
                        lcl_timing=gmytime
                        err = fstecr (lcl_timing,lcl_timing,-32,unf,Out_dateo,int(Out_deet),&
                                lcltime_step(i,j), Ptopo_npex,Ptopo_npey,1,0,0,&
                                lcltime_step(i,j),'P',lcltime_name_S(i),Out3_etik_S, &
                                'X',0,0,0,0,5,.false.)
                     endif
                     w1 = lcltime_event(i,j)
                     hyb0(1)=0.0 ; lvlel(0)=1
                     call out_fstecr(w1,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
                        lcltime_name_S(i),1.,0.,2,-1,1,lvlel(0), 1, 32,.false. )
                  end do
               endif
            end do
            if (Ptopo_myproc==0) err= fstfrm(unf)
            call out_cfile ()
         endif

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

!!$      subroutine gemtime_trace_dump
!!$      use ptopo
!!$      use step_options
!!$
!!$      character(len=1024) :: filename
!!$      type(C_PTR), dimension(10) :: array
!!$      integer :: nbeads, nbent, err
!!$      integer, dimension(:  ), pointer, contiguous :: timeT_ptr
!!$      integer, dimension(:,:), pointer, contiguous :: timeT_all
!!$!
!!$!-------------------------------------------------------------------
!!$!
!!$      nullify (timeT_ptr,timeT_all)
!!$      array(1) = C_NULL_PTR
!!$      array(1) = time_trace_get_buffer_data(gem_time_trace, nbeads, nbent, 1)
!!$      if(C_ASSOCIATED(array(1))) then
!!$         call c_f_pointer(array(1),timeT_ptr,[nbent])
!!$         allocate (timeT_all(nbent,Ptopo_numproc))
!!$         call RPN_COMM_gather(timeT_ptr,nbent,"MPI_INTEGER",timeT_all,nbent, &
!!$                              "MPI_INTEGER",0,"GRID", err)
!!$         if (Ptopo_myproc==0) then
!!$            write (filename,'("TIME_trace_",i5.5,".bin")') Step_kount
!!$            open (666,file=filename,form='UNFORMATTED')
!!$            write(666) nbent,Ptopo_npex,Ptopo_npey
!!$            write(666) timeT_all(1:nbent,1:Ptopo_numproc)
!!$            close (666)
!!$         endif
!!$         deallocate(timeT_all)
!!$         nullify (timeT_ptr,timeT_all)
!!$      endif
!!$      call time_trace_init (gem_time_trace)
!!$!
!!$!-------------------------------------------------------------------
!!$!
!!$      return
!!$      end subroutine gemtime_trace_dump

      subroutine gemlcltime_start ( mynum, myname_S, F_rset)

      character(len=*), intent(in) :: myname_S      ! name for the section
      integer         , intent(in) :: mynum, F_rset ! id and reset frequency

      integer err
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      call rpn_comm_barrier ("GRID", err)
      lcltime_name_S(mynum) = myname_S
      lcltime_tb    (mynum) = omp_get_wtime()
      lcltime_rset  (mynum) = F_rset
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemlcltime_start

      subroutine gemlcltime_stop (mynum,F_t0,F_tn)
      use step_options
      implicit none

      integer, intent(in) :: mynum     ! id
      integer, intent(in) :: F_t0,F_tn ! timestep bounds for printing

      integer knt
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      knt= Step_kount-1
      if ( (lcltime_rset(mynum)<1) .or. &
           (lcltime_last(mynum)/=Step_kount .and. mod(knt,lcltime_rset(mynum))==0) ) then
         lcltime_last(mynum)= Step_kount
         if (knt>=F_t0 .and. knt<=F_tn) then
            lcltime_nevent(mynum)= lcltime_nevent(mynum)+1
            lcltime_step (mynum,lcltime_nevent(mynum))= knt
            lcltime_event(mynum,lcltime_nevent(mynum))= lcltime_acc(mynum)
            lcltime_flag_L=.true.
         endif
         lcltime_acc(mynum)= 0.
      endif
      lcltime_acc(mynum)= lcltime_acc(mynum) + (omp_get_wtime() - lcltime_tb(mynum))
!
!-------------------------------------------------------------------
!
      return
      end subroutine gemlcltime_stop

      subroutine gemtime ( unf, from, last )
      use step_options, only : Step_maxwall
      use, intrinsic :: iso_fortran_env
      implicit none

      ! Unit file for output
      integer, intent(in) :: unf
      ! Does something special if the value is 'CURRENT TIMESTEP'
      character(len=*), intent(in) :: from
      ! If true, print the accumulated time
      logical, intent(in) :: last

      integer, external :: get_max_rss

      character date_S*8, time_S*10, jour_S*11, heure_S*10
      logical, save :: timini_L = .false.
      real          users, systs
      real, save :: user0 = 0.0, syst0 = 0.0
      real(kind=REAL64), save ::  START = -1.d0, END, avgtime(10) = 0.d0, &
                                 ACCUM_w = 0.d0, ACCUM_u = 0.d0, ACCUM_s = 0.d0
      real(kind=REAL64) omp_get_wtime

!----------------------------------------------------------------

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
                          systs-syst0, get_max_rss(), from
         ACCUM_w = ACCUM_w + END - START
         ACCUM_u = ACCUM_u + users - user0
         ACCUM_s = ACCUM_s + systs - syst0
         if (last) write(unf,1001) ACCUM_w,ACCUM_u,ACCUM_s
      else
         Step_maxwall = -1.d0
      endif

      user0 = users
      syst0 = systs
      START = END

      timini_L = .true.

!----------------------------------------------------------------

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
