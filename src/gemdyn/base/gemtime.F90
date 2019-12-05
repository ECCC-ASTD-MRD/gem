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

      type(time_context) :: gem_time_trace

      integer, parameter :: MAX_instrumented = 400
      integer, parameter :: MAX_event = 500

      logical, save :: Gem_timing_dyn_L, New_timing_dyn_L, &
                       timing_barrier_L, lcltime_flag_L

      character(len = 16), dimension(MAX_instrumented) :: nam_subr_S
      character(len = 4), dimension(MAX_instrumented) ::  lcltime_name_S

      integer, dimension(MAX_instrumented) :: timer_cnt
      integer, dimension(MAX_instrumented) :: timer_level
      integer, dimension(MAX_instrumented) :: lcltime_id
      integer, dimension(MAX_instrumented) :: lcltime_rset
      integer, dimension(MAX_instrumented) :: lcltime_nevent
      integer, dimension(MAX_instrumented) :: lcltime_last

      integer, dimension(MAX_instrumented, MAX_event) :: lcltime_step

      real(kind=REAL64) :: total_time
      real(kind=REAL64), dimension(MAX_instrumented) :: tb
      real(kind=REAL64), dimension(MAX_instrumented) :: sum_tb
      real(kind=REAL64), dimension(MAX_instrumented) :: lcltime_tb
      real(kind=REAL64), dimension(MAX_instrumented) :: lcltime_acc

      real(kind=REAL64), dimension(MAX_instrumented, MAX_event) :: lcltime_event

contains

   subroutine gemtime_init ( myproc, msg )

      character(len=*), intent(in) :: msg
      integer, intent(in) :: myproc
#include <clib_interface_mu.hf>

      character(len=16) timer_type_S,dumc_S
      real(kind=REAL64) omp_get_wtime

      !-------------------------------------------------------------------

      if (clib_getenv ('GEMDYN_TIMING', timer_type_S) < 0) then
         Gem_timing_dyn_L = .false.
         New_timing_dyn_L = .true.
         timing_barrier_L = .false.
      else
         call low2up (timer_type_S, dumc_S)
         timer_type_S = dumc_S
         Gem_timing_dyn_L = timer_type_S(1:3) == 'DYN'
         New_timing_dyn_L = timer_type_S(1:3) == 'NEW'
         timing_barrier_L = timer_type_S(4:6) == '_WB'
      endif

      if (Gem_timing_dyn_L) then
         sum_tb = 0d0
         timer_cnt = 0
         timer_level = 0
         nam_subr_S= ''
         lcltime_id = 0
         lcltime_acc = 0d0
         lcltime_name_S = ''
         lcltime_rset = -1
         lcltime_last = 0
         lcltime_nevent = 0
         lcltime_flag_L = .false.
         total_time = omp_get_wtime()
         call time_trace_init(gem_time_trace)
      else
         call timing_init2(myproc, msg)
      endif

   end subroutine gemtime_init


   subroutine gemtime_start ( id, name, level )

      character(len=*), intent(in) :: name
      integer, intent(in) :: id, level

      integer :: err
      real(kind=REAL64) :: omp_get_wtime

      if (timing_barrier_L) call rpn_comm_barrier ("GRID", err)

      if (Gem_timing_dyn_L) then
         nam_subr_S (id) = name
         timer_level(id) = level
         tb         (id) = omp_get_wtime()
         timer_cnt  (id) = timer_cnt(id) + 1
      else if(New_timing_dyn_L) then
         call timing_start2 ( id, name, level )
      endif

   end subroutine gemtime_start


   subroutine gemtime_stop (id)

      integer, intent(in) :: id

      integer :: err
      real(kind=REAL64) :: omp_get_wtime

      !-------------------------------------------------------------------

      if (timing_barrier_L) call rpn_comm_barrier ("GRID", err)

      if (Gem_timing_dyn_L) then
         sum_tb(id) = sum_tb(id) + (omp_get_wtime() - tb(id))
      else if(New_timing_dyn_L) then
         call timing_stop (id)
      endif

   end subroutine gemtime_stop

   subroutine gemtime_terminate ( myproc, msg )
      use ptopo
      use glb_ld
      use path
      use out_mod
      use out_options

      character(len=*), intent(in) :: msg
      integer, intent(in) :: myproc

#include <rmnlib_basics.hf>

      character(len = 16) :: name
      character(len = 64) :: fmt, nspace
      logical, dimension(MAX_instrumented) :: flag
      integer :: i, j, k, elem, lvl, err, unf
      integer, dimension(0:100) :: lvlel

      real(kind=REAL64) :: omp_get_wtime
      real(kind=REAL64) :: mytime(Ptopo_npex, Ptopo_npey), gmytime(Ptopo_npex, Ptopo_npey)
      real :: lcl_timing(Ptopo_npex, Ptopo_npey)
      real :: w1(l_minx:l_maxx, l_miny:l_maxy), hyb0(1)

!-------------------------------------------------------------------

      if (.not.New_timing_dyn_L .and. .not.Gem_timing_dyn_L) return

      if (Gem_timing_dyn_L) then

         if (lcltime_flag_L) then
            call out_open_file ('ti')
            call out_href('Mass_point', 1, G_ni, 1, 1, G_nj, 1)
            if (Ptopo_myproc == 0) then
               unf = 0
               err = fnom(unf, trim(Path_basedir_S) // '/lcl_timing.fst', 'STD+RND', 0 )
               err = fstouv(unf, 'RND')
            endif
            do i = 1, MAX_instrumented
               if ( lcltime_nevent(i) > 0 ) then
                  do j = 1, lcltime_nevent(i)
                     mytime = 0d0
                     mytime(Ptopo_mycol + 1, Ptopo_myrow + 1) = lcltime_event(i,j)
                     call rpn_comm_REDUCE ( mytime, gmytime, Ptopo_npex * Ptopo_npey, &
                        "MPI_DOUBLE_PRECISION", "MPI_SUM", 0, "grid", err )
                     if (Ptopo_myproc == 0) then
                        lcl_timing = gmytime
                        err = fstecr (lcl_timing, lcl_timing, -32, unf, Out_dateo, &
                           int(Out_deet), lcltime_step(i, j), Ptopo_npex, Ptopo_npey, &
                           1, 0, 0, lcltime_step(i, j), 'P', lcltime_name_S(i), &
                           Out3_etik_S, 'X', 0, 0, 0, 0, 5, .false.)
                     endif
                     w1 = lcltime_event(i, j)
                     hyb0(1) = 0.0
                     lvlel(0) = 1
                     call out_fstecr(w1, l_minx, l_maxx, l_miny, l_maxy, hyb0, &
                        lcltime_name_S(i), 1.0, 0.0, 2, -1, 1, lvlel(0), 1, 32, .false. )
                  end do
               endif
            end do
            if (Ptopo_myproc == 0) err = fstfrm(unf)
            call out_cfile
         endif

         call time_trace_dump_text(gem_time_trace, 'time_list', Ptopo_myproc)

         if (myproc /= 0) return

         print *,'___________________________________________________________'
         print *,'__________________TIMINGS ON PE #0_________________________'

         flag = .false.
         total_time = omp_get_wtime() - total_time

         do i = 1, MAX_instrumented
            lvl = 0
            elem = i
 55      if ( (trim(nam_subr_S(elem)) /= '') .and. (.not.flag(elem)) ) then


               write (nspace, '(i3)') 5 * lvl + 1
               fmt = '(f6.2,"%",' // trim(nspace) // 'x,a,1x,a,i3,a,3x,a,1pe13.6,2x,a,i8)'

               ! Fill name with dots
               do k = 1, len(name)
                  name(k:k) = '.'
               end do
               ! Right-justify the name
               name (len(name) - len(trim(nam_subr_S(elem))) + 1:len(name)) = &
                  trim(nam_subr_S(elem))

               write (output_unit, trim(fmt)) &
                  sum_tb(elem) / total_time * 100.0, name, '(', elem, ')', &
                  'Wall clock= ', sum_tb(elem), 'count= ', timer_cnt(elem)

               flag(elem) = .true.
               lvlel(lvl) = elem
 65         do j = 1, MAX_instrumented
                  if ((timer_level(j) == elem) .and. (.not.flag(j)) ) then
                     lvl = lvl + 1
                     elem = j
                     goto 55
                  endif
               end do
               lvl = lvl - 1
               if (lvl >= 0) then
                  elem = lvlel(lvl)
                  goto 65
               endif
            endif
         end do

         print *,'___________________________________________________________'

      else
         call timing_terminate2 ( myproc, msg )
      endif

   end subroutine gemtime_terminate


   subroutine gemlcltime_start ( id, name, F_rset)

      character(len=*), intent(in) :: name      ! name for the section
      integer         , intent(in) :: id, F_rset ! id and reset frequency

      integer :: err
      real(kind=REAL64) :: omp_get_wtime

      !-------------------------------------------------------------------

      call rpn_comm_barrier ("GRID", err)
      lcltime_name_S(id) = name
      lcltime_tb    (id) = omp_get_wtime()
      lcltime_rset  (id) = F_rset

   end subroutine gemlcltime_start


   subroutine gemlcltime_stop (id, F_t0, F_tn)
      use step_options
      implicit none

      integer, intent(in) :: id     ! id
      integer, intent(in) :: F_t0, F_tn ! timestep bounds for printing

      integer knt
      real(kind=REAL64) omp_get_wtime

      !-------------------------------------------------------------------

      knt= Step_kount-1
      if ( (lcltime_rset(id) < 1) .or. &
           (lcltime_last(id) /= Step_kount .and. mod(knt,lcltime_rset(id))==0) ) then
         lcltime_last(id)= Step_kount
         if (knt >= F_t0 .and. knt <= F_tn) then
            lcltime_nevent(id) = lcltime_nevent(id) + 1
            lcltime_step (id, lcltime_nevent(id)) = knt
            lcltime_event(id, lcltime_nevent(id)) = lcltime_acc(id)
            lcltime_flag_L = .true.
         endif
         lcltime_acc(id) = 0d0
      endif
      lcltime_acc(id)= lcltime_acc(id) + (omp_get_wtime() - lcltime_tb(id))

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

      character(len = 8) :: date_S
      character(len = 10) :: time_S
      character(len = 11) :: jour_S
      character(len = 10) :: heure_S
      logical, save :: timini_L = .false.
      real :: users
      real :: systs
      real, save :: user0 = 0.0, syst0 = 0.0
      real(kind=REAL64), save ::  START = -1.d0, END, avgtime(10) = 0.d0, &
                                 ACCUM_w = 0.d0, ACCUM_u = 0.d0, ACCUM_s = 0.d0
      real(kind=REAL64) :: omp_get_wtime

      !----------------------------------------------------------------

      if (unf <= 0) return

      if (START < 0.d0) START = omp_get_wtime()

      call date_and_time( date_S, time_S )
      jour_S  = date_S(1:4) // '/' // date_S(5:6) // '/' // date_S(7:8) // ' '
      heure_S = time_S(1:2) // 'h' // time_S(3:4) // 'm' // time_S(5:6) // 's,'

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
            avgtime(10)  = END - START
            Step_maxwall = sum(avgtime) / 10.d0
         endif
         write(unf, 1000) 'TIME: ' // jour_S // heure_S, END-START, users-user0, &
                          systs-syst0, get_max_rss(), from
         ACCUM_w = ACCUM_w + END - START
         ACCUM_u = ACCUM_u + users - user0
         ACCUM_s = ACCUM_s + systs - syst0
         if (last) write(unf, 1001) ACCUM_w, ACCUM_u, ACCUM_s
      else
         Step_maxwall = -1.d0
      endif

      user0 = users
      syst0 = systs
      START = END

      timini_L = .true.

 1000    format(/A,' W: ', 1pe13.6, &
                   ' U: ', 1pe13.6, &
                   ' S: ', 1pe13.6, &
                  ', Mem: ', i7, ' (Kbytes/PE) ', a)
 1001    format(/'ACCUMULATED TIME:' &
                   ' W: ', 1pe13.6, &
                   ' U: ', 1pe13.6, &
                   ' S: ', 1pe13.6)

   end subroutine gemtime

end module gem_timing
