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

module omp_timing

      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      implicit none

      public

      integer, parameter :: MAX_instrumented=400
      integer, parameter :: MAX_event=500

      character(len=16) nam_subr_S(MAX_instrumented)

      logical, save :: Gem_timing_dyn_L, New_timing_dyn_L, &
                       timing_barrier_L

      integer timer_cnt (MAX_instrumented), timer_level (MAX_instrumented)

      real(kind=REAL64) tb(MAX_instrumented), sum_tb(MAX_instrumented), total_time

contains

      subroutine gtmg_init ( myproc, msg )
      implicit none
      character(len=*), intent(in) :: msg
      integer         , intent(in) :: myproc
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
         total_time= omp_get_wtime()
      else
         call timing_init2 ( myproc, msg )
      endif

      call rpn_comm_barrier ("GRID", err)
!
!-------------------------------------------------------------------
!
      return
      end subroutine gtmg_init

      subroutine gtmg_start ( mynum, myname_S, mylevel )
      implicit none

      character(len=*), intent(in) :: myname_S
      integer         , intent(in) :: mynum,mylevel

      integer err
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (timing_barrier_L) then
!$OMP BARRIER
      endif
!$omp single
      if (timing_barrier_L) call rpn_comm_barrier ("GRID", err)
      if (Gem_timing_dyn_L) then
         nam_subr_S (mynum) = myname_S
         timer_level(mynum) = mylevel
         tb         (mynum) = omp_get_wtime()
         timer_cnt  (mynum) = timer_cnt(mynum) + 1
      else if(New_timing_dyn_L) then
         call timing_start2 ( mynum, myname_S, mylevel )
      endif
!$omp end single
!
!-------------------------------------------------------------------
!
      return
      end subroutine gtmg_start

      subroutine gtmg_stop (mynum)
      implicit none

      integer, intent(in) :: mynum

      integer err
      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (timing_barrier_L) then
!$OMP BARRIER
      endif
!$omp single
      if (timing_barrier_L) call rpn_comm_barrier ("GRID", err)
      if (Gem_timing_dyn_L) then
         sum_tb(mynum)= sum_tb (mynum) + (omp_get_wtime() - tb (mynum))
      else if(New_timing_dyn_L) then
         call timing_stop (mynum)
      endif
!$omp end single
!
!-------------------------------------------------------------------
!
      return
      end subroutine gtmg_stop

      subroutine gtmg_terminate ( myproc, msg )
      use ptopo
      use glb_ld
      use path
      use out_mod
      use out_options
      use gem_options
      use step_options
      implicit none

      character(len=*), intent(in) :: msg
      integer         , intent(in) :: myproc

      character(len=16) name
      character(len=64) fmt,nspace
      logical flag(MAX_instrumented)
      integer i,j,k,elem,lvl,lvlel(0:100)

      real(kind=REAL64) omp_get_wtime
!
!-------------------------------------------------------------------
!
      if (.not.New_timing_dyn_L .and. .not.Gem_timing_dyn_L) return


      if (Gem_timing_dyn_L) then

      if (myproc /= 0) return
      print *,'_______________________________________________________________'
      print *,'__________________OMP_TIMINGS ON PE #0_________________________'

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
      print *,'_______________________________________________________________'

      else
         call timing_terminate2 ( myproc, msg )
      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine gtmg_terminate
      
end module omp_timing
