!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

subroutine timing_terminate2(myproc, msg)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_toupper
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   character(len=*), intent(in) :: msg
   integer, intent(in) :: myproc
   !@author M. Desgagne   -- Winter 2012 --

#include <rmnlib_basics.hf>

   include "timing.cdk"

   character(len=64) :: fmt,nspace,nspace2,nspace3,tmp1_s,tmp2_s
   logical flag(MAX_instrumented)
   integer i,j,elem,lvl,lvlel(0:100),err,maxlen
   integer(INT64),dimension(MAX_instrumented) :: timer_cnt2
   real(REAL64),dimension(MAX_instrumented) :: sum_tb_mn, sum_tb_mx, sum_tb_mn2, sum_tb_mx2
   real(REAL64) :: mymax, mymin
   real :: sum_tb_per_mn, sum_tb_per_mx

   if (Timing_S=='YES') call tmg_terminate ( myproc, msg )

   maxlen = 0
   do i=1,MAX_instrumented
      mymax = maxval(sum_tb(i,2:MAX_threads))
      sum_tb_mx(i) = sum_tb(i,1) + mymax
      mymin = mymax
      do j=2,MAX_threads
         if (sum_tb(i,j) > 0.) mymin = min(mymin, sum_tb(i,j))
      enddo
      sum_tb_mn(i)   = sum_tb(i,1)    + mymin
      timer_cnt(i,1) = timer_cnt(i,1) + maxval(timer_cnt(i,2:MAX_threads))
      maxlen = max(maxlen, len_trim(nam_subr_S(i)))
   enddo

   sum_tb_mn2 = 0.D0
   sum_tb_mx2 = 0.D0
   timer_cnt2 = 0
   call rpn_comm_reduce(sum_tb_mn,    sum_tb_mn2,    MAX_instrumented, &
        RPN_COMM_REAL8,    RPN_COMM_MIN, RPN_COMM_MASTER, RPN_COMM_GRID, err)
   call rpn_comm_reduce(sum_tb_mx,    sum_tb_mx2,    MAX_instrumented, &
        RPN_COMM_REAL8,    RPN_COMM_MAX, RPN_COMM_MASTER, RPN_COMM_GRID, err)
   call rpn_comm_reduce(timer_cnt(:,1), timer_cnt2, MAX_instrumented, &
        RPN_COMM_INTEGER8, RPN_COMM_MAX, RPN_COMM_MASTER, RPN_COMM_GRID, err)

   if (myproc.ne.0) return

   !#TODO: Add Mean and Var/Std to timings

   write(6,'(a)') '________________________________________________________________________________________'
   write(6,'(a)') '|____ TIMINGS __________________________________________________________________________|'
   write(6,'(a)') '|   |                               |   Wallclock [%]| Wallclock [Sec]         |        |'
   write(6,'(a)') '| ID| NAME                          |   Min%:    Max%| MinSec     : MaxSec     |  Count |'
   write(6,'(a)') '|---|-------------------------------|----------------|-------------------------|--------|'
   flag=.false.
   mymax = maxval(sum_tb_mx2)
   write (nspace3,'(i3)') max(5, maxlen)
   do i = 1,MAX_instrumented
      lvl= 0 ; elem= i
55    if ( (trim(nam_subr_S(elem)).ne.'') .and. (.not.flag(elem)) ) then

         err = clib_toupper(nam_subr_S(elem))
         sum_tb_per_mn = real(sum_tb_mn2(elem)) / max(0.00001,real(mymax))
         sum_tb_per_mx = real(sum_tb_mx2(elem)) / max(0.00001,real(mymax))
         write (nspace,'(i3)') 3*lvl+1
         write (nspace2,'(i3)') max(1,18-(3*lvl+1))
         fmt = '(a,i3,a,'//trim(nspace)//'x,a'//trim(nspace3)//','// &
              trim(nspace2)//'x,a,f6.2,a,f6.2,a)'
         write(tmp1_S, fmt) &
              "|", elem, "|", nam_subr_S(elem), &
              "|", 100.*sum_tb_per_mn, "%: ", 100.*sum_tb_per_mx, "%|"
         fmt = '(1pe10.4,a,1pe10.4,a,i8,a)'
         write(tmp2_S, fmt) &
              sum_tb_mn2(elem), " : ", sum_tb_mx2(elem), &
              ' |', timer_cnt2(elem), '|'
         write(6, '(a,1x,a)') trim(tmp1_S),trim(tmp2_S)
         flag(elem) = .true. ; lvlel(lvl) = elem
65       do j = 1,MAX_instrumented
            if ((timer_level(j) .eq. elem) .and. (.not.flag(j)) )then
               lvl= lvl+1
               elem= j
               goto 55
            endif
         end do
         lvl= lvl - 1
         if (lvl .ge. 0) then
            elem= lvlel(lvl)
            goto 65
         endif
      endif
   enddo

   write(6,'(a)') '________________________________________________________________________________________'

   return
end subroutine timing_terminate2
