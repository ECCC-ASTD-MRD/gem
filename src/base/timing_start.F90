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

subroutine timing_start2(mynum, myname_S, mylevel)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   integer, intent(in) :: mynum,mylevel
   character(len=*), intent(in) :: myname_S
   !@author M. Desgagne   -- Winter 2012 --
   !@revision
   ! v4_40 - Desgagne - initial version
   ! v4_80 - Desgagne - introduce timer_level and timer_cnt

#include <rmnlib_basics.hf>
   include "timing.cdk"

   DOUBLE PRECISION, external :: omp_get_wtime
   integer :: mynum2

   mynum2 = max(1, min(mynum, MAX_instrumented))
   if (mynum /= mynum2) &
        print *, 'WARNING: (timing_start2) called with mnum=', mynum, ' > MAX_instrumented=', MAX_instrumented

   if (Timing_S=='YES') call tmg_start ( mynum2, myname_S )

   nam_subr_S(mynum2)   = myname_S ; timer_level(mynum2) = mylevel
   tb        (mynum2,1) = omp_get_wtime()
   timer_cnt (mynum2,1) = timer_cnt(mynum2,1) + 1

   return
end subroutine timing_start2


subroutine timing_start_omp(mynum, myname_S, mylevel)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   implicit none
!!!#include <arch_specific.hf>
   integer, intent(in) :: mynum, mylevel
   character(len=*), intent(in) :: myname_S
#include <rmnlib_basics.hf>
   include "timing.cdk"
   DOUBLE PRECISION, external :: omp_get_wtime
   integer, external :: omp_get_thread_num
   integer :: t, mynum2

   t = min(max(0, omp_get_thread_num()) + 2, MAX_threads)
   mynum2 = max(1, min(mynum, MAX_instrumented))
   if (mynum /= mynum2) then
      print *, 'WARNING: (timing_start_omp) called with mnum=', &
           mynum, ' > MAX_instrumented=', MAX_instrumented
   endif

   !$omp single
   if (Timing_S == 'YES') call tmg_start(mynum2, myname_S)
   nam_subr_S(mynum2)  = myname_S
   timer_level(mynum2) = mylevel
   !$omp end single nowait

   tb(mynum2,t) = omp_get_wtime()
   timer_cnt(mynum2,t) = timer_cnt(mynum2,t) + 1

   return
end subroutine timing_start_omp

!#TODO: MPI version
