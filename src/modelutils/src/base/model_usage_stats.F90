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

!/@*
subroutine model_usage_stats(F_msg_S,F_print_accum_L)
   implicit none
!!!#include <arch_specific.hf>
   character(len=*), intent(in) :: F_msg_S
   logical, intent(in) :: F_print_accum_L
!*@/
#include <rmn/msg.h>

   integer, external :: get_max_rss
   double precision, external :: omp_get_wtime
   
   logical, save :: timini_L = .false.
   real,    save :: user0    = 0.0
   real,    save :: syst0    = 0.0
   double precision, save :: wtime_start = -1.d0
   double precision, save :: avgtime(10) =  0.d0
   double precision, save :: accum_w = 0.d0
   double precision, save :: accum_u = 0.d0
   double precision, save :: accum_s = 0.d0

   logical :: canWrite_L
   character(len=8)  :: date_S
   character(len=10) :: time_S, heure_S
   character(len=11) :: jour_S
   character(len=256) :: tmp_S
   integer   :: msgLevelMin, msgUnit
   real      :: users, systs
   double precision :: wtime_end
   !----------------------------------------------------------------
   call msg_getInfo(canWrite_L, msgLevelMin, msgUnit, tmp_S)
   if (MSG_INFO < msgLevelMin .or. .not.canWrite_L) return

   if (wtime_start < 0.d0) wtime_start = omp_get_wtime()

   call date_and_time(date_S, time_S)
   jour_S  = date_S(1:4)//'/'//date_S(5:6)//'/'//date_S(7:8)//' '
   heure_S = time_S(1:2)//'h'//time_S(3:4)//'m'//time_S(5:6)//'s,'

#if defined (AIX)
   call setrteopts('cpu_time_type=total_usertime')
#endif
   call cpu_time(users)

#if defined (AIX)
   call setrteopts('cpu_time_type=total_systime')
#endif
   call cpu_time(systs)

   wtime_end = omp_get_wtime()

   if (timini_L) then
      if (trim(F_msg_S) =='CURRENT TIMESTEP') then
         avgtime(1:9) = avgtime(2:10)
         avgtime(10)  = wtime_end - wtime_start
!!$         Step_maxwall = sum(avgtime) / 10.d0
      endif

1000 format(/A,' W: ',1pe13.6, &
        ' U: ',1pe13.6, &
        ' S: ',1pe13.6, &
        ', Mem: ',i9,' (Kbytes/PE) ',a)
      write(msgUnit,1000) 'TIME: '//jour_S//heure_S, &
           wtime_end - wtime_start, users - user0, &
           systs - syst0, get_max_rss(), trim(F_msg_S)

      accum_w= accum_w + wtime_end - wtime_start
      accum_u= accum_u + users - user0
      accum_s= accum_s + systs - syst0

      if (F_print_accum_L) then
1001     format(/'ACCUMULATED TIME: W: ',1pe13.6, &
              ' U: ',1pe13.6, &
              ' S: ',1pe13.6)
         write(msgUnit,1001) accum_w, accum_u, accum_s
      endif
!!$   else
!!$      Step_maxwall = -1.d0
   endif

   user0 = users
   syst0 = systs
   wtime_start = wtime_end
   timini_L = .true.
   !----------------------------------------------------------------
   return
end subroutine model_usage_stats
