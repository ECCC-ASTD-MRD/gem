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

!**s/r gemtim4 - Timing routine

      subroutine gemtim4 ( unf, from, last )
      use step_options
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: from
      logical, intent(in) :: last
      integer, intent(in) :: unf

      integer, external :: get_max_rss

      character date_S*8 ,time_S*10, jour_S*11, heure_S*10
      logical, save :: timini_L = .false.
      real          users,systs
      real, save :: user0=0.0, syst0=0.0
      real*8, save ::  START=-1.d0, END, avgtime(10)=0.d0, &
                                 ACCUM_w=0.d0, ACCUM_u=0.d0, ACCUM_s=0.d0
      real*8 omp_get_wtime
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
      end
