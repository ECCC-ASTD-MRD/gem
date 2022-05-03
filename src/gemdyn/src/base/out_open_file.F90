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

!**s/r out_open_file - open an output FST file

      subroutine out_open_file ( F_prefix_S )
      use timestr_mod, only: timestr_prognum
      use step_options
      use gem_options
      use out_options
      use cstv
      use out_mod
      use out3
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in ) :: F_prefix_S

#include <rmnlib_basics.hf>

      character(len=4) :: unit_ext
      character(len=8) :: my_block
      character(len=16) :: datev,fdate
      character(len=1024) :: filen,myformat_S,my_hour
      integer prognum,err,i,indx,len0,len1
      real(kind=REAL64), parameter :: OV_day = 1.0d0/86400.0d0
      real(kind=REAL64)  dayfrac
!
!------------------------------------------------------------------
!
      if ( (Ptopo_couleur > 0) .or. (Out3_iome < 0) )  return

      call datf2p (fdate, Out3_date)

      if ( Out_unf > 0 ) return

      write(my_block,'(a,i3.3,a,i3.3)') '-',Ptopo_mycol,'-',Ptopo_myrow

      if ( trim(F_prefix_S) == 'casc' ) then

         dayfrac = dble(lctl_step-Step_delay) * Cstv_dt_8 * OV_day
         call incdatsd (datev,Step_runstrt_S,dayfrac)
         filen = trim(F_prefix_S)//'_'//trim(datev)//my_block

      else

         err= timestr_prognum (prognum, Out3_unit_S, Out3_close_interval,&
                          Out_dateo,sngl(Cstv_dt_8),lctl_step,Out_endstepno)
         unit_ext = ' '
         if (Out3_unit_S(1:3) == 'SEC') unit_ext = 's'
         if (Out3_unit_S(1:3) == 'MIN') unit_ext = 'm'
         if (Out3_unit_S(1:3) == 'DAY') unit_ext = 'd'
         if (Out3_unit_S(1:3) == 'STE') unit_ext = 'p'
         if (Out3_unit_S(1:3) == 'MON') unit_ext = 'n'

         len1 = max(3,Out3_ndigits)
         if (any(Out3_unit_s(1:3) == (/'SEC','STE'/))) len1 = max(6,len1)
         write(my_hour,'(i20)')  abs(prognum)
         indx=1
 77      i=index(my_hour(indx:20)," ")
         if ((i > 0).and.(indx <= 20)) then
            indx=indx+1
            goto 77
         end if
         len1 = max(len1,max(0,21-indx))
         len0 = len1
         if (prognum < 0) len0 = len0+1
         write(myformat_S,'(a,i1.1,a,i1.1,a)') '(a,i',len0,'.',len1,')'
         my_hour = ' '
         write(my_hour,trim(myformat_S)) '_', prognum

         filen= trim(F_prefix_S)//fdate(1:8)//fdate(10:11)// &
                     my_block//trim(my_hour)//trim(unit_ext)

      end if

      filen= trim(Out_dirname_S)//'/'//trim(filen)

      err= fnom  ( Out_unf, trim(filen), 'STD+RND', 0 )
      err= fstouv( Out_unf, 'RND' )
!
!------------------------------------------------------------------
!
      return
      end

