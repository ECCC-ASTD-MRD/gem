
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

#include <rmn/msg.h>

!/@
module cmcdate_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   implicit none
   private
   !@objective Manipulate CMC dates
   !@author Stephane Chamberland,2011-08
   !@description
   ! Public functions
   public :: cmcdate_toprint,cmcdate_fromprint,cmcdate_month, cmcdate_year,cmcdate_midmonth
   ! Public constants
   !@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

contains

   !/@*
   function cmcdate_toprint(F_date) result(F_date_S)
      integer,intent(in) :: F_date
      character(len=16) :: F_date_S
      !*@/
      integer yy,mo,dd,hh,mm,ss
      integer dat2,dat3,istat,date0
      !------------------------------------------------------------------
      F_date_S = ' '
      if (F_date <= 0) then
         write(F_date_S,"(i16)") F_date
         return
      endif
      date0 = F_date
      istat = newdate(date0,dat2,dat3,RMN_DATE_STAMP2PRINT)
      if (.not.RMN_IS_OK(istat)) return

      yy = dat2/10000
      mo = mod(dat2,10000)/100
      dd = mod(dat2,100)
      hh = dat3/1000000
      mm = mod(dat3,1000000)/10000
      ss = mod(dat3,10000)/100

      write(F_date_S,"(i4.2,i2.2,i2.2,'.',i2.2,i2.2,i2.2)") yy,mo,dd,hh,mm,ss
      !------------------------------------------------------------------
      return
   end function cmcdate_toprint


   !/@*
   function cmcdate_fromprint(F_date_S) result(F_date)
      character(len=*),intent(in) :: F_date_S
      integer :: F_date
      !*@/
      integer :: yy,mo,dd,hh,mm,ss,dat2,dat3,istat
      character(len=16) :: date_S
      !-------------------------------------------------------------------
      date_S = adjustl(F_date_S)
      F_date = RMN_ERR
      if (len_trim(F_date_S) < 8) return

      read(date_S(1:4),'(I4)',iostat=istat) yy
      if (.not.RMN_IS_OK(istat)) return

      read(date_S(5:6),'(I2)',iostat=istat) mo
      if (.not.RMN_IS_OK(istat)) return

      read(date_S(7:8),'(I2)',iostat=istat) dd
      if (.not.RMN_IS_OK(istat)) return

      hh = 0
      if (len_trim(date_S) >= 11) then
         read(date_S(10:11),'(I2)',iostat=istat) hh
      endif

      mm = 0
      if (len_trim(date_S) >= 13) then
         read(date_S(12:13),'(I2)',iostat=istat) mm
      endif

      ss = 0
      if (len_trim(date_S) >= 15) then
        read(date_S(14:15),'(I2)',iostat=istat) ss
      endif

      dat2= yy*10000 + mo*100 + dd
      dat3= hh*1000000 + mm*10000 + ss*100
      istat = newdate(F_date,dat2,dat3,RMN_DATE_PRINT2STAMP)
     !------------------------------------------------------------------
      return
   end function cmcdate_fromprint


   !/@*
   function cmcdate_month(F_date) result(F_month)
      implicit none
      integer,intent(in) :: F_date
      integer :: F_month
      !*@/
      character(len=64) :: date_S,month_S
      integer :: month,istat
      !------------------------------------------------------------------
      F_month = RMN_ERR
      if (F_date <= 0) return
      date_S = cmcdate_toprint(F_date)
      month_S = date_S(5:6)
      read(month_S,*,iostat=istat) month
      if (istat == 0) F_month = month
      !------------------------------------------------------------------
      return
   end function cmcdate_month


   !/@*
   function cmcdate_year(F_date) result(F_year)
      implicit none
      integer,intent(in) :: F_date
      integer :: F_year
      !*@/
      character(len=64) :: date_S,year_S
      integer :: year,istat
      !------------------------------------------------------------------
      F_year = RMN_ERR
      if (F_date <= 0) return
      date_S = cmcdate_toprint(F_date)
      year_S = date_S(1:4)
      read(year_S,*,iostat=istat) year
      if (istat == 0) F_year = year
      !------------------------------------------------------------------
      return
   end function cmcdate_year


   !/@*
   function cmcdate_midmonth(F_year,F_month) result(F_dateh)
      implicit none
      integer,intent(in) :: F_year,F_month
      integer :: F_dateh
      !*@/
      character(len=15) :: month_S,year_S,date1_S,date2_S
      integer :: month,year,date1,date2
      real(REAL64) :: nhours_8
      !------------------------------------------------------------------
      F_dateh = RMN_ERR
      if (F_year <= 0 .or. F_year > 3000 .or. F_month < 1 .or. F_month > 12) return
      write(month_S,'(I2.2)') F_month
      write(year_S,'(I4.4)') F_year
      date1_S = year_S(1:4)//month_S(1:2)//'01.000000'

      month = F_month + 1
      year = F_year
      if (month > 12) then
         month = 1
         year = year + 1
      endif
      write(month_S,'(I2.2)') month
      write(year_S,'(I4.4)') year
      date2_S = year_S(1:4)//month_S(1:2)//'01.000000'

      date1 = cmcdate_fromprint(date1_S)
      date2 = cmcdate_fromprint(date2_S)
      call difdatr(date1,date2,nhours_8)
      nhours_8 = abs(0.5D0 * nhours_8)
      call incdatr(F_dateh,date1,nhours_8)
      !------------------------------------------------------------------
      return
   end function cmcdate_midmonth


end module cmcdate_mod
