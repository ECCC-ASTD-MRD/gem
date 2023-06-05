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

module mu_jdate_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use iso_c_binding
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   private
   !@objective Manipulate dates from Gregorian and "cmcdate" formats to "Julian second"
   !@author Stephane Chamberland, 2016-01
   !@description
   ! Public functions
   public :: mu_set_leap_year, mu_jd2ymd, mu_ymd2jd, &
        mu_jd2ymd_leap, mu_ymd2jd_leap, &
        mu_jd2ymd_noleap, mu_ymd2jd_noleap, &
        mu_jdhms2js, mu_js2jdhms, &
        mu_ymdhms2js, mu_js2ymdhms, &
        mu_ymdhms2js_leap, mu_js2ymdhms_leap, &
        mu_ymdhms2js_noleap, mu_js2ymdhms_noleap, &
        mu_is_leapyear, mu_days_in_month_noleap, mu_days_in_month, &
        mu_day_of_year, &
        jdate_from_cmc, jdate_to_cmc, &
        jdate_from_print, jdate_to_print, &
        jdate_year, jdate_month, jdate_day_of_month, &
        jdate_day_of_year, jdate_midmonth
   ! Public constants
   integer(INT64), parameter, public :: MU_JDATE_ANY = RMN_ANY_DATE !# -1
   integer, parameter, public :: MU_JDATE_LEAP_ON = 0
   integer, parameter, public :: MU_JDATE_LEAP_IGNORED = 1
   integer, parameter, public :: MU_JDATE_PDF_LEN = 16
   integer(INT64), parameter, public :: MU_JDATE_MAX_INT = 464269103999_IDOUBLE
   character(len=MU_JDATE_PDF_LEN), parameter, public :: MU_JDATE_MAX_STR = '99991231.235959'
   integer(INT64), parameter, public :: MU_JDATE_EPOCH_INT = 0
   character(len=MU_JDATE_PDF_LEN), parameter, public :: MU_JDATE_EPOCH_STR = '-47131124.000000'
   !@/

   interface

      subroutine mu_set_leap_year(flag) bind(C, NAME='mu_set_leap_year')
         use iso_c_binding
         implicit none
         integer(C_INT), value :: flag
      end subroutine mu_set_leap_year


      subroutine mu_jd2ymd(jd, yy, mo, dd) bind(C, NAME='mu_jd2ymd')
         use iso_c_binding
         implicit none
         integer(C_LONG), value :: jd
         integer(C_INT) :: yy, mo, dd
      end subroutine mu_jd2ymd

      subroutine mu_ymd2jd(jd, yy, mo, dd) bind(C, NAME='mu_ymd2jd')
         use iso_c_binding
         implicit none
         integer(C_LONG) :: jd
         integer(C_INT), value :: yy, mo, dd
      end subroutine mu_ymd2jd


      subroutine mu_jd2ymd_leap(jd, yy, mo, dd) bind(C, NAME='mu_jd2ymd_leap')
         use iso_c_binding
         implicit none
         integer(C_LONG), value :: jd
         integer(C_INT) :: yy, mo, dd
      end subroutine mu_jd2ymd_leap

      subroutine mu_ymd2jd_leap(jd, yy, mo, dd) bind(C, NAME='mu_ymd2jd_leap')
         use iso_c_binding
         implicit none
         integer(C_LONG) :: jd
         integer(C_INT), value :: yy, mo, dd
      end subroutine mu_ymd2jd_leap


      subroutine mu_jd2ymd_noleap(jd, yy, mo, dd) bind(C, NAME='mu_jd2ymd_noleap')
         use iso_c_binding
         implicit none
         integer(C_LONG), value :: jd
         integer(C_INT) :: yy, mo, dd
      end subroutine mu_jd2ymd_noleap

      subroutine mu_ymd2jd_noleap(jd, yy, mo, dd) bind(C, NAME='mu_ymd2jd_noleap')
         use iso_c_binding
         implicit none
         integer(C_LONG) :: jd
         integer(C_INT), value :: yy, mo, dd
      end subroutine mu_ymd2jd_noleap


      subroutine mu_jdhms2js(js, jd, hh, mn, ss) bind(C, NAME='mu_jdhms2js')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG) :: js
         integer(C_LONG), value :: jd
         integer(C_INT),  value :: hh, mn, ss
      end subroutine mu_jdhms2js

      subroutine mu_js2jdhms(js, jd, hh, mn, ss) bind(C, NAME='mu_js2jdhms')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG), value :: js
         integer(C_LONG) :: jd
         integer(C_INT)  :: hh, mn, ss
      end subroutine mu_js2jdhms


      subroutine mu_ymdhms2js(js, yy, mo, dd, hh, mn, ss) &
           bind(C, NAME='mu_ymdhms2js')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG) :: js
         integer(C_INT), value :: yy, mo, dd, hh, mn, ss
      end subroutine mu_ymdhms2js

      subroutine mu_js2ymdhms(js, yy, mo, dd, hh, mn, ss) &
           bind(C, NAME='mu_js2ymdhms')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG), value :: js
         integer(C_INT) :: yy, mo, dd, hh, mn, ss
      end subroutine mu_js2ymdhms


      subroutine mu_ymdhms2js_leap(js, yy, mo, dd, hh, mn, ss) &
           bind(C, NAME='mu_ymdhms2js_leap')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG) :: js
         integer(C_INT), value :: yy, mo, dd, hh, mn, ss
      end subroutine mu_ymdhms2js_leap

      subroutine mu_js2ymdhms_leap(js, yy, mo, dd, hh, mn, ss) &
           bind(C, NAME='mu_js2ymdhms_leap')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG), value :: js
         integer(C_INT) :: yy, mo, dd, hh, mn, ss
      end subroutine mu_js2ymdhms_leap


      subroutine mu_ymdhms2js_noleap(js, yy, mo, dd, hh, mn, ss) &
           bind(C, NAME='mu_ymdhms2js_noleap')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG) :: js
         integer(C_INT), value :: yy, mo, dd, hh, mn, ss
      end subroutine mu_ymdhms2js_noleap

      subroutine mu_js2ymdhms_noleap(js, yy, mo, dd, hh, mn, ss) &
           bind(C, NAME='mu_js2ymdhms_noleap')
         use iso_c_binding
         implicit none
         integer(C_LONG_LONG), value :: js
         integer(C_INT) :: yy, mo, dd, hh, mn, ss
      end subroutine mu_js2ymdhms_noleap


      logical(C_BOOL) function mu_is_leapyear(yy) &
           bind(C, NAME='mu_is_leapyear')
         use iso_c_binding
         implicit none
         integer(C_INT), value :: yy
      end function mu_is_leapyear

      integer(C_INT) function mu_days_in_month_noleap(mo) &
           bind(C, NAME='mu_days_in_month_noleap')
         use iso_c_binding
         implicit none
         integer(C_INT), value :: mo
      end function mu_days_in_month_noleap

      integer(C_INT) function mu_days_in_month(yy, mo) &
           bind(C, NAME='mu_days_in_month')
         use iso_c_binding
         implicit none
         integer(C_INT), value :: yy, mo
      end function mu_days_in_month

      integer(C_INT) function mu_day_of_year(yy, mo, dd) &
           bind(C, NAME='mu_day_of_year')
         use iso_c_binding
         implicit none
         integer(C_INT), value :: yy, mo, dd
      end function mu_day_of_year

   end interface


contains


   !/@*
   function jdate_from_cmc(F_cmcdate) result(F_jsec_8)
      implicit none
      !@Objective Convert CMC date-time stamp to 'Julian date in seconds'
      integer, intent(in) :: F_cmcdate
      integer(INT64) :: F_jsec_8
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer :: dat2,dat3,istat,date0
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_cmcdate == RMN_ANY_DATE) then
         F_jsec_8 = MU_JDATE_ANY
         return
      endif
      if (F_cmcdate == 0) then
         F_jsec_8 = 0
         return
      endif
      F_jsec_8 = RMN_ERR
      if (F_cmcdate <= 0) return

      date0 = F_cmcdate
      istat = newdate(date0,dat2,dat3,RMN_DATE_STAMP2PRINT)
      if (.not.RMN_IS_OK(istat)) return

      yy = dat2/10000
      mo = mod(dat2,10000)/100
      dd = mod(dat2,100)
      hh = dat3/1000000
      mn = mod(dat3,1000000)/10000
      ss = mod(dat3,10000)/100

      jsec_8_c = RMN_ERR
      call mu_ymdhms2js(jsec_8_c, yy, mo, dd, hh, mn, ss)
      F_jsec_8 = jsec_8_c
      !------------------------------------------------------------------
      return
   end function jdate_from_cmc


   !/@*
   function jdate_to_cmc(F_jsec_8) result(F_cmcdate)
      implicit none
      !@Objective Convert 'Julian date in seconds' to CMC date-time stamp
      integer(INT64), intent(in) :: F_jsec_8
      integer :: F_cmcdate
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer :: dat2,dat3,istat
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_jsec_8 == MU_JDATE_ANY) then
         F_cmcdate = RMN_ANY_DATE
         return
      endif
      if (F_jsec_8 == 0) then
         F_cmcdate = 0
         return
      endif
      F_cmcdate = RMN_ERR
!!$      if (F_jsec_8 < 0) return !#jsec is <0 for noleap

      jsec_8_c = F_jsec_8
      yy=0; mo=0; dd=0; hh=0; mn=0; ss=0
      call mu_js2ymdhms(jsec_8_c, yy, mo, dd, hh, mn, ss)

      dat2 = yy*10000 + mo*100 + dd
      dat3 = hh*1000000 + mn*10000 + ss*100
      istat = newdate(F_cmcdate,dat2,dat3,RMN_DATE_PRINT2STAMP)
      if (.not.RMN_IS_OK(istat)) F_cmcdate = RMN_ERR
      !------------------------------------------------------------------
      return
   end function jdate_to_cmc


   !/@*
   function jdate_from_print(F_date_S) result(F_jsec_8)
      character(len=*),intent(in) :: F_date_S
      integer(INT64) :: F_jsec_8
      !@Objective Convert string YYYYMMDD.hhmmss to 'Julian date in seconds'
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer :: istat, a
      character(len=MU_JDATE_PDF_LEN) :: date_S
      integer(C_LONG_LONG) :: jsec_8_c
      !-------------------------------------------------------------------
      F_jsec_8 = RMN_ERR
      date_S = adjustl(F_date_S)
      if (date_S == '-1') then
         F_jsec_8 = MU_JDATE_ANY
         return
      endif
      if (len_trim(date_S) < 8) return

      a = 0
      if (date_S(1:1) == '-') then
         read(date_S(1:5),'(I5)',iostat=istat) yy
         a = 1
      else
         read(date_S(1:4),'(I4)',iostat=istat) yy
      endif
      if (.not.RMN_IS_OK(istat)) return

      read(date_S(5+a:6+a),'(I2)',iostat=istat) mo
      if (.not.RMN_IS_OK(istat)) return

      read(date_S(7+a:8+a),'(I2)',iostat=istat) dd
      if (.not.RMN_IS_OK(istat)) return

      hh = 0
      if (len_trim(date_S) >= 11+a) then
         read(date_S(10+a:11+a),'(I2)',iostat=istat) hh
      endif

      mn = 0
      if (len_trim(date_S) >= 13+a) then
         read(date_S(12+a:13+a),'(I2)',iostat=istat) mn
      endif

      ss = 0
      if (len_trim(date_S) >= 15+a) then
        read(date_S(14+a:15+a),'(I2)',iostat=istat) ss
      endif

      jsec_8_c = RMN_ERR
      call mu_ymdhms2js(jsec_8_c, yy, mo, dd, hh, mn, ss)
      F_jsec_8 = jsec_8_c
      !------------------------------------------------------------------
      return
   end function jdate_from_print


   !/@*
   function jdate_to_print(F_jsec_8) result(F_date_S)
      implicit none
      !@Objective Convert 'Julian date in seconds' to string YYYYMMDD.hhmmss
      integer(INT64), intent(in) :: F_jsec_8
      character(len=MU_JDATE_PDF_LEN) :: F_date_S
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_jsec_8 == MU_JDATE_ANY) then
         F_date_S = '-1'
         return
      endif
      F_date_S = ' '

      jsec_8_c = F_jsec_8
      yy=0; mo=0; dd=0; hh=0; mn=0; ss=0
      call mu_js2ymdhms(jsec_8_c, yy, mo, dd, hh, mn, ss)

      write(F_date_S,"(i5.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)") yy,mo,dd,hh,mn,ss
      F_date_S = trim(adjustl(F_date_S))
      !------------------------------------------------------------------
      return
   end function jdate_to_print


   !/@*
   function jdate_year(F_jsec_8) result(F_year)
      implicit none
      !@Objective Return year value from 'Julian date in seconds'
      integer(INT64), intent(in) :: F_jsec_8
      integer :: F_year
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_jsec_8 == MU_JDATE_ANY) then
         F_year = RMN_ERR
         return
      endif
      jsec_8_c = F_jsec_8
      yy=0; mo=0; dd=0; hh=0; mn=0; ss=0
      call mu_js2ymdhms(jsec_8_c, yy, mo, dd, hh, mn, ss)
      F_year = yy
      !------------------------------------------------------------------
      return
   end function jdate_year


   !/@*
   function jdate_month(F_jsec_8) result(F_month)
      implicit none
      !@Objective Return month value from 'Julian date in seconds'
      integer(INT64), intent(in) :: F_jsec_8
      integer :: F_month
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_jsec_8 == MU_JDATE_ANY) then
         F_month = RMN_ERR
         return
      endif
      jsec_8_c = F_jsec_8
      yy=0; mo=0; dd=0; hh=0; mn=0; ss=0
      call mu_js2ymdhms(jsec_8_c, yy, mo, dd, hh, mn, ss)
      F_month = mo
      !------------------------------------------------------------------
      return
   end function jdate_month


   !/@*
   function jdate_day_of_month(F_jsec_8) result(F_day)
      implicit none
      !@Objective Return day of month value from 'Julian date in seconds'
      integer(INT64), intent(in) :: F_jsec_8
      integer :: F_day
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_jsec_8 == MU_JDATE_ANY) then
         F_day = RMN_ERR
         return
      endif
      jsec_8_c = F_jsec_8
      yy=0; mo=0; dd=0; hh=0; mn=0; ss=0
      call mu_js2ymdhms(jsec_8_c, yy, mo, dd, hh, mn, ss)
      F_day = dd
      !------------------------------------------------------------------
      return
   end function jdate_day_of_month


   !/@*
   function jdate_day_of_year(F_jsec_8) result(F_day)
      implicit none
      !@Objective Return day of month value from 'Julian date in seconds'
      integer(INT64), intent(in) :: F_jsec_8
      integer :: F_day
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss
      integer(C_LONG_LONG) :: jsec_8_c
      !------------------------------------------------------------------
      if (F_jsec_8 == MU_JDATE_ANY) then
         F_day = RMN_ERR
         return
      endif
      jsec_8_c = F_jsec_8
      yy=0; mo=0; dd=0; hh=0; mn=0; ss=0
      call mu_js2ymdhms(jsec_8_c, yy, mo, dd, hh, mn, ss)
      F_day = mu_day_of_year(yy, mo, dd)
      !------------------------------------------------------------------
      return
   end function jdate_day_of_year


   !/@*
   function jdate_midmonth(F_year, F_month) result(F_jsec2_8)
      implicit none
      !@Objective Return mid month jsec value for provided year,month
      integer,intent(in) :: F_year, F_month
      integer(INT64) :: F_jsec2_8
      !*@/
      integer(C_INT) :: yy,mo,dd,hh,mn,ss, yy2,mo2
      integer(C_LONG_LONG) :: jsec_8_c, jsec2_8_c
      !------------------------------------------------------------------
      F_jsec2_8 = RMN_ERR
      if (F_month < 1 .or. F_month > 12) return

      yy=F_year; mo=F_month; dd=1; hh=0; mn=0; ss=0
      yy2=F_year; mo2=F_month+1
      if (F_month == 12) then
         yy2=F_year+1; mo2=1
      endif

      jsec_8_c = RMN_ERR ; jsec2_8_c = RMN_ERR
      call mu_ymdhms2js(jsec_8_c,  yy,  mo,  dd, hh, mn, ss)
      call mu_ymdhms2js(jsec2_8_c, yy2, mo2, dd, hh, mn, ss)

      if (jsec_8_c == RMN_ERR .or. jsec2_8_c == RMN_ERR) return

      F_jsec2_8 = jsec_8_c + (jsec2_8_c - jsec_8_c)/2
      !------------------------------------------------------------------
      return
   end function jdate_midmonth


end module mu_jdate_mod
