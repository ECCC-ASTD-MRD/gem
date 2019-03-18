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

!**s/r pdfjdate2 - number of days since "Step_runstrt_S(year) - 1"
      subroutine pdfjdate2 (jdate,yyyy,mo,dd,hh,mm,ss)
      use step_options
      implicit none
#include <arch_specific.hf>
      real*8 jdate
      integer yyyy,mo,dd,hh,mm,ss

!author
!    Michel Desgagne - RPN - ?????
!
!revision
! v3_32 - Dugas B.          - use newdate/difdatr rather than jd inline function
!
!arguments I/O
! jdate                 (O) - number of days since Jan 01, 00Z of "Step_runstrt_S(year) - 1"
! "yyyy mo dd hh mm ss" (I) - calendar date to compare with

!  calculate the number of days since the first day of the year before
!  this run was started. This can account for leap-years support, which
!  can be turned ON or OFF via calls to ACCEPT_LeapYear() and
!  Ignore_LeapYear(), respectively
!!

      integer, save :: stamp0
      integer ier, TIM1, stamp1
      integer, parameter :: TIM2=0, MOD=3
      logical, save :: done=.false.

      integer newdate
      external newdate, difdatr

      jdate = -1.0_8

      if (.not.done) then
         ! initialization of stamp0 to "Step_runstrt_S(year) - 1"
         read(Step_runstrt_S,'(I4)') TIM1
         if (TIM1 > 0) TIM1 = TIM1-1
         TIM1 = tim1*10000+0101
         ier = newdate( stamp0, TIM1,TIM2, MOD )
         if (ier /= 0) then
            print *, 'pdfjdate2: error in call to newdate(+3), tim1,tim2= ',TIM1,TIM2
            stop ' in pdfjdate2'
         end if
         done = .true.
      end if

      TIM1 = (yyyy*100+mo)*100+dd
      ier = newdate( stamp1, TIM1,TIM2, MOD )

      if (ier /= 0) then
         print *, 'pdfjdate2: error in call to newdate(+3), tim1,tim2= ',TIM1,TIM2
         stop ' in pdfjdate2'
      end if

      ! number of hours between "yyyy mo dd" and "Step_runstrt_S(year) - 1"
      call difdatr( stamp1,stamp0, jdate )

      ! add "hh mm ss" hours to jdate and convert to days
      jdate = (jdate+hh+(mm*60+ss)/3600.0_8)/24.0_8

      return
      end

