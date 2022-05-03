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

!/@*
subroutine difdatsd(diff,date1,date2)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@object find time difference in days between two dates
   !@arguments
   ! diff     (O) - ( date2 - datel ) in units of days
   ! date1    (I) - time in format "yyyymmdd.hhmmss"
   ! date2    (I) - time in format "yyyymmdd.hhmmss"
   real(REAL64), intent(out) :: diff
   character(len=16), intent(in) :: date1,date2
   !@author Michel Desgagne
   !@revision
   ! v3_32 - Dugas B.          - use newdate/difdatr rather than pdfjdate and
   !                             account for the possible three-hourly resolution
   !                             of these routines. The result is correct to the
   !                             second from year 0 to 9999 inclusive.
   !*@/
!!$   integer, external :: newdate

   integer :: diff_sec
   integer :: tim1,tim2, stamp1,stamp2, ier
   integer :: yy1,mo1,dd1,hh1,mm1,ss1, sign
   integer :: yy2,mo2,dd2,hh2,mm2,ss2
   !---------------------------------------------------------------

   call prsdate(yy1,mo1,dd1,hh1,mm1,ss1,sign,date1)
   call prsdate(yy2,mo2,dd2,hh2,mm2,ss2,sign,date2)

   ! stamp1 whithout mm1 and ss1 and a modulo(3) number of hours
   tim1 = (yy1*100+mo1)*100+dd1
   tim2 = (hh1-mod( hh1,3 ))*1000000
   hh1  =  mod( hh1,3 ) ! modulo(3) remainder
   ier  = newdate( stamp1, tim1,tim2, +3 )

   if (ier /= 0) then
      print *, 'difdatsd: error in call to newdate(+3), tim1,tim2= ',tim1,tim2
      stop ' in difdatsd'
   endif

   ! stamp2 whithout mm2 and ss2 and a modulo(3) number of hours
   tim1 = (yy2*100+mo2)*100+dd2
   tim2 = (hh2-mod( hh2,3 ))*1000000
   hh2  =  mod( hh2,3 ) ! modulo(3) remainder
   ier  = newdate( stamp2, tim1,tim2, +3 )

   if (ier /= 0) then
      print *, 'difdatsd: error in call to newdate(+3), tim1,tim2= ',tim1,tim2
      stop ' in difdatsd'
   endif

   ! calculate differences arising from mm/ss
   diff_sec = (mm2-mm1)*60+ss2-ss1
   if (diff_sec < 0) then
      diff_sec = 3600+diff_sec
      hh2 = hh2-1
   endif

   ! differences whithout the hh2/hh1 modulo(3)
   ! remainders and the mm/ss differences
   call difdatr( stamp2,stamp1, diff )

   ! add hh2/hh1 and mm/ss differences
   diff = (diff + hh2-hh1 + diff_sec / 3600.0_8 ) / 24.0_8
   !----------------------------------------------------------------
   return
end subroutine difdatsd
