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
subroutine incdatsd(newdat, olddat, dt)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@object add dt days to olddat
   !@arguments
   ! olddat   (I) - time in format "yyyymmdd.hhmmss"
   ! newdat   (O) - time in format "yyyymmdd.hhmmss"
   ! dt       (I) - time interval to add in units of days to olddat (can be negative)
   character(len=16) :: newdat,olddat
   real(REAL64) :: dt
   !@author Michel Desgagne - RPN - ?????
   !@revision
   ! v3_32 - Dugas B.          - use newdate/incdatr rather than pdfcdate and pdfjdate
   !                             and account for the possible three-hourly resolution
   !                             of these routines. The result is correct to the
   !                             second from year 0 to 9999 inclusive.
   !*@/
!!$   integer, external :: newdate

   real(REAL64) :: days,hours
   real :: hour_frac

   integer :: newyy,newmo,newdd,newhh,newmm,newss
   integer :: oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign
   integer :: tim1,tim2, oldstamp,newstamp, ier

   !---------------------------------------------------------------
   call prsdate(oldyy,oldmo,olddd,oldhh,oldmm,oldss,oldsign,olddat)

   ! newhh, newmm and newss are calculated first
   hours = oldhh + (oldmm*60+oldss)/3600.0_REAL64 ! old hours
   hours = hours + dt*24.0_REAL64 ! new hours

   ! fractionnal hours in seconds
   hour_frac = nint((hours - int(hours, INT64)) * 3600.0_REAL64)
   if (hour_frac == 3600.0) then
      hour_frac = 0.0
      hours = int(hours, INT64)+1
   endif

   days  = int(hours/24.0_REAL64, INT64) ! integer number of days to add

   if (hours < 0.0_REAL64) then
      if (days*24.0_REAL64 /= hours) days = days-1
      if (hour_frac /= 0.0) hour_frac = 3600.0+hour_frac
   endif

   newhh = hours-days*24.0_REAL64! final number of hours
   newmm = hour_frac / 60.0_REAL64 ! final number of minutes
   newss = nint(hour_frac - 60.0_REAL64 * newmm, INT64)
   if (newss == 60) then
      newss = 0 ; newmm = newmm+1
   endif

   days  = days+olddd-1 ! days to add to oldyy/oldmm/01

   ! calculate the date-time-stamp for oldyy/oldmm/01
   tim1 =  (oldyy * 100 + oldmo) * 100 + 1
   tim2 =   0

   ier = newdate(oldstamp, tim1,tim2, +3)

   if (ier /= 0) then
      print *, 'incdatsd: error in call to newdate(+3), tim1,tim2= ',tim1,tim2
      stop ' in incdatsd'
   endif

   hours = int(days, INT64)*24
   ! add days*24 hours to oldstamp
   call incdatr(newstamp, oldstamp, hours)

   ier = newdate(newstamp, tim1,tim2, -3)

   if (ier /= 0) then
      print *, 'incdatsd: error in call to newdate(-3), stamp= ',newstamp
      stop ' in incdatsd'
   endif

   newyy = mod( tim1/10000   , 10000 )
   newmo = mod( tim1/100     , 100 )
   newdd = mod( tim1         , 100 )

   write(newdat,16) newyy,newmo,newdd,newhh,newmm,newss
16 format(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
   !---------------------------------------------------------------
   return
end subroutine incdatsd

