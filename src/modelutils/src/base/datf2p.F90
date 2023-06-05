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
subroutine datf2p(mc2date,fstdate)
   implicit none
!!!#include <arch_specific.hf>
   !@object conversion of RPN datestp into mc2 date format
   !@arguments
   !   mc2date    O  date encoded in mc2 format
   !   fstdate    I  date encoded in RPN standard file format
   character(len=16), intent(out) :: mc2date
   integer :: fstdate
   !*@/
   integer, external :: newdate
   integer :: yy,mo,dd,hh,mm,ss
   integer :: dat2,dat3,err
   !-------------------------------------------------------------------
   err = newdate(fstdate,dat2,dat3,-3)
   yy = dat2/10000
   mo = mod(dat2,10000)/100
   dd = mod(dat2,100)
   hh = dat3/1000000
   mm = mod(dat3,1000000)/10000
   ss = mod(dat3,10000)/100
   write(mc2date,10) yy,mo,dd,hh,mm,ss
10 format(i4.2,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
   !-------------------------------------------------------------------
   return
end subroutine datf2p


