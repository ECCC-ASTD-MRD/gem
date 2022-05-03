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
subroutine datp2f (fstdate,mc2date)
   implicit none
!!!#include <arch_specific.hf>
   !@object conversion of RPN datestp into mc2 date format
   !@arguments
   !   mc2date    I  date encoded in mc2 format
   !   fstdate    O  date encoded in RPN standard file format
   character(len=*), intent(in) :: mc2date
   integer :: fstdate
   !*@/
   integer, external :: newdate
   integer :: yy,mo,dd,hh,mm,ss,dat2,dat3,err
   character(len=16) :: mc2date2
   character(len=4) :: cyy
   character(len=2) :: cmo,cdd,chh,cmm,css
   !-------------------------------------------------------------------
   mc2date2 = mc2date
   cyy=mc2date2(1:4)
   cmo=mc2date2(5:6)
   cdd=mc2date2(7:8)
   chh=mc2date2(10:11)
   cmm=mc2date2(12:13)
   css=mc2date2(14:15)

   read(cyy,'(I4)') yy
   read(cmo,'(I2)') mo
   read(cdd,'(I2)') dd
   read(chh,'(I2)') hh
   read(cmm,'(I2)') mm
   read(css,'(I2)') ss

   dat2= yy*10000 + mo*100 + dd
   dat3= hh*1000000 + mm*10000 + ss*100
   err = newdate(fstdate,dat2,dat3,3)
   !-------------------------------------------------------------------
   return
end subroutine datp2f


