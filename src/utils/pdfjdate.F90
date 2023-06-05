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

subroutine pdfjdate(jdate,yyyy,mo,dd,hh,mm,ss)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   real(REAL64) :: jdate
   integer yyyy,mo,dd,hh,mm,ss

   !  calculate julian calendar day
   !  see cacm letter to editor by fliegel and flandern 1968
   !  page 657

   integer jd,jyy,jmo,jdd

   jd(jyy,jmo,jdd)=jdd-32075+1461*(jyy+4800+(jmo-14)/12)/4 &
        +  367*(jmo-2-(jmo-14)/12*12)/12 - 3 &
        *((jyy+4900+(jmo-14)/12)/100)/4

   jdate = jd(yyyy,mo,dd)
   jdate = jdate + (hh*3600+mm*60+ss)/86400.0

   return
end subroutine pdfjdate

