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

subroutine pdfcdate(yyyy,mo,dd,hh,mm,ss,jdate)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   real(REAL64) :: jdate
   integer yyyy,mo,dd,hh,mm,ss,seconds

   real(REAL64) :: f,rj

   rj = int(jdate)
   f = jdate - rj
   seconds = nint(f * 86400.0)

   ss = mod(seconds, 60)
   mm = mod(seconds - ss,3600)/60


   hh = (seconds-60*mm-ss) / 3600
   if (hh.eq.24) then
      hh = 0
      seconds = seconds - 86400
      rj = rj+1.0
   endif
   mm = (seconds - hh * 3600 - ss) / 60

   call datec(int(rj),yyyy,mo,dd)

   return
end subroutine pdfcdate

