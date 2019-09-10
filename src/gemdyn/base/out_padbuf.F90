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

      subroutine out_padbuf (f,minx,maxx,miny,maxy,lnk)
      use glb_ld
      implicit none
#include <arch_specific.hf>


      integer minx,maxx,miny,maxy,lnk
      real f(minx:maxx,miny:maxy,lnk)

      integer i,j,k,iff,id,jd,jf
!
!----------------------------------------------------------------------
!
      id   = 1
      jd   = 1
      iff  = l_ni
      jf   = l_nj

      do k=1,lnk
         do i= minx , id-1
            f(i,jd:jf,k)= f(id  ,jd:jf,k)
         end do
         do i= iff+1, maxx
            f(i,jd:jf,k)= f(iff ,jd:jf,k)
         end do
      end do

      do k=1,lnk
         do j= miny, jd-1
            f(minx:maxx,j,k)= f(minx:maxx,jd,k)
         end do
         do j= jf+1, maxy
            f(minx:maxx,j,k)= f(minx:maxx,jf,k)
         end do
      end do
!
!----------------------------------------------------------------------
!
      return
      end
