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

!**s/r  tab_vec -
!
      subroutine tab_vec (tab,Minx,Maxx,Miny,Maxy,nk, &
                                vec,i0,in,j0,jn,idir)
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,nk,i0,in,j0,jn,idir
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,nk) :: tab
      real(kind=REAL64), dimension(*) :: vec(*)
!
! author    Abdessamad Qaddouri - December 2006
!
      integer iloc,i,j,k,nij
!
!     ---------------------------------------------------------------
!
      nij=(in-i0+1)*(jn-j0+1)
!
      if (idir == 1) then
         do k = 1  , nk
            iloc=(k-1)*nij
            do j=j0,jn
               do i=i0,in
                  iloc=iloc+1
                  vec(iloc) = tab(i,j,k)
               end do
            end do
         end do
      else if (idir == (-1)) then
         do k = 1  , nk
            iloc=(k-1)*nij
            do j=j0,jn
               do i=i0,in
                  iloc=iloc+1
                  tab(i,j,k)=vec(iloc)
               end do
            end do
         end do
      end if
!
!     ---------------------------------------------------------------
!
      return
      end
