!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */

      subroutine RPN_COMM_compress(array,minx,maxx,miny,maxy,nil,njl,nk)
      implicit none
      integer, intent(IN) :: minx,maxx,miny,maxy,nil,njl,nk
      integer, dimension(*), intent(INOUT) :: array
!
      integer expanded, compressed
      real *8 expanded8, compressed8
      pointer (expanded_,expanded(minx:maxx,miny:maxy,nk))
      pointer (expanded8_,expanded8(minx:maxx,miny:maxy,nk))
      pointer (compressed_,compressed(nil,njl,nk))
      pointer (compressed8_,compressed8(nil,njl,nk))
!
      integer i,j,k
!
      compressed_ = loc(array)
      expanded_ = loc(array)
!
      do k=1,nk
      do j=1,njl
!VDIR NODEP
      do i=1,nil
        compressed(i,j,k)=expanded(i,j,k)
      enddo
      enddo
      enddo
!
      entry RPN_COMM_compress8(array,minx,maxx,miny,maxy,nil,njl,nk)
!
      compressed8_ = loc(array)
      expanded8_ = loc(array)
!
      do k=1,nk
      do j=1,njl
!VDIR NODEP
      do i=1,nil
        compressed8(i,j,k)=expanded8(i,j,k)
      enddo
      enddo
      enddo
!
      return
!
!
      entry RPN_COMM_expand(array,minx,maxx,miny,maxy,nil,njl,nk)
!
      compressed_ = loc(array)
      expanded_ = loc(array)
!
      do k=nk,1,-1
      do j=njl,1,-1
!VDIR NODEP
      do i=nil,1,-1
        expanded(i,j,k)=compressed(i,j,k)
      enddo
      enddo
      enddo
!
!
      entry RPN_COMM_expand8(array,minx,maxx,miny,maxy,nil,njl,nk)
!
      compressed8_ = loc(array)
      expanded8_ = loc(array)
!
      do k=nk,1,-1
      do j=njl,1,-1
!VDIR NODEP
      do i=nil,1,-1
        expanded8(i,j,k)=compressed8(i,j,k)
      enddo
      enddo
      enddo
!
      return
      end
