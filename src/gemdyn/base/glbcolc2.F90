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

      subroutine glbcolc2 ( f2rc, g_id,g_if, g_jd,g_jf, g_kd,g_kf, &
                            f2cc,lminx,lmaxx,lminy,lmaxy,lmink,lmaxk)
      use glb_ld
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer, intent(IN) :: g_id,g_if,g_jd,g_jf,g_kd,g_kf,&
                         lminx,lmaxx,lminy,lmaxy,lmink,lmaxk
      real, intent(OUT) :: f2rc(g_id:g_if,g_jd:g_jf,g_kd:g_kf)
      real, intent(IN ) :: f2cc(lminx:lmaxx,lminy:lmaxy,lmink:lmaxk)

      integer i, j, k, iproc, tag, err, status
      integer si,sj,loindx,hiindx,loindy,hiindy
      integer len,l_id,l_if,l_jd,l_jf
      common /gatherit/ len,l_id,l_if,l_jd,l_jf
      real buf ((lmaxx-lminx+1)*(lmaxy-lminy+1)*(l_nk+1))
      data tag /210/
!
!----------------------------------------------------------------------
!
      loindx=1
      loindy=1
      hiindx=l_ni
      hiindy=l_nj
      if (l_west ) loindx = lminx
      if (l_south) loindy = lminy
      if (l_north) hiindy = lmaxy
      if (l_east ) hiindx = lmaxx
      si = Ptopo_gindx(1,Ptopo_myproc+1) - 1
      sj = Ptopo_gindx(3,Ptopo_myproc+1) - 1
      l_id = max(loindx,(g_id-si))
      l_if = min(hiindx,(g_if-si))
      l_jd = max(loindy,(g_jd-sj))
      l_jf = min(hiindy,(g_jf-sj))
      len = max(0,(l_if-l_id+1))*max(0,(l_jf-l_jd+1))*(g_kf-g_kd+1)

      if (Ptopo_myproc == 0) then

!       Copy local data (LD) segment to global field on processor 1

         if (len > 0) then
            len = 0
            do k = g_kd, g_kf
               do j = l_jd, l_jf
                  do i = l_id, l_if
                     len = len + 1
                     buf(len) = f2cc(i,j,k)
                  end do
               end do
            end do
            len = 0
            do k = g_kd, g_kf
               do j = Ptopo_gindx(3,Ptopo_myproc+1)+l_jd-1,  &
                      Ptopo_gindx(3,Ptopo_myproc+1)+l_jf-1
               do i = Ptopo_gindx(1,Ptopo_myproc+1)+l_id-1,  &
                      Ptopo_gindx(1,Ptopo_myproc+1)+l_if-1
                  len = len + 1
                  f2rc(i,j,k) = buf(len)
               end do
               end do
            end do
         end if

!       Receive the local data (LD) segments from all other processors

         do iproc = 1, Ptopo_numproc-1
            call RPN_COMM_recv ( len, 5, 'MPI_INTEGER', iproc, &
                                 tag,'GRID', status, err )
            if (len > 0) then
               call RPN_COMM_recv ( buf, len, 'MPI_REAL', iproc, &
                                 tag,'GRID', status, err )
               len = 0
               do k = g_kd, g_kf
               do j = Ptopo_gindx(3,iproc+1)+l_jd-1, Ptopo_gindx(3,iproc+1)+l_jf-1
               do i = Ptopo_gindx(1,iproc+1)+l_id-1, Ptopo_gindx(1,iproc+1)+l_if-1
                  len = len + 1
                  f2rc(i,j,k) = buf(len)
               end do
               end do
               end do
            end if
         end do

      else

!       Send local data (LD) segment to processor 1

         len = 0
         do k = g_kd, g_kf
            do j = l_jd, l_jf
            do i = l_id, l_if
               len = len + 1
               buf(len) = f2cc(i,j,k)
            end do
            end do
         end do

         call RPN_COMM_send ( len, 5, 'MPI_INTEGER', 0, tag,'GRID',err )
         if (len > 0) &
         call RPN_COMM_send ( buf, len, 'MPI_REAL', 0, tag, 'GRID',err )

      end if
!
!----------------------------------------------------------------------
!
      return
      end

