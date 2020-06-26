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

      subroutine adz_export (F_export)
      use ISO_C_BINDING
      use adz_mem
      use ptopo
      implicit none
#include <arch_specific.hf>

      type(ADZ_SLOD), intent(INOUT) :: F_export

      integer tag1, tag2, tag3,tag4, ireq, kk, ierr, cnt
      integer request(16)
!
!     ---------------------------------------------------------------
!
      tag1=14 ; tag2=15 ; tag3=16 ; tag4=17 ; ireq=0

! Receive individual number of points from each neighbor
      do kk= 1, 8
         if (Adz_exppe(kk)>=0) then
            ireq = ireq+1
            call RPN_COMM_IRecv (F_export%recv(kk), 1, 'MPI_INTEGER',&
                                 Adz_exppe(kk), tag1+Adz_exppe(kk)  ,&
                                 'GRID', request(ireq), ierr)
         endif

      end do

! Send individual number of points to each neighbor
      do kk= 1, 8
         if (Adz_exppe(kk)>=0) then
            ireq = ireq+1
            call RPN_COMM_ISend (F_export%send(kk),1, 'MPI_INTEGER',&
                                 Adz_exppe(kk), tag1+Ptopo_myproc  ,&
                                 'GRID', request(ireq), ierr )
         endif
      end do

      call RPN_COMM_waitall_nostat (ireq, request, ierr)

      ireq=0 ; cnt=0
! Receive interpolation requests from each neighbor
      do kk= 1, 8
         if ((Adz_exppe(kk)>=0).and.(F_export%recv(kk)>0)) then
            ireq = ireq+1
            call RPN_COMM_IRecv ( F_export%req(cnt+1), &
                                 3*F_export%recv(kk), 'MPI_REAL', &
                                 Adz_exppe(kk),tag2+Adz_exppe(kk),&
                                 'GRID',request(ireq),ierr )
            cnt= cnt+3*F_export%recv(kk)
         endif
      end do

! Send interpolation requests to each neighbor
      do kk= 1, 8
         if ((Adz_exppe(kk)>=0).and.(F_export%send(kk)>0)) then
            ireq = ireq+1
            call RPN_COMM_ISend ( F_export%gpos(1,KK)            ,&
                                  3*F_export%send(kk), 'MPI_REAL',&
                                  Adz_exppe(kk), tag2+Ptopo_myproc,&
                                  'GRID', request(ireq), ierr )
         endif
      end do

      F_export%COMM_handle(1) = ireq
      F_export%COMM_handle(2:ireq+1) = request(1:ireq)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_export
