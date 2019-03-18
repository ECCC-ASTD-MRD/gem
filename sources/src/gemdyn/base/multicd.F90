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

      subroutine multicD (nproc,proc,myproc,type,mycol,iwk)
      implicit none
#include <arch_specific.hf>
      integer nproc,proc(nproc),myproc,type,mycol,iwk(*)
!     include 'mpif.h'
!-----------------------------------------------------------------------
!     P-SPARSLIB ROUTINE MULTICD
!     parallel node multicoloring. This code assigns a color to
!     each processor such that no two neighboring processors are
!     assigned the same color. The algorithm used is based on a topo-
!     gical sorting of the upward directed  version of the graph
!     (i.e., graph obtained by orienting all edges from lower labels
!     to higher labels). As in level scheduling, the parallelism is of
!     the order of the diameter of the graph.
!
!     written by Y. Saad, modified by A. Malevsky, January 25, 1995
! revision:
!     Abdessamad Qaddouri: adds RPN_send and RPN_receiv
!
!-----------------------------------------------------------------------
! on entry:
!---------
! nproc   = number of processors that are adjacent to my processor
!
! proc    = list of the processors adajacent to my processor.
!
! myproc  = label of my processor
!
! type    = tag to be used for sends / receives.
!
!
! on return:
! ----------
! mycol   = color assigned to my processor
!
! work space
! -----------
! iwk     = integer whose size equal the maximum number of different
!           colors assigned to adjacent processors/
!
! NOTE: processor ID's are supposed to be >= 1 in list proc.
!-----------------------------------------------------------------------
!****Feb 1996
!     modified the  code to allow nproc = 0
      integer kol,ii,k,j,low,len,ncol,status,ierr
!
!
!     if one processor is used or no adajacent processors at all
!     return as mycol = 1

      if(nproc == 0)  then
         mycol = 1
         return
      end if
      kol=0
      len = 1
!
!     determine the processors with lower id's than mine
!
      low = 1
 1    if (proc(low) < myproc) then
         low = low+1
         if (low <= nproc) goto 1
      end if
      low = low - 1
      ncol = 0
      iwk(1) = 0
!
!     receive all colors of neighbors
!
      do 10 ii = 1, low
!         call MSG_receive(proc(ii),type,kol,len,imsg)
!      call MPI_BARRIER(MPI_COMM_WORLD,imsg)

       call RPN_COMM_recv ( kol, 1, 'MPI_INTEGER', &
        proc(ii)-1,type,'grid',status,ierr)

!
!     sorted insertion -- first find where to insert
!
         j = 1
 2       if (j <= ncol .and. iwk(j) < kol) then
            j = j+1
            goto 2
         else  if (iwk(j) == kol) then
            goto 10
         end if
         j = j-1
!
         do k= j+1,ncol,-1
            iwk(k+1) = iwk(k)
         end do
         iwk(j+1) = kol
         ncol = ncol+1
 10   continue
!
!     determine my color by searching for  a gap in iwk
!
         mycol = 1
         k = 1
 3       if (iwk(k) == mycol) then
            k = k+1
            mycol = mycol+1
            if (k <= ncol) goto 3
         end if
      do 20 ii = low+1, nproc
!         call MSG_send(proc(ii),type,mycol,len,imsg)
      call RPN_COMM_send ( mycol, 1, 'MPI_INTEGER', proc(ii)-1, &
                           type, 'grid', ierr)

 20   continue
      return
!-----------------------------------------------------------------------
!-----end-of-multicD----------------------------------------------------
      end

