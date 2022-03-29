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

!**s/r glbpos - sets up global indices for each PE

      subroutine glbpos()
      use glb_ld
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer :: i, dimens, err
      integer, dimension(6,Ptopo_numproc) :: gindx
!
!----------------------------------------------------------------------
!
      gindx = 0

      gindx(1,Ptopo_myproc+1) = l_i0
      gindx(2,Ptopo_myproc+1) = l_i0 + l_ni - 1
      gindx(3,Ptopo_myproc+1) = l_j0
      gindx(4,Ptopo_myproc+1) = l_j0 + l_nj - 1
      gindx(5,Ptopo_myproc+1) = 1
      gindx(6,Ptopo_myproc+1) = G_nk

      if ( .not. allocated(Ptopo_gindx) ) then
         allocate (Ptopo_gindx(6,Ptopo_numproc))
      end if
      if ( .not. allocated(Ptopo_gindx_alongX) ) then
         allocate (Ptopo_gindx_alongX(2,Ptopo_npex))
      end if
      if ( .not. allocated(Ptopo_gindx_alongY) ) then
         allocate (Ptopo_gindx_alongY(2,Ptopo_npey))
      end if

      Ptopo_gindx = 0
      dimens = 6*Ptopo_numproc
      call rpn_comm_ALLREDUCE (gindx,Ptopo_gindx,dimens,"MPI_INTEGER", &
                               "MPI_BOR","grid",err)

      do i=0, Ptopo_npex-1
         Ptopo_gindx_alongX(1,i+1) = Ptopo_gindx(1,Ptopo_colrow(0,i,0)+1)
         Ptopo_gindx_alongX(2,i+1) = Ptopo_gindx(2,Ptopo_colrow(0,i,0)+1)
      end do
      do i=0, Ptopo_npey-1
         Ptopo_gindx_alongY(1,i+1) = Ptopo_gindx(3,Ptopo_colrow(0,0,i)+1)
         Ptopo_gindx_alongY(2,i+1) = Ptopo_gindx(4,Ptopo_colrow(0,0,i)+1)
      end do
!
!----------------------------------------------------------------------
      return
      end

