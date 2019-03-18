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

!**s/r glbdist
!

!
      subroutine glbdist (F_2bc,bni,bnj,F_2rc,Minx,Maxx,Miny,Maxy,nk,hx,hy)
!
      use glb_ld
      implicit none
#include <arch_specific.hf>
!
      integer bni,bnj,Minx,Maxx,Miny,Maxy,nk,hx,hy
      real F_2bc(bni,bnj,nk), F_2rc(Minx:Maxx,Miny:Maxy,nk)
!
!author
!     M. Desgagne
!
!revision
! v2_00 - Desgagne M.       - initial MPI version (from MC2 v4.9)
! v2_21 - Desgagne M.       - rpn_comm stooge for glbdist
!
!object
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_2bc         I           Global array to distribute
! bni,bnj       I           Horizontal dimension of F_2bc
! F_2rc         O           Local reception array
! Minx,Maxx,Miny,Maxy      I           Horizontal dimension of F_2rc
! nk            I           Vertical dimension of F_2bc and F_2rc
!----------------------------------------------------------------
!
      integer err
!----------------------------------------------------------------------
!
      call RPN_COMM_dist (F_2bc,1,bni,1,bnj,bni,bnj,nk,0,0,1, &
                          F_2rc,Minx,Maxx,Miny,Maxy,hx,hy,G_periodx,G_periody,err)
!
!----------------------------------------------------------------------
!
      return
      end

