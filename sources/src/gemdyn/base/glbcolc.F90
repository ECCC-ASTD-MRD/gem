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

!**s/r glbcolc
!
      subroutine glbcolc (F_2rc,bni,bnj,F_2cc,Minx,Maxx,Miny,Maxy,nk)
      implicit none
#include <arch_specific.hf>

      integer bni,bnj,Minx,Maxx,Miny,Maxy,nk
      real F_2rc(bni,bnj,nk), F_2cc(Minx:Maxx,Miny:Maxy,nk)

!author
!     M. Desgagne - v. lee
!
!revision
! v2_00 - Desgagne M.       - initial MPI version (from MC2 v4.9)
! v2_21 - Desgagne M.       - rpn_comm stooge for glbcolc
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_2rc         O           Global reception array
! bni,bnj       I           Horizontal dimension of F_2rc
! F_2cc         O           Local array to collect
! Minx,Maxx,Miny,Maxy       Horizontal dimension of F_2rc
! nk            I           Vertical dimension of F_2bc and F_2rc
!----------------------------------------------------------------
!
      integer err
!
!----------------------------------------------------------------------
!
      call RPN_COMM_coll (F_2rc,1,bni,1,bnj,bni,bnj,nk,0,0,1, &
                          F_2cc,Minx,Maxx,Miny,Maxy,0,0,err)
!
!----------------------------------------------------------------------
!
      return
      end

