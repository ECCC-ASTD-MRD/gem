!---------------------------------- LICENCE BEGIN -------------------------------

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

!**s/r set_io_pes

      integer function set_io_pes ( F_comm_id, F_comm_setno, F_iome, &
                                    F_comm_io, F_iobcast, F_npes )
      use iso_c_binding
      use glb_ld
      use lun
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer F_comm_id, F_comm_setno, F_iome, F_comm_io, &
              F_iobcast, F_npes

      include "rpn_comm.inc"

      integer pe_xcoord(F_npes), pe_ycoord(F_npes), err
!
!-------------------------------------------------------------------
!
      set_io_pes = -1
      err= RPN_COMM_io_pe_valid_set (pe_xcoord,pe_ycoord,F_npes,&
                             Ptopo_npex,Ptopo_npey,Lun_debug_L,0)
      if (err < 0) return

      F_comm_id    = RPN_COMM_create_2dgrid (  G_ni, G_nj, &
                              l_minx, l_maxx, l_miny, l_maxy )
      F_comm_setno = RPN_COMM_create_io_set ( F_npes, 0 )
      F_iome       = RPN_COMM_is_io_pe      ( F_comm_setno )
      F_comm_io    = RPN_COMM_io_pe_comm    ( F_comm_setno )
      F_iobcast    = RPN_COMM_io_pe_gridid  ( F_comm_setno, 0 )

      set_io_pes= 0
!
!-------------------------------------------------------------------
!
      return
      end
