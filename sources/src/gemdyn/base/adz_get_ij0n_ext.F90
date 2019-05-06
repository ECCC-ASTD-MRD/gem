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

!**s/r adz_get_ij0n_ext: Establish scope of extended advection computational i,j domain

      subroutine adz_get_ij0n_ext (F_i0_e,F_in_e,F_j0_e,F_jn_e,F_zone)

      use glb_ld
      use HORgrid_options

      implicit none

#include <arch_specific.hf>

      integer :: F_i0_e,F_j0_e,F_in_e,F_jn_e,F_zone

      !object
      !===================================================================
      !     Establish scope of extended advection computational i,j domain
      !===================================================================

      integer :: jext,BCS_BASE

      if (F_zone==0) then !CORE

         F_i0_e = 1+pil_w
         F_in_e = l_ni-pil_e
         F_j0_e = 1+pil_s
         F_jn_e = l_nj-pil_n

      else if (F_zone==1) then !EXTENSION (CFL)

         F_i0_e = 1
         F_in_e = l_ni
         F_j0_e = 1
         F_jn_e = l_nj

         jext = Grd_maxcfl + 1

         if (l_west)  F_i0_e =    1 + pil_w - jext
         if (l_east)  F_in_e = l_ni - pil_e + jext
         if (l_south) F_j0_e =    1 + pil_s - jext
         if (l_north) F_jn_e = l_nj - pil_n + jext

      else if (F_zone==2) then !EXTENSION (BCS_BASE)

         F_i0_e = 1
         F_in_e = l_ni
         F_j0_e = 1
         F_jn_e = l_nj

         BCS_BASE= 4

         if (l_west)  F_i0_e = 1+BCS_BASE
         if (l_east)  F_in_e = l_ni-BCS_BASE
         if (l_south) F_j0_e = 1+BCS_BASE
         if (l_north) F_jn_e = l_nj-BCS_BASE

      else

         call handle_error (-1,'ADZ_GET_IJ0N_EXT','ZONE not defined')

      end if

      return

      end subroutine adz_get_ij0n_ext
