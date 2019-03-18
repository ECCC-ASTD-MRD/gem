
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

!**s/r adv_BC_LAM_zlf - Various setups to prepare Bermejo-Conde LAM with ZLF

      subroutine adv_BC_LAM_zlf (F1,F2,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,F_setup)

      use glb_ld

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy
      integer, intent(in) :: F_ni,F_nj,F_nk
      integer, intent(in) :: F_setup
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(inout) :: F1,F2

      !object
      !=========================================================
      !     Various setups to prepare Bermejo-Conde LAM with ZLF
      !=========================================================

      integer i0_2,j0_2,in_2,jn_2,i0_0,j0_0,in_0,jn_0,i0_3,j0_3,in_3,jn_3

      !In F1: Set ZERO piloting conditions outside EXTENSION (CFL)
      !-----------------------------------------------------------
      if (F_setup==1) then

         call adv_get_ij0n_ext (i0_2,in_2,j0_2,jn_2,2)

         if (l_west)  F1(1     :i0_2-1,1     :F_nj,  1:F_nk) = 0.
         if (l_east)  F1(in_2+1:F_ni,  1     :F_nj,  1:F_nk) = 0.
         if (l_south) F1(1     :F_ni,  1     :j0_2-1,1:F_nk) = 0.
         if (l_north) F1(1     :F_ni,  jn_2+1:F_nj,  1:F_nk) = 0.

      !In F2: Keep piloting conditions outside CORE of F1
      !--------------------------------------------------
      else if (F_setup==2) then

         call adv_get_ij0n_ext (i0_0,in_0,j0_0,jn_0,0)

         if (l_west)  F2(1     :i0_0-1,1     :F_nj,  1:F_nk) = F1(1     :i0_0-1,1     :F_nj,  1:F_nk)
         if (l_east)  F2(in_0+1:F_ni,  1     :F_nj,  1:F_nk) = F1(in_0+1:F_ni,  1     :F_nj,  1:F_nk)
         if (l_south) F2(1     :F_ni,  1     :j0_0-1,1:F_nk) = F1(1     :F_ni,  1     :j0_0-1,1:F_nk)
         if (l_north) F2(1     :F_ni,  jn_0+1:F_nj,  1:F_nk) = F1(1     :F_ni,  jn_0+1:F_nj,  1:F_nk)

      !In F1: Set ZERO piloting conditions outside EXTENSION (BCS_BASE)
      !----------------------------------------------------------------
      else if (F_setup==3) then

         call adv_get_ij0n_ext (i0_3,in_3,j0_3,jn_3,3)

         if (l_west)  F1(1     :i0_3-1,1     :F_nj,  1:F_nk) = 0.
         if (l_east)  F1(in_3+1:F_ni,  1     :F_nj,  1:F_nk) = 0.
         if (l_south) F1(1     :F_ni,  1     :j0_3-1,1:F_nk) = 0.
         if (l_north) F1(1     :F_ni,  jn_3+1:F_nj,  1:F_nk) = 0.

      !In F1: Reset piloting conditions outside CORE stored in F2
      !----------------------------------------------------------
      else if (F_setup==4) then

         call adv_get_ij0n_ext (i0_0,in_0,j0_0,jn_0,0)

         if (l_west)  F1(1     :i0_0-1,1     :F_nj,  1:F_nk) = F2(1     :i0_0-1,1     :F_nj,  1:F_nk)
         if (l_east)  F1(in_0+1:F_ni,  1     :F_nj,  1:F_nk) = F2(in_0+1:F_ni,  1     :F_nj,  1:F_nk)
         if (l_south) F1(1     :F_ni,  1     :j0_0-1,1:F_nk) = F2(1     :F_ni,  1     :j0_0-1,1:F_nk)
         if (l_north) F1(1     :F_ni,  jn_0+1:F_nj,  1:F_nk) = F2(1     :F_ni,  jn_0+1:F_nj,  1:F_nk)

      else

         call handle_error (-1,'ADV_BC_LAM_ZLF','SETUP not defined')

      end if

      return

      end subroutine adv_BC_LAM_zlf
