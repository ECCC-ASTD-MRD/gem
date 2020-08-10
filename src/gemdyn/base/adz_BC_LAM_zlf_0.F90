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

!**s/r adz_BC_LAM_zlf_0 - Call setups to prepare Bermejo-Conde LAM with ZLF

      subroutine adz_BC_LAM_zlf_0 (F_nptr,F_number)

      use adz_mem
      use adz_options
      use gem_timing
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_nptr,F_number

      !object
      !===========================================================================
      !     Call setups to prepare Bermejo-Conde LAM with ZLF:
      !     F_number=0: Before High-order Interpolation 
      !     F_number=1: Before adz_post_tr
      !     F_number=2:  After adz_post_tr
      !     ----------------------------------------------------------------------
      !     NOTE: When Bermejo-Conde Flux ZLF is required for at least one tracer,
      !           advections over extended zone are done for all tracers
      !===========================================================================

      integer :: n
      real, dimension(1,1,1) :: empty
!
!---------------------------------------------------------------------
!
      call gemtime_start (85, 'ZLF_0______', 36)

      !------------------------------------------------------------------
      !Prepare Bermejo-Conde LAM with ZLF before High-order Interpolation 
      !------------------------------------------------------------------
      if (F_number==0) then

         do n=1,F_nptr

            !ZERO piloting conditions of SRC outside EXTENSION (CFL)
            !-------------------------------------------------------
            call adz_BC_LAM_zlf (Adz_stack(n)%src,empty,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,1)

            !Keep piloting conditions of DST outside CORE
            !--------------------------------------------
            call adz_BC_LAM_zlf (Adz_stack(n)%dst,Adz_stack(n)%pil,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,2)

         end do

      !-----------------------------------------------------
      !Prepare Bermejo-Conde LAM with ZLF before adz_post_tr
      !-----------------------------------------------------
      else if (F_number==1) then

         do n=1,F_nptr

            !ZERO piloting conditions of DST outside EXTENSION (BCS_BASE)
            !------------------------------------------------------------
            call adz_BC_LAM_zlf (Adz_stack(n)%dst,empty,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,3)

         end do

      !----------------------------------------------------
      !Prepare Bermejo-Conde LAM with ZLF after adz_post_tr
      !----------------------------------------------------
      else if (F_number==2) then

         do n=1,F_nptr

            !Bermejo-Conde LAM Flux ZLF: Reset piloting conditions of DST outside CORE
            !-------------------------------------------------------------------------
            call adz_BC_LAM_zlf (Adz_stack(n)%dst,Adz_stack(n)%pil,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,4)

         end do

      end if

      call gemtime_stop  (85)
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_BC_LAM_zlf_0
