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

!**s/r adz_BC_LAM_zlf - Various setups to prepare Bermejo-Conde LAM with ZLF

      subroutine adz_BC_LAM_zlf_hlt (F1,F2,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,F_setup)

      use adz_mem
      use HORgrid_options

      use, intrinsic :: iso_fortran_env
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

      integer :: i,j,k,i0_e,j0_e,in_e,jn_e,i0_c,j0_c,in_c,jn_c,i0_b,j0_b,in_b,jn_b,k0,k0d2,ext,BCS_BASE
!
!---------------------------------------------------------------------
!
      k0 = Adz_k0t

      k0d2 = nint(k0/2.)

      !In F1: Set ZERO piloting conditions outside [i0_e:in_e]x[j0_e:jn_e]x[k0/2:nk]
      !-----------------------------------------------------------------------------
      if (F_setup==1) then

         ext = Grd_maxcfl + 1

         i0_e =    1 + pil_w - ext*west
         in_e = l_ni - pil_e + ext*east
         j0_e =    1 + pil_s - ext*south
         jn_e = l_nj - pil_n + ext*north

         if (l_west) then
!$omp do collapse(2)
           do i=1,i0_e-1
              do j=1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)=0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

         if (l_east) then
!$omp do collapse(2)
           do i=in_e+1,F_ni
              do j=1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)=0.
                 enddo
              enddo
           enddo 
!$omp enddo 
         endif

         if (l_south) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=1,j0_e-1
                 do k=1,F_nk 
                    F1(i,j,k)=0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

         if (l_north) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=jn_e+1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)=0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

!$omp do collapse(2)
         do i=i0_e,in_e
            do j=j0_e,jn_e
               do k=1,k0d2-1
                  F1(i,j,k)=0.
               enddo
            enddo
         enddo
!$omp enddo 

      !In F2: Keep piloting conditions outside [i0_c:in_c]x[j0_c:jn_c]x[k0:nk] of F1
      !-----------------------------------------------------------------------------
      else if (F_setup==2) then

         i0_c = 1    + pil_w
         in_c = l_ni - pil_e
         j0_c = 1    + pil_s
         jn_c = l_nj - pil_n

         if (l_west) then
!$omp do collapse(2)
           do i=1,i0_c-1
              do j=1,F_nj
                 do k=1,F_nk 
                    F2(i,j,k)= F1(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

         if (l_east) then
!$omp do collapse(2)
           do i=in_c+1,F_ni
              do j=1,F_nj
                 do k=1,F_nk 
                    F2(i,j,k)= F1(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo
         endif

         if (l_south) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=1,j0_c-1
                 do k=1,F_nk 
                    F2(i,j,k)= F1(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo
         endif

         if (l_north) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=jn_c+1,F_nj
                 do k=1,F_nk 
                    F2(i,j,k)= F1(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

!$omp do collapse(2)
         do i=i0_c,in_c
            do j=j0_c,jn_c
               do k=1,k0-1
                  F2(i,j,k)= F1(i,j,k)
               enddo
            enddo
         enddo
!$omp enddo 

      !In F1: Set ZERO piloting conditions outside EXTENSION (BCS_BASE)
      !----------------------------------------------------------------
      else if (F_setup==3) then

         BCS_BASE = 4

         i0_b =    1 + BCS_BASE*west
         in_b = l_ni - BCS_BASE*east
         j0_b =    1 + BCS_BASE*south
         jn_b = l_nj - BCS_BASE*north

         if (l_west) then
!$omp do collapse(2)
           do i=1,i0_b-1
              do j=1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)= 0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

         if (l_east) then
!$omp do collapse(2)
           do i=in_b+1,F_ni
              do j=1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)= 0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

         if (l_south) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=1,j0_b-1
                 do k=1,F_nk 
                    F1(i,j,k)= 0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

         if (l_north) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=jn_b+1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)= 0.
                 enddo
              enddo
           enddo
!$omp enddo 
         endif

      !In F1: Reset piloting conditions outside [i0_c:in_c]x[j0_c:jn_c]x[k0:nk] stored in F2
      !-------------------------------------------------------------------------------------
      else if (F_setup==4) then

         i0_c = 1    + pil_w
         in_c = l_ni - pil_e
         j0_c = 1    + pil_s
         jn_c = l_nj - pil_n

         if (l_west) then
!$omp do collapse(2)
           do i=1,i0_c-1
              do j=1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)= F2(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo
         endif

         if (l_east) then
!$omp do collapse(2)
           do i=in_c+1,F_ni
              do j=1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)= F2(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo
         endif

         if (l_south) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=1,j0_c-1
                 do k=1,F_nk 
                    F1(i,j,k)= F2(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo
         endif

         if (l_north) then
!$omp do collapse(2)
           do i=1,F_ni
              do j=jn_c+1,F_nj
                 do k=1,F_nk 
                    F1(i,j,k)= F2(i,j,k)
                 enddo
              enddo
           enddo
!$omp enddo
         endif

!$omp do collapse(2)
         do i=i0_c,in_c
            do j=j0_c,jn_c
               do k=1,k0-1
                  F1(i,j,k)= F2(i,j,k)
               enddo
            enddo
         enddo
!$omp enddo 

      else

         call gem_error (-1,'ADZ_BC_LAM_ZLF','SETUP not defined')

      end if
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_BC_LAM_zlf_hlt
