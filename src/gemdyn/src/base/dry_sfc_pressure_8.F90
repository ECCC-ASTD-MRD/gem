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

!**s/r dry_sfc_pressure_8 - Compute dry air surface pressure REAL64
!
      subroutine dry_sfc_pressure_8 (F_drysfcp0_8, presT_8, p0T_8, &
                               Minx,Maxx,Miny,Maxy,Nk,F_k0,F_timelevel_S)
      use glb_ld
      use tr3d
      use mem_tracers
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=1) :: F_timelevel_S
      integer Minx,Maxx,Miny,Maxy,Nk,F_k0
      real(kind=REAL64) F_drysfcp0_8(Minx:Maxx,Miny:Maxy),&
                             presT_8(Minx:Maxx,Miny:Maxy,Nk),&
                               p0T_8(Minx:Maxx,Miny:Maxy)
      integer i,j,k
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: sumq
      real, pointer, dimension(:,:,:)         :: tr
!     ________________________________________________________________
!
      if (F_timelevel_S == 'P') then
         call sumhydro (sumq,Minx,Maxx,Miny,Maxy,Nk,Tr3d_ntr, trt1)
         tr=>tracers_P(Tr3d_hu)%pntr
      endif
      if (F_timelevel_S == 'M') then
         call sumhydro (sumq,Minx,Maxx,Miny,Maxy,Nk,Tr3d_ntr, trt0)
         tr=>tracers_M(Tr3d_hu)%pntr
      endif

      do k=F_k0,Nk
         sumq(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,k)= &
         sumq(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,k)+ &
         tr  (1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,k)
      end do

      do j=1+pil_s,l_nj-pil_n
         F_drysfcp0_8(:,j) = 0.0d0
         do k=F_k0,Nk-1
            do i=1+pil_w,l_ni-pil_e
               F_drysfcp0_8(i,j)= F_drysfcp0_8(i,j) + &
                    (1.-sumq(i,j,k))*(presT_8(i,j,k+1) - presT_8(i,j,k))
            end do
         end do
         do i=1+pil_w,l_ni-pil_e
            F_drysfcp0_8(i,j)= F_drysfcp0_8(i,j) + &
                 (1.-sumq(i,j,Nk))*(p0T_8(i,j) - presT_8(i,j,Nk))
         end do
      end do

!     ________________________________________________________________
!
      return
      end
