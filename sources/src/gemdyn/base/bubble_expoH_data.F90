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
!**s/r bubble_expoH_data - generates initial condition for Robert's bubble
!                    experiment (Robert 1993 JAS) - EXPO height coord.
!
      subroutine bubble_expoH_data ( F_u, F_v, F_t, F_s, F_q, F_topo, F_sls,&
                               Mminx,Mmaxx,Mminy,Mmaxy,nk )
      use glb_ld
      use bubble_options
      use dyn_fisl_options
      use tdpack
      use type_mod
      use ver
      implicit none
#include <arch_specific.hf>

      integer Mminx,Mmaxx,Mminy,Mmaxy,nk
      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_t    (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_sls  (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_q    (Mminx:Mmaxx,Mminy:Mmaxy,nk+1)

      integer :: i,j,k,ii,err
      real*8 :: theta, r,rad
!
!     ---------------------------------------------------------------
!
      F_topo(:,:) = 0.0
      F_sls (:,:) = 0.0
      F_s   (:,:) = 0.0
      F_u (:,:,:) = 0.0
      F_v (:,:,:) = 0.0
      F_q (:,:,:) = 0.0
!
!---------------------------------------------------------------------
!     Initialize temperature
!---------------------------------------------------------------------
!
      if(.not.bubble_gaus_L) then

         do k=1,g_nk
            do j=1,l_nj
            do i=1,l_ni
               ii=i+l_i0-1
               theta = bubble_theta
               if ( (((ii)-bubble_ictr)**2 +((k)-bubble_kctr)**2) < bubble_rad**2 ) then
                   theta = theta + 0.5d0
               end if
               F_t(i,j,k) = theta
            end do
            end do
         end do

      else

      ! Gaussian

         do k=1,g_nk
            do j=1,l_nj
            do i=1,l_ni
               ii=i+l_i0-1
               theta = bubble_theta
               r = sqrt( ( ((ii)-bubble_ictr) * dble(bubble_dx) )**2 + ( ((k)-bubble_kctr) * dble(bubble_dz) )**2)
               rad = bubble_rad*dble(bubble_dx)
               if ( r <= rad ) then
                   theta = bubble_theta + 0.5d0
               else
                   theta = bubble_theta + 0.5d0 * exp( -(r - rad)**2 / 100.d0**2 );
               end if
               F_t(i,j,k) = theta
            end do
            end do
         end do

      end if

!
!---------------------------------------------------------------------
!     Initialize (horizontally uniform) Exner pressure
!---------------------------------------------------------------------
!
      do k=1,g_nk+1
         do j=1,l_nj
         do i=1,l_ni
               F_q(i,j,k) = 1.d0 - grav_8 / (cpd_8 * bubble_theta) * Ver_z_8%m(k)
         end do
         end do
      end do
!
      err=0
      if ((l_north).and.(l_nj-2*pil_n+1<1)) err=-1
      if ((l_east ).and.(l_ni-2*pil_e+1<1)) err=-1
      if ((l_south).and.(2*pil_s>l_nj)    ) err=-1
      if ((l_west ).and.(2*pil_w>l_ni)    ) err=-1
      call gem_error(err,'ABORT in BUBBLE_DATA',&
          'Partitionning NOT allowed for MIRROR')

 9000 format(/,'CREATING INPUT DATA FOR BUBBLE expo (height)THEORETICAL CASE' &
            /,'====================================-======================')
!
!     -----------------------------------------------------------------
!
      return
      end subroutine bubble_expoH_data
