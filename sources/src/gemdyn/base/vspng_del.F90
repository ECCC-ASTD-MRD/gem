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

!*s/r vspng_del2 - horizontal diffusion problem

      subroutine vspng_del2 ( F_sol, F_opsxp0_8, F_opsyp0_8  , &
                              F_aix_8,F_bix_8,F_cix_8,F_dix_8, &
                              F_aiy_8,F_biy_8,F_ciy_8        , &
                              Minx,Maxx,Miny,Maxy, nke, ntrp12, ntrp22, fnjb)
      use glb_ld
      use ldnh
      use trp
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy, nke, ntrp12, ntrp22, fnjb

      real   F_sol(Minx:Maxx,Miny:Maxy,*)
      real*8 F_opsxp0_8(*),F_opsyp0_8(*), &
             F_aix_8(nke,ntrp12,*),F_bix_8(nke,ntrp12,*), &
             F_cix_8(nke,ntrp12,*),F_dix_8(nke,ntrp12,*), &
             F_aiy_8(nke,ntrp22,*),F_biy_8(nke,ntrp22,*), &
             F_ciy_8(nke,ntrp22,*)
!
      integer i, j, k, cnt
      real*8 w1_8(ldnh_maxx,nke,ldnh_maxy),  &
             w2_8(ldnh_maxy,nke,ldnh_maxx),  &
             t1_8(nke,Trp_12emax,G_ni+Ptopo_npex), &
             t2_8(nke,Trp_22emax,G_nj+Ptopo_npey), &
             g1(nke*Trp_12emax,G_ni),ax(nke*Trp_12emax,G_ni), &
             cx(nke*Trp_12emax,G_ni),g2(nke*Trp_22emax,G_nj), &
             ay(nke*Trp_22emax,G_nj),cy(nke*Trp_22emax,G_nj)
!*
!     __________________________________________________________________
!
      do j = 1, l_nj
      do k = 1, nke
      do i = 1, l_ni
         w1_8(i,k,j) = dble(F_sol(i,j,k))
      end do
      end do
      end do
!
      call rpn_comm_transpose ( w1_8 , 1, ldnh_maxx, G_ni, nke, &
                                1, Trp_12emax, ldnh_maxy, t1_8,  1, 2 )
!
      cnt = 0
      do j = 1, Trp_12en
      do k = 1, nke
         cnt = cnt + 1
         do i = 1, G_ni-1
            g1(cnt,i) = F_bix_8(k,j,i)*F_opsxp0_8(i)*t1_8(k,j,i)
            ax(cnt,i) = F_aix_8(k,j,i)
            cx(cnt,i) = F_cix_8(k,j,i)
         end do
         g1(cnt,G_ni) = F_opsxp0_8(G_ni)*t1_8(k,j,G_ni)
      end do
      end do
!
      do i = 2, G_ni-1
      do k = 1, cnt
         g1(k,i) = g1(k,i) - ax(k,i)*g1(k,i-1)
      end do
      end do
      do i = G_ni-2, 1, -1
      do k = 1, cnt
         g1(k,i) = g1(k,i) - cx(k,i)*g1(k,i+1)
      end do
      end do
!
      cnt = 0
      do j = 1, Trp_12en
!VDIR NOVECTOR
      do k = 1, nke
         cnt = cnt + 1
         t1_8(k,j,G_ni) =    F_bix_8(k,j,G_ni)*g1(cnt,G_ni  ) &
                           + F_cix_8(k,j,G_ni)*g1(cnt,1     ) &
                           + F_aix_8(k,j,1   )*g1(cnt,G_ni-1)
         do i = 1, G_ni - 1
            t1_8(k,j,i) = g1(cnt,i) + F_dix_8(k,j,i)*t1_8(k,j,G_ni)
         end do
      end do
      end do
!
      call rpn_comm_transpose ( w1_8 , 1, ldnh_maxx, G_ni, nke, &
                                1, Trp_12emax, ldnh_maxy, t1_8,  -1, 2 )
!
      do j = 1, l_nj
      do k = 1, nke
      do i = 1, l_ni
         w2_8(j,k,i) = w1_8(i,k,j)
      end do
      end do
      end do
!
      call rpn_comm_transpose ( w2_8 , 1, ldnh_maxy, G_nj, nke, &
                                1, Trp_22emax, ldnh_maxx, t2_8,  2, 2 )
!
      cnt = 0
      do i = 1, Trp_22en
      do k = 1, nke
         cnt = cnt + 1
         do j = 1, fnjb
            g2 (cnt,j) = F_biy_8(k,i,j)*F_opsyp0_8(j)*t2_8(k,i,j)
            ay (cnt,j) = F_aiy_8(k,i,j)
            cy (cnt,j) = F_ciy_8(k,i,j)
         end do
      end do
      end do
!
      do j = 2, fnjb
      do k = 1, cnt
         g2 (k,j) = g2(k,j) - ay(k,j)*g2(k,j-1)
      end do
      end do
      do j = fnjb-1, 1, -1
      do k = 1, cnt
         g2 (k,j) = g2(k,j) - cy(k,j)*g2(k,j+1)
      end do
      end do
!
      cnt = 0
      do i = 1, Trp_22en
      do k = 1, nke
         cnt = cnt + 1
         do j = 1, fnjb
            t2_8(k,i,j)= g2(cnt,j)
         end do
      end do
      end do
!
      call rpn_comm_transpose ( w2_8 , 1, ldnh_maxy, G_nj, nke, &
                                1, Trp_22emax, ldnh_maxx, t2_8,  -2, 2 )
!
      do k = 1, nke
      do j = 1, l_nj
      do i = 1, l_ni
         F_sol(i,j,k) = w2_8(j,k,i)
      end do
      end do
      end do
!
!     __________________________________________________________________
!
      return
      end
