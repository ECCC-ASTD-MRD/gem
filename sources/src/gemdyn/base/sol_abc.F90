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

!**s/r sol_abc - initialisation de l'operateur unidimensionnel.
!                direction de travail: suivant la latitude.
!
      subroutine sol_abc ( F_hcon_8, F_yg_8, F_opsyp0_8, F_opsyp2_8, F_xeval_8, &
                           Minz, Maxz, F_nk, F_k0, &
                           Minx, Maxx, F_ni, F_i0, &
                           GNI, GNJ, GNK, F_ai_8, F_bi_8, F_ci_8 )
      use glb_ld
      use glb_pil
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer F_k0, F_i0, Minz, Maxz, Minx, Maxx, F_nk, F_ni, GNI,GNJ,GNK
      real(kind=REAL64)  F_hcon_8(GNK), F_yg_8(GNJ), F_opsyp0_8(3*Gnj), &
              F_opsyp2_8(3*Gnj), F_xeval_8(Gni), &
              F_ai_8(Minz:Maxz,Minx:Maxx,GNJ)  , &
              F_bi_8(Minz:Maxz,Minx:Maxx,GNJ)  , &
              F_ci_8(Minz:Maxz,Minx:Maxx,GNJ)
!
!Author
!     Abdessamad Qaddouri- JULY 1999
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_hcon_8     I    - vertical projector
! F_yg_8       I    - latitudes  of the scalar grid in radians
! F_opsyp0_8   I    - north-south projection operators
! F_opsyp2_8   I    - north-south projection operators
! F_xeval_8    I    - horizontal eigenvalues
! F_k0         I    - the global starting index for K on PEx (Trp_12sn0)
! F_i0         I    - the global starting index for I on PEy (Trp_22n0)
! Minx1        I    - minimum index on local PEx for K (trp_12smin)
! Maxx1        I    - maximum index on local PEx for K (trp_12smax)
! Minx2        I    - minimum index on local PEy for I (trp_22min)
! Maxx2        I    - maximum index on local PEy for I (trp_22max)
! F_nk         I    - number of points on local PEx for K (Trp_12sn)
! F_ni         I    - number of points on local PEy for I (Trp_22n)
! GNJ          I    - number of points on J globally (G_nj)
! F_ai_8       O    - sub-   diagonal of LU factorization
! F_bi_8       O    -        diagonal of LU factorization
! F_ci_8       O    - super- diagonal of LU factorization


      integer k, kk, i, ii, j
      integer l_pil_w,l_pil_e
      real(kind=REAL64), parameter :: zero = 0.0, one  = 1.0
      real(kind=REAL64)  a_8(F_nk,Minx:Maxx,GNJ), b_8(F_nk,Minx:Maxx,GNJ), &
              c_8(F_nk,Minx:Maxx,GNJ), di_8
!
!     ---------------------------------------------------------------
!
!  The I vector lies on the Y processor so, l_pil_w and l_pil_e will
!  represent the pilot region along I
!
      l_pil_w=0
      l_pil_e=0
      if (l_south) l_pil_w= Lam_pil_w
      if (l_north) l_pil_e= Lam_pil_e

         a_8= zero ;    b_8= one  ;    c_8= zero
      F_ai_8= zero ; F_bi_8= zero ; F_ci_8= zero

!  Preparing matrices: a_8 b_8 c_8

      do k = 1, F_nk
         kk = F_k0 + k - 1
         do i = 1+l_pil_w, F_ni-l_pil_e
            ii = i + F_i0 - 1
            do j=2+Lam_pil_s, GNJ-Lam_pil_n
               di_8 = F_opsyp0_8(GNJ+j) / cos( F_yg_8 (j) )**2
               b_8(k,i,j)= F_xeval_8(ii) * di_8 + &
                           F_hcon_8(kk)*F_opsyp0_8(GNJ+j    )+F_opsyp2_8( GNJ+j   )
               c_8(k,i,j)= F_hcon_8(kk)*F_opsyp0_8(2*GNJ+j  )+F_opsyp2_8(2*GNJ+j  )
               a_8(k,i,j)= F_hcon_8(kk)*F_opsyp0_8(2*GNJ+j-1)+F_opsyp2_8(2*GNJ+j-1)
            end do
            a_8(k,i,1+Lam_pil_s  )= zero
            c_8(k,i,GNJ-Lam_pil_n)= zero
            di_8 = F_opsyp0_8(GNJ+1+Lam_pil_s) / cos( F_yg_8 (1+Lam_pil_s) )**2
            b_8(k,i,1+Lam_pil_s  )= F_xeval_8(ii) * di_8 + &
                                    F_hcon_8(kk)*F_opsyp0_8(  GNJ+1+Lam_pil_s) + &
                                    F_opsyp2_8(  GNJ+1+Lam_pil_s)
            c_8(k,i,1+Lam_pil_s  )= F_hcon_8(kk)*F_opsyp0_8(2*GNJ+1+Lam_pil_s) + &
                                    F_opsyp2_8(2*GNJ+1+Lam_pil_s)
         end do
      end do

!  Compute matrices: F_ai_8 F_bi_8 F_ci_8

      j= 1+Lam_pil_s
      do k= 1, F_nk
      do i= 1+l_pil_w, F_ni-l_pil_e
         F_bi_8(k,i,j) = one / b_8(k,i,j)
         F_ci_8(k,i,j) = c_8(k,i,j) * F_bi_8(k,i,j)
      end do
      end do

      do j= 2+Lam_pil_s, GNJ-Lam_pil_n-1
      do k= 1, F_nk
      do i= 1+l_pil_w, F_ni-l_pil_e
         F_bi_8(k,i,j) = one / ( b_8(k,i,j)-a_8(k,i,j) * F_ci_8(k,i,j-1) )
         F_ci_8(k,i,j) = F_bi_8(k,i,j) * c_8(k,i,j)
         F_ai_8(k,i,j) = F_bi_8(k,i,j) * a_8(k,i,j)
      end do
      end do
      end do

      j= GNJ-Lam_pil_n
      do k= 1, F_nk
      do i= 1+l_pil_w, F_ni-l_pil_e
         F_bi_8(k,i,j) = one / ( b_8(k,i,j)-a_8(k,i,j) * F_ci_8(k,i,j-1) )
         F_ai_8(k,i,j) = F_bi_8(k,i,j) * a_8(k,i,j)
      end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end

