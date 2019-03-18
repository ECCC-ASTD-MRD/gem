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

!**s/r sol_mxma - parallel direct solution of an elliptic problem using MXMA

      subroutine  sol_mxma ( Sol, Rhs, Xevec                   , &
                             F_t0nis, F_t0njs, F_t0nj, F_nini  , &
                             F_t1nks, F_t1nk, F_t2nis, F_t2ni  , &
                             F_gni, F_gnj, F_gnk, F_nk         , &
                             F_npex1, F_npey1                  , &
                             F_ai,F_bi,F_ci, F_dg1,F_dg2,F_dwfft )
      use glb_ld
      use glb_pil
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer F_t0nis, F_t0njs, F_t0nj, F_t2nis, F_t2ni, F_nini
      integer F_t1nks, F_gnk,   F_nk  , F_t1nk , F_gni , F_gnj
      integer F_npex1, F_npey1

      real*8  Sol (1:F_t0nis, 1:F_t0njs, F_gnk), &
              Rhs (1:F_t0nis, 1:F_t0njs, F_gnk), &
              F_ai(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_bi(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_ci(1:F_t1nks, 1:F_t2nis, F_gnj), Xevec(F_nini)
      real*8  F_dwfft(1:F_t0njs, 1:F_t1nks, F_gni+2+F_npex1)
      real*8  F_dg1  (1:F_t0njs, 1:F_t1nks, F_gni+  F_npex1)
      real*8  F_dg2  (1:F_t1nks, 1:F_t2nis, F_gnj  +F_npey1)

!author    Abdessamad Qaddouri- July 1999

!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! Sol          O    - result of solver
! Rhs          I    - r.h.s. of elliptic equation
! Xevec        I    - eigenvectors
! F_t0nis      I    - maximum index on X for Rhs,Sol
! F_t0njs      I    - maximum index on Y for Rhs,Sol
! F_t0nj       I    - number of points on local PEy for J (ldnh_nj)
! F_t1nks      I    - maximum index on local PEx for K (trp_12smax)
! F_gnk        I    - G_nk-1 points in z direction globally (Schm_nith)
! F_t1nk       I    - number of points on local PEx for K (trp_12sn)
! F_gni        I    - number of points in x direction globally (G_ni)
! F_gnj        I    - number of points in y direction globally (G_nj)
! F_t2nis      I    - maximum index on local PEy for I (trp_22max)
! F_t2ni       I    - number of points on local PEy for I (trp_22n)
! F_t0nis1     I    - maximum index on local PEx for K (trp_12smax)
! F_t0nis2     I    - maximum index on local PEy for I (trp_22max)
! F_gnj        I    - number of points along J globally (G_nj)
! F_npex1      I    - number of processors on X
! F_npey1      I    - number of processors on Y
! F_ai         I    - sub   diagonal of LU factorization
! F_bi         I    -       diagonal of LU factorization
! F_ci         I    - super diagonal of LU factorization
! F_dg1        I    - work field
! F_dg2        I    - work field
! F_dwfft      I    - work field


      integer i,j,k,jr,l_pil_w,l_pil_e
      integer piece, p0, pn, ptotal, plon
      real*8, parameter :: zero= 0.d0, one= 1.d0
!     __________________________________________________________________
!
      l_pil_w=0
      l_pil_e=0
      if (l_south) l_pil_w= Lam_pil_w
      if (l_north) l_pil_e= Lam_pil_e
!
      call rpn_comm_transpose( Rhs, 1, F_t0nis, F_gni, (F_t0njs-1+1),&
                                       1, F_t1nks, F_gnk, F_dg1, 1,2 )

!$omp parallel private(i,j,k,p0,pn,piece,jr) &
!$omp          shared(ptotal,plon,l_pil_w,l_pil_e)

!$omp do
      do i= 1,F_gni
         F_dg1(F_t0nj+1-pil_n:F_t0njs,1:F_t1nk,i)= zero
         F_dg1(             1:pil_s  ,1:F_t1nk,i)= zero
         F_dwfft(1:F_t0njs,F_t1nk+1:F_t1nks,i)   = zero
      end do
!$omp enddo

!     projection ( wfft = x transposed * g )

!$omp do
      do k=1,F_nk
         call dgemm('N','N', (F_t0njs-1+1-pil_s-pil_n),  &
                  F_gni-Lam_pil_w-Lam_pil_e, &
                  F_gni-Lam_pil_w-Lam_pil_e, &
                  1._8, F_dg1(1+pil_s,k,1+Lam_pil_w),  &
                  (F_t0njs-1+1)* (F_t1nks-1+1), xevec, &
                  F_gni-Lam_pil_w-Lam_pil_e, &
                  0._8, F_dwfft(1+pil_s,k,1+Lam_pil_w), &
                  (F_t0njs-1+1)* (F_t1nks-1+1))
      end do
!$omp enddo

!$omp single
      call rpn_comm_transpose &
           ( F_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                        1, F_t2nis, F_gni, F_dg2, 2, 2 )
!$omp end single

      ptotal = F_t2ni-l_pil_e-l_pil_w
      plon   = (ptotal+Ptopo_npeOpenMP) / Ptopo_npeOpenMP

!$omp do
       do piece=1,Ptopo_npeOpenMP
          p0 = 1+l_pil_w + plon*(piece-1)
          pn = min(F_t2ni-l_pil_e,plon*piece+l_pil_w)

          j =1+Lam_pil_s
          do i= p0, pn
          do k= 1 , F_nk
                F_dg2(k,i,j) = F_bi(k,i,j)*F_dg2(k,i,j)
          end do
          end do

          do j =2+Lam_pil_s, F_gnj-Lam_pil_n
             jr =  j - 1
             do i= p0, pn
             do k= 1 , F_nk
                F_dg2(k,i,j) = F_bi(k,i,j)* F_dg2(k,i,j) - F_ai(k,i,j) &
                                          * F_dg2(k,i,jr)
             end do
             end do
          end do

         do j = F_gnj-1-Lam_pil_n, 1+Lam_pil_s, -1
            jr =  j + 1
            do i= p0, pn
            do k= 1 , F_nk
              F_dg2(k,i,j) = F_dg2(k,i,j) - F_ci(k,i,j) * F_dg2(k,i,jr)
            end do
            end do
         end do
      end do
!$omp enddo

!$omp single
      call rpn_comm_transpose &
           ( F_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1), &
                     1, F_t2nis, F_gni, F_dg2,- 2, 2 )
!$omp end single

!     inverse projection ( r = x * w )

!$omp do
      do k=1,F_nk
         call dgemm('N','T', (F_t0njs-1+1-pil_s-pil_n),  &
                  F_gni-Lam_pil_w-Lam_pil_e, &
                  F_gni-Lam_pil_w-Lam_pil_e, &
                  1._8, F_dwfft(1+pil_s,k,1+Lam_pil_w), &
                  (F_t0njs-1+1) * (F_t1nks-1+1),xevec , &
                  F_gni-Lam_pil_w-Lam_pil_e, &
                  0._8, F_dg1(1+pil_s,k,1+Lam_pil_w), &
                  (F_t0njs-1+1) * (F_t1nks-1+1))
      end do
!$omp end do

!$omp end parallel

      call rpn_comm_transpose( Sol, 1, F_t0nis, F_gni, (F_t0njs-1+1),&
                                    1, F_t1nks, F_gnk, F_dg1, -1, 2 )

!     __________________________________________________________________
!
      return
      end
