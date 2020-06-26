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
!**s/r sol_fft_lam - parallel direct solution of an elliptic problem
!                    for LAM grids using FFT

      subroutine sol_fft_lam ( F_sol, Rhs                      , &
                               F_t0nis, F_t0njs, F_t0nj        , &
                               F_t1nks, F_t1nk, F_t2nis, F_t2ni, &
                               F_gni, F_gnj, F_gnk, F_nk       , &
                               F_ai, F_bi, F_ci )
      use iso_c_binding
      use glb_ld
      use glb_pil
      use ptopo
      use sol
      use adz_mem
      use dyn_fisl_options
      use gem_timing
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer F_t0nis, F_t0njs, F_t0nj, F_t2nis, F_t2ni
      integer F_t1nks, F_gnk,   F_nk  , F_t1nk , F_gni, F_gnj

      real(kind=REAL64)  F_Sol (1:F_t0nis, 1:F_t0njs, F_gnk ), &
              Rhs (1:F_t0nis, 1:F_t0njs, F_gnk ), &
              F_ai(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_bi(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_ci(1:F_t1nks, 1:F_t2nis, F_gnj)

!author    Abdessamad Qaddouri- JULY 1999
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! Sol          O    - result of solver
! Rhs          I    - r.h.s. of elliptic equation
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
! F_ai         I    - sub   diagonal of LU factorization
! F_bi         I    -       diagonal of LU factorization
! F_ci         I    - super diagonal of LU factorization

      integer i, j, k, jr, n
      integer piece, p0, pn, plon, ptotal
      real(kind=REAL64), parameter :: zero = 0.d0
!     __________________________________________________________________
!
      n= (Adz_icn-1)*Schm_itnlh + Adz_itnl
      call rpn_comm_transpose ( Rhs, 1, F_t0nis, F_gni, (F_t0njs-1+1), &
                                1, F_t1nks, F_gnk, Sol_dwfft, 1, 2 )

      call time_trace_barr(gem_time_trace, 11100+n, Gem_trace_barr, &
                                      Ptopo_intracomm, MPI_BARRIER)
      ! Use OpenMP, if configured, to zero out unused portions of the
      ! transpose array
!$omp parallel private(i,j,k,jr,p0,pn,piece) &
!$omp          shared(plon,ptotal)
!$omp do
      do i= 1,F_gni
         Sol_dwfft(F_t0nj+1-pil_n:F_t0njs,        1:F_t1nk ,i)= zero
         Sol_dwfft(             1:pil_s  ,        1:F_t1nk ,i)= zero
         Sol_dwfft(             1:F_t0njs, F_t1nk+1:F_t1nks,i)= zero
      end do
!$omp enddo
!$omp end parallel

      ! The FFT routine uses internal parallelism (fftw), so it should be
      ! called outside of a parallel region

      ! Perform the real->spectral transform; the Sol_dwfft slice refers to the
      ! meaningful portion of this array.
      call execute_r2r_dft_plan(forward_plan, &
              Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
              Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)))


      ! Normalize, again with the help of OpenMP.
!$omp parallel private(i,j,k,jr,p0,pn,piece) &
!$omp          shared(plon,ptotal)
!$omp do
      do i = 0+Lam_pil_w, F_gni-1-Lam_pil_e
         do k = 1, F_nk
            do j = 1+pil_s, (F_t0njs-1+1)-pil_n
               Sol_dwfft(j,k,i+1) = Sol_pri * Sol_dwfft(j,k,i+1)
            end do
         end do
      end do
!$omp enddo

      ! Transpose in preparation for the tridiagonal solve.  Since the
      ! communication here involves MPI, it should only be called from
      ! one processor.
!$omp single
      call time_trace_barr(gem_time_trace, 11200+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)
      call rpn_comm_transpose( Sol_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                               1, F_t2nis, F_gni, Sol_dg2, 2, 2 )
      call time_trace_barr(gem_time_trace, 11300+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)
!$omp end single
      ptotal = F_t2ni-Sol_pil_e-Sol_pil_w-1
      plon   = (ptotal+Ptopo_npeOpenMP)/ Ptopo_npeOpenMP

      ! Perform the tridiagonal solve
!$omp do
      do piece=1,Ptopo_npeOpenMP
         p0 = 1+Sol_pil_w + plon*(piece-1)
         pn = min(F_t2ni-Sol_pil_e,plon*piece+Sol_pil_w)
         j =1+Lam_pil_s
         do i=p0,pn
            do k=1, F_nk
               Sol_dg2(k,i,j) = F_bi(k,i,j)*Sol_dg2(k,i,j)
            end do
         end do
         do j =2+Lam_pil_s, F_gnj-Lam_pil_n
            jr =  j - 1
            do i=p0,pn
               do k=1, F_nk
                  Sol_dg2(k,i,j) = F_bi(k,i,j)* Sol_dg2(k,i,j) - F_ai(k,i,j) &
                                            * Sol_dg2(k,i,jr)
               end do
            end do
         end do
         do j = F_gnj-1-Lam_pil_n, 1+Lam_pil_s, -1
            jr =  j + 1
            do i=p0,pn
               do k=1, F_nk
                  Sol_dg2(k,i,j) = Sol_dg2(k,i,j) - F_ci(k,i,j)*Sol_dg2(k,i,jr)
               end do
            end do
         end do
      end do
!$omp enddo
!$omp end parallel
      call time_trace_barr(gem_time_trace, 11400+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)
      call time_trace_barr(gem_time_trace, 11900+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

      ! Again transpose the data, for inversion of the FFT
      call rpn_comm_transpose( Sol_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                               1, F_t2nis, F_gni, Sol_dg2,- 2, 2 )
      call time_trace_barr(gem_time_trace, 11500+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

!     inverse projection ( r = x * w )
      call execute_r2r_dft_plan(reverse_plan, &
             Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
             Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)))
      call time_trace_barr(gem_time_trace, 11600+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

      ! And finally transpose the data into the block structure used
      ! by the rest of the model
      call rpn_comm_transpose ( F_Sol, 1, F_t0nis, F_gni, (F_t0njs-1+1), &
                                     1, F_t1nks, F_gnk,  Sol_dwfft, -1, 2)
      call time_trace_barr(gem_time_trace, 11700+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

!     __________________________________________________________________
!
      return
      end
