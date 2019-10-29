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

      subroutine sol_fft_lam ( sol, Rhs                        , &
                               F_t0nis, F_t0njs, F_t0nj        , &
                               F_t1nks, F_t1nk, F_t2nis, F_t2ni, &
                               F_gni, F_gnj, F_gnk, F_nk       , &
                               F_npex1, F_npey1                , &
                               F_ai, F_bi, F_ci  , F_dg2)
      use iso_c_binding
      use gem_fft
      use HORgrid_options
      use gem_options
      use gem_timing
      use glb_ld
      use glb_pil
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
      include 'mpif.h'
!
      integer F_t0nis, F_t0njs, F_t0nj, F_t2nis, F_t2ni
      integer F_t1nks, F_gnk,   F_nk  , F_t1nk , F_gni, F_gnj
      integer F_npex1, F_npey1

      real(kind=REAL64)  Sol (1:F_t0nis, 1:F_t0njs, F_gnk ), &
              Rhs (1:F_t0nis, 1:F_t0njs, F_gnk ), &
              F_ai(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_bi(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_ci(1:F_t1nks, 1:F_t2nis, F_gnj)
      real(kind=REAL64), allocatable, save ::  F_dwfft(:,:,:)
      real(kind=REAL64)  F_dg2  (1:F_t1nks, 1:F_t2nis, F_gnj  +F_npey1)

      ! FFTW plans
      type(dft_descriptor), save :: forward_plan, reverse_plan
      real(kind=REAL64), save :: pri ! Normalization constant
      external :: MPI_barrier

!author    Abdessamad Qaddouri- JULY 1999
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! Sol          O    - result of solver
! Rhs          I    - r.h.s. of elliptic equation
! Pri          I    - inverse projector in Fourier space
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
! F_dg2        I    - work field
! F_dwfft      I    - work field


      character(len=4), save :: type_fft
      integer i, j, k, jr
      integer, save :: l_pil_w, l_pil_e
      integer piece, p0, pn, plon, ptotal
      real(kind=REAL64), parameter :: zero = 0.d0
!     __________________________________________________________________
!
      if (.not. allocated(F_dwfft)) then
                            type_fft = 'QCOS'
         if (Grd_yinyang_L) type_fft = 'SIN'
         ! Perform initial setup for the Fourier transforms:
         ! Allocate the static array used for the transforms + transposes
         allocate(F_dwfft(1:F_t0njs, 1:F_t1nks, F_gni+2+F_npex1))

         ! Build the transform plans
         call make_r2r_dft_plan(forward_plan, & ! Plan variable
                                F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                                F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                                3, type_fft, DFT_FORWARD)
         call make_r2r_dft_plan(reverse_plan, & ! Plan variable
                                F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                                F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                                3, type_fft, DFT_BACKWARD)

      ! Get the appropriate scaling factor.  Begin with the normalization
      ! constant from the FFT, needed to make a round-trip of transforms
      ! act as the identity function:
         pri =  get_dft_norm_factor(G_ni-Lam_pil_w-Lam_pil_e,type_fft)

      ! And then modify the constant based on the discretization, essentially
      ! dividing by dx.  This calculation is modified from itf_fft_set, which
      ! is no longer used in the fftw-based interface.
         pri = pri/(G_xg_8(G_ni-Lam_pil_e+1)-G_xg_8(G_ni-Lam_pil_e))
!  The I vector lies on the Y processor so, l_pil_w and l_pil_e will
!  represent the pilot region along I

         l_pil_w=0
         l_pil_e=0
         if (l_south) l_pil_w= Lam_pil_w
         if (l_north) l_pil_e= Lam_pil_e
      end if

      call rpn_comm_transpose ( Rhs, 1, F_t0nis, F_gni, (F_t0njs-1+1), &
                                1, F_t1nks, F_gnk, F_dwfft, 1, 2 )

      ! Use OpenMP, if configured, to zero out unused portions of the
      ! transpose array
!$omp parallel private(i,j,k,jr,p0,pn,piece) &
!$omp          shared(plon,ptotal,l_pil_w,l_pil_e,pri)
!$omp do
      do i= 1,F_gni
         F_dwfft(F_t0nj+1-pil_n:F_t0njs,        1:F_t1nk ,i)= zero
         F_dwfft(             1:pil_s  ,        1:F_t1nk ,i)= zero
         F_dwfft(             1:F_t0njs, F_t1nk+1:F_t1nks,i)= zero
      end do
!$omp enddo
!$omp end parallel

      ! The FFT routine uses internal parallelism (fftw), so it should be
      ! called outside of a parallel region

      ! Perform the real->spectral transform; the F_dwfft slice refers to the
      ! meaningful portion of this array.
      call execute_r2r_dft_plan(forward_plan, &
                             F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                             F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)))

      ! Normalize, again with the help of OpenMP.
!$omp parallel private(i,j,k,jr,p0,pn,piece) &
!$omp          shared(plon,ptotal,l_pil_w,l_pil_e,pri)
!$omp do
      do i = 0+Lam_pil_w, F_gni-1-Lam_pil_e
         do k = 1, F_nk
            do j = 1+pil_s, (F_t0njs-1+1)-pil_n
               F_dwfft(j,k,i+1) = pri * F_dwfft(j,k,i+1)
            end do
         end do
      end do
!$omp enddo

      ! Transpose in preparation for the tridiagonal solve.  Since the
      ! communication here involves MPI, it should only be called from
      ! one processor.
!$omp single
      call rpn_comm_transpose( F_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                               1, F_t2nis, F_gni, F_dg2, 2, 2 )
!$omp end single

      ptotal = F_t2ni-l_pil_e-l_pil_w-1
      plon   = (ptotal+Ptopo_npeOpenMP)/ Ptopo_npeOpenMP

      ! Perform the tridiagonal solve
!$omp do
      do piece=1,Ptopo_npeOpenMP
         p0 = 1+l_pil_w + plon*(piece-1)
         pn = min(F_t2ni-l_pil_e,plon*piece+l_pil_w)
         j =1+Lam_pil_s
         do i=p0,pn
            do k=1, F_nk
               F_dg2(k,i,j) = F_bi(k,i,j)*F_dg2(k,i,j)
            end do
         end do
         do j =2+Lam_pil_s, F_gnj-Lam_pil_n
            jr =  j - 1
            do i=p0,pn
               do k=1, F_nk
                  F_dg2(k,i,j) = F_bi(k,i,j)* F_dg2(k,i,j) - F_ai(k,i,j) &
                                            * F_dg2(k,i,jr)
               end do
            end do
         end do
         do j = F_gnj-1-Lam_pil_n, 1+Lam_pil_s, -1
            jr =  j + 1
            do i=p0,pn
               do k=1, F_nk
                  F_dg2(k,i,j) = F_dg2(k,i,j) - F_ci(k,i,j)*F_dg2(k,i,jr)
               end do
            end do
         end do
      end do
!$omp enddo
!$omp end parallel

      ! Again transpose the data, for inversion of the FFT
      call rpn_comm_transpose( F_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                               1, F_t2nis, F_gni, F_dg2,- 2, 2 )

!     inverse projection ( r = x * w )
      call execute_r2r_dft_plan(reverse_plan, &
                             F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                             F_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)))

      ! And finally transpose the data into the block structure used
      ! by the rest of the model
      call rpn_comm_transpose ( Sol, 1, F_t0nis, F_gni, (F_t0njs-1+1), &
                                     1, F_t1nks, F_gnk,  F_dwfft, -1, 2)

!     __________________________________________________________________
!
      return
      end
