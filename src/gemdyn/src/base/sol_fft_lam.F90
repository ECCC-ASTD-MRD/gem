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
      use lun
      use ldnh
      use, intrinsic :: iso_fortran_env
      implicit none
#include <RPN_MPI.hf>
#include <arch_specific.hf>

      integer F_t0nis, F_t0njs, F_t0nj, F_t2nis, F_t2ni
      integer F_t1nks, F_gnk,   F_nk  , F_t1nk , F_gni, F_gnj

      real(kind=REAL64)  F_Sol (1:F_t0nis, 1:F_t0njs, F_gnk ), &
              Rhs (1:F_t0nis, 1:F_t0njs, F_gnk ), &
              F_ai(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_bi(1:F_t1nks, 1:F_t2nis, F_gnj), &
              F_ci(1:F_t1nks, 1:F_t2nis, F_gnj)

      ! Variables to manage the "stacked" x/z transpose
      logical, save :: xtranspose_setup = .false.
      integer, save :: row_communicator, column_communicator
      integer       :: proc, err


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
      call gemtime_start ( 53, 'TRP1', 24 )

      ! Update: 2020 Dec, csubich
      ! RPN_Comm_transpose inside the RPNComm library assumes that the x-transpose reduces
      ! to an even transpose of a potentially larget set of vertical levels.  E.g, 5 levels
      ! distributed over 4 processors would split as [2, 2, 1, 0] levels after the transpose,
      ! which is a truncated form of splitting 8 levels over 4 processors.  This has obvious
      ! implications for load balancing, which ultimately affect solver performance.

      ! The newer ``stacked'' transpose formulation is more flexible, but it requires some
      ! postprocessing.  This invocation is very similar to that in sol_fft_numa.

      if (.not. xtranspose_setup) then
         ! The stacked transpose function requires an intiailzation step, which we perform before
         ! the first transpose call.
         row_communicator = RPN_COMM_comm('EW')
         column_communicator = RPN_COMM_comm('NS')
         call RPN_MPI_transpose_setup(G_nk, sol_nk, row_communicator, column_communicator, err)
         if (err /= 0) then
            call gem_error(-1,'SOL_FFT_LAM','Invalid MPI transpose setup')
         end if
         xtranspose_setup = .true.
      end if

      ! Perform the stacked x/z transpose
      call RPN_MPI_ez_transpose_xz(LoC(Rhs), LoC(Sol_xpose), .true., & ! Arrays, marking the forward transpose
               ldnh_maxx*2, ldnh_maxy, sol_nk, & ! Logical bounds on the transpose, multiplying x by 2 to account for double precision
               err)

      ! Now, the MPI-tranposed array Sol_xpose must be reshaped into a contiguous internal representation.
      ! The destination of the transpose is in [j, k, i] order, and we wish to write this linearly.
      ! The i-indices of the MPI-transposed array are split amongst the (:,:,:,npex) Sol_xpose in processor order,
      ! and that forms an implied outer loop layer.
      do proc=1,Ptopo_npex
         ! Read Ptopo_gindx_alongX to determine the logical bounds given to us by processor (proc)
         do i=Ptopo_gindx_alongX(1,proc),Ptopo_gindx_alongX(2,proc)
            do k=1,sol_nk
               do j=1,ldnh_maxy
                  Sol_dwfft(j,k,i) = Sol_xpose(i-Ptopo_gindx_alongX(1,proc)+1,j,k,proc)
               end do
            end do
         end do
      end do

      call gemtime_stop (53)
      call gemtime_start ( 54, 'FFTW', 24 )

      do i= 1,F_gni
         Sol_dwfft(F_t0nj+1-pil_n:F_t0njs,        1:F_t1nk ,i)= zero
         Sol_dwfft(             1:pil_s  ,        1:F_t1nk ,i)= zero
         Sol_dwfft(             1:F_t0njs, F_t1nk+1:F_t1nks,i)= zero
      end do
      call execute_r2r_dft_plan(forward_plan, &
              Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
              Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)))

      do i = 0+Lam_pil_w, F_gni-1-Lam_pil_e
         do k = 1, F_nk
            do j = 1+pil_s, (F_t0njs-1+1)-pil_n
               Sol_dwfft(j,k,i+1) = Sol_pri * Sol_dwfft(j,k,i+1)
            end do
         end do
      end do
      call gemtime_stop (54)
      call gemtime_start ( 55, 'TRP2', 24 )
      call rpn_comm_transpose( Sol_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                               1, F_t2nis, F_gni, Sol_dg2, 2, 2 )
      call gemtime_stop (55)
      call gemtime_start ( 56, 'TRID', 24 )

      ptotal = F_t2ni-Sol_pil_e-Sol_pil_w-1
      plon   = (ptotal+Ptopo_npeOpenMP)/ Ptopo_npeOpenMP
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
      call gemtime_stop (56)

      call gemtime_start ( 55, 'TRP2', 24 )
      ! Invert the x/y transpose
      call rpn_comm_transpose( Sol_dwfft, 1, F_t0njs, F_gnj, (F_t1nks-1+1),&
                               1, F_t2nis, F_gni, Sol_dg2,- 2, 2 )
      call gemtime_stop (55)

      call gemtime_start ( 54, 'FFTW', 24 )
      call execute_r2r_dft_plan(reverse_plan, &
             Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
             Sol_dwfft((1+pil_s):(F_t0njs-pil_n),1:F_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)))
      call gemtime_stop (54)

      call gemtime_start ( 53, 'TRP1', 24 )
      ! Invert the x/z transpose, again using the "stacked" transpose routine.  First, repackage
      ! Sol_dwfft into the stacked representation:

      ! We wish to write Sol_xpose in a contiguous order, so the loops should follow
      ! its natural [i,j,k,proc] ordering.
      do proc=1,Ptopo_npex
         do k=1,sol_nk
            do j=1,ldnh_maxy
               do i=Ptopo_gindx_alongX(1,proc),Ptopo_gindx_alongX(2,proc)
                  Sol_xpose(i-Ptopo_gindx_alongX(1,proc)+1,j,k,proc) = Sol_dwfft(j,k,i)
               end do
            end do
         end do
      end do

      ! Next, call the stacked transposer with the inverse flag:
      call RPN_MPI_ez_transpose_xz(LoC(F_Sol), LoC(Sol_xpose), .false., & ! Arrays, marking the inverse transpose
               ldnh_maxx*2, ldnh_maxy, sol_nk, & ! Logical bounds on the transpose, multiplying x by 2 to account for double precision
               err) ! Error flag

      call gemtime_stop (53)

!     __________________________________________________________________
!
      return
      end
