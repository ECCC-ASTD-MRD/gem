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
!**s/r sol_fft_numa - parallel direct solution of an elliptic problem
!                     for LAM grids using FFT in shared memory segment
!                     of sockets on a machine node

      subroutine sol_fft_numa ( F_sol, F_Rhs, F_fft, F_a, F_b, F_c )
      use iso_c_binding
      use glb_ld
      use glb_pil
      use ptopo
      use ldnh
      use sol
      use numa
      use trp
      use adz_mem
      use dyn_fisl_options
      use gem_timing
      use, intrinsic :: iso_fortran_env
      implicit none
#include <RPN_MPI.hf>
#include <arch_specific.hf>

      real(kind=REAL64), dimension(ldnh_maxx,ldnh_maxy,G_nk), intent(IN ) :: F_Rhs
      real(kind=REAL64), dimension(ldnh_maxx,ldnh_maxy,G_nk), intent(OUT) :: F_Sol
      real(kind=REAL64), dimension(Sol_mink:Sol_maxk,G_ni,Sol_miny  :Sol_maxy  ), intent(IN ) :: F_a, F_b
      real(kind=REAL64), dimension(Sol_mink:Sol_maxk,G_ni,Sol_miny-1:Sol_maxy  ), intent(IN ) :: F_c
      real(kind=REAL64), dimension(Sol_mink:Sol_maxk,G_ni,Sol_miny-1:Sol_maxy+1), intent(out) :: F_fft

      integer i,j,k, i0,in, j0,jn,jr, n, err,dim,tag1,tag2, proc
      integer status(MPI_STATUS_SIZE)

      logical, save :: transpose_setup = .false.
      integer, save :: row_communicator, column_communicator
!     __________________________________________________________________
!
      n= (Adz_icn-1)*Schm_itnlh + Adz_itnl
      
! Forward transpose (resuls= G_ni in memory on all Pes)
      call gemtime_start ( 53, 'TRP1', 24 )
      ! Transpose to the 4D representation.  This representation takes the local array
      ! of shape (lni, lnj, gnk) and returns an array of shape (lni, lnj, lnk, npex),
      ! which turns the MPI portion of the transpose into the sending/receiving of single,
      ! large chunks of contiguous memory. 

      if (.not. transpose_setup) then
         ! The stacked-transpose function requires an intiailization step, so complete this
         ! once at the first call.
         row_communicator = RPN_COMM_comm('EW')
         column_communicator = RPN_COMM_comm('NS')
         call RPN_MPI_transpose_setup(G_nk, sol_nk, row_communicator, column_communicator,err)
         if (err /= 0) then
            call gem_error(-1,'SOL_FFT_NUMA','Invalid MPI transpose setup')
         end if
         transpose_setup = .true.
      end if

      ! Perform the stacked transpose

      call RPN_MPI_ez_transpose_xz(LoC(F_rhs), LoC(Sol_xpose), .true., &  ! Arrays, forward transform
               ldnh_maxx*2, ldnh_maxy, sol_nk, & ! Bounds on the transpose, x multiplied by 2 for double precision
               err)

      ! Afterwards, we need to reshape this array into a contiguous internal representation.
      ! FFTW can deal well with arrays that have non-unit strides, but this 4D form doesn't
      ! have a single, consistent stride along the first dimension.  To maximize locality
      ! during this in-core transpose, we'll go with a local representation of dimensions
      ! (gni, lnj, lnk) -- the array can be written to in a single pass.

      do k=1,sol_nk
         ! Looping from 1:ldnh_maxy covers the full in-memory extent of the array, inclusive
         ! of piloting regions
         do j=1,ldnh_maxy
            ! Since the source data is broken up by processor, we want to first loop over
            ! the processor, then over their "owned" i-regions.
            do proc=1,Ptopo_npex
               ! Grab each processor's i-range and copy it to our integrated array

               ! Using the Ptopo_gindx_alongX arary preserves generality against future changes to
               ! the domain decomposition, by not building in assumptions about how the array
               ! is intially distributed.

               do i=Ptopo_gindx_alongX(1,proc),Ptopo_gindx_alongX(2,proc)
                  Sol_dwfft(i,j,k) = Sol_xpose(i-Ptopo_gindx_alongX(1,proc)+1,j,k,proc)
               end do
            end do
         end do
      end do

      call gemtime_stop (53)
      call gemtime_start ( 54, 'FFTW', 24 )
      call time_trace_barr(gem_time_trace, 11100+n, Gem_trace_barr, &
                                      Ptopo_intracomm, MPI_BARRIER)

      ! The FFT routine uses internal parallelism (fftw), so it should be
      ! called outside of a parallel region

      ! This executes an out-of-place transform.  The source array is Sol_dwfft, which
      ! contains our MPI-transposed field, and the destination is F_fft, the shared-
      ! memory region used in the tridiagonal solver.  Sol_dwfft has memory order (j,k,i),
      ! and F_fft has memory order (k,i,j); FFTW handles the permutation internally.

      call execute_r2r_dft_plan(forward_plan, &
            Sol_dwfft((1+Lam_pil_w):(G_ni-Lam_pil_e),        & ! Source array, 1st dim = i
                      (1+pil_s):(ldnh_nj-pil_n),             & ! 2nd dim = j
                      1:sol_nk),                             & ! 3rd dim = k
            F_fft(trp_12sn0:(trp_12sn0 + sol_nk - 1),        & ! Dest array, 1st dim = k
                  (1+Lam_pil_w):(G_ni-Lam_pil_e),            & ! 2nd dim = i
                  (ldnh_j0+pil_s):(ldnh_nj+ldnh_j0-pil_n-1)))  ! 3rd dim = j

      call gemtime_stop (54)

      ! The MPI barrier acts to ensure that each process in this NUMA region has completed
      ! its transform and assigned the output to F_fft.  If this is not present, then a
      ! process could run ahead of that shared-memory assignment and operate on invalid
      ! data.
      call MPI_Barrier ( Numa_sockcomm, err )
      call gemtime_start ( 56, 'TRID', 24 )

      ! Tri-diagonal solver
      ! We only need to execute the tridiagonal solve if there are points to solve for
      ! on this socket
      if (Sol_sock_nk > 0) then

         ! Forward pass

         ! Bookkeeping to set up appropriate bounds
         dim= Sol_sock_nk*G_ni ; tag1=95 ; tag2=96
         i0= max(1+Lam_pil_w,Sol_istart)
         in= min(G_ni-Lam_pil_e,Sol_iend)
         
         ! This socket isn't the start of the process -- wait to receive partial work from
         ! the miny part
         if (Sol_miny /= 1) then
            if (Numa_sockrank == 0) &
            call MPI_RECV(F_fft(Sol_mink,1,Sol_miny-1), dim    ,&
                          MPI_DOUBLE_PRECISION, Numa_peerrank-1,&
               tag1+Numa_peerrank-1 , Numa_peercomm, status, err)
            call MPI_Barrier ( Numa_sockcomm, err )
         endif
         
         j0= max(Sol_miny,2   +Lam_pil_s)
         jn= min(Sol_maxy,G_nj-Lam_pil_n)

         ! The FFT normalization factor has not yet been applied, but we can do this alongside the
         ! application of F_b, when F_fft is read for the first time inside the tridiagonal solve.
         if (Sol_miny <= 1 + Lam_pil_s) then 
            ! The first j-column does not involve a lower-index value of j, so it must be
            ! normalized separately.
            j = 1+Lam_pil_s
            do i = i0, in
               do k=Sol_mink, Sol_maxk
                  F_fft(k,i,j) = Sol_pri*F_fft(k,i,j)*F_b(k,i,j)
               end do
            end do
         end if
         do j= j0, jn
            jr =  j - 1
            do i= i0, in
               do k= Sol_mink, Sol_maxk
                  ! Complete the first (forward substitution) step of the Thomas / ripple algorithm.
                  ! Additionally, apply normalization to the RHS F_fft(k,i,j).
                  F_fft(k,i,j) = Sol_pri*F_fft(k,i,j)*F_b(k,i,j) - F_a(k,i,j)*F_fft(k,i,jr)
               end do
            end do
         end do

         ! This socket isn't the end of the process -- send the partial work to the next socket
         if (Sol_maxy /= G_nj) then
            call MPI_Barrier ( Numa_sockcomm, err )
            if (Numa_sockrank == 0) &
            call MPI_SEND(F_fft(Sol_mink,1,jn), dim            ,&
                          MPI_DOUBLE_PRECISION, Numa_peerrank+1,&
                          tag1+Numa_peerrank, Numa_peercomm, err)
         endif

         !BACKWARD
         j0= max(Sol_miny,1   +Lam_pil_s)
         if (Sol_maxy == G_nj) then ! End of the line, so begin the backward step
            do j = G_nj-1-Lam_pil_n, j0, -1
               jr =  j + 1
               do i= i0, in
                  do k= Sol_mink, Sol_maxk
                     F_fft(k,i,j)= F_fft(k,i,j)-F_c(k,i,j)*F_fft(k,i,jr)
                  end do
               end do
            end do
         else ! Not the end of the line, so receive partial work from the next socket
            if (Numa_sockrank == 0) &
            call MPI_RECV( F_fft(Sol_mink,1,Sol_maxy+1), dim, &
                        MPI_DOUBLE_PRECISION, Numa_peerrank+1,&
              tag2+Numa_peerrank+1, Numa_peercomm, status, err)
            call MPI_Barrier ( Numa_sockcomm, err )
            
            do j = Sol_maxy, j0, -1
               jr =  j + 1
               do i= i0, in
                  do k= Sol_mink, Sol_maxk
                     F_fft(k,i,j)= F_fft(k,i,j)-F_c(k,i,j)*F_fft(k,i,jr)
                  end do
               end do
            end do
         endif
         if (Sol_miny /= 1) then ! Not the start of the line, so send partial work to the prior socket
            call MPI_Barrier ( Numa_sockcomm, err )
            if (Numa_sockrank == 0) &
            call MPI_SEND( F_fft(Sol_mink,1,j0), dim, &
                MPI_DOUBLE_PRECISION, Numa_peerrank-1,&
               tag2+Numa_peerrank, Numa_peercomm, err )
         endif

         ! End backward pass
         call MPI_Barrier ( Numa_sockcomm, err )
      endif ! End if sol_sock_nk > 0
      
      call gemtime_stop (56)


! Inverse projection ( r = x * w )
      call gemtime_start ( 54, 'FFTW', 24 )

      ! Execute the inverse FFT.  The source is now the shared-memory F_fft array, whereas
      ! the destination is Sol_dwfft.

      call execute_r2r_dft_plan(reverse_plan, &
            F_fft(trp_12sn0:(trp_12sn0 + sol_nk - 1),          & ! Source array, 1st dim = k
                    (1+Lam_pil_w):(G_ni-Lam_pil_e),            & ! 2nd dim = i
                    (ldnh_j0+pil_s):(ldnh_nj+ldnh_j0-pil_n-1)),& ! 3rd dim = j
            Sol_dwfft((1+Lam_pil_w):(G_ni-Lam_pil_e),        & ! Source array, 1st dim = i
                      (1+pil_s):(ldnh_nj-pil_n),             & ! 2nd dim = j
                      1:sol_nk))                               ! 3rd dim = k

      call gemtime_stop (54)
      call time_trace_barr(gem_time_trace, 11600+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

! Backward transpose to bring back the solution in model space
      call gemtime_start ( 53, 'TRP1', 24 )

      ! The solution is now in Sol_dwfft, but we need to create the stacked representation for the
      ! new inverse transpose.  This time, to assign memory in contiguous order we must loop over
      ! processors first:
      do proc=1,Ptopo_npex
         do k=1,sol_nk
            do j=1,ldnh_maxy
               do i=Ptopo_gindx_alongX(1,proc),Ptopo_gindx_alongX(2,proc)
                 Sol_xpose(i-Ptopo_gindx_alongX(1,proc)+1,j,k,proc) = Sol_dwfft(i,j,k)
               end do
            end do
         end do
      end do

      ! Perform the stacked inverse transpose
      call RPN_MPI_ez_transpose_xz(LoC(F_Sol), LoC(Sol_xpose), .false., &  ! Arrays, reverse transpose
               ldnh_maxx*2, ldnh_maxy, sol_nk, & ! Bounds on the transpose, x multiplied by 2 for double precision
               err)

      call gemtime_stop (53)
      call time_trace_barr(gem_time_trace, 11700+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)
!     __________________________________________________________________
!
      return
      end subroutine sol_fft_numa
