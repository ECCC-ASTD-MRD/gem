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

module sol_mem
   use, intrinsic :: iso_fortran_env
   use gem_fft
   implicit none
   public
   save

!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH THE SOLVER                                |
!                            -1                                        |
!   we solve:   (I + Ai) * Bi   * (I + Ci) * X = RHS                   |
!       with:   Ai = Sol_ai_8, Bi = Sol_bi_8 and Ci = Sol_ci_8         |
!______________________________________________________________________|
!                      |                                                 |
! NAME                 | DESCRIPTION                                     |
!----------------------|-------------------------------------------------|
! Sol_ai_8             |  sub-   diagonal of LU factorization            |
! Sol_bi_8             |         diagonal of LU factorization            |
! Sol_ci_8             |  super- diagonal of LU factorization            |
! Sol_stencil2_8,3,4,5 | stencils for Yin-Yang  (Qaddouri)               |
! Sol_stencilh         | stencils for H coordinates MATVEC               |
! sol_stencilp         | stencils for P coordinates MATVEC               |
!-------------------------------------------------------------------------

   character(len=4) :: Sol_type_fft
   integer :: Sol_pil_w,Sol_pil_e,Sol_pil_n,Sol_pil_s
   integer :: Sol_niloc,Sol_njloc,Sol_nloc,Sol_nk,Sol_sock_nk
   integer :: Sol_i0,Sol_in,Sol_j0,Sol_jn,Sol_numank
   integer :: Sol_miny,Sol_maxy,Sol_mink,Sol_maxk,Sol_k0
   integer :: Sol_ldni, Sol_istart, Sol_iend
   integer :: Sol_ii0,Sol_iin,Sol_jj0,Sol_jjn,Sol_imin,Sol_imax,Sol_jmin,Sol_jmax

   real, dimension(:,:,:), allocatable :: Sol_rhs

   real(kind=REAL64), dimension(:), allocatable :: Sol_ai_8,Sol_bi_8,Sol_ci_8
   real(kind=REAL64), dimension(:), allocatable :: Sol_stencil2_8,Sol_stencil3_8
   real(kind=REAL64), dimension(:), allocatable :: Sol_stencil4_8,Sol_stencil5_8
   real(kind=REAL64), dimension(:,:,:), allocatable :: Sol_dwfft, Sol_dg2, &
                                                       Sol_rhs_8, Sol_sol_8

   real(kind=REAL64), dimension(:,:,:,:), allocatable :: Sol_stencilp_8, Sol_stencilh_8

   !! For the one-transpose solve, the following are pointers into the allocated,
   !  shared-memory space that eliminates the y-transpose.  Defining these pointers
   !  as contiguous (and ensuring it on assignment) enables compile-time optimizations.
   !  - Sol_fft contains the transposed output of the FFT
   !  - Sol_a, Sol_b, and Sol_c contain the precomputed tridiagonal solve entries
   !  Each array is in (k,i,j) order.
   !  - Sol_xpose contains the 4D working array (lni, lnj, lnk, npex) for the improved
   !              MPI transposer
   real(kind=REAL64), dimension(:,:,:), contiguous, pointer :: Sol_a,Sol_b,Sol_c,Sol_fft
   real(kind=REAL64), dimension(:,:,:,:), contiguous, pointer :: Sol_xpose
   real(kind=REAL64), dimension(:,:,:  ), pointer :: Sol_saved=>null()

   real(kind=REAL64) :: isol_i, isol_d
   real(kind=REAL64) :: norm_residual, relative_tolerance, nu , r0, rr2, ro2, lcl_sum(2)
   real(kind=REAL64),dimension(:      ), allocatable :: gg,rot_cos, rot_sin
   real(kind=REAL64),dimension(:,:    ), allocatable :: v_lcl_sum,rr,tt, hessenberg,thread_s
   real(kind=REAL64),dimension(:,:,:  ), allocatable :: work_space,thread_s2,fdg,w2_8,w3_8
   real(kind=REAL64),dimension(:,:,:,:), allocatable :: vv, wint_8
   real(kind=REAL64),dimension(:,:,:,:), allocatable :: A1,A2,B1,B2,C1,C2
   real(kind=REAL128),dimension(:,:    ), allocatable :: thread_s128
   real             ,dimension (:,:,: ), allocatable :: fdg2

   ! Normalization factor for the FFT
   real(kind=REAL64) :: Sol_pri

   ! Opaque pointers to hold the forward and reverse FFT plans
   type(dft_descriptor) :: forward_plan, reverse_plan

end module sol_mem
