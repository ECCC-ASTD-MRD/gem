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

module sol
   use, intrinsic :: iso_fortran_env
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
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Sol_ai_8           |  sub-   diagonal of LU factorization            |
! Sol_bi_8           |         diagonal of LU factorization            |
! Sol_ci_8           |  super- diagonal of LU factorization            |
! Sol_stencil2_8,3,4,5 | stencils for Yin-Yang  (Qaddouri)             |
!----------------------------------------------------------------------

   integer :: sol_pil_w,sol_pil_e,sol_pil_n,sol_pil_s
   integer :: sol_niloc,sol_njloc,sol_nloc,sol_nk
   integer :: sol_i0,sol_in,sol_j0,sol_jn

   real, dimension(:,:,:), allocatable :: Sol_rhs

   real(kind=REAL64), dimension(:), allocatable :: Sol_ai_8,Sol_bi_8,Sol_ci_8
   real(kind=REAL64), dimension(:), allocatable :: Sol_stencil2_8,Sol_stencil3_8
   real(kind=REAL64), dimension(:), allocatable :: Sol_stencil4_8,Sol_stencil5_8

end module sol
