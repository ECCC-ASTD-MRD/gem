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

!**s/r  pre_redblack -- Red-black Gauss-Seidel preconditioning

module redblack_2d

contains

      subroutine pre_redblack2D( Lhs, Rhs, nil,njl, kidx )
      use cstv
      use geomh
      use glb_ld
      use opr
      use sol
      implicit none
#include <arch_specific.hf>
!
! Input parameters:
      integer nil, njl, &   ! Extent of the grid as seen by GMRES, containing no pilot and no halo region
              kidx          ! Index of the vertical level, used for defining the operator
      integer, parameter :: halo=1 ! Extent of the halo region
      !real*8, Dimension(nil,njl) :: Rhs(nil,njl), & ! Residual of the current problem
      !                              Lhs(nil,njl)    ! Estimated solution (preconditioner output)

      ! Use assumed-shape arrays, so no copy of the array slices is required
      real*8, Dimension(:,:) :: Rhs(:,:), & ! Residual of the current problem
                                Lhs(:,:)    ! Estimated solution (preconditioner output)

      real*8, Dimension((1-halo):(nil+halo),(1-halo):(njl+halo)) :: work ! Working array, including a halo region
!dir$ attributes align:64 :: work

!
! author   Christopher Subich - March 2017, adapted from pre_diago and mat_vecs2D by Abdessamad Qaddouri
!
!revision
! v4_5x - Subich C. - Revision for single-plane
!
!
      integer j,i,ii,jj ! Iteration and indexing variables
      integer ioffset     ! Row offset for red-black iteration
      real*8  cst,di_8    ! Parameters used for building the stencil
      real*8  stencil_pt, stencil_w, & ! Stencil elements: here, west
              stencil_e, stencil_s, & ! east, north
              stencil_n              ! south

!
!     ---------------------------------------------------------------
!

      ! Calculate the diagonal term of the matrix operator.

      ! The term corresponding to the current vertical mode, the sum of the
      ! helmholtz constant and the scaled eigenvalue of the vertical mode
      cst= (Cstv_hco1_8+Cstv_hco0_8*Opr_zeval_8(kidx))

      ! Initialize the working array
      work = 0.0

      ! Red-black Gauss-Seidel proceeds by performing diagonal preconditioning
      ! on one colour (call it red here), then solving for the other colour (black)
      ! based on the already-computed red values.  The colours are distributed in
      ! a checkerboard pattern, such that on the global array points where
      ! ii+jj = 1 mod 2 (mod(ii+jj,2)==1) are 'red' and the other points are 'black'
      do j=1, njl
         ! Convert the local j index to the global jj index
         jj=j+l_j0-1 + sol_pil_s
         di_8= Opr_opsyp0_8(G_nj+jj) / cos( G_yg_8 (jj) )**2
         ! Decide if we're iterating on i = 1,3,5 or i = 2,4,6
         ! Red points are ones where mod(ii+jj,2)==1, so to find the offset take
         ioffset = 1-mod(jj+l_i0+sol_pil_w,2)
         ! ... which corresponds to jj calculated above plus ii at the first index.

         do i=1+ioffset, nil, 2 ! NOTE! This loop steps by 2, as red-black is checkerboard
            ! Convert the local i index to the global ii index
            ii=i+l_i0-1 + sol_pil_w
            stencil_pt=cst*Opr_opsxp0_8(G_ni+ii)* &
              Opr_opsyp0_8(G_nj+jj) +Opr_opsxp2_8(G_ni+ii)*di_8+ &
              Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(G_nj+jj)
            work(i,j) = Rhs(i,j)/stencil_pt
            Lhs(i,j) = work(i,j)
         end do
      end do

      ! Now, perform a halo exchange.  Boundary black points will depend upon their out-of-domain
      ! red neighbours for computation.

      call rpn_comm_xch_halon( work, & ! Array to exchange
                               1-halo, nil + halo, 1-halo, njl + halo, & ! Lower and upper bounds of this array
                               nil, njl, 1, & ! Inner (non-halo) dimensions of the array and # of vertical levels
                               halo, halo, & ! Halo size
                               G_periodx, G_periody, & ! Grid periodicity
                               nil, 0, & ! Unused parameters -- necessary for global row accesses only
                               2) ! Data size, in units of integers.  8-byte doubles -> 2

      ! Now solve on the 'black' points.
      do j=1, njl
         ! Convert the local j index to the global jj index
         jj=j+l_j0-1 + sol_pil_s
         di_8= Opr_opsyp0_8(G_nj+jj) / cos( G_yg_8 (jj) )**2
         ! Black points are those where mod(ii+jj,2) == 0, so the computation is opposite that above
         ioffset = mod(jj+l_i0+sol_pil_w,2)

         do i=1+ioffset, nil, 2 ! NOTE! This loop steps by 2, as red-black is checkerboard
            ! Convert the local i index to the global ii index
            ii=i+l_i0-1 + sol_pil_w

            ! Compute now all 5 elements of the stencil.  Stencil_pt is the same as above
            stencil_pt=cst*Opr_opsxp0_8(G_ni+ii)* &
              Opr_opsyp0_8(G_nj+jj) +Opr_opsxp2_8(G_ni+ii)*di_8+ &
              Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(G_nj+jj)
            ! And the remainder are copied from sol_matvecs2D
            stencil_w= Opr_opsxp2_8(ii)*di_8
            stencil_e= Opr_opsxp2_8(2*G_ni+ii)*di_8
            stencil_s= Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(jj)
            stencil_n= Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(2*G_nj+jj)

            ! Solve for lhs(ii,jj), assuming the neighbour values in the work array are
            ! themselves correct lhs values and that work(i,j) is rhs(i,j)

            Lhs(i,j) = (Rhs(i,j) - stencil_w*work(i-1,j) - &
                         stencil_e*work(i+1,j) - stencil_s*work(i,j-1) - &
                         stencil_n*work(i,j+1)) / stencil_pt
         end do
      end do


!
!     ---------------------------------------------------------------
!
      return
      end subroutine pre_redblack2d

end module redblack_2d

