!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1930-2010 - Division de Recherche en Prevision Numerique
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

! pre_multigrid_jac3D.F90 -- Preconditioner for the 3D elliptic problem via
!                            multigrid coarsening along the (i,j) dimensions
!                            and line relaxation along the (k) dimension.  This
!                            approach uses a (1,1) V-cycle with an underweighted-
!                            (line) Jacobi iteration as the smoother.
!
!                            As with pre_redblack3D, this module derives the
!                            operator by using the matrix defined in the matvec
!                            module.  Additionally, this preconditioner primarily
!                            operates in single precision.
!
! Author:
!     Christopher Subich - July 2018
! Revision
!     v5.1 - Subich C. - Initial version

module multigrid_3d_jac
   implicit none
   ! With a sufficiently large halo region, the relaxation can proceed with just
   ! one boundary exchange per relaxation pass.  This comes at the expense of
   ! using more memory and redundant computations between nodes.
   integer, parameter :: HALO = 2 ! 1 is an absolute minimum

   ! Some procedures, such as the tridiagonal solve, can be explicitly vectorized
   ! even if the natural loop ordering would seem to prevent it.  The following
   ! parameter sets the target (maximal) vector length.  For single-precision
   ! floating point, reasonable options are
   !integer, parameter :: VECLEN = 4 ! Width of an SSE register
   !integer, parameter :: VECLEN = 8 ! Width of an AVX register
   integer, parameter :: VECLEN = 16 ! Width of a cache line, also allows loop unrolling

   ! Rather than go all the way to a coarsest grid, the multigrid iteration can be
   ! terminated at an early level.  There, the coarse problem is solved not by
   ! a direct or exact solver, but by a few relaxation passes.  This parameter sets
   ! the number of levels, the last of which is treated as "the" coarse level.
   integer, parameter :: MAX_LEVEL = 5
   integer :: LEVELS

   ! The multigrid V-cycle uses an underweighted (k-block) Jacobi iteration for the
   ! relaxation step.  Matlab testing suggested 0.87 was a good value for the weight
   ! parameter, but it remains to see if this can be forever fixed or needs to vary.
   real*4, parameter :: JACWEIGHT = 0.95

   ! The multigrid iteration needs the operator not just at the default, finest
   ! level, but also at a number of coarse levels along the way.  It is convenient
   ! to store these operators in advance, here using a specific derived datatype
   ! to keep both the operator and grid information together.  Each instance of
   ! this type corresponds to a multigrid level
   type matrix_info
      ! The matrix on this level, with indexing (i,j,k,STENCIL).
      real*4, dimension(:,:,:,:), allocatable :: mat
!dir$ attributes align: 64 :: mat

      ! Cached factors for the tridiagonal solve, derived via the Thomas algorithm:
      ! In this algorithm, the LHS is solved via a forward and backward substitution pass:
      ! d(k) = rhs(k)*scal(k) - m(k)*d(k-1)
      ! lhs(k) = d(k) - cp(k)*lhs(k+1)
      ! where b, m, and cp are derived from the matrix coefficients.  These can be kept
      ! between passes, saving on some repeated floating-point work
      real*4, dimension(:,:,:), allocatable :: scal, cp, m
!dir$ attributes align: 64 :: scal, cp, m

      ! The number of global points in i, j, k at this level, using the prefix ge
      ! for "global extent."  For the purposes of this module, the global extent
      ! at the finest level is not necessarily the same grid seen in GEM, but
      ! instead it is the portion of the grid touched by the linear solver;
      ! the point i=1 corresponds to (1+glb_pil_w) for example.
      integer :: ge_i, ge_j, ge_k

      ! The *local* lower and upper bounds on i and j at this level.  In contrast
      ! to the indexing convention in the rest of GEM, this module uses a global
      ! array ordering: process N will see array bounds that aren't 1->l_ni or some
      ! variation thereof.  This makes index math much easier, since indices are
      ! multiplied or divided by 2 with each level change.
      !
      ! These variables denote the area of primary responsibility for this proess,
      ! and it is exclusive  of the halo region.  The nomenclature uses the prefix
      ! l{lb,ub} for "local, lower/upper bound"
      integer :: llb_i, lub_i, llb_j, lub_j

      ! The allocated lower bounds on i and j at this level are *inclusive* of
      ! the halo region, and they may also be expanded for array alignment
      ! purposes pending vectorization tests.  These variables use the prefix
      ! "a{lb,ub}" for "allocated".  Note also that all internal arrays will use
      ! this structure, even if a halo isn't locally necessary -- this preserves
      ! alignment for vectorization
      integer :: alb_i, aub_i, alb_j, aub_j
   end type matrix_info

   ! Pre-define an array of the above type to hold the information as-calculated
   type (matrix_info) grid_info(MAX_LEVEL)

   ! Flag for whether the grid_info has been initialized
   logical :: grids_initialized = .false.

contains

   ! Public interface function, called from the solver
   subroutine pre_multigrid_jac3d(Lhs, Rhs)
      ! pre_multigrid_jac3d: estimates the left-hand-side Lhs corresponding
      ! to the given right-hand-side Rhs, by applying a single multigrid V-cycle

      ! Since this procedure acts as the interface to the linear solver, it needs
      ! to reference its grid dimensions.  The 3D Krylov solvers are called from
      ! sol_iterative3D, which defines the linear solver grid based on
      ! ldnh_{min,max}[xy], so those are used here as well.

      use ldnh, only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy ! Allocated limits for Lhs/Rhs
      use glb_ld, only: l_ni, l_nj, l_nk, l_i0, l_j0 ! Local grid offsets and extents
      use glb_pil, only: glb_pil_s, glb_pil_w ! Global pilot region sizes
      use sol, only: sol_pil_w, sol_pil_e, sol_pil_s, sol_pil_n ! Local pilot region sizes

      real*8, intent(in), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk) :: Rhs
      real*8, intent(out), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk) :: Lhs

      ! Define internal rhs/sol arrays for exchange with the v-cycle.  These arrays
      ! use the global grid numbering, but that numbering may not yet be defined
      ! if grid_info is uninitialized.  So, declare these arrays allocatable:

      real*4, allocatable, dimension(:,:,:) :: rhs_v, sol_v, sol_inc
!dir$ attributes align: 64 :: rhs_v, sol_v, sol_inc

      ! This interface routine copies between these array, so define loop variables:
      integer :: i, j, k

      ! Convenient shortcut variables:
      integer :: alb_i, aub_i, alb_j, aub_j, ge_i, ge_j, ge_k, llb_i, lub_i, llb_j, lub_j

      ! Alignment note: we do _not_ assume that Lhs and Rhs have any particular alignment;
      ! they're allocated and managed by krylov_mod

      ! Initialize the grid operators if not already done
      if (.not. grids_initialized) then
         call init_levels()
      end if

      ! Copy the top-level array bounds to local variables, for shorter naming below
      alb_i = grid_info(1)%alb_i
      aub_i = grid_info(1)%aub_i
      alb_j = grid_info(1)%alb_j
      aub_j = grid_info(1)%aub_j
      ge_i  = grid_info(1)%ge_i
      ge_j  = grid_info(1)%ge_j
      ge_k  = grid_info(1)%ge_k
      llb_i = grid_info(1)%llb_i
      lub_i = grid_info(1)%lub_i
      llb_j = grid_info(1)%llb_j
      lub_j = grid_info(1)%lub_j

      ! Now, we know the bounds are there to allocate {rhs,sol}_v


      allocate (rhs_v(alb_i:aub_i, alb_j:aub_j, ge_k))
      allocate (sol_v(alb_i:aub_i, alb_j:aub_j, ge_k))
      allocate (sol_inc(alb_i:aub_i, alb_j:aub_j, ge_k))

      ! Initialize sol_v to zero.  This should not numerically affect results,
      ! but the M*x calculations do execute multiplications that can reference
      ! outside the valid region.  0*NaN causes a floating point exception.
      sol_v = 0
      sol_inc = 0

      ! Copy the incoming RHS to the internal array; this also converts
      ! it from REAL*8 to REAL*4

      do k=1,ge_k
         rhs_v(:,alb_j:(llb_j-1),k) = 0 ! Explicitly zero unused parts of RHS
         do j=llb_j,lub_j
            rhs_v(alb_i:(llb_i-1),j,k) = 0
            do i=llb_i,lub_i
               rhs_v(i,j,k) = Rhs(i + glb_pil_w - l_i0 + 1, j + glb_pil_s - l_j0 + 1, k)
            end do
            rhs_v((lub_i+1):aub_i,j,k) = 0
         end do
         rhs_v(:,(lub_j+1):aub_j,k) = 0
      end do

      !sumsq2 = sum(real(rhs_v,8)**2)


      ! Conduct two V-cycles
      call vcycle(rhs_v, sol_v, 1, .true.) ! The first also updates rhs_v with the new residual
      call vcycle(rhs_v, sol_inc, 1, .false.) ! while the second leaves rhs_v in an undefined state

      ! Write the output LHS back to the solution

      do k=1,ge_k
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               Lhs(i,j,k) = sol_v(i - glb_pil_w + l_i0 - 1, j - glb_pil_s + l_j0 - 1, k) + &
                          sol_inc(i - glb_pil_w + l_i0 - 1, j - glb_pil_s + l_j0 - 1, k)
            end do
         end do
      end do



   end subroutine


   ! Initialize the grid_info structures, to be called once
   subroutine init_levels()
      ! Use the matrix and associated index enum from the matvec module
      use matvec, only : stencil, IDX_WEST, IDX_EAST, IDX_NORTH, IDX_SOUTH, IDX_POINT, IDX_TOP, IDX_BOTTOM

      ! Indices required to make sense of the matvec array
      !use ldnh, only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy
      use glb_ld, only: l_nk, l_i0, l_j0, & ! Vertical level count, local->global grid offsets
                        G_ni, G_nj, l_ni, l_nj ! Number of points in i/j, global/local
      use sol, only: sol_pil_n, sol_pil_e ! Local solver-specific pilot regions
      use glb_pil, only: glb_pil_s, glb_pil_w, glb_pil_e, glb_pil_n ! Global pilot region sizes

      ! Loop variables
      integer :: i, j, k, ci, cj, ilevel

      ! Variables for bounds-management at each level
      integer :: alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated bounds
                 llb_i, lub_i, llb_j, lub_j,       & ! Local bounds
                 ge_i, ge_j                          ! Global extents
      integer :: extent_i ! Working variable for i-extent
      integer :: loc_max_level, ierr ! Local variables for local level computation and MPI error

      ! Begin with the finest level, which is set based on global grid/matrix parameters
      ilevel = 1

      ge_i = G_ni - glb_pil_w - glb_pil_e ! Number of points in the global grid, i direction
      ge_j = G_nj - glb_pil_n - glb_pil_s ! Idem, j direction
      ge_k = l_nk

      llb_i = max(1,l_i0 - glb_pil_w) ! Lower bound of this process along i, in global-solver coordinates
      lub_i = l_i0 + l_ni - sol_pil_e - glb_pil_w - 1! Upper bound of this process along i
      llb_j = max(1,l_j0 - glb_pil_s) ! Lower bound along j
      lub_j = l_j0 + l_nj - sol_pil_n - glb_pil_s - 1! Upper bound along j

      ! Allocation bounds are given by the local bounds plus the halo region

      ! A vector-optimization implemented here is also to ensure that (aub_i-alb_i)
      ! is a multiple of 16 (VECLEN). This ensures that arr(i,j,k) has the same alignment
      ! characteristic (modulo 64 bytes) as arr(i,j+1,k) and arr(i,j,k+1).

      ! ifort 16.0 (the deployed version as of this writing) doesn't make the most-possible
      ! use of this information, but this adjustment does seem to improve performance in
      ! the line Jacobi relaxation by about 15%.
      alb_i = llb_i - HALO
      aub_i = lub_i + HALO
      extent_i = (aub_i-alb_i+1)
      aub_i = alb_i + VECLEN*((extent_i+VECLEN-1)/VECLEN) - 1
      alb_j = llb_j - HALO
      aub_j = lub_j + HALO

      ! Write these bounds to the grid_info structure
      grid_info(ilevel)%llb_i = llb_i
      grid_info(ilevel)%lub_i = lub_i
      grid_info(ilevel)%llb_j = llb_j
      grid_info(ilevel)%lub_j = lub_j
      grid_info(ilevel)%alb_i = alb_i
      grid_info(ilevel)%aub_i = aub_i
      grid_info(ilevel)%alb_j = alb_j
      grid_info(ilevel)%aub_j = aub_j
      grid_info(ilevel)%ge_i  = ge_i
      grid_info(ilevel)%ge_j  = ge_j
      grid_info(ilevel)%ge_k  = ge_k

      ! We can now allocate the finest-level matrix
      allocate (grid_info(ilevel)%mat(alb_i:aub_i, alb_j:aub_j, ge_k, 7))

      ! And the tridiagonal terms
      allocate (grid_info(ilevel)%scal(alb_i:aub_i, alb_j:aub_j, ge_k))
      allocate (grid_info(ilevel)%cp(alb_i:aub_i, alb_j:aub_j, ge_k))
      allocate (grid_info(ilevel)%m(alb_i:aub_i, alb_j:aub_j, ge_k))

      ! Copy the matrix from 'matrix' in the matvec module, transposing to an
      ! (i,j,k,IDX) ordering, aligning with the grid conventions in this module,
      ! and downconverting to REAL*4
      do k=1,ge_k
         do j=alb_j,(llb_j-1)
            grid_info(ilevel)%mat(:,j,k,:) = 0
         end do
         do j=llb_j,lub_j
            grid_info(ilevel)%mat(alb_i:(llb_i-1),j,k,:) = 0
            do i=llb_i,lub_i
               grid_info(ilevel)%mat(i,j,k,:) = stencil(i+glb_pil_w-l_i0+1, & ! Convert i to the matrix-local numbering
                                                        j+glb_pil_s-l_j0+1, & ! idem for j
                                                        :,k)
            end do
            grid_info(ilevel)%mat((lub_i+1):aub_i,j,k,:) = 0
         end do
         do j=(lub_j+1),aub_j
            grid_info(ilevel)%mat(:,j,k,:) = 0
         end do
      end do

      ! Perform a halo exchange of the matrix, for calculation of the coarse-grid operator
      call rpn_comm_xch_halo(grid_info(ilevel)%mat, & ! Exchange the operator
                             1-HALO, aub_i - (alb_i+HALO-1), & ! Adjusted allocated bounds along i
                             1-HALO, aub_j - (alb_j+HALO-1), & ! ... and j
                             lub_i-llb_i+1, lub_j-llb_j+1, & ! Extent of the interior (valid) region
                             ge_k*7, & ! "Flatten" the array, to treat the stencil alongside vertical levels
                             HALO, HALO, 0, 0, & ! Halo size, nonperiodic
                             0, 0) ! Unused -- related to grids that cover the poles

      ! Build the tridiagonal factors
      call build_trid_factors(grid_info(ilevel)%mat, & ! The operator
                              grid_info(ilevel)%scal, grid_info(ilevel)%cp, grid_info(ilevel)%m, & ! The L/U factors
                              alb_i, aub_i, alb_j, aub_j, & ! Allocated array bounds
                              llb_i, lub_i, llb_j, lub_j, & ! Local array bounds
                              ge_i, ge_j, ge_k) ! Global array extents

      ! Debugging printouts -- uncomment to receive a per-processor report on grid allocation
      !write(0,'("Multigrid Level ",I1," initialized, bounds ",2("[",I4,"-",I4,"-",I4,"-",I4,"] "),I3)'), ilevel, alb_i, llb_i, lub_i, aub_i, alb_j, llb_j, lub_j, aub_j, ge_k
      !write(0,'("              global extents ",4I4.3)'), ge_i, ge_j, ge_k

      ! Define the local grid bounds at each level.  Along the way, determine the maximum level we
      ! can safely support on this processor -- interior cells need at least HALO points in both
      ! i and j in order for subsequent exchanges to work properly.

      loc_max_level = MAX_LEVEL ! Start with the optimistic assumption, that this processors supports all levels
      do ilevel=2,MAX_LEVEL
         ! For coarsening, the index math is straightforward but nitpicky: every second
         ! point, beginning with (i,j)=1, is kept for the coarse grid.  Visually, this gives:

         ! Fine:    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
         ! Coarse:  1  -  2  -  3  -  4  -  5  -  6  -  7  -  8

         ! We want a processor to own a coarse-grid point if and only if
         ! it owns the fine-grid point that exactly corresponds in location.

         ! Assuming the fine-grid extent was at least 2 (it should be at least 2*HALO),
         ! we get:
         llb_i = llb_i/2+1   ! Round the lower bound up, s.t. 1->1, 2->2, 3->2, 4->3, etc
         lub_i = (lub_i+1)/2 ! Round the upper bound down, s.t. 1->1, 2->1, 3->2, 4->2, etc
         llb_j = llb_j/2+1
         lub_j = (lub_j+1)/2
         ge_i  = (ge_i+1)/2
         ge_j  = (ge_j+1)/2

         ! Allocated bounds are not coarsened from the previous level; instead they are based
         ! on local bounds plus and minus the halo region
         alb_i = llb_i-HALO
         aub_i = lub_i+HALO
         extent_i = (aub_i-alb_i+1)
         aub_i = alb_i + VECLEN*((extent_i+VECLEN-1)/VECLEN) - 1
         alb_j = llb_j-HALO
         aub_j = lub_j+HALO

         ! Write these bounds to the grid_info structure
         grid_info(ilevel)%llb_i = llb_i
         grid_info(ilevel)%lub_i = lub_i
         grid_info(ilevel)%llb_j = llb_j
         grid_info(ilevel)%lub_j = lub_j
         grid_info(ilevel)%alb_i = alb_i
         grid_info(ilevel)%aub_i = aub_i
         grid_info(ilevel)%alb_j = alb_j
         grid_info(ilevel)%aub_j = aub_j
         grid_info(ilevel)%ge_i  = ge_i
         grid_info(ilevel)%ge_j  = ge_j
         grid_info(ilevel)%ge_k  = ge_k

         ! Now, check for invalidity: an interior range of insufficient size
         if ( (llb_i > 1 .and. lub_i < ge_i .and. (lub_i-llb_i+1) < HALO) .or. &
              (llb_j > 1 .and. lub_j < ge_j .and. (lub_j-llb_j+1) < HALO) .or. &
              (lub_i < llb_i) .or. (lub_j < llb_j)) then
            loc_max_level = ilevel-1
            exit
         end if
      end do

      ! Now, use MPI_REDUCE to find the lowest supported level amongst all processors
      call RPN_COMM_allreduce(loc_max_level,LEVELS,1,"MPI_INTEGER","MPI_MIN","grid",ierr)

      ! With that set, initialize the matrix operators
      do ilevel=2,LEVELS
         ! Retrieve array bounds from the structure
         llb_i = grid_info(ilevel)%llb_i
         lub_i = grid_info(ilevel)%lub_i
         llb_j = grid_info(ilevel)%llb_j
         lub_j = grid_info(ilevel)%lub_j
         alb_i = grid_info(ilevel)%alb_i
         aub_i = grid_info(ilevel)%aub_i
         alb_j = grid_info(ilevel)%alb_j
         aub_j = grid_info(ilevel)%aub_j
         ge_i  = grid_info(ilevel)%ge_i
         ge_j  = grid_info(ilevel)%ge_j
         ge_k  = grid_info(ilevel)%ge_k

         ! Allocate the matrix operator
         allocate (grid_info(ilevel)%mat(alb_i:aub_i, alb_j:aub_j, ge_k, 7))
         ! And the tridiagonal terms
         allocate (grid_info(ilevel)%scal(alb_i:aub_i, alb_j:aub_j, ge_k))
         allocate (grid_info(ilevel)%cp(alb_i:aub_i, alb_j:aub_j, ge_k))
         allocate (grid_info(ilevel)%m(alb_i:aub_i, alb_j:aub_j, ge_k))

         ! The coarse-grid operator is defined by Mc*xc = R*M*(P*xc): the effect of
         ! applying the coarse-grid operator is the same as interpolating to the fine
         ! grid, applying the fine-grid operator, and coarsening.

         ! The V-cycle below uses full-weighting for the coarsening operator alongside
         ! linear interpolation.  Doing so for the operator, however, broadens the
         ! Laplacian from a 7-point stencil to a 27-point one -- the "point" and up/down
         ! values diffuse to the neighbouring cells.

         ! Using _subsampling_ (injection) for the coarsening operator, however, maintains
         ! a 7-point stencil.  Testing in MATLAB shows that this has a slight effect on
         ! the predicted convergence rate, but the smaller stencil should offer much
         ! greater efficiency.

         ! Iterate over the elements of the coarse matrix.  Point (ci,cj) on the coarse level
         ! corresponds to the point (2*ci-1, 2*cj-1) on the fine level,
         do k=1,ge_k
            do cj=alb_j,(llb_j-1)
               grid_info(ilevel)%mat(:,cj,k,:) = 0
            end do
            do cj=llb_j,lub_j

               if (cj < ge_j) then
                  ! General case
!dir$ ivdep
                  do ci=llb_i,min(ge_i-1,lub_i)
                     ! The point (ci,cj), when interpolated to the fine grid, has nonzero support from
                     ! i=(2*ci-2 .. 2*ci) and j=(2*cj-2 .. 2*cj), giving:

                     ! 1/4 1/2 1/4
                     ! 1/2  1  1/2
                     ! 1/4 1/2 1/4

                     ! Therefore, the coarse-level "point" term will involve contributions from all
                     ! of its horizontal components:

                     grid_info(ilevel)%mat(ci,cj,k,IDX_POINT) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_POINT)

                     ! While the other horizontal neighbours are scaled -- a coarse-grid delta function to the
                     ! west, for example, is seen at this point with a 1/2-weighting.

                     grid_info(ilevel)%mat(ci,cj,k,IDX_WEST)  = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_EAST)  = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_SOUTH) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_NORTH) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH)

                     ! And the vertical operators remain unchanged, since the interpolation acts only in the horizontal
                     ! and subsampling-coarsening does not average terms.

                     grid_info(ilevel)%mat(ci,cj,k,IDX_TOP) = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_TOP)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_BOTTOM) = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_BOTTOM)
                  end do

                  if (lub_i == ge_i) then
                     ci = lub_i
                     ! Case for max-i.  Here, the interpolation stencil is modified because the underlying operator
                     ! has a Neumann boundary condition; a coarse-grid delta expands to

                     ! 1/4 1/2 1/2
                     ! 1/2  1   1
                     ! 1/4 1/2 1/2

                     ! So the coarse-grid generation is also modified for the EAST-related terms:

                     grid_info(ilevel)%mat(ci,cj,k,IDX_POINT) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_POINT)

                     grid_info(ilevel)%mat(ci,cj,k,IDX_WEST)  = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_EAST)  = 0!  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_SOUTH) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_NORTH) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH)

                     ! And the vertical operators remain unchanged, since the interpolation acts only in the horizontal
                     ! and subsampling-coarsening does not average terms.

                     grid_info(ilevel)%mat(ci,cj,k,IDX_TOP)    = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_TOP)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_BOTTOM) = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_BOTTOM)
                  end if
               else
                  ! Here, cj == ge_j, and we cannot guarantee that no points between this location and the max-j
                  ! boundary have been dropped.  The interpolation consequently needs adjusted

                  ! General-i case:
!dir$ ivdep
                  do ci=llb_i,min(ge_i-1,lub_i)
                     ! The coarse-grid delta interpolates to:

                     ! 1/4 1/2 1/4  (ge_j-2)
                     ! 1/2  1  1/2  (ge_j-1)
                     ! 1/2  1  1/2  (ge_j)

                     ! Therefore, the coarse-level "point" term will involve contributions from all
                     ! of its horizontal components:

                     grid_info(ilevel)%mat(ci,cj,k,IDX_POINT) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_POINT)

                     ! While the other horizontal neighbours are scaled -- a coarse-grid delta function to the
                     ! west, for example, is seen at this point with a 1/2-weighting.

                     grid_info(ilevel)%mat(ci,cj,k,IDX_WEST)  = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_EAST)  = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_SOUTH) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_NORTH) = 0!  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH)

                     ! And the vertical operators remain unchanged, since the interpolation acts only in the horizontal
                     ! and subsampling-coarsening does not average terms.

                     grid_info(ilevel)%mat(ci,cj,k,IDX_TOP) = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_TOP)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_BOTTOM) = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_BOTTOM)
                  end do

                  if (lub_i == ge_i) then
                     ci = lub_i
                     ! Case for max-i, modifying things in both i and j.  The coarse-grid delta interpolates to:

                     ! 1/4 1/2 1/2
                     ! 1/2  1   1
                     ! 1/2  1   1

                     ! So the coarse-grid generation is also modified for the EAST-related terms:

                     grid_info(ilevel)%mat(ci,cj,k,IDX_POINT) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST) + &
                                                                0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH) + &
                                                                  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_POINT)

                     grid_info(ilevel)%mat(ci,cj,k,IDX_WEST)  = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_WEST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_EAST)  = 0!  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_EAST)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_SOUTH) = 0.5*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_SOUTH)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_NORTH) = 0!  1*grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_NORTH)

                     ! And the vertical operators remain unchanged, since the interpolation acts only in the horizontal
                     ! and subsampling-coarsening does not average terms.

                     grid_info(ilevel)%mat(ci,cj,k,IDX_TOP)    = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_TOP)
                     grid_info(ilevel)%mat(ci,cj,k,IDX_BOTTOM) = grid_info(ilevel-1)%mat(2*ci-1,2*cj-1,k,IDX_BOTTOM)
                  end if ! Endif lub_i == ge_i
               end if ! Endif j == ge_j
            end do ! End do j
         end do ! End do k

         ! Perform a halo exchange of the matrix
         call rpn_comm_xch_halo(grid_info(ilevel)%mat, & ! Exchange the operator
                                1-HALO, aub_i - (alb_i+HALO-1), & ! Adjusted allocated bounds along i
                                1-HALO, aub_j - (alb_j+HALO-1), & ! ... and j
                                lub_i-llb_i+1, lub_j-llb_j+1, & ! Extent of the interior (valid) region
                                ge_k*7, & ! "Flatten" the array, to treat the stencil alongside vertical levels
                                HALO, HALO, 0, 0, & ! Halo size, nonperiodic
                                0, 0) ! Unused -- related to grids that cover the poles

         ! Allocation debug printouts
         !write(0,'("Multigrid Level ",I1," initialized, bounds ",2("[",I4,"-",I4,"-",I4,"-",I4,"] "),I3)'), ilevel, alb_i, llb_i, lub_i, aub_i, alb_j, llb_j, lub_j, aub_j, ge_k
         !write(0,'("              global extents ",4I4.3)'), ge_i, ge_j, ge_k

         ! Build the tridiagonal factors
         call build_trid_factors(grid_info(ilevel)%mat, & ! The operator
                                 grid_info(ilevel)%scal, grid_info(ilevel)%cp, grid_info(ilevel)%m, & ! The L/U factors
                                 alb_i, aub_i, alb_j, aub_j, & ! Allocated array bounds
                                 llb_i, lub_i, llb_j, lub_j, & ! Local array bounds
                                 ge_i, ge_j, ge_k) ! Global array extents


      end do ! Do over levels

      grids_initialized = .true.

   end subroutine

   ! Build LU factors for tridiagonal solve
   subroutine build_trid_factors(mat,scal,cp,m,alb_i,aub_i,alb_j,aub_j,llb_i,lub_i,llb_j,lub_j,ge_i,ge_j,ge_k)
      use matvec, only : IDX_POINT,IDX_BOTTOM,IDX_TOP

      ! Array bounds, both allocation and local:
      integer, intent(in) :: alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated bounds, inclusive of halo
                             llb_i, lub_i, llb_j, lub_j,       & ! Local bounds, exclusive of halo
                             ge_i, ge_j                          ! Global extents over all processors
      ! Matrix
      real*4, intent(in) :: mat(alb_i:aub_i,alb_j:aub_j,ge_k,7)

      ! Storage arrays for the tridiagonal factors
      real*4, intent(out), dimension(alb_i:aub_i,alb_j:aub_j,ge_k) :: scal, cp, m

      ! Loop variables
      integer :: i, j, k

      ! The vertical line-solve in this multigrid method is based on the Thomas algorithm
      ! (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm), which implements
      ! Gaussian substitution without pivoting.  This algorithm is reasonably fast, and
      ! pivoting is unnecessary because the vertical-only operator is diagonally dominant.

      ! Here, we precompute the factors derived by this algorithm, so that the tridiagonal
      ! solve reduces to a forward and backward substitution pass requiring only multiplication
      ! and addition.  When executing the Thomas algorithm in full at each V-cycle, the
      ! line solve takes the lion's share of CPU time.

      ! Zero the factor arrays outside of the written range
      scal(:,alb_j:max(0,llb_j-HALO-1),:) = 0
      cp(:,alb_j:max(0,llb_j-HALO-1),:) = 0
      m(:,alb_j:max(0,llb_j-HALO-1),:) = 0
      do j=max(1,llb_j-HALO),min(ge_j,lub_j+HALO)
         scal(alb_i:max(0,llb_i-HALO-1),j,:) = 0
         cp(alb_i:max(0,llb_i-HALO-1),j,:) = 0
         m(alb_i:max(0,llb_i-HALO-1),j,:) = 0
         do i=max(1,llb_i-HALO),min(ge_i,lub_i+HALO)
            ! To translate from the linked reference:
            ! b(k) = mat(i,j,k,IDX_POINT) -- matrix diagonal
            ! c(k) = mat(i,j,k,IDX_BOTTOM) -- matrix superdiagonal
            ! a(k) = mat(i,j,k,IDX_TOP) -- matrix subdiagonal

            ! Calculations from the forward sweep:
            ! scal(k) = 1/(b(k) - a(k)*cp(k-1))
            ! cp(k) = c(k)*scal(k)
            ! m(k) =  a(k)*scal(k)

            ! k = 1
            scal(i,j,1) = 1.0/mat(i,j,1,IDX_POINT)
            cp(i,j,1)   = mat(i,j,1,IDX_BOTTOM)*scal(i,j,1)

            do k=2,ge_k
               scal(i,j,k) = 1.0/(mat(i,j,k,IDX_POINT)-mat(i,j,k,IDX_TOP)*cp(i,j,k-1))
               cp(i,j,k)   = mat(i,j,k,IDX_BOTTOM)*scal(i,j,k)
               m(i,j,k)    = mat(i,j,k,IDX_TOP)*scal(i,j,k)
            end do
         end do
         scal(min(ge_i+1,lub_i+HALO+1):aub_i,j,:) = 0
         cp(min(ge_i+1,lub_i+HALO+1):aub_i,j,:) = 0
         m(min(ge_i+1,lub_i+HALO+1):aub_i,j,:) = 0
      end do
      scal(:,min(ge_j+1,lub_j+HALO+1):aub_j,:) = 0
      cp(:,min(ge_j+1,lub_j+HALO+1):aub_j,:) = 0
      m(:,min(ge_j+1,lub_j+HALO+1):aub_j,:) = 0
   end subroutine ! build_trid_factors

   ! V-cycle management
   recursive subroutine vcycle(rhs_v, sol_v, mylevel,update_rhs)
      ! Perform a single multigrid V-cycle

      ! The RHS array is passed with intent(inout) to allow for in-place updates; we
      ! have already allocated it to contain the proper halo region

      ! On input: rhs_v contains valid right-hand-side information on the region
      ! (llb_i:lub_i,llb_j:lub_j,:) and zeroes outside this region
      ! On output: if update_rhs is .true., the remaining residual rhs_v - M*sol_v;
      !            if update_rhs is .false., undefined (an intermediate result)
      integer, intent(in) :: mylevel
      real*4, contiguous, intent(inout) :: rhs_v(grid_info(mylevel)%alb_i:,grid_info(mylevel)%alb_j:,:)

      logical, intent(in) :: update_rhs

      ! On output, sol_v contains a valid left-hand-side estimate on the region
      ! (llb_i:lub_i,llb_j:lub_j,:), and undefined outside this region
      real*4, contiguous, intent(out)   :: sol_v(grid_info(mylevel)%alb_i:,grid_info(mylevel)%alb_j:,:)


      ! Fine-grid working array, for storing the solution increment
      real*4 :: sol_inc(grid_info(mylevel)%alb_i : grid_info(mylevel)%aub_i, &
                        grid_info(mylevel)%alb_j : grid_info(mylevel)%aub_j, &
                        grid_info(mylevel)%ge_k)
!dir$ attributes align: 64 :: sol_inc

      ! Coarse-grid working arrays, defined as allocatable because they are unused
      ! at the coarsest level
      real*4, allocatable :: rhs_cg(:,:,:), sol_cg(:,:,:)
!dir$ attributes align: 64 :: rhs_cg, sol_cg

      ! Loop and bounds variables
      integer :: j, k, ci, cj                  ! Loop variables
      integer :: alb_i, aub_i, alb_j, aub_j, ge_k ! Allocation variables
      integer :: llb_i, lub_i, llb_j, lub_j       ! Local lower/upper bounds
      integer :: ge_i, ge_j                       ! Global array extents in i, j
      integer :: alb_ci, aub_ci, alb_cj, aub_cj   ! Allocation variables for the coarse grid
      integer :: llb_ci, lub_ci, llb_cj, lub_cj   ! Lower/upper bounds for the coarse grid
      integer :: ge_ci, ge_cj                     ! Global array extents in i, j for coarse grid

      ! Alignment directive for vectorization
!dir$ assume_aligned rhs_v:64, sol_v:64

      sol_inc = 0

      ! Mnemonics for frequently-used grid information
      alb_i = grid_info(mylevel)%alb_i
      aub_i = grid_info(mylevel)%aub_i
      alb_j = grid_info(mylevel)%alb_j
      aub_j = grid_info(mylevel)%aub_j
      ge_i  = grid_info(mylevel)%ge_i
      ge_j  = grid_info(mylevel)%ge_j
      ge_k  = grid_info(mylevel)%ge_k
      llb_i = grid_info(mylevel)%llb_i
      lub_i = grid_info(mylevel)%lub_i
      llb_j = grid_info(mylevel)%llb_j
      lub_j = grid_info(mylevel)%lub_j

      ! Step zero: exchange the halo region

      ! rpn_comm_xch_halo is... awkward.  It makes the implicit assumption that
      ! the input array has a lower bound of 1-HALO along i and j, so we need
      ! to adjust the input bounds to match this requirement.  That might not
      ! meet the allocated bounds of the array, so first test this:
      if (alb_i /= llb_i - HALO) then
         STOP 'alb_i != llb_i - HALO' !, alb_i, llb_i, HALO, mylevel
      end if
      if (alb_j /= llb_j - HALO) then
         STOP 'alb_j != llb_j - HALO' !, alb_j, llb_j, HALO, mylevel
      end if
      call rpn_comm_xch_halo(rhs_v, & ! Exchange the rhs_v array
                             1-HALO, aub_i - (alb_i+HALO-1), & ! Adjusted allocated bounds along i
                             1-HALO, aub_j - (alb_j+HALO-1), & ! ... and j
                             lub_i-llb_i+1, lub_j-llb_j+1, ge_k, & ! Extents of the interior region
                             HALO, HALO, 0, 0, & ! Halo size, nonperiodic
                             0, 0) ! Unused -- related to grids that cover the poles

      ! Step one: perform a Jacobi pass of the full region

      ! This call directly writes the solution to sol_v.  Other calls will need to use sol_inc,
      ! but since we know this is the first relaxation pass we can avoid the extra copy.

      call line_jacobi(rhs_v, sol_v, JACWEIGHT, & ! RHS, solution, and underrelaxation weight
                       grid_info(mylevel)%scal, & ! Matrix tridiagonal factors
                       grid_info(mylevel)%cp, &   !
                       grid_info(mylevel)%m, &    !
                       alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated array bounds
                       max(1,llb_i-HALO), min(ge_i,lub_i+HALO), & ! Array subset i, confined to the global domain
                       max(1,llb_j-HALO), min(ge_j,lub_j+HALO)) ! Array subset j

      ! With a halo size of 2, this means that we have a valid estimate for x from lb-2 to ub+2
      ! Step two: update the RHS, based on the horizontal part of the operator only

      call rhs_update_ij(sol_v, rhs_v, grid_info(mylevel)%mat, JACWEIGHT, & ! Perform rhs = rhs - mat*sol
                         alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated bounds
                         max(1,llb_i-(HALO-1)), min(ge_i,lub_i+(HALO-1)), & ! Valid subset i
                         max(1,llb_j-(HALO-1)), min(ge_j,lub_j+(HALO-1)))   ! Valid subset j

      ! With a halo size of 2, this means we have a valid calculation for f from lb-1 to ub+1

      ! Step three: the recursive call
      if (mylevel < LEVELS) then ! We're not on the finest level, so make the recursive call
         ! Define coarse-grid bounds
         alb_ci = grid_info(mylevel+1)%alb_i
         aub_ci = grid_info(mylevel+1)%aub_i
         alb_cj = grid_info(mylevel+1)%alb_j
         aub_cj = grid_info(mylevel+1)%aub_j
         llb_ci = grid_info(mylevel+1)%llb_i
         lub_ci = grid_info(mylevel+1)%lub_i
         llb_cj = grid_info(mylevel+1)%llb_j
         lub_cj = grid_info(mylevel+1)%lub_j
         ge_ci  = grid_info(mylevel+1)%ge_i
         ge_cj  = grid_info(mylevel+1)%ge_j

         ! Allocate the coarse-grid LHS and RHS arrays
         allocate(sol_cg(alb_ci:aub_ci,alb_cj:aub_cj,ge_k))
         allocate(rhs_cg(alb_ci:aub_ci,alb_cj:aub_cj,ge_k))

         ! Initialize the coarse-grid array; failing to do so can cause floating-
         ! point exceptions when M*x tries to multiply junk data by 0
         sol_cg = 0

         ! For coarsening, the index math is straightforward but nitpicky: every second
         ! point, beginning with (i,j)=1, is kept for the coarse grid.  Visually, this gives:

         ! Fine:    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
         ! Coarse:  1  -  2  -  3  -  4  -  5  -  6  -  7  -  8

         ! So ci = (i+1)/2. and i = 2*ci - 1


         ! Step 3.1: Coarsen rhs_cg = R*rhs
         do k = 1, ge_k
            ! Zero the unused portions of rhs_cg
            rhs_cg(:,alb_cj:(llb_cj-1),k) = 0
            do cj = llb_cj, lub_cj
               rhs_cg(alb_ci:(llb_ci-1),cj,k) = 0
               do ci = llb_ci, lub_ci
                  ! Apply full weighting, which is the stencil [0.25 0.5 0.25] * [0.25 0.5 0.25]^T
                  rhs_cg(ci,cj,k) = 0.0625*rhs_v(2*ci-2,2*cj-2,k) + & ! (-1,-1)
                                    0.125 *rhs_v(2*ci-1,2*cj-2,k) + & ! ( 0,-1)
                                    0.0625*rhs_v(2*ci  ,2*cj-2,k) + & ! (+1,-1)
                                    0.125 *rhs_v(2*ci-2,2*cj-1,k) + & ! (-1, 0)
                                    0.25  *rhs_v(2*ci-1,2*cj-1,k) + & ! ( 0, 0)
                                    0.125 *rhs_v(2*ci  ,2*cj-1,k) + & ! (+1, 0)
                                    0.0625*rhs_v(2*ci-2,2*cj  ,k) + & ! (-1,+1)
                                    0.125 *rhs_v(2*ci-1,2*cj  ,k) + & ! ( 0,+1)
                                    0.0625*rhs_v(2*ci  ,2*cj  ,k)     ! (+1,+1)
               end do
               rhs_cg((lub_ci+1):aub_ci,cj,k) = 0
            end do
            rhs_cg(:,(lub_cj+1):aub_cj,k) = 0
         end do

         ! Step 3.2: Invoke vcycle

         call vcycle(rhs_cg,sol_cg,mylevel+1,.false.)

         ! Step 3.3: Halo exchange of coarse-grid results

         ! The objective here is to fill the halo on coarse-grid-x such that
         ! we can interpolate a fine-grid-x correction valid with a halo size
         ! of 1.  (That then lets us update RHS on the interior and compute
         ! the final Jacobi correction).

         ! Structurally, this means we only need a coarse-halo of 1.  However,
         ! rpn_comm_xch_halo doesn't appear to give us that flexibility, so
         ! we exchange a full halo instead.

         call rpn_comm_xch_halo(sol_cg, & ! Exchange the coarse-grid solution
                              1-HALO, aub_ci - (alb_ci+HALO-1), & ! Adjusted allocated coarse bounds on i
                              1-HALO, aub_cj - (alb_cj+HALO-1), & ! Adjusted allocated coarse bounds on j
                              lub_ci - llb_ci + 1, lub_cj - llb_cj + 1, ge_k, & ! Extent of the interior
                              HALO, HALO, 0, 0, & ! Halo size, nonperiodic
                              0, 0) ! Unused

         ! Step 3.4: Interpolate coarse-grid results to fine grid

         ! The interpolation process is a straightforward linear interpolation, the
         ! transpose of the full weighting restriction operator above.  For unit-stride
         ! access, we'll use a modest degree of loop unrolling

         do k=1, ge_k
            do j=max(1,llb_j-1), min(ge_j,lub_j+1)
               if (mod(j,2) == 1 .or. j == ge_j) then
                  ! If j is odd, then the interpolation along the j-direction
                  ! is an identity -- the row corresponds to one kept on the
                  ! coarse grid, and the interpolation stencil uses either one
                  ! or two points in i.

                  ! This also applies to points along the max-j boundary, even
                  ! if it is not a row kept on the coarse grid.  This is a
                  ! consequence of the implied Neumann boundary condition
                  ! for the linear operator.

                  ! Note: these loops can write past [llb-1,lub+1] because of
                  ! the index math involved.  This should not be a problem provided
                  ! HALO is at least 2, but care should be taken with a smaller
                  ! halo region.

                  !do ci=max(1,llb_ci-1),min(ge_ci-1,lub_ci+1)

                  ! Loop over the coarse array, and write the points in the fine
                  ! array corresponding to the direct-injection point and the
                  ! point to the right.  This relies on (ci) and (ci+1).

                  do ci=max(1,(llb_i-1+1)/2), & ! Convert the fine point with halo into a coarse index to remain in-bounds
                        min(ge_ci-1,(lub_i+1+1)/2)
                     ! The left point on the fine grid has a coarse-grid equivalent
                     sol_inc(2*ci-1,j,k) = sol_cg(ci,(j+1)/2,k)

                     ! The right point is interpolated between adjacent CG points
                     sol_inc(2*ci  ,j,k) = 0.5*sol_cg(ci  ,(j+1)/2,k) + &
                                           0.5*sol_cg(ci+1,(j+1)/2,k)
                  end do

                  ! The rightmost point on the global grid is special, for the
                  ! same reasons as max-j
                  if (lub_ci == ge_ci) then
                     sol_inc(2*ge_ci-1,j,k) = sol_cg(ge_ci,(j+1)/2,k)
                     sol_inc(2*ge_ci  ,j,k) = sol_cg(ge_ci,(j+1)/2,k)
                  end if
               else
                  ! j is even and not at the very top, so we need to interpolate
                  ! between coarse-j rows.
                  do ci=max(1,(llb_i-1+1)/2), & ! Convert the fine point with halo into a coarse index to remain in-bounds
                        min(ge_ci-1,(lub_i+1+1)/2)
                     ! The left point on the fine grid interpolates between coarse-j
                     ! levels, but has a matching coarse-i point
                     sol_inc(2*ci-1,j,k) = 0.5 *sol_cg(ci  ,j/2  ,k) + &
                                           0.5 *sol_cg(ci  ,j/2+1,k)

                     ! The right point on the fine grid interpolates from four
                     ! coarse points
                     sol_inc(2*ci  ,j,k) = 0.25*sol_cg(ci  ,j/2  ,k) + &
                                           0.25*sol_cg(ci+1,j/2  ,k) + &
                                           0.25*sol_cg(ci  ,j/2+1,k) + &
                                           0.25*sol_cg(ci+1,j/2+1,k)
                  end do

                  ! The global-max-i test is a special case, as above
                  if (lub_ci == ge_ci) then
                     sol_inc(2*ge_ci-1,j,k) = 0.5*sol_cg(ge_ci,j/2  ,k) + &
                                              0.5*sol_cg(ge_ci,j/2+1,k)
                     sol_inc(2*ge_ci  ,j,k) = 0.5*sol_cg(ge_ci,j/2  ,k) + &
                                              0.5*sol_cg(ge_ci,j/2+1,k)
                  end if
               end if
            end do ! j
         end do ! k

         ! Step 3.5: Update RHS based on coarse-grid correction

         call rhs_update_ijk(sol_inc,rhs_v,grid_info(mylevel)%mat, & ! rhs = rhs - mat*sol_inc
                             alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated bounds
                             llb_i, lub_i, llb_j, lub_j) ! Valid subset, no halo region

         ! Add the solution increment to sol_v

         sol_v(llb_i:lub_i,llb_j:lub_j,:) =   sol_v(llb_i:lub_i,llb_j:lub_j,:) + &
                                            sol_inc(llb_i:lub_i,llb_j:lub_j,:)

      end if

      ! Step four: Perform a Jacobi pass of the restricted (no halo) region
      call line_jacobi(rhs_v, sol_inc, JACWEIGHT, & ! RHS, solution, and underrelaxation weight
                       grid_info(mylevel)%scal, & ! Matrix tridiagonal factors
                       grid_info(mylevel)%cp, &   !
                       grid_info(mylevel)%m, &    !
                       alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated array bounds
                       llb_i, lub_i, llb_j, lub_j) ! Valid subset, no halo region

      if (update_rhs) then
         call rpn_comm_xch_halo(sol_inc, & ! Exchange the solution increment
                                1-HALO, aub_i - (alb_i+HALO-1), & ! Adjusted allocated bounds along i
                                1-HALO, aub_j - (alb_j+HALO-1), & ! ... and j
                                lub_i-llb_i+1, lub_j-llb_j+1, ge_k, & ! Extents of the interior region
                                HALO, HALO, 0, 0, & ! Halo size, nonperiodic
                                0, 0) ! Unused -- related to grids that cover the poles
         call rhs_update_ij(sol_inc, rhs_v, grid_info(mylevel)%mat, JACWEIGHT, & ! Perform rhs = rhs - mat*sol_inc
                            alb_i, aub_i, alb_j, aub_j, ge_k, & ! Allocated bounds
                            llb_i, lub_i, llb_j, lub_j) ! Interior region
      end if

      ! Add the solution increment to sol_v
      sol_v(llb_i:lub_i,llb_j:lub_j,:) =   sol_v(llb_i:lub_i,llb_j:lub_j,:) + &
                                         sol_inc(llb_i:lub_i,llb_j:lub_j,:)

   end subroutine

   ! Update rhs = rhs - M*sol after a Jacobi pass
   subroutine rhs_update_ij(sol, rhs, mat, weight, &
                           alb_i, aub_i, alb_j, aub_j, ge_k, &
                           lb_i, ub_i, lb_j, ub_j)
      ! Get matrix index references from matvec, only those used here
      use matvec, only : IDX_WEST, IDX_EAST, IDX_NORTH, IDX_SOUTH

      ! Parameters
      integer, intent(in) :: alb_i, aub_i, alb_j, aub_j, ge_k ! Allocated array bounds
      integer, intent(in) :: lb_i, ub_i, lb_j, ub_j ! Array subset

      real*4, intent(in) :: weight ! Underrelaxation weight used for the Jacobi pass

      real*4, intent(inout) :: rhs(alb_i:aub_i, alb_j:aub_j, ge_k) ! RHS
      real*4, intent(in)    :: sol(alb_i:aub_i, alb_j:aub_j, ge_k) ! Solution, for M*sol
      real*4, intent(in)    :: mat(alb_i:aub_i, alb_j:aub_j, ge_k, 7) ! Matrix

      ! Loop variables
      integer :: i, j, k

      ! Array alignment specifiers
!dir$ assume_aligned rhs:64, sol:64, mat:64

      ! After a Jacobi pass, we don't need to perform a full matrix-vector multiplication
      ! to update the RHS.

      ! The full M*u = f problem splits as (Mv + Mh)*u = f.  A full-weighted Jacobi relaxation
      ! solves u' = (Mv)^(-1) f, so f' = f - M*u' turns into f' = -Mh*u', requiring
      ! computation of only the horizontal terms (and considerable memory-bandwidth savings).

      ! In general, we apply an underweighted Jacobi relaxation, so that u' = w*(Mv)^(-1) f,
      ! so we need to account for that in the update.

      do k=1,ge_k
         do j=lb_j, ub_j
            do i = lb_i, ub_i
               rhs(i,j,k) = (1-weight)*rhs(i,j,k) - &
                              mat(i,j,k,IDX_WEST)*sol(i-1,j,k) - &
                              mat(i,j,k,IDX_EAST)*sol(i+1,j,k) - &
                              mat(i,j,k,IDX_SOUTH)*sol(i,j-1,k) - &
                              mat(i,j,k,IDX_NORTH)*sol(i,j+1,k)
            end do
         end do
      end do
   end subroutine

   ! Update rhs = rhs - M*sol, using all terms
   subroutine rhs_update_ijk(sol, rhs, mat, &
                             alb_i, aub_i, alb_j, aub_j, ge_k, &
                             lb_i, ub_i, lb_j, ub_j)
      ! Get matrix index references from matvec, only those used here
      use matvec, only : IDX_WEST, IDX_EAST, IDX_NORTH, IDX_SOUTH, IDX_POINT, IDX_TOP, IDX_BOTTOM

      ! Parameters
      integer, intent(in) :: alb_i, aub_i, alb_j, aub_j, ge_k ! Allocated array bounds
      integer, intent(in) :: lb_i, ub_i, lb_j, ub_j ! Array subset

      real*4, intent(inout) :: rhs(alb_i:aub_i, alb_j:aub_j, ge_k) ! RHS
      real*4, intent(in)    :: sol(alb_i:aub_i, alb_j:aub_j, ge_k) ! Solution, for M*sol
      real*4, intent(in)    :: mat(alb_i:aub_i, alb_j:aub_j, ge_k, 7) ! Matrix

      ! Loop variables
      integer :: i, j, k

      ! Array alignment specifiers
!dir$ assume_aligned rhs:64, sol:64, mat:64

      k = 1
      do j=lb_j, ub_j
         do i = lb_i, ub_i
            rhs(i,j,k) = rhs(i,j,k) - mat(i,j,k,IDX_WEST)  *sol(i-1,j,k) - &
                                      mat(i,j,k,IDX_EAST)  *sol(i+1,j,k) - &
                                      mat(i,j,k,IDX_SOUTH) *sol(i,j-1,k) - &
                                      mat(i,j,k,IDX_NORTH) *sol(i,j+1,k) - &
                                      mat(i,j,k,IDX_BOTTOM)*sol(i,j,k+1) - &
                                      mat(i,j,k,IDX_POINT) *sol(i,j,k)
         end do
      end do
      do k=2,ge_k-1
         do j=lb_j, ub_j
            do i = lb_i, ub_i
               rhs(i,j,k) = rhs(i,j,k) - mat(i,j,k,IDX_WEST)  *sol(i-1,j,k) - &
                                         mat(i,j,k,IDX_EAST)  *sol(i+1,j,k) - &
                                         mat(i,j,k,IDX_SOUTH) *sol(i,j-1,k) - &
                                         mat(i,j,k,IDX_NORTH) *sol(i,j+1,k) - &
                                         mat(i,j,k,IDX_TOP)   *sol(i,j,k-1) - &
                                         mat(i,j,k,IDX_BOTTOM)*sol(i,j,k+1) - &
                                         mat(i,j,k,IDX_POINT) *sol(i,j,k)
            end do
         end do
      end do
      k = ge_k
      do j=lb_j, ub_j
         do i = lb_i, ub_i
            rhs(i,j,k) = rhs(i,j,k) - mat(i,j,k,IDX_WEST)  *sol(i-1,j,k) - &
                                      mat(i,j,k,IDX_EAST)  *sol(i+1,j,k) - &
                                      mat(i,j,k,IDX_SOUTH) *sol(i,j-1,k) - &
                                      mat(i,j,k,IDX_NORTH) *sol(i,j+1,k) - &
                                      mat(i,j,k,IDX_TOP)   *sol(i,j,k-1) - &
                                      mat(i,j,k,IDX_POINT) *sol(i,j,k)
         end do
      end do
   end subroutine

   ! Perform one Jacobi pass on a subset of the arrays
   subroutine line_jacobi(rhs, sol, weight, scal, cp, m, &
                          alb_i, aub_i, alb_j, aub_j, ge_k, &
                          lb_i, ub_i, lb_j, ub_j)
      ! Get matrix index references from matvec, only those needed here
      use matvec, only : IDX_POINT, IDX_TOP, IDX_BOTTOM

      ! Parameters
      integer, intent(in) :: alb_i, aub_i, alb_j, aub_j, ge_k ! Allocated array bounds
      integer, intent(in) :: lb_i, ub_i, lb_j, ub_j ! Solved array subset

      real*4, intent(in) :: weight ! Underrelaxation weight

      real*4, intent(in) :: rhs(alb_i:aub_i, alb_j:aub_j, ge_k) ! Input RHS
      real*4, intent(out) :: sol(alb_i:aub_i, alb_j:aub_j, ge_k) ! Output LHS
      real*4, intent(in) :: scal(alb_i:aub_i, alb_j:aub_j, ge_k) ! Input matrix LU factors
      real*4, intent(in) :: cp(alb_i:aub_i, alb_j:aub_j, ge_k)   !
      real*4, intent(in) :: m(alb_i:aub_i, alb_j:aub_j, ge_k)    !

      ! Loop variables
      integer :: i, j, k
      integer :: si ! Slice length along i, for vectorization

      ! Working arrays for the Thomas algorithm; these are overwritten as
      ! part of the tridiagonal solve
      real*4, dimension(VECLEN,ge_k) :: ds

      ! Alignment directives
!dir$ attributes align: 64 :: ds
!dir$ assume_aligned rhs:64, sol:64, scal:64, cp:64, m:64

      ! This is the most work-intensive loop of the entire V-cycle, so
      ! we'll go to lengths to ensure it completely vectorizes.
      do j=lb_j, ub_j

         ! Stride the i-loop by VECLEN.  The line solve has irresolvable dependencies
         ! between k-levels, but multiple i-levels can be processed at the same time.
         ! The automatic vectorizer, however, only (easily) works on inner loops, so
         ! this is written as loop-over-k, then loop-over-stride.
         do i=lb_i,ub_i,VECLEN
            si = min(VECLEN, 1 + ub_i - i)

            ! Forward substitution pass
            ds(1:si,1) = rhs(i:(i+si-1),j,1)*scal(i:(i+si-1),j,1)
            do k=2,ge_k
               ds(1:si,k) = rhs(i:(i+si-1),j,k)*scal(i:(i+si-1),j,k) - &
                            m(i:(i+si-1),j,k)*ds(1:si,k-1)
            end do

            ! Backwards substitution: ds(:,ge_k) already contains the solution
            sol(i:(i+si-1),j,ge_k) = weight*ds(1:si,ge_k)
            do k=ge_k-1,1,-1
               ! Write directly to the sol vector, to avoid a ds->sol copy
               sol(i:(i+si-1),j,k) = weight*(ds(1:si,k)) - cp(i:(i+si-1),j,k)*sol(i:(i+si-1),j,k+1)
            end do

         end do
      end do
   end subroutine
end module multigrid_3d_jac
