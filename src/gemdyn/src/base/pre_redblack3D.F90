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

! pre_redblack3D -- Preconditioner for the 3D elliptic problem, via implicit
!                   treatment of vertical variations and Red/Black Gauss-
!                   Seidel in the horizontal.  Uses the matrix operations
!                   defined in the matvec module
module redblack_3d

contains

   subroutine pre_redblack3D ( Lhs, Rhs )

      use matvec ! Elliptic operator, plus indexing parameters
      ! Import grid allocation dimensions from ldnh
      use ldnh, only   : Ldnh_minx, Ldnh_maxx, Ldnh_miny, Ldnh_maxy
      ! Import logical grid information from grd_ld
      use glb_ld, only : l_ni, l_nj, l_nk, l_i0, l_j0, &
                         G_periodx, G_periody
      ! Import solver pilot region from sol
      use sol, only    : sol_pil_e, sol_pil_w, sol_pil_n, sol_pil_s
      use glb_pil

      use, intrinsic :: iso_fortran_env
      implicit none

      ! Input (Rhs) and output (Lhs) arrays have dimensions given by the ldnh variables
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(in)  :: Rhs
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(out) :: Lhs

      ! The Red/Black Gauss-seidel algorithm performs line relaxation on the red points first,
      ! then it uses that approximate solution to give an enhanced solution at the black points.
      ! This requires exchange of the red points across processors.  By defining separate working
      ! arrays for x(red) and f(black), we minimize data exchange and enhance contiguous data access
      real, dimension(:,:,:), allocatable, save :: lhs_red_halo, rhs_black
!dir$ attributes align: 64:: lhs_red_halo, rhs_black
      integer, save :: cl_ni

      ! Similarly, it is very convenient to have local copies of the matrix, similarly-ordered
      ! such that entries for red and black points are separated
      real, save, dimension(:,:,:,:), allocatable :: mat_red, mat_black
!dir$ attributes align: 64:: mat_red, mat_black
      logical, save :: matrix_init = .false.

      ! The vertical line solve proceeds via the Thomas algorithm, which is a simplified forward/
      ! backward substitution.  This destructively modifies the coefficients of the tridiagonal
      ! matrix and the right-hand side, so we want arrays to hold temporary copies.
      integer, parameter :: STRIP_SIZE = 16
      real, dimension(STRIP_SIZE,l_nk) :: as, bs, cs, ds ! Arrays for "strip-mined" Thomas Algorithm
      integer :: strip


      ! Finally, we need indexing variables
      integer :: i, j, k ! Loop variables

      ! This code ultimately uses five coordinate frames, an unfortunate number called for
      ! because the red/black colouring is global and we also want to take advantage of
      ! monocoloured arrays for data locality.  These reference frames are (with unused
      ! variables commented out):

      !integer :: li      ! Local I-index, into Lhs and Rhs arrays
      !integer :: gi      ! Global I-index, into the global domain inclusive of the pilot region
      integer :: gi_0    ! Global index of the leftmost array point
      integer :: di      ! Global ``in-domain'' I-index, which excludes glb_pil_w
      integer :: di_0    ! Domain I-index of the leftmost point
      !integer :: cgi     ! Global colour I-index, into the imaginary global monochromatic array
      integer :: cgi_0   ! Global colour I-index of the leftmost point
      integer :: ci      ! Local colour I-index
      integer :: cimin, cimax ! Colour I-indices to control iteration

      ! Some of these are also necessary for the j-index
      !integer :: gj      ! Global J-index
      integer:: gj_0 ! Global J-index of the southmost point
      integer :: dj      ! In-domain J-index
      integer :: dj_0 ! In-domain J-index of the southmost point

      integer :: ioffset ! Offset for red/black indexing

      gi_0 = l_i0
      di_0 = gi_0 - glb_pil_w
      cgi_0 = (di_0 + 1)/2
      gj_0 = l_j0
      dj_0 = gj_0 - glb_pil_s

      ! As general commentary, this function preserves the (x,y,z) logical array ordering
      ! used elsewhere in GEM dynamics, but this might not be the most appropriate.  Coupling
      ! is tightest within a z-line, so a (z, x, y) internal array ordering might make more sense.
      ! This should be explored only after profiling and/or examining the failures of auto-
      ! vectorization.

      if (.not.matrix_init) then
         matrix_init = .true.

         ! Split the matrix into red and black parts, for better data locality in the subsequent loops.
         ! Additionally, transpose it such that all coefficients of a type (such as POINT, or WEST, or TOP)
         ! are together in memory, again for locality
         cl_ni = (di_0 + (l_ni - 1) + 1)/2 - cgi_0 + 1
         !cl_ni_alloc = 8*((cl_ni)/8+1)
         allocate (mat_red(1:cl_ni,1:l_nj,1:l_nk,7))
         allocate (mat_black(1:cl_ni,1:l_nj,1:l_nk,7))

         do k=1,l_nk
            do j=1+sol_pil_s,l_nj-sol_pil_n
               dj  = j + dj_0 - 1
               do i=1+sol_pil_w,l_ni-sol_pil_e
                  di = i + di_0 - 1
                  ci = (di + 1)/2 - cgi_0 + 1

                  if (mod(di+dj,2) == 0) then ! Black point
                     mat_black(ci,j,k,:) = stencil(i,j,:,k)
                  else ! Red point
                     mat_red(ci,j,k,:) =  stencil(i,j,:,k)
                  end if
               end do
            end do
         end do
         ! Allocate the colour arrays.  For reduced data storage, this code separates
         ! lhs_red (red points computed after the first step) and rhs_black (the rhs
         ! of black points updated with new estimates for red), which involve
         ! approximately Ni/2 points along the first dimension.  However, this depends
         ! on how the grid splits between processors, so compute this in the global
         ! reference frame
         allocate (lhs_red_halo(0:(cl_ni+1), 0:(l_nj+1), l_nk))
         allocate (rhs_black(cl_ni,1:l_nj,l_nk))

         lhs_red_halo = 0
         rhs_black = 0

      end if

      ! Iterative over ``red'' points, defined by mod(ii+jj,2) == 1 where ii and jj are GLOBAL indices
      do j = 1+sol_pil_s, l_nj-sol_pil_n
         dj = j + dj_0 - 1 ! Domain j-index of this column
         di = sol_pil_w + di_0 ! Domain i-index of the leftmost point of the inner iteration

         ! The point (di,dj) is a red point if (di + dj) % 2 = 1, so compute an offset in I (possibly 0)
         ! such that 1+sol_pil_w+ioffset is the first (local) red point of the row
         ioffset = 1 - mod(dj+(1+sol_pil_w+di_0-1),2)

         do i = 1+sol_pil_w+ioffset, l_ni-sol_pil_e,2*STRIP_SIZE
            ! Implement the Thomas algorithm: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
            ! This should be numerically stable because the diagonal entry of the 3D operator (IDX_POINT)
            ! contains stabilizing terms from the timestep and from the horizontal Laplacian

            di = i + di_0 - 1             ! Domain I-index of this point
            ci = (di + 1) / 2 - cgi_0 + 1 ! I-index of this point in the red-only array


            ! For better vectorization, implement this in a ``strip-mined'' way.  The algorithm as-written
            ! cannot be vectorized over a line because of data dependencies, but we -can- handle multiple
            ! lines at the same time.

            strip = min(STRIP_SIZE,1+(l_ni - sol_pil_e - i)/2)

            ! Copy variables to a, b, c, d
            ! As a lexical trap, the vertical matrix entries count positive as downward.  That is,
            ! at level k, IDX_BOTTOM refers to k+1 and IDX_TOP refers to k-1
            ds(1:strip,:) = Rhs(i:(i+2*strip-2):2,j,:)               ! Right-hand side
            cs(1:strip,:) = mat_red(ci:(ci+strip-1),j,:,IDX_BOTTOM) ! Superdiagonal
            bs(1:strip,:) = mat_red(ci:(ci+strip-1),j,:,IDX_POINT) ! Main diagonal
            as(1:strip,:) = mat_red(ci:(ci+strip-1),j,:,IDX_TOP) ! Subdiagonal

            ! Forward sweep, k=1 special case
            cs(1:strip,1) = cs(1:strip,1) / bs(1:strip,1)
            ds(1:strip,1) = ds(1:strip,1) / bs(1:strip,1)
            do k = 2, l_nk ! Forward sweep
               cs(1:strip,k) = cs(1:strip,k) / (bs(1:strip,k) - as(1:strip,k)*cs(1:strip,k-1))
               ds(1:strip,k) = (ds(1:strip,k) - as(1:strip,k)*ds(1:strip,k-1)) / (bs(1:strip,k) - as(1:strip,k)*cs(1:strip,k-1))
            end do

            ! Backwards sweep, to store the solution in d
            ! d(k) already contains the solution
            do k = l_nk-1, 1, -1
               ds(1:strip,k) = ds(1:strip,k) - cs(1:strip,k)*ds(1:strip,k+1)
            end do

            ! Write d (containing the solution) to Lhs and Lhs_halo, the latter
            ! because boundary values will need to be exchanged
            !Lhs_halo(i:(i+2*strip-2):2,j,:) = ds(1:strip,:)
            lhs_red_halo(ci:(ci+strip-1),j,:) = ds(1:strip,:)
            Lhs(i:(i+2*strip-2):2,j,:) = ds(1:strip,:)

         end do
      end do

      ! Perform a halo exchange
      call rpn_comm_xch_halon (lhs_red_halo, 0, cl_ni+1, 0, l_nj+1, cl_ni, l_nj, l_nk, &
                               1, 1, G_periodx, G_periody, cl_ni, 0, 1)



      ! Update the RHS at 'black' points based on the newly-computed LHS at
      ! 'red' points'
      do k= 1, l_nk
         do j = 1+sol_pil_s, l_nj-sol_pil_n
            dj = j + dj_0 - 1 ! Domain J-indix of this column

            if (mod(dj,2) == 0) then
               ! If the J-column is even, then (ODD,j) is a red point and (EVEN,j) is a black point.
               ! The conversion from global-index to colour-index divides by two, so its inverse might
               ! be offset.  The conversion is:

               ! di(red)   = 2*cgi - 1
               ! di(black) = 2*cgi

               ! Now, we need to know which local i-index corresponds to the first black point.  The offset
               ! is 1 if the domain-index of 1+sol_pil_w (what would be the start of the iteration) is odd,
               ! so...
               ioffset = mod(1+sol_pil_w+di_0-1,2)

               ! Finally, we can define iteration bounds on the colour array, by transforming to the global
               ! colour index [(i + di_0-1)/2] and back
               cimin = (1+sol_pil_w+ioffset+di_0)/2 - cgi_0 + 1
               ! Calculate the maximum by counting the number of black points along the line
               cimax = ((l_ni-sol_pil_e) - (1+sol_pil_w+ioffset))/2 + cimin
               do ci = cimin, cimax
                  i = 2*(ci + cgi_0 - 1) - di_0 + 1
                  rhs_black(ci,j,k) = Rhs(i,j,k) - mat_black(ci,j,k,IDX_WEST)*lhs_red_halo(ci,j,k) &
                                                 - mat_black(ci,j,k,IDX_EAST)*lhs_red_halo(ci+1,j,k) &
                                                 - mat_black(ci,j,k,IDX_SOUTH)*lhs_red_halo(ci,j-1,k) &
                                                 - mat_black(ci,j,k,IDX_NORTH)*lhs_red_halo(ci,j+1,k)
               end do

            else
               ! Otherwise, the J-column is odd, and (EVEN,j) is a red point and (ODD,j) is a black point.
               ! The conversion from global-index to colour-index is now:

               ! di(red) = 2*cgi
               ! di(black) = 2*cgi - 1

               ! The offset for our first i-index black point is the inverse of the above, or
               ioffset = 1-mod(1+sol_pil_w+di_0-1,2)

               ! And the iteration bounds on the colour array are
               cimin = (1+sol_pil_w+ioffset+di_0)/2 - cgi_0 + 1
               cimax = ((l_ni-sol_pil_e) - (1+sol_pil_w+ioffset))/2 + cimin
               do ci = cimin, cimax
                  i = 2*(ci + cgi_0 - 1) - di_0
                  rhs_black(ci,j,k) = Rhs(i,j,k) - mat_black(ci,j,k,IDX_WEST)*lhs_red_halo(ci-1,j,k)  &
                                                 - mat_black(ci,j,k,IDX_EAST)*lhs_red_halo(ci,j,k)  &
                                                 - mat_black(ci,j,k,IDX_SOUTH)*lhs_red_halo(ci,j-1,k) &
                                                 - mat_black(ci,j,k,IDX_NORTH)*lhs_red_halo(ci,j+1,k)
               end do
            end if

         end do
      end do

      ! Now, iterate over ``black'' points
      do j = 1+sol_pil_s, l_nj-sol_pil_n
         dj = j + l_j0 - 1 - glb_pil_s ! Global j-index of this column, with respect to the piloting region
         di = sol_pil_w + l_i0 - glb_pil_w ! Global i-index of the first (westmost) row
         ! The point (di,dj) is a black point if (di + dj) % 2 = 0, so compute an offset in I (possibly 0)
         ! such that 1+sol_pil_w+ioffset is the first (local) black point of the row
         ioffset = mod(di+dj,2)

         do i = 1+sol_pil_w+ioffset, l_ni-sol_pil_e,2*STRIP_SIZE
            ! Again, implement the Thomas algorithm
            di = i + di_0 - 1             ! Domain I-index of this point
            ci = (di + 1) / 2 - cgi_0 + 1 ! I-index of this point in the black-only array
            strip = min(STRIP_SIZE,1+(l_ni - sol_pil_e - i)/2)

            ! Copy variables to a, b, c, d
            ds(1:strip,:) = rhs_black(ci:(ci+strip-1),j,:)
            cs(1:strip,:) = mat_black(ci:(ci+strip-1),j,:,IDX_BOTTOM) ! Superdiagonal
            bs(1:strip,:) = mat_black(ci:(ci+strip-1),j,:,IDX_POINT) ! Main diagonal
            as(1:strip,:) = mat_black(ci:(ci+strip-1),j,:,IDX_TOP) ! Subdiagonal

            ! Forward sweep, k=1 special case
            cs(1:strip,1) = cs(1:strip,1) / bs(1:strip,1)
            ds(1:strip,1) = ds(1:strip,1) / bs(1:strip,1)
            do k = 2, l_nk ! Forward sweep
               cs(1:strip,k) = cs(1:strip,k) / (bs(1:strip,k) - as(1:strip,k)*cs(1:strip,k-1))
               ds(1:strip,k) = (ds(1:strip,k) - as(1:strip,k)*ds(1:strip,k-1)) / (bs(1:strip,k) - as(1:strip,k)*cs(1:strip,k-1))
            end do

            ! Backwards sweep, to store the solution in d
            ! d(k) already contains the solution
            do k = l_nk-1, 1, -1
               ds(1:strip,k) = ds(1:strip,k) - cs(1:strip,k)*ds(1:strip,k+1)
            end do

            ! Write d (containing the solution) to Lhs
            Lhs(i:(i+2*strip-2):2,j,:) = ds(1:strip,:)
         end do
      end do
      !STOP 1


   end subroutine pre_redblack3D
end module redblack_3d
