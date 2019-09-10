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

!** fgmres - Flexible generalized minimum residual method (with restarts).
!
      subroutine sol_fgmres2d (solution, matvec, rhs_b, tolerance, maxinner, maxouter, nbiter, conv, level)
      use dyn_fisl_options
      use glb_ld
      use ldnh
      use prec
      use sol
      use redblack_2d

      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: level

      ! Initial guess on input, approximate solution on output
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(inout) :: solution

      ! A matrix-vector product routine (A.*v).
      interface
         subroutine matvec(v, prod, level)
            use ldnh, only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy
      use, intrinsic :: iso_fortran_env
            implicit none
            integer, intent(in) :: level
            real(kind=REAL64), dimension (ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(in) :: v
            real(kind=REAL64), dimension (ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(out) :: prod
         end subroutine
      end interface

      ! Right hand side of the linear system.
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(in) :: rhs_b

      ! Tolerance to achieve. The algorithm terminates when the relative
      ! residual is below tolerance.
      real(kind=REAL64), intent(in) :: tolerance

      ! Restarts the method every maxinner inner iterations.
      integer, intent(in) :: maxinner

      ! Specifies the maximum number of outer iterations.
      ! Iteration will stop after maxinner*maxouter steps
      ! even if the specified tolerance has not been achieved.
      integer, intent(in) :: maxouter

      ! Total number of iteration performed
      integer, intent(out) :: nbiter

      real(kind=REAL64), intent(out) :: conv

      ! References
      !
      ! C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, 1995
      ! (https://www.siam.org/books/textbooks/fr16_book.pdf)
      !
      ! Y. Saad, Iterative Methods for Sparse Linear Systems. SIAM, 2003.
      ! (http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf)
      !
      ! Katarzyna Świrydowicz, Julien Langou, Shreyas Ananthan, Ulrike Meier Yang, Stephen Thomas
      ! Low synchronization Gram-Schmidt and GMRES algorithms
      ! To appear in Numerical Linear Algebra with Applications (TODO: update reference when paper is published)
      !
      integer :: i, j, k, k1, ii, jj, ierr

      integer :: initer, outiter, nextit, it
      real(kind=REAL64) :: relative_tolerance, r0
      real(kind=REAL64) :: norm_residual, norm_b, nu, t
      real(kind=REAL64), dimension(maxinner+1, maxinner) :: hessenberg

      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy) :: work_space
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, maxinner) :: ww
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, maxinner+1) :: vv
      real(kind=REAL64), dimension(maxinner+1, maxinner+1) :: rr, tt

      real(kind=REAL64) :: local_dot
      real(kind=REAL64), dimension(:,:), allocatable :: v_local_prod, v_prod

      real(kind=REAL64), dimension(maxinner+1) :: rot_cos, rot_sin, gg

      logical almost_zero

      integer i0, in, j0, jn
      integer niloc,njloc

      niloc = (l_ni-pil_e)-(1+pil_w)+1
      njloc = (l_nj-pil_n)-(1+pil_s)+1

      ! Here we go !

      i0 = 1  + sol_pil_w
      in = l_ni - sol_pil_e
      j0 = 1  + sol_pil_s
      jn = l_nj - sol_pil_n

      outiter = 0
      nbiter = 0

      conv = 0.d0

      ! Residual of the initial iterate
      call matvec(solution, work_space, level)

      ! Compute ||b*b|| to determine the required error for convergence
      local_dot = 0.0d0
      do j=j0,jn
!DIR$ SIMD
         do i=i0,in
            local_dot = local_dot + (rhs_b(i,j)*rhs_b(i,j))
         end do
      end do
      call RPN_COMM_allreduce(local_dot, norm_b, 1, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)

      r0 = sqrt(norm_b)

      ! Scale tolerance according to the norm of b
      relative_tolerance = tolerance * r0

      do j=j0,jn
         do i=i0,in
            vv(i,j,1) = rhs_b(i,j) - work_space(i,j)
         end do
      end do

      do

         local_dot = 0.0d0
         do j=j0,jn
!DIR$ SIMD
            do i=i0,in
               local_dot = local_dot + (vv(i, j, 1) * vv(i, j, 1))
            end do
         end do

         call RPN_COMM_allreduce(local_dot, norm_residual, 1, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)
         norm_residual = sqrt(norm_residual)

         ! Current guess is a good enough solution
         if (norm_residual < relative_tolerance) then
            !call RPN_COMM_Finalize(i)
            !stop 'hammer time'
            return
         end if

         nu = 1.0d0 / norm_residual
         do j=j0,jn
            do i=i0,in
               vv(i,j,1) = vv(i,j,1) * nu
            end do
         end do

         conv = norm_residual / r0

         ! initialize 1-st term of rhs of hessenberg system.
         gg(1) = norm_residual
         gg(2:) = 0.d0

         tt = 0.d0
         rr = 0.d0

         initer = 0

         do

            nbiter = nbiter + 1
            initer = initer + 1
            nextit = initer + 1

            if (sol2D_precond_S == 'JACOBI')   then
               call pre_jacobi2D ( work_space(i0:in,j0:jn), &
                                   vv(i0:in,j0:jn, initer), &
                                   Prec_xevec_8, niloc, njloc,&
                                   Prec_ai_8, Prec_bi_8, Prec_ci_8, level )
            elseif (sol2D_precond_S == 'REDBLACK') then
               call pre_redblack2D (work_space(i0:in,j0:jn), &
                                    vv(i0:in,j0:jn,initer), &
                                    niloc, njloc, level)
            else
               work_space(i0:in,j0:jn) = vv(i0:in,j0:jn, initer)
            end if

            ww(i0:in,j0:jn, initer) = work_space(i0:in,j0:jn)


            call matvec ( work_space, vv(:,:,nextit), level )


            ! Modified Gram-Schmidt from Świrydowicz et al. (2018)

            allocate( v_local_prod(initer,2), v_prod(initer,2) )
            v_local_prod = 0.d0 ; v_prod = 0.d0

            do it=1,initer
               do j=j0,jn
!DIR$ SIMD
                  do i=i0,in
                     v_local_prod(it,1) = v_local_prod(it,1) + ( vv(i, j, it) * vv(i, j, initer) )
                     v_local_prod(it,2) = v_local_prod(it,2) + ( vv(i, j, it) * vv(i, j, nextit) )
                  end do
               end do
            end do

            call RPN_COMM_allreduce(v_local_prod, v_prod, initer*2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)

            tt(1:initer-1,initer) = v_prod(1:initer-1,1)
            rr(initer,initer)     = v_prod(initer,1)
            rr(1:initer,nextit)   = v_prod(1:initer,2)

            deallocate ( v_local_prod, v_prod)

            rr(initer,initer) = sqrt( rr(initer,initer) )
            vv(i0:in, j0:jn, initer)  = vv(i0:in, j0:jn, initer) / rr(initer,initer)

            rr(initer,nextit) = rr(initer,nextit) / rr(initer,initer)

            if (nextit > 2) then
               tt(1:initer-1, initer) = tt(1:initer-1, initer) / rr(initer,initer)
            end if

            tt(initer,initer) = 1.d0
            tt(1:initer-1, initer) = - matmul( tt(1:initer-1, 1:initer-1), tt(1:initer-1, initer) )
            rr(1:initer,nextit) = matmul( transpose(tt(1:initer, 1:initer)), rr(1:initer,nextit) )

            do it=1,initer
               do j=j0,jn
                  do i=i0,in
                     vv(i, j, nextit) = vv(i, j, nextit) - vv(i, j, it) * rr(it,nextit)
                  end do
               end do
            end do

            local_dot = 0.d0
            do j=j0,jn
!DIR$ SIMD
               do i=i0,in
                  local_dot = local_dot + (vv(i, j, nextit) * vv(i, j, nextit))
               end do
            end do

            call RPN_COMM_allreduce(local_dot,nu,1,"MPI_double_precision","MPI_sum","MULTIGRID",ierr)
            rr(nextit,nextit) = sqrt(nu)

            ! Watch out for happy breakdown
            if (.not. almost_zero( rr(nextit,nextit) ) ) then
               nu = 1.d0 / rr(nextit,nextit)
               do j=j0,jn
                  do i=i0,in
                     vv(i, j, nextit) = vv(i, j, nextit) * nu
                  end do
               end do
            end if

            hessenberg(1:nextit,initer) = rr(1:nextit,nextit)

            ! Form and store the information for the new Givens rotation
            if (initer > 1) then
               do k=2,initer
                  k1 = k-1
                  t = hessenberg(k1,initer)
                  hessenberg(k1,initer) = rot_cos(k1)*t + rot_sin(k1)*hessenberg(k,initer)
                  hessenberg(k,initer) = -rot_sin(k1)*t + rot_cos(k1)*hessenberg(k,initer)
               end do

            end if

            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
            if (.not. almost_zero(nu)) then
               nu = 1.d0 / nu
               rot_cos(initer) = hessenberg(initer,initer) * nu
               rot_sin(initer) = hessenberg(nextit,initer) * nu

               gg(nextit) = -rot_sin(initer) * gg(initer)
               gg(initer) =  rot_cos(initer) * gg(initer)

               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
            end if

            norm_residual = abs(gg(nextit))

            conv = norm_residual / r0

            if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            end if

         end do

         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve upper triangular system
         gg(initer) = gg(initer) / hessenberg(initer,initer)
         do ii=2,initer
            k  = initer - ii + 1
            k1 = k + 1
            t  = gg(k)
!$DIR SIMD
            do j=k1,initer
               t = t - hessenberg(k,j) * gg(j)
            end do
            gg(k) = t / hessenberg(k,k)
         end do

         ! Form linear combination to get solution.
         do it=1,initer
            t = gg(it)

            do j=j0,jn
               do i=i0,in
                  solution(i, j) = solution(i, j) + t * ww(i, j, it)
               end do
            end do

         end do

         outiter = outiter + 1

         if (norm_residual <= relative_tolerance .or. outiter > maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.
         do it=1,initer
            jj = nextit - it + 1
            gg(jj-1) = -rot_sin(jj-1) * gg(jj)
            gg(jj)   =  rot_cos(jj-1) * gg(jj)
         end do

         do it=1,nextit
            t = gg(it)
            if (it == 1) then
               t = t - 1.d0
            end if

            do j=j0,jn
!DIR$ SIMD
               do i=i0,in
                  vv(i, j, 1) = vv(i, j, 1) + t * vv(i, j, it)
               end do
            end do

         end do

      end do

      return
      end
