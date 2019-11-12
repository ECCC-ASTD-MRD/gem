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
      subroutine sol_fgmres3d_H (solution, rhs_b, tolerance, maxinner, maxouter, nbiter, conv)
      use dyn_fisl_options
      use glb_ld
      use ldnh
      use sol
      use ptopo
      use prec             ! Qaddouri-- Blockwise Red/Black Gauss-Seidel and jacobi preconditioners
      use redblack_3d      ! csubich -- Red/Black (z-column) preconditioner
      use multigrid_3d_jac ! csubich -- Multigrid (column relaxation) preconditioner
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      ! Initial guess on input, approximate solution on output
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(inout) :: solution

      ! Right hand side of the linear system.
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(in) :: rhs_b

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
      real(kind=REAL64) :: norm_residual, norm_Ax_2, b_Ax, hegedus_scaling, nu, t
      real(kind=REAL64), dimension(maxinner+1, maxinner) :: hessenberg

      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk) :: work_space
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk, maxinner) :: ww
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk, maxinner+1) :: vv
      real(kind=REAL64), dimension(maxinner+1, maxinner+1) :: rr, tt

      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk) :: wint_8
      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk) :: wint_82

      real(kind=REAL64) :: local_dot
      real(kind=REAL64), dimension(3) :: initial_dot, initial_local_dot
      real(kind=REAL64), dimension(:,:), allocatable :: v_local_prod, v_prod

      real(kind=REAL64), dimension(maxinner+1) :: rot_cos, rot_sin, gg

      logical almost_zero

      integer i0, in, j0, jn
      integer niloc,njloc

      niloc = (l_ni-pil_e)-(1+pil_w)+1
      njloc = (l_nj-pil_n)-(1+pil_s)+1
      wint_82 = 0.d0

      ! Here we go !

      i0 = 1  + sol_pil_w
      in = l_ni - sol_pil_e
      j0 = 1  + sol_pil_s
      jn = l_nj - sol_pil_n

      outiter = 0
      nbiter = 0

      conv = 0.d0

      vv = 0.d0

      ! Residual of the initial iterate
      call mat_vecs3D_H3 ( solution, work_space, ldnh_minx,ldnh_maxx, ldnh_miny,ldnh_maxy,l_ni,l_nj,l_nk)

!      call matvec(solution, work_space)

      ! Index 1 : Compute ||b*b|| to determine the required error for convergence
      ! Index 2 : Compute b^T*Ax for Hegedus trick
      ! Index 3 : Compute ||Ax||^2 for Hegedus trick
      initial_local_dot = 0.0d0
      do k=1,l_nk
         do j=j0,jn
            do i=i0,in
               initial_local_dot(1) = initial_local_dot(1) + (rhs_b(i,j,k)*rhs_b(i,j,k))
               initial_local_dot(2) = initial_local_dot(2) + (rhs_b(i,j,k)*work_space(i,j,k))
               initial_local_dot(3) = initial_local_dot(3) + (work_space(i,j,k)*work_space(i,j,k))
            end do
         end do
      end do

      call RPN_COMM_allreduce(initial_local_dot, initial_dot, 3, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)

      r0        = sqrt( initial_dot(1) )
      b_Ax      = initial_dot(2)
      norm_Ax_2 = initial_dot(3)

      ! Scale tolerance according to the norm of b
      relative_tolerance = tolerance * r0

      ! Rescale initial guess appropriately (Hegedüs trick)
      if ( .not. almost_zero( norm_Ax_2 ) ) then
         hegedus_scaling = b_Ax / norm_Ax_2
      else
         hegedus_scaling = 1.d0
      end if

      do k=1,l_nk
         do j=j0,jn
            do i=i0,in
               solution(i,j,k) = solution(i,j,k) * hegedus_scaling

               vv(i,j,k,1) = rhs_b(i,j,k) - work_space(i,j,k) * hegedus_scaling
            end do
         end do
      end do

      do

         local_dot = 0.0d0
         do k=1,l_nk
            do j=j0,jn
!DIR$ SIMD
               do i=i0,in
                  local_dot = local_dot + (vv(i, j, k, 1) * vv(i, j, k, 1))
               end do
            end do
         end do

         call RPN_COMM_allreduce(local_dot, norm_residual, 1, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)
         norm_residual = sqrt(norm_residual)

         ! Current guess is a good enough solution
         if (norm_residual < relative_tolerance) then
            return
         end if

         nu = 1.0d0 / norm_residual
         do k=1,l_nk
            do j=j0,jn
               do i=i0,in
                  vv(i,j,k,1) = vv(i,j,k,1) * nu
               end do
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

            select case(sol3D_precond_S)
               case ('JACOBI') ! Blockwise Jacobi preconditioner
                  wint_8(i0:in,j0:jn,:)=vv(i0:in,j0:jn,:,initer)

                  call pre_jacobi3D ( wint_82(i0:in,j0:jn,:), &
                                      vv(i0:in,j0:jn,:,initer), &
                                      Prec_xevec_8, niloc, njloc, l_nk,&
                                      Prec_ai_8, Prec_bi_8, Prec_ci_8 )
                  work_space(i0:in,j0:jn,:)=wint_82(i0:in,j0:jn,:)
               case ('GAUSS') ! Blockwise Red/Black Gauss-Seidel preconditioner
                  wint_8(i0:in,j0:jn,:)=vv(i0:in,j0:jn,:,initer)
                  call RB_BGAUSS(wint_82,wint_8,l_minx,l_maxx,l_miny,l_maxy,l_nk)
                  work_space(i0:in,j0:jn,:)=wint_82(i0:in,j0:jn,:)
               case ('REDBLACK') ! Column-wise Red/Black Gauss-Seidel preconditioner
                  ! As an interface note, pre_redblack3D and pre_multigrid_jac3d both take
                  ! full arrays as references; they handle the borders internally.
                  call pre_redblack3D (work_space, vv(:,:,:,initer))
               case ('MULTIGRID') ! Column relaxation multigrid preconditioner
                  call pre_multigrid_jac3d (work_space, vv(:,:,:,initer))
               case default
                  work_space(i0:in,j0:jn,:) = vv(i0:in,j0:jn,:,initer)
            end select

            ww(i0:in,j0:jn,:,initer) = work_space(i0:in,j0:jn,:)

!            call matvec ( work_space, vv(:,:,:,nextit) )
             call mat_vecs3D_H3 (  work_space,vv(:,:,:,nextit), ldnh_minx,ldnh_maxx, ldnh_miny,ldnh_maxy,l_ni,l_nj,l_nk)



            ! Modified Gram-Schmidt from Świrydowicz et al. (2018)

            ! TODO : avoid memory allocation
            allocate( v_local_prod(initer,2), v_prod(initer,2) )
            v_local_prod = 0.d0 ; v_prod = 0.d0

            do it=1,initer
               do k=1,l_nk
                  do j=j0,jn
!DIR$ SIMD
                     do i=i0,in
                        v_local_prod(it,1) = v_local_prod(it,1) + ( vv(i, j, k, it) * vv(i, j, k, initer) )
                        v_local_prod(it,2) = v_local_prod(it,2) + ( vv(i, j, k, it) * vv(i, j, k, nextit) )
                     end do
                  end do
               end do
            end do

            call RPN_COMM_allreduce(v_local_prod, v_prod, initer*2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)

            tt(1:initer-1,initer) = v_prod(1:initer-1,1)
            rr(initer,initer)     = v_prod(initer,1)
            rr(1:initer,nextit)   = v_prod(1:initer,2)

            deallocate ( v_local_prod, v_prod)

            rr(initer,initer) = sqrt( rr(initer,initer) )
            vv(i0:in, j0:jn, 1:l_nk, initer)  = vv(i0:in, j0:jn, 1:l_nk, initer) / rr(initer,initer)

            rr(initer,nextit) = rr(initer,nextit) / rr(initer,initer)

            if (nextit > 2) then
               tt(1:initer-1, initer) = tt(1:initer-1, initer) / rr(initer,initer)
            end if

            tt(initer,initer) = 1.d0
            tt(1:initer-1, initer) = - matmul( tt(1:initer-1, 1:initer-1), tt(1:initer-1, initer) )
            rr(1:initer,nextit) = matmul( transpose(tt(1:initer, 1:initer)), rr(1:initer,nextit) )

            do it=1,initer
               do k=1,l_nk
                  do j=j0,jn
                     do i=i0,in
                        vv(i, j, k, nextit) = vv(i, j, k, nextit) - vv(i, j, k, it) * rr(it,nextit)
                     end do
                  end do
               end do
            end do

            local_dot = 0.d0
            do k=1,l_nk
               do j=j0,jn
!DIR$ SIMD
                  do i=i0,in
                     local_dot = local_dot + (vv(i, j, k, nextit) * vv(i, j, k, nextit))
                  end do
               end do
            end do

            call RPN_COMM_allreduce(local_dot,nu,1,"MPI_double_precision","MPI_sum","MULTIGRID",ierr)
            rr(nextit,nextit) = sqrt(nu)

            ! Watch out for happy breakdown
            if ( .not. almost_zero( rr(nextit,nextit) ) ) then
               nu = 1.d0 / rr(nextit,nextit)

               do k=1,l_nk
                  do j=j0,jn
                     do i=i0,in
                        vv(i, j, k, nextit) = vv(i, j, k, nextit) * nu
                     end do
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
            if ( .not. almost_zero(nu) ) then
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
            do j=k1,initer
               t = t - hessenberg(k,j) * gg(j)
            end do
            gg(k) = t / hessenberg(k,k)
         end do

         ! Form linear combination to get solution.
         do it=1,initer
            t = gg(it)

            do k=1,l_nk
               do j=j0,jn
                  do i=i0,in
                     solution(i, j, k) = solution(i, j, k) + t * ww(i, j, k, it)
                  end do
               end do
            end do

         end do

         outiter = outiter + 1

         if (norm_residual <= relative_tolerance .or. outiter >= maxouter) then
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

            do k=1,l_nk
               do j=j0,jn
!DIR$ SIMD
                  do i=i0,in
                     vv(i, j, k, 1) = vv(i, j, k, 1) + t * vv(i, j, k, it)
                  end do
               end do
            end do

         end do

      end do

      return
      end
