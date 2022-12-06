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
!

! TODO : give some love to this s/r ...

   integer function fbicgstab3d(x, matvec, b, x0, ni, nj, nk, &
                                minx, maxx, miny, maxy, i0, il, j0, jl, &
                                tolerance, maxiter, precond_S, conv) result(retval)
      use glb_ld
      use lun
      use prec
      use, intrinsic :: iso_fortran_env
      implicit none
      ! Solves a linear system using a Bi-Conjugate Gradient STABilised
      !
      ! References
      ! C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, 1995
      ! (https://www.siam.org/books/textbooks/fr16_book.pdf)
      !
      integer, intent(in) :: ni, nj, nk, minx, maxx, miny, maxy, i0, il, j0, jl
      !
      ! The converged solution.
      real(kind=REAL64), dimension (minx:maxx, miny:maxy, nk), intent(out) :: x
      !
      ! Initial guess for the solution
      real(kind=REAL64), dimension (minx:maxx, miny:maxy, nk), intent(in) :: x0
      !
      ! Right hand side of the linear system.
      real(kind=REAL64), dimension (minx:maxx, miny:maxy, nk), intent(in) :: b
      !
      ! A matrix-vector product routine (A.*v).
      interface
         subroutine matvec(v, prod)
            use glb_ld, only: l_nk
            use ldnh, only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy
      use, intrinsic :: iso_fortran_env
            implicit none
            real(kind=REAL64), dimension (ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(in) :: v
            real(kind=REAL64), dimension (ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(out) :: prod
         end subroutine
      end interface
      !
      ! Tolerance to achieve. The algorithm terminates when either the relative
      ! or the absolute residual is below tolerance.
      real(kind=REAL64), intent(in) :: tolerance
      !
      ! Maximum number of iterations. Iteration will stop after maxiter steps
      ! even if the specified tolerance has not been achieved.
      integer, intent(in) :: maxiter
      ! Preconditioner
      character(len=*), intent(in) :: precond_S
      !
      real(kind=REAL64), intent(out) :: conv

      integer :: iter, i, j, k
      real(kind=REAL64) :: relative_tolerance, alpha, beta, tau, omega
      real(kind=REAL64) :: rho_old, rho, norm_residual, r0
      logical :: force_restart

      real(kind=REAL64), dimension(minx:maxx, miny:maxy, nk) :: residual, hatr0 ,res_matvec
      real(kind=REAL64), dimension(minx:maxx, miny:maxy, nk) :: vv, pp, pp_prec, ss, ss_prec, tt

      logical almost_zero
      real(kind=REAL64) dist_dotproduct

      ! Here we go !

      retval = 0

      x(:, :, :) = x0(:, :, :)

      residual(:, :, :) = 0.0d0
      ss(:, :, :) = 0.0d0
      tt(:, :, :) = 0.0d0
      vv(:, :, :) = 0.0d0
      pp(:, :, :) = 0.0d0

      ! Residual of the initial iterate
      call matvec(x, res_matvec)

      do k=1,nk
         do j=j0,jl
            do i=i0,il
               residual(i, j, k) = b(i, j, k) - res_matvec(i, j, k)
            end do
         end do
      end do

      force_restart = .false.
      hatr0 = residual
      rho_old = 1.d0
      alpha = 1.d0
      omega = 1.d0
      rho = dist_dotproduct(hatr0, residual, minx, maxx, miny, maxy, i0, il, j0, jl, nk)
      norm_residual = sqrt(rho)
      r0 = norm_residual

      ! Scale tolerance acording to the values of the residual
      relative_tolerance = tolerance * norm_residual

      iter = 0
      do while ((norm_residual > relative_tolerance) .and. (iter < maxiter))
         iter = iter + 1

         if (almost_zero(omega)) then
            write(lun_out, *) 'WARNING : BiCGSTAB breakdown (omega=0)'
            force_restart = .true.

            beta = 0.0d0
         else
            beta = (rho / rho_old) * (alpha / omega)
         end if

         do k=1,nk
            do j=j0,jl
               do i=i0,il
                  pp(i, j, k) = residual(i, j, k) + beta * (pp(i, j, k) - omega * vv(i, j, k))
               end do
            end do
         end do

         select case(precond_S)
            case ('JACOBI')
               call pre_jacobi3D (pp_prec(i0:il,j0:jl,:), pp(i0:il,j0:jl,:), Prec_xevec_8, &
                                  (Ni-pil_e)-(1+pil_w)+1, (Nj-pil_n)-(1+pil_s)+1, nk, &
                                  Prec_ai_8,Prec_bi_8,Prec_ci_8)
            case default
               pp_prec(i0:il,j0:jl,:) = pp(i0:il,j0:jl,:)
         end select

         ! Compute search direction
         call matvec(pp_prec, vv)
         tau = dist_dotproduct(hatr0, vv, minx, maxx, miny, maxy, i0, il, j0, jl, nk)

         if (almost_zero(tau)) then
            write(lun_out, *) 'WARNING : BiCGSTAB breakdown (tau=0)'
            force_restart = .true.

            alpha = 0.0d0
         else
            alpha = rho / tau
         end if

         do k=1,nk
            do j=j0,jl
               do i=i0,il
                  ss(i, j, k) = residual(i, j, k) - alpha * vv(i, j, k)
               end do
            end do
         end do

         select case(precond_S)
            case ('JACOBI')
                  call pre_jacobi3D (ss_prec(i0:il,j0:jl,:), ss(i0:il,j0:jl,:), Prec_xevec_8, &
                                     (Ni-pil_e)-(1+pil_w)+1, (Nj-pil_n)-(1+pil_s)+1, nk, &
                                     Prec_ai_8,Prec_bi_8,Prec_ci_8)
            case default
               ss_prec(i0:il,j0:jl,:) = ss(i0:il,j0:jl,:)
         end select

         call matvec(ss_prec, tt)
         tau = dist_dotproduct(tt, tt, minx, maxx, miny, maxy, i0, il, j0, jl, nk)

         if (almost_zero(tau)) then
            write(lun_out, *) 'WARNING : BiCGSTAB breakdown (tt=0)'
            force_restart = .true.

            omega = 0.0d0
         else
            omega = dist_dotproduct(tt, ss, minx, maxx, miny, maxy, i0, il, j0, jl, nk) / tau
         end if

         rho_old = rho
         rho = -omega * (dist_dotproduct(hatr0, tt, minx, maxx, miny, maxy, i0, il, j0, jl, nk))

         do k=1,nk
            do j=j0,jl
               do i=i0,il
                  ! Update the solution and the residual vectors
                  x(i, j, k) = x(i, j, k) + alpha * pp_prec(i, j, k) + omega * ss_prec(i, j, k)
                  residual(i, j, k) = ss(i, j, k) - omega * tt(i, j, k)
               end do
            end do
         end do

         norm_residual = sqrt(dist_dotproduct(residual, residual, minx, maxx, miny, maxy, &
                                               i0, il, j0, jl, nk))

         if (force_restart) then
            call matvec(x, res_matvec)

            do k=1,nk
               do j=j0,jl
                  do i=i0,il
                     residual(i, j, k) = b(i, j, k) - res_matvec(i, j, k)
                  end do
               end do
            end do

            hatr0 = residual
            rho_old = 1.d0
            alpha = 1.d0
            omega = 1.d0
            rho = dist_dotproduct(hatr0, residual, minx, maxx, miny, maxy, i0, il, j0, jl, nk)
            norm_residual = sqrt(rho)

            vv(:, :, :) = 0.d0
            pp(:, :, :) = 0.d0

            force_restart = .false.
         end if
      end do

       conv = norm_residual / r0
      retval = iter

   end function fbicgstab3d
