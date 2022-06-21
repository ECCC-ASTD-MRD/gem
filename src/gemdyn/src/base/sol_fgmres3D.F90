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
      subroutine sol_fgmres3d (solution, rhs_b, tolerance, maxinner, maxouter, nbiter, conv)
      use dyn_fisl_options
      use dynkernel_options
      use glb_ld
      use ldnh
      use sol
      use lam_options
      use omp_lib
      use mem_tstp
      use prec
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      include 'mpif.h'
      include 'rpn_comm.inc'
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
      ! K. Swirydowicz, J. Langou, S. Ananthan, U. Yang, and S. Thomas:
      ! Low synchronization Gram–Schmidt and GMRES algorithms. United States: N. p., 2020. Web. doi:10.1002/nla.2343.
      ! (https://onlinelibrary.wiley.com/doi/10.1002/nla.2343)

      integer :: i, j, k, k0, k1, ii, jj, ierr

      integer :: initer, outiter, it
      real(kind=REAL64) ::   t
      real(kind=REAL64) :: local_dot(2), glb_dot(2), l_avg_8(2), r_avg_8(2)
      real(kind=REAL64), dimension(maxinner+1,2) :: v_local_prod, v_prod, v_avg_8

      real(kind=REAL64) :: rrp

      logical almost_zero

      integer i0, in, j0, jn
      integer niloc,njloc, nii, njj, nk

      nii = Sol_iin-Sol_ii0+1
      njj = Sol_jjn-Sol_jj0+1

      niloc = (l_ni-pil_e)-(1+pil_w)+1
      njloc = (l_nj-pil_n)-(1+pil_s)+1
      k0=1+Lam_gbpil_T
      nk=l_nk-k0+1

      i0 = 1  + sol_pil_w
      in = l_ni - sol_pil_e
      j0 = 1  + sol_pil_s
      jn = l_nj - sol_pil_n

      outiter = 0
      nbiter = 0

      conv = 0.d0

!$omp do
            do it=1,maxinner+1
               do k=1,l_nk
                  do j=sol_jmin,sol_jmax
                     do i=sol_imin,sol_imax
                        vv(i,j,k,it) = 0.d0
                        work_space(i,j,k) = 0.d0
                     end do
                  end do
               end do
            enddo
!$omp end do
!$omp do
            do it=1,maxinner+1
               do k=1,l_nk
                  do j=sol_jj0,sol_jjn
                     do i=sol_ii0,sol_iin
                        wint_8(i,j,k,it) = 0.d0
                     end do
                  end do
               end do
            enddo
!$omp end do

      ! Residual of the initial iterate
      call matvec3D ( solution, work_space,ldnh_minx,ldnh_maxx, ldnh_miny,ldnh_maxy, &
                      sol_imin,sol_imax,sol_jmin,sol_jmax)

      !  Compute ||b*b|| to determine the required error for convergence

      local_dot(1)=0.d0
!$omp do
      do k=k0,l_nk
         do j=j0,jn
            do i=i0,in
               vv(i,j,k,1) = rhs_b(i,j,k) - work_space(i,j,k)
            end do
         end do
      end do
!$omp enddo

!$omp do
      do k=k0,l_nk
         do j=j0,jn
            do i=i0,in
               local_dot(1) = local_dot(1) + (rhs_b(i,j,k)*rhs_b(i,j,k))
            end do
         end do
      end do
!$omp enddo
      thread_s(1,OMP_get_thread_num()) = local_dot(1)
!$OMP BARRIER

      do
      local_dot(2)=0.d0
!$omp single
         do k=k0,l_nk
            do j=j0,jn
               do i=i0,in
                  local_dot(2) = local_dot(2) + (vv(i, j, k, 1) * vv(i, j, k, 1))
               end do
            end do
         end do
!$omp end single
         thread_s(2,OMP_get_thread_num()) = local_dot(2)
!$OMP BARRIER

!$omp single
         l_avg_8(1) = sum(thread_s(1,:))
         l_avg_8(2) = sum(thread_s(2,:))

         call RPN_COMM_allreduce(l_avg_8, glb_dot, 2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)
         r0 = sqrt(glb_dot(1))
         norm_residual = sqrt(glb_dot(2))
!$omp end single
!$OMP BARRIER

         ! Scale tolerance according to the norm of b;
         relative_tolerance = tolerance * r0

         ! Current guess is a good enough solution
         if (norm_residual < relative_tolerance) then
            return
         end if

         nu = 1.0d0 / norm_residual
         conv = norm_residual / r0

!$omp do
         do k=k0,l_nk
            do j=j0,jn
               do i=i0,in
                  vv(i,j,k,1) = vv(i,j,k,1) * nu
               end do
            end do
         end do
!$omp enddo

         ! initialize 1-st term of rhs of hessenberg system.
         gg(1) =  norm_residual
         gg(2:) = 0.d0

         tt = 0.
         rr = 0.

         do initer=1,maxinner

            nbiter = nbiter + 1

            ! Preconditionning: Restricted Additive Schwarz (RAS)
!$omp single
            call rpn_comm_xch_halo_8 (vv(:,:,k0:l_nk,initer),Sol_imin,Sol_imax,Sol_jmin,Sol_jmax, &
                                      l_ni,l_nj,nk,ovlpx,ovlpy,G_periodx,G_periody,l_ni,0)
!$omp end single
            call pre_jacobi3D2 (wint_8(Sol_ii0:Sol_iin,Sol_jj0:Sol_jjn,:,initer), &
                                    vv(Sol_ii0:Sol_iin,Sol_jj0:Sol_jjn,:,initer), nii,njj,l_nk)

            !Matrix-Vector product
            call matvec3D (wint_8(:,:,:,initer),vv(:,:,:,initer+1),Sol_ii0,Sol_iin,Sol_jj0,Sol_jjn, &
                           Sol_imin,Sol_imax,Sol_jmin,Sol_jmax)

            ! Modified Gram-Schmidt from Świrydowicz et al. (2018)
                  v_local_prod  = 0.d0
                  v_lcl_sum = 0.d0
                  v_prod   = 0.d0

!$omp do
            do it=1,initer+1
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        v_local_prod(it,1) = v_local_prod(it,1) + ( vv(i, j, k, it) * vv(i, j, k, initer) )
                        v_local_prod(it,2) = v_local_prod(it,2) + ( vv(i, j, k, it) * vv(i, j, k, initer+1) )
                     end do
                  end do
               end do
            enddo
!$omp enddo

            do it =1,initer+1
               thread_s2(1,it,OMP_get_thread_num()) = v_local_prod(it,1)
               thread_s2(2,it,OMP_get_thread_num()) = v_local_prod(it,2)
            enddo
!$omp BARRIER

!$omp single
            do it =1,initer+1
               v_avg_8(it,1) = sum(thread_s2(1,it,:))
            enddo
            do it =1,initer+1
               v_avg_8(it,2) = sum(thread_s2(2,it,:))
            enddo
            call RPN_COMM_allreduce(v_avg_8(1:initer+1,:), v_prod(1:initer+1,:), (initer+1)*2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)

            do it=1,initer-1
               tt(it,initer) = v_prod(it,1)
            enddo

            do it=1,initer
               rr(it,initer+1) = v_prod(it,2)
            enddo

            rr(initer,initer) = v_prod(initer,1)
            rr(initer,initer) = sqrt( rr(initer,initer) )
            rr(initer,initer+1) = rr(initer,initer+1) / rr(initer,initer)
            ro2 = v_prod(initer+1,2)
!$omp end single
!$omp BARRIER



!$omp do
            do k=k0,l_nk
               do j=j0,jn
                  do i=i0,in
                     vv(i, j, k, initer)  = vv(i, j, k, initer) / rr(initer,initer)
                  enddo
               enddo
            enddo
!$omp enddo


!$omp single
            if (initer > 1) then
               do it=1,initer-1
                  tt(it, initer) = tt(it, initer) / rr(initer,initer)
               enddo
            end if

            tt(initer,initer) = 1.d0
            tt(1:initer-1, initer) = - matmul( tt(1:initer-1, 1:initer-1), tt(1:initer-1, initer) )
            rr(1:initer,initer+1) = matmul( transpose(tt(1:initer, 1:initer)), rr(1:initer,initer+1) )
!$omp end single

            ! Compute rr2=rr(:,initer+1)*rr(:,initer+1) needed in the computation of vv(:,:,:,initer+1) estimated norm
            rrp=0.d0 ; rr2=0.d0
!$omp do
            do i=1,initer
                rrp = rrp + rr(i,initer+1)*rr(i,initer+1)
            enddo
!$omp enddo
            thread_s(4,OMP_get_thread_num()) = rrp
!$omp BARRIER
!$omp single
            rr2  = sum(thread_s(4,:))
!$omp end single

            do it=1,initer
!$omp do
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        vv(i, j, k, initer+1) = vv(i, j, k, initer+1) - vv(i, j, k, it) * rr(it,initer+1)
                     end do
                  end do
               end do
!$omp enddo
            end do

            ! Compute estimated norm of vv(:,:,:,initer+1): from GHYSELS et al. (Journal of Scientific Computing 2013)
              if( (ro2-rr2) > 0.) then
!$omp single
               nu=sqrt(ro2-rr2)
               rr(initer+1,initer+1) =  nu
!$omp end single
            ! Or in case (ro2-rr2)<= 0 compute exact norm
            else

               local_dot(1) = 0.d0; lcl_sum=0.d0
!$omp do
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        local_dot(1) = local_dot(1) + (vv(i, j, k, initer+1) * vv(i, j, k, initer+1))
                     end do
                  end do
               end do
!$omp enddo
            thread_s(3,OMP_get_thread_num()) = local_dot(1)
!$OMP BARRIER

!$omp single
            r_avg_8 = sum(thread_s(3,:))
            call RPN_COMM_allreduce(r_avg_8,nu,1,"MPI_double_precision","MPI_sum","MULTIGRID",ierr)
            rr(initer+1,initer+1) = sqrt(nu)
!$omp end single
!$OMP BARRIER
            endif

            if ( .not. almost_zero( rr(initer+1,initer+1) ) ) then
               nu = 1.d0 / rr(initer+1,initer+1)
!$omp do
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        vv(i, j, k, initer+1) = vv(i, j, k, initer+1) * nu
                     end do
                  end do
               end do
!$omp enddo

            end if

!$omp single
            do it=1,initer+1
               hessenberg(it,initer) = rr(it,initer+1)
            enddo

            ! Form and store the information for the new Givens rotation
            if (initer > 1) then
               do k=2,initer
                  k1 = k-1
                  t = hessenberg(k1,initer)
                  hessenberg(k1,initer) = rot_cos(k1)*t + rot_sin(k1)*hessenberg(k,initer)
                  hessenberg(k,initer) = -rot_sin(k1)*t + rot_cos(k1)*hessenberg(k,initer)
               end do
            end if

            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(initer+1,initer)**2)

            if ( .not. almost_zero(nu) ) then
               nu = 1.d0 / nu
               rot_cos(initer) = hessenberg(initer,initer) * nu
               rot_sin(initer) = hessenberg(initer+1,initer) * nu

               gg(initer+1) = -rot_sin(initer) * gg(initer)
               gg(initer) =  rot_cos(initer) * gg(initer)

               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(initer+1,initer)

            end if
!$omp end single

            norm_residual = abs(gg(initer+1))

            conv = norm_residual / r0
            if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            endif
         end do


         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve upper triangular system
!$omp single
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
!$omp end single

         ! Form linear combination to get solution.
         do it=1,initer

!$omp do
            do k=k0,l_nk
               do j=j0,jn
!DIR$ SIMD
                  do i=i0,in
                     solution(i, j, k) = solution(i, j, k) + gg(it)* wint_8(i, j, k, it)
                  end do
               end do
            end do
!$omp enddo

         end do

         outiter = outiter + 1


         if (norm_residual <= relative_tolerance .or. outiter >= maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.
!$omp  single
         do it=1,initer
            jj = initer+1 - it + 1
            gg(jj-1) = -rot_sin(jj-1) * gg(jj)
            gg(jj)   =  rot_cos(jj-1) * gg(jj)
         end do
!$omp end single

         do it=1,initer+1
            t = gg(it)
            if (it == 1) then
               t = t - 1.d0
            end if

!$omp do
            do k=k0,l_nk
               do j=j0,jn
!DIR$ SIMD
                  do i=i0,in
                     vv(i, j, k, 1) = vv(i, j, k, 1) + t * vv(i, j, k, it)
                  end do
               end do
            end do
!$omp enddo

         end do

     end do

      return
      end
