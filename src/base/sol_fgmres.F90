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

!** fgmres - Flexible generalized minimum residual method(with restarts).

      subroutine sol_fgmres ( solution, rhs_b, F_print_L )
      use dynkernel_options
      use glb_ld
      use ldnh
      use lun
      use sol_mem
      use sol_options
      use lam_options
      use omp_lib
      use mem_tstp
      use prec
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      include 'mpif.h'
      include 'rpn_comm.inc'
      ! Initial guess on input, approximate solution on output
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(inout) :: solution

      ! Right hand side of the linear system.
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(in) :: rhs_b

      logical, intent(in) :: F_print_L

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

      logical almost_zero

      integer :: i, j, k, k1, ii, jj, ierr
      integer :: initer, outiter, it, nbiter

      real(kind=REAL64) :: t, conv, local_dot(2), glb_dot(2), l_avg_8(2), r_avg_8(2)
      real(kind=REAL64), dimension(sol_im+1,2) :: v_local_prod, v_prod, v_avg_8
      real(kind=REAL64) :: rrp, wrr2, wro2, wnu , wr0, residual, Rel_tolerance
!
!     ---------------------------------------------------------------
!
      outiter= 0 ; nbiter= 0 ; conv= 0.d0

      ! Residual of the initial iterate
      call matvec_hlt ( solution  , ldnh_minx,ldnh_maxx,ldnh_miny,ldnh_maxy,&
                        work_space, sol_imin,sol_imax,sol_jmin,sol_jmax, l_nk )

      !  Compute ||b*b|| to determine the required error for convergence

      local_dot(1)=0.d0
!!$omp master
!$omp do collapse(2)
      do k=Sol_k0,l_nk
         do j=Sol_j0,Sol_jn
            do i=Sol_i0,Sol_in
               vv(i,j,k,1)  = rhs_b(i,j,k) - work_space(i,j,k)
               local_dot(1) = local_dot(1) + (rhs_b(i,j,k)*rhs_b(i,j,k))
            end do
         end do
      end do
!$omp enddo
!!$omp end master

      do ! MAIN OUTER ITERATION
         local_dot(2)=0.d0
!!$omp master
!$omp do collapse(2)
         do k=Sol_k0,l_nk
            do j=Sol_j0,Sol_jn
               do i=Sol_i0,Sol_in
                  local_dot(2) = local_dot(2) + (vv(i, j, k, 1) * vv(i, j, k, 1))
               end do
            end do
         end do
!$omp enddo
!!$omp end master
         thread_s(1:2,OMP_get_thread_num()) = local_dot(1:2)
!$OMP BARRIER

!$omp single
         l_avg_8(1) = sum(thread_s(1,:))
         l_avg_8(2) = sum(thread_s(2,:))
         call RPN_COMM_allreduce(l_avg_8, glb_dot, 2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)
         wr0      = sqrt(glb_dot(1))
         residual = sqrt(glb_dot(2))

!$omp end single copyprivate(wr0, residual)

!         wr0      = sqrt(glb_dot(1))
!         residual = sqrt(glb_dot(2))

         ! Scale tolerance according to the norm of b;
         Rel_tolerance = sol_fgm_eps * wr0
         conv = residual / wr0
         if ((OMP_get_thread_num() == 0).and.F_print_L) &
         write(Lun_out, "(3x,'FGMRES convergence at iteration', i4,' relative residual=',1pe14.7)") nbiter, conv

         ! Current guess is a good enough solution
         if (residual < Rel_tolerance) return

         wnu = 1.0d0 / residual

!$omp do collapse(2)
         do k=Sol_k0,l_nk
            do j=Sol_j0,Sol_jn
               do i=Sol_i0,Sol_in
                  vv(i,j,k,1) = vv(i,j,k,1) * wnu
               end do
            end do
         end do
!$omp enddo

         ! initialize 1-st term of rhs of hessenberg system.
         gg(1 ) = residual
         gg(2:) = 0.d0

         tt = 0. ; rr = 0.

         do initer=1,sol_im

            nbiter = nbiter + 1

            ! Preconditionning: Restricted Additive Schwarz (RAS)
!$omp single
            call rpn_comm_xch_halo_8 (vv(sol_imin,sol_jmin,Sol_k0,initer),Sol_imin,Sol_imax,Sol_jmin,Sol_jmax, &
                                      l_ni,l_nj,Sol_nk,Sol_ovlpx,Sol_ovlpy,G_periodx,G_periody,l_ni,0)
!$omp end single
            call pre_jacobi_hlt (wint_8(Sol_ii0:Sol_iin,Sol_jj0:Sol_jjn,:,initer), &
                                     vv(Sol_ii0:Sol_iin,Sol_jj0:Sol_jjn,:,initer), sol_niloc,sol_njloc,l_nk)

            !Matrix-Vector product
            call matvec_hlt (wint_8(sol_ii0,sol_jj0,1,initer),Sol_ii0,Sol_iin,Sol_jj0,Sol_jjn,&
                             vv(sol_imin,sol_jmin,1,initer+1),Sol_imin,Sol_imax,Sol_jmin,Sol_jmax,l_nk)

            ! Modified Gram-Schmidt from Świrydowicz et al. (2018)
            v_local_prod  = 0.d0
            v_lcl_sum     = 0.d0
            v_prod        = 0.d0

            do it=1,initer+1
!!$omp master
!$omp do collapse(2)
               do k=Sol_k0,l_nk
                  do j=Sol_j0,Sol_jn
                     do i=Sol_i0,Sol_in
                        v_local_prod(it,1) = v_local_prod(it,1) + ( vv(i, j, k, it) * vv(i, j, k, initer  ) )
                        v_local_prod(it,2) = v_local_prod(it,2) + ( vv(i, j, k, it) * vv(i, j, k, initer+1) )
                     end do
                  end do
               end do
!$omp enddo
!!$omp end master
               thread_s2(1:2,it,OMP_get_thread_num()) = v_local_prod(it,1:2)
            enddo
!$omp BARRIER

!$omp single
            do it =1,initer+1
               v_avg_8(it,1) = sum(thread_s2(1,it,:))
               v_avg_8(it,2) = sum(thread_s2(2,it,:))
            enddo
            call RPN_COMM_allreduce(v_avg_8(1:initer+1,:), v_prod(1:initer+1,:), (initer+1)*2,&
                                    "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)
            do it=1,initer-1
               tt(it,initer) = v_prod(it,1)
            enddo
            do it=1,initer
               rr(it,initer+1) = v_prod(it,2)
            enddo
            rr(initer,initer) = v_prod(initer,1)
            rr(initer,initer) = sqrt( rr(initer,initer) )
            rr(initer,initer+1) = rr(initer,initer+1) / rr(initer,initer)

            wro2 = v_prod(initer+1,2)
!$omp end single copyprivate(wro2)

!$omp do collapse(2)
            do k=Sol_k0,l_nk
               do j=Sol_j0,Sol_jn
                  do i=Sol_i0,Sol_in
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

            ! Compute wrr2=rr(:,initer+1)*rr(:,initer+1) needed in the computation of vv(:,:,:,initer+1) estimated norm
            rrp=0.d0 ; wrr2=0.d0
!$omp do
            do i=1,initer
               rrp = rrp + rr(i,initer+1)*rr(i,initer+1)
            enddo
!$omp enddo
            thread_s(4,OMP_get_thread_num()) = rrp
!$omp BARRIER

            wrr2  = sum(thread_s(4,:))

            do it=1,initer
!$omp do collapse(2)
               do k=Sol_k0,l_nk
                  do j=Sol_j0,Sol_jn
                     do i=Sol_i0,Sol_in
                        vv(i, j, k, initer+1) = vv(i, j, k, initer+1) - vv(i, j, k, it) * rr(it,initer+1)
                     end do
                  end do
               end do
!$omp enddo
            end do

            ! Compute estimated norm of vv(:,:,:,initer+1): from GHYSELS et al. (Journal of Scientific Computing 2013)
            if( (wro2-wrr2) > 0.) then
!$omp single
               wnu=sqrt(wro2-wrr2)
               rr(initer+1,initer+1) =  wnu
!$omp end single
            ! Or in case (wro2-wrr2)<= 0 compute exact norm
            else
               
               local_dot(1) = 0.d0; lcl_sum=0.d0
!!$omp master
!$omp do collapse(2)
               do k=Sol_k0,l_nk
                  do j=Sol_j0,Sol_jn
                     do i=Sol_i0,Sol_in
                        local_dot(1) = local_dot(1) + (vv(i, j, k, initer+1) * vv(i, j, k, initer+1))
                     end do
                  end do
               end do
!$omp enddo
!!$omp end master
               thread_s(3,OMP_get_thread_num()) = local_dot(1)
!$OMP BARRIER

!$omp single
               r_avg_8 = sum(thread_s(3,:))
               call RPN_COMM_allreduce(r_avg_8,wnu,1,"MPI_double_precision","MPI_sum","MULTIGRID",ierr)
               rr(initer+1,initer+1) = sqrt(wnu)
!$omp end single
            endif

            if ( .not. almost_zero( rr(initer+1,initer+1) ) ) then
               wnu = 1.d0 / rr(initer+1,initer+1)
!$omp do collapse(2)
               do k=Sol_k0,l_nk
                  do j=Sol_j0,Sol_jn
                     do i=Sol_i0,Sol_in
                        vv(i, j, k, initer+1) = vv(i, j, k, initer+1) * wnu
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

            wnu = sqrt(hessenberg(initer,initer)**2 + hessenberg(initer+1,initer)**2)

            if ( .not. almost_zero(wnu) ) then
               wnu = 1.d0 / wnu
               rot_cos(initer) = hessenberg(initer,initer) * wnu
               rot_sin(initer) = hessenberg(initer+1,initer) * wnu

               gg(initer+1) = -rot_sin(initer) * gg(initer)
               gg(initer) =  rot_cos(initer) * gg(initer)

               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(initer+1,initer)
            end if
!$omp end single

            residual = abs(gg(initer+1))

            conv = residual / wr0
            if ((OMP_get_thread_num() == 0).and.F_print_L) &
            write(Lun_out, "(3x,'FGMRES convergence at iteration', i4,' relative residual=',1pe14.7)") nbiter, conv
            if ((initer >= sol_im) .or. (residual <= Rel_tolerance)) exit
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

!$omp do collapse(2)
            do k=Sol_k0,l_nk
               do j=Sol_j0,Sol_jn
!DIR$ SIMD
                  do i=Sol_i0,Sol_in
                     solution(i, j, k) = solution(i, j, k) + gg(it)* wint_8(i, j, k, it)
                  end do
               end do
            end do
!$omp enddo

         end do

         outiter = outiter + 1

         if (residual <= Rel_tolerance .or. outiter >= sol_fgm_maxits) return

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

!$omp do collapse(2)
            do k=Sol_k0,l_nk
               do j=Sol_j0,Sol_jn
!DIR$ SIMD
                  do i=Sol_i0,Sol_in
                     vv(i, j, k, 1) = vv(i, j, k, 1) + t * vv(i, j, k, it)
                  end do
               end do
            end do
!$omp enddo
         end do

      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine sol_fgmres
