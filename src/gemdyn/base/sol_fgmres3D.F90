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
      use lam_options
      use ldnh
      use sol
      use opr
      use ptopo
      use gem_timing
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

!Authors: A. Qdddouri & R. Aider  (2021)

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

      integer :: i, j, k, k0, k1, nk, ii, jj, ierr
      logical, save :: first_time = .true.

      integer :: initer, outiter, nextit, it
      real(kind=REAL64) :: relative_tolerance, r0
      real(kind=REAL64) :: norm_residual, nu, t
      real(kind=REAL64), dimension(maxinner+1, maxinner) :: hessenberg

      real(kind=REAL64), dimension(1-ovlpx:l_ni+ovlpx,1-ovlpy:l_nj+ovlpy, l_nk) :: work_space
      real(kind=REAL64), dimension(1-ovlpx:l_ni+ovlpx,1-ovlpy:l_nj+ovlpy, l_nk, maxinner+1) :: vv
      real(kind=REAL64), dimension(1-ovlpx:l_ni+ovlpx,1-ovlpy:l_nj+ovlpy, l_nk,maxinner) :: wint_82
      real(kind=REAL64), dimension(maxinner+1, maxinner+1) :: rr, tt

      real(kind=REAL64) :: local_dot(2), glb_dot(2)
      real(kind=REAL64), dimension(:,:), allocatable :: v_local_prod, v_prod

      real(kind=REAL64), dimension(maxinner+1) :: rot_cos, rot_sin, gg

      logical almost_zero

      real(kind=REAL64) :: ro2,rr2, temp
      integer i0, in, j0, jn, ii0, iin, jj0, jjn
      integer imin, imax, jmin, jmax
      integer niloc,njloc, ni_val, nj_val

      k0=1+Lam_gbpil_T
      nk=l_nk-k0+1

      niloc = (l_ni-pil_e)-(1+pil_w)+1
      njloc = (l_nj-pil_n)-(1+pil_s)+1

      i0 = 1  + sol_pil_w
      in = l_ni - sol_pil_e
      j0 = 1  + sol_pil_s
      jn = l_nj - sol_pil_n

      ii0  = 1    - ovlpx
      iin  = l_ni + ovlpx
      jj0  = 1    - ovlpy
      jjn  = l_nj + ovlpy

      imin = ii0
      imax = iin
      jmin = jj0
      jmax = jjn

      if (Ptopo_mycol==1)  ii0  = 0
      if (Ptopo_mycol.eq.Ptopo_npex-2)  iin = l_ni+1
      if (Ptopo_myrow==1)  jj0  = 0
      if (Ptopo_myrow.eq.Ptopo_npey-2)  jjn = l_nj+1

      if (l_east)  iin = l_ni-sol_pil_e
      if (l_west)  ii0 = 1+sol_pil_w
      if (l_north) jjn = l_nj-sol_pil_n
      if (l_south) jj0 = 1+sol_pil_s
      ni_val = iin-ii0+1
      nj_val = jjn-jj0+1

      outiter = 0
      nbiter = 0

      conv = 0.d0

      vv = 0.d0
      wint_82 = 0.d0

      ! Residual of the initial iterate
      ! if(first_time) then
      !   call matvec3d_init()
      !   first_time=.false.
      !endif

      call matvec3D ( solution, work_space,ldnh_minx,ldnh_maxx, ldnh_miny,ldnh_maxy)

      ! Index 1 : Compute ||b*b|| to determine the required error for convergence

      local_dot(1) = 0.0d0
      do k=k0,l_nk
         do j=j0,jn
            do i=i0,in
               local_dot(1) = local_dot(1) + (rhs_b(i,j,k)*rhs_b(i,j,k))
            end do
         end do
      end do
      do k=k0,l_nk
         do j=j0,jn
            do i=i0,in
               vv(i,j,k,1) = rhs_b(i,j,k) - work_space(i,j,k)
            end do
         end do
      end do

      do

         local_dot(2) = 0.0d0
         do k=k0,l_nk
            do j=j0,jn
               do i=i0,in
                  local_dot(2) = local_dot(2) + (vv(i, j, k, 1) * vv(i, j, k, 1))
               end do
            end do
         end do

         call RPN_COMM_allreduce(local_dot, glb_dot, 2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)


         r0 = sqrt(glb_dot(1))
         norm_residual = sqrt(glb_dot(2))

         ! Scale tolerance according to the norm of b
         relative_tolerance = tolerance * r0

         ! Current guess is a good enough solution
         if (norm_residual < relative_tolerance) then
            return
         end if

         nu = 1.0d0 / norm_residual
         do k=k0,l_nk
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
! Preconditionning RAS
             call rpn_comm_xch_halo_8 (vv(:,:,k0:l_nk,initer) ,   &
                                      imin, imax, jmin, jmax, l_ni, l_nj,nk, &
                                      ovlpx, ovlpy, G_periodx, G_periody, l_ni,0 )

            call pre_jacobi3D2 (wint_82(ii0:iin,jj0:jjn,:,initer),vv(ii0:iin,jj0:jjn,:,initer), &
                                        ii0,iin,jj0,jjn,ni_val,nj_val,l_nk)
!!

            call matvec3D (wint_82(:,:,:,initer),vv(:,:,:,nextit),imin,imax,jmin,jmax)

!            call mat_vecs3D_H3 ( wint_82(:,:,:,initer),vv(:,:,:,nextit),imin,imax, jmin, jmax ,l_ni,l_nj,l_nk)

            ! Modified Gram-Schmidt from Świrydowicz et al. (2018)

            ! TODO : avoid memory allocation
            allocate( v_local_prod(initer+1,3), v_prod(initer+1,3) )
            v_local_prod = 0.d0 ; v_prod = 0.d0

            do it=1,initer+1
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        v_local_prod(it,1) = v_local_prod(it,1) + ( vv(i, j, k, it) * vv(i, j, k, initer) )
                        v_local_prod(it,2) = v_local_prod(it,2) + ( vv(i, j, k, it) * vv(i, j, k, nextit) )
                     end do
                  end do
               end do
            end do

            call RPN_COMM_allreduce(v_local_prod, v_prod, nextit*2, "MPI_double_precision", "MPI_sum", "MULTIGRID", ierr)

            tt(1:initer-1,initer) = v_prod(1:initer-1,1)
            rr(initer,initer)     = v_prod(initer,1)
            rr(1:initer,nextit)   = v_prod(1:initer,2)
            ro2=v_prod(nextit,2)

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

            rr2=0.d0
            do i=1,initer
                rr2=rr2 + rr(i,nextit)*rr(i,nextit)
            enddo

            do it=1,initer
               temp=rr(it,nextit)
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        vv(i, j, k, nextit) = vv(i, j, k, nextit) - vv(i, j, k, it) * temp
                     end do
                  end do
               end do
            end do

            ! Compute estimated norm nu=sqrt(ro2-rr2)
            if( (ro2-rr2) > 0.) then
               nu=sqrt(ro2-rr2)
               rr(nextit,nextit) =  nu
            ! Compute exact norm
            else
               local_dot = 0.d0
               do k=k0,l_nk
                  do j=j0,jn
                     do i=i0,in
                        local_dot = local_dot + (vv(i, j, k, nextit) * vv(i, j, k, nextit))
                     end do
                  end do
               end do
            call RPN_COMM_allreduce(local_dot,nu,1,"MPI_double_precision","MPI_sum","MULTIGRID",ierr)
            rr(nextit,nextit) = sqrt(nu)
            endif

            if ( .not. almost_zero( rr(nextit,nextit) ) ) then
               nu = 1.d0 / rr(nextit,nextit)

               do k=k0,l_nk
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

            do k=k0,l_nk
               do j=j0,jn
!DIR$ SIMD
                  do i=i0,in
                     solution(i, j, k) = solution(i, j, k) + t * wint_82(i, j, k, it)
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

            do k=k0,l_nk
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
