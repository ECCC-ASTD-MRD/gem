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
!**s/r sol_2d - Elliptic solver based on vertical decomposition leading
!               to F_nk 2D horizontal elliptic problems to solve.
!
      subroutine sol_direct ( F_rhs_8, F_sol_8, F_ni, F_nj, F_nk, &
                              F_print_L, F_offi, F_offj )
      use gem_options
      use HORgrid_options
      use dyn_fisl_options
      use glb_ld
      use lun
      use ldnh
      use sol_options
      use sol_mem
      use opr
      use ptopo
      use trp
      use adz_mem
      use gem_timing
      use, intrinsic :: iso_fortran_env
      implicit none

      logical, intent(in) :: F_print_L
      integer, intent(in) :: F_ni,F_nj,F_nk,F_offi, F_offj
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(in) :: F_rhs_8
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(inout) :: F_sol_8
!
      integer i,j,k,ni,nij,iter,n
      real linfini
      real(kind=REAL64), dimension (ldnh_maxx,ldnh_maxy,l_nk) :: wk3
      real(kind=REAL64), dimension (l_minx:l_maxx,l_miny:l_maxy,l_nk) :: yyrhs
!
!     ---------------------------------------------------------------
!
      ni  = ldnh_ni-pil_w-pil_e
      nij = (ldnh_maxy-ldnh_miny+1)*(ldnh_maxx-ldnh_minx+1)

      n= (Adz_icn-1)*Schm_itnlh + Adz_itnl
      call time_trace_barr(gem_time_trace, 1, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

!$omp parallel private (i,j,k) shared (F_offi,F_offj,ni,nij,Sol_rhs_8)
!$omp do
      do j=1+pil_s, ldnh_nj-pil_n
         call dgemm ('N','N', ni, G_nk, G_nk, 1.d0, F_rhs_8(1+pil_w,j,1), &
                     nij, Opr_lzevec_8, G_nk, 0.d0, Sol_rhs_8(1+pil_w,j,1), nij)
         do k=1,Schm_nith
            do i = 1+pil_w, ldnh_ni-pil_e
               Sol_rhs_8(i,j,k) = Opr_opsxp0_8(G_ni+F_offi+i) * &
                                  Opr_opsyp0_8(G_nj+F_offj+j) * Sol_rhs_8(i,j,k)
            end do
         end do
      end do
!$omp enddo
!$omp end parallel
      call time_trace_barr(gem_time_trace, 11000+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

      if (Grd_yinyang_L) then

         wk3 = Sol_rhs_8

         do iter=1, Sol_yyg_maxits

            if (sol_one_transpose_L) then
               call sol_fft_numa ( Sol_sol_8, Sol_rhs_8, &
                                   Sol_fft,Sol_a,Sol_b,Sol_c )
            else
               call sol_fft_lam ( Sol_sol_8, Sol_rhs_8, &
                               ldnh_maxx, ldnh_maxy, ldnh_nj,           &
                               trp_12smax, trp_12sn, trp_22max, trp_22n,&
                               G_ni, G_nj, G_nk, sol_nk,                &
                               Sol_ai_8, Sol_bi_8, Sol_ci_8 )
            endif
                            
            Sol_rhs_8 = wk3
            yyrhs(1:l_ni,1:l_nj,:) = -Sol_sol_8(1:l_ni,1:l_nj,:)
            call rpn_comm_xch_halo_8 (yyrhs,l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

            call yyg_rhs_xchng (Sol_rhs_8, yyrhs, ldnh_minx, ldnh_maxx   ,&
                                ldnh_miny, ldnh_maxy, l_minx,l_maxx,&
                                l_miny,l_maxy, l_nk, iter, linfini)

            if (Lun_debug_L.and.F_print_L) then
               write(Lun_out,1001) linfini,iter
            end if

            if ((iter > 1) .and. (linfini < Sol_yyg_eps)) then
               exit
            end if
         end do

         if (F_print_L) then
            write(Lun_out,1002) linfini,iter

            if (linfini > Sol_yyg_eps) then
               write(Lun_out,9001) Sol_yyg_eps
            end if
         end if

      else
         
         if (sol_one_transpose_L) then
            !if (Lun_out>0) print*, "SOL_FFT_NUMA"
            call sol_fft_numa ( Sol_sol_8, Sol_rhs_8, &
                    Sol_fft,Sol_a,Sol_b,Sol_c )
         else
            !if (Lun_out>0) print*, "SOL_FFT_LAM"
            call sol_fft_lam ( Sol_sol_8, Sol_rhs_8, &
                 ldnh_maxx, ldnh_maxy, ldnh_nj,           &
                 trp_12smax, trp_12sn, trp_22max, trp_22n,&
                 G_ni, G_nj, G_nk, sol_nk,                &
                 Sol_ai_8, Sol_bi_8, Sol_ci_8 )
         endif

      end if

!$omp parallel private (j) shared (g_nk, Sol_sol_8)
!$omp do
      do j=1+pil_s, ldnh_nj-pil_n
         call dgemm ('N','T', ni, G_nk, G_nk, 1.d0, Sol_sol_8(1+pil_w,j,1), &
                     nij, Opr_zevec_8, G_nk, 0.d0, F_sol_8(1+pil_w,j,1), nij)
      end do
!$omp enddo
!$omp end parallel

      call time_trace_barr(gem_time_trace, 11800+n, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

 1001 format (3x,'Iterative YYG    solver convergence criteria: ',1pe14.7,' at iteration', i3)
 1002 format (3x,'Final YYG    solver convergence criteria: ',1pe14.7,' at iteration', i3)
 9001 format (3x,'WARNING: iterative YYG solver DID NOT converge to requested criteria:: ',1pe14.7)
!
!     ---------------------------------------------------------------
!
      return
      end
