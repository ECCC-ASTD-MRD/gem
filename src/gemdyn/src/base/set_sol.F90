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

!**s/r set_sol - Computes matrices a,b,c for the elliptic solver

      subroutine set_sol
      use HORgrid_options
      use glb_pil
      use lam_options
      use dyn_fisl_options
      use glb_ld
      use cstv
      use lun
      use sol
      use ldnh
      use ptopo
      use gmm_itf_mod
      use opr
      use prec
      use trp
      use omp_lib
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      type(gmm_metadata) :: savemeta
      integer k,ni,nj,istat
      real(kind=REAL64) yg_8(G_nj), wk(G_nk)
!     __________________________________________________________________
!
      ! sol_nk formerly removed Lam_gbpil_T points if piloting from the top was
      ! enabled.  However, all other aspects of the solver are defined on the
      ! full vertical grid, so we should define sol_nk as a split of G_nk points.
      sol_nk = trp_12sn

      if ( Sol_type_S(1:9) == 'ITERATIVE' ) then

         if (Sol_type_S(11:12) == '2D') then
            if (Lun_out > 0) write (Lun_out,1001) trim(sol2D_precond_S)
         else
            ni=ldnh_maxx-ldnh_minx+1
            nj=ldnh_maxy-ldnh_miny+1
            call gmm_build_meta4D (savemeta,  1,ni,0,0,ni,&
                                              1,nj,0,0,nj,&
                                              1,l_nk,0,0,l_nk,&
                                    0,0,0,0,0,0,GMM_NULL_FLAGS)
            istat= gmm_create('SOL_SAVED',Sol_saved,savemeta, GMM_FLAG_RSTR+GMM_FLAG_IZER)

            if (Lun_out > 0) write (Lun_out,1002) trim(Sol3D_krylov_S), trim(sol3D_precond_S)

         end if

         do k= 1, G_nk
            wk(k) = (Cstv_hco1_8+Cstv_hco0_8*Opr_zeval_8(k))
         end do

         sol_pil_w= pil_w ; sol_pil_e= pil_e
         sol_pil_s= pil_s ; sol_pil_n= pil_n

!! Bloc extension == restrictive Jacobi preconditionner
         if (allocated(Prec_xevec_8)) deallocate (Prec_xevec_8)
         if (allocated(Prec_xeval_8)) deallocate (Prec_xeval_8)
         if (allocated(Prec_ai_8)) deallocate (Prec_ai_8)
         if (allocated(Prec_bi_8)) deallocate (Prec_bi_8)
         if (allocated(Prec_invbi_8)) deallocate (Prec_invbi_8)
         if (allocated(Prec_ci_8)) deallocate (Prec_ci_8)

         select case(sol3D_precond_S)
           case ('RAS') ! Restrictive Additive Schwarz preconditioner
              Sol_ii0  = 1    - ovlpx
              Sol_iin  = l_ni + ovlpx
              Sol_jj0  = 1    - ovlpy
              Sol_jjn  = l_nj + ovlpy

              Sol_imin = Sol_ii0
              Sol_imax = Sol_iin
              Sol_jmin = Sol_jj0
              Sol_jmax = Sol_jjn

              if (Ptopo_mycol==1)  Sol_ii0  = 0
              if (Ptopo_mycol == Ptopo_npex-2)  Sol_iin = l_ni+1
              if (Ptopo_myrow==1)  Sol_jj0  = 0
              if (Ptopo_myrow == Ptopo_npey-2)  Sol_jjn = l_nj+1

              if (l_west)  Sol_ii0 = 1 + sol_pil_w
              if (l_east)  Sol_iin = l_ni - sol_pil_e
              if (l_south) Sol_jj0 = 1 + sol_pil_s
              if (l_north) Sol_jjn = l_nj - sol_pil_n

              sol_niloc=Sol_iin-Sol_ii0+1
              sol_njloc=Sol_jjn-Sol_jj0+1
              sol_nloc = sol_niloc*sol_njloc*Schm_nith

              allocate (Prec_xevec_8(sol_niloc,sol_niloc)    ,&
                        Prec_xeval_8(sol_niloc), &
                         Prec_ai_8(sol_niloc,sol_njloc,G_nk),&
                         Prec_bi_8(sol_niloc,sol_njloc,G_nk),&
                         Prec_invbi_8(sol_niloc,sol_njloc,G_nk),&
                         Prec_ci_8(sol_niloc,sol_njloc,G_nk))

              call eigenabc_local2 (Prec_xeval_8,Prec_xevec_8,Prec_ai_8,&
                             Prec_bi_8,Prec_invbi_8,Prec_ci_8,&
                             sol_niloc,sol_njloc,Schm_nith,wk)
           case default
              sol_niloc= (l_ni-pil_e)-(1+pil_w)+1
              sol_njloc= (l_nj-pil_n)-(1+pil_s)+1
              sol_nloc = sol_niloc*sol_njloc*Schm_nith
              allocate (Prec_xevec_8(sol_niloc,sol_niloc)    ,&
                        Prec_xeval_8(sol_niloc),&
                        Prec_ai_8(sol_niloc,sol_njloc,G_nk),&
                        Prec_bi_8(sol_niloc,sol_njloc,G_nk),&
                        Prec_ci_8(sol_niloc,sol_njloc,G_nk))
              call eigenabc_local (Prec_xeval_8,Prec_xevec_8,Prec_ai_8,&
                             Prec_bi_8,Prec_ci_8,l_ni,l_nj      ,&
                              sol_niloc,sol_njloc,Schm_nith,l_i0,l_j0,wk)
         end select

         allocate (gg(1:sol_im+1),rot_cos(1:sol_im+1), rot_sin(1:sol_im+1))
         allocate (v_lcl_sum(1:sol_im+1,1:2),rr(1:sol_im+1,1:sol_im+1),&
                   tt(1:sol_im+1,1:sol_im+1), hessenberg(1:sol_im+1, 1:sol_im))
         allocate (work_space(sol_imin:sol_imax,sol_jmin:sol_jmax,1:l_nk)   ,&
                   vv(sol_imin:sol_imax,sol_jmin:sol_jmax,1:l_nk,1:sol_im+1),&
                   wint_8(sol_ii0:sol_iin,sol_jj0:sol_jjn,1:l_nk,1:sol_im+1))
         allocate (thread_s(1:4,0:OMP_get_max_threads()-1),&
                   thread_s2(1:2,1:sol_im+1,0:OMP_get_max_threads()-1))
         allocate (A1(l_minx:l_maxx, l_miny:l_maxy,1:l_nk,1:15),&
                   A2(l_minx:l_maxx, l_miny:l_maxy,1:l_nk,1:15),&
                   B1(l_minx:l_maxx, l_miny:l_maxy,1:l_nk,1:15),&
                   B2(l_minx:l_maxx, l_miny:l_maxy,1:l_nk,1:15),&
                   C1(l_minx:l_maxx, l_miny:l_maxy,1:l_nk,1:15),&
                   C2(l_minx:l_maxx, l_miny:l_maxy,1:l_nk,1:15))
         ni= Sol_iin-Sol_ii0+1
         nj= Sol_jjn-Sol_jj0+1
         allocate (fdg(ni,nj,l_nk),w2_8(ni,nj,l_nk),w3_8(ni,nj,l_nk))
         allocate (fdg2(l_minx:l_maxx,l_miny:l_maxy,l_nk+1)) ; fdg2=0.

      else                      ! Using the direct solver

         do k= 1, G_nk
            wk(k)= Cstv_hco0_8*Opr_zeval_8(k)
         end do
         ! The one-transpose solver calculates the tridiagonal coefficients in
         ! sol_prepabc rather than sol_abc, so select the right function
         if ((sol_one_transpose_L)) then
            ! If sol_sock_nk is 0, which will happen if ptopo_x > nk, then
            ! some processors will be idle during the tridiagonal solve, and
            ! they have no need of coefficients for the solve
            if (Sol_sock_nk>0)  then
               call sol_prepabc (Sol_a,Sol_b,Sol_c,wk,G_nk)
            endif
         else ! Using the two-transpose solver
               yg_8(1:G_nj)=G_yg_8(1:G_nj)
               call sol_abc ( wk,yg_8,Opr_opsyp0_8, &
                              Opr_opsyp2_8,Opr_xeval_8 , &
                   trp_12smin, trp_12smax,  sol_nk, trp_12sn0, &
                   trp_22min , trp_22max , trp_22n, trp_22n0 , &
                   G_ni,G_nj,G_nk, Sol_ai_8, Sol_bi_8, Sol_ci_8 )
                if (allocated(Sol_dg2)) deallocate(Sol_dg2)
                allocate ( Sol_dg2(1:trp_12smax,1:trp_22max,G_nj+Ptopo_npey))
          endif

         ! The correct FFT solver type depends on the fundamental grid.  The LAM
         ! solver has implicit neumann boundary conditions, whereas the yin-yang
         ! solver has dirichlet boundary conditions.  The former naturally results
         ! in a cosine expansion, the latter results in a sine expansion
         Sol_type_fft = 'QCOS'
         if (Grd_yinyang_L) Sol_type_fft = 'SIN'

         if (allocated(Sol_rhs_8)) deallocate(Sol_rhs_8)
         if (allocated(Sol_sol_8)) deallocate(Sol_sol_8)
         allocate( Sol_rhs_8(ldnh_maxx,ldnh_maxy,l_nk) ,&
                   Sol_sol_8(ldnh_maxx,ldnh_maxy,l_nk) )

         if (allocated(Sol_dwfft)) deallocate(Sol_dwfft)
         if (associated(Sol_xpose)) then
            deallocate(Sol_xpose)
            nullify(Sol_xpose)
         endif
         allocate(Sol_xpose(1:ldnh_maxx, 1:ldnh_maxy, 1:sol_nk, 1:ptopo_npex))

         if (sol_one_transpose_L) then
            ! If we use the one-transpose solver, then the FFT step also incorporates an in-memory
            ! transpose, such that the output of the forward transform is directly written to the
            ! shared-memory array.  This in-memory transpose is specified with an additional,
            ! optional parameter to make_r2r_dft_plan.

            ! With the "stacked" transposer, the transposed data is in the (i,j,k) memory order,
            ! but after the transform we want this data to be in the shared-memory array which
            ! has (k,i,j) memory order.
            allocate (Sol_dwfft(G_ni+2+Ptopo_npex,1:ldnh_maxy,1:trp_12smax))
            call make_r2r_dft_plan(forward_plan,                  & ! Plan variable
               Sol_dwfft((1+Lam_pil_w):(G_ni-Lam_pil_e),          & ! Source array, 1st dim = i
                         (1+pil_s):(ldnh_nj-pil_n),               & ! 2nd dim = j
                         1:sol_nk),                               & ! 3rd dim = k
               Sol_fft(trp_12sn0:(trp_12sn0 + sol_nk - 1),        & ! Dest array, 1st dim = k
                       (1+Lam_pil_w):(G_ni-Lam_pil_e),            & ! 2nd dim = i
                       (ldnh_j0+pil_s):(ldnh_nj+ldnh_j0-pil_n-1)),& ! 3rd dim = j
               1,                                                 & ! Transform dimension in input array
               Sol_type_fft, DFT_FORWARD,                         & ! Transform type
               [2,3,1])                                             ! Permutation array

            call make_r2r_dft_plan(reverse_plan,                  & ! Plan variable
               Sol_fft(trp_12sn0:(trp_12sn0 + sol_nk - 1),        & ! Source array, 1st dim = k
                       (1+Lam_pil_w):(G_ni-Lam_pil_e),            & ! 2nd dim = i
                       (ldnh_j0+pil_s):(ldnh_nj+ldnh_j0-pil_n-1)),& ! 3rd dim = j
               Sol_dwfft((1+Lam_pil_w):(G_ni-Lam_pil_e),          & ! Dest array, 1st dim = i
                         (1+pil_s):(ldnh_nj-pil_n),               & ! 2nd dim = j
                         1:sol_nk),                               & ! 3rd dim = k
               2, Sol_type_fft, DFT_BACKWARD, [3,1,2]) ! Note different permutation and transform dimension

            Sol_xpose = 0.0d0 ! Sol_xpose is only allocated with the one-transpose solver

         else
            ! Otherwise, the transform is in-place and does not involve a permutation.
            allocate( Sol_dwfft(1:ldnh_maxy , 1:trp_12smax, G_ni+2+Ptopo_npex))
            call make_r2r_dft_plan(forward_plan, & ! Plan variable
                Sol_dwfft((1+pil_s):(ldnh_maxy-pil_n),1:sol_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                Sol_dwfft((1+pil_s):(ldnh_maxy-pil_n),1:sol_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                                   3, Sol_type_fft, DFT_FORWARD)
            call make_r2r_dft_plan(reverse_plan, & ! Plan variable
                Sol_dwfft((1+pil_s):(ldnh_maxy-pil_n),1:sol_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                Sol_dwfft((1+pil_s):(ldnh_maxy-pil_n),1:sol_nk,(1+Lam_pil_w):(G_ni-Lam_pil_e)), &
                                   3, Sol_type_fft, DFT_BACKWARD)
         endif
         Sol_dwfft= 0.d0 ; Sol_rhs_8= 0.d0 ; Sol_sol_8= 0.d0
         Sol_pri= get_dft_norm_factor(G_ni-Lam_pil_w-Lam_pil_e,Sol_type_fft)
         Sol_pri= Sol_pri/(G_xg_8(G_ni-Lam_pil_e+1)-G_xg_8(G_ni-Lam_pil_e))
         Sol_pil_w=0 ; Sol_pil_e=0
         if (l_south) Sol_pil_w= Lam_pil_w
         if (l_north) Sol_pil_e= Lam_pil_e

      end if


 1001 format(/,'WILL USE FGMRES 2D ITERATIVE SOLVER WITH ',a, &
               ' PRECONDITIONNER' &
             /,'=============================================')
 1002 format(/,'WILL USE ',a,' 3D ITERATIVE SOLVER WITH ',a, &
               ' PRECONDITIONNER' &
             /,'=============================================')
!
!     ---------------------------------------------------------------
!
      return
      end
