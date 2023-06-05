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

!**s/r tstpdyn -  Performs a dynamical timestep of the model

      subroutine fislh_tstpdyn (F_icn)
      use adz_mem
      use cstv
      use glb_ld
      use glb_pil
      use gmm_vt0
      use mem_tstp
      use lam_options
      use dynkernel_options
      use dyn_fisl_options
      use HORgrid_options
      use ldnh
      use yyg_param
      use omp_timing
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: F_icn

      logical :: print_conv
      integer i0, in, j0, jn, k0, k0t, ni, nj, iln, icln
      real(kind=REAL64) :: dt_8
!     
!     ---------------------------------------------------------------
!
      i0= 1   +pil_w
      in= l_ni-pil_e
      j0= 1   +pil_s
      jn= l_nj-pil_n
      k0= 1+Lam_gbpil_T
      k0t=k0
      if (Schm_opentop_L) k0t=k0-1
      ni = ldnh_maxx-ldnh_minx+1
      nj = ldnh_maxy-ldnh_miny+1
      
      call gtmg_start (20, 'TSTPDYN', 10)

      dt_8= Cstv_dt_8
      
      if ( F_icn == 1 ) then  ! Compute RHS

!        Compute the right-hand sides of the governing equations

         call gtmg_start (21, 'RHS', 20)
         call fislh_rhs ( dt_8 )
         call gtmg_stop (21)
!$omp single
         call rpn_comm_xch_halo (orhsu_ext, Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,&
                    l_ni,l_nj, 6*l_nk, Adz_halox,Adz_haloy, .false.,.false., l_ni,0)
!$omp end single

! Compute time interpolation of Lateral BCs for LAM configurations
         if ( .not. Grd_yinyang_L ) call nest_bcs (dt_8, rhsu, rhsv, &
                                  l_minx,l_maxx,l_miny,l_maxy, l_nk)
      end if
                           
!$omp single
      call rpn_comm_xch_halo( ut0,l_minx,l_maxx,l_miny,l_maxy,&
                             l_ni,l_nj,2*l_nk,G_halox,G_haloy,&
                             G_periodx,G_periody,G_niu,0)
!$omp end single

!     Perform Semi-Lagrangian advection

      call gtmg_start (22, 'ADVECTION', 20)
      call adz_main ()
      call gtmg_stop (22)

      call gtmg_start (23, 'PRE', 20)
      if ( F_icn == 1 ) call oro_adj ()
      
!     Combine some rhs to obtain the linear part
!     of the right-hand side of the elliptic problem

      call fislh_pre (dt_8, i0, j0, k0, in, jn, k0t )
      call gtmg_stop (23)

      do iln=1,Schm_itnlh

!        Compute non-linear components and combine them
!        to obtain final right-hand side of the elliptic problem
         icln=F_icn*iln
         if ( .not. Grd_yinyang_L ) icln=icln+1

!$omp single
         if(icln > 1) then
            call rpn_comm_xch_halo( ut0, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj,&
                           4*G_nk+1,G_halox, G_haloy, G_periodx, G_periody, l_ni, 0 )
         end if
!$omp end single

         call gtmg_start (24, 'NLI', 20)
         call fislh_nli (dt_8, i0, j0, k0, in, jn, k0t)
         call gtmg_stop (24)

!        Solve the elliptic problem
         print_conv = (iln   == Schm_itnlh ) .and. &
                      (F_icn == Schm_itcn  ) .and. &
                      (Ptopo_couleur == 0  ) .and. &
                      (Lun_out > 0)

         call gtmg_start (25, 'SOL', 20)
         call sol_main (rhs_sol,lhs_sol,ni,nj, l_nk, print_conv)
         call gtmg_stop (25)

!        Back subtitution
         call gtmg_start (26, 'BAC', 20)
         call fislh_bac (dt_8, i0, j0, k0, in, jn, k0t)
         call gtmg_stop (26)
!$omp single
         if (Grd_yinyang_L) then
            call yyg_xchng_vec_uv2uv (ut0, vt0,&
                                      l_minx,l_maxx,l_miny,l_maxy,G_nk)

            call yyg_xchng (tt0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng (zdt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng (st0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            1,    .false., 'CUBIC', .false.)
            call yyg_xchng (qt0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk+1, .false., 'CUBIC', .false.)
        end if
!$omp end single

      end do

!$omp single
      if (Grd_yinyang_L) then
         call yyg_xchng (wt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .false., 'CUBIC', .false.)
      end if
!$omp end single
                      
      call gtmg_stop (20)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine fislh_tstpdyn
