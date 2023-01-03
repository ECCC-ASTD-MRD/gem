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

      subroutine tstpdyn (F_icn)
      use cstv
      use glb_ld
      use glb_pil
      use gmm_vt1
      use gmm_vt0
      use mem_nest
      use mem_tstp
      use adz_mem
      use gmm_geof
      use gmm_pw
      use gem_options
      use step_options
      use lam_options
      use dynkernel_options
      use dyn_fisl_options
      use HORgrid_options
      use ldnh
      use lun
      use tdpack
      use yyg_param
      use omp_timing
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(IN) :: F_icn

      logical :: print_conv
      integer i0, in, j0, jn, k0, ni, nj, iln, icln
!
!     ---------------------------------------------------------------
!
! Avoid non-associated pointers
      if (Grd_yinyang_L) then
         nest_t => ut1
         nest_q => ut1
      end if

      i0= 1   +pil_w
      in= l_ni-pil_e
      j0= 1   +pil_s
      jn= l_nj-pil_n
      k0= 1+Lam_gbpil_T

      if (Lun_debug_L) write(Lun_out,1000)

      if ( F_icn == 1 ) then       ! Compute RHS

         call gtmg_start ( 20, 'RHS', 10 )

!        Compute the right-hand sides of the governing equations
         call fislp_rhs ( orhsu, orhsv, orhsc, orhst, orhsw, orhsf,&
                          ut1, vt1, wt1               ,&
                          tt1, st1, zdt1, qt1, sls, fis0     ,&
                          l_minx,l_maxx,l_miny,l_maxy, l_nk )
         Adz_cnt_traj= 0 ; Adz_cnt_int= 0
         call rpn_comm_xch_halox( orhsu, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, Adz_halox,Adz_haloy,.false.,.false.,&
              orhsu_ext,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,l_ni,0)
         call rpn_comm_xch_halox( orhsv, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, Adz_halox,Adz_haloy,.false.,.false.,&
              orhsv_ext,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,l_ni,0)
         call rpn_comm_xch_halox( orhst, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, Adz_halox,Adz_haloy,.false.,.false.,&
              orhst_ext,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,l_ni,0)
         call rpn_comm_xch_halox( orhsc, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, Adz_halox,Adz_haloy,.false.,.false.,&
              orhsc_ext,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,l_ni,0)
         call rpn_comm_xch_halox( orhsf, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, Adz_halox,Adz_haloy,.false.,.false.,&
              orhsf_ext,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,l_ni,0)
         if (.not.Dynamics_hydro_L) then
         call rpn_comm_xch_halox( orhsw, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,l_nk, Adz_halox,Adz_haloy,.false.,.false.,&
              orhsw_ext,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,l_ni,0)
         endif
         call gtmg_stop (20)

         call firstguess ()

! Perform time interpolation of Lateral BCs for LAM configurations

         if ( .not. Grd_yinyang_L ) call nest_bcs (Cstv_dt_8,rhsu, rhsv, &
                                      l_minx,l_maxx,l_miny,l_maxy, l_nk)

      end if

      call rpn_comm_xch_halo( ut0,l_minx,l_maxx,l_miny,l_maxy,&
                              l_niu,l_nj,l_nk,G_halox,G_haloy,&
                              G_periodx,G_periody,G_niu,0)
      call rpn_comm_xch_halo( vt0,l_minx,l_maxx,l_miny,l_maxy,&
                              l_ni,l_njv,l_nk,G_halox,G_haloy,&
                              G_periodx,G_periody,G_ni,0)

! Perform Semi-Lagrangian advection

      call gtmg_start (21, 'ADZ_MAIN', 10)
      Adz_icn= F_icn
!$omp parallel
      call adz_main ()
!$omp end parallel

      call gtmg_stop(21)

      call gtmg_start (22, 'PRE', 10)

      if ( F_icn == 1 ) call oro_adj ()

!     Combine some rhs to obtain the linear part
!     of the right-hand side of the elliptic problem
      call pre ( rhsu, rhsv, fis0, rhsc, rhst, rhsw, rhsf ,&
                 rhsb, nest_t, l_minx,l_maxx,l_miny,l_maxy,&
                 i0, j0, in, jn, k0, l_nk )

      call gtmg_stop (22)

      if ( Lun_debug_L ) write (Lun_out,1005) Schm_itnlh

      ni = ldnh_maxx-ldnh_minx+1
      nj = ldnh_maxy-ldnh_miny+1

      do iln=1,Schm_itnlh

         call gtmg_start ( 23, 'NLI', 10 )

!        Compute non-linear components and combine them
!        to obtain final right-hand side of the elliptic problem
         icln=F_icn*iln ; Adz_itnl= iln
         if ( .not. Grd_yinyang_L ) icln=icln+1
         call nli (nl_u, nl_v, nl_t, nl_c, nl_w, nl_f          ,&
                   ut0, vt0, tt0, st0, zdt0, qt0, rhs_sol, rhsc,&
                   sls, fis0, nl_b, l_minx,l_maxx,l_miny,l_maxy,&
                   l_nk, ni, nj, i0, j0, in, jn, k0, icln)

         call gtmg_stop (23)

         call gtmg_start ( 24, 'SOL', 10 )

!        Solve the elliptic problem
         print_conv = (iln   == Schm_itnlh ) .and. &
                      (F_icn == Schm_itcn  ) .and. &
                      (Ptopo_couleur == 0  ) .and. &
                      (Lun_out > 0)
!$omp parallel
         call sol_main (rhs_sol,lhs_sol,ni,nj, l_nk, print_conv)
!$omp end parallel

         call gtmg_stop (24)

         call gtmg_start ( 25, 'BAC', 10 )

!        Back subtitution
         call  bac (lhs_sol, sls, fis0                        ,&
                    ut0, vt0, wt0, tt0, st0, zdt0, qt0, nest_q,&
                    rhsu, rhsv, rhst, rhsw, rhsf, rhsb        ,&
                    nl_u, nl_v, nl_t, nl_w, nl_f, nl_b        ,&
                    l_minx, l_maxx, l_miny, l_maxy            ,&
                    ni,nj,l_nk,i0, j0, k0, in, jn)

         call gtmg_stop (25)

         if (Grd_yinyang_L) then
            call yyg_xchng_vec_uv2uv (ut0, vt0,&
                                      l_minx,l_maxx,l_miny,l_maxy,G_nk)

            call yyg_int_xch_scal (tt0 , G_nk, .false., 'CUBIC', .false.)
            call yyg_int_xch_scal (zdt0, G_nk, .false., 'CUBIC', .false.)
            call yyg_int_xch_scal (st0 , 1,    .false., 'CUBIC', .false.)
            if (.not.Dynamics_hydro_L) then
               call yyg_int_xch_scal (qt0 , G_nk+1, .false., 'CUBIC', .false.)
            end if
        end if

      end do

      if (Grd_yinyang_L) then
         call yyg_int_xch_scal (wt0, G_nk, .false., 'CUBIC', .false.)
      end if

!     ---------------------------------------------------------------
!
 1000 format( &
       3X,'PERFORM A DYNAMICAL STEP: (S/R TSTPDYN)', &
      /3X,'========================================',/)
 1005 format( &
       3X,'ITERATING SCHM_ITNLH=',I3,' TIMES TO SOLVE NON-LINEAR ', &
          'HELMHOLTZ PROBLEM')

      return
      end
