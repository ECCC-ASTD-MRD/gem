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
      use glb_ld
      use glb_pil
      use gmm_vt1
      use gmm_vt0
      use gmm_nest
      use mem_tstp
      use gmm_geof
      use gmm_pw
      use gmm_vth
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
      use gem_timing
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(IN) :: F_icn

      logical :: print_conv
      integer i0, in, j0, jn, k0, ni, nj, iln, icln
!
!     ---------------------------------------------------------------
!
      call rpn_comm_xch_halo( ut1 , l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj ,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( vt1 , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_njv,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( tt1 , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( st1 , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,1   , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      if (.not.Dynamics_hydro_L) then
         call rpn_comm_xch_halo(  qt1, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1,&
                                  G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      end if

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

         call gemtime_start ( 20, 'RHS', 10 )

!        Compute the right-hand sides of the governing equations
         call rhs ( orhsu, orhsv, orhsc, orhst, orhsw, orhsf,&
                    ut1, vt1, wt1               ,&
                    tt1, st1, zdt1, qt1, sls, fis0     ,&
                    l_minx,l_maxx,l_miny,l_maxy, l_nk )

         call gemtime_stop (20)

         call firstguess ()

! Perform time interpolation of Lateral BCs for LAM configurations

         if ( .not. Grd_yinyang_L ) call nest_bcs (rhsu, rhsv, &
                              l_minx,l_maxx,l_miny,l_maxy, l_nk)

      end if

!     Perform Semi-Lagrangian advection

      call gemtime_start (21, 'ADZ_MAIN', 10)

      call adz_main ( orhsu, rhsu, orhsv, rhsv, orhsc, rhsc, orhst,&
                      rhst, orhsf, rhsf, orhsw, rhsw              ,&
                      l_minx,l_maxx,l_miny,l_maxy, l_nk )

      call gemtime_stop(21)

      call gemtime_start (22, 'PRE', 10)

      if ( F_icn == 1 ) then
         if ( .not. Grd_yinyang_L .and. .not. Lam_ctebcs_L) then
            fis0(1:l_ni,1:l_nj)= nest_fullme(1:l_ni,1:l_nj)
            call rpn_comm_xch_halo (fis0,l_minx,l_maxx,l_miny,l_maxy,&
               l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         else
            if (Vtopo_L .and. (Lctl_step >= Vtopo_start)) then
               call var_topo (fis0, real(Lctl_step),&
                              l_minx,l_maxx,l_miny,l_maxy)
               if (Grd_yinyang_L) then
                  call yyg_xchng (fis0, l_minx,l_maxx,l_miny,l_maxy, &
                                  l_ni,l_nj, 1, .false., 'CUBIC', .true.)
               else
               call rpn_comm_xch_halo (fis0,l_minx,l_maxx,l_miny,l_maxy,&
                  l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
               end if
            end if
         end if
      end if

!     Combine some rhs to obtain the linear part
!     of the right-hand side of the elliptic problem
      call pre (rhsu, rhsv, fis0, rhsc, rhst, &
                rhsw, rhsf, rhsb, nest_t, l_minx,l_maxx,l_miny,l_maxy,&
                i0, j0, in, jn, k0, l_nk)

      call gemtime_stop (22)

      if ( Lun_debug_L ) write (Lun_out,1005) Schm_itnlh

      ni = ldnh_maxx-ldnh_minx+1
      nj = ldnh_maxy-ldnh_miny+1

      do iln=1,Schm_itnlh

         call gemtime_start ( 23, 'NLI', 10 )

!        Compute non-linear components and combine them
!        to obtain final right-hand side of the elliptic problem
         icln=F_icn*iln
         if ( .not. Grd_yinyang_L ) icln=icln+1
         call nli (nl_u, nl_v, nl_t, nl_c, nl_w, nl_f          ,&
                   ut0, vt0, tt0, st0, zdt0, qt0, rhs_sol, rhsc,&
                   sls, fis0, nl_b, l_minx,l_maxx,l_miny,l_maxy,&
                   l_nk, ni, nj, i0, j0, in, jn, k0, icln)

         call gemtime_stop (23)

         call gemtime_start ( 24, 'SOL', 10 )

!        Solve the elliptic problem
         print_conv = (iln   == Schm_itnlh ) .and. &
                      (F_icn == Schm_itcn  ) .and. &
                      (Ptopo_couleur == 0  ) .and. &
                      (Lun_out > 0)

         call sol_main (rhs_sol,lhs_sol,ni,nj, l_nk, print_conv)

         call gemtime_stop (24)

         call gemtime_start ( 25, 'BAC', 10 )

!        Back subtitution
         call  bac (lhs_sol, sls, fis0                        ,&
                    ut0, vt0, wt0, tt0, st0, zdt0, qt0, nest_q,&
                    rhsu, rhsv, rhst, rhsw, rhsf, rhsb        ,&
                    nl_u, nl_v, nl_t, nl_w, nl_f, nl_b        ,&
                    l_minx, l_maxx, l_miny, l_maxy            ,&
                    ni,nj,l_nk,i0, j0, k0, in, jn)

         call gemtime_stop (25)

         if (Grd_yinyang_L) then
            call yyg_xchng_vec_uv2uv (ut0, vt0,&
                                      l_minx,l_maxx,l_miny,l_maxy,G_nk)

            call yyg_xchng (tt0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng (zdt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng (st0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            1,    .false., 'CUBIC', .false.)
            if (.not.Dynamics_hydro_L) then
               call yyg_xchng (qt0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                               G_nk+1, .false., 'CUBIC', .false.)
            end if
        end if

      end do

      if (Grd_yinyang_L) then
         call yyg_xchng (wt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .false., 'CUBIC', .false.)
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
