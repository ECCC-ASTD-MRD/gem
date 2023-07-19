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

      subroutine adz_tracers_interp ()
      use ens_options
      use adz_interp_hlt_mod
      use mem_tracers
      use omp_timing
      implicit none

      integer ::  n, deb
!
!     ---------------------------------------------------------------
!
      call gtmg_start (33, 'ADV_tracers', 10)

      if (Adz_verbose>0) call stat_mass_tracers_hlt (1,"BEFORE ADVECTION")

      if (Tr3d_ntrTRICUB_NT>0) then
         call gtmg_start (34, 'TRICUB_NT', 33)
         deb= Tr3d_debTRICUB_NT
!$omp single
         do n=1, Tr3d_ntrTRICUB_NT ! Tricubic NO post treatment
            pnt_stack(n)%src => tracers_P(deb+n-1)%pntr
            pnt_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
!$omp end single
         call adz_tricub_hlt ( pnt_stack, Tr3d_ntrTRICUB_NT ,&
                  Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                  Adz_j0,Adz_jn,Adz_k0t,F_ext_L=.false. )
         call gtmg_stop (34)
      end if

      if (Tr3d_ntrBICHQV_NT>0) then
         call gtmg_start (35, 'BICHQV_NT', 33)
         deb= Tr3d_debBICHQV_NT
!$omp single
         do n=1, Tr3d_ntrBICHQV_NT ! BicubicH+QuinticV NO post treatment
            pnt_stack(n)%src => tracers_P(deb+n-1)%pntr
            pnt_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
!$omp end single
         call adz_tricub_hlt ( pnt_stack, Tr3d_ntrBICHQV_NT ,&
                  Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                  Adz_j0,Adz_jn,Adz_k0t,F_ext_L=.false.,F_QV_L=.true. )
         call gtmg_stop (35)
      end if

      if (Tr3d_ntrTRICUB_WP>0) then
         call gtmg_start (37, 'TRICUB_WP', 33)
         deb= Tr3d_debTRICUB_WP
!$omp single
         Adz_post => Adz_post_3CWP
         Adz_flux => Adz_flux_3CWP
         do n=1, Tr3d_ntrTRICUB_WP
            pnt_stack(n)%src => tracers_P(deb+n-1)%pntr
            pnt_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
!$omp end single
         if (.not.Adz_BC_LAM_zlf_L) then
            call adz_tricub_hlt ( pnt_stack, Tr3d_ntrTRICUB_WP ,&
                  Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                  Adz_j0,Adz_jn,Adz_k0t,F_ext_L=.false.,F_post=Tr_3CWP)
         else
!$omp single
            do n=1, Tr3d_ntrTRICUB_WP
               pnt_stack(n)%pil => tracers_B(deb+n-1)%pntr
            end do
!$omp end single
            call adz_BC_LAM_zlf_0_hlt (pnt_stack,Tr3d_ntrTRICUB_WP,0)
            call adz_tricub_hlt ( pnt_stack, Tr3d_ntrTRICUB_WP ,&
                  Adz_pb,Adz_cpntr_t,Adz_num_b,Adz_i0b,Adz_inb,&
                  Adz_j0b,Adz_jnb,1,F_ext_L=.false.,F_post=Tr_3CWP)
            call adz_BC_LAM_zlf_0_hlt (pnt_stack,Tr3d_ntrTRICUB_WP,1)
         end if
         call gtmg_stop (37)
      end if

      if (Tr3d_ntrBICHQV_WP>0) then
         call gtmg_start (39, 'ADZ_BICHQV_WP', 33)
         deb= Tr3d_debBICHQV_WP
!$omp single
         Adz_post => Adz_post_BQWP
         Adz_flux => Adz_flux_BQWP
         do n=1, Tr3d_ntrBICHQV_WP
            pnt_stack(n)%src => tracers_P(deb+n-1)%pntr
            pnt_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
!$omp end single
         if (.not.Adz_BC_LAM_zlf_L) then
            call adz_tricub_hlt ( pnt_stack, Tr3d_ntrBICHQV_WP,&
                    Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                    Adz_j0,Adz_jn,Adz_k0t, F_ext_L=.false.    ,&
                    F_QV_L=.true.,F_post=Tr_BQWP  )
         else
!$omp single
            do n=1, Tr3d_ntrBICHQV_WP
               pnt_stack(n)%pil => tracers_B(deb+n-1)%pntr
            end do
!$omp end single
            call adz_BC_LAM_zlf_0_hlt (pnt_stack,Tr3d_ntrBICHQV_WP,0)
            call adz_tricub_hlt ( pnt_stack, Tr3d_ntrBICHQV_WP,&
                  Adz_pb,Adz_cpntr_t,Adz_num_b,Adz_i0b,Adz_inb,&
                  Adz_j0b,Adz_jnb,1, F_ext_L=.false.          ,&
                  F_QV_L=.true.,F_post=Tr_BQWP  )
            call adz_BC_LAM_zlf_0_hlt (pnt_stack,Tr3d_ntrBICHQV_WP,1)
         end if
         call gtmg_stop (39)
      end if

      call gtmg_stop (33)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_tracers_interp
