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

      subroutine adz_tracers (F_before_psadj_L)
      use ens_options
      use adz_interp_rhs_mod
      use adz_interp_hlt_mod
      use mem_tracers
      use omp_timing
      implicit none

      logical, intent(in) :: F_before_psadj_L

      integer ::  n, deb
      logical ::  after_psadj_L!, sto_phy_L
!
!     ---------------------------------------------------------------
!
      after_psadj_L = .not.F_before_psadj_L
!     sto_phy_L= associated(mcrhsint)
      
      call gtmg_start (33, 'ADZ_TRACERS', 10)

!$omp single
      if (Adz_verbose>0 .and. F_before_psadj_L) call stat_mass_tracers (1,"BEFORE ADVECTION")
!$omp end single

      if (Tr3d_ntrTRICUB_NT>0 .and. F_before_psadj_L) then
         call gtmg_start (34, 'TRICUB_NT', 33)
         deb= Tr3d_debTRICUB_NT
         do n=1, Tr3d_ntrTRICUB_NT ! Tricubic NO post treatment
            Adz_stack(n)%src => tracers_P(deb+n-1)%pntr
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
         call adz_tricub_hlt ( Adz_stack, Tr3d_ntrTRICUB_NT ,&
                  Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                  Adz_j0,Adz_jn,Adz_k0t,F_ext_L=.false. )
         call gtmg_stop (34)
      end if

      if (Tr3d_ntrBICHQV_NT>0 .and. F_before_psadj_L) then
         call gtmg_start (35, 'BICHQV_NT', 33)
         deb= Tr3d_debBICHQV_NT
         do n=1, Tr3d_ntrBICHQV_NT ! BicubicH+QuinticV NO post treatment
            Adz_stack(n)%src => tracers_P(deb+n-1)%pntr
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
         call adz_tricub_hlt ( Adz_stack, Tr3d_ntrBICHQV_NT ,&
                  Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                  Adz_j0,Adz_jn,Adz_k0t,F_ext_L=.false.,F_QV_L=.true. )
         call gtmg_stop (35)
      end if

      !Resetting done at each timestep before calling adz_post_tr
      !----------------------------------------------------------
      if (max(Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)>0 .and. after_psadj_L) call set_post_tr_hlt ()

      if (Tr3d_ntrTRICUB_WP>0) then
         call gtmg_start (37, 'TRICUB_WP', 33)
         Adz_post => Adz_post_3CWP
         Adz_flux => Adz_flux_3CWP
         deb= Tr3d_debTRICUB_WP
         do n=1, Tr3d_ntrTRICUB_WP
            Adz_stack(n)%src => tracers_P(deb+n-1)%pntr
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
         if (.not.Adz_BC_LAM_zlf_L .and. F_before_psadj_L) then
            call adz_tricub_hlt ( Adz_stack, Tr3d_ntrTRICUB_WP ,&
                  Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                  Adz_j0,Adz_jn,Adz_k0t,F_ext_L=.false.,F_post=Tr_3CWP)
         else if (F_before_psadj_L) then
            do n=1, Tr3d_ntrTRICUB_WP
               Adz_stack(n)%pil => tracers_B(deb+n-1)%pntr
            end do
            call adz_BC_LAM_zlf_0_hlt (Tr3d_ntrTRICUB_WP,0)
            call adz_tricub_hlt ( Adz_stack, Tr3d_ntrTRICUB_WP ,&
                  Adz_pb,Adz_cpntr_t,Adz_num_b,Adz_i0b,Adz_inb,&
                  Adz_j0b,Adz_jnb,1,F_ext_L=.false.,F_post=Tr_3CWP)
            call adz_BC_LAM_zlf_0_hlt (Tr3d_ntrTRICUB_WP,1)
         end if
         call gtmg_stop (37)

         !Apply ILMC shape-preserving for Tr_3CWP
         !---------------------------------------
!$omp single
         if (after_psadj_L) call adz_post_tr (1)
!$omp end single

      end if

      if (Tr3d_ntrBICHQV_WP>0) then
         call gtmg_start (39, 'ADZ_BICHQV_WP', 33)
         Adz_post => Adz_post_BQWP
         Adz_flux => Adz_flux_BQWP
         deb= Tr3d_debBICHQV_WP
         do n=1, Tr3d_ntrBICHQV_WP
            Adz_stack(n)%src => tracers_P(deb+n-1)%pntr
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
         if (.not.Adz_BC_LAM_zlf_L .and. F_before_psadj_L) then
            call adz_tricub_hlt ( Adz_stack, Tr3d_ntrBICHQV_WP,&
                    Adz_pt,Adz_cpntr_t,Adz_num_t,Adz_i0,Adz_in,&
                    Adz_j0,Adz_jn,Adz_k0t, F_ext_L=.false.    ,&
                    F_QV_L=.true.,F_post=Tr_BQWP  )
         else if (F_before_psadj_L) then
            do n=1, Tr3d_ntrBICHQV_WP
               Adz_stack(n)%pil => tracers_B(deb+n-1)%pntr
            end do
            call adz_BC_LAM_zlf_0_hlt (Tr3d_ntrBICHQV_WP,0)
            call adz_tricub_hlt ( Adz_stack, Tr3d_ntrBICHQV_WP,&
                  Adz_pb,Adz_cpntr_t,Adz_num_b,Adz_i0b,Adz_inb,&
                  Adz_j0b,Adz_jnb,1, F_ext_L=.false.          ,&
                  F_QV_L=.true.,F_post=Tr_BQWP  )
            call adz_BC_LAM_zlf_0_hlt (Tr3d_ntrBICHQV_WP,1)
         end if
         call gtmg_stop (39)

         !Apply ILMC shape-preserving for Tr_BQWP
         !---------------------------------------
!$omp single
         if (after_psadj_L) call adz_post_tr (2)
!$omp end single

      end if

      !Apply Bermejo-Conde mass-fixer for all tracers in Adz_bc
      !--------------------------------------------------------
!$omp single
      if (max(Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)>0 .and. after_psadj_L) call adz_post_tr (0)
!$omp end single

      if (Tr3d_ntrTRICUB_WP>0 .and. Adz_BC_LAM_zlf_L .and. after_psadj_L) then
         deb= Tr3d_debTRICUB_WP
         do n=1, Tr3d_ntrTRICUB_WP
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
            Adz_stack(n)%pil => tracers_B(deb+n-1)%pntr
         end do
         call adz_BC_LAM_zlf_0_hlt (Tr3d_ntrTRICUB_WP,2)
      end if

      if (Tr3d_ntrBICHQV_WP>0 .and. Adz_BC_LAM_zlf_L .and. after_psadj_L) then
         deb= Tr3d_debBICHQV_WP
         do n=1, Tr3d_ntrBICHQV_WP
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
            Adz_stack(n)%pil => tracers_B(deb+n-1)%pntr
         end do
         call adz_BC_LAM_zlf_0_hlt (Tr3d_ntrBICHQV_WP,2)
      end if

!$omp single
      if (Adz_verbose>0 .and. after_psadj_L) call stat_mass_tracers (0,"AFTER ADVECTION")
!$omp end single

      call gtmg_stop (33)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_tracers
