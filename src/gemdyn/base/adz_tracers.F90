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

      subroutine adz_tracers ()
      use adv_pos
      use adz_mem
      use adz_options
      use adz_interp_rhs_mod
      use gem_timing
      use gmm_itf_mod
      use tr3d
      implicit none
#include <arch_specific.hf>

      integer ::  n, err
!
!     ---------------------------------------------------------------
!
      call gemtime_start (32, 'ADZ_TRNT_3C', 21)
      do n=1, Tr3d_ntrTRICUB_NT ! Tricubic NO post treatment
         err= gmm_get ( 'TR/'//trim(Tr3d_TRICUB_NT_S(n))//':P' ,&
                        Adz_stack(n)%src )
         err= gmm_get ( 'TR/'//trim(Tr3d_TRICUB_NT_S(n))//':M' ,&
                        Adz_stack(n)%dst )
      end do
      call adz_tricub_rhs ( Adz_stack, Tr3d_ntrTRICUB_NT       ,&
           Adz_pt(1,Adz_i0,Adz_j0,Adz_k0),Adz_cpntr_t,Adz_num_q,&
           Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )
      call gemtime_stop (32)

      if (Tr3d_ntrBICHQV_NT>0) then
         call gemtime_start (33, 'ADZ_TRNT_2CQV', 21)
         do n=1, Tr3d_ntrBICHQV_NT ! BicubicH+QuinticV NO post treatment
            err= gmm_get( 'TR/'//trim(Tr3d_BICHQV_NT_S(n))//':P' ,&
                          Adz_stack(n)%src )
            err= gmm_get( 'TR/'//trim(Tr3d_BICHQV_NT_S(n))//':M' ,&
                           Adz_stack(n)%dst )
         end do
         call adz_bicubHQV_rhs ( Adz_stack, Tr3d_ntrBICHQV_NT,&
              pxt,pyt,pzt,Adz_num_q                          ,&
              Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )
         call gemtime_stop (33)
      end if

      if (max(Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)>0) call set_post_tr ()

      if (Tr3d_ntrTRICUB_WP>0) then
         call gemtime_start (34, 'ADZ_TRWP_3C', 21)
         do n=1, Tr3d_ntrTRICUB_WP
            err= gmm_get( 'TR/'//trim(Tr3d_TRICUB_WP_S(n))//':P' ,&
                          Adz_stack(n)%src )
            err= gmm_get( 'TR/'//trim(Tr3d_TRICUB_WP_S(n))//':M' ,&
                           Adz_stack(n)%dst )
         end do
         call adz_BC_LAM_zlf_0 (Tr3d_TRICUB_WP_S,Tr3d_TRICUB_WP,Tr3d_ntrTRICUB_WP,1)
         if (.not.Adz_BC_LAM_zlf_L) then
            call adz_tricub_rhs ( Adz_stack,Tr3d_ntrTRICUB_WP         ,&
                  Adz_pt(1,Adz_i0,Adz_j0,Adz_k0),Adz_cpntr_t,Adz_num_q,&
                  Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0,F_post=Tr3d_TRICUB_WP)
         else
            call adz_tricub_rhs ( Adz_stack,Tr3d_ntrTRICUB_WP           ,&
                  Adz_pb(1,Adz_i0b,Adz_j0b,Adz_k0),Adz_cpntr_t,Adz_num_b,&
                  Adz_i0b,Adz_inb,Adz_j0b,Adz_jnb,Adz_k0,F_post=Tr3d_TRICUB_WP)
         end if
         call gemtime_stop (34)

         call gemtime_start (36, 'ADZ_TR_POST', 21)
         call adz_post_tr (Tr3d_TRICUB_WP_S,Tr3d_TRICUB_WP,Tr3d_ntrTRICUB_WP)
         call adz_BC_LAM_zlf_0 (Tr3d_TRICUB_WP_S,Tr3d_TRICUB_WP,Tr3d_ntrTRICUB_WP,2)
         call gemtime_stop (36)
      end if

      if (Tr3d_ntrBICHQV_WP>0) then
         call gemtime_start (35, 'ADZ_TRWP_2CQV', 21)
         do n=1, Tr3d_ntrBICHQV_WP
            err= gmm_get( 'TR/'//trim(Tr3d_BICHQV_WP_S(n))//':P' ,&
                          Adz_stack(n)%src )
            err= gmm_get( 'TR/'//trim(Tr3d_BICHQV_WP_S(n))//':M' ,&
                          Adz_stack(n)%dst )
         end do
         call adz_BC_LAM_zlf_0 (Tr3d_BICHQV_WP_S,Tr3d_BICHQV_WP,Tr3d_ntrBICHQV_WP,1)
         if (.not.Adz_BC_LAM_zlf_L) then
            call adz_bicubHQV_rhs ( Adz_stack,Tr3d_ntrBICHQV_WP,&
                  pxt,pyt,pzt,Adz_num_q                        ,&
                  Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0,F_post=Tr3d_BICHQV_WP)
         else
            call adz_bicubHQV_rhs ( Adz_stack,Tr3d_ntrBICHQV_WP,&
                  pxt,pyt,pzt,Adz_num_b                        ,&
                  Adz_i0b,Adz_inb,Adz_j0b,Adz_jnb,Adz_k0,F_post=Tr3d_BICHQV_WP)
         end if
         call gemtime_stop (35)

         call gemtime_start (36, 'ADZ_TR_POST', 21)
         call adz_post_tr (Tr3d_BICHQV_WP_S,Tr3d_BICHQV_WP,Tr3d_ntrBICHQV_WP)
         call adz_BC_LAM_zlf_0 (Tr3d_BICHQV_WP_S,Tr3d_BICHQV_WP,Tr3d_ntrBICHQV_WP,2)
         call gemtime_stop (36)
      end if

      if ( (max(Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)>0) .and. &
           (Adz_verbose>0) ) call stat_mass_tracers (0,"AFTER ADVECTION")
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_tracers
