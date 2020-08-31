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

      subroutine adz_tracers_HY ()
      use adv_pos
      use adz_interp_rhs_mod
      use mem_tracers
      implicit none
#include <arch_specific.hf>

      !object
      !===========================================================================
      !     Dryair Conservation: Advect Humidity/Hydrometeors Tracers before PSADJ
      !===========================================================================

      integer ::  n, deb
!
!     ---------------------------------------------------------------
!
      if (Tr3d_ntrTRICUB_HY>0) then
         deb= Tr3d_debTRICUB_HY
         do n=1, Tr3d_ntrTRICUB_HY
            Adz_stack(n)%src => tracers_P(deb+n-1)%pntr
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
         call adz_tricub_rhs ( Adz_stack, Tr3d_ntrTRICUB_HY    ,&
           Adz_pt(1,Adz_i0,Adz_j0,Adz_k0),Adz_cpntr_t,Adz_num_q,&
           Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )
      end if

      if (Tr3d_ntrBICHQV_HY>0) then
         deb= Tr3d_debBICHQV_HY
         do n=1, Tr3d_ntrBICHQV_HY
            Adz_stack(n)%src => tracers_P(deb+n-1)%pntr
            Adz_stack(n)%dst => tracers_M(deb+n-1)%pntr
         end do
         call adz_bicubHQV_rhs ( Adz_stack, Tr3d_ntrBICHQV_HY,&
                                 pxt,pyt,pzt,Adz_num_q       ,&
                           Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )
      end if
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_tracers_HY
