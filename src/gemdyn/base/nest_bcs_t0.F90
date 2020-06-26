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

!**s/r nest_bcs_t0 -
!
      subroutine nest_bcs_t0 ()
      use dynkernel_options
      use dyn_fisl_options
      use lam_options
      use mem_nest
      use mem_tracers
      use gmm_vt0
      use glb_ld
      use tr3d
      implicit none
#include <arch_specific.hf>

      integer i,j,n,deb
      logical :: using_qt0
!
!----------------------------------------------------------------------
!
      using_qt0 = ( .not.Dynamics_hydro_L ) .or. (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H')

      if (l_north) then
         ut0 (1:l_niu,l_nj-pil_n+1:l_nj ,1:G_nk) = nest_u (1:l_niu,l_nj-pil_n+1:l_nj ,1:G_nk)
         vt0 (1:l_ni ,l_nj-pil_n  :l_njv,1:G_nk) = nest_v (1:l_ni ,l_nj-pil_n  :l_njv,1:G_nk)
         tt0 (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk) = nest_t (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk)
         st0 (1:l_ni ,l_nj-pil_n+1:l_nj) = nest_s (1:l_ni,l_nj-pil_n+1:l_nj)
         wt0 (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk) = nest_w (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk)
         zdt0(1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk) = nest_zd(1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk)
         if ( using_qt0 ) then
            qt0 (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk+1) = nest_q (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk+1)
         end if
      end if

      if (l_south) then
         ut0 (1:l_niu,1:pil_s ,1:G_nk) = nest_u (1:l_niu,1:pil_s ,1:G_nk)
         vt0 (1:l_ni ,1:pil_s ,1:G_nk) = nest_v (1:l_ni ,1:pil_s ,1:G_nk)
         tt0 (1:l_ni ,1:pil_s ,1:G_nk) = nest_t (1:l_ni ,1:pil_s ,1:G_nk)
         st0 (1:l_ni ,1:pil_s) = nest_s (1:l_ni,1:pil_s)
         wt0 (1:l_ni ,1:pil_s ,1:G_nk) = nest_w (1:l_ni ,1:pil_s ,1:G_nk)
         zdt0(1:l_ni ,1:pil_s ,1:G_nk) = nest_zd(1:l_ni ,1:pil_s ,1:G_nk)
         if ( using_qt0 ) then
            qt0 (1:l_ni ,1:pil_s ,1:G_nk+1) = nest_q (1:l_ni ,1:pil_s ,1:G_nk+1)
         end if
      end if

      if (l_east) then
         ut0 (l_ni-pil_e  :l_niu,1:l_nj ,1:G_nk) = nest_u (l_ni-pil_e  :l_niu,1:l_nj ,1:G_nk)
         vt0 (l_ni-pil_e+1:l_ni ,1:l_njv,1:G_nk) = nest_v (l_ni-pil_e+1:l_ni ,1:l_njv,1:G_nk)
         tt0 (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk) = nest_t (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk)
         st0 (l_ni-pil_e+1:l_ni ,1:l_nj) = nest_s (l_ni-pil_e+1:l_ni,1:l_nj)
         wt0 (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk) = nest_w (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk)
         zdt0(l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk) = nest_zd(l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk)
         if ( using_qt0 ) then
            qt0 (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk+1) = nest_q (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk+1)
         end if
      end if

      if (l_west) then
         ut0 (1:pil_w, 1:l_nj , 1:G_nk) = nest_u (1:pil_w, 1:l_nj , 1:G_nk)
         vt0 (1:pil_w, 1:l_njv, 1:G_nk) = nest_v (1:pil_w, 1:l_njv, 1:G_nk)
         tt0 (1:pil_w, 1:l_nj , 1:G_nk) = nest_t (1:pil_w, 1:l_nj , 1:G_nk)
         st0 (1:pil_w, 1:l_nj) = nest_s (1:pil_w,1:l_nj)
         wt0 (1:pil_w, 1:l_nj , 1:G_nk) = nest_w (1:pil_w, 1:l_nj , 1:G_nk)
         zdt0(1:pil_w, 1:l_nj , 1:G_nk) = nest_zd(1:pil_w, 1:l_nj , 1:G_nk)
         if ( using_qt0 ) then
            qt0 (1:pil_w, 1:l_nj, 1:G_nk+1) = nest_q (1:pil_w, 1:l_nj, 1:G_nk+1)
         end if
      end if

      if (Schm_opentop_L) then
         ut0 (1:l_niu,1:l_nj ,1:Lam_gbpil_t) = nest_u (1:l_niu,1:l_nj ,1:Lam_gbpil_t)
         vt0 (1:l_ni ,1:l_njv,1:Lam_gbpil_t) = nest_v (1:l_ni ,1:l_njv,1:Lam_gbpil_t)
         tt0 (1:l_ni ,1:l_nj ,1:Lam_gbpil_t-1) = nest_t (1:l_ni ,1:l_nj ,1:Lam_gbpil_t-1)
         wt0 (1:l_ni ,1:l_nj ,1:Lam_gbpil_t-1) = nest_w (1:l_ni ,1:l_nj ,1:Lam_gbpil_t-1)
         zdt0(1:l_ni ,1:l_nj ,1:Lam_gbpil_t-1) = nest_zd(1:l_ni ,1:l_nj ,1:Lam_gbpil_t-1)
         if ( using_qt0 ) then
            qt0 (1:l_ni,1:l_nj,1:Lam_gbpil_t)     = nest_q (1:l_ni,1:l_nj,1:Lam_gbpil_t)
         end if
      end if

      if (Lam_toptt_L) then
!        Pilot the temperature for the whole top level
         do j=1,l_nj
         do i=1,l_ni
            tt0(i,j,1) = nest_t(i,j,1)
         end do
         end do
      end if

      do n=1,Tr3d_ntr
         deb= (n-1)*G_nk + 1

         if (l_north) &
         tracers_M(n)%pntr (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk) = &
         nest_tr(1:l_ni,l_nj-pil_n+1:l_nj,deb:deb+G_nk-1)

         if (l_east ) &
         tracers_M(n)%pntr (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk) = &
         nest_tr(l_ni-pil_e+1:l_ni,1:l_nj,deb:deb+G_nk-1)

         if (l_south) &
         tracers_M(n)%pntr (1:l_ni ,1:pil_s ,1:G_nk) = &
         nest_tr (1:l_ni,1:pil_s,deb:deb+G_nk-1)

         if (l_west ) &
         tracers_M(n)%pntr (1:pil_w ,1:l_nj ,1:G_nk) = &
         nest_tr (1:pil_w,1:l_nj,deb:deb+G_nk-1)

         if (Schm_opentop_L) then
            tracers_M(n)%pntr (1:l_ni, 1:l_nj, 1:Lam_gbpil_t-1) = &
            nest_tr (1:l_ni,1:l_nj,deb:deb+Lam_gbpil_t-2)
         end if

      enddo
!
!----------------------------------------------------------------------
!
      return
      end
