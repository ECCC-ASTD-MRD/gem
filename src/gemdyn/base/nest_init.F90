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
!**s/r nest_init -- Initializes nesting data for LAM configuration

      subroutine nest_init ()
      use lam_options
      use glb_ld
      use gmm_vt1
      use mem_nest
      use gmm_geof
      use lun
      use step_options
      use tr3d
      use mem_tracers
      implicit none
#include <arch_specific.hf>

      integer n,yy,mo,dd,hh,mm,ss,dum,deb
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write(Lun_out,1000)

      if (Lam_ctebcs_L) then
!     LAM with same (constant) pilot conditions
         nest_u  = ut1
         nest_v  = vt1
         nest_t  = tt1
         nest_s  = st1
         nest_w  = wt1
         nest_q  = qt1
         nest_zd = zdt1
         nest_fullme = fis0

         do n=1,Tr3d_ntr
            deb = (n-1) * l_nk
            nest_tr(l_minx:l_maxx,l_miny:l_maxy,deb+1:deb+l_nk) = tracers_P(n)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)
         end do

      else
!     ordinary LAM with future pilot conditions
         nest_u_fin  = ut1
         nest_v_fin  = vt1
         nest_t_fin  = tt1
         nest_s_fin  = st1
         nest_w_fin  = wt1
         nest_q_fin  = qt1
         nest_zd_fin = zdt1
         nest_fullme_fin = fis0
         do n=1,Tr3d_ntr
            deb = (n-1) * l_nk + 1
            nest_tr_fin(l_minx:l_maxx,l_miny:l_maxy,deb:deb+l_nk-1) = tracers_P(n)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)
         end do

      endif

      Lam_previous_S= ''
      Lam_current_S = Step_runstrt_S
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,Lam_current_S)
      call pdfjdate2 (Lam_tfin, yy,mo,dd,hh,mm,ss)
      Lam_tdeb      = Lam_tfin
!
!     ---------------------------------------------------------------
!
 1000 format(3X,'NESTING INITIALIZATION (NEST_INIT)')
      return
      end
