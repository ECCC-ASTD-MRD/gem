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

!**s/r gem_ctrl - initiate the forward integration of the model

      subroutine gem_ctrl ()
      use dynkernel_options
      use dyn_fisl_options
      use gem_options
      use HORgrid_options
      use lam_options
      use step_options
      use init_options
      use glb_ld
      use gmm_geof
      use mem_nest
      use lun
      use rstr
      use gem_timing
      use tr3d
      implicit none

!object
!     Beginning of the integration. This subroutine
!     reads the data and performs initialization if required.
!     It then initiates the forward intergration of the model.

      logical :: rstrt_L= .false.
!     
!     ---------------------------------------------------------------
!
      call gemtime ( Lun_out, 'GEM_CTRL: START', .false. )

      if ( .not. Rstri_rstn_L ) then
         call indata()
      else
         if (Dynamics_hauteur_L) &
         call vertical_metric (GVM, fis0, sls, l_minx,l_maxx,l_miny,l_maxy)
         if ( .not. Grd_yinyang_L .and. Lam_ctebcs_L ) then
            call nest_indata  (nest_u, nest_v , nest_w, nest_t   ,&
                               nest_q, nest_zd, nest_s, nest_tr  ,&
                               nest_fullme,.false.,Step_runstrt_S,&
                               l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr)
         endif
      endif
      call glbstat ( fis0,'ME',"indata",l_minx,l_maxx,l_miny,l_maxy,1,1,&
                     1-G_halox,G_ni+G_halox,1-G_haloy,G_nj+G_haloy,1,1 )
      if (Schm_sleve_L) &
      call glbstat ( orols,'MELS',"indata",l_minx,l_maxx,l_miny,l_maxy, &
                 1,1, 1-G_halox,G_ni+G_halox,1-G_haloy,G_nj+G_haloy,1,1 )
                   
      call gemtime ( Lun_out, 'GEM_CTRL: INIT COMPLETED', .false. )
      call gemtime_stop ( 2 )

      if (  Init_mode_L ) call initial (rstrt_L)

      if ( .not.rstrt_L ) call gem_run (rstrt_L)

      if (Lun_out > 0) write(Lun_out,3000) Lctl_step

 3000 format(/,'GEM_CTRL: END OF CURRENT TIME SLICE AT TIMESTEP',I8, &
             /,'===================================================')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine gem_ctrl
