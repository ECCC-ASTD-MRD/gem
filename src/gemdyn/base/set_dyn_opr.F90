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

!**   s/r set_dyn_opr - initialize operators and some constant parameters

      subroutine set_dyn_opr()
      use adz_mem
      use mem_tstp
      use tr3d
      use dynkernel_options
      use geomh
      use glb_ld
      implicit none

      integer :: dim,ntr
!
!     ---------------------------------------------------------------
!
!     Initialize horizontal diffusion package

      call hzd_exp_set()

!     Initialize DCMIP vertical diffusion package

      call dcmip_vrd_set()

!     Initialize common block for vertical sponge

      call vspng_set()

!     Initialize common block for equatorial sponge

      call eqspng_set()

      if ( Dynamics_FISL_L ) call adz_set () 

      call grid_area_mask (geomh_area_8,geomh_mask_8,geomh_area_mask_8,l_ni,l_nj)

      call adz_check_tracers()

      ! Common memory arena for temporary computations within OMP parallel regions
      if (Dynamics_hauteur_L) then
         dim= max(Adz_nij,(l_maxx-l_minx+1)*(l_maxy-l_miny+1)) * G_nk
         ntr= max(3, Tr3d_ntrTRICUB_NT, Tr3d_ntrTRICUB_WP, Tr3d_ntrBICHQV_NT, Tr3d_ntrBICHQV_WP)
         allocate ( WS1_8 (dim), WS1 (ntr*dim))
      endif
!     
!     ---------------------------------------------------------------
!
      return
      end
