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

!**s/r set_gmm
!
      subroutine set_gmm
      use gem_options
      use glb_ld
      use gmm_itf_mod
      use rstr
      use var_gmm
      implicit none
#include <arch_specific.hf>
!
!author
!     Michel Desgagne  -  February 2010
!
!revision
! v4_12 - Desgagne M.       - initial version
!

      integer :: istat
!-------------------------------------------------------------------
!
! Establish meta_data template variables for gmm_create purpose
!
      call gmm_build_meta1D(metasc, &
                            1,1,0,0,1, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta1d, &
                            1,l_ni*l_nj*l_nk,0,0,l_ni*l_nj*l_nk, &
                            0,0,0,0,0, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta2d, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta3D(meta3d_nk, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            1,l_nk,0,0,l_nk, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta3D(meta3d_nk1, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            1,l_nk+1,0,0,l_nk+1, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta3D(meta3d_0nk, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            0,l_nk,0,0,l_nk+1, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta3D(meta3d_0nk1, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            0,l_nk+1,0,0,l_nk+2, &
                            0,GMM_NULL_FLAGS)

      if (Rstri_rstn_L) istat = gmm_checkpoint_all(GMM_READ_CKPT)
!
!-------------------------------------------------------------------
!
      return
      end
