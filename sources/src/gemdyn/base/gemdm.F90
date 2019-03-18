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

!**s/r gemdm - Main entry point for the GEMDM model
!
! Formal scientific documentation for this version can be found at
! $gemdyn/share/doc/GEM4.4.pdf
!
      subroutine gemdm
      implicit none
#include <arch_specific.hf>
!
!     ---------------------------------------------------------------
!
! Initialize: Domain, MPI, processor topology and ptopo.cdk

      call init_component

! Establish: model configuration, domain decomposition
!            and model geometry

      call set_world_view

! Initialize the ensemble prevision system

      call itf_ens_init

! Initialize the physics parameterization package

      call itf_phy_init

! Initialize tracers

      call tracers

! Setup main memory

      call main_gmm_storage

      call set_dyn_opr

! Run GEM

      call gem_ctrl

! Terminate

      call stop_world_view
!
!     ---------------------------------------------------------------
!
      return
      end
