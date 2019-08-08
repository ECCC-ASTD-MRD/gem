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
      use dynkernel_options
      use geomh
      use glb_ld
      implicit none
#include <arch_specific.hf>
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

      if ( Dynamics_FISL_L ) then
         call adv_setgrid()
         call adv_param()
         call adz_set ()
      end if

      call grid_area_mask (geomh_area_8, geomh_mask_8, l_ni, l_nj)

      call adz_check_tracers()
!
!     ---------------------------------------------------------------
!
      return
      end
