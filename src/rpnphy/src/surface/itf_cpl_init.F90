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

function itf_cpl_init(F_path_S, F_print_L, F_unout, F_dateo, F_dt) result(F_istat)
   use cpl_itf, only: cpl_init
   use phygridmap, only: drv_glb_ni, drv_glb_nj, drv_lcl_ni, drv_lcl_nj, &
        phy_lcl_i0, phy_lcl_j0, phy_lcl_in, phy_lcl_jn, phydim_nk
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>

   character(len=*), intent(in) :: F_path_S
   logical, intent(in)          :: F_print_L
   integer, intent(in)          :: F_dateo, F_unout
   real, intent(in)             :: F_dt
   !@return
   integer :: F_istat
   !@object     interface to cpl_init
   !@authors    Francois Roy -- spring 2014
   !@revision
   ! v4_70 - Roy, F.  - initial version

   integer :: istat
   !---------------------------------------------------------------
   F_istat = 0
   if (.not.cplocn) return

   istat = cpl_init(F_path_S, F_print_L, F_unout, F_dateo, F_dt, &
        drv_glb_ni, drv_glb_nj, drv_lcl_ni, drv_lcl_nj, &
        phy_lcl_i0, phy_lcl_j0, phy_lcl_in, phy_lcl_jn, &
        phydim_nk, z0mtype, z0ttype, Z0TLAT)

   if (istat < 0) then
      F_istat = -1
   else
      cplocn = (istat == 1)
      F_istat = 0
   endif
   !---------------------------------------------------------------
   return
end function itf_cpl_init
