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

!/@*
function itf_phy_prefold_opr (F_data, F_name_S, F_horiz_interp_S, &
             F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn) result(F_istat)
      use HORgrid_options
      use glb_ld
   implicit none

   !@objective - Pre-folding operations for data on the physics grid
   !@arguments
   character(len=*),intent(in) :: F_name_S,F_horiz_interp_S
   integer,intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn
   real,intent(inout) :: F_data(F_minx:F_maxx,F_miny:F_maxy,F_k0:F_kn)
   !@return
   integer :: F_istat
   !@author Stephane Chamberland, 2012-04
   !@revisions
   !*@/

#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   real, dimension(l_ni,l_nj,F_k0:F_kn) :: data_dyngrid

   !---------------------------------------------------------------------
   call msg(MSG_DEBUG,'[BEGIN] itf_phy_prefold_opr')
   F_istat = RMN_OK

   if (Grd_yinyang_L .and. F_horiz_interp_S /= ' ') then
       data_dyngrid = 0.
       data_dyngrid(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:) = F_data
       call yyg_scalgeo(data_dyngrid,1,l_ni,1,l_nj,F_kn-F_k0+1,F_horiz_interp_S)
       F_data = data_dyngrid(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)
   end if

   call msg(MSG_DEBUG,'[END] itf_phy_prefold_opr')
   !---------------------------------------------------------------------

   return
end function itf_phy_prefold_opr
