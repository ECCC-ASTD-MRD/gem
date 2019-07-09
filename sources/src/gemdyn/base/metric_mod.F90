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

module metric
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   real, dimension(:,:,:), pointer :: zmom, ztht, lg_pstar

   real, dimension(:,:,:), pointer :: mc_Ix, mc_Iy, mc_Iz
   real, dimension(:,:,:), pointer :: mc_Jx,  mc_Jy,  mc_iJz
   real, dimension(:,:,:), pointer :: mc_logJz

   real(kind=REAL64), dimension(:,:)  , pointer :: mc_css_H_8, mc_alfas_H_8, mc_betas_H_8, mc_cssp_H_8

end module metric
