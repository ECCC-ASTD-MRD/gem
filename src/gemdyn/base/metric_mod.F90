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

   type Vmetric
      real(kind=REAL64), dimension(:,:,:), allocatable :: zmom_8, ztht_8, lg_pstar_8
      real(kind=REAL64), dimension(:,:,:), allocatable :: zmom_u, ztht_u, zmom_v, ztht_v
      
      real(kind=REAL64), dimension(:,:,:), allocatable :: mc_Ix_8, mc_Iy_8, mc_Iz_8
      real(kind=REAL64), dimension(:,:,:), allocatable :: mc_Jx_8,  mc_Jy_8,  mc_iJz_8
      real(kind=REAL64), dimension(:,:,:), allocatable :: mc_logJz_8
      
      real(kind=REAL64), dimension(:,:)  , allocatable :: mc_css_H_8, mc_alfas_H_8, mc_betas_H_8, mc_cssp_H_8
      real(kind=REAL64), dimension(:,:)  , allocatable :: mc_cst_8, mc_alfat_8, mc_cstp_8
   end type Vmetric
      
   type(Vmetric) :: GVM
      
end module metric
