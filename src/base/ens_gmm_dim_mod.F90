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

module ens_gmm_dim
   use rmn_gmm
   implicit none
   public
   save

   type(gmm_metadata) :: meta3d_sh2 ! sans hallo
   type(gmm_metadata) :: meta3d_ar_p,meta2d_ar_s,meta3d_ai_p,meta2d_ai_s
   type(gmm_metadata) :: meta3d_br_p,meta2d_br_s,meta3d_bi_p,meta2d_bi_s
   type(gmm_metadata) :: meta2d_dum,meta3d_pls,meta4d_plp

end module ens_gmm_dim
