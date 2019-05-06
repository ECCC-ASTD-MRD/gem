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

module tr3d
   implicit none
   public
   save

   !-------------------------------------------------------------------------------------
   !Tr3d_mono = Type of monotonicity
   !          0 -> None
   !          1 -> Clipping QMSL : Bermejo-Staniforth,1992,MWR,120,2622-2631
   !          2 -> ILMC          : Sorenson et al.,2013,GMD,6,1029-1042
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   !Tr3d_mass = Type of mass conservation
   !          0                    -> None
   !          1                    -> Bermejo-Conde LEGACY=T : Bermejo-Conde,2002,MWR,130,423-430
   !          100+weight*10+pexp_n -> Bermejo-Conde LEGACY=F : Bermejo-Conde,2002,MWR,130,423-430
   !-------------------------------------------------------------------------------------
   !   weight = Weight in Bermejo-Conde
   !          1 -> Additive
   !          2 -> Multiplicative
   !          3 -> Additive+Factor pr_k/pr_s
   !   pexp_n = Rank in Adv_BC_pexp_list(PEXP_LIST_MAX) to choose P exponent REAL
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   !Tr3d_intp = NONE   : H:NONE  V:NONE (No advection)
   !Tr3d_intp = CUBIC  : H:Cubic V:Cubic
   !Tr3d_intp = QUINTIC: H:Cubic V:Quintic
   !-------------------------------------------------------------------------------------

   integer, parameter :: MAXTR3D = 250

   character(len=4) :: Tr3d_name_S(MAXTR3D), NTR_Tr3d_name_S(MAXTR3D)
   character(len=4) :: Tr3d_NT_S(MAXTR3D),Tr3d_M1_S(MAXTR3D),Tr3d_MC_S(MAXTR3D)
   character(len=12):: Tr3d_intp  (MAXTR3D)
   logical :: Tr3d_wload(MAXTR3D), Tr3d_hzd (MAXTR3D)
   integer :: Tr3d_mono (MAXTR3D), Tr3d_mass(MAXTR3D)
   integer :: Tr3d_ntr, Tr3d_ntrNT, Tr3d_ntrM1, Tr3d_ntrMC, NTR_Tr3d_ntr
   real    :: Tr3d_vmin(MAXTR3D),Tr3d_vmax(MAXTR3D)

end module tr3d
