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
   use, intrinsic :: iso_fortran_env
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
   !          0           -> None
   !          1           -> Bermejo-Conde LEGACY=T : Bermejo-Conde,2002,MWR,130,423-430
   ! 100+weight*10+pexp_n -> Bermejo-Conde LEGACY=F : Bermejo-Conde,2002,MWR,130,423-430
   !-------------------------------------------------------------------------------------
   !   weight = Weight in Bermejo-Conde
   !          1 -> Additive
   !          2 -> Multiplicative
   !          3 -> Additive+Factor pr_k/pr_s
   !   pexp_n = Rank in Adv_BC_pexp_list(PEXP_LIST_MAX) to choose P exponent REAL
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   !Tr3d_intp = 'NONE'    : Hor:NONE  Ver:NONE   (No advection)
   !Tr3d_intp = 'TRICUB'  : Hor:Cubic Ver:Cubic  (tricubic)
   !Tr3d_intp = 'BICUBH_QV: Hor:Cubic Ver:Quintic
   !-------------------------------------------------------------------------------------

   integer, parameter :: MAXTR3D = 250

   character(len=4) :: Tr3d_name_S(MAXTR3D), NTR_Tr3d_name_S(MAXTR3D)
   character(len=12):: Tr3d_intp  (MAXTR3D)
   logical :: Tr3d_wload(MAXTR3D), Tr3d_hzd (MAXTR3D)
   integer :: Tr3d_mono (MAXTR3D), Tr3d_mass(MAXTR3D)
   integer :: Tr3d_ntr, NTR_Tr3d_ntr, Tr3d_hu
   integer :: Tr3d_ntrTRICUB_NT, Tr3d_ntrTRICUB_WP
   integer :: Tr3d_ntrBICHQV_NT, Tr3d_ntrBICHQV_WP
   integer :: Tr3d_debTRICUB_NT, Tr3d_debTRICUB_WP
   integer :: Tr3d_debBICHQV_NT, Tr3d_debBICHQV_WP
   real    :: Tr3d_vmin(MAXTR3D),Tr3d_vmax(MAXTR3D)

   real(kind=REAL64), dimension(:), allocatable :: BCMD_3C, BCMD_BQ

   type :: meta_tracers
      character(len=4 ) :: name
      character(len=12) :: intp
      logical :: wload,hzd
      integer :: mono,mass
      real    :: vmin,vmax
      real(kind=REAL64) :: BC_mass_deficit
   end type meta_tracers

   type(meta_tracers), dimension(:), pointer :: Tr_3CNT, Tr_3CWP,&
                                                Tr_BQNT, Tr_BQWP

end module tr3d
