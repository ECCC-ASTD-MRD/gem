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

module adv_interp
   implicit none
   public
   save

!----------------------------------------------------------------------|
! LOCALISATION PARAMETERS                                              |
!----------------------------------------------------------------------|
! adv_ovdx_8   | \                                                     |
! adv_ovdy_8   |   inverse of shortest grid distance in x,y,z          |
! adv_ovdz_8   | /                                                     |
! adv_x00_8    | \                                                     |
! adv_y00_8    | / reference coordinate for localisation in x,y,z      |
! adv_lcx      | \                                                     |
! adv_lcy      |   fine grid to variable grid equivalence used for     |
! adv_lcz      | / fast localisation                                   |
!----------------------------------------------------------------------|
! PRECOMPUTED INTERPOLATION PARAMETERS                                 |
!----------------------------------------------------------------------|
! adv_dlx_8    |
! adv_dly_8    |
! adv_dlz_8    |
! adv_bsx_8    |
! adv_bsy_8    |
! adv_bsz_8    |
! adv_dix_8    |
! adv_diy_8    |
! adv_x_8      | Coefficients for linear interpolation in x,y and z
! adv_y_8      |
! adv_zbc_8    |
! adv_xabcd_8  | Coefficients for Lagrange 3D in x,y and z
! adv_xbacd_8  |
! adv_xcabd_8  |
! adv_xdabc_8  |
! adv_yabcd_8  |
! adv_ybacd_8  |
! adv_ycabd_8  |
! adv_ydabc_8  |
! adv_zabcd_8  |
! adv_zbacd_8  |
! adv_zcabd_8  |
! adv_zdabc_8  |
! adv_zxabcde_8| Coefficients for Quintic Lagrange in z
! adv_zaxbcde_8|
! adv_zbxacde_8|
! adv_zcxabde_8|
! adv_zdxabce_8|
! adv_zexabcd_8|
!----------------------------------------------------------------------|

   type :: ADV_TYPE_V_R8
      sequence
      real*8, dimension(:), pointer, contiguous :: t,m,x
   end type ADV_TYPE_V_R8

   type :: ADV_TYPE_V_I
      sequence
      integer, dimension(:), pointer, contiguous :: t,m,x
   end type ADV_TYPE_V_I

   type :: ADV_TYPE_SV_R8
      sequence
      real*8, dimension(:), pointer, contiguous :: t,m,x
   end type ADV_TYPE_SV_R8

   type :: ADV_TYPE_SV_I
      sequence
      integer, dimension(:), pointer, contiguous :: t,m,x
   end type ADV_TYPE_SV_I


   integer,dimension(:),allocatable :: adv_lcx,adv_lcy
   integer :: pnz

   real*8 ::  adv_ovdx_8, adv_ovdy_8, adv_ovdz_8
   real*8 ::  adv_x00_8,  adv_y00_8
   real*8 ::  adv_x_8, adv_y_8
   real*8 ::  adv_xabcd_8, adv_xbacd_8, adv_xcabd_8,adv_xdabc_8, adv_xbc_8
   real*8 ::  adv_yabcd_8, adv_ybacd_8, adv_ycabd_8,adv_ydabc_8, adv_ybc_8

   real*8, dimension(:), allocatable :: adv_dlx_8, adv_dly_8
   real*8, dimension(:), allocatable :: adv_bsx_8, adv_bsy_8

   real*8, dimension(:), allocatable :: adv_diz_8

   type(ADV_TYPE_V_R8)  :: adv_dlz_8
   type(ADV_TYPE_SV_R8) :: adv_zbc_8
   type(ADV_TYPE_SV_R8) :: adv_zabcd_8, adv_zbacd_8, adv_zcabd_8, adv_zdabc_8

   type(ADV_TYPE_SV_R8) :: adv_zxabcde_8, adv_zaxbcde_8, adv_zbxacde_8, &
                           adv_zcxabde_8, adv_zdxabce_8, adv_zexabcd_8

   type(ADV_TYPE_SV_R8) :: adv_bsz_8
   type(ADV_TYPE_SV_I)  :: adv_lcz

end module adv_interp
