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
module hzd_mod
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

!______________________________________________________________________
!                                                                      |
!  PROJECTION OPERATORS FOR HOR. DIFFUSION  (initialized in hzd_set))  |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Hzd_geom*          | Hor. diffu. in the dynamics (rhs)               |
!----------------------------------------------------------------------
   real(kind=REAL64), dimension(:,:), pointer, contiguous :: Hzd_geom_q => null()
   real(kind=REAL64), dimension(:,:), pointer, contiguous :: Hzd_geom_u => null()
   real(kind=REAL64), dimension(:,:), pointer, contiguous :: Hzd_geom_v => null()
   real(kind=REAL64), dimension(:), allocatable   :: Hzd_smago_lnrM_8, Hzd_smago_lnrT_8
   real(kind=REAL64) , dimension(:), allocatable :: Hzd_coef_8,Hzd_coef_8_tr,Hzd_coef_8_theta
   
   real(kind=REAL64), dimension (:,:,:,:), allocatable ::stencil_V
   real(kind=REAL64), dimension (:,:,:), allocatable ::stencil_V1,stencil_V2 
 
   real(kind=REAL64), dimension (:,:,:,:), allocatable :: skpu, sku, skpv, skv 
   real(kind=REAL64), dimension (:,:,:), allocatable :: xfactu,xfactv,xfacth,xfactz 
   real(kind=REAL64), dimension (:,:,:,:), allocatable :: Jy, Jyp,Jizx, Jizxm,Jizpx, &
                                                          Jx, Jxp,Jzv,Jzvm
   real(kind=REAL64), dimension (:,:,:,:), allocatable :: Jzpt, Jzt, Jxt, Jyt
   real(kind=REAL64), dimension (:,:,:,:), allocatable :: Jm, Jzz
   real(kind=REAL64), dimension (:,:,:), allocatable :: a_u, b_u, c_u, W_u, &
                                                        a_v, b_v, c_v, W_v, &
                                                        a_th, b_th, c_th, W_th, &
                                                        a_zdt, b_zdt, c_zdt, W_zdt
   real(kind=REAL64) :: zfact 
   real(kind=REAL64) ru,rv

   integer      Hzd_niter,Hzd_niter_tr,Hzd_niter_theta
   integer      Hzd_hyb_top, hzd_hyb_bot

   real, dimension(:,:,:), allocatable :: air_density,air_density_U,air_density_V,air_density_m


end module hzd_mod
