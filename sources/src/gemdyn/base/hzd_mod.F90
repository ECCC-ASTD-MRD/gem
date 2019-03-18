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
   real*8, dimension(:,:), pointer, contiguous :: Hzd_geom_q => null()
   real*8, dimension(:,:), pointer, contiguous :: Hzd_geom_u => null()
   real*8, dimension(:,:), pointer, contiguous :: Hzd_geom_v => null()
   real*8, dimension(:), allocatable   :: Hzd_smago_lnrM_8, Hzd_smago_lnrT_8

   integer      Hzd_niter,Hzd_niter_tr,Hzd_niter_theta
   real*8 , dimension(:), allocatable :: Hzd_coef_8,Hzd_coef_8_tr,Hzd_coef_8_theta

end module hzd_mod
