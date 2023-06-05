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
module hvdif_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Horizontal diffusion if activated, will be applied to the following
   !# variables  : Horizontal winds, ZDot, W, tracers (_tr var)

   !# Background 2 delta-x removal ratio - range(0.0-1.0)
   real :: Hzd_lnR = -1.
   namelist /hvdif  / Hzd_lnr
   namelist /hvdif_p/ Hzd_lnr

   !# Theta 2 delta-x removal ratio - range(0.0-1.0).
   real :: Hzd_lnr_theta = -1.
   namelist /hvdif  / Hzd_lnr_theta
   namelist /hvdif_p/ Hzd_lnr_theta

   !# Tracers 2 delta-x removal ratio - range(0.0-1.0)
   real :: Hzd_lnR_tr = -1.
   namelist /hvdif  / Hzd_lnr_tr
   namelist /hvdif_p/ Hzd_lnr_tr

   !# Order of the background diffusion operator
   !# 2, 4, 6, 8
   integer :: Hzd_pwr = -1
   namelist /hvdif  / Hzd_pwr
   namelist /hvdif_p/ Hzd_pwr

   !# Order of the background diffusion operator on theta
   !# 2, 4, 6, 8
   integer :: Hzd_pwr_theta = -1
   namelist /hvdif  / Hzd_pwr_theta
   namelist /hvdif_p/ Hzd_pwr_theta

   !# Order of the background diffusion operator on tracers
   integer :: Hzd_pwr_tr = -1
   namelist /hvdif  / Hzd_pwr_tr
   namelist /hvdif_p/ Hzd_pwr_tr

   !# Main Smagorinsky control parameter (usual range 0.1-0.3)
   real :: Hzd_smago_param= -1.
   namelist /hvdif  / Hzd_smago_param
   namelist /hvdif_p/ Hzd_smago_param

   !# Apply Smago diffusion on theta using Hzd_smago_param/Hzd_smago_prandtl parameter
   real :: Hzd_smago_prandtl = -1.
   namelist /hvdif  / Hzd_smago_prandtl
   namelist /hvdif_p/ Hzd_smago_prandtl

   !# Apply Smago diffusion on HU using Hzd_smago_param/Hzd_smago_prandtl_hu parameter
   real :: Hzd_smago_prandtl_hu = -1.
   namelist /hvdif  / Hzd_smago_prandtl_hu
   namelist /hvdif_p/ Hzd_smago_prandtl_hu

   !# Coefficient of background diffusion added to the coefficient computed using the
   !# Smagorinsky approach. Two options are available.
   !# (i)  Fixed background diffusion: Just assign a single value, e.g., 0.2, and the resultant
   !#      diffusion will be similar to standard del-2 diffusion with Hzd_lnr=0.2.
   !# (ii) Vertically variable diffusion: Takes three values. The first element of the array determines
   !#      the constant value of background diffusion coeff. below Hzd_smago_lev(1). The second
   !#      the value at Hzd_smago_lev(2). The third element determines the maximum coefficient
   !#      element represents at the model top. Two ramps of COS^2-type are used between
   !#      Hzd_smago_lev(1) and Hzd_smago_lev (2), and between Hzd_smago_lev(2) and the model lid.
   real, dimension(3):: Hzd_smago_lnr = [ -1., -1., -1. ]
   namelist /hvdif  / Hzd_smago_lnr
   namelist /hvdif_p/ Hzd_smago_lnr

   !# The levels (bot,top) in the hybrid coordinate where the background diffusion
   !# coefficient varies between the value defined by Hzd_smago_lnr(1) and Hzd_smago_lnr(2).
   !# In GEM-P: the values should be in the units of hyb, e.g. Hzd_smago_lev = 0.7, 0.4
   !# In GEM-H: the values should be in the units of hyb_H, e.g. Hzd_smago_lev = 3000., 7500.

   real, dimension(2):: Hzd_smago_lev = [ -1., -1. ]
   namelist /hvdif  / Hzd_smago_lev
   namelist /hvdif_p/ Hzd_smago_lev

   !# If TRUE then background diffusion is applied to THETA and HU.
   logical :: Hzd_smago_theta_base_L = .false.
   namelist /hvdif  / Hzd_smago_theta_base_L
   namelist /hvdif_p/ Hzd_smago_theta_base_L

   !# Frictional heating is considered when Hzd_smago_fric_heat>0.
   real :: Hzd_smago_fric_heat = 0.
   namelist /hvdif  / Hzd_smago_fric_heat
   namelist /hvdif_p/ Hzd_smago_fric_heat

   !# Coefficients that multiply KM to simulate sponge layer near the top
   !# of the model. Warning! if this parameter is used, the EPONGE in the
   !# physics namelist should be removed.
   real, dimension(1000) :: Eq_sponge = 0.
   namelist /hvdif/ Eq_sponge

   !# Latitudinal ramping of equatorial sponge
   logical :: Eq_ramp_L = .false.
   namelist /hvdif  / Eq_ramp_L
   namelist /hvdif_p/ Eq_ramp_L

   !# The following variables are to help compute latitudinal modulation of
   !# vertical diffusion coefficient on momentum by fitting a cubic btwn
   !# values P_lmvd_weigh_low_lat and P_lmvd_weigh_high_lat at latitudes
   !# P_lmvd_low_lat and P_lmvd_high_lat

   !# Multiplication factor of P_pbl_spng at latitude P_lmvd_high_lat
   real :: P_lmvd_weigh_high_lat = 1.0
   namelist /hvdif  / P_lmvd_weigh_high_lat
   namelist /hvdif_p/ P_lmvd_weigh_high_lat

   !# Multiplication factor of P_pbl_spng at latitude P_lmvd_low_lat
   real :: P_lmvd_weigh_low_lat = 1.0
   namelist /hvdif  / P_lmvd_weigh_low_lat
   namelist /hvdif_p/ P_lmvd_weigh_low_lat

   !# Latitude at which the multiplication factor becomes P_lmvd_weigh_high_lat
   real :: P_lmvd_high_lat = 30.0
   namelist /hvdif  / P_lmvd_high_lat
   namelist /hvdif_p/ P_lmvd_high_lat

   !# latitude at which the multiplication factor becomes P_lmvd_weigh_low_lat
   real :: P_lmvd_low_lat =  5.0
   namelist /hvdif  / P_lmvd_low_lat
   namelist /hvdif_p/ P_lmvd_low_lat

   !# Vspng:
   !# Vertical sponge if activated, will be applied to the following
   !# variables: Horizontal Wind,  Temperature (Top level), Zdot, W

   !# Top coefficient for del-2 diffusion (m2/s)
   real :: Vspng_coeftop = -1.
   namelist /hvdif  / Vspng_coeftop
   namelist /hvdif_p/ Vspng_coeftop

   !# Number of levels from the top of the model
   integer :: Vspng_nk = 0
   namelist /hvdif  / Vspng_nk
   namelist /hvdif_p/ Vspng_nk

   !# True-> Riley diffusion on vertical motion on Vspng_nk levels
   logical :: Vspng_riley_L= .false.
   namelist /hvdif  / Vspng_riley_L
   namelist /hvdif_p/ Vspng_riley_L

   integer Eq_nlev
   integer Vspng_niter
   real, dimension(:,:), allocatable :: eponmod
   real, dimension(:  ), allocatable :: coef,cm,cp
   real(kind=REAL64), dimension(:), allocatable :: Vspng_coef_8

contains

!**s/r hvdif_nml - Read namelist hvdif

      integer function hvdif_nml (F_unf)
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         hvdif_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=hvdif_p)
            if ( F_unf == -2 ) write (Lun_out,nml=hvdif)
         end if
         return
      end if

      hvdif_nml= -1 ; nml_must= .false. ; nml_S= 'hvdif'

      rewind(F_unf)
      read (F_unf, nml=hvdif, end= 1001, err=1003)
      hvdif_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         hvdif_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (hvdif_nml < 0 ) return
      if ((Lun_out>=0).and.(hvdif_nml==0)) write (Lun_out, 6004) trim(nml_S)
      hvdif_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function hvdif_nml

   function hvdif_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function hvdif_options_init

end module hvdif_options
