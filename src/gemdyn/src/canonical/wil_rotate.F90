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

!**s/r wil_rotate - Compute rotated coordinates rotlon,rotlat for a rotation by angle alpha

      subroutine wil_rotate (F_rlon_8,F_rlat_8,F_rotlon_8,F_rotlat_8,F_kind)

      use wil_options
      use tdpack

      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64) F_rlon_8, F_rlat_8, F_rotlon_8, F_rotlat_8

      integer F_kind

      !object
      !======================================================================================
      !     This subroutine computes the rotated coordinates
      !     rotlon,rotlat for a rotation by angle alpha
      !
      !     Based on Ritchie(1987),MWR,115,608-619 and Williamson et al.,1992,JCP,102,211-224
      !
      !     F_kind=0: From No rotation to    rotation
      !     F_kind=1: From    rotation to No rotation
      !======================================================================================

      !---------------------------------------------------------------

      real(kind=REAL64) theta_8,lambda_8,test_8,test1_8,test2_8,alpha_8

      real(kind=REAL64), parameter :: EPS_8 = 1.0d-06

      !---------------------------------------------------------------

      !--------------------------------------------------------------------------
      !Williamson_lon_pole_r_8 and Williamson_lon_pole_r_8 are Longitude,Latitude
      !in the no rotated system of the north pole P' of the rotated system
      !--------------------------------------------------------------------------

      theta_8  = Williamson_lat_pole_r_8
      lambda_8 = Williamson_lon_pole_r_8
      alpha_8  = Williamson_alpha

      !-----------
      !No rotation
      !-----------
      if (alpha_8==0.0d0) then

         F_rotlon_8 = F_rlon_8
         F_rotlat_8 = F_rlat_8

      !-----------------------
      !Rotation by angle alpha
      !-----------------------
      else

         !Latitude
         !--------
         if (F_kind==0) then
            test_8 = sin(F_rlat_8)*sin(theta_8)+cos(F_rlat_8)*cos(F_rlon_8-lambda_8)*cos(theta_8)
         else
            test_8 = sin(F_rlat_8)*sin(theta_8)-cos(F_rlat_8)*cos(F_rlon_8         )*cos(theta_8)
         endif

         if (test_8>1.0d0) then
            F_rotlat_8 =  pi_8/2.0
         elseif (test_8<-1.0d0) then
            F_rotlat_8 = -pi_8/2.0
         else
            F_rotlat_8 = asin(test_8)
         endif

         !Longitude
         !---------
         test_8 = cos(F_rotlat_8)

         if (test_8==0.0d0) then
            F_rotlon_8 = 0.0d0
         else

            if (F_kind==0) then
               test_8 = sin(F_rlon_8-lambda_8)*cos(F_rlat_8)/test_8
            else
               test_8 = sin(F_rlon_8         )*cos(F_rlat_8)/test_8
            endif

            if (test_8>1.0d0) then
               F_rotlon_8 =  pi_8/2.0
            elseif (test_8<-1.0d0) then
               F_rotlon_8 = -pi_8/2.0
            else
               if (F_kind==0) then
                  F_rotlon_8 = asin(test_8)
               else
                  F_rotlon_8 = asin(test_8) + lambda_8
               endif
            endif

         endif

         !Adjust for correct branch of inverse sine
         !-----------------------------------------
         if (F_kind==0) then
            test1_8 = sin(F_rlat_8)
            test2_8 = sin(F_rotlat_8)*sin(theta_8)-cos(F_rotlat_8)*cos(theta_8)*cos(F_rotlon_8)
            if (abs(test1_8-test2_8) > EPS_8) then
               F_rotlon_8 = pi_8 - F_rotlon_8
            endif
         else
            test1_8 = sin(F_rlat_8)
            test2_8 = sin(F_rotlat_8)*sin(theta_8)+cos(F_rotlat_8)*cos(theta_8)*cos(F_rotlon_8-lambda_8)
            if (abs(test1_8-test2_8) > EPS_8) then
               F_rotlon_8 = pi_8 - (F_rotlon_8-lambda_8) + lambda_8
            endif
         endif

      endif

      !---------------------------------------------------------------

      return
      end
