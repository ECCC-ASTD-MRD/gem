!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------


subroutine prphynml()
   implicit none
!!!#include <arch_specific.hf>
   !@author V.Lee - October 20,2008
   !----------------------------------------------------------------------
   print *,'creating phydict.nml'
   open(6, file="phydict.nml", access='SEQUENTIAL', form='FORMATTED')
   call nml_phy()
   call nml_convection()
   call nml_surface()
   call nml_series()
   close(6)
   !-------------------------------------------------------------------
   return
end subroutine prphynml

subroutine nml_phy()
   use phy_options
   implicit none
!!!#include <arch_specific.hf>
   integer :: istat
   istat = phy_options_init()
   write(6, nml=physics_cfgs)
   return
end subroutine nml_phy

subroutine nml_surface
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>
   integer :: istat
   istat = sfc_options_init()
   write(6, nml=surface_cfgs)
   return
end subroutine nml_surface

subroutine nml_convection
   use cnv_options
   implicit none
!!!#include <arch_specific.hf>
   integer :: istat
   istat = cnv_options_init()
   write(6, nml=convection_cfgs)
   return
end subroutine nml_convection

subroutine nml_series()
   use series_options
   implicit none
!!!#include <arch_specific.hf>
   integer :: istat
   istat = series_options_init()
   write(6, nml=series)
   return
end subroutine nml_series
