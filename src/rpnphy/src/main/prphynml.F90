

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
