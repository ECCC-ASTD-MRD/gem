!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : chm_consphychm_mod.ftn90
! Creation       : H. Landry, Mai 2008
! Description    : Modules defining physical and chemical constants as
!                  provided in the file "constante" in the GEM directory
!
! Extra info     :
!
!============================================================================

module chm_consphychm_mod
   use tdpack_const
   save
!
!     declaration of chemical constants
!
   real(kind=4), parameter :: avno      = 0.6022000000000e+24 ! avogadro's num        atoms mol-1
   real(kind=4), parameter :: kboltz    = 0.1381000000000e-15 ! boltzmann constant    (in cgs)
   real(kind=4), parameter :: rg        = 0.8202714038983e-04 ! gas constant          atm -> pa
   real(kind=4), parameter :: rgasi     = 0.8314000000000e+01 ! gas constant (J K-1 mol-1)
   real(kind=4), parameter :: mwt_air   = 0.2897000000000e+02 ! mol. wgt. of air
   real(kind=4), parameter :: mwt_h2o   = 18.015              ! mol. wgt. of water
   real(kind=4), parameter :: rho_h2o   = 0.1000000000000e+07 ! density of water      g/m**3
   real(kind=4), parameter :: consth    = (mwt_air / mwt_h2o) * rho_h2o ! convert (kg H2O/kg air) to ppmv
   real(kind=4), parameter :: du_o3     = 0.2687000000000e+20 ! Loschmidth number(molec/cm3): O3 (molec/cm2) / du_o3 (molec/cm3) * 1E3 (DU/cm) = O3 (du)
   real(kind=4), parameter :: rad2deg   = 180.0 / pi          ! convert radians to degree

end module chm_consphychm_mod


