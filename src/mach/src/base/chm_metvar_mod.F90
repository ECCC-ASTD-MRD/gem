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
! Fichier/File   : chm_metvar_mod.ftn90
! Creation       : H. Landry, Oct. 2011
! Description    : Defines tables where the meteorological variables are
!                  dumped into so that they can be easily be accessed from
!                  within the chemistry.  It also obfuscates the buses below
!                  chm_exe.
!
!============================================================================

module chm_metvar_mod
   save
   
! 2D array
! Stuff that comes from the dynamic bus
   integer(kind=4), parameter :: MV2D_PPLUS    = 1  ! Surface pressure

! Stuff that comes from the permanent bus
   integer(kind=4), parameter :: MV2D_DXDY     = 2  ! Area of grid square
   integer(kind=4), parameter :: MV2D_TSURF    = 3  ! Area averaged surface temperature
   integer(kind=4), parameter :: MV2D_WSDIAG   = 4  ! Screen level wind speed
   integer(kind=4), parameter :: MV2D_TDIAG    = 5  ! Screen level specific humidity
   integer(kind=4), parameter :: MV2D_QDIAG    = 6  ! Screen level specific humidity
   integer(kind=4), parameter :: MV2D_GLSEA    = 7  ! Sea ice fraction
   integer(kind=4), parameter :: MV2D_SNODP    = 8  ! Snow Depth
   integer(kind=4), parameter :: MV2D_H        = 9  ! Boundary layer height
   integer(kind=4), parameter :: MV2D_DLAT     =10  ! Latitude
   integer(kind=4), parameter :: MV2D_DLON     =11  ! Longitude
   integer(kind=4), parameter :: MV2D_FLUSOLIS =12  ! VIS flux towards ground
   integer(kind=4), parameter :: MV2D_MT       =13  ! Topography
   integer(kind=4), parameter :: MV2D_ILMO     =14  ! Aggregated Inverse of Monin-Obukhov length
   integer(kind=4), parameter :: MV2D_WSOIL    =15  ! Soil volumetric water content
   integer(kind=4), parameter :: MV2D_MG       =16  ! land-sea mask (0-1 fraction)
   integer(kind=4), parameter :: MV2D_AL5      =17  ! Albedo-ALVIS visible surface total (5th dimension)

! Stuff that comes from the volatile bus
   integer(kind=4), parameter :: MV2D_UE       = 18 ! Friction velocity u* from surface momentum flux
   integer(kind=4), parameter :: MV2D_CANG     = 19 ! Cosine of solar zenith
   integer(kind=4), parameter :: MV2D_RAINRATE = 20 ! liquid precip. rate
   integer(kind=4), parameter :: MV2D_SNOF     = 21 ! Snow fraction
   integer(kind=4), parameter :: MV2D_WDDIAG   = 22 ! Screen level wind direction

   integer(kind=4), parameter :: SIZE_MV2D     = 22
!
! 3D array
! Stuff that comes (or evalated from fields) from the dynamic bus
   integer(kind=4), parameter :: MV3D_TPLUS    = 1 ! Temperature
   integer(kind=4), parameter :: MV3D_WS       = 2 ! Wind speed at t+dt(dynamics)
   integer(kind=4), parameter :: MV3D_RHO      = 3 ! Air density
   integer(kind=4), parameter :: MV3D_HUPLUS   = 4 ! Specific humidty at t+dt(dynamics)
   integer(kind=4), parameter :: MV3D_QCPLUS   = 5 ! Cloud water/ice at t+dt(dynamics)
   integer(kind=4), parameter :: MV3D_SIGM     = 6 ! local sigma values for momentum levels
   integer(kind=4), parameter :: MV3D_SIGT     = 7 ! local sigma values for thermodynamic levels
   integer(kind=4), parameter :: MV3D_WPLUS    = 8 ! dz/dt (vertical velocity)

! Stuff that comes from the permanent bus
   integer(kind=4), parameter :: MV3D_FTOT     = 9 ! Total cloud (0-1)

! Stuff that comes from the volatile bus
   integer(kind=4), parameter :: MV3D_ZMOM     = 10 ! Geop. heights on momentum levels
   integer(kind=4), parameter :: MV3D_ZPLUS    = 11 ! Geop. heights on thermodynamic levels
   integer(kind=4), parameter :: MV3D_KT       = 12 ! Thermal vertical diffusion coefficient

! Stuff that comes from both volatile and permanent bus
   integer(kind=4), parameter :: MV3D_RNFLX    = 13 ! Liquid precipitation flux (stratiform + convective)
   integer(kind=4), parameter :: MV3D_SNOFLX   = 14 ! Solid precipitation flux (stratiform + convective)
   integer(kind=4), parameter :: MV3D_PEVP     = 15 ! Evap. of precip (stratiform + convective)
   integer(kind=4), parameter :: MV3D_PPRO     = 16 ! Cloud to rain collection tendency (stratiform + convective)
! Used in MESSY rad. transfer
   integer(kind=4), parameter :: MV3D_CLDRAD   = 17 ! CLOUD fraction for radiation from volatile bus (see. cloud in "prep_cw_rad.F90")
   integer(kind=4), parameter :: MV3D_LWC      = 18 ! CLOUD/WATER CONTENT AT TIME T (GEM:LWC/QD) from perm. bus
   
! Used for the indirect feedback effect   
   integer(kind=4), parameter :: MV3D_NCPLUS   = 19 ! Cloud droplets number concentration at t+dt(dynamics)

! LINOZ Ozone (from the dynamic bus)
   integer(kind=4), parameter :: MV3D_O3L      = 20 ! LINOZ Ozone (in ug/kg)

   integer(kind=4), parameter :: SIZE_MV3D     = 20
!  
end module chm_metvar_mod
