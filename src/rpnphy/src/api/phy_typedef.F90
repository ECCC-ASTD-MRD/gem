module phy_typedef
   use phy_status, only: PHY_ERROR, PHY_NONE, PHY_OK
   use phymem, only: phymeta, PHY_NAMELEN, PHY_MAXFLAGS, PHY_STAG_SFC, PHY_STAG_MOM, PHY_STAG_THERMO, PHY_STAG_ENERGY, PHY_NPATH_DEFAULT, PHY_BPATH_DEFAULT
   implicit none
   private
   public :: PHY_MAXNAMELENGTH, PHY_MAXFLAGS, phymeta, NPATH_DEFAULT, BPATH_DEFAULT, PHY_ERROR, PHY_NONE, PHY_OK

   integer, parameter :: PHY_MAXNAMELENGTH = PHY_NAMELEN
   
   character(len=*), parameter :: NPATH_DEFAULT = PHY_NPATH_DEFAULT
   character(len=*), parameter :: BPATH_DEFAULT = PHY_BPATH_DEFAULT

   ! iname : The input FST name of the field (referred to as 'I' throughout the API)
   ! oname : The output FST name of the field (referred to as 'O' throughout the API).
   ! vname : The internal variable name used in the physics (referred to as 'V' throughout the API).
   ! bus   : The name of the bus on which the field resides ('D' dynamics; 'P' permanent/physics; 'V' volatile; 'E' entry).
   ! init  : Boolean to define whether or not the field is initialized to 0 on the bus.
   ! stag  : Level staggering information (dynamics or energy, not Charney-Phillips staggering).
   ! esp   : The bus space allocated for this field.
   ! fmul  : The number of arbitrarily-defined levels in the field.
   ! nk    : The number of atmospheric levels.
   ! mosaic: The number of surface sub-types applicable to the field.
   ! wload : A flag defining whether or not the field should be used in the density calculation (water loaded).
   ! hzd   : tracer attribute, a flag defining whether or not the perform horizontal diffusion
   ! monot : tracer attribute, a flag defining whether or not the use monoton interpolation in the advection
   ! massc : tracer attribute, a flag defining whether or not the use a mass conservation scheme in the advection
   ! n     : A length three vector that defines the extent of each dimension of the field (shape).
   ! flags : list of keywords associated with the field

contains

   !# This s/r is needed for the compiler to produce a .o of this file (make requirement)
   subroutine phy_typedef_dummy()
      print *,'dummy'
   end subroutine phy_typedef_dummy

end module phy_typedef
