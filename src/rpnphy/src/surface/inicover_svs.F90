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

module inicover_svs_mod
   implicit none
   private
   
   public :: inicover_svs
   
contains

subroutine inicover_svs(pvars, kount, ni)
   use mu_jdate_mod, only: jdate_day_of_year
   use sfc_options
   use sfcbus_mod
   use phymem, only: phyvar
   use lookup4ccilc_we_mod, only: lookup4ccilc_we
   use aggcovernat_mod, only: aggcovernat
   use aggveghigh_mod, only: aggveghigh
   use aggveglow_mod, only: aggveglow
   use interpveg_mod, only: interpveg
   use veglowhigh_mod, only: veglowhigh
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   type(phyvar), pointer, contiguous :: pvars(:)
   integer, intent(in) :: ni, kount

   !@Author S. Belair et  M. Abrahamowicz (Jan 2016)
   !@Revision
   !@Object Initialize vegetation fields for SVS scheme
   !@Arguments
   !       - Input/Ouput -
   ! pvars    list of all phy vars (meta + slab data)
   !       - Input -
   ! kount    current timestep number
   ! ni       horizontal slice dimension
   !
   !@Notes    inisurf has been split in two subroutines:
   !          inisurf and inicover. the former calls the latter.
   !
   !     the geophysical fields determined from vegetation
   !     are done so using the following classification:
   !
   !     Class       Vegetation type
   !     =====       ===============
   !       1         (salt) water
   !       2         ice
   !       3         inland lake
   !       4         evergreen needleleaf trees
   !       5         evergreen broadleaf trees
   !       6         deciduous needleleaf trees
   !       7         deciduous broadleaf trees
   !       8         tropical broadleaf trees
   !       9         drought deciduous trees
   !       10        evergreen broadleaf shrub
   !       11        deciduous shrubs
   !       12        thorn shrubs
   !       13        short grass and forbs
   !       14        long grass
   !       15        crops
   !       16        rice
   !       17        sugar
   !       18        maize
   !       19        cotton
   !       20        irrigated crops
   !       21        urban
   !       22        tundra
   !       23        swamp
   !       24        desert
   !       25        mixed wood forests
   !       26        mixed shrubs

   !********************************************************************
   ! Tables for the veg characteristics for each veg type
   !********************************************************************
      REAL ALDAT(NCLASS), D2DAT(NCLASS), D50DAT(NCLASS), D95DAT(NCLASS)
      REAL RSMINDAT(NCLASS), GEXPDAT(NCLASS)
      REAL LAIDAT(NCLASS), VEGDAT(NCLASS),EMISDAT(NCLASS) 
      REAL CVDAT(NCLASS), RGLDAT(NCLASS), GAMMADAT(NCLASS)
      REAL Z0MDAT(NCLASS), MAXPDAT(NCLASS)
!
      DATA ALDAT/ &
                    0.13   , 0.70   , 0.13   , 0.14   , 0.12   , &
                    0.14   , 0.18   , 0.13   , 0.17   , 0.14   , &
                    0.18   , 0.19   , 0.20   , 0.19   , 0.20   , & 
                    0.21   , 0.18   , 0.18   , 0.25   , 0.18   , & 
                    0.12   , 0.17   , 0.12   , 0.30   , 0.15   , &
                    0.15   / 
!    
      DATA D2DAT/    &
                    0.0    , 0.0    , 0.0    , 2.0    , 2.0    , &
                    1.0    , 2.0    , 2.0    , 2.0    , 2.0    , & 
                    2.0    , 2.0    , 1.5    , 2.0    , 2.0    , & 
                    1.0    , 1.0    , 1.5    , 2.0    , 1.5    , & 
                    1.0    , 1.0    , 2.0    , 2.0    , 2.0    , & 
                    2.0    / 

!    
      DATA D50DAT/    &
                    0.0    , 0.0    , 0.0    , 0.2    , 0.2    , &
                    0.2    , 0.2    , 0.2    , 0.2    , 0.2    , & 
                    0.2    , 0.3    , 0.5    , 0.2    , 0.5    , & 
                    0.2    , 0.2    , 0.2    , 0.2    , 0.2    , & 
                    0.2    , 0.1    , 0.15   , 0.75   , 0.2    , & 
                    0.5    / 
!    
      DATA D95DAT/    &
                    0.0    , 0.0    , 0.0    , 0.9    , 0.9    , &
                    0.9    , 0.9    , 1.2    , 0.9    , 0.9    , & 
                    0.9    , 1.5    , 1.5    , 1.2    , 1.5    , & 
                    0.9    , 0.9    , 0.9    , 0.9    , 0.9    , & 
                    0.9    , 0.3    , 0.7    , 2.0    , 0.9    , & 
                    1.5    / 
!    
      DATA RSMINDAT/    &
                    500.   , 500.   , 500.   , 250.   , 250.   , &
                    250.   , 250.   , 250.   , 250.   , 150.   , & 
                    150.   , 150.   , 100.   , 100.   , 100.   , & 
                    100.   , 100.   , 100.   , 100.   , 150.   , & 
                    150.   , 150.   , 150.   , 500.   , 250.   , & 
                    250.   / 
      DATA LAIDAT/ &
                    0.00   , 0.00   , 0.00   , 5.00   , 6.00   , & 
                   -99.    , -99.   , 6.00   , 4.00   , 3.00   , & 
                   -99.    , 3.00   , 1.00   , -99.   , -99.   , &
                   -99.    , -99.   , -99.   , -99.   , 1.00   , & 
                    1.00   , -99.   , 4.00   , 0.00   , -99.   , & 
                   -99.    / 
      DATA MAXPDAT/ &
                    0.00   , 0.00   , 0.00   , 0.10   , 0.10   , &
                    0.10   , 0.10   , 0.10   , 0.10   , 0.05   , &
                    0.05   , 0.05   , 0.05   , 0.05   , 0.05   , &
                    0.05   , 0.05   , 0.05   , 0.05   , 0.005  , &
                    0.005  , 0.05   , 0.05   , 0.05   , 0.10   , &
                    0.10   /
      DATA VEGDAT/ &
                    0.00   , 0.00   , 0.00   , 0.90   , 0.99   , & 
                    0.90   , 0.90   , 0.99   , 0.90   , 0.50   , & 
                    0.50   , 0.50   , 0.7    , 0.30   , -99.   , & 
                    -99.   , -99.   , -99.   , -99.   , 0.85   , & 
                    0.80   , 0.50   , 0.60   , 0.00   , 0.90   , & 
                    0.70   / 
      DATA CVDAT/  &   
                    2.0E-5 , 2.0E-5 , 2.0E-5 , 1.0E-5 , 1.0E-5 , & 
                    1.0E-5 , 1.0E-5 , 1.0E-5 , 1.0E-5 , 2.0E-5 , & 
                    2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
                    2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
                    2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
                    1.5E-5 / 
      DATA RGLDAT/  &  
                    100.   , 100.   , 100.   , 30.    , 30.    , & 
                    30.    , 30.    , 30.    , 30.    , 100.   , & 
                    100.   , 100.   , 100.   , 100.   , 100.   , & 
                    100.   , 100.   , 100.   , 100.   , 100.   , & 
                    100.   , 100.   , 100.   , 100.   , 100.   , & 
                    100.   / 
      DATA GAMMADAT/ & 
                     0.    , 0.     , 0.     , 0.04   , 0.04   , & 
                     0.04  , 0.04   , 0.04   , 0.04   , 0.     , & 
                     0.    , 0.     , 0.     , 0.     , 0.     , & 
                     0.    , 0.     , 0.     , 0.     , 0.     , & 
                     0.04  , 0.     , 0.     , 0.     , 0.     , & 
                     0.04  / 
!
      DATA Z0MDAT / &
                    0.001  , 0.001  , 0.001  , 1.75   , 2.0    , &
                    1.0    , 2.0    , 3.0    , 0.8    , 0.1    , &
                    0.2    , 0.2    , 0.1    , 0.1    , 0.15   , &
                    0.15   , 0.35   , 0.25   , 0.10   , 0.25   , &
                    5.0    , 0.1    , 0.1    , 0.1    , 1.75   , &
                    0.5    / 
!
      DATA EMISDAT/ & 
                    0.991  , 1.000  , 0.991  , 0.996  , 0.996  , & 
                    0.990  , 0.990  , 0.996  , 0.990  , 0.954  , & 
                    0.954  , 0.954  , 0.993  , 0.993  , 0.981  , &
                    0.981  , 0.981  , 0.981  , 0.981  , 0.981  , &
                    1.000  , 0.992  , 0.995  , 0.941  , 0.993  , & 
                    0.993  /
!
      DATA GEXPDAT/ & 
                    2.  , 2.  , 2.  , 2.  , 2.  , & 
                    2.  , 2.  , 2.  , 2.  , 2.  , & 
                    2.  , 2.  , 2.  , 2.  , 2.  , &
                    2.  , 2.  , 2.  , 2.  , 2.  , &
                    2.  , 2.  , 2.  , 2.  , 2.  , & 
                    2.  /
!
!

   !********************************************************************
   !                tables describing the annual evolution of veg fields
   !********************************************************************

   real, save :: d50veg15(13),d95veg15(13),d50veg16(13),d95veg16(13)
   
   data d50veg15/ &
        0.15   , 0.15   , 0.15  , 0.15  , 0.15   , &
        0.15   , 0.30   , 0.30  , 0.30  , 0.15   , &
        0.15   , 0.15   , 0.15                      /
                    
   data  d95veg15/ &
        1.0    , 1.0    , 1.0   , 1.0   , 1.0    , &
        1.0    , 2.5    , 2.5   , 2.5   , 1.0    , &
        1.0    , 1.0    , 1.0                       /
             
   data d50veg16/ &
        0.15   , 0.15   , 0.15  , 0.15  , 0.15   , &
        0.15   , 0.30   , 0.30  , 0.30  , 0.15   , &
        0.15   , 0.15   , 0.15                      /
                    
   data  d95veg16/ &
        1.0    , 1.0    , 1.0   , 1.0   , 1.0    , &
        1.0    , 2.5    , 2.5   , 2.5   , 1.0    , &
        1.0    , 1.0    , 1.0                       /

   real, save :: vegcrops(13),vegdat14(13),vegdat16(13),vegdat17(13)

   data vegcrops/ &
        0.05   , 0.05   , 0.05  , 0.10  , 0.20   , &
        0.40   , 0.60   , 0.6   , 0.6   , 0.05   , &
        0.05   , 0.05   , 0.05                      /
        
   !Meme chose que vegcrops par default
   data vegdat14/ &
        0.05   , 0.05   , 0.05  , 0.10  , 0.20   , &
        0.40   , 0.60   , 0.6   , 0.6   , 0.05   , &
        0.05   , 0.05   , 0.05                      /

   !Meme chose que vegcrops par default
   data vegdat16/ &
        0.05   , 0.05   , 0.05  , 0.10  , 0.20   , &
        0.40   , 0.60   , 0.6   , 0.6   , 0.05   , &
        0.05   , 0.05   , 0.05                      /

   !Meme chose que vegcrops par default
   data vegdat17/ &
        0.05   , 0.05   , 0.05  , 0.10  , 0.20   , &
        0.40   , 0.60   , 0.6   , 0.6   , 0.05   , &
        0.05   , 0.05   , 0.05                      /


   real, save :: lai6(13), lai7(13), lai11(13), lai14(13), lai15(13), &
        lai16(13), lai17(13), lai18(13), lai19(13), lai22(13), &
        lai25(13), lai26(13)

   DATA LAI6 / & 
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   ,  &
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   ,  & 
        1.0   , 0.1   , 0.1                      /
   DATA LAI7 /  &
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   ,  & 
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   ,  &
        1.0   , 0.1   , 0.1                      /
   DATA LAI11/  &
        0.5   , 0.5   , 1.0   , 1.0   , 1.5   ,  & 
        2.0   , 3.0   , 3.0   , 2.0   , 1.5   ,  & 
        1.0   , 0.5   , 0.5                      /
   DATA LAI14/    &
        0.5   , 0.5   , 0.5   , 0.5   , 0.5   ,  &
        0.5   , 1.0   , 2.0   , 2.0   , 1.5   ,  &
        1.0   , 1.0   , 0.5                      /
   data lai15/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   ,  &
        0.1   , 0.1   , 0.1                      /
   DATA LAI16/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  &
        2.5   , 4.0   , 5.0   , 6.0   , 0.1   ,  &
        0.1   , 0.1   , 0.1                      /
   DATA LAI17/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  & 
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   ,  & 
        0.1   , 0.1   , 0.1                      /
   DATA LAI18/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   ,  &
        0.1   , 0.1   , 0.1                      /
   DATA LAI19/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  & 
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   ,  & 
        0.1   , 0.1   , 0.1                      /
   DATA LAI22/  &
        1.0   , 1.0   , 0.5   , 0.1   , 0.1   ,  & 
        0.1   , 0.1   , 1.0   , 2.0   , 1.5   ,  & 
        1.5   , 1.0   , 1.0                      /
   DATA LAI25/  &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   ,  &
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   ,  & 
        3.0   , 3.0   , 3.0                      /
   DATA LAI26/  &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   ,  & 
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   ,  & 
        3.0   , 3.0   , 3.0                      /

   !********************************************************************
   integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: use value from my_jdate_mod

   integer :: i,k
   integer(INT64) :: delti64
   real :: julien, juliens

   real, dimension(nclass) :: laidatdn, laidatds, logz0mloc
   real, dimension(nclass) :: vegdatdn, vegdatds
   
   real, pointer, dimension (:) :: zdlat,zlaideci, zz0mvh,zz0mvl
   real, pointer, dimension (:,:) :: zlaivf26, zvegf
   !--- numeric min for z0m
   real numin_z0m
   data numin_z0m /1.E-4/

! This macro only gets pointer address to pass as an argument to calling function --- to MANIPULATE VARIABLE, USE MACROS  MK... below
#define PTR1D(NAME2) pvars(vd%NAME2%idxv)%data(:)
! Assign local variables below
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni) => pvars(vd%NAME2%idxv)%data(:)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => pvars(vd%NAME2%idxv)%data(:)

      MKPTR1D(zdlat,dlat)
      MKPTR1D(zlaideci,laideci)
      MKPTR1D(zz0mvh,z0mvh)
      MKPTR1D(zz0mvl,z0mvl)

      MKPTR2D(zlaivf26,laivf26)
      MKPTR2D(zvegf,vegf)

!    Default look-up for vf_type .eq. CCILC_WE:
!    This have being tweaked to work with geophy file produced
!    using CCILC_WE and tested just for North America (HRDPS)

      if (vf_type .eq. 'CCILC_WE') &
      call lookup4ccilc_we(ALDAT, D2DAT, D50DAT, D95DAT, VEGDAT, Z0MDAT, GEXPDAT, &
                           d50veg15,d95veg15,d50veg16,d95veg16,vegcrops,vegdat14,vegdat16,vegdat17, &
                           lai11, lai14, lai15, lai16, lai17, nclass)

      ! Read some of the look up tables from namelist
      if (svs_read_aldat)    aldat    = svs_aldat
      if (svs_read_d2dat)    d2dat    = svs_d2dat
      if (svs_read_d50dat)   d50dat   = svs_d50dat
      if (svs_read_d95dat)   d95dat   = svs_d95dat
      if (svs_read_vegdat)   vegdat   = svs_vegdat
      if (svs_read_z0mdat)   z0mdat   = svs_z0mdat
      if (svs_read_gexpveg)  gexpdat   = svs_gexpveg
      
      ! Read look-up monthly climatology
      
      if (svs_read_d50veg15) d50veg15 = svs_d50veg15
      if (svs_read_d50veg16) d50veg16 = svs_d50veg16
      if (svs_read_d95veg15) d95veg15 = svs_d95veg15
      if (svs_read_d95veg16) d95veg16 = svs_d95veg16
      
      if (svs_read_vegcrops) vegcrops = svs_vegcrops
      if (svs_read_vegdat14) vegdat14 = svs_vegdat14
      if (svs_read_vegdat16) vegdat16 = svs_vegdat16
      if (svs_read_vegdat17) vegdat17 = svs_vegdat17
      
      if (svs_read_lai11) lai11    = svs_lai11
      if (svs_read_lai14) lai14    = svs_lai14
      if (svs_read_lai15) lai15    = svs_lai15
      if (svs_read_lai16) lai16    = svs_lai16
      if (svs_read_lai17) lai17    = svs_lai17

      if (svs_urban_params) then ! modify urban surface parameters set above
         cvdat(21)   = 0.3E-5
         z0mdat(21)  = 1.0
         emisdat(21) = 0.950
      endif
      
      ! Initialize arrays that differ for N and S hemispheres
      do i=1,nclass
         laidatdn(i)  = laidat(i)
         laidatds(i)  = laidat(i)
         vegdatdn(i)  = vegdat(i)
         vegdatds(i)  = vegdat(i)
       ! for z0m combine the logs 
         logz0mloc(i) = log(z0mdat(i))
      end do

     ! Determine the current julian day
      delti64 = int (delti64)
      julien = real(jdate_day_of_year(jdateo + kount*delti64 + MU_JDATE_HALFDAY))
      
      ! Use a monthly climatology of rooting depth for veg15 and veg16
      ! Only applies to class 15 et 16 and the same values are used in both hemispheres!!!
      ! TODO: generalize to other crop classes and both hemispheres.
      if (svs_read_d50veg15 .or. vf_type .eq. 'CCILC_WE') then
         d50dat(15) = interpveg(julien, d50veg15)
      endif

      if (svs_read_d50veg16 .or. vf_type .eq. 'CCILC_WE') then
         d50dat(16) = interpveg(julien, d50veg16)
      endif

      if (svs_read_d95veg15 .or. vf_type .eq. 'CCILC_WE') then
         d95dat(15) = interpveg(julien, d95veg15)
      endif

      if (svs_read_d95veg16 .or. vf_type .eq. 'CCILC_WE') then
         d95dat(16) = interpveg(julien, d95veg16)
      endif


      ! Fill the laidatd and vegdatd fields for
      ! land use classes varying with seasons
      ! (i.e., replace the -99 values in the table
      ! with temporal interpolations from the tables above)

      ! tables for northern hemisphere
      laidatdn( 6)  = interpveg(julien , lai6 )
      laidatdn( 7)  = interpveg(julien , lai7 )
      laidatdn(11)  = interpveg(julien , lai11)
      laidatdn(14)  = interpveg(julien , lai14)
      laidatdn(15)  = interpveg(julien , lai15)
      laidatdn(16)  = interpveg(julien , lai16)
      laidatdn(17)  = interpveg(julien , lai17)
      laidatdn(18)  = interpveg(julien , lai18)
      laidatdn(19)  = interpveg(julien , lai19)
      laidatdn(22)  = interpveg(julien , lai22)
      laidatdn(25)  = interpveg(julien , lai25)
      laidatdn(26)  = interpveg(julien , lai26)

      if (svs_read_vegdat14 .or. vf_type .eq. 'CCILC_WE') &
               vegdatdn(14)  = interpveg(julien , vegdat14)
      vegdatdn(15)  = interpveg(julien , vegcrops)
      vegdatdn(16)  = interpveg(julien , vegcrops)
      if (svs_read_vegdat16 .or. vf_type .eq. 'CCILC_WE') &
               vegdatdn(16)  = interpveg(julien , vegdat16)
      vegdatdn(17)  = interpveg(julien , vegcrops)
      if (svs_read_vegdat17 .or. vf_type .eq. 'CCILC_WE') &
               vegdatdn(17)  = interpveg(julien , vegdat17)
      vegdatdn(18)  = interpveg(julien , vegcrops)
      vegdatdn(19)  = interpveg(julien , vegcrops)

      !  tables for southern hermisphere
      juliens = julien  - 183
      if (juliens < 0.) juliens = juliens + 366.

      laidatds( 6)  = interpveg(juliens, lai6 )
      laidatds( 7)  = interpveg(juliens, lai7 )
      laidatds(11)  = interpveg(juliens, lai11)
      laidatds(14)  = interpveg(juliens, lai14)
      laidatds(15)  = interpveg(juliens, lai15)
      laidatds(16)  = interpveg(juliens, lai16)
      laidatds(17)  = interpveg(juliens, lai17)
      laidatds(18)  = interpveg(juliens, lai18)
      laidatds(19)  = interpveg(juliens, lai19)
      laidatds(22)  = interpveg(juliens, lai22)
      laidatds(25)  = interpveg(juliens, lai25)
      laidatds(26)  = interpveg(juliens, lai26)

      if (svs_read_vegdat14 .or. vf_type .eq. 'CCILC_WE') &
               vegdatds(14)  = interpveg(juliens, vegdat14)
      vegdatds(15)  = interpveg(juliens, vegcrops)
      vegdatds(16)  = interpveg(juliens, vegcrops)
      if (svs_read_vegdat16 .or. vf_type .eq. 'CCILC_WE') &
               vegdatds(16)  = interpveg(juliens , vegdat16)
      vegdatds(17)  = interpveg(juliens, vegcrops)
      if (svs_read_vegdat17 .or. vf_type .eq. 'CCILC_WE') &
               vegdatds(17)  = interpveg(juliens , vegdat17)
      vegdatds(18)  = interpveg(juliens, vegcrops)
      vegdatds(19)  = interpveg(juliens, vegcrops)

      do i=1,ni
        do k=1, nclass
           if(zdlat(i).ge.0.0) then
!             northern hemisphere
            zlaivf26(i,k) = zvegf(i,k) * laidatdn(k)
!            zlaictem(i,k) = laidatdn(k)
           else
!             southern hemisphere
            zlaivf26(i,k) = zvegf(i,k) * laidatds(k)
!            zlaictem(i,k) = laidatds(k)
           endif
        enddo
      enddo

      call aggcovernat(PTR1D(vegf), d2dat, d2dat , PTR1D(rootdp), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), d50dat, d50dat , PTR1D(d50), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), d95dat, d95dat , PTR1D(d95), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), maxpdat, maxpdat , PTR1D(maxpond), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), gexpdat, gexpdat , PTR1D(gexp), &
           PTR1D(dlat), ni, nclass)

      do i=1,ni
           if(zdlat(i).ge.0.0) then
           !  northern hemisphere
              zlaideci(i) = laidatdn( 6 )
           else
!             southern hemisphere
              zlaideci(i) = laidatds( 6 )
           endif
        enddo
!
!    Agg Fields for HIGH VEGETATION TYPES
!

      call aggveghigh(PTR1D(vegf), laidatdn, laidatds, PTR1D(laivh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf), aldat, aldat, PTR1D(alvh), &
           PTR1D(dlat), ni, nclass)
      if( .not. read_emis ) &
           call aggveghigh(PTR1D(vegf), emisdat, emisdat, PTR1D(emisvh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf), rsmindat, rsmindat, PTR1D(stomrvh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf), cvdat, cvdat, PTR1D(cvh), &
           PTR1D(dlat), ni, nclass)
!
!         aggregate logs ...
!
      if ( .not. read_z0vh ) then

         call aggveghigh(PTR1D(vegf), logz0mloc, logz0mloc, PTR1D(z0mvh), &
              PTR1D(dlat), ni, nclass)
      
         ! reverse log operation to get final z0h
         DO i=1,ni
            zz0mvh(i)= max( exp(zz0mvh(i)) , numin_z0m )
         ENDDO
         
      endif


      call aggveghigh(PTR1D(vegf), rgldat, rgldat, PTR1D(rglvh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf), gammadat , gammadat, PTR1D(gamvh), &
           PTR1D(dlat), ni, nclass)
!
!    Agg Fields for LOW VEGETATION TYPES
!

      call aggveglow(PTR1D(vegf), laidatdn, laidatds, PTR1D(laivl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf), aldat, aldat, PTR1D(alvl), &
           PTR1D(dlat), ni, nclass)
      if( .not. read_emis ) &
           call aggveglow(PTR1D(vegf), emisdat, emisdat, PTR1D(emisvl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf), rsmindat, rsmindat, PTR1D(stomrvl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf), cvdat, cvdat, PTR1D(cvl), &
           PTR1D(dlat), ni, nclass)
!
!         aggregate logs ...
!
      call aggveglow(PTR1D(vegf), logz0mloc, logz0mloc, PTR1D(z0mvl), &
           PTR1D(dlat), ni, nclass)

      ! reverse log operation to get final z0h
      DO i=1,ni
         zz0mvl(i)= max( exp(zz0mvl(i)) , numin_z0m )
      ENDDO

      call aggveglow(PTR1D(vegf), rgldat, rgldat, PTR1D(rglvl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf), gammadat , gammadat, PTR1D(gamvl), &
           PTR1D(dlat), ni, nclass)
!
!    Compute LOW and HIGH vegetation fractions as well as 
!    DECIDUOUS and EVERGREEN (trees) vegetation fraction
!
      call veglowhigh(PTR1D(vegf), vegdatdn , vegdatds, PTR1D(vegl), &
           PTR1D(vegh), PTR1D(deciduous), PTR1D(evergreen), PTR1D(impervu), &
           PTR1D(dlat), PTR1D(agrifrac), ni, nclass )

   return
 end subroutine inicover_svs

end module inicover_svs_mod
