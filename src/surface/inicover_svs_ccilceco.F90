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

subroutine inicover_svs_ccilceco(kount, ni, trnch)
   use mu_jdate_mod, only: jdate_day_of_year
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   integer ni, kount, trnch

   !@Author S. Belair et  M. Abrahamowicz (Jan 2016)
   !@Revision
   !@Object Initialize vegetation fields for SVS scheme
   !@Arguments
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
   !     Class       Vegetation type --- ***NEW in 2019 (Maria A.)
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
   !       14        ***TEMPERATE long grass
   !       15        ***TEMPERATE crops
   !       16        ***TROPICAL and SUBTROPICAL long grass
   !       17        ***TROPICAL and SUBTROPICAL crops
   !       18        ***TROPICAL and SUBTROPICAL mixed forest/shrubs
   !       19        ***Evergreen needleleafs in savannas
   !       20        irrigated crops
   !       21        urban
   !       22        tundra
   !       23        swamp
   !       24        desert
   !       25        ***TEMPERATE mixed wood forests
   !       26        ***TEMPERATE mixed shrubs

   !********************************************************************
   ! Tables for the veg characteristics for each veg type
   !********************************************************************
      REAL ALDAT(NCLASS), D2DAT(NCLASS), D50DAT(NCLASS), D95DAT(NCLASS)
      REAL RSMINDAT(NCLASS)
      REAL LAIDAT(NCLASS), VEGDAT(NCLASS),EMISDAT(NCLASS) 
      REAL CVDAT(NCLASS), RGLDAT(NCLASS), GAMMADAT(NCLASS)
      REAL Z0MDAT(NCLASS)
!
      DATA ALDAT/ &
                     0.13   , 0.70   , 0.13   , 0.11   , 0.12   , &
                     0.11   , 0.14   , 0.13   , 0.17   , 0.13   , &
                     0.14   , 0.19   , 0.20   , 0.18   , 0.16   , & 
                     0.17   , 0.17   , 0.16   , 0.11   , 0.18   , & 
                     0.12   , 0.14   , 0.12   , 0.26   , 0.12   , &
                     0.13   / 
!    
      DATA D2DAT/    &
                    0.0    , 0.0    , 0.0    , 2.0    , 2.0    , &
                    1.0    , 2.0    , 2.0    , 2.0    , 2.0    , & 
                    2.0    , 2.0    , 1.5    , 3.0    , 3.0    , & 
                    2.0    , 2.0    , 2.0    , 3.0    , 1.5    , & 
                    1.0    , 1.0    , 2.0    , 2.0    , 2.0    , & 
                    2.0    /
      
!    
      DATA D50DAT/    &
                    0.0    , 0.0    , 0.0    , 0.2    , 0.2    , &
                    0.2    , 0.2    , 0.2    , 0.2    , 0.2    , & 
                    0.2    , 0.3    , 0.5    , 1.0    , 1.0    , & 
                    0.2    , 0.5    , 0.2    , 0.75   , 0.2    , & 
                    0.2    , 0.1    , 0.15   , 0.75   , 0.2    , & 
                    0.5    / 
!
      DATA D95DAT/    &
                    0.0    , 0.0    , 0.0    , 0.9    , 0.9    , &
                    0.9    , 0.9    , 1.2    , 0.9    , 0.9    , & 
                    0.9    , 1.5    , 1.5    , 2.4    , 2.9    , & 
                    1.2    , 1.5    , 0.9    , 2.5    , 0.9    , & 
                    0.9    , 0.3    , 0.7    , 2.0    , 0.9    , & 
                    1.5    / 

      !
      DATA LAIDAT/ &
                     0.00   , 0.00   , 0.00   , 4.00   , 6.00   , & 
                    -99.    , -99.   , 6.00   , -99.   , 2.00   , & 
                    -99.    , 0.75   , 1.00   , -99.   , -99.   , &
                    -99.    , -99.   , 4.00   , 4.00   , 1.00   , & 
                     1.00   , -99.   , 4.00   , 0.00   , -99.   , & 
                    -99.    / 
      DATA VEGDAT/ &
                     0.00   , 0.00   , 0.00   , 0.90   , 0.99   , & 
                     0.90   , 0.90   , 0.99   , 0.50   , 0.50   , & 
                     0.50   , 0.30   , 0.7    , -99.   , -99.   , & 
                     0.60   , -99.   , 0.5    , 0.70   , 0.85   , & 
                     0.80   , 0.50   , 0.70   , 0.00   , 0.90   , & 
                     0.70   / 

! 1 for VH, 2 for VL
      DATA CVDAT/  &   
                     2.0E-5 , 2.0E-5 , 2.0E-5 , 1.0E-5 , 1.0E-5 , & 
                     1.0E-5 , 1.0E-5 , 1.0E-5 , 1.0E-5 , 2.0E-5 , & 
                     2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
                     2.0E-5 , 2.0E-5 , 1.0E-5 , 1.0E-5 , 2.0E-5 , & 
                     2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
                     1.5E-5 / 
!
      DATA Z0MDAT / &
                    0.001  , 0.001  , 0.001  , 1.75   , 2.0    , &
                    1.0    , 2.0    , 3.0    , 0.8    , 0.1    , &
                    0.2    , 0.2    , 0.1    , 0.15   , 0.15   , &
                    0.15   , 0.15   , 1.75   , 1.75   , 0.25   , &
                    5.0    , 0.1    , 0.1    , 0.1    , 1.75   , &
                    0.5    / 
! READ IN --- DO NOT RE-SPECIFY....
      DATA EMISDAT/ & 
                     0.991  , 1.000  , 0.991  , 0.996  , 0.996  , & 
                     0.990  , 0.990  , 0.996  , 0.990  , 0.954  , & 
                     0.954  , 0.954  , 0.993  , 0.993  , 0.981  , &
                     0.981  , 0.981  , 0.981  , 0.996  , 0.981  , &
                     1.000  , 0.992  , 0.995  , 0.941  , 0.993  , & 
                     0.993  /

!NOT USED WITH PHOTO:

    
      DATA RSMINDAT/    &
                     500.   , 500.   , 500.   , 250.   , 250.   , &
                     250.   , 250.   , 250.   , 250.   , 150.   , & 
                     150.   , 150.   ,  100.   , 100.   ,  100.   , & 
                     100.   ,  100.   ,  100.   , 250.   , 150.   , & 
                     150.   , 150.   , 150.   , 500.   , 250.   , & 
                     250.   / 


      DATA RGLDAT/  &  
                     100.   , 100.   , 100.   , 30.    , 30.    , & 
                     30.    , 30.    , 30.    , 30.    , 100.   , & 
                     100.   , 100.   , 100.   , 100.   , 100.   , & 
                     100.   , 100.   , 100.   , 30.   , 100.   , & 
                     100.   , 100.   , 100.   , 100.   , 100.   , & 
                     100.   / 
      DATA GAMMADAT/ & 
                     0.    , 0.     , 0.     , 0.04   , 0.04   , & 
                     0.04  , 0.04   , 0.04   , 0.04   , 0.     , & 
                     0.    , 0.     , 0.     , 0.     , 0.     , & 
                     0.    , 0.     , 0.     , 0.04   , 0.     , & 
                     0.04  , 0.     , 0.     , 0.     , 0.     , & 
                     0.04  / 
!
!

   !********************************************************************
   !                tables describing the annual evolution of veg fields
   !********************************************************************

      real, save :: vegdat14(13),vegdat15(13), &
           vegdat15_NHwintercrops(13), vegdat15_SA(13), &
           vegdat15_australia(13),vegdat16(13),vegdat17(13), &
           vegdat17_SH(13),vegdat14_south49(13)


!Temperate long grass north/south of 49N/49S
   data vegdat14/ &
        0.10   , 0.10   , 0.10  , 0.10  , 0.10   , &
        0.30   , 0.40   , 0.50  , 0.40  , 0.20   , &
        0.15   , 0.15   , 0.10                      /


  !Temperate long grass south/north of 49N/49S
  data vegdat14_south49/ &
        0.90   , 0.90   , 0.90  , 0.90  , 0.90   , &
        0.90   , 0.90   , 0.90  , 0.90  , 0.90   , &
        0.90   , 0.90   , 0.90                      /


   
! Temperate crops
   data vegdat15/ &
        0.10   , 0.10   , 0.10  , 0.10   , 0.20  , &
        0.60   , 0.70   , 0.70   , 0.5   , 0.15  , &
        0.15   , 0.10   , 0.10                      /

   ! Temperate crops with "NH winter" crops -- below 44N
   ! also boosted summer 
   data vegdat15_NHwintercrops/ &
        0.90   , 0.90   , 0.90  , 0.90   , 0.90  , &
        0.90   , 0.90   , 0.90  , 0.90   , 0.90 , &
        0.90   , 0.90   , 0.90                      /
! Temperate crops in SH, in South America patch
   data vegdat15_SA/  &
        0.80   , 0.80   , 0.80  , 0.80   , 0.80  , &
        0.80   , 0.80   , 0.80  , 0.80   , 0.80 , &
        0.80   , 0.80   , 0.80 /
   ! Temperate crops in SH, Australia
   data vegdat15_australia/ &
        0.70   , 0.70   , 0.60  , 0.20   , 0.20  , &
        0.60   , 0.70   , 0.70   , 0.5   , 0.15  , &
        0.15   , 0.70   , 0.70                      /  
   
! Tropical/subtropical grass/ savanna
   data vegdat16/ &
        0.90   , 0.90   , 0.90   , 0.90   ,0.90  , &
        0.90   , 0.90   , 0.90   , 0.90  , 0.90  , &
        0.90   , 0.90   , 0.90                      /

   ! Tropical/subtropical crops
   ! NH
   data vegdat17/ &
        0.9    , 0.9  , 0.9   , 0.7  , 0.5  , &
        0.65   , 0.80   , 0.80   , 0.80  , 0.50  , &
        0.90   , 0.9   , 0.9                      /

   data vegdat17_SH/ &
        0.80   , 0.80   , 0.60   , 0.15  , 0.15  , &
        0.15   , 0.20   , 0.50   , 0.60  , 0.50  , &
        0.30   , 0.60   , 0.80                     /

   real, save :: lai6(13), lai7(13), lai9(13),lai11(13), lai14(13), lai15(13), &
        lai16(13), lai17(13), lai22(13), &
        lai25(13), lai26(13)

   DATA LAI6 / & 
        0.5  , 0.5  , 0.5  , 0.5  , 0.5   ,  &
        2.0  , 4.0  , 3.0  , 1.0  , 0.5   , & 
        0.5  , 0.5  , 0.5                     /
   DATA LAI7 /  &
        0.5   , 0.5   , 0.5   , 1.0   , 3.0   , & 
        4.0   , 5.0   , 5.0   , 3.0   , 1.0   ,  &
        0.5   , 0.5  , 0.5                     /
   DATA LAI9 /  &
        1.0   , 0.5   , 0.5   , 1.0   , 1.0   , & 
        1.0   , 3.0   , 4.0   , 4.0   , 2.0   ,  &
        2.0   , 2.0  ,  1.0                     /
   DATA LAI11/  &
        0.5   , 0.5   , 1.0   , 1.0   , 1.5   , & 
        2.0   , 2.0   , 2.0   , 2.0   , 1.5   , & 
        1.0   , 0.5   , 0.5                      /
   DATA LAI14/    &
        1.0   , 1.0   , 1.0   , 1.0   , 1.0   ,  &
        1.0   , 1.0   , 1.0   , 1.0   , 1.0   ,  &
        1.0   , 1.0   , 1.0                      /
   data lai15/ &
        0.2   , 0.2   , 0.2   , 0.5   , 1.0   , &
        2.0   , 3.0   , 3.5   , 2.0   , 0.5   , &
        0.2   , 0.2   , 0.2                      /
   DATA LAI16/  &
        1.0   , 1.0  , 1.0   , 1.0   , 1.0   ,  &
        1.0    , 1.0   , 1.0   , 1.0   , 1.0   ,  &
        1.0   , 1.0   , 1.0                      /
   DATA LAI17/  &
        1.0   , 1.0   , 1.0   , 1.0   , 1.0   , & 
        1.0   , 1.0   , 1.0   , 1.0   , 1.0   , & 
        1.0   , 1.0   , 1.0                      /
   DATA LAI22/  &
        0.1   , 0.1   , 0.1   , 0.1   , 0.1   , & 
        0.2   , 0.3  ,  0.3   , 0.3   , 0.2   , & 
        0.2   , 0.1   , 0.1                      /
   DATA LAI25/  &
        1.0   , 1.0   , 1.2   , 1.5   , 2.0   ,  &
        5.0   , 5.5   , 5.0   , 2.0   , 1.5   , & 
        1.2   , 1.0   , 1.0                      /
   DATA LAI26/  &
        0.6   , 0.6   , 0.6   , 0.6   , 0.6   , & 
        1.0   , 1.5   , 2.0   , 1.2   , 1.0   , & 
        0.6   , 0.6   , 0.6                      /

   !********************************************************************
   integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: use value from my_jdate_mod
   real, external :: interpveg

   integer :: i,k
   real :: julien, juliens

   real, dimension(nclass) :: laidatdn, laidatds, logz0mloc
   real, dimension(nclass) :: vegdatdn, vegdatds
   real, dimension(ni) :: mask_central_usa, mask_SA_vf15, mask_australia, lat, lon

   real, pointer, dimension (:) :: zdlat,zdlon, zlaideci, zz0mvh,zz0mvl
   real, pointer, dimension (:,:) :: zlaivf26, zvegdati, zvegf, zvegf_evol
   !--- numeric min for z0m
   real numin_z0m
   data numin_z0m /1.E-4/
   real rad2deg
   !convert radians to degrees
   rad2deg=180./acos(-1.)
   
! This macro only gets pointer address to pass as an argument to calling function --- to MANIPULATE VARIABLE, USE MACROS  MK... below
#define PTR1D(NAME2) busptr(vd%NAME2%i)%ptr(1,trnch)
! Assign local variables below
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

      MKPTR1D(zdlat,dlat)
      MKPTR1D(zlaideci,laideci)
      MKPTR1D(zdlon,dlon)
      MKPTR1D(zz0mvh,z0mvh)
      MKPTR1D(zz0mvl,z0mvl)

      MKPTR2D(zlaivf26,laivf26)
      MKPTR2D(zvegdati,vegdati)
      MKPTR2D(zvegf,vegf)
      MKPTR2D(zvegf_evol,vegf_evol)

      
     ! Determine the current julian day
      julien = real(jdate_day_of_year(jdateo + kount*int(delt) + MU_JDATE_HALFDAY))

      do i=1,nclass
         ! initialize arrays that differ for N and S hemispheres
         laidatdn(i)  = laidat(i)
         laidatds(i)  = laidat(i)
         vegdatdn(i)  = vegdat(i)
         vegdatds(i)  = vegdat(i)
       ! for z0m combine the logs 
         logz0mloc(i) = log(z0mdat(i))
      end do

      ! Fill the laidatd and vegdatd fields for
      ! land use classes varying with seasons
      ! (i.e., replace the -99 values in the table
      ! with temporal interpolations from the tables above)

      ! tables for northern hemisphere

      laidatdn( 6)  = interpveg(julien , lai6 )
      laidatdn( 7)  = interpveg(julien , lai7 )
      laidatdn( 9)  = interpveg(julien , lai9 )
      laidatdn(11)  = interpveg(julien , lai11)
      laidatdn(14)  = interpveg(julien , lai14)
      laidatdn(15)  = interpveg(julien , lai15)
      laidatdn(16)  = interpveg(julien , lai16)
      laidatdn(17)  = interpveg(julien , lai17)
      laidatdn(22)  = interpveg(julien , lai22)
      laidatdn(25)  = interpveg(julien , lai25)
      laidatdn(26)  = interpveg(julien , lai26)

      vegdatdn(14)  = interpveg(julien , vegdat14)
      vegdatdn(15)  = interpveg(julien , vegdat15)
      vegdatdn(16)  = interpveg(julien , vegdat16)
      vegdatdn(17)  = interpveg(julien , vegdat17)


      !  tables for southern hermisphere
      juliens = julien  - 183
      if (juliens < 0.) juliens = juliens + 366.

      laidatds( 6)  = interpveg(juliens, lai6 )
      laidatds( 7)  = interpveg(juliens, lai7 )
      laidatds( 9)  = interpveg(juliens, lai9 )
      laidatds(11)  = interpveg(juliens, lai11)
      laidatds(14)  = interpveg(juliens, lai14)
      laidatds(15)  = interpveg(juliens, lai15)
      laidatds(16)  = interpveg(juliens, lai16)
      laidatds(17)  = interpveg(juliens, lai17)      
      laidatds(22)  = interpveg(juliens, lai22)
      laidatds(25)  = interpveg(juliens, lai25)
      laidatds(26)  = interpveg(juliens, lai26)

  
      vegdatds(14)  = interpveg(juliens , vegdat14)
      vegdatds(15)  = interpveg(juliens , vegdat15)
      vegdatds(16)  = interpveg(juliens , vegdat16)
      vegdatds(17)  = interpveg(juliens , vegdat17)

      ! create masks

      do i=1,ni
         lat(i) = rad2deg * zdlat(i)         
         lon(i) = rad2deg * zdlon(i)          
         if (lon(i).gt.180.0) then
            lon(i)=lon(i)-360.
         endif
         
          if(lat(i).ge.36.0.and.lat(i).le.49.0.and.&
              lon(i).ge.-108.5.and.lon(i).le.-88.8) then
             mask_central_usa(i) =1.0
          else
             mask_central_usa(i) =0.0
          endif

          if(lat(i).gt.-40.0.and.lat(i).le.-30..and.&
              lon(i).ge.-70.0.and.lon(i).le.-50.0) then
             mask_SA_vf15(i) =1.0
          else
             mask_SA_vf15(i) =0.0
          endif

          if(lat(i).gt.-45.0.and.lat(i).le.-10..and.&
              lon(i).ge.112.0.and.lon(i).le.154.0) then
             mask_Australia(i) =1.0
          else
             mask_Australia(i) =0.0
          endif
          
       enddo


      ! calculate actual vegetation fraction, to be used in
      ! weighted means below
       ! i.e., want vegf x vegdat

       !first calculate vegdat:

       ! Classes 1-13 -- simple interpolation
       do k=1,13
          do i=1,ni
             if (zvegf(i,k).gt.0.0) then
                if(lat(i).ge.0.0) then
                   zvegdati(i,k) = vegdatdn(k)
                else
                  zvegdati(i,k) = vegdatds(k)
               endif
            else
               zvegdati(i,k) = 0.0
            endif
          enddo
       enddo

       ! Class 14 -- Temperate long grass: Use 2 tables:
       ! Default one, and boosted one in [49N,49S] latitudes
       ! except for central USA
       k=14
       do i=1,ni
          if(zvegf(i,k).gt.0.0) then
             if(abs(lat(i)).le.49.0.and.mask_central_usa(i).eq.0.0) then
                ! use boosted table
                vegdatds(14)  = interpveg(juliens , vegdat14_south49)
                vegdatdn(14)  = interpveg(julien , vegdat14_south49)
             else
                vegdatds(14)  = interpveg(juliens , vegdat14)
                vegdatdn(14)  = interpveg(julien , vegdat14)
             endif
          
             if(lat(i).ge.0.0) then
                zvegdati(i,k) = vegdatdn(k)
             else
                zvegdati(i,k) = vegdatds(k)
             endif
          else
             zvegdati(i,k) = 0.0
          endif
       enddo

       ! Class 15 -- Temperate crops
       ! Several look-up tables
       k=15
       do i=1,ni
          if(zvegf(i,k).gt.0.0) then
             if(lat(i).gt.0.0) then
                !NH
                if(lat(i).le.44.0.and.mask_central_usa(i).eq.0.0) then
                ! Southern latitudes in NH + not central USA-- use winter crop LT
                   zvegdati(i,k) = interpveg(julien, vegdat15_NHwintercrops)
                else
                   zvegdati(i,k) = interpveg(julien, vegdat15)
                endif
             else
                ! SH
                ! patch in southern America that should probably
                ! be reclassified as grassland -- boost otherwise,default
                if(mask_SA_vf15(i).eq.1.0) then
                   zvegdati(i,k) = interpveg(juliens, vegdat15_SA)
                else if (mask_australia(i).eq.1.0) then
                   zvegdati(i,k) = interpveg(juliens, vegdat15_australia)
                else
                   zvegdati(i,k)  = interpveg(juliens , vegdat15)
                endif               
             endif
          else
             zvegdati(i,k) = 0.0
          endif
       enddo


       ! Classes 16-26 -- simple interpolation
       ! except 17, use 2 different LT for NH and SH
       ! swamp 23 -- boost 90 under 49N
       do k=16,26
          if(k.ne.17) then
             do i=1,ni
                if (zvegf(i,k).gt.0.0) then
                   if(lat(i).ge.0.0) then
                      zvegdati(i,k) = vegdatdn(k)
                   else
                      zvegdati(i,k) = vegdatds(k)
                   endif
                else
                   zvegdati(i,k) = 0.0
                endif
             enddo
          else if(k.eq.17) then
             do i=1,ni
               if (zvegf(i,k).gt.0.0) then
                   if(lat(i).ge.0.0) then
                      zvegdati(i,k) = interpveg(julien, vegdat17)
                   else
                      zvegdati(i,k) = interpveg(juliens, vegdat17_SH)
                   endif
                else
                   zvegdati(i,k) = 0.0
                endif
             enddo
          endif
          if(k.eq.23) then
             do i=1,ni
                if(lat(i).le.49.0) zvegdati(i,k)=0.9
             enddo
          endif
       enddo

       ! compute time evolving veg. fraction

       do k=1,nclass
          do i=1,ni
             zvegf_evol(i,k) = zvegf(i,k) * zvegdati(i,k)       
         enddo
      enddo
      
      ! in photo code, expect LAI to be already weighted by vegfrac (see photo code)
      do k=1, nclass
         do i=1,ni
            if(zdlat(i).ge.0.0) then
               !             northern hemisphere
               zlaivf26(i,k) = zvegf_evol(i,k) * laidatdn(k)
            else
               !             southern hemisphere
               zlaivf26(i,k) = zvegf_evol(i,k) * laidatds(k)
           endif
        enddo
      enddo

      call aggcovernat(PTR1D(vegf_evol), d2dat, d2dat , PTR1D(rootdp), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf_evol), d50dat, d50dat , PTR1D(d50), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf_evol), d95dat, d95dat , PTR1D(d95), &
           PTR1D(dlat), ni, nclass)

      do i=1,ni
           if(zdlat(i).ge.0.0) then
           !             northern hemisphere
              zlaideci(i) = laidatdn( 6 )
           else
!             southern hemisphere
              zlaideci(i) = laidatds( 6 )
           endif
        enddo
!
!    Agg Fields for HIGH VEGETATION TYPES
!

      call aggveghigh(PTR1D(vegf_evol), laidatdn, laidatds, PTR1D(laivh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf_evol), aldat, aldat, PTR1D(alvh), &
           PTR1D(dlat), ni, nclass)
      if( .not. read_emis ) &
           call aggveghigh(PTR1D(vegf_evol), emisdat, emisdat, PTR1D(emisvh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf_evol), rsmindat, rsmindat, PTR1D(stomrvh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf_evol), cvdat, cvdat, PTR1D(cvh), &
           PTR1D(dlat), ni, nclass)
!
!         aggregate logs ...
!
      if ( .not. read_z0vh ) then

         call aggveghigh(PTR1D(vegf_evol), logz0mloc, logz0mloc, PTR1D(z0mvh), &
              PTR1D(dlat), ni, nclass)
      
         ! reverse log operation to get final z0h
         DO i=1,ni
            zz0mvh(i)= max( exp(zz0mvh(i)) , numin_z0m )
         ENDDO
         
      endif


      call aggveghigh(PTR1D(vegf_evol), rgldat, rgldat, PTR1D(rglvh), &
           PTR1D(dlat), ni, nclass)
      call aggveghigh(PTR1D(vegf_evol), gammadat , gammadat, PTR1D(gamvh), &
           PTR1D(dlat), ni, nclass)
!
!    Agg Fields for LOW VEGETATION TYPES
!

      call aggveglow(PTR1D(vegf_evol), laidatdn, laidatds, PTR1D(laivl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf_evol), aldat, aldat, PTR1D(alvl), &
           PTR1D(dlat), ni, nclass)
      if( .not. read_emis ) &
           call aggveglow(PTR1D(vegf_evol), emisdat, emisdat, PTR1D(emisvl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf_evol), rsmindat, rsmindat, PTR1D(stomrvl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf_evol), cvdat, cvdat, PTR1D(cvl), &
           PTR1D(dlat), ni, nclass)
!
!         aggregate logs ...
!
      call aggveglow(PTR1D(vegf_evol), logz0mloc, logz0mloc, PTR1D(z0mvl), &
           PTR1D(dlat), ni, nclass)

      ! reverse log operation to get final z0h
      DO i=1,ni
         zz0mvl(i)= max( exp(zz0mvl(i)) , numin_z0m )
      ENDDO

      call aggveglow(PTR1D(vegf_evol), rgldat, rgldat, PTR1D(rglvl), &
           PTR1D(dlat), ni, nclass)
      call aggveglow(PTR1D(vegf_evol), gammadat , gammadat, PTR1D(gamvl), &
           PTR1D(dlat), ni, nclass)
!
!    Compute LOW and HIGH vegetation fractions as well as 
!    DECIDUOUS and EVERGREEN (trees) vegetation fraction

      !    Contrary to all weighted means above, veglowhigh requires the VF fraction and not
      !   VF x vegdat... this calculation is done in the subroutine below.
!
      call veglowhigh_ccilceco(PTR1D(vegf), PTR1D(vegf_evol), PTR1D(vegl), &
           PTR1D(vegh), PTR1D(deciduous), PTR1D(evergreen), PTR1D(impervu), &
           ni, nclass )

!      call veglowhigh(PTR1D(vegf), vegdatdn , vegdatds, bogus1, &
!           bogus2, PTR1D(deciduous), PTR1D(evergreen), PTR1D(impervu), &
!           PTR1D(dlat), ni, nclass )

   return
 end subroutine inicover_svs_ccilceco
