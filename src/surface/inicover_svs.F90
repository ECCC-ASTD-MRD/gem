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

subroutine inicover_svs(kount, ni, trnch)
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
      REAL RSMINDAT(NCLASS)
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
                     150.   , 150.   ,  100.   , 100.   ,  100.   , & 
                     100.   ,  100.   ,  100.   , 100.   , 150.   , & 
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
                     0.05   , 0.05   , 0.05   , 0.05   , 0.05  , &
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
!

   !********************************************************************
   !                tables describing the annual evolution of veg fields
   !********************************************************************

   real, save :: vegcrops(13)

   data vegcrops/ &
        0.05   , 0.05   , 0.05   , 0.10   , 0.20   , &
        0.40   , 0.60   , 0.6   , 0.6   , 0.05   , &
        0.05   , 0.05   , 0.05                      /


  ! DATA VEGCROPS/  &
  !                  0.05   , 0.05   , 0.95   , 0.95   , 0.95   , &
  !                  0.95   , 0.95   , 0.80   , 0.90   , 0.05   , & 
   !                 0.05   , 0.05   , 0.05                      /


   real, save :: lai6(13), lai7(13), lai11(13), lai14(13), lai15(13), &
        lai16(13), lai17(13), lai18(13), lai19(13), lai22(13), &
        lai25(13), lai26(13)

   DATA LAI6 / & 
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   ,  &
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   , & 
        1.0   , 0.1   , 0.1                      /
   DATA LAI7 /  &
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   , & 
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   ,  &
        1.0   , 0.1   , 0.1                      /
   DATA LAI11/  &
        0.5   , 0.5   , 1.0   , 1.0   , 1.5   , & 
        2.0   , 3.0   , 3.0   , 2.0   , 1.5   , & 
        1.0   , 0.5   , 0.5                      /
   DATA LAI14/    &
        0.5   , 0.5   , 0.5   , 0.5   , 0.5   ,  &
        0.5   , 1.0   , 2.0   , 2.0   , 1.5   ,  &
        1.0   , 1.0   , 0.5                      /
   ! DATA LAI15/  &
   !      0.1   , 0.1   , 0.1   , 0.1   , 0.3   , & 
   !      0.5   , 0.5   , 3.5   , 4.0   , 0.1   , & 
   !      0.1   , 0.1   , 0.1                      /
   data lai15/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   DATA LAI16/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  &
        2.5   , 4.0   , 5.0   , 6.0   , 0.1   ,  &
        0.1   , 0.1   , 0.1                      /
   DATA LAI17/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , & 
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   , & 
        0.1   , 0.1   , 0.1                      /
   DATA LAI18/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   ,  &
        0.1   , 0.1   , 0.1                      /
   DATA LAI19/  &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , & 
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   , & 
        0.1   , 0.1   , 0.1                      /
   DATA LAI22/  &
        1.0   , 1.0   , 0.5   , 0.1   , 0.1   , & 
        0.1   , 0.1   , 1.0   , 2.0   , 1.5   , & 
        1.5   , 1.0   , 1.0                      /
   DATA LAI25/  &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   ,  &
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   , & 
        3.0   , 3.0   , 3.0                      /
   DATA LAI26/  &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   , & 
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   , & 
        3.0   , 3.0   , 3.0                      /

   !********************************************************************
   integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: use value from my_jdate_mod
   real, external :: interpveg

   integer :: i,k
   real :: julien, juliens

   real, dimension(nclass) :: laidatdn, laidatds, logz0mloc
   real, dimension(nclass) :: vegdatdn, vegdatds
   
   real, pointer, dimension (:) :: zdlat,zlaideci, zz0mvh,zz0mvl
   real, pointer, dimension (:,:) :: zlaivf26, zvegf
   !--- numeric min for z0m
   real numin_z0m
   data numin_z0m /1.E-4/

! This macro only gets pointer address to pass as an argument to calling function --- to MANIPULATE VARIABLE, USE MACROS  MK... below
#define PTR1D(NAME2) busptr(vd%NAME2%i)%ptr(1,trnch)
! Assign local variables below
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

      MKPTR1D(zdlat,dlat)
      MKPTR1D(zlaideci,laideci)
      MKPTR1D(zz0mvh,z0mvh)
      MKPTR1D(zz0mvl,z0mvl)

      MKPTR2D(zlaivf26,laivf26)
      MKPTR2D(zvegf,vegf)

      if (svs_urban_params) then ! modify urban surface parameters set above
         cvdat(21)   = 0.3E-5
         z0mdat(21)  = 1.0
         emisdat(21) = 0.950
      endif
      
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

      vegdatdn(15)  = interpveg(julien , vegcrops)
      vegdatdn(16)  = interpveg(julien , vegcrops)
      vegdatdn(17)  = interpveg(julien , vegcrops)
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

      vegdatds(15)  = interpveg(juliens, vegcrops)
      vegdatds(16)  = interpveg(juliens, vegcrops)
      vegdatds(17)  = interpveg(juliens, vegcrops)
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
           PTR1D(dlat), ni, nclass )

!      call veglowhigh(PTR1D(vegf), vegdatdn , vegdatds, bogus1, &
!           bogus2, PTR1D(deciduous), PTR1D(evergreen), PTR1D(impervu), &
!           PTR1D(dlat), ni, nclass )

   return
 end subroutine inicover_svs
