!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END --------------------------

subroutine inicover2(kount, ni, trnch)
   use mu_jdate_mod, only: jdate_day_of_year
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ni, kount, trnch

   !@Author Bernard Bilodeau and Stephane Belair (May 2000)
   !@Revision
   ! 001      see version 5.5.0 for previous history
   !@Object Initialize vegetation fields for the surface schemes
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

   real aldat(nclass), d2dat(nclass), rsminxdat(nclass)
   real laidat(nclass), vegdat(nclass)
   real cvdat(nclass), rgldat(nclass), gammadat(nclass)
 
   data aldat/ &
        0.13   , 0.70   , 0.13   , 0.14   , 0.12   , &
        0.14   , 0.18   , 0.13   , 0.17   , 0.14   , &
        0.18   , 0.19   , 0.20   , 0.19   , 0.20   , &
        0.21   , 0.18   , 0.18   , 0.25   , 0.18   , &
        0.12   , 0.17   , 0.12   , 0.30   , 0.15   , &
        0.15   /

   data d2dat/ &
        1.0    , 1.0    , 1.0    , 3.0    , 3.0    , &
        1.0    , 3.0    , 5.0    , 5.0    , 2.0    , &
        2.0    , 2.0    , 1.5    , 2.0    , 2.0    , &
        1.2    , 1.0    , 1.5    , 2.0    , 1.5    , &
        1.0    , 1.0    , 2.0    , 1.0    , 2.0    , &
        2.0    /

   data rsminxdat/ &
        500.   , 500.   , 500.   , 250.   , 250.   , &
        250.   , 250.   , 250.   , 250.   , 150.   , &
        150.   , 150.   ,  40.   ,  40.   ,  40.   , &
        40.   ,  40.   ,  40.   ,  40.   , 150.   , &
        150.   , 150.   , 150.   , 500.   , 250.   , &
        150.   /
   data laidat/ &
        0.00   , 0.00   , 0.00   , 5.00   , 6.00   , &
        -99.    , -99.   , 6.00   , 4.00   , 3.00   , &
        -99.    , 3.00   , 1.00   , -99.   , -99.   , &
        -99.    , -99.   , -99.   , -99.   , 1.00   , &
        1.00   , -99.   , 4.00   , 0.00   , -99.   , &
        -99.    /
   data vegdat/ &
        0.00   , 0.00   , 0.00   , 0.90   , 0.99   , &
        0.90   , 0.90   , 0.99   , 0.90   , 0.50   , &
        0.50   , 0.50   , 0.85   , 0.30   , -99.   , &
        -99.   , -99.   , -99.   , -99.   , 0.85   , &
        0.10   , 0.50   , 0.60   , 0.00   , 0.90   , &
        0.90   /
   data cvdat/ &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 1.0E-5 , 1.0E-5 , &
        1.0E-5 , 1.0E-5 , 1.0E-5 , 1.0E-5 , 2.0E-5 , &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , &
        2.0E-5 /
   data rgldat/ &
        100.   , 100.   , 100.   , 30.    , 30.    , &
        30.    , 30.    , 30.    , 30.    , 100.   , &
        100.   , 100.   , 100.   , 100.   , 100.   , &
        100.   , 100.   , 100.   , 100.   , 100.   , &
        100.   , 100.   , 100.   , 100.   , 100.   , &
        100.   /
   data gammadat/ &
        0.    , 0.     , 0.     , 0.04   , 0.04   , &
        0.04  , 0.04   , 0.04   , 0.04   , 0.     , &
        0.    , 0.     , 0.     , 0.     , 0.     , &
        0.    , 0.     , 0.     , 0.     , 0.     , &
        0.    , 0.     , 0.     , 0.     , 0.     , &
        0.    /

   !********************************************************************
   !                tables describing the annual evolution of veg fields
   !********************************************************************

   real, save :: vegcrops(13)

   data vegcrops/ &
        0.05   , 0.05   , 0.05   , 0.10   , 0.20   , &
        0.40   , 0.80   , 0.80   , 0.90   , 0.05   , &
        0.05   , 0.05   , 0.05                      /

   real, save :: lai6(13), lai7(13), lai11(13), lai14(13), lai15(13), &
        lai16(13), lai17(13), lai18(13), lai19(13), lai22(13), &
        lai25(13), lai26(13)

   data lai6 / &
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   , &
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   , &
        1.0   , 0.1   , 0.1                      /
   data lai7 / &
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   , &
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   , &
        1.0   , 0.1   , 0.1                      /
   data lai11/ &
        0.5   , 0.5   , 1.0   , 1.0   , 1.5   , &
        2.0   , 3.0   , 3.0   , 2.0   , 1.5   , &
        1.0   , 0.5   , 0.5                      /
   data lai14/ &
        0.5   , 0.5   , 0.5   , 0.5   , 0.5   , &
        0.5   , 1.0   , 2.0   , 2.0   , 1.5   , &
        1.0   , 1.0   , 0.5                      /
   data lai15/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai16/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.5   , 4.0   , 5.0   , 6.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai17/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai18/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai19/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai22/ &
        1.0   , 1.0   , 0.5   , 0.1   , 0.1   , &
        0.1   , 0.1   , 1.0   , 2.0   , 1.5   , &
        1.5   , 1.0   , 1.0                      /
   data lai25/ &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   , &
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   , &
        3.0   , 3.0   , 3.0                      /
   data lai26/ &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   , &
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   , &
        3.0   , 3.0   , 3.0                      /

   !********************************************************************

   integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: use value from my_jdate_mod
   real, external :: interpveg

   integer :: i
   real :: julien, juliens

   real, dimension(nclass) :: aldatd, cvdatd, d2datd, gammadatd, &
        laidatdn, laidatds, rgldatd, rsmindatd, &
        vegdatdn, vegdatds

   IF_ISBA: if (schmsol == 'ISBA') then

      ! Determine the current julian day
      julien = real(jdate_day_of_year(jdateo + kount*int(delt) + MU_JDATE_HALFDAY))

      ! Do the aggregation
      do i=1,nclass
         aldatd(i)    = aldat(i)
         d2datd(i)    = d2dat(i)
         rsmindatd(i) = rsminxdat(i)
         laidatdn(i)  = laidat(i)
         laidatds(i)  = laidat(i)
         vegdatdn(i)  = vegdat(i)
         vegdatds(i)  = vegdat(i)
         cvdatd(i)    = cvdat(i)
         rgldatd(i)   = rgldat(i)
         gammadatd(i) = gammadat(i)
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

#define PTR1D(NAME2) busptr(vd%NAME2%i)%ptr(1,trnch)

      call aggcovernat(PTR1D(vegf), laidatdn, laidatds, PTR1D(lai), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), vegdatdn, vegdatds, PTR1D(vegfrac), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), aldatd, aldatd, PTR1D(alveg), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), d2datd, d2datd , PTR1D(rootdp), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), rsmindatd, rsmindatd, PTR1D(stomr), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), cvdatd, cvdatd, PTR1D(cveg), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), rgldatd, rgldatd, PTR1D(rgl), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), gammadatd , gammadatd, PTR1D(gamveg), &
           PTR1D(dlat), ni, nclass)

   endif IF_ISBA

   return
end subroutine inicover2
