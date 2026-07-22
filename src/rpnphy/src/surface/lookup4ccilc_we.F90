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

module lookup4ccilc_we_mod
  implicit none
  public
contains

subroutine lookup4ccilc_we(ALDAT, D2DAT, D50DAT, D95DAT, VEGDAT, Z0MDAT, GEXPDAT,  &
                           d50veg15,d95veg15,d50veg16,d95veg16,vegcrops,vegdat14,vegdat16,vegdat17, &
                           lai11, lai14, lai15, lai16, lai17, nclass)
   implicit none
!!!#include <arch_specific.hf>

   !@Arguments
   !            - Input -
   ! NCLASS     Number of natural landuse classes
   !            - Output -
   ! All the look-up tables declared in inicovert_svs

   integer nclass
   
   REAL ALDAT(NCLASS), D2DAT(NCLASS), D50DAT(NCLASS), D95DAT(NCLASS)
   REAL VEGDAT(NCLASS), Z0MDAT(NCLASS), GEXPDAT(NCLASS)
   
   real d50veg15(13), d95veg15(13), d50veg16(13), d95veg16(13)
   real vegcrops(13),vegdat14(13),vegdat16(13),vegdat17(13)
   real lai11(13), lai14(13), lai15(13), lai16(13), lai17(13)
        
   !@Author Nicolas Gasset
   !@Revision
   !
   !@Object
   !     Define default look-up tables for vf_type = CCILC_WE
   !     Note there is a redifnition/refinement of four VF!
   !     - VF(13) becomes North American grassland west
   !     - VF(15) becomes North American crops west
   !     - VF(16) becomes North American crops east
   !     - VF(17) becomes North American grassland east
   !     Various look-up table are thus ajusted to be representative 
   !     fr this new definition.
   !     This have being tweaked to work with geophy file produced
   !     using CCILC_WE and tested just for North America (HRDPS).
             
!   
       ALDAT = (/  &
                     0.13   , 0.70   , 0.13   , 0.14   , 0.12   , &
                     0.14   , 0.18   , 0.13   , 0.17   , 0.14   , &
                     0.18   , 0.19   , 0.20   , 0.19   , 0.20   , & 
                     0.20   , 0.20   , 0.18   , 0.25   , 0.18   , & 
                     0.12   , 0.17   , 0.12   , 0.30   , 0.15   , &
                     0.15   /)
!    
       D2DAT = (/  &
                     0.0    , 0.0    , 0.0    , 2.0    , 2.0    , &
                     1.0    , 2.0    , 2.0    , 2.0    , 2.0    , & 
                     2.0    , 2.0    , 1.5    , 2.0    , 2.0    , & 
                     2.0    , 1.5    , 1.5    , 2.0    , 1.5    , & 
                     1.0    , 1.0    , 2.0    , 2.0    , 2.0    , & 
                     2.0    /)
 
!     
       D50DAT = (/  &
                     0.0    , 0.0    , 0.0    , 0.2    , 0.2    , &
                     0.2    , 0.3    , 0.2    , 0.2    , 0.2    , & 
                     0.2    , 0.3    , 0.2    , 0.3    , 0.15   , & 
                     0.15   , 0.3    , 0.2    , 0.2    , 0.2    , & 
                     0.2    , 0.1    , 0.15   , 0.75   , 0.3    , & 
                     0.5    /)
!      
       D95DAT = (/  &
                     0.0    , 0.0    , 0.0    , 0.9    , 0.9    , &
                     0.9    , 2.5    , 1.2    , 0.9    , 0.9    , & 
                     0.9    , 1.5    , 0.9    , 2.5    , 1.0    , & 
                     1.0    , 2.5    , 0.9    , 0.9    , 0.9    , & 
                     0.9    , 0.3    , 0.7    , 2.0    , 2.5    , & 
                     1.5    /)
      GEXPDAT  = (/  &
                    2.  , 2.  , 2.  , 2.  , 2.  , & 
                    2.  , 4.  , 2.  , 2.  , 2.  , & 
                    2.  , 2.  , 2.  , 4.  , 4.  , &
                    4.  , 4.  , 2.  , 2.  , 2.  , &
                    2.  , 2.  , 4.  , 2.  , 4.  , & 
                    2.  /)
!    
!       RSMINDAT = (/  &
!                      500.   , 500.   , 500.   , 250.   , 250.   , &
!                      250.   , 250.   , 250.   , 250.   , 150.   , & 
!                      150.   , 150.   ,  100.  , 100.   ,  100.  , & 
!                      100.   ,  100.  ,  100.  , 100.   , 150.   , & 
!                      150.   , 150.   , 150.   , 500.   , 250.   , & 
!                      250.   /) 
!       LAIDAT = (/  &
!                      0.00   , 0.00   , 0.00   , 5.00   , 6.00   , & 
!                     -99.    , -99.   , 6.00   , 4.00   , 3.00   , & 
!                     -99.    , 3.00   , 1.00   , -99.   , -99.   , &
!                     -99.    , -99.   , -99.   , -99.   , 1.00   , & 
!                      1.00   , -99.   , 4.00   , 0.00   , -99.   , & 
!                     -99.    /) 
!       MAXPDAT = (/  &
!                      0.00   , 0.00   , 0.00   , 0.10   , 0.10   , &
!                      0.10   , 0.10   , 0.10   , 0.10   , 0.05   , &
!                      0.05   , 0.05   , 0.05   , 0.05   , 0.05   , &
!                      0.05   , 0.05   , 0.05   , 0.05   , 0.005  , &
!                      0.005  , 0.05   , 0.05   , 0.05   , 0.10   , &
!                      0.10   /)
!
       VEGDAT = (/  &
                     0.00   , 0.00   , 0.00   , 0.90   , 0.99   , & 
                     0.90   , 0.90   , 0.99   , 0.90   , 0.50   , & 
                     0.40   , 0.50   , 0.70   , -99.   , -99.   , & 
                     -99.   , -99.   , -99.   , -99.   , 0.85   , & 
                     0.80   , 0.50   , 0.60   , 0.00   , 0.90   , & 
                     0.90   /)
!
!       CVDAT = (/  &
!                      2.0E-5 , 2.0E-5 , 2.0E-5 , 1.0E-5 , 1.0E-5 , & 
!                      1.0E-5 , 1.0E-5 , 1.0E-5 , 1.0E-5 , 2.0E-5 , & 
!                      2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
!                      2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
!                      2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , & 
!                      1.5E-5 /) 
!       RGLDAT = (/  &
!                      100.   , 100.   , 100.   , 30.    , 30.    , & 
!                      30.    , 30.    , 30.    , 30.    , 100.   , & 
!                      100.   , 100.   , 100.   , 100.   , 100.   , & 
!                      100.   , 100.   , 100.   , 100.   , 100.   , & 
!                      100.   , 100.   , 100.   , 100.   , 100.   , & 
!                      100.   /) 
!       GAMMADAT = (/  &
!                       0.    , 0.     , 0.     , 0.04   , 0.04   , & 
!                       0.04  , 0.04   , 0.04   , 0.04   , 0.     , & 
!                       0.    , 0.     , 0.     , 0.     , 0.     , & 
!                       0.    , 0.     , 0.     , 0.     , 0.     , & 
!                       0.04  , 0.     , 0.     , 0.     , 0.     , & 
!                       0.04  /)
!
       Z0MDAT = (/  &
                     0.001  , 0.001  , 0.001  , 1.75   , 2.0    , &
                     1.0    , 2.0    , 3.0    , 0.8    , 0.1    , &
                     0.1    , 0.2    , 0.05   , 0.2    , 0.1    , &
                     0.15   , 0.15   , 0.25   , 0.1    , 0.25   , &
                     5.0    , 0.01   , 0.1    , 0.1    , 1.75   , &
                     0.5    /)
!
!        EMISDAT = (/  & 
!                      0.991  , 1.000  , 0.991  , 0.996  , 0.996  , & 
!                      0.990  , 0.990  , 0.996  , 0.990  , 0.954  , & 
!                      0.954  , 0.954  , 0.993  , 0.993  , 0.981  , &
!                      0.981  , 0.993  , 0.981  , 0.981  , 0.981  , &
!                      1.000  , 0.992  , 0.995  , 0.941  , 0.993  , & 
!                      0.993  /)
 
 ! Monthly climatology for the various quatities of lookup tables defined above
   
   d50veg15 = (/  &
         0.15   , 0.15   , 0.15  , 0.15  , 0.15   , &
         0.15   , 0.30   , 0.30  , 0.30  , 0.15   , &
         0.15   , 0.15   , 0.15                   /)
                    
   d95veg15 = (/  &
         1.0    , 1.0    , 1.0   , 1.0   , 1.0    , &
         1.0    , 2.5    , 2.5   , 2.5   , 1.0    , &
         1.0    , 1.0    , 1.0                    /)

   d50veg16 = (/  &
         0.15   , 0.15   , 0.15  , 0.15  , 0.15   , &
         0.15   , 0.30   , 0.30  , 0.30  , 0.15   , &
         0.15   , 0.15   , 0.15                   /)
                    
   d95veg16 = (/  &
         1.0    , 1.0    , 1.0   , 1.0   , 1.0    , &
         1.0    , 2.5    , 2.5   , 2.5   , 1.0    , &
         1.0    , 1.0    , 1.0                    /)

    vegcrops = (/  &
         0.05   , 0.05   , 0.05  , 0.20  , 0.40   , &
         0.60   , 0.80   , 0.80  , 0.60  , 0.40   , &
         0.20   , 0.05   , 0.05                   /)
                
    vegdat14 = (/  &
         0.30   , 0.30   , 0.30  , 0.40  , 0.50   , &
         0.70   , 0.90   , 0.90  , 0.80  , 0.50   , &
         0.40   , 0.30   , 0.30                   /)
 
    vegdat16 = (/  &
         0.05   , 0.05   , 0.05  , 0.20  , 0.40   , &
         0.60   , 0.90   , 0.90  , 0.70  , 0.40   , &
         0.20   , 0.05   , 0.05                   /)
 
    vegdat17 = (/  &
         0.40   , 0.40   , 0.50  , 0.60  , 0.70   , &
         0.80   , 0.90   , 0.90  , 0.80  , 0.70   , &
         0.60   , 0.40   , 0.40                   /)
 

    LAI11 = (/  &
         1.0   , 1.0   , 1.0   , 1.0   , 1.0   , & 
         1.0   , 1.0   , 1.0   , 1.0   , 1.0   , & 
         1.0   , 1.0   , 1.0                      /)

    LAI14 = (/  &
         1.0   , 1.0   , 1.0   , 1.5   , 2.0   ,  &
         3.0   , 4.5   , 4.5   , 3.5   , 2.0   ,  &
         1.0   , 1.0   , 1.0                      /)

    LAI15 = (/  &
         0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
         2.0   , 3.5   , 3.0   , 1.0   , 0.1   , &
         0.1   , 0.1   , 0.1                      /)
         
    LAI16 = (/  &
         0.1   , 0.1   , 0.1   , 0.5   , 1.0   ,  &
         2.0   , 5.0   , 5.0   , 2.5   , 0.5   ,  &
         0.1   , 0.1   , 0.1                      /)
         
    LAI17 = (/  & 
         1.0   , 1.0   , 1.0   , 1.5   , 2.0   , & 
         2.5   , 3.5   , 3.5   , 3.0   , 2.0   , & 
         1.0   , 1.0   , 1.0                      /)

   return
 end subroutine lookup4ccilc_we
end module lookup4ccilc_we_mod
