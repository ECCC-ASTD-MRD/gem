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
!**S/P TRACEDATA - TRACE GAS CONCENTRATIONS
!
      BLOCK DATA CCC_TRACEDATA
!
      implicit none
!
!Authors
!
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001
!
!Object
!
!        Input trace gas concentrations in unit ppmv,
!
!Arguments
!
!Implicites
!
#include "ccc_tracegases.cdk"
!
!
!    GHG values provided by CCCma
!
!     DATA   CO2_PPM,  CH4_PPM,  N2O_PPM,   F11_PPM,  F12_PPM,
!    1     / 350.E-6,  1.75E-6,  0.28E-6,   0.280E-9, 0.530E-9 /,
!    2       F113_PPM, F114_PPM
!    3     / 0.050E-9, 0.030E-9 /
!
!    GHG-2006-WMO (note: only the mixing ratios provided below are used in the code)
!
      DATA   CO2_PPM,  CH4_PPM,  N2O_PPM,   F11_PPM,  F12_PPM &
           / 380.E-6,  1.783E-6, 0.3186E-6, 0.280E-9, 0.530E-9 /, &
             F113_PPM, F114_PPM &
           / 0.050E-9, 0.030E-9 /
!
!----------------------------------------------------------------------C
!     INPUT TRACE GAS CONCENTRATIONS IN UNIT PPMV,                     C
!     PARTS PER MILLION BY VOLUME, TRANSFORM TO MASS MIXING RATIO.     C
!     THE SAME AS WATER VAPOR AND OZONE.                               C
!     1.5188126 = 44.    / 28.97                                       C
!     0.5522955 = 16.    / 28.97                                       C
!     1.5188126 = 44.    / 28.97                                       C
!     O2 INPUT AS A CONSTANT, UNIT MIXING RATIO BY MASS                C
!     4.7418019 = 137.37 / 28.97                                       C
!     4.1736279 = 120.91 / 28.97                                       C
!     5.2440456 = 151.92 / 28.97                                       C
!     5.8301691 = 168.90 / 28.97                                       C
!     28.97 MOLECULAR WEIGHT OF AIR, E-06 PER MILLION                  C
!----------------------------------------------------------------------C
!
!
!    GHG values provided by CCCma
!
!
!     DATA   RMCO2,        RMCH4,        RMN2O,         RMO2,
!    1     / 5.3158441E-4, 9.6651712E-7, 4.2526752E-7,  2.315E-1      /,
!    2       RMF11,        RMF12,        RMF113,        RMF114
!    3       1.3277045E-9, 2.2120227E-9, 2.6220228E-10, 1.7490507E-10 /
!
!
!    GHG-2006-WMO
!
      DATA   RMCO2,        RMCH4,        RMN2O,         RMO2 &
           / 5.7714879E-4, 9.8474E-7,    4.8405E-7,     2.315E-1      /, &
             RMF11,        RMF12,        RMF113,        RMF114 &
           / 1.3277045E-9, 2.2120227E-9, 2.6220228E-10, 1.7490507E-10 /
!
      END BLOCK DATA CCC_TRACEDATA

