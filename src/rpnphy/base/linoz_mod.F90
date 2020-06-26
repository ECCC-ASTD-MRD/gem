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
!-------------------------------------- LICENCE END ---------------------------

module linoz_mod
   implicit none
   public
   save

!     Linoz coefficients  

   integer, parameter ::  jplat  =18
   integer, parameter ::  jplev  =25
   integer, parameter ::  jptime =12
   integer, parameter ::  jpcoeff=14

   real, dimension(jplat,jplev,jptime,jpcoeff)  :: coeff_table

   real, parameter :: p_linoz   = 0.1000000000000e+05 !pressure level for stratospheric linoz ozone Pa  (100 hPa)
   real, parameter :: hu_linoz  = 0.1000000000000e+02 !water vapour level for Linoz stratospheric ozone ppmv (10 ppmv)
   real, parameter :: ptop_clim = 0.4000000000000e+05 !pressure level to relax Linoz O3 to Fortuin-Kelder O3 climatology Pa (400 hPa)
   real, parameter :: ptop_haloe= 0.5000000000000e+02 !pressure level to replace Linoz O3 with Fortuin-Kelder O3 climatology (using HALO obs) Pa (0.5 hPa)
   real, parameter :: mwt_air   = 0.2897000000000e+02 !mol. wgt. of air     (g/mol)
   real, parameter :: mwt_o3    = 0.4800000000000e+02 !mol. weight of ozone (g/mol)
   real, parameter :: mwt_ch4   = 0.1604000000000e+02 !mol. weight of ch4   (g/mol)
   real, parameter :: mwt_n2o   = 0.4401300000000e+02 !mol. weight of n2o   (g/mol)
   real, parameter :: mwt_f11   = 0.1373700000000e+03 !mol. weight of cfc11 (g/mol)
   real, parameter :: mwt_f12   = 0.1209100000000e+03 !mol. weight of cfc12 (g/mol)
   real, parameter :: avno      = 0.6022000000000e+24 ! avogadro's num        atoms mol-1
   real, parameter :: du_o3     = 0.2687000000000e+20 !1 DU = 2.687 E-19 molecules of O3 per square metre
   real, parameter :: gravity   = 0.9806160000000e+1  !m s-2; gravitational acceleration

   real, parameter :: consth    = 0.1607549380000e+07 ! const for conv hum  1.0e+06 / (18.015/28.96)
   real, parameter :: tau_linoz = 2.*24.*3600.        ! 2 days relaxation time (sec)
   real, parameter :: tau_haloe = 3.*3600.            ! 3 hours relaxation time (sec)
   real, parameter :: QepsO3    = 1.e-10              ! mole /mole
   real, parameter :: Qeps      = 1.e-20              ! mole /mole
   
   logical, parameter :: local_dbg = .true.

end module linoz_mod
