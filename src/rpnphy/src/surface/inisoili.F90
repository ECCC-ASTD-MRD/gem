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

subroutine inisoili2(ni, trnch)
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>

   integer ni, trnch

   !@Author Stephane Belair (February 1999)
   !@Object Initialize the soil properties from the sand and clay
   !         fraction for 5 layers of the soil
   !@Arguments
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   integer :: i
   real, pointer, dimension(:) :: zacoef, zbcoef, zc1sat, zc2ref, zc3ref, &
        zcgsat, zclay, zpcoef, zsand, zwfc, &
        zwsat, zwwilt

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zacoef, acoef)
   MKPTR1D(zbcoef, bcoef)
   MKPTR1D(zc1sat, c1sat)
   MKPTR1D(zc2ref, c2ref)
   MKPTR1D(zc3ref, c3ref)
   MKPTR1D(zcgsat, cgsat)
   MKPTR1D(zclay,  clay)
   MKPTR1D(zpcoef, pcoef)
   MKPTR1D(zsand,  sand)
   MKPTR1D(zwfc,   wfc)
   MKPTR1D(zwsat,  wsat)
   MKPTR1D(zwwilt, wwilt)

   do i=1,ni
      zwsat  (i)  =  0.001*( -1.08*zsand(i) + 494.305 )
      zwwilt (i)  =  37.1342e-3*sqrt(max(1.,zclay(i)))
      zwfc   (i)  =  89.0467e-3 * max(1.,zclay(i))**0.3496
      zbcoef (i)  =  0.137 * zclay(i)  + 3.501
      zcgsat (i)  =  -1.557e-2 * zsand(i) &
           -  1.441e-2 * zclay(i) + 4.7021
      zcgsat (i)  = 1.e-6 * zcgsat(i)
      zc1sat (i)  = 0.01*( 5.58*zclay(i) + 84.88 )
      zc2ref (i)  = 13.815 * max(1.,zclay(i))**(-0.954)
      zc3ref (i)  = 5.327 * max(1.,zclay(i))**(-1.043)
      zacoef (i)  = 732.42e-3 * max(1.,zclay(i))**(-0.539)
      zpcoef (i)  = 0.134 * zclay(i) + 3.4
   end do

   return
end subroutine inisoili2
