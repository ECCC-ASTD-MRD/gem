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

subroutine inisoili_svs(ni, trnch)
   use sfcbus_mod
   use svs_configs
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>

   integer ni, trnch

   !@Author  Maria Abrahamowicz, Stephane Belair , Vincent Fortin (20xx)
   !@Object  Compute soil properties for given soil texture. Compute these properties on 
   !         native levels of database providing soil texture, and then map properties 
   !         unto SVS levels.
   !@Arguments
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   integer :: i, k, kk
   REAL b, usb, fb, crit1_wfcint, crit2_wfcint
   
   ! "geo" variables are on the levels of the geophysical soil texture datbase
   REAL, dimension(ni,nl_stp) :: wsat_geo, wwilt_geo, wfc_geo, b_geo, psisat_geo, &
           ksat_geo, wfcint_geo, fb_geo
   real, pointer, dimension(:) :: zcgsat, zgrkef, zdraindens, zslop

   ! variables on the levels of SVS
   real, pointer, dimension(:,:) :: zbcoef, zclay, zfbcof, zksat, zpsisat, zsand, zwfc, zwfcint, zwsat, zwwilt 

  
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zcgsat, cgsat)
   MKPTR1D(zdraindens, draindens)
   MKPTR1D(zgrkef, grkef)
   MKPTR1D(zslop, slop)

   MKPTR2D(zbcoef, bcoef)
   MKPTR2D(zclay, clay)
   MKPTR2D(zfbcof, fbcof)
   MKPTR2D(zksat, ksat)
   MKPTR2D(zpsisat , psisat)
   MKPTR2D(zsand, sand)
   MKPTR2D(zwfc, wfc)
   MKPTR2D(zwfcint, wfcint)
   MKPTR2D(zwsat, wsat)
   MKPTR2D(zwwilt , wwilt)
   

   ! calculate soil parameters on native GEO layers, and then map them unto model layers. 
   ! calculate weights to be used in phybusinit.... because here... we are
   ! re-doing the calculation for each row of domain...
   ! but the weights are the same !

   ! Compute averaged values for sand and clay for each layer when SOILGRID is used
   if ( soiltext == "SOILGRIDS" ) then   
      do i=1,ni
         do k=1,(nl_stp-1)
            zsand(i,k) = (zsand(i,k)+zsand(i,k+1))/2.
            zclay(i,k) = (zclay(i,k)+zclay(i,k+1))/2.
         enddo
      enddo
   endif

   !     Computer soil properties for GEO layers
   do i=1,ni
      do k=1,nl_stp
         wsat_geo  (i,k)  =  -0.00126   * zsand(i,k) + 0.489
         wwilt_geo (i,k)  =  37.1342e-3 * sqrt(max(1.,zclay(i,k)))
         wfc_geo   (i,k)  =  89.0467e-3 * max(1.,zclay(i,k))**0.3496
         psisat_geo(i,k)  =  0.01 * ( 10.0**(-0.0131 * zsand(i,k) + 1.88) )
         ksat_geo  (i,k)  =  ( 10.0**(0.0153 * zsand(i,k) - 0.884) ) * 7.0556E-6

         b                 =  0.137 * zclay(i,k)  + 3.501
         b_geo     (i,k)  =  b
         usb               =  1./b
         fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
         fb_geo(i,k)      =  fb
         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation

         crit1_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*(wsat_geo(i,k)/wwilt_geo(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*fb**b

         if (abs(zslop(i)).gt.crit1_wfcint) then
            wfcint_geo(i,k) = wwilt_geo(i,k)        
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            wfcint_geo(i,k) = wsat_geo(i,k) 
         elseif (zslop(i).ne.0.0) then
            wfcint_geo(i,k) = wsat_geo(i,k) * fb * &
                 ( psisat_geo(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            wfcint_geo(i,k) = wfc_geo(i,k)
         endif
         
      enddo
   enddo
   ! "Map" GEO soil properties unto model soil layers
   Do i = 1 , ni
      Do k = 1, nl_svs
         do kk = 1 , nl_stp
            
            zwsat  (i,k)  = zwsat  (i,k) + wsat_geo  (i,kk)  * weights( k , kk)
            zwwilt (i,k)  = zwwilt (i,k) + wwilt_geo (i,kk)  * weights( k , kk)
            
            zwfc   (i,k)  = zwfc   (i,k) + wfc_geo   (i,kk)  * weights( k , kk)
            zbcoef (i,k)  = zbcoef (i,k) + b_geo     (i,kk)  * weights( k , kk)
            zfbcof (i,k)  = zfbcof (i,k) + fb_geo    (i,kk)  * weights( k , kk)
            zpsisat(i,k)  = zpsisat(i,k) + psisat_geo(i,kk)  * weights( k , kk)
            zksat  (i,k)  = zksat  (i,k) + ksat_geo  (i,kk)  * weights( k , kk)
            zwfcint(i,k)  = zwfcint(i,k) + wfcint_geo(i,kk)  * weights( k , kk)
            
         enddo
      enddo
      ! compute thermal coeff. 
      ! for 1st model layer only --- here simply use 1st GEO soil texture !!! Do not map !
      zcgsat (i)  = ( -1.557e-2 * zsand(i,1) &
           -  1.441e-2 * zclay(i,1) + 4.7021 ) * 1.E-6 
      
      ! Compute effective parameter for watdrain
      zgrkef(i)   = 2.* zdraindens(i) * zslop(i)

   enddo

   return
 end subroutine inisoili_svs
