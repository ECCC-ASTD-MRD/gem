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

module inisoili_svs_mod
   implicit none
   private
   
   public :: inisoili_svs
   
contains
   
subroutine inisoili_svs(pvars, ni)
   use tdpack
   use sfcbus_mod
   use svs_configs
   use sfc_options
   use phymem, only: phyvar
   implicit none
!!!#include <arch_specific.hf>

   type(phyvar), pointer, contiguous :: pvars(:)
   integer, intent(in) :: ni

   !@Author  Maria Abrahamowicz, Stephane Belair , Vincent Fortin (20xx)
   !@Object  Compute soil properties for given soil texture. Compute these properties on 
   !         native levels of database providing soil texture, and then map properties 
   !         unto SVS levels.
   !@Arguments
   !             - Input/Ouput -
   ! pvars       list of all phy vars (meta + slab data)
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   integer :: i, k, kk, jj
   REAL b, usb, fb, crit1_wfcint, crit2_wfcint, ts
   
   ! "geo" variables are on the levels of the geophysical soil texture datbase
   REAL, dimension(ni,nl_stp) :: cgsat_geo, wsat_geo, wwilt_geo, wfc_geo, b_geo, psisat_geo, &
           ksat_geo, wfcint_geo, fb_geo, quartz_geo,rhosoil_geo,conddry_geo, condminfac_geo, condsld_geo , wunfrz_geo
   real, pointer, dimension(:) :: zgrkef, zdraindens, zslop, zagrifrac

   ! variables on the levels of SVS
   real, pointer, dimension(:,:) :: zbcoef, zclay, zfbcof, zcgsat, zksat, zksatnat, zpsisat, zsand, zwfc, zwfcint, zwsat, &
                               zwwilt, zconddry, zcondminfac, zcondsld , zquartz, zrhosoil, zwunfrz

   ! SVS multiplying coefficient to adjust ksat in agricultural areas
   real, pointer, dimension(:) :: zkasmod_a

   ! Variables used to compute soil properties using USDA2006 method
   real OM, S, C     ! organic, sand and clay content (by weight)
   real theta_33t    ! initial estimate of wfc
   real theta_1500t  ! initial estimate of wwilt
   real theta_S33t   ! initial estimate of wsat - wfc
   real theta_S33    ! 2nd estimate of wsat - wfc

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni) => pvars(vd%NAME2%idxv)%data(:)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => pvars(vd%NAME2%idxv)%data(:)

   MKPTR1D(zdraindens, draindens)
   MKPTR1D(zgrkef, grkef)
   MKPTR1D(zslop, slop)
   MKPTR1D(zagrifrac, agrifrac)
   MKPTR1D(zkasmod_a, kasmod_a)

   MKPTR2D(zbcoef, bcoef)
   MKPTR2D(zcgsat, cgsat)
   MKPTR2D(zclay, clay)
   MKPTR2D(zfbcof, fbcof)
   MKPTR2D(zksat, ksat)
   MKPTR2D(zksatnat, ksatnat)
   MKPTR2D(zpsisat , psisat)
   MKPTR2D(zsand, sand)
   MKPTR2D(zwfc, wfc)
   MKPTR2D(zwfcint, wfcint)
   MKPTR2D(zwsat, wsat)
   MKPTR2D(zwunfrz, wunfrz)
   MKPTR2D(zwwilt , wwilt)
   MKPTR2D(zconddry , conddry)
   MKPTR2D(zcondminfac , condminfac)
   MKPTR2D(zcondsld , condsld)
   MKPTR2D(zrhosoil , rhosoil)
   MKPTR2D(zquartz , quartz)

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
   if (svs_soiltext2prop == "SURFEXV8") then
      do i=1,ni
         do k=1,nl_stp
            wsat_geo  (i,k)  =  -0.00126   * zsand(i,k) + 0.489
            wwilt_geo (i,k)  =  37.1342e-3 * sqrt(max(1.,zclay(i,k)))
            wfc_geo   (i,k)  =  89.0467e-3 * max(1.,zclay(i,k))**0.3496
         enddo
      enddo
   elseif (svs_soiltext2prop == "USDA2006") then
      do i=1,ni
         do k=1,nl_stp
            C = max(0.,min(100.,zclay(i,k)))/100.
            S = max(0.,min(100.-zclay(i,k),zsand(i,k)))/100.
            ! Organic content is identical everywhere except in deserts where it is set to zero
            ! TODO: read from soil texture database
            if (S .gt. 0.85 .and. C .lt. 0.10) then
               ! Sandy soil: set organic content to zero
               OM = 0.
            else
               ! Ensure that fractional weight of sand, clay and organic content combined is less than one
               OM = min(DEFAULT_ORGANIC_CONTENT/100., 1. - S - C)
            endif
            theta_33t = -0.251 * S + 0.195 * C + 0.011 * OM &
               + 0.006 * S * OM - 0.027 * C * OM + 0.452 * S * C + 0.299
            wfc_geo(i,k) = theta_33t + 1.283 * (theta_33t ** 2) - 0.374 * theta_33t - 0.015
            theta_1500t = -0.024 * S + 0.487 * C + 0.006 * OM &
               + 0.005 * S * OM - 0.013 * C * OM + 0.068 * S * C + 0.031
            wwilt_geo(i,k) = theta_1500t + 0.14 * theta_1500t - 0.02
            theta_S33t = 0.278 * S + 0.034 * C + 0.022 * OM &
               - 0.018 * S * OM - 0.027 * C * OM - 0.584 * S * C + 0.078
            theta_S33 = theta_S33t + 0.636 * theta_S33t - 0.107
            wsat_geo(i,k) = wfc_geo(i,k) + theta_S33 - 0.097 * S + 0.043
            wsat_geo(i,k) = max(CRITWATER,min(1.-CRITWATER,wsat_geo(i,k)))
            wfc_geo(i,k) = max(CRITWATER,min(wsat_geo(i,k),wfc_geo(i,k)))
            wwilt_geo(i,k) = max(CRITWATER,min(wfc_geo(i,k),wwilt_geo(i,k)))
         enddo
      enddo
   endif

   do i=1,ni
      do k=1,nl_stp
         cgsat_geo (i,k)  = ( -1.557e-2 * zsand(i,k) &
                             -  1.441e-2 * zclay(i,k) + 4.7021 ) * 1.E-6 
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
         
       ! Compute soil thermal properties for soil freezing

!       Quartz content (ref : NL95 & PL98)):
        quartz_geo(i,k)  = 0.038 + 0.0095*zsand(i,k)

!       Soil dry density (PL98):
        rhosoil_geo(i,k) = (1.0-wsat_geo(i,k))*2700.

!       Soil solid conductivity:
        if (quartz_geo(i,k).gt.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (2.0**(1.0-quartz_geo(i,k)))
        endif
        if (quartz_geo(i,k).le.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (3.0**(1.0-quartz_geo(i,k)))
        endif

!       Soil dry conductivity:
        conddry_geo(i,k) = (0.135*rhosoil_geo(i,k) + 64.7) / &
                        (2700. - 0.947*rhosoil_geo(i,k))

!       factor used to compute the thermal conductivity of soil mineral (only relevant if soil_cond = TIAN2016)
        condminfac_geo(i,k) = (0.182*zsand(i,k)/100 + 0.00775*zclay(i,k)/100 + 0.0534*(100 - zsand(i,k) - zclay(i,k))/100)

        !Use the physical model from Tian et al. (2016) [https://doi.org/10.1111/ejss.12366]
        if (soil_cond == 'TIAN2016') then    
            condsld_geo(i,k) = 7.7**(zsand(i,k)/100)*1.93**(zclay(i,k)/100)*2.74**((100 - zsand(i,k) - zclay(i,k))/100)
        endif
        
!       Unfrozen residual water content obtained from Niu and Yang (2006)
!       Average value between -10 and -2 deg C
        wunfrz_geo(i,k) = 0. 
        do jj = 0,4
            ts =263.15+ 2.*jj
            wunfrz_geo(i,k)     =  wunfrz_geo(i,k)+ wsat_geo(i,k)*(CHLF*(ts-273.15)/(ts*(-1.0*psisat_geo(i,k))*9.81))**(-1.0*usb)
        enddo
        wunfrz_geo(i,k) =  wunfrz_geo(i,k)/5.

      enddo
   enddo
   ! "Map" GEO soil properties unto model soil layers
   Do i = 1 , ni
      Do k = 1, nl_svs
         do kk = 1 , nl_stp
            
            zwsat  (i,k)  = zwsat  (i,k) + wsat_geo  (i,kk)  * weights( k , kk)
            zwwilt (i,k)  = zwwilt (i,k) + wwilt_geo (i,kk)  * weights( k , kk)
            zcgsat (i,k)  = zcgsat (i,k) + cgsat_geo (i,kk)  * weights( k , kk)
            zwfc   (i,k)  = zwfc   (i,k) + wfc_geo   (i,kk)  * weights( k , kk)
            zbcoef (i,k)  = zbcoef (i,k) + b_geo     (i,kk)  * weights( k , kk)
            zfbcof (i,k)  = zfbcof (i,k) + fb_geo    (i,kk)  * weights( k , kk)
            zpsisat(i,k)  = zpsisat(i,k) + psisat_geo(i,kk)  * weights( k , kk)
            zksat  (i,k)  = zksat  (i,k) + ksat_geo  (i,kk)  * weights( k , kk)
            zwfcint(i,k)  = zwfcint(i,k) + wfcint_geo(i,kk)  * weights( k , kk)
            zconddry  (i,k)  = zconddry  (i,k) + conddry_geo  (i,kk)  * weights( k , kk)
            zcondminfac (i,k)  = zcondminfac (i,k) + condminfac_geo (i,kk)  * weights( k , kk)
            zcondsld  (i,k)  = zcondsld  (i,k) + condsld_geo  (i,kk)  * weights( k , kk)
            zquartz   (i,k)  = zquartz   (i,k) + quartz_geo   (i,kk)  * weights( k , kk)
            zrhosoil  (i,k)  = zrhosoil  (i,k) + rhosoil_geo  (i,kk)  * weights( k , kk)  
            zwunfrz   (i,k)  = zwunfrz   (i,k) + wunfrz_geo   (i,kk)  * weights( k , kk)
            
         enddo

         ! Modify ksat with multiplying coefficient in agricultural areas, to represent the effect
         ! of ploughing, which generally affects soils down to 20cm, so the first 3 soil layers
         ! The KASMOD_A coefficient only applies in agricultural areas so weigthed average is computed
         IF(svs_tdrains_plough) then
           zksatnat(i,k) = zksat(i,k)
           IF(k.le.kplough)THEN
             zksat(i,k) = zksat(i,k) * ( 1.0 - zagrifrac(i) + zkasmod_a(i) * zagrifrac(i) )
           ENDIF
         ENDIF

      enddo
      ! -- vfo001 ---
      ! zcgsat is now computed for each layer but for
      ! retrocompatibility with IC4 version of SVS in NSRPS that
      ! relies on GSDE database we re-compute zcgsat for the first
      ! level from the mapped soil texture. This way results are
      ! neutral for NSRPS IC4. This code has no impact on SVS
      ! configs that use a 5cm depth for the first layer and use
      ! the SOILGRIDS database (e.g. MoSA IC5). After IC5 (assuming
      ! no OPS system based on SVS uses GSDE anymore) this line of
      ! code that re-computes zcgsat(i,1) should be removed.
      ! --- vfo001 ---
      ! compute thermal coeff.
      ! for 1st model layer only --- here simply use 1st GEO soil texture !!! Do not map !
      zcgsat (i,1)  = ( -1.557e-2 * zsand(i,1) &
           -  1.441e-2 * zclay(i,1) + 4.7021 ) * 1.E-6
      ! Compute effective parameter for watdrain
      zgrkef(i)   = 2.* zdraindens(i) * zslop(i)

   enddo

   return
 end subroutine inisoili_svs

end module inisoili_svs_mod
