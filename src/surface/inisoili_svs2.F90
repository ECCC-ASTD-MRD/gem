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

module inisoili_svs2_mod
   implicit none
   private
   
   public :: inisoili_svs2
   
contains
   
subroutine inisoili_svs2(pvars, ni)
   !use tdpack  - M.A. Don't think we need here
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

   integer :: i, k, kk
   REAL b, usb, fb, crit1_wfcint, crit2_wfcint, soillayer_depth, beta_soc
   REAL :: depth_sapric = 1., & ! Depth (m) where the soil properties reach sapric values (Decharm et al. 2016)
           depth_fibric = 0.01 , &! Depth (m) where the soil properties start to depart from fibric values (Decharm et al. 2016)
           rho_sc_max = 130.   ! Maximum soil carbon density (Lawrence and Slater, 2008)
   
   ! "geo" variables are on the levels of the geophysical soil texture datbase
   REAL, dimension(ni,nl_stp) :: wsat_geo, wwilt_geo, wfc_geo, b_geo, psisat_geo, &
           ksat_geo, wfcint_geo, fb_geo ,quartz_geo,rhosoil_geo,conddry_geo,condsld_geo, &
           soilhcapz_dry_geo, fsoc_geo
   ! 100% soil organic content variables 
   REAL, dimension(ni,nl_svs) :: wsat_soc, wwilt_soc, wfc_soc, b_soc, psisat_soc, &
           ksat_soc, wfcint_soc, fb_soc ,conddry_soc,condsld_soc, hcap_soc
   real, pointer, dimension(:) :: zcgsat, zgrkef, zdraindens, zslop

   ! variables on the levels of SVS
   real, pointer, dimension(:,:) :: zbcoef, zclay, zfbcof, zksat, zpsisat, zsand, zwfc, zwfcint, zwsat, zwwilt, zfsoc 
   real, pointer, dimension(:,:) :: zconddry, zcondsld, zquartz, zrhosoil, zsoilhcapz_dry, zgravel, zbulksoil, zoc


  
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni) => pvars(vd%NAME2%idxv)%data(:)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => pvars(vd%NAME2%idxv)%data(:)

   MKPTR1D(zcgsat, cgsat)
   MKPTR1D(zdraindens, draindens)
   MKPTR1D(zgrkef, grkef)
   MKPTR1D(zslop, slop)


   MKPTR2D(zbcoef, bcoef)
   MKPTR2D(zbulksoil , bulksoil)
   MKPTR2D(zclay, clay)
   MKPTR2D(zconddry, conddry)
   MKPTR2D(zcondsld, condsld)  
   MKPTR2D(zfbcof, fbcof)
   MKPTR2D(zfsoc , fsoc)! fraction of soil organic content
   MKPTR2D(zgravel , gravel)
   MKPTR2D(zksat, ksat)
   MKPTR2D(zoc , oc)
   MKPTR2D(zpsisat , psisat)
   MKPTR2D(zquartz, quartz)  
   MKPTR2D(zrhosoil, rhosoil)  
   MKPTR2D(zsand, sand)
   MKPTR2D(zwfc, wfc)
   MKPTR2D(zwfcint, wfcint)
   MKPTR2D(zwsat, wsat)
   MKPTR2D(zwwilt , wwilt)
   MKPTR2D(zsoilhcapz_dry , soilhcapz_dry)


   !call subroutine to compute layer thicknesses
   call layer_thickness()

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
   if (read_oc) then
      do i=1,ni
         do k=1,nl_stp
           ! From GSDE doc: zoc*0.01 is in [% of weight], zbulksoil*0.01 in [g/cm3] and zgravel in [% of volume]
           fsoc_geo(i,k) = zoc(i,k) * 0.01 / 100. * zbulksoil(i,k) * 0.01 *1000. * (1.-zgravel(i,k)/100.) / rho_sc_max
           fsoc_geo(i,k) = max(min(fsoc_geo(i,k), 1.),0.)
         enddo
      enddo
   else
      do i=1,ni
         do k=1,nl_stp
            fsoc_geo(i,k) = 0.
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
         
         ! Compute soil thermal properties for heat diffusion

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

!       Soil heat capacity
        if ((zsand(i,k)+zclay(i,k)) .gt. epsilon_svs) then
            ! From Lawrence and Slater 2008
            soilhcapz_dry_geo(i,k) = (2.128*zsand(i,k) + 2.385*zclay(i,k))/(zsand(i,k)+zclay(i,k))*1.E6
        else    
            ! Previous formulation used in SVS - constant value
            soilhcapz_dry_geo(i,k) = 2700. * 733.
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
            zconddry  (i,k)  = zconddry  (i,k) + conddry_geo  (i,kk)  * weights( k , kk)
            zcondsld  (i,k)  = zcondsld  (i,k) + condsld_geo  (i,kk)  * weights( k , kk)
            zquartz   (i,k)  = zquartz   (i,k) + quartz_geo   (i,kk)  * weights( k , kk)
            zrhosoil  (i,k)  = zrhosoil  (i,k) + rhosoil_geo  (i,kk)  * weights( k , kk)
            zsoilhcapz_dry  (i,k)  = zsoilhcapz_dry  (i,k) + soilhcapz_dry_geo  (i,kk)  * weights( k , kk)
            zfsoc  (i,k)  = zfsoc  (i,k) + fsoc_geo  (i,kk)  * weights( k , kk)
            
         enddo
      enddo
      ! compute thermal coeff. 
      ! for 1st model layer only --- here simply use 1st GEO soil texture !!! Do not map !
      zcgsat (i)  = ( -1.557e-2 * zsand(i,1) &
           -  1.441e-2 * zclay(i,1) + 4.7021 ) * 1.E-6 
      
      ! Compute effective parameter for watdrain
      zgrkef(i)   = 2.* zdraindens(i) * zslop(i)

   enddo

   !     Computer 100% soil organic content properties for typical peat soil profile (Decharme et al. 2016)
   if (read_oc) then
   do i=1,ni
      do k=1,nl_svs

         if (k .EQ. 1) then
            soillayer_depth = delz(k) * 0.5
         else
                soillayer_depth = soillayer_depth + (delz(k) + delz(k-1)) * 0.5
         endif

         beta_soc = log(0.845/0.93)/log(depth_sapric/depth_fibric) ! See Table 1 in Decharme et al. 2016
         wsat_soc  (i,k)  =  0.93 * (soillayer_depth/depth_fibric)**(beta_soc)
         wsat_soc     (i,k)  =  min(max( wsat_soc(i,k), 0.845), 0.93) ! Bound the value to min and max

         beta_soc = log(0.222/0.073)/log(depth_sapric/depth_fibric)
         wwilt_soc (i,k)  =  0.073 * (soillayer_depth/depth_fibric)**(beta_soc)
         wwilt_soc     (i,k)  =  min(max( wwilt_soc(i,k), 0.073), 0.222)

         beta_soc = log(0.719/0.369)/log(depth_sapric/depth_fibric)
         wfc_soc   (i,k)  =  0.369 * (soillayer_depth/depth_fibric)**(beta_soc)
         wfc_soc     (i,k)  =  min(max( wfc_soc(i,k), 0.369), 0.719)

         beta_soc = log(0.0101/0.0103)/log(depth_sapric/depth_fibric)
         psisat_soc(i,k)  =  0.0103 * (soillayer_depth/depth_fibric)**(beta_soc)
         psisat_soc     (i,k)  =  min(max( psisat_soc(i,k), 0.0101), 0.0103)

         beta_soc = log((1.E-7)/(2.8E-4))/log(depth_sapric/depth_fibric)
         ksat_soc  (i,k)  = 2.8E-4 * (soillayer_depth/depth_fibric)**(beta_soc)
         ksat_soc     (i,k)  =  min(max( ksat_soc(i,k), 1.E-7), 2.8E-4)

         beta_soc = log(12./2.7)/log(depth_sapric/depth_fibric)
         b_soc     (i,k)  =  2.7 * (soillayer_depth/depth_fibric)**(beta_soc)
         b_soc     (i,k)  =  min(max( b_soc(i,k), 2.7), 12.)

         hcap_soc(i,k) = 2.5E6

         condsld_soc(i,k) =  0.25 

         conddry_soc(i,k) = 0.05   

      enddo
   enddo
   else
        wsat_soc    (:,:)  = 0.
        wwilt_soc   (:,:)  = 0.
        wfc_soc     (:,:)  = 0.
        psisat_soc  (:,:)  = 0.
        ksat_soc    (:,:)  = 0.
        b_soc       (:,:)  = 0.
        hcap_soc    (:,:)  = 0.
        condsld_soc (:,:)  = 0.
        conddry_soc (:,:)  = 0.
   endif

   ! model soil layers incluting soil organic content
   Do i = 1 , ni
      Do k = 1, nl_svs

         if (read_oc) then
            ! Arithmetic mean
            zwsat  (i,k)  = (1.-zfsoc(i,k)) * zwsat  (i,k) + zfsoc(i,k) * wsat_soc  (i,k)
            zwwilt (i,k)  = (1.-zfsoc(i,k)) * zwwilt  (i,k) + zfsoc(i,k) * wwilt_soc  (i,k)
            zwfc   (i,k)  = (1.-zfsoc(i,k)) * zwfc  (i,k) + zfsoc(i,k) * wfc_soc  (i,k)
            zpsisat(i,k)  = (1.-zfsoc(i,k)) * zpsisat  (i,k) + zfsoc(i,k) * psisat_soc  (i,k)
            zsoilhcapz_dry (i,k) = (1.-zfsoc(i,k)) * zsoilhcapz_dry  (i,k) + zfsoc(i,k) * hcap_soc  (i,k)

            ! Geometric mean
            zksat  (i,k)  = zksat(i,k)**(1.-zfsoc(i,k)) * ksat_soc(i,k)**zfsoc(i,k)
            zconddry  (i,k)  = zconddry(i,k)**(1.-zfsoc(i,k)) * conddry_soc(i,k)**zfsoc(i,k)
            zcondsld  (i,k)  = zcondsld(i,k)**(1.-zfsoc(i,k)) * condsld_soc(i,k)**zfsoc(i,k)

            ! bcoeff and fbcoeff, Arithmetic mean
            b                 =  (1.-zfsoc(i,k)) * zbcoef  (i,k) + zfsoc(i,k) * b_soc  (i,k)
         else
            b                 =  zbcoef  (i,k)
         endif
         
         zbcoef     (i,k)  =  b
         usb               =  1./b
         fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
         zfbcof(i,k)      =  fb
         
         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation
         
         crit1_wfcint   = 2.*zdraindens(i)*zpsisat(i,k)*(zwsat(i,k)/zwwilt(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*zpsisat(i,k)*fb**b
         
         if (abs(zslop(i)).gt.crit1_wfcint) then
            zwfcint(i,k) = zwwilt(i,k)        
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            zwfcint(i,k) = zwsat(i,k) 
         elseif (zslop(i).ne.0.0) then
            zwfcint(i,k) = zwsat(i,k) * fb * &
                 ( zpsisat(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            zwfcint(i,k) = zwfc(i,k)
         endif

      enddo
      
   enddo

   return
 end subroutine inisoili_svs2

end module inisoili_svs2_mod
