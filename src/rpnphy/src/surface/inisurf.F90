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

!/@*
subroutine inisurf4(kount, ni, nk, trnch)
   use tdpack_const, only: TRPL, TCDK, RAUW
   use sfc_options
   use sfcbus_mod
   use svs_configs
   implicit none
!!!#include <arch_specific.hf>
#include <rmn/msg.h>
#include <rmnlib_basics.hf>
   !@Object Transfer and initialize geophysical fields for the
   !        surface schemes
   !@Arguments
   !       - Input/Ouput -
   ! f        field for permanent physics variables
   ! fsiz     dimension of f
   ! e        field for entry variables
   ! esiz     dimension of e
   ! ni       horizontal dimension

   integer ni, nk, kount, trnch

   !@Author Stephane Belair (February 1999)
   !@NOTE: This subroutine expects snow depth in cm.
   !       The snow depth is converted in metre (in this s/r)
   !       when the 'entry variables' are transfered to the
   !       permanent variables.
   !@Revisions
   ! 001      M. Mackay   (Sep 2022)  - CSLM added
   !*@/

   include "isbapar.cdk"
   include "sfcinput.cdk"
   !include "cslm.cdk"

   real, parameter :: z0ice = 0.001
   real, parameter :: z0sea = 0.001
   real, parameter :: Z0MLAND_MIN_SVS = 0.1

   real, save :: almin  = 0.50
   real, save :: tauf   = 0.24
   real, save :: tauday = 24.

   real    :: tempsum
   integer :: i, k, nk1

   real, pointer, dimension(:) :: &
        zdrainaf, zemisr, zemistg, zemistgen, zglacier, zglsea, &
        zglsea0, zicedp, ziceline, zlakefr, zlhtg, zmg, zml, zresa, zresagr, zresavg, &
        zresasa, zresasv, zriverfr, zslop, zsnoal, zsnoalen, zsnoagen, zsnodpl, zsnoden, &
        zsnoma, zsnoro, zsnvden, zsnvdp, zsnvma, ztsrad, ztwater, &
        zz0en, zz0mland, zz0mlanden, zz0mvh, zz0mvhen, zz0veg, zz0tveg

   real, pointer, dimension(:) :: &
        zrofinlakaf, zlfxi , zlfxo , zevlak

   
   real, pointer, dimension(:,:) :: &
        zalvis, zclay, zclayen, zsand, zsanden, zsnodp, &
        ztglacier, ztmice, ztmoins, ztsoil, zvegf, zz0, zz0t
  

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   !-------------------------------------------------------------

   !  Nothing to process if no fields were read
   if (phyinread_n == 0) return

   MKPTR1D(zdrainaf,drainaf)
   MKPTR1D(zemisr,emisr)
   MKPTR1D(zemistg,emistg)
   MKPTR1D(zemistgen,emistgen)
   MKPTR1D(zglacier,glacier)
   MKPTR1D(zglsea,glsea)
   MKPTR1D(zglsea0,glsea0)
   MKPTR1D(zicedp,icedp)
   MKPTR1D(ziceline,iceline)
   MKPTR1D(zlakefr,lakefr)
   MKPTR1D(zlhtg,lhtg)
   MKPTR1D(zmg,mg)
   MKPTR1D(zml,ml)
   MKPTR1D(zresa,resa)
   MKPTR1D(zresagr,resagr)
   MKPTR1D(zresavg,resavg)
   MKPTR1D(zresasa,resasa)
   MKPTR1D(zresasv,resasv)
   MKPTR1D(zriverfr, riverfr)
   MKPTR1D(zslop,slop)
   MKPTR1D(zsnoal,snoal)
   MKPTR1D(zsnoalen,snoalen)
   MKPTR1D(zsnoagen,snoagen)
   MKPTR1D(zsnoden,snoden)
   MKPTR1D(zsnodpl,snodpl) 
   MKPTR1D(zsnoma,snoma)
   MKPTR1D(zsnoro,snoro)
   MKPTR1D(zsnvden,snvden)
   MKPTR1D(zsnvdp,snvdp) 
   MKPTR1D(zsnvma,snvma)
   MKPTR1D(ztsrad,tsrad)
   MKPTR1D(ztwater,twater)
   MKPTR1D(zz0en,z0en)
   MKPTR1D(zz0mland,z0mland)
   MKPTR1D(zz0mlanden,z0mlanden)
   MKPTR1D(zz0mvh,z0mvh)
   MKPTR1D(zz0mvhen,z0mvhen)
   MKPTR1D(zz0veg,z0veg)
   MKPTR1D(zz0tveg,z0tveg)

   MKPTR1D(zrofinlakaf,rofinlakaf)
   MKPTR1D(zlfxi,lfxi)
   MKPTR1D(zlfxo,lfxo)
   MKPTR1D(zevlak,evlak)
   

   MKPTR2D(zalvis,alvis)
   MKPTR2D(zclay,clay)
   MKPTR2D(zclayen,clayen)
   MKPTR2D(zsand,sand)
   MKPTR2D(zsanden,sanden)
   MKPTR2D(zsnodp,snodp)
   MKPTR2D(ztglacier,tglacier)
   MKPTR2D(ztmice,tmice)
   MKPTR2D(ztmoins,tmoins)
   MKPTR2D(ztsoil,tsoil)
   MKPTR2D(zvegf,vegf)
   MKPTR2D(zz0,z0)
   MKPTR2D(zz0t,z0t)



   ! Several treatments on geophysical fields valid for isba
   ! the water temperature (tm) is decreased for points where the
   ! filtering of mountains lead to an icrease of the water level
   ! (old subroutine modtmtp of gem's dynamic library)

   ! Other consistency tests ...

   if (any('snodp' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do k=1,nsurf
         do i=1,ni
            zsnodp(i,k) = max( 0., zsnodp(i,k))
         end do
      end do
   endif

   if (any('tglacier' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         ztglacier(i,1) = min( trpl, ztglacier(i,1))
         ztglacier(i,2) = min( trpl, ztglacier(i,2))
      end do
   endif

   !========================================================================
   !          for variables common to all surface schemes
   !========================================================================
   !
   !
!VDIR NODEP
   
   if (any('alvis' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         nk1 = size(zalvis,2)
         zalvis(i,indx_soil   ) = zalvis(i,nk1)
         zalvis(i,indx_glacier) = zalvis(i,nk1)
         zalvis(i,indx_water  ) = zalvis(i,nk1)
         zalvis(i,indx_ice    ) = zalvis(i,nk1)
         zalvis(i,indx_agrege ) = zalvis(i,nk1)
         if (schmurb.ne.'NIL') then
            zalvis(i,indx_urb ) = zalvis(i,nk1)
         endif
         if (schmlake.ne.'NIL') then
            zalvis(i,indx_lake ) = zalvis(i,nk1)
         endif
         if (schmriver.ne.'NIL') then
            zalvis(i,indx_river ) = zalvis(i,nk1)
         endif
      enddo
   endif

   if (kount == 0 .and. .not.any('emisr' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         zemisr(i) = 1.
      enddo
   endif

   !       --- snodp deja en metres
   if (any('snodp' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         zsnodp(i,indx_water  ) = 0.0
      enddo
   endif

   if (any('tsoil' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         ztsrad(i) = ztsoil(i,1)
      enddo
   endif

   if (z0veg_only) then
      do i=1,ni
         if (schmsol == 'ISBA')then
            zz0en(i) = zz0veg(i)
         else if (schmsol == 'SVS') then
            zz0en(i) = max(zz0mlanden(i), Z0MLAND_MIN_SVS)
         endif
      enddo
   endif

   if (z0veg_only .or. any('z0en' == phyinread_list_s(1:phyinread_n))) then
      do i=1, ni
         zz0 (i,indx_soil   ) = max(zz0en(i),z0min)
         if (z0veg_only) then
            zz0 (i,indx_glacier) = Z0GLA
         else
            zz0 (i,indx_glacier) = max(zz0en(i),Z0GLA)
         endif
         zz0 (i,indx_water  ) = z0sea
         zz0 (i,indx_ice    ) = z0ice
         zz0 (i,indx_agrege ) = max(zz0en(i),z0min)
         if (schmlake.ne.'NIL') then
            zz0(i,indx_lake ) = z0sea
         endif
         if (schmriver.ne.'NIL') then
            zz0(i,indx_river ) = z0sea
         endif

         zz0t(i,indx_soil   ) = max(zz0en(i),z0min)
         if (z0veg_only)  then
            zz0t(i,indx_glacier) = Z0GLA
         else
            zz0t(i,indx_glacier) = max(zz0en(i),Z0GLA)
         endif
         zz0t(i,indx_water  ) = z0sea
         zz0t(i,indx_ice    ) = z0ice
         zz0t(i,indx_agrege ) = max(zz0en(i),z0min)
         if (schmlake.ne.'NIL') then
            zz0t(i,indx_lake ) = z0sea
         endif
         if (schmriver.ne.'NIL') then
            zz0t(i,indx_river ) = z0sea
         endif
      enddo
   endif
      
   if (any('z0veg' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         zz0veg (i) = max(zz0veg(i),z0min)
         zz0tveg(i) = max(zz0veg(i),z0min)
      enddo
   endif
   
   if (any('glsea0' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         zglsea (i) = zglsea0(i)
      enddo
   endif

   !  Mask for the lakes, lake fraction and river fraction
   !  set lake fraction to zero if no lake scheme 
   !  Also set river fraction to 0.0 in all cases for now
   !  leave mask for lake untouched because used even when no lake
   !  scheme
   if (any('vegf' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         zml(i) = zvegf(i,3)
          if (schmlake.ne.'NIL') then
             zlakefr(i)  = zvegf(i,3)
             zriverfr(i) = 0.0
          else            
             zlakefr(i) = 0.0
             zriverfr(i) = 0.0
          endif
       enddo
    endif

    if (kount == 0 .and. .not.icelac) then
       do i=1,ni
          ziceline(i) = 1.
       enddo
    endif

    if (any('tmice' == phyinread_list_s(1:phyinread_n))) then
       do k=1,nl
         do i=1,ni
            ztmice(i,k) = min(tcdk, ztmice(i,k))
         end do
      end do
   endif

   !========================================================================
   !                             for lakes only
   !========================================================================

   call lacs4(climat, ni, trnch)

   !========================================================================
   !     Special cases

   if (any('icedp' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         !           no snow allowed in the absence of marine ice
         if (zicedp(i).lt.himin) then
            zsnodp(i,indx_ice) = 0.0
         endif
      end do
   endif

   !=========================================================================
   !                                      FOR ISBA ... FOR ISBA ... FOR ISBA
   !=========================================================================

   IF_ISBA: if (schmsol == 'ISBA') then

      if (kount == 0) zresa(1:ni) = 50.

      ! Special operations for the snow variables
      !
      ! Careful here about the units:
      ! "snoro" is the relative density of snow, 
      !         i.e., rho_ice / rho_water (no units)
      ! "snoma" is the snow water equivalent in mm (i.e., kg / m2)
      ! "snoal" is the snow albedo determined from the snow age
      !
      ! Note that "snoag" is in hours ... (tauday also)

!VDIR NODEP
      do i=1,ni
         if (any('snoro' == phyinread_list_s(1:phyinread_n))) then
            zsnoro(i) = max(100.,zsnoro(i)) / rauw
         endif
         if (any('snoro' == phyinread_list_s(1:phyinread_n)) .or. &
              any('snodp' == phyinread_list_s(1:phyinread_n))) then
            zsnoma(i) = rauw * zsnoro(i) * zsnodp(i,indx_soil)
         endif
      end do

      ! For the albedo, there are two possibilities:
      !
      ! 1) if switch "snoalb_anl" is true, then the "i6"
      !    record in the starting standard file (snoalen) contains the snow albedo
      !
      ! 2) if switch "snoalb_anl" is false, then we use the snow age (snoagen)
      !    to derive the snow albedo
      !
      IF_SNO_ALB: if (snoalb_anl) then

         if (any('snoalen' == phyinread_list_s(1:phyinread_n))) then
            do i=1,ni
               zsnoal(i)  =  zsnoalen(i)
            end do
         endif

      else

         ! snow albedo is determined from the snow age according to two different
         ! expressions depending if the snow pack is melting or not

         if (any('snoagen' == phyinread_list_s(1:phyinread_n)) .or. &
              any('snoalen' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
            do i=1,ni
               if (ztmoins(i,nk).lt.trpl) then
                  zsnoal(i)  = ansmax - todry*zsnoagen(i)/tauday
               else
                  zsnoal(i)  = (ansmax-almin) * &
                       exp( -tauf*zsnoagen(i)/tauday ) &
                       + almin
               end if
               zsnoal(i)  = max( zsnoal(i) , almin )
               zsnoal(i)  = min( zsnoal(i) , ansmax )
            end do
         endif

      end if IF_SNO_ALB

      !  Initialize the parameters that depend on vegetation

      if (any('vegf' == phyinread_list_s(1:phyinread_n)) .or. &
           (kntveg > 0 .and. mod(kount,kntveg) == 0)) then
         call inicover2(kount, ni, trnch)
      endif

      ! Sand and clay fractions of the soil are taken as simple averages
      ! of the first 3 layers

!VDIR NODEP
      do i=1,ni
         if (any('sand' == phyinread_list_s(1:phyinread_n))) then
            zsand(i,1) = (zsand(i,1) + zsand(i,2) + zsand(i,3)) / 3.
         endif
         if (any('clay' == phyinread_list_s(1:phyinread_n))) then
            zclay(i,1) = (zclay(i,1) + zclay(i,2) + zclay(i,3)) / 3.
         endif
      end do

      ! Make sure the entry fields are coherent ...

      call coherence3(ni, trnch)

      ! Initialize the soil characteristics using the soil texture

      if (any('clay' == phyinread_list_s(1:phyinread_n)) .or. &
           any('sand' == phyinread_list_s(1:phyinread_n))) then
         call inisoili2(ni, trnch)
      endif

   end if IF_ISBA
!=========================================================================
!                                      FOR SVS  ... FOR SVS  ... FOR SVS 
!=========================================================================
!
!
   IF_SVS: IF (schmsol.EQ.'SVS') THEN
!
!VDIR NODEP
      if (kount == 0) then
         do i=1,ni
            zresagr(i)         = 100.
            zresavg(i)         = 50.
            zresasa(i)         = 100.
            zresasv(i)         = 100.
         end do
      endif
        
      ! snow depth, density + mass calc.
      if (any('snodpl' == phyinread_list_s(1:phyinread_n)) .or. &
           any('snoden' == phyinread_list_s(1:phyinread_n))) then
         do i=1,ni        
            zsnoma(i)  = zsnodpl(i) * zsnoden(i)
         end do
      endif
            
      ! snow depth, density + mass calc.
      if (any('snvdp' == phyinread_list_s(1:phyinread_n)) .or. &
           any('snvden' == phyinread_list_s(1:phyinread_n))) then
         do i=1,ni  
            zsnvma(i)  = zsnvdp(i) * zsnvden(i)
         end do
      endif

      !
      ! For the ALBEDO, for SVS, only one possibility
      !
      ! 1) if switch "snoalb_anl" is true, then read-in 
      ! snow albedo normally (snoal, snval)
      ! 2) OTHERWISE ABORT
      !
      !
      if (.not.snoalb_anl) then

         call msg(MSG_ERROR, '(inisurf) SVS requires snoalb_anl=.TRUE.; snow albedo must be read at entry')
         return

      endif

      if (any('slop' == phyinread_list_s(1:phyinread_n))) then
         ! DDeacu: Ensure that slope is positive and set its minimum value     
         ! max. angle for slope is 45 degrees.       
         do i=1,ni
            if ( zmg(i).gt.critmask ) then
               zslop(i)  = min ( max( abs( zslop(i) ) , 5.e-03 ) , 1.0 ) 
            else
               zslop(i)  = 0.0
            endif
         end do
      endif

      if ( read_emis ) then
        if (any('emistgen' == phyinread_list_s(1:phyinread_n))) then
           do i=1,ni
              zemistg(i)         = zemistgen(i)
           end do
        endif
     endif


     if (svs_dynamic_z0h .or. svs_local_z0m ) then
     
        if (any('z0mlanden' == phyinread_list_s(1:phyinread_n))) then
           ! impose minimum of 0.1 for z0mloc
           do i=1,ni
              zz0mland(i)         = max( zz0mlanden(i) , 0.1)
           enddo
           
           if ( svs_local_z0m ) then
              do i=1,ni
                 zz0 (i,indx_soil   ) = zz0mland(i)
                 ! use z0h=z0m/5. as first approximation
                 zz0t(i,indx_soil   ) = zz0mland(i) * 0.2
              end do
           endif
        endif
     endif
     
     if ( read_z0vh ) then
        if (any('z0mvhen' == phyinread_list_s(1:phyinread_n))) then
            ! impose minimum of 0.1 for z0mvh
           do i=1,ni
              zz0mvh(i)         = max( zz0mvhen(i) , 0.1)
           end do
        endif
     endif

!
!                          Initialize the parameters that depend
!                          on vegetation
!
   
      if (any('vegf' == phyinread_list_s(1:phyinread_n)) .or. &
           (kntveg > 0 .and. mod(kount,kntveg) == 0)) then
         if(vf_type == "CCILCECO") then            
            call inicover_svs_ccilceco(kount, ni, trnch)
         else
            call inicover_svs(kount, ni, trnch)
         endif
      endif
!
!
!
!
!                           Sand and clay fractions 
!
!VDIR NODEP
      soil_data: if ( soiltext == "GSDE" .or. soiltext == "SLC" &
           .or. soiltext == "SOILGRIDS" ) then 

         if (any('sanden' == phyinread_list_s(1:phyinread_n))) then
            do k=1,nl_stp
               do i=1,ni
                  zsand(i,k) = zsanden(i,k)
               end do
            end do
         endif


         if (any('clayen' == phyinread_list_s(1:phyinread_n))) then
            do k=1,nl_stp
               do i=1,ni
                  zclay(i,k) = zclayen(i,k)
               end do
            end do
         endif

         clay_n_sand:if (any('clayen' == phyinread_list_s(1:phyinread_n)) .or. &
              any('sanden' == phyinread_list_s(1:phyinread_n))) then

            do k=1,nl_stp
               do i=1,ni
                  watmask2: if (zmg(i).lt.critmask) then
                     ! OVER WATER...
                     zsand  (i,k)    = 0.0
                     zclay  (i,k)    = 0.0
                  else
                     ! OVER LAND
                     
                     if (zsand(i,k)+zclay(i,k).lt.critexture) then
                        !                If no sand and clay component
                        !                attribute to these points characteristics
                        !                of typical loamy soils
                        zsand(i,k) = 35.
                        zclay(i,k) = 35.
                     else 
                        !                 Minimum of 1% of sand and clay 
                        zsand(i,k) =  max( zsand(i,k) , 1.0) 
                        
                        zclay(i,k) =  max( zclay(i,k) , 1.0)
                        
                        if ( zsand(i,k)+zclay(i,k).gt.100 ) then
                           ! reduce sand & clay  percentage proportionally 
                           tempsum= zsand(i,k) + zclay(i,k)
                           zsand(i,k) = zsand(i,k)/tempsum * 100.
                           zclay(i,k) = zclay(i,k)/tempsum * 100.
                        endif
                     endif
                  endif watmask2              
               enddo
            enddo
            ! initialize soil characteristics 
            call inisoili_svs( ni, trnch )
         endif clay_n_sand
      endif soil_data

      ! Make sure the entry fields are coherent ...
      call coherence3(ni, trnch)
   

     END IF IF_SVS



   !========================================================================
   !                             for TEB only
   !========================================================================

   ! Note that TEB variables do not support reading for kount>0:  phyincread_list_s
   !  would need to be processed within initown() to implement this support.
   if (kount == 0 .and. schmurb == 'TEB') &
        call initown3(ni, trnch)

   return
end subroutine inisurf4
