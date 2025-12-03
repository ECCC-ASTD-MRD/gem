
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
!** S/P SVS
!
subroutine svs2(BUS, BUSSIZ, PTSURF, PTSURFSIZ, DT, KOUNT, TRNCH, N, M, NK)
   use, intrinsic :: iso_fortran_env, only: INT64
   use phy_status, only: phy_error_L, physeterror
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
   use sfcbus_mod
   use sfc_options, only: atm_external, atm_tplus, radslope, jdateo, &
        use_photo, nclass, zu, zt, sl_Lmin_soil, VAMIN, svs_local_z0m, &
        vf_type, nsl, lunique_profile_svs2, lsnow_interception_svs2,  &
        cano_ref_forcing, lwater_ponding_svs,critwater, z0snow_svs2, &
         lsfclayer_crocus_svs2
   use svs_configs
   use suncos, only: suncos2

   use tdpack
   USE MODD_CSTS,     ONLY : XRHOLW

   implicit none
!!!#include <arch_specific.hf>
!
!Authors
!          V. Vionnet, N. Leroux, N. Gauthier, M. Abrahamowicz (2017-2024) 
!Revisions
!
!Object
!         Land surface scheme SVS2 based on the SVS model
!
!Arguments
!
!               - Input/Output -
! BUS           bus of surface variables
!
!               - Input -
! BUSSIZ        size of the surface bus
! PTSURF        surface pointers
! PTSURFSIZ     dimension of ptsurf
! KOUNT         number of timestep
! TRNCH         row number
! DT            timestep
! N             running length
! M             horizontal dimension
! NK            vertical dimension
!
!
!

   integer BUSSIZ, N, NK, KOUNT, TRNCH
   real DT
   real,target :: bus(bussiz)
   integer PTSURFSIZ
   integer PTSURF(PTSURFSIZ)

   integer SURFLEN

! WARNING !!!! x in bus(x(varname,1,1)) is defined in the line below
! it is now case sensitive
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1
! so defined X() also to make it NOT case sensitive
#define X(fptr,fj,fk) x(fptr,fj,fk)

   integer, parameter :: INDX_SFC = INDX_SOIL

   logical, parameter :: TDIAGLIM_FALSE = .false.



!     
! LOCAL ARRAYS defined for variables passed to 
! explicit interface in sl_prelim, sl_sfcmod ... need to pass arrays not address of first
! element, so use:
! bus(x(varname,i,k) :)        instead of 
! bus(x(varname,i,k)  )
! PASSING BUSES WILL NOT WORK FOR EXPLICIT INTERFACE... DIMENSION of VARIABLES
! DEFINED LOCALLY based on size of first variable... which in this case is WHOLE! BUS

   real,pointer,dimension(:) :: hu
   real,pointer,dimension(:) :: ps
   real,pointer,dimension(:) :: tt
   real,pointer,dimension(:) :: uu
   real,pointer,dimension(:) :: vv
   real,pointer,dimension(:) :: z0h
   real,pointer,dimension(:) :: z0m
   real,pointer,dimension(:) :: z0mland   
   real,pointer,dimension(:) :: zdlat
   real,pointer,dimension(:) :: zemisr
   real,pointer,dimension(:) :: zfcor
   real,pointer,dimension(:) :: zqdiag
   real,pointer,dimension(:) :: zqdiagtyp
   real,pointer,dimension(:) :: zqsurf
   real,pointer,dimension(:) :: zsnodp
   real,pointer,dimension(:) :: ztdiag
   real,pointer,dimension(:) :: ztdiagtyp
   real,pointer,dimension(:) :: zthetaa
   real,pointer,dimension(:) :: ztsa
   real,pointer,dimension(:) :: zudiag
   real,pointer,dimension(:) :: zudiagtyp
   real,pointer,dimension(:) :: zvdiag
   real,pointer,dimension(:) :: zvdiagtyp
   real,pointer,dimension(:) :: zzusl
   real,pointer,dimension(:) :: zztsl

   real,pointer,dimension(:) :: zslop
   real,pointer,dimension(:) :: wsatur1
   real,pointer,dimension(:) :: isoil1
   real,pointer,dimension(:) :: wsoil1
   real,pointer,dimension(:) :: zwatpond
   real,pointer,dimension(:) :: zmaxpond
   real,pointer,dimension(:) :: zvegh
   real,pointer,dimension(:) :: zvegl 
!
!

!******************************************************
!     LOCAL ARRAYS  --- ALPHABETICAL ORDER
!******************************************************
!

   integer i,j,m, masklat50(n)

   real,dimension(n) :: alva, cg, cvpa, del_vl, del_vh,  dwaterdt
   real,dimension(n) :: eva, gamva
   real,dimension(n) :: leff, lesnofrac, lesvnofrac, rainrate_mm, rainrate_mm_veg
   real,dimension(n) :: hrsurf, hrsurfgv, leslnofrac, lesvlnofrac
   real,dimension(n) :: rgla, rhoa, snowrate_mm,snowrate_mm_veg, stom_rs, stomra, rpp
   real,dimension(n) :: suncosa, sunother1, sunother2, sunother3
   real,dimension(n) :: sunother4, trad, tva, vdir, vmod, vmod_lmin, wrmax_vl, wrmax_vh, wveglt, wveght
   real,dimension(n) :: wsaturc1
! 
   real, dimension(n,nl_svs) :: isoilt, wsoilt
!
!  for SVS2 only (to reorganise and clean some var)
   real,dimension(n) :: pct, pz0avg_snow,pz0loc_snow, pz0h_snow, pzenith
   real,dimension(n) :: pgfluxsnow
   real,dimension(n) :: pg ! Water flux to the soil column [kg/m2/dt]
   real,dimension(n) :: pforest
   real,dimension(n) :: pgfluxsnow_v,pforest_v,phvegapol_v

   real,dimension(n) :: ptvege  ! Average skin temperature of the vegetation (low and high veg)

   real,dimension(n) :: prg_veg    ! Surface incoming shortwave radiation under high vegetation
   real,dimension(n) :: prat_veg   ! Surface incoming longwave radiation under high vegetation
   real,dimension(n) :: pwind_drift_open ! Wind speed for snowdrift routine in open terrain
   real,dimension(n) :: pwind_top  ! Wind speed at canopy top
   real,dimension(n) :: punload_open ! Unload term in open terrain (set to zero)
   real,dimension(n) :: punload_forest  ! Unload term in forested terrain (computed if snow interception is simulated)
   real,dimension(n) :: puref_veg  ! Forcing height for wind under high vegetation
   real,dimension(n) :: ptref_veg  ! Forcing height for temperature/humidity under high vegetation
   real,dimension(n) :: PZ0HVH  ! Canopy roughness length for heat
   real,dimension(n) :: zz0nat, zz0hnat ! Local variables for grid box average roughness length
   real,dimension(n) :: phm_can ! Heat mass for the high vegetation layer (J K-1 m-2)
   real,dimension(n) :: pscap ! Vegetation layer snow capacities (kg m-2)
   real,dimension(n) :: pfcans ! Canopy layer snowcover fractions from FSM2
   real,dimension(n) :: pres_snca ! Resistance for intercepted snow in high canopy
   real, dimension(n) ::  eg_grid! evaporation rate over bare ground and bare ground below high veg (grid box average) [kg/m2/s]
   real, dimension(n) ::  HVSN_VH !Halstead coefficient of the high vegetation canopy accounting for intercepted snow

     ! NL_SVS VARIABLES
   real, dimension(n,nl_svs) ::  pd_g, pdzg
   real,dimension(n,nl_svs) :: psoil_temp_vgh  ! Soil temperature at the bottom of the snowpack under high vegetation
!

   real, dimension(n,nl_svs) :: wsoiltt
   real, dimension(n,nl_svs) ::  etr_grid ! Evapotranspiration from each layer (grid box average) [m/s]
   real, dimension(n,nl_svs) :: wft, wftv, wftg, wdttv, wdttg
   real, dimension(n,nl_svs) :: delwatgr, delwatvg, delicegr, delicevg





!******************************************************
!
      real,pointer,dimension(:) :: zfsolis
!     
      integer yy, mo, dd, hh, mn, sec
      REAL HZ, HZ0, JULIEN, pond_infilt

      integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200    
!
!     In the offline mode the t-step 0 is (correctly) not performed
      if (atm_external .and. kount == 0) return
!
      SURFLEN = M

! assign pointers
      z0h      (1:n) => bus( x(z0t,1,indx_sfc)   : )
      if ( svs_local_z0m ) then
         z0m      (1:n) => bus( x(z0mland,1,1)       : )
      else
         z0m      (1:n) => bus( x(z0,1,indx_sfc)   : )
      endif
      z0mland  (1:n) => bus( x(z0mland,1,1)      : )
      zdlat    (1:n) => bus( x(dlat,1,1)         : )
      zemisr   (1:n) => bus( x(emisr,1,1)         : )      
      zfcor    (1:n) => bus( x(fcor,1,1)         : )
      zqdiag   (1:n) => bus( x(qdiag,1,1)        : )
      zqdiagtyp(1:n) => bus( x(qdiagtyp,1,indx_sfc) : )
      zqsurf   (1:n) => bus( x(qsurf,1,indx_sfc) : )
      zsnodp   (1:n) => bus( x(snodp,1,indx_sfc) : )
      ztdiag   (1:n) => bus( x(tdiag,1,1)        : )
      ztdiagtyp(1:n) => bus( x(tdiagtyp,1,indx_sfc) : )
      ztsa     (1:n) => bus( x(tsa,1,1)          : )     
      zudiag   (1:n) => bus( x(udiag,1,1)        : )
      zudiagtyp(1:n) => bus( x(udiagtyp,1,indx_sfc) : )
      zvdiag   (1:n) => bus( x(vdiag,1,1)        : )
      zvdiagtyp(1:n) => bus( x(vdiagtyp,1,indx_sfc) : )
      zzusl    (1:n) => bus( x(zusl,1,1)         : )
      zztsl    (1:n) => bus( x(ztsl,1,1)         : )

      wsatur1  (1:n) => bus( x(wsat,1,1)         : )
      isoil1   (1:n) => bus( x(isoil,1,1)        : )
      zwatpond (1:n) => bus( x(watpond,1,1)      : )
      zmaxpond (1:n) => bus( x(maxpond,1,1)      : )
      wsoil1   (1:n) => bus( x(wsoil,1,1)        : )
      zslop    (1:n) => bus( x(slop,1,1)        : )
      zvegh    (1:n) => bus( x(vegh,1,1)        : )
      zvegl    (1:n) => bus( x(vegl,1,1)        : )

      if (atm_tplus) then
         hu       (1:n) => bus( x(huplus,1,nk)      : )
         ps       (1:n) => bus( x(pplus,1,1)        : )
         tt       (1:n) => bus( x(tplus,1,nk)       : )
         zthetaa  (1:n) => bus( x(thetaap,1,1)      : )
         uu       (1:n) => bus( x(uplus,1,nk)       : )
         vv       (1:n) => bus( x(vplus,1,nk)       : )
      else
         hu       (1:n) => bus( x(humoins,1,nk)     : )
         ps       (1:n) => bus( x(pmoins,1,1)       : )
         zthetaa  (1:n) => bus( x(thetaa,1,1)       : )
         tt       (1:n) => bus( x(tmoins,1,nk)      : )
         uu       (1:n) => bus( x(umoins,1,nk)      : )
         vv       (1:n) => bus( x(vmoins,1,nk)      : )
      endif




!  
!
      IF (RADSLOPE) THEN
         zFSOLIS(1:n)   => bus( x(fluslop,1,1)      : )
      ELSE
         zFSOLIS(1:n)   => bus( x(flusolis,1,1)     : )
      ENDIF

     
      ! CONVERT RAINRATE AND SNOWRATE FROM M/S TO MM/S TO MATCH UNITS OF
      ! OTHER WATER FLUXES (EVAPORATION etc.)
      
      DO I=1,N
          rainrate_mm(i) = bus(x(rainrate,i,1)) * M_TO_MM
          snowrate_mm(i) = bus(x(snowrate,i,1)) * M_TO_MM
      ENDDO


      ! Calculate greenwich hour 
      call mu_js2ymdhms(jdateo, yy, mo, dd, hh, mn, sec)
      hz0 = hh + float(mn)/60. + float(sec)/3600.
      hz = amod(hz0+ (float(kount)*dt)/3600., 24.)
      
      !Determine the current julian day
      julien = real(jdate_day_of_year(jdateo + kount*int(dt) + MU_JDATE_HALFDAY))
      !Get local solar angle
      call suncos2(suncosa,sunother1,sunother2,sunother3,sunother4,n, &
                   bus(x(dlat,1,1)),bus(x(dlon,1,1)),hz,julien,.false.)
!

      ! Calculate mask for VF26 , |LAT|<=50 have masklat50=1, otherwise masklat50=0

      DO I=1,N

         if( abs (   bus(x(DLAT,I,1)) * 180./acos(-1.)  ) .le.50.) then
            masklat50(i)=1
         else
            masklat50(i)=0
         endif


      ENDDO

      IF(lwater_ponding_svs .and. kount==1) THEN
          DO I=1,N
!           EG: Adjust max. ponding depth according to bare ground fraction: consider 10mm over bare ground
	    zmaxpond(I) = zmaxpond(I) * (zvegh(I)+zvegl(I)) + 0.01 * (1.-zvegh(I)-zvegl(I))
!	    EG: Adjust max. ponding depth according to slope
            zmaxpond(I) = max(0.0,zmaxpond(I)*(1.0E-10)**zslop(I))
         END DO
      ENDIF


      IF(KOUNT.EQ.1) then

         ! ---------------- Initialize variables for canopy scheme  --------------------
         DO I=1,N 
            bus(x(ESNC,I,1))=0.
            bus(x(ESNCAF,I,1))=0.
         END DO

         ! ---------------- Initialize variables for ES and Crocus snowpack schemes--------------------

         DO I=1,N
            PGFLUXSNOW(I)=0.0
            IF(bus(x(SNOMA_SVS,I,1))>0.) THEN
                bus(x(SNOAL,I,1))=0.8
            ELSE
                bus(x(SNOAL,I,1))=0.1
            ENDIF

            PGFLUXSNOW_V(I)=0.0
            IF(bus(x(SNOMA_SVS,I,1))>0.) THEN
                bus(x(SNVAL,I,1))=0.8

            ELSE
                bus(x(SNVAL,I,1))=0.1
            ENDIF
         END DO
      ENDIF


! ---------------- For Crocus and ES scheme--------------------

     ! Option for roughness lenghts for snow
      DO I=1,N
         IF(Z0SNOW_SVS2=='SURFEX') THEN 
            ! Option used in SURFEX: Use local value for momentum and heat
            ! Note that local value for heat is different in SURFEX and SVS1
            PZ0H_SNOW(I)     = Z0HSNOW_CRO ! For heat  transfer
            PZ0AVG_SNOW(I)   = Z0MSNOW_CRO ! For momentum transfer
            PZ0LOC_SNOW(I)   = Z0MSNOW_CRO !For wind calculation below canopy
         ELSE IF(Z0SNOW_SVS2=='SVS1') THEN
            ! Option used in SVS1: Use grid-averaged value for momentum and local value for heat
            ! Note that local value for heat is different in SURFEX and SVS1
            PZ0H_SNOW(I)     =  Z0HSNOW  ! For heat  transfer
            PZ0AVG_SNOW(I)   =  Z0M(I)  ! For momentum transfer
            PZ0LOC_SNOW(I)   =  Z0HSNOW/Z0M_TO_Z0H  !For wind calculation below canopy
         ELSE IF(Z0SNOW_SVS2=='HYB') THEN
            ! Option used in SVS1: Use grid-averaged value for momentum and local value for heat
            ! Note that local value for heat is different in SURFEX and SVS1
            PZ0H_SNOW(I)     =  Z0HSNOW_CRO  ! For heat  transfer
            PZ0AVG_SNOW(I)   =  Z0M(I)  ! For momentum transfer
            PZ0LOC_SNOW(I)   =  Z0HSNOW_CRO/Z0M_TO_Z0H  !For wind calculation below canopy            
         ENDIF
      ENDDO


      DO I=1,N
             bus(x(RSNOWSA,I,1)) = 0.
             bus(x(RSNOWSV,I,1)) = 0.

             PZENITH(I) =  ACOS(SUNCOSA(I))
             PFOREST(I)=0.
             PFOREST_V(I)=1.
             PHVEGAPOL_V(I) = 0. ! Effect of basal vegetation on snowpack properties are not taken into account in high vegetation. 
             PFCANS(I) = 0.
             PSCAP(I) = 0.
             PRES_SNCA(I) = 0.

            ! TO BE CHECKED======================
            PCT(I)= 1.E-4
            !PCT(I)= 1./(30000.)
            PZ0HVH(I) = Z0M_TO_Z0H * BUS(x(Z0MVH  ,1,1)) ! Z0M_TO_Z0H = 0.2 from svs_configs,


            DO J=1,NL_SVS
               PD_G(I,J)=DL_SVS(J)
               IF(J == 1) THEN
                  PDZG(i,j) = DL_SVS(J)
               ELSE
                  PDZG(i,j) = DL_SVS(J) - DL_SVS(J-1)
               ENDIF
           ENDDO
       ENDDO

! Compute snow diagnostics for some inputs
!
      do I=1,N
!        total snow mass
         bus(x(SNOMA,I,1)) = 0.
         bus(x(SNVMA,I,1)) = 0.
!        total snow depths
         bus(x(SNODPL,I,1))  = 0.
         bus(x(SNVDP ,I,1))  = 0.

         do J=1,NSL
            bus(x(SNOMA ,I,1)) = bus(x(SNOMA ,I,1)) + bus(x(SNOMA_SVS ,I,J))
            bus(x(SNVMA ,I,1)) = bus(x(SNVMA ,I,1)) + bus(x(SNOMAV_SVS,I,J))
            bus(x(SNODPL,I,1)) = bus(x(SNODPL,I,1)) + bus(x(SNOMA_SVS ,I,J))/bus(x(SNODEN_SVS ,I,J))
            bus(x(SNVDP ,I,1)) = bus(x(SNVDP ,I,1)) + bus(x(SNOMAV_SVS,I,J))/bus(x(SNODENV_SVS,I,J))
         enddo
!        mean snow densities (absolute and relative)
         if ( bus(x(SNODPL,I,1)) .gt.0.0) then
            bus(x(SNODEN,I,1)) = bus(x(SNOMA ,I,1))/bus(x(SNODPL, I,1))
            bus(x(SNORO ,I,1)) = bus(x(SNODEN,I,1))/1000.
         else
            bus(x(SNODEN,I,1)) = 50.0
            bus(x(SNORO ,I,1)) = 0.05
         endif
         if ( bus(x(SNVDP,I,1)) .gt.0.0) then
            bus(x(SNVDEN,I,1)) = bus(x(SNVMA ,I,1))/bus(x(SNVDP,I,1))
            bus(x(SNVRO ,I,1)) = bus(x(SNVDEN,I,1))/1000.
         else
            bus(x(SNVDEN,I,1)) = 50.0
            bus(x(SNVRO ,I,1)) = 0.05
         endif
      enddo

!
!******************************************************************
!                  SVS SUBROUTINES START HERE
!******************************************************************
!

!   2 possible approaches for flux and coeff. calculations... IMPOSE minimum wind or minimum Monin?Obukhov Length
!   For minimum Monin?Obukhov, atm wind will be modified internally to insure coupling and desired minimum value
!   In this case, min_wind_speed is really numeric = VAMIN


      if (sl_Lmin_soil > 0.) then
         ! option using minimun Monin?Obukhov Length ( vmod=max(uv,vamin) )
         ! impose minimum wind = VAMIN 
         i = sl_prelim(tt,hu,uu,vv,ps,zzusl,VMOD,VDIR,TVA,RHOA,min_wind_speed=VAMIN)
      else
         ! option using minimum wind speed vmod=max(uv,2.5)
         i = sl_prelim(tt,hu,uu,vv,ps,zzusl,VMOD,VDIR,TVA,RHOA, &
              min_wind_speed=2.5,min_wind_reduc='linear')
      endif


      if (i /= SL_OK) then
         call physeterror('svs', 'error returned by sl_prelim()')
         return
      endif 

!
! EG CODE RELATED TO PONDING OPTION
      IF (lwater_ponding_svs) then
          DO I=1,N
             
             wsaturc1(I)= max((wsatur1(I)-isoil1(I)-0.00001), CRITWATER)
             if (wsoil1(I).lt.wsaturc1(I).and.zwatpond(I).gt.0.0) then 
               ! Amount of ponded water that infiltrates in the soil      
               pond_infilt = min(zwatpond(I),(wsaturc1(I)-wsoil1(I))*dl_svs(1))
               ! Update near-surface soil moisture and amount of ponded water
               wsoil1(I) = wsoil1(I) + pond_infilt / dl_svs(1)
               zwatpond(I) = zwatpond(I) - pond_infilt
             endif
          ENDDO
      ELSE
         DO I=1,N
            zwatpond(I)= 0.0
         ENDDO
      END IF
! END EG CODE RELATED TO PONDING OPTION


      CALL SOILI_SVS2( BUS(x(WSOIL ,1,1)), &
           BUS(x(ISOIL  ,1,1)), &
           BUS(x(SNOMA  ,1,1)), BUS(x(SNVMA  ,1,1)), &
           BUS(x(SNORO  ,1,1)), BUS(x(SNVRO  ,1,1)), &
           BUS(x(VEGH   ,1,1)), &
           BUS(x(VEGL   ,1,1)), BUS(x(CGSAT  ,1,1)), &
           BUS(x(WSAT   ,1,1)), BUS(x(WWILT  ,1,1)), &
           BUS(x(BCOEF  ,1,1)), &
           BUS(x(CVH    ,1,1)), BUS(x(CVL    ,1,1)), &
           BUS(x(ALVH   ,1,1)), BUS(x(ALVL   ,1,1)), &
           BUS(x(EMISVH ,1,1)), BUS(x(EMISVL ,1,1)), &
           BUS(x(EMISTG ,1,1)), &
           BUS(x(RGLVH  ,1,1)), BUS(x(RGLVL  ,1,1)), &
           BUS(x(STOMRVH,1,1)), BUS(x(STOMRVL,1,1)), &
           BUS(x(GAMVH  ,1,1)), BUS(x(GAMVL  ,1,1)), &
           BUS(x(LAIVH  ,1,1)), BUS(x(LAIVL  ,1,1)),  &
           BUS(x(Z0MVH  ,1,1)),   &
           BUS(x(Z0MVL  ,1,1)), z0m, &
           BUS(x(CLAY   ,1,1)), BUS(x(SAND   ,1,1)), &
           BUS(x(DECIDUOUS,1,1)),BUS(x(EVERGREEN,1,1)), &
           BUS(x(LAIDECI,1,1)), BUS(x(SOILHCAPZ_DRY,1,1)), bus(x(VGH_DENS,1,1)),   &
           BUS(x(SVS_WTA,1,1)),BUS(x(SVS_WTG,1,1)), CG, &
           BUS(x(SOILHCAPZ,1,1)),BUS(x(SOILCONDZ,1,1)),  &
           BUS(x(PSNGRVL,1,1)),  &
           BUS(x(Z0T  ,1,indx_soil)),  &
           BUS(x(ALGR,1,1)),BUS(x(EMISGR,1,1)), &
           BUS(x(ALGRV,1,1)),BUS(x(EMISGRV,1,1)), &
           BUS(x(PSNVH  ,1,1)), BUS(x(PSNVHA ,1,1)), &
           BUS(x(PSURFVHA ,1,1)),         &
           ALVA, BUS(x(LAIVA  ,1,1)), CVPA, EVA, BUS(x(Z0HA ,1,1)),&
           BUS(x(Z0MVG,1,1)), RGLA, STOMRA,   &
           GAMVA,bus(x(CONDSLD    ,1,1)) , bus(x(CONDDRY   ,1,1)), &
           N)
!
      ! Update vegetation temperature with average of low and high vegetation.
      ! VV TO BE MODIFIED: Intermediate step during developement.
      !
      DO I=1,N
          IF (BUS(x(VEGL,I,1)) + BUS(x(VEGH,I,1))  .GT. 0) THEN
              PTVEGE(I)   =  (BUS(x(VEGL,I,1)) *bus(x(TVEGEL,I,1)) + BUS(x(VEGH,I,1)) *bus(x(TVEGEH,I,1)) )/ &
                                                 (BUS(x(VEGL,I,1)) + BUS(x(VEGH,I,1)))
          ELSE ! There is no low or high veg, here we are just putting a value in PTVEGE equal to the skin temp of bg
              PTVEGE(I)   =  bus(x(TGROUND,I,1))
          ENDIF
      ENDDO
!
      CALL VEGI_SVS2 ( zfsolis,   &
           tt                  , PTVEGE,   &
           hu                  , ps                  ,   &
           BUS(x(WSOIL ,1,1)),  &
           RGLA                ,  &
           bus(x(LAIVA  ,1,1))     , bus(x(LAIVH   ,1,1)),   &
           bus(x(LAIVL   ,1,1)), BUS(x(HVEGLPOL,1,1)),  STOMRA,     &
           GAMVA, bus(x(WWILT   ,1,1)),      &
           bus(x(WFC     ,1,1)), SUNCOSA,     &
           bus(x(ROOTDP     ,1,1)),  bus(x(D50   ,1,1)),    &
           bus(x(D95   ,1,1)),  BUS(x(PSNGRVL,1,1)), &
           BUS(x(VEGH   ,1,1)), BUS(x(VEGL   ,1,1)), &
           BUS(x(Z0MVH  ,1,1)), bus(x(VGH_HEIGHT   ,1,1)),bus(x(VGH_DENS,1,1)), &
           bus(x(SNCMA     ,1,1)), bus(x(WVEG_VH,1,1)),bus(x(RST     ,1,1)),     &
           bus(x(SKYVIEW ,1,1)), bus(x(SKYVIEWA ,1,1)), &
           bus(x(VEGTRANS,1,1)), bus(x(VEGTRANSA,1,1)),   &
           bus(x(frootd   ,1,1)), bus(x(acroot ,1,1)), WRMAX_VL, &
           WRMAX_VH, PHM_CAN,BUS(x(HVEGAPOL,1,1)), PSCAP, N)

      If ( (.not.atm_external) .AND. (kount.EQ.0) ) then
         ! GEM first timestep
         ! long-term ... define default value for rcctem in inisurf
         DO I=1,N
            STOM_RS(I) = bus(x(RST,I,1))
         ENDDO
      else if (atm_external .and. KOUNT.EQ.1) then
         !SPS first timestep 
         DO I=1,N
            STOM_RS(I) = bus(x(RST,I,1))
         ENDDO
      else
         IF( USE_PHOTO ) THEN
            DO I=1,N
               STOM_RS(I) =  bus(x(RCCTEM,I,1))
            END DO
         ELSE
            DO I=1,N
               STOM_RS(I) = bus(x(RST,I,1))
            END DO
         ENDIF
      endif


!
!      Effect of high vegetation on met forcing
!
     IF (CANO_REF_FORCING .EQ. 'FOR') THEN ! Forcing in the forest, do not modify them
        DO I=1,N
            PUREF_VEG(I) = bus(x(zusl,I,1))
            PTREF_VEG(I) =  bus(x(ztsl,I,1))
            bus(x(TCA,I,1))  = TT(I)  !
            bus(x(QCA,I,1))  = HU(I)  ! Air specific humidity in the canopy
            bus(x(VCA,I,1))  = VMOD(I)
            PWIND_TOP(I) = VMOD(I)
            bus(x(VCA_DRIFT,I,1)) = VMOD(I) ! Wind speed forcing used for the snow drift routine
            bus(x(SWCA,I,1))    = zfsolis(I)
            bus(x(LWCA,I,1))  = bus(x(FDSI,I,1))
         ENDDO
      ELSE ! ABV or O2F
         CALL  CANOPY_MET_SVS2(tt, hu, vmod, zfsolis, bus(x(FDSI,1,1)), &
                         bus(x(TVEGEH,1,1)),bus(x(zusl,1,1)),  bus(x(ztsl,1,1)),  &
                         SUNCOSA, bus(x(VGH_HEIGHT,1,1)),bus(x(VGH_DENS,1,1)), &
                         BUS(x(Z0MVH  ,1,1)),PZ0LOC_SNOW, BUS(x(SVS_WTG,1,1)), &
                     	 BUS(x(LAIVH  ,1,1)),bus(x(SKYVIEW,1,1)),BUS(x(EMISVH ,1,1)) , bus(x(SWCA,1,1)), bus(x(LWCA,1,1)), &
                         bus(x(VCA,1,1)), bus(x(TCA,1,1)), bus(x(QCA,1,1)) , PWIND_TOP,  &
                         PUREF_VEG,PTREF_VEG,  bus(x(VCA_DRIFT,1,1)), N)
      ENDIF

      IF(LSNOW_INTERCEPTION_SVS2) THEN

         CALL SNOW_INTERCEPTION_SVS2(DT,bus(x(TVEGEH,1,1)), tt, hu, ps, PWIND_TOP,zfsolis,RHOA,     &
                           rainrate_mm,snowrate_mm, bus(x(SNCMA     ,1,1)), wrmax_vh, bus(x(SKYVIEW,1,1)),&
                           bus(x(ESNC     ,1,1)), bus(x(ESNCAF     ,1,1)),  BUS(x(LAIVH  ,1,1)),   &
                           BUS(x(SVS_WTG,1,1)),PHM_CAN, BUS(x(VGH_DENS   ,1,1)), PSCAP,   &
                           bus(x(wveg_vh  ,1,1)), rainrate_mm_veg, snowrate_mm_veg, punload_forest,      &
                           PFCANS, N)

      ELSE
         DO I=1,N
            ! Rainfall and snowfall rate below high-vegetation are not impacted by the presence of high-vegetation

            rainrate_mm_veg(i) = rainrate_mm(i)
            snowrate_mm_veg(i) = snowrate_mm(i)
         ENDDO

      ENDIF


      ! Store rainfall and snowfall rate below high vegetation (in m) to be consistent with rainrate and snowrate in the bus
      ! Unloading is added for mass conservation. 
      DO I=1,N
         ! Set snowfall and rainfall rate below vegetation to zero if no high vegetation is present. 
         ! This is used to make sure that Crocus is not called to simulate snowpack evolution below high vegetation 
         ! when high vegetation is not present in a grid cell. 
         IF ( BUS(x(VEGH,I,1)) .LT. EPSILON_SVS) THEN
             rainrate_mm_veg(i) = 0. 
             snowrate_mm_veg(i) = 0. 
         ENDIF
         bus(x(rainrate_vgh,i,1))  = rainrate_mm_veg(i)/1000.
         bus(x(snowrate_vgh,i,1))  = snowrate_mm_veg(i)/1000. + punload_forest(i)/1000.
      ENDDO
!
      CALL DRAG_SVS2 ( bus(x(TGROUND,1,1)),bus(x(TGROUNDV,1,1))  , &
           bus(x(TVEGEL,1,1)), bus(x(TVEGEH,1,1)), bus(x(TSNOWV_SVS,1,1)), &
           bus(x(TSNOW_SVS,1,1)),bus(x(WSOIL ,1,1)) ,  &
           bus(x(WVEG_VL,1,1)),bus(x(WVEG_VH,1,1)),  zthetaa,  &
           VMOD, VDIR, hu, RHOA,    &
           ps, STOM_RS,   &
           z0m, z0mland, bus(x(Z0MVG,1,1)), bus(x(WFC,1,1)),      &
           bus(x(WSAT,1,1)),  bus(x(CLAY,1,1)), bus(x(SAND,1,1)), &
           bus(x(LAIVL,1,1)),bus(x(LAIVH,1,1)), WRMAX_VL, WRMAX_VH, &
           bus(x(zusl,1,1)), bus(x(ztsl,1,1)),    &
           bus(x (DLAT,1,1)), bus(x(PSNVH ,1,1)),&
           bus(x(FCOR,1,1)),bus(x(Z0HA ,1,1)), BUS(x(SVS_WTG,1,1)), &
           bus(x(VGH_DENS,1,1)), BUS(x(Z0MVH  ,1,1)),  BUS(x(Z0MVL  ,1,1)), PZ0LOC_SNOW, PZ0H_SNOW, &
           bus(x(VGH_HEIGHT   ,1,1)),BUS(x(LAIVH  ,1,1)), bus(x(VCA,1,1)),PFCANS,bus(x(SNCMA,1,1)),  &
           bus(x(RESAGR,1,1)),bus(x(RESAGRV,1,1)), &
           bus(x(RESA_VL,1,1)),bus(x(RESA_VH,1,1)), pres_snca, bus(x(RESASA,1,1)), bus(x(RESASV,1,1)), &
           bus(x(HUSURF,1,1)),bus(x(HUSURFGV,1,1)),   &
           HRSURF, HRSURFGV,      &
           bus(x(HV_VL,1,1)),bus(x(HV_VH,1,1)), HVSN_VH, DEL_VL, DEL_VH,     &
           bus(x(Z0HBG,1,1)), bus(x(Z0HVL,1,1)), bus(x(Z0HVH,1,1)), bus(x(Z0HGV,1,1)), &
            N )
      if (phy_error_L) return




!     Snow over bare/low ground

      ! Compute wind speed for snow drift effect and unloading term 
      do I=1,N
            PWIND_DRIFT_OPEN(I) = VMOD(I)
            PUNLOAD_OPEN(I) = 0. ! No unloading in open terrain
      enddo

      CALL SNOW_SVS2(   bus(x(SNOMA_SVS,1,1)), bus(x(TSNOW_SVS,1,1)), bus(x(WSNOW_SVS,1,1)),    &
                         bus(x(SNODEN_SVS,1,1)),  bus(x(SNOAL,1,1)),bus(x(SNOAGE_SVS,1,1)),    &
                         bus(x(SNODIAMOPT_SVS,1,1)), bus(x(SNOSPHERI_SVS,1,1)),bus(x(SNOHIST_SVS,1,1)),   &
                         DT, bus(x(TPSOIL    ,1,1)) ,  PCT, bus(x(SOILHCAPZ,1,1)), bus(x(SOILCONDZ,1,1)),                 &
                         ps,tt,zfsolis,     &
                         hu, VMOD, PWIND_DRIFT_OPEN, &
                         bus(x(FDSI,1,1)),         &
                         RAINRATE_MM, SNOWRATE_MM,PUNLOAD_OPEN,bus(x(RESASA,1,1)),                   &
                         RHOA, bus(x(zusl,1,1)),  bus(x(ztsl,1,1)),             &
                         BUS(X(ALGR,1,1)), PD_G, PDZG,                          &
                         bus(x(RSNOWSA,1,1)), bus(x(GFLUXSA,1,1)),bus(x(RNETSA,1,1)),bus(x(HFLUXSA,1,1)) , &
                         PGFLUXSNOW,bus(x(SWNETSA,1,1)), bus(x(LWNETSA,1,1)), bus(x(SUBLDRIFTA,1,1)), &
                         bus(x(HPSA ,1,1)),  bus(x(PSNGRVL ,1,1)), PZ0AVG_SNOW,PZ0LOC_SNOW,PZ0H_SNOW, &
                         LESNOFRAC, LESLNOFRAC, bus(x(ESA,1,1)), PZENITH, &
                         bus(x (DLAT,1,1)), bus(x (DLON,1,1)),PFOREST,bus(x(SNOTYPE_SVS,1,1)),  &
                         BUS(x(HVEGAPOL,1,1)),BUS(x(AGINGCOEF,1,1)),N, NL_SVS)
      if (phy_error_L) return




! Define temperature use as a lower boundary condition for the snowpack below high vegetation
      DO I=1,N
          DO J=1,NL_SVS
             PSOIL_TEMP_VGH(I,J) = bus(x(TPSOIL,I,J))
          ENDDO
      ENDDO

!
!     Snow under high veg  as in SVS1

      !
      ! The effect of low basal vegetation is also considered for snow in forested environment
      DO I = 1,N
          PHVEGAPOL_V(I) = BUS(x(HVEGAPOL  ,I,1))
      ENDDO

      CALL SNOW_SVS2(   bus(x(SNOMAV_SVS,1,1)), bus(x(TSNOWV_SVS,1,1)), bus(x(WSNOWV_SVS,1,1)),    &
                             bus(x(SNODENV_SVS,1,1)), bus(x(SNVAL,1,1)),bus(x(SNOAGEV_SVS,1,1)),    &
                             bus(x(SNODIAMOPTV_SVS,1,1)), bus(x(SNOSPHERIV_SVS,1,1)),bus(x(SNOHISTV_SVS,1,1)),   &
                             DT,PSOIL_TEMP_VGH, PCT, bus(x(SOILHCAPZ,1,1)), bus(x(SOILCONDZ,1,1)),               &
                             ps, bus(x(TCA,1,1)),bus(x(SWCA,1,1)),     &
                             bus(x(QCA,1,1)), bus(x(VCA,1,1)), bus(x(VCA_DRIFT,1,1)), &
                             bus(x(LWCA,1,1)),         &
                             RAINRATE_MM_VEG, SNOWRATE_MM_VEG,PUNLOAD_FOREST, bus(x(RESASV,1,1)),              &
                             RHOA,  PUREF_VEG,   PTREF_VEG,            &
                             BUS(X(ALGR,1,1)), PD_G, PDZG,                          &
                             bus(x(RSNOWSV,1,1)), bus(x(GFLUXSV,1,1)),bus(x(RNETSV,1,1)) , bus(x(HFLUXSV ,1,1)), &
                             PGFLUXSNOW_V,bus(x(SWNETSV,1,1)),bus(x(LWNETSV,1,1)),bus(x(SUBLDRIFTV,1,1)), &
                             bus(x(HPSV ,1,1)),bus(x(PSNVH ,1,1)), PZ0AVG_SNOW,PZ0LOC_SNOW,PZ0H_SNOW,  &
                             LESVNOFRAC, LESVLNOFRAC, bus(x(ESV,1,1)),PZENITH, &
                             bus(x (DLAT,1,1)), bus(x (DLON,1,1)), PFOREST_V,bus(x(SNOTYPEV_SVS,1,1)), &
                             PHVEGAPOL_V,BUS(x(AGINGCOEF,1,1)), N, NL_SVS)


      if (phy_error_L) return

! Compute snow diagnostics for hydro and outputs
!
      DO I=1,N
!        total snow mass
         bus(x(SNOMA,I,1))  = 0.
         bus(x(SNVMA,I,1)) = 0
!        total snow depth
         bus(x(SNODPL,I,1))  = 0.
         bus(x(SNVDP ,I,1))   = 0.
!        total snow liquid water content
         bus(x(WSNOW,I,1))  = 0.
         bus(x(WSNV ,I,1))   = 0.

         DO J=1,NSL
            bus(x(SNOMA,I,1))   =  bus(x(SNOMA,I,1)) + bus(x(SNOMA_SVS ,I,J))
            bus(x(SNVMA,I,1))   =  bus(x(SNVMA,I,1)) + bus(x(SNOMAV_SVS,I,J))
            bus(x(SNODPL,I,1))  =  bus(x(SNODPL,I,1)) + bus(x(SNOMA_SVS ,I,J))/bus(x(SNODEN_SVS ,I,J))
            bus(x(SNVDP ,I,1))  =  bus(x(SNVDP ,I,1)) + bus(x(SNOMAV_SVS,I,J))/bus(x(SNODENV_SVS,I,J))
            bus(x(WSNOW,I,1))   =  bus(x(WSNOW,I,1)) + bus(x(WSNOW_SVS ,I,J))*1000.
            bus(x(WSNV,I,1))   =  bus(x(WSNV,I,1)) + bus(x(WSNOWV_SVS ,I,J))*1000.
         ENDDO

!        Cumulated liquid water runoff leaving the snowpack
         bus(x(RSNOWS_ACC,I,1)) = bus(x(RSNOWS_ACC,I,1)) + bus(x(RSNOWSA,I,1))*DT
         bus(x(RSNOWSV_ACC,I,1)) = bus(x(RSNOWSV_ACC,I,1)) + bus(x(RSNOWSV,I,1))*DT

         ! Total latent heat flux from snow
         bus(x(LFLUXSA ,I,1))  = LESNOFRAC(I) + LESLNOFRAC(I)
         bus(x(LFLUXSV ,I,1))  = LESVNOFRAC(I) + LESVLNOFRAC(I)

      ENDDO

!

      CALL EBUDGET_SVS2(bus(x(TSA ,1,1)),  &
                  bus(x(WSOIL     ,1,1)) , bus(x(ISOIL,1,1)),  &
                  bus(x(TGROUND   ,1,1)) , bus(x(TGROUNDV,1,1)),  &
                  bus(x(TVEGEL    ,1,1)) , bus(x(TVEGEH  ,1,1)) ,    &
                  bus(x(TPSOIL    ,1,1)) ,    &
                  bus(x(TPERM     ,1,1)) , bus(x(GFLUXSA,1,1)), bus(x(GFLUXSV,1,1)), &
                  DT                     , VMOD, VDIR, bus(x(DLAT,1,1)),     &
                  zfsolis, bus(x(SWCA,1,1)),ALVA ,bus(x(laiva,1,1)),         &
                  GAMVA , BUS(x(ALVL,1,1)), &
                  BUS(x(ALVH,1,1)), BUS(x(ALGR,1,1)), BUS(x(EMISGR,1,1)),    &
                  BUS(x(ALGRV,1,1)), BUS(x(EMISGRV,1,1)),    &
                  bus(x(FDSI       ,1,1)), bus(x(LWCA,1,1)), zthetaa ,    &
                  bus(x(FCOR       ,1,1)), bus(x(zusl,1,1)),    &
                  bus(x(ztsl       ,1,1)), hu, &
                  ps, RHOA, BUS(x(SVS_WTA,1,1)), BUS(x(SVS_WTG,1,1)),  bus(x(VGH_DENS,1,1)), &
                  z0m, z0mland , bus(x(Z0T,1,indx_soil)),&
                  HRSURF,HRSURFGV,       &
                  bus(x(HV_VL,1,1)) , bus(x(HV_VH,1,1)), HVSN_VH, DEL_VL, DEL_VH, STOM_RS ,&
                  CG,CVPA,BUS(x(EMISVL ,1,1)), BUS(x(EMISVH ,1,1)) ,  BUS(x(SKINCOND_VL ,1,1)),  &
                  bus(x(RESAGR,1,1)), bus(x(RESA_VL,1,1)),bus(x(RESA_VH,1,1)),   &
                  bus(x(RESASA,1,1)), bus(x(RESASV,1,1)) ,bus(x(RESAGRV,1,1)),pres_snca, &
                  bus(x(RNETSA     ,1,1)) , bus(x(HFLUXSA,1,1)),   &
                  LESLNOFRAC, LESNOFRAC        , bus(x(ESA,1,1)), bus(x(SUBLDRIFTA,1,1)),  &
                  bus(x(SNOAL      ,1,1)) ,  bus(x(RSNOWSA,1,1)),   &
                  bus(x(TSNOW_SVS  ,1,1)) ,    &
                  bus(x(RNETSV     ,1,1)) , bus(x(HFLUXSV ,1,1)),   &
                  LESVLNOFRAC, LESVNOFRAC              , bus(x(ESV,1,1)),  bus(x(SUBLDRIFTV,1,1)),  &
                  bus(x(SNVAL      ,1,1)) ,    &
                  bus(x(TSNOWV_SVS ,1,1)) , PHM_CAN,  bus(x(SNCMA     ,1,1)), &
                  bus(x(VGH_HEIGHT   ,1,1)),  &
                  bus(x(SKYVIEW   ,1,1)), bus(x(SKYVIEWA   ,1,1)),  PFCANS, &
                  bus(x(SOILHCAPZ ,1,1)) ,bus(x(SOILCONDZ,1,1)),   &
                  rainrate_mm,bus(x(WVEG_VL,1,1)),bus(x(WVEG_VH,1,1)), &
                  bus(x(snoma,1,1)), bus(x(snvma,1,1)),&
                  bus(x(VEGTRANSA  ,1,1)) , bus(x(ALVIS,1,indx_soil)),     &
                  bus(x(RNET_S     ,1,1)),    &
                  bus(x(FC  ,1,indx_soil)), bus(x(FV  ,1,indx_soil)),   &
                  bus(x(LEG        ,1,1)) , bus(x(LEVL  ,1,1)), bus(x(LEVH ,1,1)),    &
                  bus(x(LES        ,1,1)) , bus(x(LESV   ,1,1)),    &
                  bus(x(LEGV       ,1,1)) ,  &
                  bus(x(LER_VL        ,1,1)) , bus(x(LETR_VL       ,1,1)) ,   &
                  bus(x(LER_VH        ,1,1)) , bus(x(LETR_VH       ,1,1)) ,   &
                  bus(x(EG            ,1,1)) , bus(x(EGV            ,1,1)) ,   &
                  bus(x(ER_VL         ,1,1)) , bus(x(ETR_VL    ,1,1)),    &
                  bus(x(ER_VH         ,1,1)) , bus(x(ETR_VH    ,1,1)),  bus(x(ESNC     ,1,1)), &
                  bus(x(FL         ,1,1)),  bus(x(EFLUX      ,1,1)) ,    &
                  bus(x(BM         ,1,1)) , bus(x(FQ   ,1,1)),    &
                  bus(x(bt, 1,indx_soil)) , bus(x(RESAEF,1,1)),   &
                  LEFF                    ,    &
                  bus(x(FTEMP,1,indx_soil)), BUS(x(FVAP,1,indx_soil)),   &
                  bus(x(qsurf,1,indx_soil)), bus(x(frv ,1,indx_soil)),   &
                  bus(x(ALFAT      ,1,1)) , bus(x(ALFAQ      ,1,1)) ,    &
                  bus(x(ilmo  ,1,indx_soil)), bus(x(hst  ,1,indx_soil)), &
                  TRAD, N,   &
                  bus(x(QVEG ,1,1)), bus(x(QGV   ,1,1)), bus(x(QGR   ,1,1)), &
                  RPP, bus(x(Z0HA ,1,1)))


      ! Update vegetation temperature with average of low and high vegetation.
      ! Update aerodynamical resistance  with average of low and high vegetation.
      ! VV TO BE MODIFIED: Intermediate step during developement.
      !
      DO I=1,N
          IF (BUS(x(VEGL,I,1)) + BUS(x(VEGH,I,1))  .GT. 0) THEN
              PTVEGE(I)   =  (BUS(x(VEGL,I,1)) *bus(x(TVEGEL,I,1)) + BUS(x(VEGH,I,1)) *bus(x(TVEGEH,I,1)) )/ &
                                                 (BUS(x(VEGL,I,1)) + BUS(x(VEGH,I,1)))
              BUS(x(RESAVG,I,1))   =  (BUS(x(VEGL,I,1)) *bus(x(RESA_VL,I,1)) + BUS(x(VEGH,I,1)) *bus(x(RESA_VH,I,1)) )/ &
                                                 (BUS(x(VEGL,I,1)) + BUS(x(VEGH,I,1)))
          ELSE ! There is no low or high veg, here we are just putting a value in PTVEGE equal to the skin temp of bg
              PTVEGE(I)   =  bus(x(TGROUND,I,1))
              BUS(x(RESAVG,I,1))  = BUS(x(RESAGR,I,1))
          ENDIF
      ENDDO


      if (phy_error_L) return
!
!

      CALL WATSURF_BUDGET_SVS2 ( DT,      &
           bus(x(ESNC     ,1,1)), bus(x(ESNCAF     ,1,1)), &
           bus(x(eg      ,1,1)), bus(x(egv      ,1,1)),   &
           bus(x(er_vl      ,1,1)),                       &
           bus(x(er_vh   ,1,1)),bus(x(etr_vl      ,1,1)), &
           bus(x(etr_vh  ,1,1)), rainrate_mm, rainrate_mm_veg ,&
           bus(x(rsnowsa ,1,1)), bus(x(rsnowsv ,1,1)),&
           bus(x(svs_wta ,1,1)),&
           bus(x(svs_wtg ,1,1)), bus(x(acroot  ,1,1)),&
           wrmax_vl,wrmax_vh,  &
           bus(x(snoma   ,1,1)), bus(x(snvma   ,1,1)),&
           bus(x(SNCMA     ,1,1)), &
           bus(x(wveg_vl ,1,1)),bus(x(wveg_vh  ,1,1)),&
           wveglt, wveght                            ,&
           PG, ETR_GRID, eg_grid, &
           N)

      CALL HYDRO_SVS2 ( DT,      &
           bus(x(impervu ,1,1)), PG,  &
           ETR_GRID,  eg_grid,  bus(x(wsat    ,1,1)),&
           bus(x(ksat    ,1,1)), bus(x(psisat  ,1,1)),&
           bus(x(bcoef   ,1,1)), bus(x(fbcof   ,1,1)),&
           bus(x(wfcint  ,1,1)), bus(x(grkef   ,1,1)),&
           bus(x(wsoil   ,1,1)), wsoilt              ,&
           bus(x(isoil   ,1,1)), isoilt              ,&
           bus(x(ksatc   ,1,1)), bus(x(khc     ,1,1)),&
           bus(x(psi     ,1,1)), bus(x(grksat  ,1,1)),&
           bus(x(wfcdp   ,1,1)), bus(x(watflow ,1,1)),&
           bus(x(latflw  ,1,1)),bus(x(runofftot ,1,indx_soil)), &
           bus(x(watpond ,1,1)),bus(x(maxpond ,1,1)), &
           N)

      IF( USE_PHOTO ) THEN

         if(vf_type == "CCILCECO") then
            CALL PHTSYN_SVS_CCILCECO( BUS(x(LAIVF26,1,1))  , BUS(x(VEGF_EVOL   ,1,1)), &
                        PTVEGE  , ps, &
                        BUS(x(RESAVG ,1,1))  , hu, &
                        zFSOLIS              , BUS(x(WSOIL ,1,1)), &
                        BUS(x(FROOTD ,1,1))  , SUNCOSA            , &
                        BUS(x(WFC    ,1,1))  , BUS(x(WWILT  ,1,1)), &
                        MASKLAT50            , BUS(x(VGCTEM ,1,1))  , &
                        BUS(x(LAICTEM,1,1))  ,                      &
                        BUS(x(RCCTEM ,1,1))  , BUS(x(CO2I1  ,1,1)), &
                        BUS(x(AVG_GWSOL,1,1)), &
                        NCLASS, N)

         else

            ! WARNING:
            ! USING VEGF in call below
            ! SHould probably use VEGF_EVOL
            !

            CALL PHTSYN_SVS2 ( BUS(x(LAIVF26,1,1))  , BUS(x(VEGF,1,1)), &
                        PTVEGE  , ps, &
                        BUS(x(RESAVG ,1,1))  , hu, &
                        zFSOLIS              , BUS(x(WSOIL ,1,1)), &
                        BUS(x(FROOTD ,1,1))  , SUNCOSA            , &
                        BUS(x(WFC    ,1,1))  , BUS(x(WWILT  ,1,1)), &
                        MASKLAT50            , BUS(x(VGCTEM ,1,1))  , &
                        BUS(x(LAICTEM,1,1))  ,                      &
                        BUS(x(RCCTEM ,1,1))  , BUS(x(CO2I1  ,1,1)), &
                        BUS(x(AVG_GWSOL,1,1)), &
                        NCLASS, N)

         endif
      ENDIF
!

!
!     Phase change for the soil column
!
      CALL PHASE_CHANGES (DT, bus(x(LAIVA  ,1,1)), BUS(x(SOILHCAPZ,1,1)) , &
                      bus(x(WSAT   ,1,1)), bus(x(PSISAT  ,1,1)), bus(x(BCOEF  ,1,1)), &
                      bus(x(TPSOIL ,1,1)), bus(x(ISOIL  ,1,1)), &
                      wsoilt, WFTG, WDTTG, DELWATGR, DELICEGR     , &
                      bus(x(FROOTD ,1,1)), N                   , &
                      bus(x(PHASEF ,1,1)), bus(x(PHASEM ,1,1)) , &
                      bus(x(DELTAT ,1,1)), bus(x(APPHEATCAP ,1,1)), bus(x(TMAX ,1,1)) )

       ! Update the soil liquid water and ice content after phase changes
       DO I=1,N
         DO J=1,NL_SVS
            WSOILT(I,J) = WDTTG(I,J)
            ISOILT(I,J) = WFTG(I,J)
          END DO
       END DO
!
!     Update prognostic variable in SVS2
!
      CALL UPDATE_SVS2 ( WSOILT, ISOILT, WVEGLT,WVEGHT,   &
                       bus(x(WSOIL   ,1,1)), bus(x(ISOIL   ,1,1)),  &
                       bus(x(WVEG_VL ,1,1)), bus(x(WVEG_VH ,1,1)),  &
                       bus(x(WSOILM  , 1,1)), &
                       N )
!
  
      !# Compute values at the diagnostic level

      ! for now z0m with orography
      ! compute z0h in soili or ebudget
      
      i = sl_sfclayer(zthetaa,hu,vmod,vdir,zzusl,zztsl,ztsa,zqsurf, &
           z0m,z0h,zdlat,zfcor,L_min=sl_Lmin_soil,spdlim=vmod_lmin, &
           hghtm_diag=zu,hghtt_diag=zt,t_diag=ztdiag,q_diag=zqdiag, &
           u_diag=zudiag,v_diag=zvdiag,tdiaglim=TDIAGLIM_FALSE) 
      
      if (i /= SL_OK) then
         call physeterror('svs', 'error 2 returned by sl_sfclayer()')
         return
      endif

      if (sl_Lmin_soil > 0.) then
         ! re-scale diagnostic winds 
         zudiag = zudiag * vmod / vmod_lmin
         zvdiag = zvdiag * vmod / vmod_lmin
      endif

   !# Fill surface type-specific diagnostic values
   zqdiagtyp = zqdiag
   ztdiagtyp = ztdiag
   zudiagtyp = zudiag
   zvdiagtyp = zvdiag

      do i=1,n
!
!
        ! TO DO - INITIALIZE EMISR WITH MEAN SURFACE EMISSIVITY
        !bus(x(emisr  ,i,1        )) = bus(x(emis ,i,1        ))   
        bus(x(tsurf  ,i,indx_sfc )) = bus(x(tsa  ,i,1        ))
        bus(x(tsrad  ,i,1        )) = TRAD(i)
!
!       CALCULATE LAND-ATMOSPHERE OUTCOMING WATER FLUX
        BUS(x(WFLUX,I,1)) = RHOA(I)*BUS(x(EFLUX,I,1))
        BUS(x(ACCEVAP,I,1)) = BUS(x(ACCEVAP,I,1)) + BUS(x(WFLUX,I,1)) * DT
!
!       CALCULATE MEAN SNOW DEPTH FOR ESTHETIC PURPOSE ONLY
        zsnodp(i) = bus(x(VEGH,i,1)) * bus(x(SNVDP,i,1)) + (1. -  bus(x(VEGH,i,1))) * bus(x(SNODPL,i,1))
      end do
!
!     FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
      CALL FILLAGG ( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_SOIL,  &  
                    SURFLEN )
!



      RETURN
    END subroutine svs2
