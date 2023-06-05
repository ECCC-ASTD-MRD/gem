
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
subroutine svs(BUS, BUSSIZ, PTSURF, PTSURFSIZ, DT, KOUNT, TRNCH, N, M, NK)
   use, intrinsic :: iso_fortran_env, only: INT64
   use phy_status, only: phy_error_L
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
   use sfcbus_mod
   use sfc_options, only: atm_external, atm_tplus, radslope, jdateo, &
        use_photo, nclass, zu, zt, sl_Lmin_soil, VAMIN, svs_local_z0m, &
        vf_type,lsoil_freezing_svs1,lwater_ponding_svs1,critwater
   use svs_configs
   implicit none
!!!#include <arch_specific.hf>
!
!Author
!          S. Belair (January 1997), M. Abrahamowicz, S.Z. Husain (2012) 
!Revisions
!
! 001      Rewrite ISBA into SVS (multi-budget, multi-layer isba)
!              Add sea ice surface
!Object
!          Multitasking of the surface scheme SVS
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

   integer i,m, masklat50(n)

   real,dimension(n) :: alva, cg, cvpa, del, dwaterdt
   real,dimension(n) :: esnofrac, esvnofrac, eva, gamva, hrsurf
   real,dimension(n) :: leff, lesnofrac, lesvnofrac, rainrate_mm
   real,dimension(n) :: rgla, rhoa, snowrate_mm, stom_rs, stomra
   real,dimension(n) :: suncosa, sunother1, sunother2, sunother3
   real,dimension(n) :: sunother4, trad, tva, vdir, vmod, vmod_lmin, wrmax, wvegt
   real,dimension(n) :: wsaturc1
! 
   real, dimension(n,nl_svs) :: isoilt, wsoilt
!
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

      IF(lwater_ponding_svs1 .and. kount==1) THEN
          DO I=1,N
!           EG: Adjust max. ponding depth according to bare ground fraction: consider 10mm over bare ground
	    zmaxpond(I) = zmaxpond(I) * (zvegh(I)+zvegl(I)) + 0.01 * (1.-zvegh(I)-zvegl(I))
!	    EG: Adjust max. ponding depth according to slope
            zmaxpond(I) = max(0.0,zmaxpond(I)*(1.0E-10)**zslop(I))
         END DO
      ENDIF
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
!
! EG CODE RELATED TO PONDING OPTION
      IF (lwater_ponding_svs1) then
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

      CALL SOILI_SVS( BUS(x(WSOIL ,1,1)), &
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
           BUS(x(LAIDECI,1,1)),   &
           BUS(x(SVS_WTA,1,1)), CG, &
           BUS(x(PSNGRVL,1,1)),  &  
           BUS(x(Z0T  ,1,indx_soil)),  & 
           BUS(x(ALGR,1,1)),BUS(x(EMISGR,1,1)), &  
           BUS(x(PSNVH  ,1,1)), BUS(x(PSNVHA ,1,1)), &  
           ALVA, BUS(x(LAIVA  ,1,1)), CVPA, EVA, BUS(x(Z0HA ,1,1)),&
           BUS(x(Z0MVG,1,1)), RGLA, STOMRA,   &
           GAMVA, N,    &
           BUS(x(SOILHCAPZ,1,1)), BUS(x(SOILCONDZ,1,1)), BUS(x(CONDDRY   ,1,1)), BUS(x(CONDSLD  ,1,1)) )
!
!     

      CALL VEGI_SVS ( zfsolis,   &
           tt                  , bus(x(TVEGE   ,1,1)),   &  
           hu                  , ps                  ,   &  
           BUS(x(WSOIL ,1,1)),  &
           RGLA                ,  &   
           bus(x(LAIVA  ,1,1))     , bus(x(LAIVH   ,1,1)),   &  
           STOMRA,     &
           GAMVA, bus(x(WWILT   ,1,1)),      &
           bus(x(WFC     ,1,1)), SUNCOSA,     &
           bus(x(ROOTDP     ,1,1)),  bus(x(D50   ,1,1)),    &
           bus(x(D95   ,1,1)),  BUS(x(PSNGRVL,1,1)), &
           BUS(x(VEGH   ,1,1)), BUS(x(VEGL   ,1,1)), &
           bus(x(RST     ,1,1)),     &
           bus(x(SKYVIEW ,1,1)), bus(x(VEGTRANS,1,1)),   &   
           bus(x(frootd   ,1,1)), bus(x(acroot ,1,1)), WRMAX, N)

      IF(KOUNT.EQ.1) then
         DO I=1,N
            STOM_RS(I) = bus(x(RST,I,1))
         ENDDO
         ! long-term ... define default value for rcctem in inisurf
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
      CALL DRAG_SVS ( bus(x(TGROUND,1,1)), bus(x(TVEGE,1,1)),  &   
           bus(x(WSOIL ,1,1)) ,  &   
           bus(x(WVEG   ,1,1)), zthetaa,  &   
           VMOD, VDIR, hu,     &
           ps, STOM_RS,   &  
           z0m, z0mland, bus(x(Z0MVG,1,1)), bus(x(WFC,1,1)),      &
           bus(x(WSAT,1,1)),  bus(x(CLAY,1,1)), bus(x(SAND,1,1)), &
           bus(x(LAIVA,1,1)), WRMAX, bus(x(zusl,1,1)), bus(x(ztsl,1,1)),    & 
           bus(x (DLAT,1,1)), &
           bus(x(FCOR,1,1)),bus(x(Z0HA ,1,1)), &  
           bus(x(RESAGR,1,1)), bus(x(RESAVG,1,1)), &    
           bus(x(HUSURF,1,1)),   &  
           HRSURF,      &
           bus(x(HV,1,1)), DEL,     &
           bus(x(Z0HBG,1,1)), bus(x(Z0HVG,1,1)), &
           N )  
      if (phy_error_L) return


      CALL SNOW_ALONE ( bus(x(TSNOW,1,1)), bus(x(TSNOW,1,2)),  &     
                     bus(x(SNORO,1,1)),    &   
                     bus(x(SNOAL,1,1)), bus(x(WSNOW,1,1)),  &     
                     bus(x(SNODPL,1,1)),   &   
                     bus(x(SNOMA,1,1)),   &   
                     ps , VMOD, VDIR,  RHOA,   &    
                     zthetaa,   &   
                     zfsolis, bus(x(FDSI,1,1)),  &    
                     hu, rainrate_mm,  &      
                     snowrate_mm,  tt,  &     
                     bus(x(TDIAG,1,1)),     &   
                     bus(x(UDIAG,1,1)), bus(x(VDIAG,1,1)),  &     
                     bus(x(TSNAVG,1,1)),&
                     bus(x(MELTS,1,1)), bus(x(MELTSR,1,1)),   &    
                     bus(x(RNETSA,1,1)), bus(x(HFLUXSA,1,1)),  &     
                     LESNOFRAC, ESNOFRAC,   &      
                     bus(x(RSNOWSA,1,1)), bus(x(SNODEN,1,1)), bus(x(RESASA,1,1)),  &    
                     DT, z0m, z0mland, bus(x(FCOR,1,1)),  &     
                     bus(x(zusl,1,1)), bus(x(ztsl,1,1)), &
                     bus(x (DLAT,1,1)), bus(x(PSNGRVL ,1,1)), N)     
      if (phy_error_L) return
!    
!
      CALL SNOW_VEG ( bus(x(TSNOWVEG  ,1,1)), bus(x(TSNOWVEG,1,2)),  &  
                     bus(x(SNVRO     ,1,1)),   & 
                     bus(x(SNVAL     ,1,1)), bus(x(WSNV    ,1,1)),   & 
                     bus(x(SNVDP     ,1,1)),    &
                     bus(x(SNVMA     ,1,1)),    &
                     ps , VMOD, VDIR, RHOA,   &  
                     zthetaa,    &
                     zfsolis, bus(x(FDSI    ,1,1)),  &  
                     hu , rainrate_mm,  &  
                     snowrate_mm, tt,  &  
                     bus(x(TDIAG     ,1,1)),     &
                     bus(x(UDIAG     ,1,1)), bus(x(VDIAG   ,1,1)),   & 
                     BUS(x(VEGH   ,1,1)) ,   & 
                     bus(x(SKYVIEW   ,1,1)), bus(x(VEGTRANS,1,1)),   & 
                     bus(x(TVEGE     ,1,2)), bus(x(EMISVH ,1,1)),   &
                     bus(x(WVEG      ,1,1)), bus(x(TSVAVG,1,1)), &
                     bus(x(RNETSV    ,1,1)), bus(x(HFLUXSV ,1,1)),  &  
                     LESVNOFRAC            , ESVNOFRAC,     & 
                     bus(x(RSNOWSV   ,1,1)), bus(x(SNVDEN  ,1,1)), bus(x(RESASV,1,1)), &  
                     DT,    &
                     z0m, z0mland, &
                     bus(x(FCOR   ,1,1)), bus(x (DLAT,1,1)),  & 
                     bus(x(zusl,1,1)), bus(x(ztsl,1,1)), bus(x(PSNVHA ,1,1)), N)
      if (phy_error_L) return

!
      CALL EBUDGET_SVS(bus(x(TSA ,1,1)),  &  
                  bus(x(WSOIL     ,1,1)) , bus(x(ISOIL,1,1)),  &   
                  bus(x(TGROUND    ,1,1)) , bus(x(TGROUND,1,2)),   & 
                  bus(x(TVEGE      ,1,1)) , bus(x(TVEGE,1,2)),   &   
                  DT                     , VMOD, VDIR, bus(x(DLAT,1,1)),   &   
                  zfsolis ,ALVA ,bus(x(laiva,1,1)),GAMVA ,     & 
                  BUS(x(ALGR,1,1))        , BUS(x(EMISGR,1,1)),    & 
                  bus(x(FDSI       ,1,1)) , zthetaa ,    &   
                  bus(x(FCOR       ,1,1)) , bus(x(zusl,1,1)),    &  
                  bus(x(ztsl       ,1,1)) , hu, &
                  ps, RHOA, BUS(x(SVS_WTA,1,1)), &
                  z0m, z0mland , bus(x(Z0T,1,indx_soil)),&
                  HRSURF,       & 
                  bus(x(HV         ,1,1)) , DEL, STOM_RS ,& 
                  CG,CVPA,EVA,bus(x(PSNGRVL    ,1,1)) ,    &    
                  bus(x(RESAGR,1,1)), bus(x(RESAVG,1,1)),   &        
                  bus(x(RESASA,1,1)), bus(x(RESASV,1,1)), &
                  bus(x(RNETSA     ,1,1)) , bus(x(HFLUXSA,1,1)),   &   
                  LESNOFRAC               , ESNOFRAC,   &   
                  bus(x(SNOAL      ,1,1)) ,    &  
                  bus(x(TSNOW      ,1,1)) ,    &  
                  bus(x(RNETSV     ,1,1)) , bus(x(HFLUXSV ,1,1)),   &   
                  LESVNOFRAC              , ESVNOFRAC,    &    
                  bus(x(SNVAL      ,1,1)) ,    &  
                  bus(x(TSNOWVEG   ,1,1)) ,   &   
                  bus(x(VEGH       ,1,1)) , bus(x(VEGL   ,1,1)),   &   
                  bus(x(PSNVH      ,1,1)) ,    &   
                  bus(x(PSNVHA     ,1,1)),  &   
                  rainrate_mm,bus(x(WVEG   ,1,1)),bus(x(snoma,1,1)),&
                  bus(x(snvma,1,1)),&
                  bus(x(ALVIS,1,indx_soil)),     & 
                  bus(x(RNET_S     ,1,1)),    &   
                  bus(x(FC  ,1,indx_soil)), bus(x(FV  ,1,indx_soil)),   &    
                  bus(x(LEG        ,1,1)) , bus(x(LEV  ,1,1)),    &   
                  bus(x(LES        ,1,1)) , bus(x(LESV   ,1,1)),    &  
                  bus(x(LER        ,1,1)) , bus(x(LETR       ,1,1)) ,   &  
                  bus(x(EG         ,1,1)) ,   &    
                  bus(x(ER         ,1,1)) , bus(x(ETR    ,1,1)),    &  
                  bus(x(FL         ,1,1)),  bus(x(EFLUX      ,1,1)) ,    &  
                  bus(x(BM         ,1,1)) , bus(x(FQ   ,1,1)),    &  
                  bus(x(bt, 1,indx_soil)) , bus(x(RESAEF,1,1)),   &  
                  LEFF                    , DWATERDT,     & 
                  bus(x(FTEMP,1,indx_soil)), BUS(x(FVAP,1,indx_soil)),   &   
                  bus(x(qsurf,1,indx_soil)), bus(x(frv ,1,indx_soil)),   &   
                  bus(x(ALFAT      ,1,1)) , bus(x(ALFAQ      ,1,1)) ,    &  
                  bus(x(ilmo  ,1,indx_soil)), bus(x(hst  ,1,indx_soil)), &   
                  TRAD, N )
      if (phy_error_L) return
!
!
      IF(LSOIL_FREEZING_SVS1) THEN

              ! Optional activation of soil freezing

              CALL SOIL_FREEZING(DT, &
                bus(x(tpsoil   ,1,1)),bus(x(vegl    ,1,1)),&
                bus(x(vegh    ,1,1)), bus(x(psngrvl ,1,1)),&
                bus(x(psnvha  ,1,1)), bus(x(soilcondz,1,1)), &
                bus( x(soilhcapz,1,1)), &
                bus(x(tground, 1,1)), &
                bus(x(wsoil   ,1,1)) , bus(x(isoil   ,1,1)), &
                bus(x(snoro   ,1,1)) , bus(x(snodpl   ,1,1)), &
                bus(x(tsnow   ,1,2)) ,  &
                bus(x(snvro   ,1,1)) , bus(x(snvdp   ,1,1)), &
                bus(x(tsnowveg   ,1,2)) ,bus(x(tperm, 1,1)),   &
                bus(x(wunfrz, 1,1)), &
                N )
      ENDIF

      CALL HYDRO_SVS ( DT,      & 
           bus(x(eg      ,1,1)), bus(x(er      ,1,1)),&
           bus(x(etr     ,1,1)), rainrate_mm         ,&
           bus(x(rsnowsa ,1,1)), bus(x(rsnowsv ,1,1)),&
           bus(x(impervu ,1,1)), bus(x(vegl    ,1,1)),&
           bus(x(vegh    ,1,1)), bus(x(psngrvl ,1,1)),&
           bus(x(psnvha  ,1,1)), bus(x(acroot  ,1,1)),&
           wrmax,                bus(x(wsat    ,1,1)),&
           bus(x(ksat    ,1,1)), bus(x(psisat  ,1,1)),&
           bus(x(bcoef   ,1,1)), bus(x(fbcof   ,1,1)),&
           bus(x(wfcint  ,1,1)), bus(x(grkef   ,1,1)),&
           bus(x(snoma   ,1,1)), bus(x(snvma   ,1,1)),&
           bus(x(wveg    ,1,1)), wvegt               ,&
           bus(x(wsoil   ,1,1)), wsoilt              ,&
           bus(x(isoil   ,1,1)), isoilt              ,&
           bus(x(ksatc   ,1,1)), bus(x(khc     ,1,1)),&
           bus(x(psi     ,1,1)), bus(x(grksat  ,1,1)),&
           bus(x(wfcdp   ,1,1)), bus(x(watflow ,1,1)),&
           bus(x(latflw  ,1,1)), &
           bus(x(runofftot ,1,indx_soil)), N, bus(x(watpond ,1,1)), &
           bus(x(maxpond ,1,1)))



      IF( USE_PHOTO ) THEN

         if(vf_type == "CCILCECO") then            
            CALL PHTSYN_SVS_CCILCECO( BUS(x(LAIVF26,1,1))  , BUS(x(VEGF_EVOL   ,1,1)), &
                        BUS(x(TVEGE  ,1,1))  , ps, &
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

            CALL PHTSYN_SVS ( BUS(x(LAIVF26,1,1))  , BUS(x(VEGF,1,1)), &
                        BUS(x(TVEGE  ,1,1))  , ps, &
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
      CALL UPDATE_SVS ( WSOILT, ISOILT, WVEGT, &
                       bus(x(WSOIL   ,1,1)), bus(x(ISOIL   ,1,1)),  &
                       bus(x(WVEG    ,1,1)), bus(x(WSOILM  ,1,1)), &
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
    END subroutine svs
