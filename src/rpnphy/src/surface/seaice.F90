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
subroutine seaice3(BUS, BUSSIZ, PTSURF, PTSURFSIZ, lcl_indx, &
     N, M, NK)
   use tdpack
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use cpl_itf, only: cpl_update
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
   !@Object The multi-level model calculates the temperature profile across a
   !          snow/ice slab and the surface fluxes of heat, moisture, and momentum
   !          over floating ice-covered surfaces (sea/lake).
   !@Arguments
   !           - Input/Output -
   ! BUS       bus of surface variables
   !           - Input -
   ! BUSSIZ    size of the surface bus
   ! PTSURF    surface pointers
   ! PTSURFSIZ dimension of ptsurf
   ! lcl_indx
   ! N         running length
   ! M         horizontal dimension
   ! NK        vertical dimension

   integer BUSSIZ, N, M, NK
   integer PTSURFSIZ
   integer PTSURF(PTSURFSIZ), lcl_indx(2,n)
   real,target :: bus(bussiz)

   !@Author J. Mailhot (April 1999)
   !Revisions
   ! 001     M. Carrera and V. Fortin (Nov 2007) - Compute total runoff
   !                as precipitation - evaporation (no storage)
   !@Notes
   !          One-dimensional thermodynamic sea ice model:
   !          - based on modified version of nl-layer model of Semtner (1976, JPO 6, 379-389).
   !          - includes snow cover on the top of ice, heat conduction through snow and ice,
   !            thermal inertia of snow and ice layers, and penetrating solar radiation
   !          - includes a parameterization of albedo, conductivity and heat capacity
   !            (cf. Ebert and Curry 1993, JGR 98, 10085-10109; Flato and Brown 1996,
   !             JGR 101, 25767-25777)
   !          - effects of open leads in ice are considered (using climatological values
   !            of lead fraction to modify the analysed ice fraction)
   !          - vertical discretization similar to Flato and Brown (1996)
   !            but a flux boundary condition applied at the surface

   !          When the option ICEMELT = .TRUE. the model allows melting/growth of snow depth
   !          and ice thickness, and open water in summer with a crude oceanic mixed layer.

   !          The routine works with an arbitrary number of levels, but needs at least 2
   !          (NL .GE. 2).
   !*@/

   integer SURFLEN
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   integer, parameter :: INDX_SFC = INDX_ICE

   real, parameter :: zt_rho = 1.5      !Height used to compute air density
   logical         :: cplupd
   real, parameter :: delh_tresh=1.e-10 !Treshold to compute ctu in ocean (ice) coupled mode
   logical, parameter :: SEAICE_TDIAGLIM=.false.
   real delh, delq

   integer,dimension(n) :: nsnow

   real,dimension(n) :: dqsat, qsat, rhoa, scr1, scr2, scr3, scr4
   real,dimension(n) :: scr5,  scr6,  scr7, scr8, scr9, scr10,scr11
   real,dimension(n) :: tb,    vmod,  vdir, zsnodp_m,  vmod0
   real,dimension(n) :: my_ta,my_qa
   real,dimension(n) :: zu10, zusr   ! wind at 10m and at sensor level
   real,dimension(n) :: zref_sw_surf, zemit_lw_surf, zzenith
   real,dimension(n) :: zusurfzt, zvsurfzt, zqd

   real,dimension(n,nl) :: a, b, c, cap, cond, d, dz, sour, tp, z

   real,pointer,dimension(:) :: albsfc, cmu, ctu, fc_ice
   real,pointer,dimension(:) :: fsol, fv_ice, zemisr
   real,pointer,dimension(:) :: hice, hst_ice, hu, ilmo_ice
   real,pointer,dimension(:) :: ps, qsice, th, ts, tt, uu, vv
   real,pointer,dimension(:) :: z0h, z0m, zalfaq, zalfat
   real,pointer,dimension(:) :: zdlat, zfcor, zfdsi, zftemp, zfvap
   real,pointer,dimension(:) :: zml, zqdiag, zrainrate, zrunofftot
   real,pointer,dimension(:) :: zsnodp, zsnowrate, ztdiag
   real,pointer,dimension(:) :: ztsrad, zudiag, zvdiag, zfrv, zzusl, zztsl
   real,pointer,dimension(:) :: zfsd, zfsf, zcoszeni
   real,pointer,dimension(:) :: zutcisun, zutcishade, zwbgtsun, zwbgtshade
   real,pointer,dimension(:) :: zradsun, zradshade, ztglbsun, ztglbshade
   real,pointer,dimension(:) :: ztwetb, zq1, zq2, zq3, zq4, zq5, zq6, zq7
   real,pointer,dimension(:) :: zqdiagtyp, ztdiagtyp, zudiagtyp, zvdiagtyp
   real,pointer,dimension(:) :: zqdiagtypv, ztdiagtypv, zudiagtypv, zvdiagtypv

   real,pointer,dimension(:,:) :: t

   integer I, K
   real BETA1, SC

   real, save :: CON1,CON2,CON3,CON4,CON5,CON6
   real, save :: CON7,CON8,CON9,CON10,CON11
   real, save :: FI0,CONDFI,TFRZW,TMELI,TMELS
   real, save :: ALBOW,ALBDI,ALBMI,ALBDS,ALBMS,EMISW
   real, save :: COEFCOND,COEFHCAP,COEFEXT
   real, save :: ROICE,ROSNOW(2),ROSWTR
   real, save :: Z0W
   real, save :: Z0LAKEICE
   real, save :: BASEHF
   real, save :: HCAPI,HCAPW,VHFICE,VHFBAS,VHFSNO
   real, save :: HSMIN,DMIX
   real :: EMISNO,EMISI

   data   CON1     , CON2   , CON3  , CON4    / &
        2.845E-6 , 2.7E-4 , 233.0 , 0.2   /
   data   CON5     , CON6   , CON7  , CON8  / &
        92.88    , 7.364  , 3.2   , 14.24 /
   data   CON9     , CON10 , CON11 / &
        19.39    , 0.1   , 0.44  /
   data  TFRZW , TMELI , TMELS /  271.2 , 273.05 , 273.15  /
   data  ALBOW   ,  ALBDI  ,  ALBMI ,  ALBDS  ,  ALBMS / &
        0.08    ,  0.57   ,  0.50  ,  0.83   ,  0.77  /
   data  FI0   ,  CONDFI , COEFCOND , COEFHCAP , COEFEXT / &
        0.17  ,  2.034  , 0.1172   , 1.715E+7 , 1.5     /
   data  EMISW / 0.97 /
   data  ROICE ,  ROSWTR  / 913.0 ,  1025.0  /
   data  ROSNOW / 330.0 , 450.0 /
   data  Z0W  / 3.2E-5 /
   data  Z0LAKEICE / 1.6E-4 /
   data  BASEHF / 2.0 /
   data  HCAPI    , HCAPW    , VHFICE   , VHFBAS   , VHFSNO   / &
        2.062E+3 , 4.088E+3 , 3.014E+8 , 2.679E+8 , 1.097E+8 /
   data  HSMIN , DMIX / 0.010 , 30.0 /

   !------------------------------------------------------------------------
   SURFLEN = M

   albsfc   (1:n) => bus( x(alvis,1,indx_sfc) : )
   cmu      (1:n) => bus( x(bm,1,1)           : )
   ctu      (1:n) => bus( x(bt,1,indx_sfc)    : )
   fc_ice   (1:n) => bus( x(fc,1,indx_sfc)    : )
   fsol     (1:n) => bus( x(flusolis,1,1)     : )
   fv_ice   (1:n) => bus( x(fv,1,indx_sfc)    : )
   hice     (1:n) => bus( x(icedp,1,1)        : )
   hst_ice  (1:n) => bus( x(hst,1,indx_sfc)   : )
   ilmo_ice (1:n) => bus( x(ilmo,1,indx_sfc)  : )
   qsice    (1:n) => bus( x(qsurf,1,indx_sfc) : )
   ts       (1:n) => bus( x(tsurf,1,indx_sfc) : )
   z0h      (1:n) => bus( x(z0t,1,indx_sfc)   : )
   z0m      (1:n) => bus( x(z0,1,indx_sfc)    : )
   zalfaq   (1:n) => bus( x(alfaq,1,1)        : )
   zalfat   (1:n) => bus( x(alfat,1,1)        : )
   zcoszeni (1:n) => bus( x(cang,1,1)         : )
   zdlat    (1:n) => bus( x(dlat,1,1)         : )
   zemisr   (1:n) => bus( x(emisr,1,1)        : )
   zfcor    (1:n) => bus( x(fcor,1,1)         : )
   zfdsi    (1:n) => bus( x(fdsi,1,1)         : )
   zfrv     (1:n) => bus( x(frv,1,indx_sfc)   : )
   zftemp   (1:n) => bus( x(ftemp,1,indx_sfc) : )
   zfvap    (1:n) => bus( x(fvap,1,indx_sfc)  : )
   zml      (1:n) => bus( x(ml,1,1)           : )
   zrainrate(1:n) => bus( x(rainrate,1,1)     : )
   zrunofftot(1:n) => bus( x(runofftot,1,indx_sfc) : )
   zsnodp   (1:n) => bus( x(snodp,1,indx_sfc) : )
   zsnowrate(1:n) => bus( x(snowrate,1,1)     : )
   ztsrad   (1:n) => bus( x(tsrad,1,1)        : )
   zudiag   (1:n) => bus( x(udiag,1,1)        : )
   zudiagtyp(1:n) => bus( x(udiagtyp,1,indx_sfc) : )
   zudiagtypv(1:n) => bus( x(udiagtypv,1,indx_sfc) : )
   zvdiag   (1:n) => bus( x(vdiag,1,1)        : )
   zvdiagtyp(1:n) => bus( x(vdiagtyp,1,indx_sfc) : )
   zvdiagtypv(1:n) => bus( x(vdiagtypv,1,indx_sfc) : )
   ztdiag   (1:n) => bus( x(tdiag,1,1)        : )
   ztdiagtyp(1:n) => bus( x(tdiagtyp,1,indx_sfc) : )
   ztdiagtypv(1:n) => bus( x(tdiagtypv,1,indx_sfc) : )
   zqdiag   (1:n) => bus( x(qdiag,1,1)        : )
   zqdiagtyp(1:n) => bus( x(qdiagtyp,1,indx_sfc) : )
   zqdiagtypv(1:n) => bus( x(qdiagtypv,1,indx_sfc) : )
   zzusl    (1:n) => bus( x(zusl,1,1)         : )
   zztsl    (1:n) => bus( x(ztsl,1,1)         : )

   zfsd (1:n)       => bus( x(fsd,1,1)         : )
   zfsf (1:n)       => bus( x(fsf,1,1)         : )
   zutcisun (1:n)   => bus( x(yutcisun,1,indx_sfc)    : )
   zutcishade (1:n) => bus( x(yutcishade,1,indx_sfc)  : )
   zwbgtsun (1:n)   => bus( x(ywbgtsun,1,indx_sfc)    : )
   zwbgtshade (1:n) => bus( x(ywbgtshade,1,indx_sfc)  : )
   zradsun (1:n)    => bus( x(yradsun,1,indx_sfc)     : )
   zradshade (1:n)  => bus( x(yradshade,1,indx_sfc)   : )
   ztglbsun (1:n)   => bus( x(ytglbsun,1,indx_sfc)    : )
   ztglbshade (1:n) => bus( x(ytglbshade,1,indx_sfc)  : )
   ztwetb (1:n)     => bus( x(ytwetb,1,indx_sfc)      : )
   zq1 (1:n)        => bus( x(yq1,1,indx_sfc)         : )
   zq2 (1:n)        => bus( x(yq2,1,indx_sfc)         : )
   zq3 (1:n)        => bus( x(yq3,1,indx_sfc)         : )
   zq4 (1:n)        => bus( x(yq4,1,indx_sfc)         : )
   zq5 (1:n)        => bus( x(yq5,1,indx_sfc)         : )
   zq6 (1:n)        => bus( x(yq6,1,indx_sfc)         : )
   zq7 (1:n)        => bus( x(yq7,1,indx_sfc)         : )

   t(1:n,1:nl)    => bus( x(tmice,1,1)        : )

   if (atm_tplus) then
      hu       (1:n) => bus( x(huplus,1,nk)      : )
      ps       (1:n) => bus( x(pplus,1,1)        : )
      th       (1:n) => bus( x(thetaap,1,1)      : )
      tt       (1:n) => bus( x(tplus,1,nk)       : )
      uu       (1:n) => bus( x(uplus,1,nk)       : )
      vv       (1:n) => bus( x(vplus,1,nk)       : )
   else
      hu       (1:n) => bus( x(humoins,1,nk)     : )
      ps       (1:n) => bus( x(pmoins,1,1)       : )
      th       (1:n) => bus( x(thetaa,1,1)       : )
      tt       (1:n) => bus( x(tmoins,1,nk)      : )
      uu       (1:n) => bus( x(umoins,1,nk)      : )
      vv       (1:n) => bus( x(vmoins,1,nk)      : )
   endif
   
   ! Test on minimum value for number of layers
   if( NL .lt. 2 ) then
      call physeterror('seaice', 'NL < 2')
      return
   endif

   ! Fully-implicit time scheme
   BETA1 = 1.0
   SC = (1.0 - BETA1)*DELT

!*     1.     Preliminaries
!      --------------------

      i = sl_prelim(tt,hu,uu,vv,ps,zzusl,spd_air=vmod0,dir_air=vdir,rho_air=rhoa, &
           min_wind_speed=sqrt(vamin))
      if (i /= SL_OK) then
         call physeterror('seaice', 'error returned by sl_prelim()')
         return
      endif

      ! Set emissivity values
      select case (snow_emiss)
      case DEFAULT
         EMISNO = snow_emiss_const
      end select
      select case (ice_emiss)
      case DEFAULT
         EMISI = ice_emiss_const
      end select

      do I=1,N


!        Modifications for the sea ice surface:
!            - under-ice seawater temperature
         TB(I) = TFRZW

         ! Sea and lake ice roughness
         if ( ZML(I) .gt. 0.01 ) then
            ! non-zero lake fraction
            Z0M(I) = Z0LAKEICE
         else
            Z0M(I) = Z0SEAICE
         endif

         ! Ice thickness and icemelt
         if( ICEMELT ) then
!           Roughness lengths for the surface (ice/water)
            if( HICE(I) .lt. HIMIN ) then
               Z0M(I) = Z0W
!              remove snow if ice is too thin
               ZSNODP(I) = 0.0
            endif
         else
!         minimum ice thickness
            HICE(I) = max ( HICE(I) , HIMIN )
         endif

         Z0H(I) = Z0M(I)
!       Damp snow depth with h0*tanh(ff/h0)

         if (LIMSNODP)then
            ZSNODP_M(I) = SNOH0*tanh(ZSNODP(I)/SNOH0)
         else
            ZSNODP_M(I) = ZSNODP(I)
         endif

      end do

!       2.     Initialization of temperature profiles
!       ---------------------------------------------


!                               Temperature at mid-layers
      do K=1,NL-1
        do I=1,N
          TP(I,K) = 0.5*( T(I,K) + T(I,K+1) )
        end do
      end do

      do I=1,N
        TP(I,NL) = 0.5*( T(I,NL) + TB(I) )
      end do

!       3.     Calculate the drag and heat coefficients
!       -----------------------------------------------

      do I=1,N
!                               Saturated specific humidity at surface
        TS(I)    = T(I,1)
        QSAT(I)  = FOQST ( TS  (I), PS(I) )
        DQSAT(I) = FODQS ( QSAT(I), TS(I) )
        QSICE(I) = QSAT(I)
      end do

      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
           coefm=cmu,coeft=ctu,flux_t=zftemp,flux_q=zfvap,ilmo=ilmo_ice,ue=zfrv, &
           h=hst_ice,L_min=sl_Lmin_seaice,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('seaice', 'error returned by sl_sfclayer()')
         return
      endif



!       4.     Parameterizations (albedo, conductivity, capacity,...)
!       -------------------------------------------------------------


      do I=1,N
!                               - surface albedo (function of surface type and temperature)
        if( TS(I) .lt. TMELS ) then
          SCR2(I) = ALBDS
        else
          SCR2(I) = ALBMS
        endif

        if( TS(I) .lt. TMELI ) then
          SCR3(I) = min(ALBDI,max(ALBOW,0.08+CON11*HICE(I)**0.28))
          SCR9(I) = CHLC+CHLF
        else
          SCR3(I) = ALBMI
          SCR9(I) = CHLC
        endif

      end do

      do I=1,N
!                               - emissivity and fraction of penetrating solar radiation
!                                 (function of surface type snow/ice/water)
        if( HICE(I) .ge. HIMIN ) then
          if( ZSNODP(I) .gt. CON10 ) then
            ALBSFC(I) = SCR2(I)
            ZEMISR(I) = EMISNO
            SCR6(I) = 1.0
          else
            ALBSFC(I) = min( SCR2(I) , SCR3(I)+ZSNODP(I)* &
                            (SCR2(I)-SCR3(I))/CON10 )
            ZEMISR(I) = EMISI
            SCR6(I) = 1.0 - FI0
            SCR6(I) = min( 1.0 , SCR6(I)+ZSNODP(I)* &
                            (1.0-SCR6(I))/CON10 )
          endif
        else
          ALBSFC(I) = ALBOW
          ZEMISR(I) = EMISW
          SCR6(I) = 1.0
          SCR9(I) = CHLC
        endif
!                               - bulk salinity of the ice
        SCR8(I) = max( CON7 , CON8-CON9*HICE(I) )

      end do


      do I=1,N
       if( HICE(I) .ge. HIMIN ) then
!                                                       snow layer
        if( ZSNODP(I) .ge. HSMIN ) then
          NSNOW(I) = 1
          DZ(I,1) = ZSNODP_M(I)
          SCR4(I) = ROSNOW(1)
          if( TP(I,1) .ge. TMELS ) SCR4(I) = ROSNOW(2)
          COND(I,1) = CON1*SCR4(I)**2+CON2*2.**((TP(I,1)-CON3)*CON4)
          CAP(I,1) = SCR4(I)*(CON5+CON6*TP(I,1))
        else
          NSNOW(I) = 0
          DZ(I,1) = HICE(I)/FLOAT(NL)
        endif

       else
         NSNOW(I) = 0
         DZ(I,1) = DMIX
         CAP(I,1) = DMIX*ROSWTR*HCAPW
       endif
      end do

!                                                      ice slab
      do K=2,NL
        do I=1,N
          if( HICE(I) .ge. HIMIN ) then
            SCR3(I) = ROICE*HCAPI
            DZ(I,K) = HICE(I)/FLOAT(NL-NSNOW(I))
            COND(I,K) = CONDFI+COEFCOND*SCR8(I)/min(TP(I,K)-TMELI,-0.1)
            COND(I,K) = max( 0.2*CONDFI , COND(I,K) )

            CAP(I,K) = SCR3(I)+COEFHCAP*SCR8(I)/ &
                     min( T(I,K)-TMELI , -0.1 )**2
            CAP(I,K) = min( 100.0*SCR3(I) , CAP(I,K) )
          else
            NSNOW(I) = 0
            DZ(I,K) = DMIX
            CAP(I,K) = DMIX*ROSWTR*HCAPW
          endif
        end do
      end do

      do I=1,N
        if( HICE(I) .ge. HIMIN .and. NSNOW(I).eq.0 ) then
          SCR3(I) = ROICE*HCAPI
          DZ(I,1) = HICE(I)/FLOAT(NL-NSNOW(I))
          COND(I,1) = CONDFI+COEFCOND*SCR8(I)/min(TP(I,1)-TMELI,-0.1)
          COND(I,1) = max( 0.2*CONDFI , COND(I,1) )

          CAP(I,1) = SCR3(I)+COEFHCAP*SCR8(I)/ &
                   min( T(I,1)-TMELI , -0.1 )**2
          CAP(I,1) = min( 100.0*SCR3(I) , CAP(I,1) )
        endif
      end do

      do I=1,N
       if( HICE(I) .ge. HIMIN ) then
        if( NSNOW(I) .eq. 0 ) then
          CAP(I,1) = SCR3(I)+COEFHCAP*SCR8(I)/ &
                     min( 0.75*T(I,1)+0.25*T(I,2)-TMELI , -0.1 )**2
          CAP(I,1) = min( 100.0*SCR3(I) , CAP(I,1) )
        else
          if( NL .gt. 2) then
            CAP(I,2) = SCR3(I)+COEFHCAP*SCR8(I)/ &
                     min( 0.75*T(I,2)+0.25*T(I,3)-TMELI , -0.1 )**2
          else
            CAP(I,2) = SCR3(I)+COEFHCAP*SCR8(I)/ &
                     min( 0.75*T(I,2)+0.25*TB(I)-TMELI , -0.1 )**2
          endif
          CAP(I,2) = min( 100.0*SCR3(I) , CAP(I,2) )
        endif
!                                              add "thin" snow layer
        if( NSNOW(I) .eq. 0 ) then
          SCR4(I) = ROSNOW(1)
          if( T(I,1) .ge. TMELS ) SCR4(I) = ROSNOW(2)
          SCR4(I) = CON1*SCR4(I)**2+CON2*2.**((T(I,1)-CON3)*CON4)
          COND(I,1) = COND(I,1)/ &
                      (1.0+ZSNODP(I)*COND(I,1)/(DZ(I,1)*SCR4(I)))
        endif

       endif
      end do

      do I=1,N
!                                              penetrating solar radiation
        SCR7(I) = (1.0-SCR6(I))*FSOL(I)*(1.-ALBSFC(I))
        if( NSNOW(I).eq.1 ) then
          Z(I,1) = 0.0
          Z(I,2) = 0.5*DZ(I,2)
        else
          Z(I,1) = 0.5*DZ(I,1)
          Z(I,2) = Z(I,1)+0.5*(DZ(I,2)+DZ(I,1))
        endif
        SOUR(I,1) = SCR7(I)*(1.-exp(-COEFEXT*Z(I,1)))
        SOUR(I,2) = SCR7(I)*(exp(-COEFEXT*Z(I,1)) &
                            -exp(-COEFEXT*Z(I,2)))

      end do

      if( NL .gt. 2) then
        do K=3,NL
          do I=1,N
            Z(I,K) = Z(I,K-1)+0.5*(DZ(I,K)+DZ(I,K-1))
            SOUR(I,K) = SCR7(I)*(exp(-COEFEXT*Z(I,K-1)) &
                                -exp(-COEFEXT*Z(I,K)))
          end do
        end do
        do I=1,N
          Z(I,NL) = Z(I,NL)+0.5*DZ(I,NL)
          SOUR(I,NL) = SCR7(I)*(exp(-COEFEXT*Z(I,NL-1)) &
                              -exp(-COEFEXT*Z(I,NL)))
        end do
      else
        do I=1,N
          Z(I,NL) = Z(I,NL)+0.5*DZ(I,NL)
          SOUR(I,NL) = SCR7(I)*(1.0-exp(-COEFEXT*Z(I,NL)))
        end do
      endif



!       5.     Compute the temperature profile
!       --------------------------------------

!                                              linearized terms in surface heat budget
      do I=1,N

        SCR1(I) = 4. * ZEMISR(I) * STEFAN * TS(I)**3 &
           +  RHOA(I) * CTU(I) * (DQSAT(I) * SCR9(I) + CPD)


        SCR2(I) = 3. * ZEMISR(I) * STEFAN  * TS(I)**4 &
           + RHOA(I) * CTU(I) * DQSAT(I) * SCR9(I) * TS(I)

        SCR3(I) = RHOA(I)*CTU(I)*(CPD*TH(I)-SCR9(I)*(QSAT(I)-HU(I))) &
           + SCR6(I)*FSOL(I)*(1.-ALBSFC(I)) &
           + ZEMISR(I)*ZFDSI(I)

      end do

!                                              setup tridiagonal terms A, B, C
!                                              and right-hand-side term D

      if( NL.gt.3 ) then
        do K=3,NL-1
          do I=1,N
            if( HICE(I) .ge. HIMIN ) &
              CAP(I,K) = 0.5*( DZ(I,K)+DZ(I,K-1) )*CAP(I,K)
          end do
        end do
      endif

      do I=1,N
        if( HICE(I) .ge. HIMIN ) then
          SOUR(I,2) = SOUR(I,2)+SOUR(I,1)
          if( NSNOW(I) .eq. 0 ) then
            CAP(I,2) = ( 0.5*DZ(I,2)+DZ(I,1) )*CAP(I,2)
!                                              add "thin" snow layer
            SCR4(I) = ROSNOW(1)
            if( T(I,1) .ge. TMELS ) SCR4(I) = ROSNOW(2)
            SCR4(I) = SCR4(I)*(CON5+CON6*T(I,1))
            CAP(I,2) = CAP(I,2) + ZSNODP(I)*SCR4(I)
            CAP(I,NL) = 0.5*( DZ(I,NL-1)+DZ(I,NL) )*CAP(I,NL)
          else
            CAP(I,2) = 0.5*DZ(I,2)*CAP(I,2)+DZ(I,1)*CAP(I,1)
            CAP(I,NL) = 0.5*( DZ(I,NL-1)+DZ(I,NL) )*CAP(I,NL)
          endif
        endif
      end do

      do K=2,NL
        do I=1,N
          if( HICE(I) .ge. HIMIN ) then
            A(I,K) = (COND(I,K-1)/DZ(I,K-1))/CAP(I,K)
            C(I,K) = (COND(I,K)/DZ(I,K))/CAP(I,K)
            B(I,K) = -A(I,K)-C(I,K)
            D(I,K) = DELT*SOUR(I,K)/CAP(I,K) + T(I,K)
          endif
        end do
      end do

      do I=1,N
        if( HICE(I) .ge. HIMIN ) then
          A(I,1) = 0.0
          B(I,1) = 0.0
          C(I,1) = 0.0
          D(I,1) = 0.0
          C(I,NL) = 0.0
        endif
      end do

      do K=1,NL
        do I=1,N
          if( HICE(I) .lt. HIMIN ) then
            A(I,K) = 0.0
            B(I,K) = 0.0
            C(I,K) = 0.0
            D(I,K) = 0.0
          endif
        end do
      end do

      do I=1,N
          if( HICE(I) .lt. HIMIN ) &
                 D(I,1) = T(I,1)
      end do



!                                             back to full levels
      do K=1,NL
        do I=1,N
          TP(I,K) = T(I,K)
        end do
      end do

      call DIFUVD1(D, SC, A, B, C, TP, D, N, N, NL)

      do K=1,NL
        do I=1,N
          A(I,K) =    -BETA1*DELT*A(I,K)
          B(I,K) = 1.0-BETA1*DELT*B(I,K)
          C(I,K) =    -BETA1*DELT*C(I,K)
        end do
      end do

      do I=1,N
!                                              add upper boundary condition
        if( HICE(I) .ge. HIMIN ) then
          C(I,1) = -BETA1
          B(I,1) = -C(I,1)+DZ(I,1)*SCR1(I)/COND(I,1)
          D(I,1) = DZ(I,1)*(SCR2(I)+SCR3(I))/COND(I,1) &
                   +(1.0-BETA1)*(T(I,2)-T(I,1))
!                                              add lower boundary condition
          D(I,NL) = D(I,NL)+DELT*TB(I)*(COND(I,NL)/ &
                            DZ(I,NL))/CAP(I,NL)
        endif
      end do

      do K=1,NL
        do I=1,N
          if( HICE(I) .lt. HIMIN ) then
            A(I,K) = -1.0
            B(I,K) =  1.0
          endif
        end do
      end do

      do I=1,N
        if( HICE(I) .lt. HIMIN ) then
          A(I,1) =  0.0
          B(I,1) =  1.0
          C(I,1) =  0.0
          B(I,1) = B(I,1)+DELT*SCR1(I)/CAP(I,1)
          D(I,1) = D(I,1)+DELT*(SCR2(I)+SCR3(I))/CAP(I,1)
        endif
      end do


      call DIFUVD2(TP, A, B, C, D, D, N, N, NL)

!                                              prevent temperatures
!                                              from exceeding melting temperature
      do I=1,N
        if( HICE(I) .ge. HIMIN ) then
          if( ZSNODP(I) .ge. HSMIN ) then
            SCR4(I) = TMELS
          else
            SCR4(I) = TMELI
          endif
        endif
      end do

      do K=1,NL
        do I=1,N
          if( HICE(I) .ge. HIMIN ) &
            TP(I,K) = min ( TP(I,K) , SCR4(I) )
        end do
      end do

      do I=1,N
        if( HICE(I) .lt. HIMIN ) &
          SCR4(I) = TB(I)
      end do


!                                              surface and snow/ice temperatures
      do K=1,NL
        do I=1,N
          T(I,K) = TP(I,K)
        end do
      end do

      do I=1,N
        TS(I) = TP(I,1)
!                                              saturated specific humidity at surface
        QSICE(I) = FOQST( TS(I), PS(I) )
      end do


!       6.     Melting and growth
!       -------------------------
!                                              melting at the surface (snow/ice)
!                                              growth/melting of ice at the lower boundary
      if(ICEMELT) then

        do I=1,N

!                               recompute surface albedo with new TS

          if( TS(I) .lt. TMELS ) then
            SCR2(I) = ALBDS
          else
            SCR2(I) = ALBMS
          endif

          if( TS(I) .lt. TMELI ) then
            SCR3(I) = min(ALBDI,max(ALBOW,0.08+CON11*HICE(I)**0.28))
            SCR9(I) = CHLC+CHLF
          else
            SCR3(I) = ALBMI
            SCR9(I) = CHLC
          endif

          if( HICE(I) .ge. HIMIN ) then
            if( ZSNODP(I) .gt. CON10 ) then
              ALBSFC(I) = SCR2(I)
            else
              ALBSFC(I) = min( SCR2(I) , SCR3(I)+ZSNODP(I)* &
                              (SCR2(I)-SCR3(I))/CON10 )
            endif
          else
            ALBSFC(I) = ALBOW
            SCR9(I) = CHLC
          endif

        end do

        do I=1,N
!                               compute heat conduction flux in upper layer
          if( HICE(I) .ge. HIMIN ) then
            if( ZSNODP(I) .ge. HSMIN ) then
              SCR7(I) = ROSNOW(1)
              SCR5(I) = 0.5*(T(I,1)+T(I,2))
              if( SCR5(I) .ge. TMELS ) SCR7(I) = ROSNOW(2)
              COND(I,1) = CON1*SCR7(I)**2 &
                        +CON2*2.**((SCR5(I)-CON3)*CON4)
            else
              SCR7(I) = 0.5*(T(I,1)+T(I,2))
              COND(I,1) = CONDFI+COEFCOND*SCR8(I) &
                          /min(SCR7(I)-TMELI,-0.1)
              COND(I,1) = max( 0.2*CONDFI , COND(I,1) )
!                                              add "thin" snow layer
              SCR5(I) = ROSNOW(1)
              if( T(I,1) .ge. TMELS ) SCR5(I) = ROSNOW(2)
              SCR5(I) = CON1*SCR5(I)**2+CON2*2.**((T(I,1)-CON3)*CON4)
              COND(I,1) = COND(I,1)/ &
                      (1.0+ZSNODP(I)*COND(I,1)/(DZ(I,1)*SCR5(I)))
            endif
            SCR11(I) = COND(I,1)*(T(I,2)-T(I,1))/DZ(I,1)
          else
            SCR11(I) = 0.0
          endif

!                               compute HA = H + LE + FWnet + FLnet
          SCR5(I) = RHOA(I)*CTU(I)* &
                    ( CPD*(TS(I)-TH(I))+SCR9(I)*(QSICE(I)-HU(I)) ) &
           - SCR6(I)*FSOL(I)*(1.-ALBSFC(I)) &
                   + ZEMISR(I)*( STEFAN*TS(I)**4 - ZFDSI(I) )

          SCR3(I) = SCR5(I) - SCR11(I)

        end do


        do I=1,N
!                               criteria for melting at the surface
          SCR10(I) = 0.0
          if( TS(I).ge.SCR4(I) .and. &
              SCR5(I).le.0.0 .and. SCR3(I).lt.0.0 ) then
!                               melt available snow...
            ZSNODP(I) = ZSNODP(I) + DELT*SCR3(I)/VHFSNO
!                               ...then ice...
            if( ZSNODP(I).lt.0.0 ) then
              HICE(I) = HICE(I) + ZSNODP(I)*VHFSNO/VHFICE
              ZSNODP(I) = 0.0
!                               ...then warm oceanic mixed layer
              if( HICE(I).lt.0.0 ) then
                TS(I) = TS(I) - HICE(I)*VHFICE/CAP(I,1)
                HICE(I) = 0.0
                SCR10(I) = 1.0
              endif
            endif

          endif
        end do
!                               add growth/melt of ice at lower boundary
        do I=1,N
          if( HICE(I) .ge. HIMIN ) then

!                               compute heat conduction flux in lower layer
            SCR11(I) = COND(I,NL)*(TB(I)-T(I,NL))/DZ(I,NL)
!                               compute absorption of penetrating solar radiation
            SCR8(I) = FSOL(I)*(1.-ALBSFC(I))*(1.-SCR6(I)) &
                      *exp(-COEFEXT*Z(I,NL))
            HICE(I) = max ( 0.0 , HICE(I) + &
                      DELT*( SCR11(I)-SCR8(I)-BASEHF )/VHFBAS )
!                                if ice is too thin, switch to oceanic mixed layer
            if( HICE(I).lt.HIMIN ) then
              TS(I) = TFRZW
              ZSNODP(I) = 0.0
              SCR10 (I) = 1.0
            endif

          else
            if( SCR5(I).gt.0.0 ) then
!                               cool oceanic mixed layer...
              TS(I) = TS(I) - DELT*SCR5(I)/CAP(I,1)
!                               ...then form new ice
              if( TS(I).le.TFRZW ) then
                HICE(I) = HICE(I) + (TFRZW-TS(I))*CAP(I,1)/VHFBAS
                TS(I) = TFRZW
              endif
              SCR10(I) = 1.0
            endif

          endif
        end do
!                               ...and adjust final temperature profile

        do K=1,NL
          do I=1,N
            if( SCR10(I).eq.1.0 ) then
              T(I,K) = TS(I)
              TP(I,K) = TS(I)
            endif
          end do
        end do
!                               add snow fall (units changed from
!                               water equivalent to snow equivalent)
        do I=1,N
          if( TS(I).lt.TMELI .and. HICE(I).gt.HIMIN ) then
            ZSNODP(I) = ZSNODP(I) + DELT*ZSNOWRATE(I)*(1000./ROSNOW(1))
          endif
        end do

!                               transition periods
        do I=1,N
          if( NSNOW(I).eq.0 .and. ZSNODP(I).ge.HSMIN ) then
            T(I,1) = TP(I,1)
            T(I,2) = TP(I,1)
            Z(I,1) = 0.0
          endif
        end do
        do K=2,NL
          do I=1,N
            if( NSNOW(I).eq.0 .and. ZSNODP(I).ge.HSMIN ) &
              Z(I,K) = Z(I,K-1)+DZ(I,K-1)
          end do
        end do

        do I=1,N
          if( NSNOW(I).eq.0 .and. ZSNODP(I).ge.HSMIN ) then
            DZ(I,1) = 0.0
            DZ(I,2) = 0.0
          endif
        end do

        if( NL.gt.2 ) then
          do K=3,NL
            do I=1,N
              if( NSNOW(I).eq.0 .and. ZSNODP(I).ge.HSMIN ) then
                DZ(I,K) = DZ(I,K-1)+HICE(I)/FLOAT(NL-1)
                T(I,K) = TP(I,K-1)+((DZ(I,K)-Z(I,K-1)) &
                         /(Z(I,K)-Z(I,K-1)))*(TP(I,K)-TP(I,K-1))
              endif
            end do
          end do
        endif

        do I=1,N
          if( NSNOW(I).eq.1 .and. ZSNODP(I).lt.HSMIN ) then
            T(I,1) = TP(I,2)
            Z(I,1) = 0.0
            Z(I,2) = 0.0
          endif
        end do
        if( NL.gt.2 ) then
          do K=3,NL
            do I=1,N
              if( NSNOW(I).eq.1 .and. ZSNODP(I).lt.HSMIN ) &
                Z(I,K) = Z(I,K-1)+DZ(I,K-1)
            end do
          end do
        endif
        do I=1,N
          if( NSNOW(I).eq.1 .and. ZSNODP(I).lt.HSMIN ) then
            DZ(I,1) = 0.0
            DZ(I,2) = HICE(I)/FLOAT(NL)
            T(I,2) = TP(I,2)+((DZ(I,2)-Z(I,2)) &
                     /(Z(I,3)-Z(I,2)))*(TP(I,3)-TP(I,2))
          endif
        end do
        if( NL.gt.3 ) then
          do K=3,NL-1
            do I=1,N
              if( NSNOW(I).eq.1 .and. ZSNODP(I).lt.HSMIN ) then
                DZ(I,K) = DZ(I,K-1)+HICE(I)/FLOAT(NL)
                T(I,K) = TP(I,K)+((DZ(I,K)-Z(I,K)) &
                       /(Z(I,K+1)-Z(I,K)))*(TP(I,K+1)-TP(I,K))
              endif
            end do
          end do
        endif
        do I=1,N
          if( NSNOW(I).eq.1 .and. ZSNODP(I).lt.HSMIN ) then
            DZ(I,NL) = DZ(I,NL-1)+HICE(I)/FLOAT(NL)
            T(I,NL) = TP(I,NL)+((DZ(I,NL)-Z(I,NL)) &
                  /(max(HICE(I),HIMIN)/FLOAT(NL-1)))*(TB(I)-TP(I,NL))
          endif
        end do

      endif

!       7.     Update fluxes
!       ------------------------------------------------

!     Update TS and QSICE
      do I=1,N
        TS   (I) = T(I,1)
        QSICE(I) = FOQST( TS(I), PS(I) )
      end do

      ! Compute diagnostic quantities at 1.5m
      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
           hghtt_diag=zt_rho,t_diag=my_ta,q_diag=my_qa,tdiaglim=SEAICE_TDIAGLIM, &
           L_min=sl_Lmin_seaice,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('seaice', 'error 2 returned by sl_sfclayer()')
         return
      endif

      do I=1,N
!                               remove snow if ice is too thin
        if( HICE(I).lt.HIMIN ) ZSNODP(I) = 0.0

        ZTSRAD   (I) = TS(I)

        ZALFAT   (I) = - CTU(I) * ( TS    (I) - TH(I) )
        ZALFAQ   (I) = - CTU(I) * ( QSICE (I) - HU(I) )
        if (.not.IMPFLX) CTU (I) = 0.
        RHOA     (I) = PS(I)/(RGASD * my_ta(I)*(1.+DELTA*my_qa(I)))
        ZRUNOFFTOT(I) = (1000.*(ZRAINRATE(I) +ZSNOWRATE(I)) &
             + RHOA(I)*ZALFAQ(I)) * DELT
        FC_ICE(I)    = -CPD *RHOA(I)*ZALFAT(I)
        FV_ICE(I)    = -(CHLC+CHLF)*RHOA(I)*ZALFAQ(I)

        if (IMPFLX) then
          ZALFAT   (I) = - CTU(I) *  TS(I)
          ZALFAQ   (I) = - CTU(I) *  QSICE(I)
        endif

      end do

      ! Estimate diagnostic quantities at user-specified level
      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
           hghtm_diag=zu,hghtt_diag=zt,t_diag=ztdiag,q_diag=zqdiag,u_diag=zudiag, &
           v_diag=zvdiag,tdiaglim=SEAICE_TDIAGLIM,L_min=sl_Lmin_seaice,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('seaice', 'error 3 returned by sl_sfclayer()')
         return
      endif
      if (sl_Lmin_seaice > 0.) then
         zudiag = zudiag * vmod0 / vmod
         zvdiag = zvdiag * vmod0 / vmod
      endif
      
      ! Fill surface type-specific diagnostic values
      zqdiagtyp = zqdiag
      ztdiagtyp = ztdiag
      zudiagtyp = zudiag
      zvdiagtyp = zvdiag
      zqdiagtypv = zqdiag
      ztdiagtypv = ztdiag
      zudiagtypv = zudiag
      zvdiagtypv = zvdiag

      if (cplocn) then
         ! Update with fluxes and diagnostic variables from ocean model
         cplupd=.false.
         call cpl_update (vmod(1:n), 'UVI', lcl_indx, n, u=UU, v=VV, cplu=cplupd)
         if (cplupd) then
           call cpl_update (FC_ICE(1:n), 'SHI' , lcl_indx, n)
           call cpl_update (FV_ICE(1:n), 'LHI' , lcl_indx, n)
           call cpl_update (ZUDIAG(1:n), 'ZUI' , lcl_indx, n)
           call cpl_update (ZVDIAG(1:n), 'ZVI' , lcl_indx, n)
           call cpl_update (ZTDIAG(1:n), 'ZTI' , lcl_indx, n)
           call cpl_update (ZQDIAG(1:n), 'ZQI' , lcl_indx, n)
           call cpl_update (TS    (1:n), 'I7I' , lcl_indx, n)
           call cpl_update (ZTSRAD(1:n),   'T4I' , lcl_indx, n)
           call cpl_update (QSICE (1:n),   'QSI' , lcl_indx, n)
           call cpl_update (ILMO_ICE(1:n), 'ILI' , lcl_indx, n)
           call cpl_update (Z0M   (1:n),   'ZMI' , lcl_indx, n)
           call cpl_update (Z0H   (1:n),   'ZHI' , lcl_indx, n)
           
           ! Derives other variables from cpl_update output
           ! assuming FC_ICE and FV_ICE were computed with
           ! RHOA at the same level
           if (zt_rho == zt) then
             my_ta(1:n)=ztdiag(1:n)
             my_qa(1:n)=zqdiag(1:n)
           else
             call physeterror('seaice', 'inconsistent density level')
             return
           endif
           do I=1,N
              RHOA  (I)= PS(I)/(RGASD * my_ta(I)*(1.+DELTA*my_qa(I)))
              ZALFAT(I)= FC_ICE(I)/(-CPD *RHOA(I))
              ZALFAQ(I)= FV_ICE(I)/(-(CHLC+CHLF)*RHOA(I))
              ZFTEMP(I) = -ZALFAT(I)
              ZFVAP(I)  = -ZALFAQ(I)
              T   (I,1)= TS(I)
              if (IMPFLX) then
                ! CTU consistent with ice-ocean model fluxes (as possible)
                ! Uncertainties in CTU from interpolation/agregation
                ! will be compensated by ZALFAT and ZALFAQ
                delh=TS(I)-TH(I) ; delq=QSICE(I)-HU(I)
                CTU(I) = 0.5*( -ZALFAT(I)/sign(max(abs(delh),delh_tresh),delh) &
                               -ZALFAQ(I)/sign(max(abs(delq),delh_tresh),delq) )
                CTU(I) = max(0.,CTU(I)) ! CTU<0 should never occur except under vicinity
                                        ! of null fluxes agregation/interpolation anyway
                                        ! e.g. sign(ZALFAT) /= sign(delh) or sign(ZALFAQ) /= sign(delq)
                ! ZALFAT and ZALFAQ consistent with FC_ICE and FV_ICE
                ZALFAT(I) = ZALFAT(I) - CTU(I)*TH(I)
                ZALFAQ(I) = ZALFAQ(I) - CTU(I)*HU(I)
              else
                CTU(I) = 0. !Redundant but less dangerous
              endif
           end do
           ! Diagnostic ustar based on agregated ice model stress (tau=rho*ustar**2)
           call cpl_update (ZFRV  (1:n), 'FRI' , lcl_indx, n, rho=RHOA, vmod=VMOD, cmu=CMU)
           ! Updated consistent CMU needed by implicit scheme using the following relation
           ! ==> CM=ustar/vmod (ustar**2=tau/rho=CM*CM*vmod**2)
         endif
      endif

      !--------------------------------------
      !   8.     Heat Stress Indices
      !------------------------------------
      !#TODO: at least 4 times identical code in surface... separeted s/r to call
      IF_TERMAL_STRESS: if (thermal_stress) then

      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
           hghtm_diag=zt,hghtt_diag=zt,u_diag=zusurfzt, &
           v_diag=zvsurfzt,tdiaglim=SEAICE_TDIAGLIM,L_min=sl_Lmin_seaice,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('seaice', 'error 3 returned by sl_sfclayer()')
         return
      endif
      if (sl_Lmin_seaice > 0.) then
         zusurfzt = zusurfzt * vmod0 / vmod
         zvsurfzt = zvsurfzt * vmod0 / vmod
      endif
               
      do I=1,N

         if (abs(zzusl(i)-zu) <= 2.0) then
            zu10(i) = sqrt(uu(i)**2+vv(i)**2)
         else
            zu10(i) = sqrt(zudiag(i)**2+zvdiag(i)**2)
         endif

         ! wind  at SensoR level zubos at z=zt
         if( (abs(zusurfzt(i)) >= 0.1) .and. (abs(zvsurfzt(i)) >= 0.1)) then
         zusr(i) = sqrt( zusurfzt(i)**2 + zvsurfzt(i)**2)
         else
         zusr(i) = zu10(i)
         endif

             zqd(i) = max( ZQDIAG(i) , 1.e-6)

            zref_sw_surf(i) = albsfc(i) * fsol(i)
            zemit_lw_surf(i) = (1. -zemisr(i)) * zfdsi(i) + zemisr(i)*stefan*ztsrad(i)**4

         zzenith(i) = acos(zcoszeni(i))
         if (fsol(i) > 0.0) then
            zzenith(i) = min(zzenith(i), pi/2.)
         else
            zzenith(i) = max(zzenith(i), pi/2.)
         endif

      end do

         call SURF_THERMAL_STRESS(ZTDIAG, zqd,            &
              ZU10,ZUSR,  ps,                             &
              ZFSD, ZFSF, ZFDSI, ZZENITH,                 &
              ZREF_SW_SURF,ZEMIT_LW_SURF,                 &
              Zutcisun ,Zutcishade,                       &
              zwbgtsun, zwbgtshade,                       &
              zradsun, zradshade,                         &
              ztglbsun, ztglbshade, ztwetb,               &
              ZQ1, ZQ2, ZQ3, ZQ4, ZQ5,                    &
              ZQ6,ZQ7, N)
      endif IF_TERMAL_STRESS

!     FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
      call FILLAGG ( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_ICE, &
                     SURFLEN )

   return
end subroutine seaice3
