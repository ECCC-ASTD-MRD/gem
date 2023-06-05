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
subroutine glaciers2(BUS, BUSSIZ, PTSURF, PTSURFSIZ, N, M, NK)
   use tdpack
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
   !@Object Calculate the surface temperature (and specific humidity) and
   !          the surface fluxes of heat, moisture, and momentum over
   !          continental glaciers and ice sheets.
   !@Arguments
   !               - Input/Output -
   ! BUS           bus of surface variables
   !               - Input -
   ! BUSSIZ        size of the surface bus
   ! PTSURF        surface pointers
   ! PTSURFSIZ     dimension of ptsurf
   ! N             running length
   ! M             horizontal dimension
   ! NK            vertical dimension

   integer BUSSIZ, N, M, NK
   integer PTSURFSIZ
   integer PTSURF(PTSURFSIZ)
   real,target :: bus(bussiz)

   !@Author J. Mailhot (April 1999)
   !@Notes
   !          One-dimensional thermodynamic model to treat glaciers and ice sheets:
   !          - includes snow cover on the top of ice, heat conduction through snow and ice,
   !            thermal inertia of snow and ice layers, and penetrating solar radiation
   !          - includes a parameterization of albedo, conductivity and heat capacity
   !
   ! Note:     - this subroutine expects snow depth in metre 
   !*@/
   !Revisions
   ! 001    M. Carrera and V. Fortin (Nov 2007) - Compute total runoff

   integer SURFLEN
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   integer, parameter :: INDX_SFC = INDX_GLACIER
   real, parameter :: ZT_RHO = 1.5              !Height used to compute air density (m)

   real, save :: CON10
   real, save :: FI0,CONDFI,TMELICE,TMELSNO
   real, save :: ALBDI,ALBMI,ALBDS,ALBMS,COEFEXT
   real, save :: ROICE,ROSNOW(2)
   real, save :: HCAPI,HFICE,VHFSNO
   real, save :: HCAPS,KSDS,KSMS
   real :: EMISICE, EMISNOW

   real,dimension(n) :: a,    b,     c,     c2,   ct,   dqsat, rhoa
   real,dimension(n) :: scr1, scr2,  scr3,  scr4, scr5, scr6, scr7, scr8
   real,dimension(n) :: scr10, scr11, t2, vmod, vdir, zsnodp_m, vmod0
   real, dimension(n) :: zu10, zusr    ! wind at 10m and at the sensor level
   real, dimension(n) :: zref_sw_surf, zemit_lw_surf, zzenith
   real, dimension(n) :: my_ta, my_qa
   real, dimension(n) :: zusurfzt, zvsurfzt, zqd

   real,pointer,dimension(:) :: albsfc, cmu, ctu, fc_glac
   real,pointer,dimension(:) :: fsol, fv_glac, zemisr
   real,pointer,dimension(:) :: hst_glac, hu, ilmo_glac
   real,pointer,dimension(:) :: ps, qsice, th, snorate, tdeep, ts, tt, uu, vv
   real,pointer,dimension(:) :: z0h, z0m
   real,pointer,dimension(:) :: zalfaq, zalfat, zdlat, zfcor, zfdsi
   real,pointer,dimension(:) :: zftemp, zfvap, zqdiag, zrunofftot, zsnodp, ztdiag
   real,pointer,dimension(:) :: ztsurf, ztsrad, zudiag, zvdiag
   real,pointer,dimension(:) :: zfrv, zzusl, zztsl
   real,pointer,dimension(:) :: zqdiagtyp, ztdiagtyp, zudiagtyp, zvdiagtyp
   real,pointer,dimension(:) :: zqdiagtypv, ztdiagtypv, zudiagtypv, zvdiagtypv
   real,pointer,dimension(:) :: zfsd, zfsf, zcoszeni
   real,pointer,dimension(:) :: zutcisun, zutcishade, zwbgtsun, zwbgtshade
   real,pointer,dimension(:) :: zradsun, zradshade, ztglbsun, ztglbshade
   real,pointer,dimension(:) :: ztwetb, zq1, zq2, zq3, zq4, zq5, zq6, zq7

   logical, parameter :: GLACIER_TDIAGLIM=.false.

   integer I
   real DAY, DEEPFRAC

   data  CON10 / &
         0.1 /
   data  FI0   ,  CONDFI   / &
         0.17  ,  2.034    /
   data  TMELICE , TMELSNO  /  273.05 , 273.15  /
   data  ALBDI   ,  ALBMI ,  ALBDS  ,  ALBMS / &
         0.57    ,  0.50  ,  0.83   ,  0.77  /
   data  COEFEXT / &
         1.5     /
   data  ROICE  / &
         913.0  /
   data  ROSNOW  / 330.0 , 450.0 /
   data  HCAPI    , HFICE   , VHFSNO   / &
         2.062E+3 , 3.34E+5 , 1.097E+8 /
   data  HCAPS   , KSDS  , KSMS   / &
         2.04E+3 , 0.325 , 0.665 /
 
   !------------------------------------------------------------------------
   SURFLEN = M
   DAY = 86400.
   !# Fraction of diurnal temperature cycle
   !# (determines depth of deep in-glacier temperature)
   DEEPFRAC = exp(-1.0)

   albsfc   (1:n) => bus( x(alvis,1,indx_sfc) : )
   cmu      (1:n) => bus( x(bm,1,1)           : )
   ctu      (1:n) => bus( x(bt,1,indx_sfc)    : )
   fc_glac  (1:n) => bus( x(fc,1,indx_sfc)    : )
   if (radslope)  then
      fsol  (1:n) => bus( x(fluslop,1,1)      : )
   else
      fsol  (1:n) => bus( x(flusolis,1,1)     : )
   endif
   fv_glac  (1:n) => bus( x(fv,1,indx_sfc)    : )
   hst_glac (1:n) => bus( x(hst,1,indx_sfc)   : )
   ilmo_glac(1:n) => bus( x(ilmo,1,indx_sfc)  : )
   qsice    (1:n) => bus( x(qsurf,1,indx_sfc) : )
   snorate  (1:n) => bus( x(snowrate,1,1)     : )
   tdeep    (1:n) => bus( x(tglacier,1,2)     : )
   ts       (1:n) => bus( x(tglacier,1,1)     : )
   z0h      (1:n) => bus( x(z0t,1,indx_sfc)   : )
   z0m      (1:n) => bus( x(z0,1,indx_sfc)    : )
   zalfaq   (1:n) => bus( x(alfaq,1,1)        : )
   zalfat   (1:n) => bus( x(alfat,1,1)        : )
   zcoszeni (1:n) => bus( x(cang,1,1)        : )
   zdlat    (1:n) => bus( x(dlat,1,1)         : )
   zemisr   (1:n) => bus( x(emisr,1,1)        : )
   zfcor    (1:n) => bus( x(fcor,1,1)         : )
   zfdsi    (1:n) => bus( x(fdsi,1,1)         : )
   zfrv     (1:n) => bus( x(frv,1,indx_sfc)   : )
   zftemp   (1:n) => bus( x(ftemp,1,indx_sfc) : )
   zfvap    (1:n) => bus( x(fvap,1,indx_sfc)  : )
   zrunofftot(1:n) => bus( x(runofftot,1,indx_sfc) : )
   zsnodp   (1:n) => bus( x(snodp,1,indx_sfc) : )
   ztsrad   (1:n) => bus( x(tsrad,1,1)        : )
   ztsurf   (1:n) => bus( x(tsurf,1,indx_sfc) : )
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


   
!*     1.     Preliminaries
!      --------------------

      i = sl_prelim(tt,hu,uu,vv,ps,zzusl,spd_air=vmod0,dir_air=vdir,rho_air=rhoa, &
           min_wind_speed=sqrt(vamin))
      if (i /= SL_OK) then
         call physeterror('glaciers', 'error returned by sl_prelim()')
         return
      endif

      
      ! Set emissivity values
      select case (snow_emiss)
      case DEFAULT
         EMISNOW = snow_emiss_const
      end select
      select case (ice_emiss)
      case DEFAULT
         EMISICE = ice_emiss_const
      end select

      do I=1,N

!                           - in-glacier temperature

        T2(I) = TDEEP(I)

!                           Thermal roughness length for
!                           the ice surface


        Z0H(I) = max(   0.2 * Z0M(I)  ,  Z0GLA   )
        Z0H(I) = min(    Z0H(I)  ,  0.2   )


!                           Saturated specific humidity above
!                           ice surface

        QSICE(I) = FOQST( TS(I), PS(I) )
        DQSAT(I) = FODQS ( QSICE(I), TS(I) )

!       Damp snow depth with h0*tanh(ff/h0)

      if (LIMSNODP) then
        ZSNODP_M(I) = SNOH0*tanh(ZSNODP(I)/SNOH0)
      else
       ZSNODP_M(I)=ZSNODP(I)
      endif

      end do


!       2.     Calculate the drag and heat coefficients
!       -----------------------------------------------

      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
           coefm=cmu,coeft=ctu,flux_t=zftemp,flux_q=zfvap,ilmo=ilmo_glac,ue=zfrv, &
           h=hst_glac,L_min=sl_Lmin_glacier,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('glaciers', 'error returned by sl_sfclayer()')
         return
      endif 

!       3.     Parameterizations (albedo, emissivity, conductivity,...)
!       ---------------------------------------------------------------




!                           surface albedo (function of surface
!                           type and temperature)


      do I=1,N

        if( TS(I) .lt. TMELSNO ) then
          SCR2(I) = ALBDS
        else
          SCR2(I) = ALBMS
        endif

        if( TS(I) .lt. TMELICE ) then
          SCR3(I) = ALBDI
        else
          SCR3(I) = ALBMI
        endif

        if( ZSNODP(I) .gt. CON10 ) then
          ALBSFC(I) = SCR2(I)
        else
          ALBSFC(I) = min( SCR2(I) , SCR3(I)+ZSNODP(I)* &
                          (SCR2(I)-SCR3(I))/CON10 )
        endif

      end do

!                               conductivity and capacity (function of surface temperature)
      do I=1,N
!                                            snow layer
        if( TS(I) .lt. TMELSNO ) then
          SCR2(I) = ROSNOW(1)*HCAPS
          SCR4(I) = KSDS
        else
          SCR2(I) = ROSNOW(2)*HCAPS
          SCR4(I) = KSMS
        endif
!                                            ice layer
        SCR3(I) = ROICE*HCAPI
        SCR5(I) = CONDFI
!                                            depths D1 and D2
        SCR7(I) = -sqrt((SCR4(I)/SCR2(I))*DAY/PI)*ALOG(DEEPFRAC)
        SCR8(I) = -sqrt((SCR5(I)/SCR3(I))*DAY/PI)*ALOG(DEEPFRAC)
!                                            equivalent ice depth and
!                                            limit if Tp is in snow layer
        if( ZSNODP_M(I) .lt. SCR7(I) ) then
          SCR6(I) = SCR8(I)*( 1. -ZSNODP_M(I)/SCR7(I) )
          SCR10(I) = ZSNODP_M(I)
        else
          SCR6(I) = 0.0
          SCR10(I) = SCR7(I)
        endif
!                                            compute K/H
        SCR11(I) = SCR4(I)/( SCR10(I) + SCR6(I)*SCR4(I)/SCR5(I) )

      end do

      do I=1,N
!                               heat capacity of combined snow/ice layers H*C
        SCR5(I) = 0.5*( SCR10(I)*SCR2(I)*( 1.+SCR11(I)*SCR6(I)/SCR5(I) ) &
                + SCR3(I)*SCR11(I)*SCR6(I)**2/SCR5(I) )
        CT(I) = 1./SCR5(I)
        C2(I) = CT(I)*SCR11(I)

      end do

!                               emissivity and penetration of solar radiation
!                               (function of surface type snow/ice)
      do I=1,N

        if( ZSNODP(I) .gt. CON10 ) then
          ZEMISR(I) = EMISNOW
          SCR1(I) = 1.0
        else
          ZEMISR(I) = EMISICE
          SCR1(I) = 1.0 - FI0
          SCR1(I) = min( 1.0 , SCR1(I)+ZSNODP(I)* &
                        (1.0-SCR1(I))/CON10 )
        endif

        SCR2(I) = 1.-exp(-COEFEXT*SCR6(I))

      end do

!                                              terms A, B, and C for the
!                                              calculation of TS at time 'T+DT'
      do I=1,N

        A(I) = CT(I) * ( 4. * ZEMISR(I) * STEFAN * TS(I)**3 &
           +  RHOA(I) * CTU(I) * (DQSAT(I) * (CHLC+CHLF) + CPD) ) &
           + C2(I)

        B(I) = CT(I) * ( 3. * ZEMISR(I) * STEFAN  * TS(I)**3 &
           + RHOA(I) * CTU(I) * DQSAT(I) * (CHLC+CHLF) )

        C(I) = C2(I) * T2(I) + CT(I) * &
           (  RHOA(I)*CTU(I) * (CPD*TH(I) - (CHLC+CHLF)*(QSICE(I)-HU(I))) &
           +  FSOL(I)*(1.-ALBSFC(I)) * (SCR1(I)+(1.-SCR1(I))*SCR2(I)) &
           + ZEMISR(I)*ZFDSI(I) )

      end do


!       4.     Calculate the surface temperature from
!              force-restore equation at the snow or ice surface
!       -----------------------------------------------------------------

!                                       Energy balance for surface temperature
!  d TS / dt = -2*HA/(C*H)-2*K(TS-TP)/(C*H*H)
!                                       after linearization
!  d TS / dt = - A * TS+ + B * TS- + C
!                                       and deep temperature
!  d TP / dt = - ( TP+ -  TS+ ) / DAY
!                                       can be solved analytically
        do I=1,N
           TS(I) = TS(I) + ( TS(I)*(B(I)-A(I)) + C(I) ) &
                 * ( 1.-exp(-DELT*A(I)) ) / A(I)
           TDEEP(I) = (TDEEP(I)+TS(I)*DELT/DAY)/(1.0+DELT/DAY)
        end do


      do I=1,N
!                                              prevent surface and deep temperatures
!                                              from exceeding melting temperature
        if( ZSNODP(I) .gt. 0.0 ) then
          SCR4(I) = TMELSNO
        else
          SCR4(I) = TMELICE
        endif

        TS(I) = min ( TS(I) , SCR4(I) )
        TDEEP(I) = min ( TDEEP(I) , SCR4(I) )

!                               saturated specific humidity at
!                               surface

        QSICE(I) = FOQST( TS(I), PS(I) )

      end do

      ! Estimate diagnostic-level quantities
      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
           hghtm_diag=zu,hghtt_diag=zt,t_diag=ztdiag,q_diag=zqdiag,u_diag=zudiag, &
           v_diag=zvdiag,tdiaglim=GLACIER_TDIAGLIM,L_min=sl_Lmin_glacier,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('glaciers', 'error 2 returned by sl_sfclayer()')
         return
      endif
      if (sl_Lmin_glacier > 0.) then
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

!       5.     Melting and growth
!       -------------------------
!                                              melting of snow at the surface
!                                              accumulation of snow fall
      if(ICEMELT) then

        do I=1,N

!                               recompute surface albedo with new TS

          if( TS(I) .lt. TMELSNO ) then
            SCR2(I) = ALBDS
          else
            SCR2(I) = ALBMS
          endif

          if( TS(I) .lt. TMELICE ) then
            SCR3(I) = ALBDI
          else
            SCR3(I) = ALBMI
          endif

          if( ZSNODP(I) .gt. CON10 ) then
            ALBSFC(I) = SCR2(I)
          else
            ALBSFC(I) = min( SCR2(I) , SCR3(I)+ZSNODP(I)* &
                           (SCR2(I)-SCR3(I))/CON10 )
          endif

!                               compute K/H * (TB - TS)
          SCR11(I) = SCR11(I)*( TDEEP(I) - TS(I) )

!                               compute HA = H + LE + FWnet + FLnet
          SCR5(I) = RHOA(I)*CTU(I)* &
                    ( CPD*(TS(I)-TH(I))+(CHLC+CHLF)*(QSICE(I)-HU(I)) ) &
                   - FSOL(I)*(1.-ALBSFC(I))*SCR1(I) &
                   + ZEMISR(I)*( STEFAN*TS(I)**4 - ZFDSI(I) )

          SCR3(I) = SCR5(I) - SCR11(I)

!                               criteria for melting at the surface
          ZRUNOFFTOT(I) = 0.0
          if( TS(I).ge.SCR4(I) .and. &
              SCR5(I).le.0.0 .and. SCR3(I).lt.0.0 ) then

!                               melt available snow...
            ZSNODP(I) = max( 0.0, ZSNODP(I) + DELT*SCR3(I)/VHFSNO )
!                               (in units of water equivalent - kg/m2 or mm)
            ZRUNOFFTOT(I) = - DELT*SCR3(I)/HFICE
!
          endif
!                               add snow fall (units changed from
!                               water equivalent to snow equivalent)
          if( TS(I).lt.TMELICE ) then
            ZSNODP(I) = ZSNODP(I) + DELT*SNORATE(I)*(1000./ROSNOW(1))
          endif

        end do

      endif


!       6.     The fluxes
!       -----------------

      ! Compute diagnostic quantities at 1.5m
      i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor,  &
           hghtt_diag=ZT_RHO,t_diag=my_ta,q_diag=my_qa,tdiaglim=GLACIER_TDIAGLIM, &
           L_min=sl_Lmin_glacier,spdlim=vmod)
      if (i /= SL_OK) then
         call physeterror('glaciers', 'error 3 returned by sl_sfclayer()')
         return
      endif

!VDIR NODEP
      do I=1,N

        ZTSURF    (I) = TS (I)
        ZTSRAD    (I) = TS (I)

        ZALFAT    (I) = - CTU(I) * ( TS   (I)-TH(I) )
        ZALFAQ    (I) = - CTU(I) * ( QSICE(I)-HU(I) )
        if (.not.IMPFLX) CTU (I) = 0.
        RHOA      (I) = PS(I)/(RGASD * my_ta(i)*(1.+DELTA*my_qa(i)))
        FC_GLAC(I)    = -CPD *RHOA(I)*ZALFAT(I)
        FV_GLAC(I)    = -(CHLC+CHLF)*RHOA(I)*ZALFAQ(I)

        if (IMPFLX) then
          ZALFAT    (I) = - CTU(I) *  TS(I)
          ZALFAQ    (I) = - CTU(I) *  QSICE(I)
        endif
!***

      end do
 
      !------------------------------------
      !   7.     Heat Stress Indices
      !------------------------------------
      !#TODO: at least 4 times identical code in surface... separeted s/r to call
      IF_THERMAL_STRESS: if (thermal_stress) then

         i = sl_sfclayer(th,hu,vmod0,vdir,zzusl,zztsl,ts,qsice,z0m,z0h,zdlat,zfcor, &
              hghtm_diag=zt,hghtt_diag=zt,u_diag=zusurfzt, &
              v_diag=zvsurfzt,tdiaglim=GLACIER_TDIAGLIM,L_min=sl_Lmin_glacier,spdlim=vmod)
         if (i /= SL_OK) then
            call physeterror('glaciers', 'error 3 returned by sl_sfclayer()')
            return
         endif
         zusurfzt = zusurfzt * vmod0 / vmod
         zvsurfzt = zvsurfzt * vmod0 / vmod
         
         do i=1,N

            if (abs(zzusl(i)-zu) <= 2.0) then
               zu10(i) = sqrt(uu(i)**2+vv(i)**2)
            else
               zu10(i) = sqrt(zudiag(i)**2+zvdiag(i)**2)
            endif

            ! wind  at SensoR level zusr at z=zt
            if( (abs(zusurfzt(i)) >= 0.1) .and. (abs(zvsurfzt(i)) >= 0.1)) then
               zusr(i) = sqrt( zusurfzt(i)**2 + zvsurfzt(i)**2)
            else
               zusr(i) = zu10(i)
            endif

            zqd(i) =  max(ZQDIAG(i) , 1.e-6)

            zref_sw_surf(i) = albsfc(i) * fsol(i)
            zemit_lw_surf(i) = (1. -zemisr(i)) * zfdsi(i) + zemisr(i)*stefan*ztsurf(i)**4

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

      endif IF_THERMAL_STRESS

      !     FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
      call FILLAGG( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_GLACIER, &
           SURFLEN)

   return
end subroutine glaciers2
