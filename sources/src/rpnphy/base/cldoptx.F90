!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

!/@*
subroutine cldoptx6(LWC,LWP,IWP,neb,T,sig,ps, &
     lat,mg,ml,m,lmx,nk, &
     pbl,ipbl,dz,sdz,eneb,opdepth,asymg, &
     topthw,topthi, ctp,ctt, &
     omegav,tauae,ioptix)
   use tdpack, only: GRAV, PI, RGASD
   implicit none
!!!#include <arch_specific.hf>
   !@OBJECT
   !     computes optical parameters as input to visible and infrared
   !             radiation also includes aerosol parameterization
   !             Optical parameters refer to entire VIS or IR spectrum
   !             but could be extended to several bands matching those
   !             of the radiation codes.
   !@ARGUMENTS
   !          - Output -
   ! ENEB     cloud amount times emissivity in each layer (0. to 1.)
   !          (effective nebulosity to be used by IR code)
   ! OPDEPTH  layer visible optical depth (dimensionless)
   ! ASYMG    layer visible asymmetry factor (G in literature, 0. to 1. )
   ! OMEGAV   layer visible single scattering albedo (0. to 1.)
   ! TAUAE    layer aerosol optical depth for VIS code
   ! IPBL     closest model level matching PBL (LMX)
   ! TOPTHW   total integrated optical thickness of water (from TLWP)
   ! TOPTHI   total integrated optical thickness of ice   (from TIWP)
   !
   !          -Output -
   ! LWC      TOTAL (liquid and solid) cloud water content for
   !          CONDS and OLDSUND schemes (cw_rad=0).
   !          Units : Kg water/Kg air (caution: not in Kg/m3) (LMX,NK)
   !
   !          -Input -
   ! LWC      * TOTAL cloud water content for NEWSUND, CONSUN, EXMO
   !            and WARM K-Y condensation schemes (cw_rad=1);
   !          * LIQUID water content for MIXED PHASE and
   !            COLD K-Y schemes (cw_rad=2).
   !
   !          -Input -
   ! NEB      layer cloud amount (0. to 1.) (LMX,NK)
   ! T        layer temperature (K) (M,NK)
   ! SIG      sigma levels (0. to 1.) (LMX,NK; local sigma)
   ! LAT      latitude in radians (-pi/2 to +pi/2) (LMX)
   ! PBL      height of planetary boundary layer in meters (LMX)
   ! DZ       work array, on output geometrical thickness (LMX,NK)
   ! SDZ      work array (LMX)
   ! MG       ground cover (ocean=0.0,land <= 1.)  (LMX)
   ! ML       fraction of lakes (0.-1.) (LMX)
   ! PS       surface pressure (N/m2) (LMX)
   ! LMX      number of profiles to compute
   ! M        first dimension of temperature (usually LMX)
   ! NK       number of layers
   ! IOPTIX   parameterizations for cloud optical properties
   !          = 1 for simpler condensation schemes
   !          = 2 for microphysics schemes

   integer lmx,m,nk,ioptix
   real ipbl(lmx)
   real LWC(LMX,nk), LWP(LMX,nk), IWP(LMX,nk)
   real neb(LMX,nk), t(m,nk), sig(LMX,nk)
   real ps(LMX),lat(LMX),eneb(LMX,nk),mg(LMX),ml(LMX)
   real opdepth(LMX,nk),asymg(LMX,nk),omegav(LMX,nk)
   real tauae(LMX,nk,5),pbl(LMX),dz(LMX,nk),sdz(LMX)
   real topthw(lmx),topthi(lmx)
   real ctp(lmx),ctt(lmx)

   !@AUTHOR L. Garand and H. Barker (April 1995)
   !@REVISION
   ! 001 R. Sarrazin and L. Garand (May 1996) - Correct bug for omegav
   !                                            and change tuneopi
   ! 002 N. Brunet (Oct 96) Correct bug for mg
   ! 003 C. Beaudoin (Jan 98) Eliminate fictitious stratospheric clouds
   !                          above 50 mb for CONDS condensation option
   ! 004 B. Bilodeau and L. Garand (Aug 1999) - Add IWC for interaction
   !                                           with microphysics schemes
   ! 005 B. Dugas (April 1999) Never any clouds above 70 mb, but this only
   !                           when the new input parameter climat is true
   ! 006 A. Methot and L. Garand (Jun 2000) - introduce a maximum in the
   !                                         total optical depth
   ! 007 A. Methot (May 2000) - modify effective radius relationship
   !                           with lwc
   ! 008 A. Methot and L. Garand (Jun 2000) - introduce Hu and Stamnes
   !                                         parameters for liquid water
   ! 009 A. Methot and Mailhot (Jun 2000) - introduce Fu & Liou parameters
   !                                       for ice
   ! 010 B. Bilodeau (Mar 2001) - Old cldotpx code as option
   ! 011 B. Bilodeau (Nov 2001) - Tuneopi = 0.8 for new optical properties
   ! 012 B. Bilodeau (Nov 2002) - Back to old optical properties for ice
   !                              Lakes treated as land
   ! 013 B. Bilodeau, P. Vaillancourt and A. Glazer (Dec 2002)
   !                            - Calculate ctp and ctt
   ! 014 B. Dugas - Rename CLIMAT to NOSTRLWC
   !
   ! 015 M. Lepine  (March 2003) -  CVMG... Replacements
   ! 016 D. Talbot (June 2003) - IBM conversion
   !                - calls to vsexp routine (from massvp4 library)
   !                - calls to exponen4 (to calculate power function '**')
   !                - divisions replaced by reciprocals
   ! 017 B. Bilodeau (Aug 2003) - exponen4 replaced by vspow
   ! 018 B. Bilodeau (Feb 2004) - take ml into account for aerosols
   !                              when ioptix == OPT_OPTIX_OLD
   ! 038      L. Spacek (Aug 2004) - cloud clean-up, cloud and liquid/ice water
   !                                 calculations moved out of cldoptx4
   !                                 neb= cloud even for IOPTIX=OPT_OPTIX_OLD
   !                                 remove NOSTRLWC
   !*@/

   integer, parameter :: OPT_OPTIX_OLD = 1
   integer, parameter :: OPT_OPTIX_NEW = 2

   real, dimension(LMX,NK  ) :: CLOUD
   real, dimension(LMX,NK  ) :: DP
   real, dimension(LMX,NK  ) :: TRANSMISSINT
   logical, dimension(LMX     ) :: TOP

   integer i,k
   integer ire(m,nk)
   real tuneopw
   real rei, rec_rei, ki, tuneopi, omi, gi, ssai
   real dp1,dp2,dp3
   real elsa, emiss
   real third,rec_grav,rec_180,rec_rgasd
   real ct,aero,rlat,eps
   real aird(m,nk),rew(m,nk),rec_cdd(m,nk),kw(m,nk)
   real kwf_ire(m,nk),kvf_ire(m,nk),ssf_ire(m,nk),gwf_ire(m,nk)
   real ei(m,nk),omw(m,nk),ssaw(m,nk),gw(m,nk)
   real ew(m,nk)
   real kwf(3,3), kvf(3,3), ssf(3,3), gwf(3,3)


   ! LWP,IWP are liquid, ice water paths in g/m2
   data third/0.3333333/
   ! diffusivity factor of Elsasser
   data elsa/1.66/
   data eps/1.e-10/
   save third,elsa,eps

   rec_grav=1./GRAV
   rec_rgasd=1./RGASD

   IF_OPTIX_NEW: if (IOPTIX == OPT_OPTIX_NEW) then

      ! ************************************************************
      !           LIQUID WATER parameterization after Hu and Stamnes
      ! -----------------------------------------------------------
      !
      !        Parameters for the relationship between the
      !        diffusivity factor and equivalent radius
      !        at 11um wavelenght ( from table 4)
      !        After Hu and Stamnes 1993, JCAM p 728-742
      !                              2.5 < radius < 12.5 um
      kwf(1,1) =   -3.88e-05
      kwf(2,1) =    5.24
      kwf(3,1) =  140.
      !                             12.5 <= radius <= 30 um
      kwf(1,2) =  794.
      kwf(2,2) =    -.148
      kwf(3,2) = -423.
      !                             30.0 <  radius <  60 um
      kwf(1,3) = 1680.
      kwf(2,3) =    -.966
      kwf(3,3) =   -5.84

      !        Parameters for the relationship between the
      !        optical thickness  and equivalent radius
      !        at .719 um wavelenght ( from table 1)
      !        After Hu and Stamnes 1993, JCAM p 728-742
      !                              2.5 < radius < 12.5 um
      kvf(1,1) = 1810.
      kvf(2,1) =   -1.08
      kvf(3,1) =    6.85
      !                             12.5 <= radius <= 30 um
      kvf(1,2) = 1700.
      kvf(2,2) =   -1.04
      kvf(3,2) =    1.04
      !                             30.0 <  radius <  60 um
      kvf(1,3) =  978.
      kvf(2,3) =    -.816
      kvf(3,3) =   -9.89

      !        Parameters for the relationship between the
      !        asymmetry factor  and equivalent radius
      !        at .719 um wavelenght ( from table 2)
      !        After Hu and Stamnes 1993, JCAM p 728-742
      !                              2.5 < radius < 12.5 um
      ssf(1,1) =    9.95e-7
      ssf(2,1) =    -.856
      ssf(3,1) =   -4.37e-7
      !                             12.5 <= radius <= 30 um
      ssf(1,2) =    1.88e-7
      ssf(2,2) =    1.32
      ssf(3,2) =    3.08e-6
      !                             30.0 <  radius <  60 um
      ssf(1,3) =    2.03e-5
      ssf(2,3) =    -.332
      ssf(3,3) =   -4.32e-5

      !        Parameters for the relationship between the
      !        single scatering albedo and equivalent radius
      !        at .719 um wavelenght ( from their table 3)
      !        After Hu and Stamnes 1993, JCAM p 728-742
      !                              2.5 < radius < 12.5 um
      gwf(1,1) = -.141
      gwf(2,1) = -.694
      gwf(3,1) =  .889
      !                             12.5 <= radius <= 30 um
      gwf(1,2) = -.157
      gwf(2,2) = -.782
      gwf(3,2) =  .886
      !                             30.0 <  radius <  60 um
      gwf(1,3) = -.214
      gwf(2,3) = -.916
      gwf(3,3) =  .885

   endif IF_OPTIX_NEW

! ************************************************************
!           MISCELLANEOUS
! -----------------------------------------------------------
!
! tuning for optical thickness in visible
! (ONLY VIS IS AFFECTED)
! this tuning affects outgoing and surface
! radiation; has little impact on atmosperic heating

      if (ioptix == OPT_OPTIX_OLD) then
         tuneopw = 0.30
         tuneopi = 0.30
      else if (ioptix == OPT_OPTIX_NEW) then
         tuneopw = 0.80
         tuneopi = 0.30
      endif

! initialize output fields to zero

      do i=1,lmx
         topthw(i) = 0.0
         topthi(i) = 0.0
      end do
!
! ************************************************************
!                        PRELIMINARY WORK
! -----------------------------------------------------------
!
      do k=1,nk
      do I=1,lmx
         dp1=0.5*(sig(i,min(k+1,nk))-sig(i,max(k-1,1)))
         dp2=0.5*(sig(i,1)+sig(i,2))
         dp3=0.5*(1.-sig(i,nk))
         if (k .eq. 1) then
            dp(i,k) = dp2
         else if (k .eq. nk) then
            dp(i,k) = dp3
         else
            dp(i,k) = dp1
         endif

         dp(i,k)=max(dp(i,k)*ps(i),0.)

         cloud(i,k)=neb(i,k)
      end do
      end do
!
! *****************************************************************
!     MAIN LOOP FOR MICROPHYSICS SCHEMES ("NEW" OPTICAL PROPERTIES)
! -----------------------------------------------------------------
!
      if (IOPTIX == OPT_OPTIX_NEW) then

      do k=1,nk
      do I=1,lmx
!                        EQUIVALENT SIZE of PARTICLES
!
!------------->    determines equivalent size of WATER particles
!                  set number of drops per cm**3 100 for water
!                  and 500 for land

            if (mg(i) .le. 0.5 .and. ml(i) .le. 0.5) then
!              cdd(i,k) = 100.
               rec_cdd(i,k) = 1. / 100.
            else
!              cdd(i,k) = 500.
               rec_cdd(i,k) = 1. / 500.
            endif
!
!                  aird is air density in kg/m3
!
            aird(i,k) = sig(i,k)*ps(i)/T(i,k)*REC_RGASD
!
      end do
      end do
!
      call VSPOWN1  (REW,(1.+LWC*1.e4)*LWC*aird*rec_cdd, &
                      third, nk*lmx)
!
      do k=1,nk
      do I=1,lmx
!           REW(i,k) = 3000. * ( (1.+LWCM(i,k)*1.e4)*
!    +            LWCM(i,k)*frac(i,k)*aird(i,k)*rec_cdd(i,k))**third
            REW(i,k) = 3000. * REW(i,k)
            REW(i,k) = max(  2.5,REW(i,k))
            REW(i,k) = min( 60., REW(i,k))
!
!                  determines array index for given REW

                                      ire(i,k) = 2
            if ( rew(i,k) .lt. 12.5 ) ire(i,k) = 1
            if ( rew(i,k) .gt. 30.0 ) ire(i,k) = 3
!           avant d'appeler vspownn il faut definir
            kwf_ire(i,k) = kwf(2,ire(i,k))
            kvf_ire(i,k) = kvf(2,ire(i,k))
            ssf_ire(i,k) = ssf(2,ire(i,k))
            gwf_ire(i,k) = gwf(2,ire(i,k))
!
      end do
      end do
!
      call VSPOWNN  (kw,rew,kwf_ire,nk*lmx)
!
      do k=1,nk
      do I=1,lmx
!
!                  water diffusivity after Hu and Stammes, JCAM 1993
!
!           kw      = ( kwf(1,ire(i,k))*( rew(i,k)**kwf(2,ire(i,k)) ) + kwf(3,ire(i,k)) )*1.e-3
            kw(i,k) = ( kwf(1,ire(i,k)) * kw(i,k)                     + kwf(3,ire(i,k)) )*1.e-3
!
      end do
      end do

            REI = 15.
            REC_REI = 1. / 15.
            KI = .0003 + 1.290 * REC_REI
            call VSEXP (EI,-elsa*ki*IWP,nk*lmx)
            call VSEXP (ENEB,-elsa*ki*IWP-kw*LWP,nk*lmx)
!
!           on elimine une boucle en faisant ceci plus tot
            call VSPOWNN  (omw,rew,kvf_ire,nk*lmx)
            call VSPOWNN  (ssaw,rew,ssf_ire,nk*lmx)
            call VSPOWNN  (gw,rew,gwf_ire,nk*lmx)
            SSAI = 1.0 - 1.295E-2 - 1.321E-4 *REI
!
!
      do k=1,nk
      do I=1,lmx
!
!
!                  emissivity of ice after Ebert and Curry, JGR-1992,p3833
!                  using 11.1 micron parameters as representative
!                  of entire spectrum assume equivalent ice particles
!                  radius of 25 micron will affect both VIS and
!                  IR radiation
!
!           REI = 15.
!           KI = .0003 + 1.290 / REI
!           EI = 1. -exp(-elsa*ki *IWP(i,k))
            EI(i,k) = 1. - EI(i,k)
!
!
!                  compute combined ice/water cloud emissivity assuming
!                  the transmission is the product of the ice and water
!                  phase transmissions
!
!                  cloud amount is considered outside of the current
!                  main loop due to potential optical depth correction
!
!                  cloud emissivity temporarly in ENEB

!           ENEB(i,k) = 1. - exp( - elsa*ki *IWP(i,k) - kw(i,k) * LWP(i,k) )
            ENEB(i,k) = 1. -   ENEB(i,k)
!
!     end do
!     end do
!
!
!     do k=1,nk
!     do I=1,lmx
!                        OPTICAL THICKNESS
!
!                   water optical thickness: Hu and Stammes, JCAM 1993
!                                                for .719 um
!
!           omw     =LWP(i,k) * ( kvf(1,ire(i,k))*( rew(i,k)**kvf(2,ire(i,k)) ) + kvf(3,ire(i,k)) )*1.e-3
            omw(i,k)=LWP(i,k) * ( kvf(1,ire(i,k))*( omw(i,k)                  ) + kvf(3,ire(i,k)) )*1.e-3
            omw(i,k)=omw(i,k)*tuneopw

            OMI = min(IWP(i,k)*(3.448E-3+2.431*REC_REI) * tuneopi, 25.)

            OPDEPTH(i,k)= max(omw(i,k) + omi,1.e-10)
!
!                 save integrated quantities for output purposes
!
         topthw(i) = topthw(i) + omw(i,k)
         topthi(i) = topthi(i) + omi
!
!
!                        SINGLE SCATTERING ALBEDO
!
!                 note that this parameter is very close to one for
!                 most of VIS but lowers in near infrared; proper
!                 spectral weighting is important because small changes
!                 will affect substantially the solar heating outgoing
!                 radiance or planetary albedo; It has a smaller influence
!                 on the surface flux. A SSA of unity will create division
!                 by zero in solar code.

!           SSAI = 1.0 - 1.295E-2 - 1.321E-4 *REI
!

!                 WATER  Single scattering Albedo: Hu and Stammes, JCAM 1993
!                                                for .719 um

!           SSAW     = 1.- ( ssf(1,ire(i,k))*( rew(i,k)**ssf(2,ire(i,k)) ) + ssf(3,ire(i,k)) )
            SSAW(i,k)= 1.- ( ssf(1,ire(i,k))*( SSAW(i,k)                 ) + ssf(3,ire(i,k)) )
!
!                  weighting by optical depth
!
         OMEGAV(i,k)= (OMW(i,k) * SSAW(i,k) + OMI* SSAI)/ (OMW(i,k)+OMI+eps)
         OMEGAV(i,k)=max(OMEGAV(i,k),0.9850)
         OMEGAV(i,k)=min(OMEGAV(i,k),0.999999)
!
!                  Ice   Asymmetry factor: Ebert and Curry JGR 1992 Eq.8
!                        for 25 micron particles and a weighting for the
!                        five bands given in their Table 2
!
            GI = min(0.777 + 5.897E-4 * REI, 0.9999)
!
!
!                  Water  Asymetry factor: Hu and Stammes, JCAM 1993
!                                                for .719 um

!           GW      = ( gwf(1,ire(i,k))*( rew(i,k)**gwf(2,ire(i,k)) ) + gwf(3,ire(i,k)) )
            GW(i,k) = ( gwf(1,ire(i,k))*( GW(i,k)                   ) + gwf(3,ire(i,k)) )
!
!                  weighting by SSA * opt depth
!
         ASYMG(i,k) = (SSAI*OMI*GI + SSAW(i,k)*OMW(i,k)*GW(i,k))/(SSAI*OMI+SSAW(i,k)*OMW(i,k)+eps)
!
         asymg(i,k) = max(asymg(i,k), 0.75 )
!
         asymg(i,k) = min(asymg(i,k),0.9999)
!
!                  geometrical thickness
!
         dz(i,k)= dp(i,k)/aird(i,k)*rec_grav
!
!
      end do
      end do

      else if (IOPTIX == OPT_OPTIX_OLD) then
!
! ****************************************************************
!     MAIN LOOP FOR CONVENTIONAL CONDENSATION SCHEMES
!     ("OLD" OPTICAL PROPERTIES)
! ----------------------------------------------------------------
!
            REI = 15.
            REC_REI = 1. / 15.
            KI = .0003 + 1.290 * REC_REI
!
            call VSEXP (EW,-0.087*elsa*LWP,nk*lmx)
            call VSEXP (EI,-elsa*ki*IWP,nk*lmx)
!
      do k=1,nk
      do I=1,lmx
!
!                  emissivity of water after Stephens, JAS 78
!                  we take a 0.144 factor as average of his
!                  downward and upward emissivity which divided
!                  by diffusivity factor elsa yields 0.087
!
!           EW =      1. -exp(-0.087*elsa* LWP(i,k))
            EW(i,k) = 1. - EW(i,k)
!
!                  emissivity of ice after Ebert and Curry, JGR-1992,p3833
!                  using 11.1 micron parameters as representative
!                  of entire spectrum assume equivalent ice particles
!                  radius of 25 micron will affect both VIS and
!                  IR radiation
!
!           EI(i,k) = 1. -exp(-elsa*ki *IWP(i,k))
            EI(i,k) = 1. - EI(i,k)
!
!                  compute combined ice/water cloud emissivity assuming
!                  the transmission is the product of the ice and water
!                  phase transmissions
!
            EMISS = 1. - (1.-EI(i,k))* (1.-EW(i,k))
!
!                   effective cloud cover
!
            ENEB(i,k)= cloud(i,k)*emiss
!
!                 optical thickness at nadir is computed
!                 set number of drops per cm**3 100 for water
!                 and 500 for land
!
!           The following line is commented out to ensure bitwise
!           validation of the GEMDM global model. It should be
!           activated in a future version of the code.
            if (mg(i).le.0.5) then
!              cdd(i,k) = 100.
               rec_cdd(i,k) = 1. / 100.
            else
!              cdd(i,k) = 500.
               rec_cdd(i,k) = 1. / 500.
            endif
!
!                 aird is air density in kg/m3
            aird(i,k) = sig(i,k)*ps(i)/T(i,k)*REC_RGASD
!
      end do
      end do
!
!		  check for small negative values in the LWC
!		  variable and clip quietly if appropriate.
!
            if (minval(LWC) < 1.e-10) then
              if (minval(LWC) < -1.e-10) then
                call physeterror('cldoptx', 'large negative LWC')
                return
              else
                LWC = max(LWC,1.e-10)
              endif
            endif
!
            call VSPOWN1  (REW,LWC*aird*rec_cdd,third,nk*lmx)
!
      do k=1,nk
      do I=1,lmx
!
!                 this parameterization from H. Barker, based
!                 on aircraft data range 4-17 micron is that specified
!                 by Slingo for parameterizations
!
!           REW(i,k) = max(4., 754.6 * (LWC(i,k)*aird(i,k)*rec_cdd(i,k))**third)
            REW(i,k) = max(4., 754.6 * REW(i,k))
            REW(i,k) = min(REW(i,k),17.)

!
!                 slingo JAS 89 p 1420 for weighted average of
!                 bands 1-5 (all spectrum)
!
            OMW(i,k) =  LWP(i,k)* (2.622E-2 + 1.356/REW(i,k) )
!
!                 follows Ebert and Curry, 1992
!                 no variation as function of band
!                 ice optical depth limited to 25; water to 125.
!
            omw(i,k)= min(omw(i,k)*tuneopw,125.)
            OMI = min(IWP(i,k)*(3.448E-3+2.431*REC_REI) * tuneopi, 25.)
            OPDEPTH(i,k)= max(omw(i,k) + omi,1.e-10)
!
!
!                 save integrated quantities for output purposes
!
         topthw(i) = topthw(i) + omw(i,k)
         topthi(i) = topthi(i) + omi
!
            SSAI = 1.0 - 1.295E-2 - 1.321E-4 *REI
!           Slingo JAS 1989 for weighted average bands 1-4 p 1420
            SSAW(i,k) = 1.0 - 6.814E-3 - 4.842E-4 *REW(i,k)
!
!
!                  weighting by optical depth
!
         OMEGAV(i,k)= (OMW(i,k) * SSAW(i,k) + OMI* SSAI)/ (OMW(i,k)+OMI+eps)
         OMEGAV(i,k)=max(OMEGAV(i,k),0.9850)
         OMEGAV(i,k)=min(OMEGAV(i,k),0.9999)
!
!
!                  Ice   Asymmetry factor: Ebert and Curry JGR 1992 Eq.8
!                        for 25 micron particles and a weighting for the
!                        five bands given in their Table 2
!
            GI = min(0.777 + 5.897E-4 * REI, 0.9999)
!
!                  Water  Asymetry factor: Slingo 1989
!
            GW(i,k) = min(0.804 + 3.850E-3*REW(i,k), 0.9999)
!
!                  weighting by SSA * opt depth
!
         ASYMG(i,k) = (SSAI*OMI*GI + SSAW(i,k)*OMW(i,k)*GW(i,k))/(SSAI*OMI+SSAW(i,k)*OMW(i,k)+eps)
!
            asymg(i,k) = max(asymg(i,k),GI)
!
         asymg(i,k) = min(asymg(i,k),0.9999)
!
!                  geometrical thickness
!
         dz(i,k)= dp(i,k)/aird(i,k)*rec_grav
!
      end do
      end do
!
!
      endif
!
!
!
! ************************************************************
!                        END OF MAIN LOOPS
! -----------------------------------------------------------

!
       if (IOPTIX == OPT_OPTIX_NEW) then
!
! ******************************************************************
!                  CORRECTION FOR MAXIMUM TOTAL OPTICAL DEPTH &
!                  FINAL CLOUD*EMISSIVITY CALCULATIONS
! ------------------------------------------------------------------
!
!                  temporary use of ipbl as a work field.
!                  ipbl is 1. when there is no correction
!                  otherwise ipbl is a number between zero and one
!
         do i=1, lmx
            ipbl(i)=min( 1., 20./( max(topthi(i)+topthw(i), 1.e-10) ) )
         enddo
!                  only optical depth and cloud emissivity
!                  need to be re-scaled
!
         do k=1, nk
         do i=1, lmx
            OPDEPTH(i,k) = ipbl(i) * OPDEPTH(i,k)
            eneb(i,k) = neb(i,k) * &
                        ( 1. - ( 1.-eneb(i,k) )**ipbl(i) )
         enddo
         enddo
!
!
      endif
!
!     Diagnostics : cloud top pressure (ctp) and temperature (ctt)
!
      do i=1,lmx
         ctp (i)   = 110000.
         ctt (i)   = 310.
         top(i) = .true.
         transmissint(i,1) = 1. - eneb(i,1)
         if ( (1.-transmissint(i,1)) .gt. 0.99 .and. top(i) ) then
            ctp(i) = sig(i,1)*ps(i)
            ctt(i) = t(i,1)
            top(i) = .false.
         end if
      end do

      do k=2,nk
         do i=1,lmx
            transmissint(i,k) = transmissint(i,k-1) * (1.-eneb(i,k))
            if ( (1.-transmissint(i,k)) .gt. 0.99 .and. top(i) ) then
               ctp(i) = sig(i,k)*ps(i)
               ctt(i) = t(i,k)
               top(i) = .false.
            end if
         end do
      end do
!
!
! ******************************************************************
!                          AEROSOLS LOADING
! ------------------------------------------------------------------
!
        do 10 I =1,lmx
        sdz(i) =  0.
        ipbl(i) = 0.
 10     continue
!
        do 5 k=nk,1,-1
        do 6 I=1,lmx
           if ( int(ipbl(i)).eq.0 ) then
!                   pbl heiht recomputed as sum of layer thicknesses
             sdz(i)=sdz(i)+dz(i,k)
!                                              level closest to pbl
             if (sdz(i) .gt. pbl(i)) ipbl(i) = float(k)
           endif
 6      continue
 5      continue
!
!                  distributing aerosols
!                  optical thickness higher over land;
!                  decreases with latitude
!                  See Toon and Pollack, J. Appl. Meteor, 1976, p.235
!                  aero being total optical thickness, we distribute it
!                  within PBL in proportion with geometrical thickness
!                  above pbl aerosol optical depth is kept negligible
!
        ct=2./PI
        REC_180=1./180.
        do 3 k=1,nk
        do 4 i=1,lmx
             tauae(i,k,1)=1.e-10
             tauae(i,k,2)=1.e-10
             rlat = lat(i)*PI*REC_180
             if ( k.ge.int(ipbl(i)) ) then
!
                if ( mg(i).ge.0.5.or.ml(i).ge.0.5) then
!
!                  over land
!
                      aero = 0.25 - 0.2*ct * abs(lat(i))
                      tauae(i,k,1)= max(aero * dz(i,k)/sdz(i), 1.e-10)
                else
!
!                  over ocean
!
                      aero = 0.13 - 0.1*ct * abs(lat(i))
                      tauae(i,k,2)= max(aero  *dz(i,k)/sdz(i), 1.e-10)
                endif
             endif
!
!                  other types of aerosols set to negligible
!
             tauae(i,k,3)=1.e-10
             tauae(i,k,4)=1.e-10
             tauae(i,k,5)=1.e-10
 4      continue
 3      continue
        return
        end
