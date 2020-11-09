module my_dmom_mod
  use tdpack
  implicit none
!!!#include <arch_specific.hf>

  private
  public :: mydmom_main

  contains

!_______________________________________________________________________________________!
 SUBROUTINE MYDMOM_MAIN(WZ,T,Q,QC,QR,QI,QN,QG,QH,NC,NR,NY,NN,NG,NH,PS,TM,QM,QCM,        &
     QRM,QIM,QNM,QGM,QHM,NCM,NRM,NYM,NNM,NGM,NHM,PSM,S,RT_rn1,RT_rn2,RT_fr1,RT_fr2,     &
     RT_sn1,RT_sn2,RT_sn3,RT_pe1,RT_pe2,RT_peL,RT_snd,GZ,T_TEND,Q_TEND,QCTEND,QRTEND,   &
     QITEND,QNTEND,QGTEND,QHTEND,NCTEND,NRTEND,NYTEND,NNTEND,NGTEND,NHTEND,dt,NI,N,NK,  &
     J,KOUNT,CCNtype,precipDiag_ON,sedi_ON,warmphase_ON,autoconv_ON,icephase_ON,snow_ON,&
     initN,dblMom_c,dblMom_r,dblMom_i,dblMom_s,dblMom_g,dblMom_h,Dm_c,Dm_r,Dm_i,Dm_s,   &
     Dm_g,Dm_h,ZET,ZEC,SLW,VIS,VIS1,VIS2,VIS3,h_CB,h_ML1,h_ML2,h_SN,SS01,SS02,SS03,SS04,&
     SS05,SS06,SS07,SS08,SS09,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,   &
     Tc_3comp,PRFLXKY,SWFLXKY,F12KY,FEVPKY,cloud_bin,reff_c,reff_i1,       &
     reff_i2,reff_i3, reff_i4)

  use my_fncs_mod
  use my_sedi_mod

  implicit none

!CALLING PARAMETERS:
  integer,               intent(in)    :: NI,NK,N,J,KOUNT,CCNtype
  real,                  intent(in)    :: dt
  real, dimension(:),    intent(in)    :: PS,PSM,Tc_3comp
  real, dimension(:),    intent(out)   :: h_CB,h_ML1,h_ML2,h_SN
  real, dimension(:),    intent(out)   :: RT_rn1,RT_rn2,RT_fr1,RT_fr2,RT_sn1,RT_sn2,   &
                                          RT_sn3,RT_pe1,RT_pe2,RT_peL,ZEC,RT_snd
  real, dimension(:,:),  intent(in)    :: WZ,S,GZ
  real, dimension(:,:),  intent(inout) :: T,Q,QC,QR,QI,QN,QG,QH,NC,NR,NY,NN,NG,NH,     &
        TM,QM,QCM,QRM,QIM,QNM,QGM,QHM,NCM,NRM,NYM,NNM,NGM,NHM
  real, dimension(:,:),  intent(out)   :: T_TEND,QCTEND,QRTEND,QITEND,QNTEND,          &
        QGTEND,QHTEND,Q_TEND,NCTEND,NRTEND,NYTEND,NNTEND,NGTEND,NHTEND,ZET,Dm_c,       &
        Dm_r,Dm_i,Dm_s,Dm_g,Dm_h,SLW,VIS,VIS1,VIS2,VIS3,SS01,SS02,SS03,SS04,SS05,SS06, &
        SS07,SS08,SS09,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20,reff_c,  &
        reff_i1,reff_i2,reff_i3,reff_i4
  real, dimension(NI,NK),intent(out)   :: cloud_bin             ! cloud fraction (0 or 1)     
  logical,               intent(in)    :: dblMom_c,dblMom_r,dblMom_i,dblMom_s,         &
        dblMom_g,dblMom_h,precipDiag_ON,sedi_ON,icephase_ON,snow_ON,warmphase_ON,      &
        autoconv_ON,initN
  real, dimension(size(QC,dim=1),size(QC,dim=2)),   intent(out) :: prflxky,swflxky
  real, dimension(size(QC,dim=1),size(QC,dim=2)),   intent(out) :: f12ky,fevpky

!_______________________________________________________________________________________
!                                                                                       !
!                    Milbrandt-Yau Multimoment Bulk Microphysics Scheme                 !
!                              - double-moment version   -                              !
!_______________________________________________________________________________________!
!  Package version:   2.18.5              (internal bookkeeping)                        !
!  Last modified  :   2015-03-09                                                        !
!_______________________________________________________________________________________!
!
!  Author:
!       J. Milbrandt, McGill University (August 2004)
!
!  Major revisions:
!
!  001  J. Milbrandt  (Dec 2006) - Converted the full Milbrandt-Yau (2005) multimoment
!        (RPN)                     scheme to an efficient fixed-dispersion double-moment
!                                  version
!  002  J. Milbrandt  (Mar 2007) - Added options for single-moment/double-moment for
!                                  each hydrometeor category
!  003  J. Milbrandt  (Feb 2008) - Modified single-moment version for use in GEM-LAM-2.5
!  004  J. Milbrandt  (Nov 2008) - Modified double-moment version for use in 2010 Vancouver
!                                  Olympics GEM-LAM configuration
!  005  J. Milbrandt  (Aug 2009) - Modified (dmom) for PHY_v5.0.4, for use in V2010 system:
!                                  + reduced ice/snow capacitance to C=0.25D (from C=0.5D)
!                                  + added diagnostic fields (VIS, levels, etc.)
!                                  + added constraints to snow size distribution (No_s and
!                                    LAMDA_s limits, plus changed m-D parameters
!                                  + modified solid-to-liquid ratio calculation, based on
!                                    volume flux (and other changes)
!                                  + added back sedimentation of ice category
!                                  + modified condition for conversion of graupel to hail
!                                  + corrected bug it diagnostic "ice pellets" vs. "hail"
!                                  + minor bug corrections (uninitialized values, etc.)
!  006  J. Milbrandt  (Jan 2011) - Bug fixes and minor code clean-up from PHY_v5.1.3 version
!                                  + corrected latent heat constants in thermodynamic functions
!                                    (ABi and ABw) for sublimation and evaporation
!                                  + properly initialized variables No_g and No_h
!                                  + changed max ice crystal size (fallspeed) to 5 mm (2 m s-1)
!                                  + imposed maximum ice number concentration of 1.e+7 m-3
!                                  + removed unused supersaturation reduction
!  007  J. Milbrandt  (Nov 2012) - Bug fixes:
!                                  + suppress riming of snow, graupel, or hail if number
!                                    concentration is completed depleted
!                                  + removed extra QMLrg reduction in 'Adjustment of S/S'
!                                  + correction to limits in calculation of QCLsr and NCLir
!                                  + initialization of N{X} and N{X}M if dblMom_{x}=.false.
!  008  J. Milbrandt  (Jun 2013)  - Removed vertical motion unit conversion (now passed in as W)
!                                   Added limit on supersaturation w.r.t. ice for ice nucleation
!                                   Added radar reflectivity contribution from cloud water
!  009  J. Milbrandt  (Apr 2014)  - Bug fix to calculation of Dho (prevents divide-by-zero)
!  010  S. Gravel     (Oct 2014)  - Additional outputs for chemistry module
!  011  J. Milbrandt  (Mar 2015)  - Added user-specified temperature condition for 3-component freezing
!
!  Object:
!          Computes changes to the temperature, water vapor mixing ratio, and the
!          mixing ratios and total number concentrations of six hydrometeor species
!          resulting from cloud microphysical interactions at saturated grid points.
!          Liquid and solid surface precipitation rates from sedimenting hydrometeor
!          categories are also computed.
!
!          This subroutine and the associated modules form the single/double-moment
!          switchable verion of the multimoment bulk microphysics package, the full
!          version of which is described in the references below.
!
!  References:  Milbrandt and Yau, (2005a), J. Atmos. Sci., vol.62, 3051-3064
!               --------- and ---, (2005b), J. Atmos. Sci., vol.62, 3065-3081
!               (and references therein)
!
!  Please report bugs to:  jason.milbrandt@ec.gc.ca
!_______________________________________________________________________________________!
!
! Arguments:         Description:                                         Units:
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!            - Input -
!
! NI                 number of x-dir points (in local subdomain)
! NK                 number of vertical levels
! N                  not used (to be removed)
! J                  y-dir index (local subdomain)
! KOUNT              current model time step number
! dt                 model time step                                      [s]
! CCNtype            switch for airmass type
!                      1 = maritime                   --> N_c =  80 cm-3  (1-moment cloud)
!                      2 = continental 1              --> N_c = 200 cm-3     "       "
!                      3 = continental 2  (polluted)  --> N_c = 500 cm-3     "       "
!                      4 = land-sea-mask-dependent (TBA)
! WZ                 vertical velocity                                    [m s-1]
! S                  sigma (=p/p_sfc)
! GZ                 geopotential
! dblMom_(x)         logical switch for double(T)-single(F)-moment for category (x)
! precipDiag_ON      logical switch, .F. to suppress calc. of sfc precip types
! sedi_ON            logical switch, .F. to suppress sedimentation
! warmphase_ON       logical switch, .F. to suppress warm-phase (Part II)
! autoconv_ON        logical switch, .F. to supppress autoconversion (cld->rn)
! icephase_ON        logical switch, .F. to suppress ice-phase (Part I)
! snow_ON            logical switch, .F. to suppress snow initiation
! Tc_3comp           temperature (C) above which 3-component freezing is suppressed for
!                     Tc_3comp(1): ice  + rain
!                     Tc_3comp(2): snow + rain
!                     Tc_3comp(3): grpl + rain
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!            - Input/Output -
!
! T                  air temperature at time (t*)                         [K]
! TM                 air temperature at time (t-dt)                       [K]
! Q                  water vapor mixing ratio at (t*)                     [kg kg-1]
! QM                 water vapor mixing ratio at (t-dt)                   [kg kg-1]
! PS                 surface pressure at time (t*)                        [Pa]
! PSM                surface pressure at time (t-dt)                      [Pa]
!
!  For x = (C,R,I,N,G,H):  C = cloud
!                          R = rain
!                          I = ice (pristine) [except 'NY', not 'NI']
!                          N = snow
!                          G = graupel
!                          H = hail
!
! Q(x)               mixing ratio for hydrometeor x at (t*)               [kg kg-1]
! Q(x)M              mixing ratio for hydrometeor x at (t-dt)             [kg kg-1]
! N(x)               total number concentration for hydrometeor x  (t*)   [m-3]
! N(x)M              total number concentration for hydrometeor x  (t-dt) [m-3]
!
! Note:  The arrays "VM" (e.g. variables TM,QM,QCM etc.) are declared as INTENT(INOUT)
!        such that their values are modified in the code [VM = 0.5*(VM + V)].
!        This is to approxiate the values at time level (t), which are needed by
!        this routine but are unavailable to the PHYSICS.  The new values are discared
!        by the calling routine ('vkuocon6.ftn').  However, care should be taken with
!        interfacing with other modelling systems.  For GEM/MC2, it does not matter if
!        VM is modified since the calling module passes back only the tendencies
!        (VTEND) to the model.

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!            - Output -
!
! Q_TEND             tendency for water vapor mixing ratio                [kg kg-1 s-1]
! T_TEND             tendency for air temperature                         [K s-1]
! Q(x)TEND           tendency for mixing ratio for hydrometeor x          [kg kg-1 s-1]
! N(x)TEND           tendency for number concentration for hydrometeor x  [m-3 s-1]
! Dm_(x)             mean-mass diameter for hydrometeor x                 [m]
! H_CB               height of cloud base                                 [m]
! h_ML1              height of first melting level from ground            [m]
! h_ML2              height of first melting level from top               [m]
! h_SN               height of snow level                                 [m]
! RT_rn1             precipitation rate (at sfc) of liquid rain           [m+3 m-2 s-1]
! RT_rn2             precipitation rate (at sfc) of liquid drizzle        [m+3 m-2 s-1]
! RT_fr1             precipitation rate (at sfc) of freezing rain         [m+3 m-2 s-1]
! RT_fr2             precipitation rate (at sfc) of freezing drizzle      [m+3 m-2 s-1]
! RT_sn1             precipitation rate (at sfc) of ice crystals (liq-eq) [m+3 m-2 s-1]
! RT_sn2             precipitation rate (at sfc) of snow    (liq-equiv)   [m+3 m-2 s-1]
! RT_sn3             precipitation rate (at sfc) of graupel (liq-equiv)   [m+3 m-2 s-1]
! RT_snd             precipitation rate (at sfc) of snow    (frozen)      [m+3 m-2 s-1]
! RT_pe1             precipitation rate (at sfc) of ice pellets (liq-eq)  [m+3 m-2 s-1]
! RT_pe2             precipitation rate (at sfc) of hail (total; liq-eq)  [m+3 m-2 s-1]
! RT_peL             precipitation rate (at sfc) of hail (large only)     [m+3 m-2 s-1]
! SSxx               S/S terms (for testing purposes)
! SLW                supercooled liquid water content                     [kg m-3]
! VIS                visibility resulting from fog, rain, snow            [m]
! VIS1               visibility component through liquid cloud (fog)      [m]
! VIS2               visibility component through rain                    [m]
! VIS3               visibility component through snow                    [m]
! ZET                total equivalent radar reflectivity                  [dBZ]
! ZEC                composite (column-max) of ZET                        [dBZ]
!_______________________________________________________________________________________!


!LOCAL VARIABLES:

 !Variables to count active grid points:
  logical :: log1,log2,log3,log4,doneK,rainPresent,calcDiag,CB_found,ML_found,      &
             SN_found
  logical, dimension(size(QC,dim=1),size(QC,dim=2)) :: activePoint
  integer, dimension(size(QC,dim=1)) :: ktop_sedi
  integer :: i,k,niter,ll,start

  real    :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,                    &
       VDmax,NNUmax,X,D,DEL,QREVP,NuDEPSOR,NuCONTA,NuCONTB,NuCONTC,iMUkin,          &
       NuCONT,GG,Na,Tcc,F1,F2,Kdiff,PSIa,Kn,source,sink,sour,ratio,qvs0,            &
       DELqvs,ft,esi,Si,Simax,Vq,Vn,Vz,LAMr,No_r_DM,No_i,No_s,No_g,No_h,            &
       iABi,ABw,VENTr,VENTs,VENTg,VENTi,VENTh,Cdiff,Ka,MUdyn,MUkin,DEo,             &
       gam,ScTHRD,Tc,mi,ff,Ec,Ntr,Dho,DMrain,Ech,DMice,DMsnow,DMgrpl,DMhail,        &
       ssat,Swmax,dey,Esh,Eii,Eis,Ess,Eig,Eih,FRAC,JJ,Dirg,Dirh,Dsrs,Dsrg,Dsrh,     &
       Dgrg,Dgrh,SIGc,L,TAU,DrAUT,DrINIT,Di,Ds,Dg,Dh,qFact,nFact,Ki,Rz,             &
       vr0,vi0,vs0,vg0,vh0,Dc,Dr,QCLcs,QCLrs,QCLis,QCLcg,QCLrg,QCLig,QCNig,NCNig,   &
       QCLch,QCLrh,QCLsh,QMLir,QMLsr,QMLgr,QMLhr,QCLih,QVDvg,QVDvh,QSHhr,           &
       QFZci,QNUvi,QVDvi,QCNis,QCNis1,QCNis2,QCLir,QCLri,QCNsg,QCLsr,QCNgh,         &
       QCLgr,QHwet,QVDvs,QFZrh,QIMsi,QIMgi,NMLhr,NCNgh,NVDvh,NCLir,NCLri,NCLrh,     &
       NCLch,NCLsr,NCLirg,NCLirh,NrFZrh,NhFZrh,NCLsrs,NCLsrg,NCLsrh,NCLgrg,         &
       NCLgrh,NVDvg,NMLgr,NiCNis,NsCNis,NVDvs,NMLsr,NCLsh,NCLss,NNUvi,NFZci,NVDvi,  &
       NCLis,NCLig,NCLih,NMLir,NCLrs,NCNsg,NCLcs,NCLcg,NIMsi,NIMgi,NCLgr,NCLrg,     &
       NSHhr,RCAUTR,RCACCR,CCACCR,CCSCOC,CCAUTR,CRSCOR,ALFx,des_pmlt,Ecs,des,ides,  &
       LAMx,iLAMx,iLAMxB0,Dx,ffx,iLAMc,iNCM,iNRM,iNYM,iNNM,iNGM,iLAMs_D3,           &
       iLAMg,iLAMg2,iLAMgB0,iLAMgB1,iLAMgB2,iLAMh,iLAMhB0,iLAMhB1,iLAMhB2,iNHM,     &
       iLAMi,iLAMi2,iLAMi3,iLAMi4,iLAMi5,iLAMiB0,iLAMiB1,iLAMiB2,iLAMr6,iLAMh2,     &
       iLAMs,iLAMs2,iLAMsB0,iLAMsB1,iLAMsB2,iLAMr,iLAMr2,iLAMr3,iLAMr4,iLAMr5,      &
       iLAMc2,iLAMc3,iLAMc4,iLAMc5,iLAMc6,iQCM,iQRM,iQIM,iQNM,iQGM,iQHM,iEih,iEsh,  &
       N_c,N_r,N_i,N_s,N_g,N_h,fluxV_i,fluxV_g,fluxV_s,rhos_mlt,fracLiq,maxQlClip,  &
       maxQiClip,iNC,iNY,iNN,iNG,iNH,cloudliq,cloudice,cloudsnow

 !Variables that only need to be calulated on the first step (and saved):
  real, save :: idt,iMUc,cmr,cmi,cms,cmg,cmh,icmr,icmi,icmg,icms,icmh,idew,idei,    &
       ideh,ideg,GC1,imso,icexc9,cexr1,cexr2,cexr3,No_s_SM,No_r,idms,imgo,icexs2,   &
       cexr4,cexr5,cexr6,cexr9,icexr9,ckQr1,ckQr2,ckQr3,ckQi1,ckQi2,ckQi3,ckQi4,    &
       icexi9,ckQs1,ckQs2,cexs1,cexs2,ckQg1,ckQg2,ckQg4,ckQh1,ckQh2,ckQh4,GR37,dms, &
       LCP,LFP,LSP,ck5,ck6,PI2,PIov4,PIov6,CHLS,iCHLF,cxr,cxi,Gzr,Gzi,Gzs,Gzg,Gzh,  &
       N_c_SM,iGC1,GC2,GC3,GC4,GC5,iGC5,GC6,GC7,GC8,GC11,GC12,GC13,GC14,iGR34,mso,  &
       GC15,GR1,GR3,GR13,GR14,GR15,GR17,GR31,iGR31,GR32,GR33,GR34,GR35,GR36,GI4,    &
       GI6,GI20,GI21,GI22,GI31,GI32,GI33,GI34,GI35,iGI31,GI11,GI36,GI37,GI40,iGG34, &
       GS09,GS11,GS12,GS13,iGS20,GS31,iGS31,GS32,GS33,GS34,GS35,GS36,GS40,iGS40,    &
       GS50,GG09,GG11,GG12,GG13,GG31,iGG31,GG32,GG33,GG34,GG35,GG36,GG40,iGG99,GH09,&
       GH11,GH12,GH13,GH31,GH32,GH33,GH40,GR50,GG50,iGH34,GH50,iGH99,iGH31,iGS34,   &
       iGS20_D3,GS40_D3,cms_D3,eds,fds,rfact_FvFm,Gzc

!Size distribution parameters:
  real, parameter :: MUc      =  3.    !shape parameter for cloud
  real, parameter :: alpha_c  =  1.    !shape parameter for cloud
  real, parameter :: alpha_r  =  0.    !shape parameter for rain
  real, parameter :: alpha_i  =  0.    !shape parameter for ice
  real, parameter :: alpha_s  =  0.    !shape parameter for snow
  real, parameter :: alpha_g  =  0.    !shape parameter for graupel
  real, parameter :: alpha_h  =  0.    !shape parameter for hail
  real, parameter :: No_s_max =  1.e+8 !max. allowable intercept for snow [m-4]
  real, parameter :: lamdas_min= 500.  !min. allowable LAMDA_s [m-1]

 !For single-moment:
  real, parameter :: No_r_SM  =  1.e+7  !intercept parameter for rain    [m-4]
  real, parameter :: No_g_SM  =  4.e+6  !intercept parameter for graupel [m-4]
  real, parameter :: No_h_SM  =  1.e+5  !intercept parameter for hail    [m-4]
  !note: No_s = f(T), rather than a fixed value
  !------------------------------------!
  ! Symbol convention: (dist. params.) ! MY05: Milbrandt & Yau, 2005a,b (JAS)
  !       MY05    F94       CP00       ! F94:  Ferrier, 1994            (JAS)
  !       ------  --------  ------     ! CP00: Cohard & Pinty, 2000a,b  (QJGR)
  !       ALFx    ALPHAx    MUx-1      !
  !       MUx     (1)       ALPHAx     !
  !       ALFx+1  ALPHAx+1  MUx        !
  !------------------------------------!
  !  Note: The symbols for MU and ALPHA are REVERSED from that of CP2000a,b
  !        Explicit appearance of MUr = 1. has been removed.

  ! Fallspeed parameters:
  real, parameter :: afr=  149.100,  bfr= 0.5000   !Tripoloi and Cotton (1980)
  real, parameter :: afi=   71.340,  bfi= 0.6635   !Ferrier (1994)
  real, parameter :: afs=   11.720,  bfs= 0.4100   !Locatelli and Hobbs (1974)
  real, parameter :: afg=   19.300,  bfg= 0.3700   !Ferrier (1994)
  real, parameter :: afh=  206.890,  bfh= 0.6384   !Ferrier (1994)
 !options:
 !real, parameter :: afs=    8.996,  bfs= 0.4200   !Ferrier (1994)
 !real, parameter :: afg=   6.4800,  bfg= 0.2400   !LH74 (grpl-like snow of lump type)

  real, parameter :: epsQ  = 1.e-14   !kg kg-1, min. allowable mixing ratio
  real, parameter :: epsN  = 1.e-3    !m-3,     min. allowable number concentration
  real, parameter :: epsQ2 = 1.e-6    !kg kg-1, mixing ratio threshold for diagnostics
  real, parameter :: epsVIS= 1.       !m,       min. allowable visibility

  real, parameter :: iLAMmin1= 1.e-6  !min. iLAMx (prevents underflow in Nox and VENTx calcs)
  real, parameter :: iLAMmin2= 1.e-10 !min. iLAMx (prevents underflow in Nox and VENTx calcs)
  real, parameter :: eps   = 1.e-32
  real, parameter :: k1    = 0.001
  real, parameter :: k2    = 0.0005
  real, parameter :: k3    = 2.54
  real, parameter :: CPW   = 4218., CPI_=2093.

  real, parameter :: deg   =  400., mgo= 1.6e-10
  real, parameter :: deh   =  900.
  real, parameter :: dei   =  500., mio=1.e-12, Nti0=1.e3
  real, parameter :: dew   = 1000.
  real, parameter :: desFix=  100.  !used for snowSpherical = .true.
  real, parameter :: desMax=  500.
  real, parameter :: Dso   =  125.e-6  ![m]; embryo snow diameter (mean-volume particle)
  real, parameter :: dmr   = 3., dmi= 3., dmg= 3., dmh= 3.

  ! NOTE: VxMAX below are the max.allowable mass-weighted fallspeeds for sedimentation.
  !       Thus, Vx corresponds to DxMAX (at sea-level) times the max. density factor, GAM.
  !       [GAMmax=sqrt(DEo/DEmin)=sqrt(1.25/0.4)~2.]  e.g. VrMAX = 2.*8.m/s = 16.m/s
  real, parameter :: DrMax=  5.e-3,   VrMax= 16.,   epsQr_sedi= 1.e-8
  real, parameter :: DiMax=  5.e-3,   ViMax=  2.,   epsQi_sedi= 1.e-10
  real, parameter :: DsMax=  5.e-3,   VsMax=  2.,   epsQs_sedi= 1.e-8
  real, parameter :: DgMax= 50.e-3,   VgMax=  8.,   epsQg_sedi= 1.e-8
  real, parameter :: DhMax= 80.e-3,   VhMax= 25.,   epsQh_sedi= 1.e-10

  real, parameter :: thrd    = 1./3.
  real, parameter :: sixth   = 0.5*thrd
  real, parameter :: Ecg     = 1., Ers= 1., Eci= 1. !collection efficiencies, Exy, between categories x and y
  real, parameter :: Eri     = 1., Erg= 1., Erh= 1.
  real, parameter :: Xdisp   = 0.25               !dispersion of the fall velocity of ice
  real, parameter :: aa11    = 9.44e15, aa22= 5.78e3, Rh= 41.e-6
  real, parameter :: Avx     = 0.78, Bvx= 0.30    !ventilation coefficients [F94 (B.36)]
  real, parameter :: Abigg   = 0.66, Bbigg= 100.  !parameters in probabilistic freezing
  real, parameter :: fdielec     = 4.464          !ratio of dielectric factor, |K|w**2/|K|i**2
  real, parameter :: zfact       = 1.e+18         !conversion factor for m-3 to mm2 m-6 for Ze
  real, parameter :: minZET      = -99.           ![dBZ] min threshold for ZET
  real, parameter :: maxVIS      = 99.e+3         ![m] max. allowable VIS (visibility)
  real, parameter :: Drshed      = 0.001          ![m] mean diam. of drop shed during wet growth
  real, parameter :: SIGcTHRS    = 15.e-6         !threshold cld std.dev. before autoconversion
  real, parameter :: KK1         = 3.03e3         !parameter in Long (1974) kernel
  real, parameter :: KK2         = 2.59e15        !parameter in Long (1974) kernel
  real, parameter :: Dhh         = 82.e-6         ![m] diameter that rain hump first appears
  real, parameter :: gzMax_sedi  = 180000.        !GZ value below which sedimentation is computed
  real, parameter :: Dr_large    = 200.e-6        ![m] size threshold to distinguish rain/drizzle for precip rates
  real, parameter :: Ds_large    = 200.e-6        ![m] size threshold to distinguish snow/snow-grains for precip rates
  real, parameter :: Dh_large    = 1.0e-2         ![m] size threshold for "large" hail precipitation rate
  real, parameter :: Dh_min      = 5.0e-3         ![m] size threhsold for below which hail converts to graupel
  real, parameter :: Dr_3cmpThrs = 2.5e-3         ![m] size threshold for hail production from 3-comp freezing
  real, parameter :: w_CNgh      = 3.             ![m s-1] vertical motion  threshold for CNgh
  real, parameter :: r_CNgh      = 0.05           !Dg/Dho ratio threshold for CNgh
  real, parameter :: Tc_FZrh     = -10.           !temp-threshold (C) for FZrh
  real, parameter :: CNsgThres   = 1.0            !threshold for CLcs/VDvs ratio for CNsg
  real, parameter :: capFact_i   = 0.5            !capacitace factor for ice  (C= 0.5*D*capFact_i)
  real, parameter :: capFact_s   = 0.5            !capacitace factor for snow (C= 0.5*D*capFact_s)
  real, parameter :: noVal_h_XX  = -1.            !non-value indicator for h_CB, h_ML1, h_ML2, h_SN
  real, parameter :: minSnowSize = 1.e-4          ![m] snow size threshold to compute h_SN
  real, parameter :: Fv_Dsmin    = 125.e-6        ![m] min snow size to compute volume flux
  real, parameter :: Fv_Dsmax    = 0.008          ![m] max snow size to compute volume flux
  real, parameter :: Ni_max      = 1.e+7          ![m-3] max ice crystal concentration
  real, parameter :: maxTclip    = 10.            ![k] max allowable drop in T from Q clip
  real, parameter :: satw_peak   = 1.01           !assumed max. peak saturation w.r.t. water (for calc of Simax)

  ! Constants used for contact ice nucleation:
  real, parameter :: LAMa0  = 6.6e-8     ![m] mean free path at T0 and p0 [W95_eqn58]
  real, parameter :: T0     = 293.15     ![K] ref. temp.
  real, parameter :: p0     = 101325.    ![Pa] ref. pres.
  real, parameter :: Ra     = 1.e-6      ![m] aerosol (IN) radius         [M92 p.713; W95_eqn60]
  real, parameter :: kBoltz = 1.381e-23  !Boltzmann's constant
  real, parameter :: KAPa   = 5.39e5     !aerosol thermal conductivity

 !Test switches:
  logical, parameter :: iceDep_ON     = .true.  !.false. to suppress depositional growth of ice
  logical, parameter :: grpl_ON       = .true.  !.false. to suppress graupel initiation
  logical, parameter :: hail_ON       = .true.  !.false. to suppress hail initiation
  logical, parameter :: rainAccr_ON   = .true.  ! rain accretion and self-collection ON/OFF
  logical, parameter :: snowSpherical = .false. !.true.: m(D)=(pi/6)*const_des*D^3 | .false.: m(D)= 0.069*D^2
  integer, parameter :: primIceNucl   = 1       !1= Meyers+contact ;  2= Cooper
  real,    parameter :: outfreq       =  60.    !frequency to compute output diagnostics [s]

!Passed as physics namelist parameters:
! logical, parameter :: precipDiag_ON = .true.  !.false. to suppress calc. of sfc precip types
! logical, parameter :: sedi_ON       = .true.  !.false. to suppress sedimentation
! logical, parameter :: warmphase_ON  = .true.  !.false. to suppress warm-phase (Part II)
! logical, parameter :: autoconv_ON   = .true.  ! autoconversion ON/OFF
! logical, parameter :: icephase_ON   = .true.  !.false. to suppress ice-phase (Part I)
! logical, parameter :: snow_ON       = .true.  !.false. to suppress snow initiation
! logical, parameter :: initN         = .true.  !.true.  to initialize Nx of Qx>0 and Nx=0

  real, dimension(size(QC,dim=1),size(QC,dim=2)) :: DE,iDE,DP,QSS,QSW,QSI,DZ,RHOQX,VQQ,      &
        gamfact,gamfact_r,massFlux3D_r,massFlux3D_s,massFlux3D_i,massFlux3D_g,massFlux3D_h
  real, dimension(size(QC,dim=1))                :: fluxM_r,fluxM_i,fluxM_s,fluxM_g,fluxM_h, &
        HPS,dum
  integer, dimension(size(QC,dim=1))             :: activeColumn

  !==================================================================================!

  !----------------------------------------------------------------------------------!
  !                      PART 1:   Prelimiary Calculations                           !
  !----------------------------------------------------------------------------------!

  !---
  if (.not.dblMom_c) then
     NCM= 0.;  NC = 0.
  endif
  if (.not.dblMom_r) then
     NRM= 0.;  NR = 0.
  endif
  if (.not.dblMom_i) then
     NYM= 0.;  NY = 0.
  endif
  if (.not.dblMom_s) then
     NNM= 0.;  NN = 0.
  endif
  if (.not.dblMom_g) then
     NGM= 0.;  NG = 0.
  endif
  if (.not.dblMom_h) then
     NHM= 0.;  NH = 0.
  endif
  !===

 !-------------
 !Convert N from #/kg to #/m3:
  do k= 1,nk
    do i= 1,ni
      tmp1= S(i,k)*PSM(i)/(RGASD*TM(i,k))  !air density at time (t-1)
      tmp2= S(i,k)*PS(i)/(RGASD*T(i,k))    !air density at time (*)

      NCM(i,k)= NCM(i,k)*tmp1;   NC(i,k)= NC(i,k)*tmp2
      NRM(i,k)= NRM(i,k)*tmp1;   NR(i,k)= NR(i,k)*tmp2
      NYM(i,k)= NYM(i,k)*tmp1;   NY(i,k)= NY(i,k)*tmp2
      NNM(i,k)= NNM(i,k)*tmp1;   NN(i,k)= NN(i,k)*tmp2
      NGM(i,k)= NGM(i,k)*tmp1;   NG(i,k)= NG(i,k)*tmp2
      NHM(i,k)= NHM(i,k)*tmp1;   NH(i,k)= NH(i,k)*tmp2
    enddo
  enddo
 !=============

  ! The SSxx arrays are for passed to the volatile bus for output as 3-D diagnostic
  ! output variables, for testing purposes.  For example, to output the
  ! instantanous value of the deposition rate, add 'SS01(i,k) = QVDvi'  in the
  ! appropriate place.  It can then be output as a 3-D physics variable by adding
  ! it to the sortie_p list in 'outcfgs.out'

  SS01= 0.; SS02= 0.; SS03= 0.; SS04= 0.; SS05= 0.; SS06= 0.; SS07= 0.; SS08= 0.
  SS09= 0.; SS10= 0.; SS11= 0.; SS12= 0.; SS13= 0.; SS14= 0.; SS15= 0.; SS16= 0.
  SS17= 0.; SS18= 0.; SS19= 0.; SS20= 0.

 !Determine the upper-most level in each column to which to compute sedimentation:
  ktop_sedi= 0
  do i=1,ni
     do k=1,nk
       ktop_sedi(i)= k
       if (GZ(i,k)<gzMax_sedi) exit
     enddo
  enddo

 !Compute diagnostic values only every 'outfreq' minutes:
 !calcDiag= (mod(DT*float(KOUNT),outfreq)==0.)
  calcDiag = .true.  !compute diagnostics every step (for time-series output)

!####  These need only to be computed once per model integration:
!      (note:  These variables must be declared with the SAVE attribute)

! if (KOUNT==0) then
!*** For restarts, these values are not saved.  Therefore, the condition statement
!    must be modified to something like: IF (MOD(Step_rsti,KOUNT).eq.0) THEN
!    in order that these be computed only on the first step of a given restart.
!    (...to be done.  For now, changing condition to IF(TRUE) to compute at each step.)


  if (.TRUE.) then

   PI2    = PI*2.
   PIov4  = 0.25*PI
   PIov6  = PI*sixth
   CHLS   = CHLC+CHLF  !J k-1; latent heat of sublimation
   LCP    = CHLC/CPD
   LFP    = CHLF/CPD
   iCHLF  = 1./CHLF
   LSP    = LCP+LFP
   ck5    = 4098.170*LCP
   ck6    = 5806.485*LSP
   idt    = 1./dt
   imgo   = 1./mgo
   idew   = 1./dew
   idei   = 1./dei
   ideg   = 1./deg
   ideh   = 1./deh
   maxQlClip = maxTclip/LCP  !kg kg-1; max liquid mass allowed to clip to vapor
   maxQiClip = maxTclip/LSP  !kg kg-1; max frozen mass allowed to clip to vapor

   !Constants based on size distribution parameters:

   ! Mass parameters [ m(D) = cD^d ]
   cmr    = PIov6*dew;  icmr= 1./cmr
   cmi    = 440.;       icmi= 1./cmi
   cmg    = PIov6*deg;  icmg= 1./cmg
   cmh    = PIov6*deh;  icmh= 1./cmh

   cms_D3 = PIov6*desFix !used for snowSpherical = .T. or .F.
   if (snowSpherical) then
      cms = cms_D3
      dms = 3.
   else
!     cms = 0.0690;  dms = 2.000   !Cox, 1988 (QJRMS)
      cms = 0.1597;  dms = 2.078   !Brandes et al., 2007 (JAMC)
   endif
   icms   = 1./cms
   idms   = 1./dms
   mso    = cms*Dso**dms
   imso   = 1./mso
  !bulk density parameters: [rho(D) = eds*D^fds]
  !  These are implied by the mass-diameter parameters, by computing the bulk
  !  density of a sphere with the equaivalent mass.
  !  e.g. m(D) = cD^d = (pi/6)rhoD^3 and solve for rho(D)
   eds    = cms/PIov6
   fds    = dms-3.
   if (fds/=-1. .and..not.snowSpherical) GS50= gamma(1.+fds+alpha_s)

   ! Cloud:
   iMUc   =  1./MUc
   GC1    =  gamma(alpha_c+1.0)
   iGC1   = 1./GC1
   GC2    =  gamma(alpha_c+1.+3.0*iMUc)  !i.e. gamma(alf + 4)
   GC3    =  gamma(alpha_c+1.+6.0*iMUc)  !i.e. gamma(alf + 7)
   GC4    =  gamma(alpha_c+1.+9.0*iMUc)  !i.e. gamma(alf + 10)
   GC11   =  gamma(1.0*iMUc+1.0+alpha_c)
   GC12   =  gamma(2.0*iMUc+1.0+alpha_c)
   GC5    =  gamma(1.0+alpha_c)
   iGC5   = 1./GC5
   GC6    =  gamma(1.0+alpha_c+1.0*iMUc)
   GC7    =  gamma(1.0+alpha_c+2.0*iMUc)
   GC8    =  gamma(1.0+alpha_c+3.0*iMUc)
   GC13   =  gamma(3.0*iMUc+1.0+alpha_c)
   GC14   =  gamma(4.0*iMUc+1.0+alpha_c)
   GC15   =  gamma(5.0*iMUc+1.0+alpha_c)
   icexc9 =  1./(GC2*iGC1*PIov6*dew)
  !specify cloud droplet number concentration [m-3] based on 'CCNtype' (1-moment):
   if     (CCNtype==1) then
      N_c_SM =  0.8e+8          !maritime
   elseif (CCNtype==2) then
      N_c_SM =  2.0e+8          !continental 1
   elseif (CCNtype==3) then
      N_c_SM =  5.0e+8          !continental 2 (polluted)
   else
      N_c_SM =  2.0e+8          !default (cont1), if 'CCNtype' specified incorrectly
   endif

   ! Rain:
   cexr1  = 1.+alpha_r+dmr+bfr
   cexr2  = 1.+alpha_r+dmr
   GR17   = gamma(2.5+alpha_r+0.5*bfr)
   GR31   = gamma(1.+alpha_r)
   iGR31  = 1./GR31
   GR32   = gamma(2.+alpha_r)
   GR33   = gamma(3.+alpha_r)
   GR34   = gamma(4.+alpha_r)
   iGR34  = 1./GR34
   GR35   = gamma(5.+alpha_r)
   GR36   = gamma(6.+alpha_r)
   GR37   = gamma(7.+alpha_r)
   GR50   = (No_r_SM*GR31)**0.75  !for 1-moment or Nr-initialization
   cexr5  = 2.+alpha_r
   cexr6  = 2.5+alpha_r+0.5*bfr
   cexr9  = cmr*GR34*iGR31;    icexr9= 1./cexr9
   cexr3  = 1.+bfr+alpha_r
   cexr4  = 1.+alpha_r
   ckQr1  = afr*gamma(1.+alpha_r+dmr+bfr)/gamma(1.+alpha_r+dmr)
   ckQr2  = afr*gamma(1.+alpha_r+bfr)*GR31
   ckQr3  = afr*gamma(7.+alpha_r+bfr)/GR37
   if (.not.dblMom_r) then
      No_r = No_r_SM
   endif

   ! Ice:
   GI4    = gamma(alpha_i+dmi+bfi)
   GI6    = gamma(2.5+bfi*0.5+alpha_i)
   GI11   = gamma(1.+bfi+alpha_i)
   GI20   = gamma(0.+bfi+1.+alpha_i)
   GI21   = gamma(1.+bfi+1.+alpha_i)
   GI22   = gamma(2.+bfi+1.+alpha_i)
   GI31   = gamma(1.+alpha_i)
   iGI31  = 1./GI31
   GI32   = gamma(2.+alpha_i)
   GI33   = gamma(3.+alpha_i)
   GI34   = gamma(4.+alpha_i)
   GI35   = gamma(5.+alpha_i)
   GI36   = gamma(6.+alpha_i)
   GI40   = gamma(1.+alpha_i+dmi)
   icexi9 = 1./(cmi*gamma(1.+alpha_i+dmi)*iGI31)
   ckQi1  = afi*gamma(1.+alpha_i+dmi+bfi)/GI40
   ckQi2  = afi*GI11*iGI31
   ckQi4  = 1./(cmi*GI40*iGI31)

   ! Snow:
   cexs1  = 2.5+0.5*bfs+alpha_s
   cexs2  = 1.+alpha_s+dms
   icexs2 = 1./cexs2
   GS09   = gamma(2.5+bfs*0.5+alpha_s)
   GS11   = gamma(1.+bfs+alpha_s)
   GS12   = gamma(2.+bfs+alpha_s)
   GS13   = gamma(3.+bfs+alpha_s)
   GS31   = gamma(1.+alpha_s)
   iGS31  = 1./GS31
   GS32   = gamma(2.+alpha_s)
   GS33   = gamma(3.+alpha_s)
   GS34   = gamma(4.+alpha_s)
   iGS34  = 1./GS34
   GS35   = gamma(5.+alpha_s)
   GS36   = gamma(6.+alpha_s)
   GS40   = gamma(1.+alpha_s+dms)
   iGS40  = 1./GS40
   iGS20  = 1./(GS40*iGS31*cms)
   ckQs1  = afs*gamma(1.+alpha_s+dms+bfs)*iGS40
   ckQs2  = afs*GS11*iGS31
   GS40_D3 = gamma(1.+alpha_s+3.)
   iGS20_D3= 1./(GS40_D3*iGS31*cms_D3)
   rfact_FvFm= PIov6*icms*gamma(4.+bfs+alpha_s)/gamma(1.+dms+bfs+alpha_s)

   ! Graupel:
   GG09   = gamma(2.5+0.5*bfg+alpha_g)
   GG11   = gamma(1.+bfg+alpha_g)
   GG12   = gamma(2.+bfg+alpha_g)
   GG13   = gamma(3.+bfg+alpha_g)
   GG31   = gamma(1.+alpha_g)
   iGG31  = 1./GG31
   GG32   = gamma(2.+alpha_g)
   GG33   = gamma(3.+alpha_g)
   GG34   = gamma(4.+alpha_g)
   iGG34  = 1./GG34
   GG35   = gamma(5.+alpha_g)
   GG36   = gamma(6.+alpha_g)
   GG40   = gamma(1.+alpha_g+dmg)
   iGG99  = 1./(GG40*iGG31*cmg)
   GG50   = (No_g_SM*GG31)**0.75     !for 1-moment only
   ckQg1  = afg*gamma(1.+alpha_g+dmg+bfg)/GG40
   ckQg2  = afg*GG11*iGG31
   ckQg4  = 1./(cmg*GG40*iGG31)

   ! Hail:
   GH09   = gamma(2.5+bfh*0.5+alpha_h)
   GH11   = gamma(1.+bfh+alpha_h)
   GH12   = gamma(2.+bfh+alpha_h)
   GH13   = gamma(3.+bfh+alpha_h)
   GH31   = gamma(1.+alpha_h)
   iGH31  = 1./GH31
   GH32   = gamma(2.+alpha_h)
   GH33   = gamma(3.+alpha_h)
   iGH34  = 1./gamma(4.+alpha_h)
   GH40   = gamma(1.+alpha_h+dmh)
   iGH99  = 1./(GH40*iGH31*cmh)
   GH50   = (No_h_SM*GH31)**0.75     !for 1-moment only
   ckQh1  = afh*gamma(1.+alpha_h+dmh+bfh)/GH40
   ckQh2  = afh*GH11*iGH31
   ckQh4  = 1./(cmh*GH40*iGH31)

  endif  !if (KOUNT=0)
!####

!=======================================================================================!

  !Compute thickness of layers for sedimentation calcuation:
  ! (note; 'GZ' passed in is geopotential, not geopotential height)
  tmp1= 1./GRAV
  do k=2,nk
     DZ(:,k)= (GZ(:,k-1)-GZ(:,k))*tmp1
  enddo
  DZ(:,1)= DZ(:,2)

! Temporarily store arrays at time (t*) in order to compute (at the end of subroutine)
! the final VXTEND as VXTEND = ( VX{t+1} - VX{t*} )/dt :
  T_TEND = T ;  Q_TEND = Q
  QCTEND = QC;  QRTEND = QR;  QITEND = QI;  QNTEND = QN;  QGTEND = QG;  QHTEND = QH
  NCTEND = NC;  NRTEND = NR;  NYTEND = NY;  NNTEND = NN;  NGTEND = NG;  NHTEND = NH

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Initialize Nx if Qx>0 and Nx=0:  (for nesting from 1-moment to 2-moment):
  IF (initN) THEN
     do k= 1,nk
        do i= 1,ni
           tmp1= S(i,k)*PSM(i)/(RGASD*TM(i,k))  !air density at time (t-1)
           tmp2= S(i,k)*PS(i)/(RGASD*T(i,k))    !air density at time (*)

        !cloud:
           if (QCM(i,k)>epsQ .and. NCM(i,k)<epsN)                                        &
              NCM(i,k)= N_c_SM
           if (QC(i,k)>epsQ  .and. NC(i,k)<epsN)                                         &
              NC(i,k) = N_c_SM
        !rain
           if (QRM(i,k)>epsQ .and. NRM(i,k)<epsN)                                        &
              NRM(i,k)= (No_r_SM*GR31)**(3./(4.+alpha_r))*(GR31*iGR34*tmp1*QRM(i,k)*     &
                        icmr)**((1.+alpha_r)/(4.+alpha_r))
           if (QR(i,k)>epsQ  .and. NR(i,k)<epsN)                                         &
              NR(i,k)= (No_r_SM*GR31)**(3./(4.+alpha_r))*(GR31*iGR34*tmp2*QR(i,k)*       &
                       icmr)**((1.+alpha_r)/(4.+alpha_r))
        !ice:
           if (QIM(i,k)>epsQ .and. NYM(i,k)<epsN)                                        &
              NYM(i,k)= N_Cooper(TRPL,TM(i,k))
           if (QI(i,k)>epsQ  .and. NY(i,k)<epsN)                                         &
              NY(i,k)= N_Cooper(TRPL,T(i,k))
        !snow:
           if (QNM(i,k)>epsQ .and. NNM(i,k)<epsN) then
              No_s= Nos_Thompson(TRPL,TM(i,k))
              NNM(i,k)= (No_s*GS31)**(dms*icexs2)*(GS31*iGS40*icms*tmp1*QNM(i,k))**      &
                        ((1.+alpha_s)*icexs2)
           endif
           if (QN(i,k)>epsQ  .and. NN(i,k)<epsN)  then
              No_s= Nos_Thompson(TRPL,T(i,k))
              NN(i,k)= (No_s*GS31)**(dms*icexs2)*(GS31*iGS40*icms*tmp2*QN(i,k))**        &
                       ((1.+alpha_s)*icexs2)
           endif
        !grpl:
           if (QGM(i,k)>epsQ .and. NGM(i,k)<epsN)                                        &
              NGM(i,k)= (No_g_SM*GG31)**(3./(4.+alpha_g))*(GG31*iGG34*tmp1*QGM(i,k)*     &
                        icmg)**((1.+alpha_g)/(4.+alpha_g))
           if (QG(i,k)>epsQ  .and. NG(i,k)<epsN)                                         &
              NG(i,k)= (No_g_SM*GG31)**(3./(4.+alpha_g))*(GG31*iGG34*tmp2*QG(i,k)*       &
                   icmg)**((1.+alpha_g)/(4.+alpha_g))
        !hail:
           if (QHM(i,k)>epsQ .and. NHM(i,k)<epsN)                                        &
              NHM(i,k)= (No_h_SM*GH31)**(3./(4.+alpha_h))*(GH31*iGH34*tmp1*QHM(i,k)*     &
                        icmh)**((1.+alpha_h)/(4.+alpha_h))
           if (QH(i,k)>epsQ  .and. NH(i,k)<epsN)                                         &
              NH(i,k)= (No_h_SM*GH31)**(3./(4.+alpha_h))*(GH31*iGH34*tmp2*QH(i,k)*       &
                        icmh)**((1.+alpha_h)/(4.+alpha_h))

        enddo !i-loop
     enddo    !k-loop
  ENDIF  !N-initialization

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! Clip all moments to zero if one or more corresponding category moments are less than
!  the minimum allowable value:
! (Note: Clipped mass is added back to water vapor field to conserve total mass)
  do k= 1,nk
     do i= 1,ni

       IF (dblMom_c) THEN
         if(QC(i,k)<epsQ .or. NC(i,k)<epsN)    then
            Q(i,k) = Q(i,k) + QC(i,k)
            QC(i,k)= 0.;   NC(i,k)= 0.
         endif
         if(QCM(i,k)<epsQ .or. NCM(i,k)<epsN)  then
            QM(i,k) = QM(i,k) + QCM(i,k)
            QCM(i,k)= 0.;  NCM(i,k)= 0.
         endif
       ELSE
         if(QC(i,k)<epsQ)    then
            Q(i,k) = Q(i,k) + QC(i,k)
            QC(i,k)= 0.
         endif
         if(QCM(i,k)<epsQ)  then
            QM(i,k) = QM(i,k) + QCM(i,k)
            QCM(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_r) THEN
         if (QR(i,k)<epsQ .or. NR(i,k)<epsN)   then
            Q(i,k) = Q(i,k) + QR(i,k)
            QR(i,k)= 0.;  NR(i,k)= 0.
         endif
         if (QRM(i,k)<epsQ .or. NRM(i,k)<epsN) then
            QM(i,k) = QM(i,k) + QRM(i,k)
            QRM(i,k)= 0.;  NRM(i,k)= 0.
         endif
       ELSE
         if (QR(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QR(i,k)
            QR(i,k)= 0.
         endif
         if (QRM(i,k)<epsQ) then
            QM(i,k) = QM(i,k) + QRM(i,k)
            QRM(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_i) THEN
         if (QI(i,k)<epsQ .or. NY(i,k)<epsN)   then
            Q(i,k) = Q(i,k) + QI(i,k)
            QI(i,k)= 0.;  NY(i,k)= 0.
         endif
         if (QIM(i,k)<epsQ .or. NYM(i,k)<epsN) then
            QM(i,k) = QM(i,k) + QIM(i,k)
            QIM(i,k)= 0.;  NYM(i,k)= 0.
         endif
       ELSE
         if (QI(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QI(i,k)
            QI(i,k)= 0.
         endif
         if (QIM(i,k)<epsQ) then
            QM(i,k) = QM(i,k) + QIM(i,k)
            QIM(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_s) THEN
         if (QN(i,k)<epsQ .or. NN(i,k)<epsN)   then
            Q(i,k) = Q(i,k) + QN(i,k)
            QN(i,k)= 0.;  NN(i,k)= 0.
         endif
         if (QNM(i,k)<epsQ .or. NNM(i,k)<epsN) then
            QM(i,k) = QM(i,k) + QNM(i,k)
            QNM(i,k)= 0.;  NNM(i,k)= 0.
         endif
       ELSE
         if (QN(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QN(i,k)
            QN(i,k)= 0.
         endif
         if (QNM(i,k)<epsQ) then
            QM(i,k) = QM(i,k) + QNM(i,k)
            QNM(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_g) THEN
         if (QG(i,k)<epsQ .or. NG(i,k)<epsN)   then
            Q(i,k) = Q(i,k) + QG(i,k)
            QG(i,k)= 0.;  NG(i,k)= 0.
         endif
         if (QGM(i,k)<epsQ  .or. NGM(i,k)<epsN) then
            QM(i,k) = QM(i,k) + QGM(i,k)
            QGM(i,k)= 0.;  NGM(i,k)= 0.
         endif
       ELSE
         if (QG(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QG(i,k)
            QG(i,k)= 0.
         endif
         if (QGM(i,k)<epsQ) then
            QM(i,k) = QM(i,k) + QGM(i,k)
            QGM(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_h) THEN
         if (QH(i,k)<epsQ .or. NH(i,k)<epsN)   then
            Q(i,k) = Q(i,k) + QH(i,k)
            QH(i,k)= 0.;  NH(i,k)= 0.
         endif
         if (QHM(i,k)<epsQ .or. NHM(i,k)<epsN) then
            QM(i,k) = QM(i,k) + QHM(i,k)
            QHM(i,k)= 0.;  NHM(i,k)= 0.
         endif
       ELSE
         if (QH(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QH(i,k)
            QH(i,k)= 0.
         endif
         if (QHM(i,k)<epsQ) then
            QM(i,k) = QM(i,k) + QHM(i,k)
            QHM(i,k)= 0.
         endif
       ENDIF

    enddo  !i-loop
  enddo    !k-loop;    (clipping)
  QM = max(QM,0.)
  Q  = max(Q ,0.)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

! Approximate values at time {t}:
!  [ ave. of values at {*} (advected, but no physics tendency added) and {t-dt} ]:
  HPS= 0.5*(PSM+PS);   TM = 0.5*(TM + T);   QM = 0.5*(QM + Q)
  QCM= 0.5*(QCM+QC);   QRM= 0.5*(QRM+QR);   QIM= 0.5*(QIM+QI)
  QNM= 0.5*(QNM+QN);   QGM= 0.5*(QGM+QG);   QHM= 0.5*(QHM+QH)

  if (dblMom_c) NCM= 0.5*(NCM+NC)
  if (dblMom_r) NRM= 0.5*(NRM+NR)
  if (dblMom_i) NYM= 0.5*(NYM+NY)
  if (dblMom_s) NNM= 0.5*(NNM+NN)
  if (dblMom_g) NGM= 0.5*(NGM+NG)
  if (dblMom_h) NHM= 0.5*(NHM+NH)

  do k=1,nk
     do i=1,ni
      !Saturation mixing ratios:  [used in both Parts I and II]
        QSW(i,k)= sngl(FOQSA(TM(i,k),HPS(i)*S(i,k)))      !wrt. liquid water at (t)
        QSS(i,k)= sngl(FOQST( T(i,k), PS(i)*S(i,k)))      !wrt. ice surface  at (*)
        QSI(i,k)= sngl(FOQST(TM(i,k),HPS(i)*S(i,k)))      !wrt. ice surface  at (t)
      !Air density at time (t)
        DE(i,k) = S(i,k)*HPS(i)/(RGASD*TM(i,k))           !air density at time  (t)
        iDE(i,k)= 1./DE(i,k)
     enddo
  enddo

  do i= 1,ni
    !Air-density factor: (for fall velocity computations)
     DEo           = DE(i,nk)
     gamfact(i,:)  = sqrt(DEo/(DE(i,:)))
     gamfact_r(i,:)= sqrt( 1./(DE(i,:)))
  enddo

  !-- Diagnostics for chemistry module:
  PRFLXKY = 0.
  SWFLXKY = 0.
  F12KY   = 0.
  FEVPKY  = 0.
  !==

  !----------------------------------------------------------------------------------!
  !                 End of Preliminary Calculation section (Part 1)                  !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !                      PART 2: Cold Microphysics Processes                         !
  !----------------------------------------------------------------------------------!

! Determine the active grid points (i.e. those which scheme should treat):
  activePoint = .false.
  DO k=2,nk
     DO i=1,ni
        log1= ((QIM(i,k)+QGM(i,k)+QNM(i,k)+QHM(i,k))<epsQ)  !no solid  (i,g,s,h)
        log2= ((QCM(i,k)+QRM(i,k))                  <epsQ)  !no liquid (c,r)
        log3= ((TM(i,k)>TRPL) .and. log1)                   !T>0C & no i,g,s,h
        log4= log1.and.log2.and.(QM(i,k)<QSI(i,k))          !no sol. or liq.; subsat(i)
        if (.not.( log3 .or. log4 ) .and. icephase_ON) then
          activePoint(i,k)= .true.
        endif
     ENDDO
  ENDDO

    ! Size distribution parameters:
    !  Note: + 'thrd' should actually be '1/dmx'(but dmx=3 for all categories x)
    !        + If Qx=0, LAMx etc. are never be used in any calculations
    !          (If Qc=0, CLcy etc. will never be calculated. iLAMx is set to 0
    !           to avoid possible problems due to bugs.)

  DO k= 2,nk                       !Main loop for Part 2
    DO i= 1,ni
      IF (activePoint(i,k)) THEN


       Tc= TM(i,k)-TRPL
!       if (Tc<-120. .or. Tc>50.)   &
!        print*, '***WARNING*** -- In MICROPHYSICS --  Ambient Temp.(C):',Tc
       Cdiff = (2.2157e-5+0.0155e-5*Tc)*1.e5/(S(i,k)*HPS(i))
       MUdyn = 1.72e-5*(393./(TM(i,k)+120.))*(TM(i,k)/TRPL)**1.5 !RYp.102
       MUkin = MUdyn*iDE(i,k)
       iMUkin= 1./MUkin
       ScTHRD= (MUkin/Cdiff)**thrd       ! i.e. Sc^(1/3)
       Ka    = 2.3971e-2 + 0.0078e-2*Tc                                   !therm.cond.(air)
       Kdiff = (9.1018e-11*TM(i,k)*TM(i,k)+8.8197e-8*TM(i,k)-(1.0654e-5)) !therm.diff.(air)
       gam   = gamfact(i,k)

      !Collection efficiencies:
      !note:  Ecg=Eri=Ers=Erg=Erh=1. (constant parameters)
      !       - Ecs is computed in CLcs section
      !       - Ech is computed in CLch section
       Eis   = min(0.05*exp(0.1*Tc),1.)     !Ferrier, 1995 (Table 1)
       Eig   = min(0.01*exp(0.1*Tc),1.)     !dry (Eig=1.0 for wet growth)
       Eii   = 0.1*Eis
       Ess   = Eis;   Eih = Eig;   Esh = Eig
       iEih  = 1./Eih
       iEsh  = 1./Esh

       qvs0  = sngl(FOQSA(TRPL,HPS(i)*S(i,k)))      !sat.mix.ratio at 0C
       DELqvs= qvs0-(QM(i,k))

    ! Cloud:
       if (QCM(i,k)>epsQ) then
          if (.not. dblMom_c) NCM(i,k)= N_c_SM
          iQCM   = 1./QCM(i,k)
          iNCM   = 1./NCM(i,k)
          Dc     = Dm_x(DE(i,k),QCM(i,k),iNCM,icmr,thrd)

          iLAMc  = iLAMDA_x(DE(i,k),QCM(i,k),iNCM,icexc9,thrd)
          iLAMc2 = iLAMc *iLAMc
          iLAMc3 = iLAMc2*iLAMc
          iLAMc4 = iLAMc2*iLAMc2
          iLAMc5 = iLAMc3*iLAMc2
       else
          Dc     = 0.;   iLAMc3= 0.
          iLAMc  = 0.;   iLAMc4= 0.
          iLAMc2 = 0.;   iLAMc5= 0.
       endif

    ! Rain:
       if (QRM(i,k)>epsQ) then
          if (.not. dblMom_r) NRM(i,k)= GR50*sqrt(sqrt(GR31*iGR34*DE(i,k)*QRM(i,k)*icmr))
          iQRM   = 1./QRM(i,k)
          iNRM   = 1./NRM(i,k)
          Dr     = Dm_x(DE(i,k),QRM(i,k),iNRM,icmr,thrd)
          iLAMr  = max( iLAMmin1, iLAMDA_x(DE(i,k),QRM(i,k),iNRM,icexr9,thrd) )
          tmp1   = 1./iLAMr
          iLAMr2 = iLAMr *iLAMr
          iLAMr3 = iLAMr2*iLAMr
          iLAMr4 = iLAMr2*iLAMr2
          iLAMr5 = iLAMr3*iLAMr2
          if (Dr>40.e-6) then
             vr0 = gamfact_r(i,k)*ckQr1*iLAMr**bfr
          else
             vr0 = 0.
          endif
       else
          iLAMr  = 0.;  Dr    = 0.;  vr0   = 0.
          iLAMr2 = 0.;  iLAMr3= 0.;  iLAMr4= 0.;  iLAMr5 = 0.
       endif

    ! Ice:
       if (QIM(i,k)>epsQ) then
          if (.not. dblMom_i) NYM(i,k)= N_Cooper(TRPL,TM(i,k))

          iQIM   = 1./QIM(i,k)
          iNYM   = 1./NYM(i,k)
          iLAMi  = max( iLAMmin2, iLAMDA_x(DE(i,k),QIM(i,k),iNYM,icexi9,thrd) )
          iLAMi2 = iLAMi *iLAMi
          iLAMi3 = iLAMi2*iLAMi
          iLAMi4 = iLAMi2*iLAMi2
          iLAMi5 = iLAMi3*iLAMi2
          iLAMiB0= iLAMi**(bfi)
          iLAMiB1= iLAMi**(bfi+1.)
          iLAMiB2= iLAMi**(bfi+2.)
          vi0    = gamfact(i,k)*ckQi1*iLAMiB0
          Di     = Dm_x(DE(i,k),QIM(i,k),iNYM,icmi,thrd)
       else
          iLAMi  = 0.;  vi0    = 0.;  Di     = 0.
          iLAMi2 = 0.;  iLAMi3 = 0.;  iLAMi4 = 0.;  iLAMi5= 0.
          iLAMiB0= 0.;  iLAMiB1= 0.;  iLAMiB2= 0.
       endif

    ! Snow:
       if (QNM(i,k)>epsQ) then
          if (.not.dblMom_s) then
             No_s_SM = Nos_Thompson(TRPL,TM(i,k))
             NNM(i,k)= (No_s*GS31)**(dms*icexs2)*(GS31*iGS40*icms*DE(i,k)*QNM(i,k))**    &
                       ((1.+alpha_s)*icexs2)
          endif
          iQNM   = 1./QNM(i,k)
          iNNM   = 1./NNM(i,k)
          iLAMs  = max( iLAMmin2, iLAMDA_x(DE(i,k),QNM(i,k),iNNM,iGS20,idms) )
          iLAMs_D3= max(iLAMmin2, iLAMDA_x(DE(i,k),QNM(i,k),iNNM,iGS20_D3,thrd) )
          iLAMs2 = iLAMs*iLAMs
          iLAMsB0= iLAMs**(bfs)
          iLAMsB1= iLAMs**(bfs+1.)
          iLAMsB2= iLAMs**(bfs+2.)
          vs0    = gamfact(i,k)*ckQs1*iLAMsB0
          Ds     = min(DsMax, Dm_x(DE(i,k),QNM(i,k),iNNM,icms,idms))
          if (snowSpherical) then
             des = desFix
          else
             des = des_OF_Ds(Ds,desMax,eds,fds)
          endif
         !!-- generalized equations (any alpha_s):
         !    No_s  = (NNM(i,k))*iGS31/iLAMs**(1.+alpha_s)
         !    VENTs = Avx*GS32*iLAMs**(2.+alpha_s)+Bvx*ScTHRD*sqrt(gam*afs*iMUkin)*      &
         !!--         GS09*iLAMs**(2.5+0.5*bfs+alpha_s)
         !The following equations for No_s and VENTs is based on m(D)=(pi/6)*100.*D**3 for snow.
         !  Strict application of m(D)=c*D**2 would require re-derivation using implied
         !  definition of D as the MAXIMUM DIMENSION of an ellipsoid, rather than a sphere.
         !  For simplicity, the m-D^3 relation is applied -- used for VDvs and MLsr only.
         if (dblMom_s) then
           !No_s= NNM(i,k)*iGS31/iLAMs     !optimized for alpha_s=0
            No_s= NNM(i,k)*iGS31/iLAMs_D3  !based on m-D^3 (consistent with VENTs, below)
         else
            No_s= No_s_SM
         endif
         VENTs= Avx*GS32*iLAMs_D3**2. + Bvx*ScTHRD*sqrt(gamfact(i,k)*afs*iMUkin)*GS09*   &
                iLAMs_D3**cexs1
       else
          iLAMs  = 0.;  vs0    = 0.;  Ds     = 0.;  iLAMs2= 0.
          iLAMsB0= 0.;  iLAMsB1= 0.;  iLAMsB1= 0.
          des    = desFix !used for 3-component freezing if QNM=0 (even for snowSpherical=.F.)
       endif
       ides  = 1./des


    ! Graupel:
       if (QGM(i,k)>epsQ) then
          if (.not.dblMom_g) NGM(i,k)= GG50*sqrt(sqrt(GG31*GG34*DE(i,k)*QGM(i,k)*icmg))
          iQGM   = 1./QGM(i,k)
          iNGM   = 1./NGM(i,k)
          iLAMg  = max( iLAMmin1, iLAMDA_x(DE(i,k),QGM(i,k),iNGM,iGG99,thrd) )
          iLAMg2 = iLAMg *iLAMg
          iLAMgB0= iLAMg**(bfg)
          iLAMgB1= iLAMg**(bfg+1.)
          iLAMgB2= iLAMg**(bfg+2.)
          if (dblMom_g) then
            !No_g = (NGM(i,k))*iGG31/iLAMg**(1.+alpha_g)
             No_g= NGM(i,k)*iGG31/iLAMg     !optimized for alpha_g=0
          else
             No_g= No_g_SM
          endif
          vg0    = gamfact(i,k)*ckQg1*iLAMgB0
          Dg     = Dm_x(DE(i,k),QGM(i,k),iNGM,icmg,thrd)
       else
          iLAMg  = 0.;  vg0    = 0.;  Dg     = 0.;  No_g   = 0.
          iLAMg2 = 0.;  iLAMgB0= 0.;  iLAMgB1= 0.;  iLAMgB1= 0.
       endif

    ! Hail:
       if (QHM(i,k)>epsQ) then
          if (.not.dblMom_h) NHM(i,k)= GH50*sqrt(sqrt(GH31*iGH34*DE(i,k)*QHM(i,k)*icmh))
          iQHM   = 1./QHM(i,k)
          iNHM   = 1./NHM(i,k)
          iLAMh  = max( iLAMmin1, iLAMDA_x(DE(i,k),QHM(i,k),iNHM,iGH99,thrd) )
          iLAMh2 = iLAMh*iLAMh
          iLAMhB0= iLAMh**(bfh)
          iLAMhB1= iLAMh**(bfh+1.)
          iLAMhB2= iLAMh**(bfh+2.)
          if (dblMom_h) then
               No_h= NHM(i,k)*iGH31/iLAMh**(1.+alpha_h)
          else
               No_h= No_h_SM
          endif
          vh0    = gamfact(i,k)*ckQh1*iLAMhB0
          Dh     = Dm_x(DE(i,k),QHM(i,k),iNHM,icmh,thrd)
       else
          iLAMh  = 0.;  vh0    = 0.;  Dh     = 0.;  No_h= 0.
          iLAMhB0= 0.;  iLAMhB1= 0.;  iLAMhB1= 0.
       endif
!------

 !Calculating ice-phase source/sink terms:

 ! Initialize all source terms to zero:
       QNUvi=0.;  QVDvi=0.;  QVDvs=0.;  QVDvg=0.;  QVDvh=0.;   QCNig=0.
       QCLcs=0.;  QCLcg=0.;  QCLch=0.;  QFZci=0.;  QCLri=0.;   QMLsr=0.
       QCLrs=0.;  QCLrg=0.;  QMLgr=0.;  QCLrh=0.;  QMLhr=0.;   QFZrh=0.
       QMLir=0.;  QCLsr=0.;  QCLsh=0.;  QCLgr=0.;  QCNgh=0.
       QCNis=0.;  QCLir=0.;  QCLis=0.;  QCLih=0.
       QIMsi=0.;  QIMgi=0.;  QCNsg=0.;  QHwet=0.

       NCLcs= 0.; NCLcg=0.;  NCLch=0.;  NFZci=0.;  NMLhr=0.;   NCNig=0.
       NCLri= 0.; NCLrs=0.;  NCLrg=0.;  NCLrh=0.;  NMLsr=0.;   NMLgr=0.
       NMLir= 0.; NSHhr=0.;  NNUvi=0.;  NVDvi=0.;  NVDvh=0.;   QCLig=0.
       NCLir= 0.; NCLis=0.;  NCLig=0.;  NCLih=0.;  NIMsi=0.;   NIMgi=0.
       NiCNis=0.; NsCNis=0.; NVDvs=0.;  NCNsg=0.;  NCLgr=0.;   NCLsrh=0.
       NCLss= 0.; NCLsr=0.;  NCLsh=0.;  NCLsrs=0.; NCLgrg=0.;  NCNgh=0.
       NVDvg= 0.; NCLirg=0.; NCLsrg=0.; NCLgrh=0.; NrFZrh=0.;  NhFZrh=0.
       NCLirh=0.

       Dirg=0.; Dirh=0.; Dsrs= 0.; Dsrg= 0.; Dsrh= 0.; Dgrg=0.; Dgrh=0.

   !-------------------------------------------------------------------------------------------!

           ! COLLECTION by snow, graupel, hail:
           !  (i.e. wet or dry ice-categories [=> excludes ice crystals])

           ! Collection by SNOW:
       if (QNM(i,k)>epsQ) then
          ! cloud:
          if (QCM(i,k)>epsQ) then

            !Approximation of Ecs based on Pruppacher & Klett (1997) Fig. 14-11
             Ecs= min(Dc,30.e-6)*3.333e+4*sqrt(min(Ds,1.e-3)*1.e+3)
             QCLcs= dt*gam*afs*cmr*Ecs*PIov4*iDE(i,k)*(NCM(i,k)*NNM(i,k))*iGC5*iGS31*    &
                    (GC13*GS13*iLAMc3*iLAMsB2+2.*GC14*GS12*iLAMc4*iLAMsB1+GC15*GS11*     &
                    iLAMc5*iLAMsB0)

             NCLcs= dt*gam*afs*PIov4*Ecs*(NCM(i,k)*NNM(i,k))*iGC5*iGS31*(GC5*GS13*       &
                    iLAMsB2+2.*GC11*GS12*iLAMc*iLAMsB1+GC12*GS11*iLAMc2*iLAMsB0)

            !continuous collection: (alternative; gives values ~0.95 of SCE [above])
            !QCLcs= dt*gam*Ecs*PIov4*afs*QCM(i,k)*NNM(i,k)*iLAMs**(2.+bfs)*GS13*iGS31
            !NCLcs= QCLcs*NCM(i,k)/QCM(i,k)

            !Correction factor for non-spherical snow [D = maximum dimension] which
            !changes projected area:   [assumption: A=0.50*D**2 (vs. A=(PI/4)*D**2)]
            ! note: Strictly speaking, this correction should only be applied to
            !       continuous growth approximation for cloud.  [factor = 0.50/(pi/4)]
             if (.not. snowSpherical) then
                tmp1 = 0.6366      !factor = 0.50/(pi/4)
                QCLcs= tmp1*QCLcs
                NCLcs= tmp1*NCLcs
             endif

             QCLcs= min(QCLcs, QCM(i,k))
             NCLcs= min(NCLcs, NCM(i,k))
          else
             QCLcs= 0.;   NCLcs= 0.
          endif

          ! ice:
          if (QIM(i,k)>epsQ) then
             tmp1= vs0-vi0
             tmp3= sqrt(tmp1*tmp1+0.04*vs0*vi0)

             QCLis= dt*cmi*iDE(i,k)*PI*6.*Eis*(NYM(i,k)*NNM(i,k))*tmp3*iGI31*iGS31*(0.5* &
                    iLAMs2*iLAMi3+2.*iLAMs*iLAMi4+5.*iLAMi5)

             NCLis= dt*PIov4*Eis*(NYM(i,k)*NNM(i,k))*GI31*GS31*tmp3*(GI33*GS31*iLAMi2+   &
                    2.*GI32*GS32*iLAMi*iLAMs+GI31*GS33*iLAMs2)

             QCLis= min(QCLis, (QIM(i,k)))
             NCLis= min(QCLis*(NYM(i,k)*iQIM), NCLis)
          else
             QCLis= 0.;   NCLis= 0.
          endif

          if (dblMom_s) then
             !snow: (i.e. self-collection [aggregation])
             NCLss= dt*0.93952*Ess*(DE(i,k)*(QNM(i,k)))**((2.+bfs)*thrd)*(NNM(i,k))**    &
                    ((4.-bfs)*thrd)
               !Note: 0.91226 = I(bfs)*afs*PI^((1-bfs)/3)*des^((-2-bfs)/3); I(bfs=0.41)=1138
               !      0.93952 = I(bfs)*afs*PI^((1-bfs)/3)*des^((-2-bfs)/3); I(bfs=0.42)=1172
               !      [interpolated from 3rd-order polynomial approx. of values given in RRB98;
               !       see eqn(A.35)]
             NCLss= min(NCLss, 0.5*(NNM(i,k)))
          endif

       else
          QCLcs= 0.;   NCLcs= 0.;   QCLis= 0.;   NCLis= 0.;  NCLss= 0.
       endif

       ! Collection by GRAUPEL:
       if (QGM(i,k)>epsQ) then

          ! cloud:
          if (QCM(i,k)>epsQ) then

             QCLcg= dt*gam*afg*cmr*Ecg*PIov4*iDE(i,k)*(NCM(i,k)*NGM(i,k))*iGC5*iGG31*    &
                    (GC13*GG13*iLAMc3*iLAMgB2+ 2.*GC14*GG12*iLAMc4*iLAMgB1+GC15*GG11*    &
                    iLAMc5*iLAMgB0)

             NCLcg= dt*gam*afg*PIov4*Ecg*(NCM(i,k)*NGM(i,k))*iGC5*iGG31*(GC5*GG13*       &
                    iLAMgB2+2.*GC11*GG12*iLAMc*iLAMgB1+GC12*GG11*iLAMc2*iLAMgB0)

             QCLcg= min(QCLcg, (QCM(i,k)))
             NCLcg= min(NCLcg, (NCM(i,k)))
          else
             QCLcg= 0.;   NCLcg= 0.
          endif

          ! ice:
          if (QIM(i,k)>epsQ) then
             tmp1= vg0-vi0
             tmp3= sqrt(tmp1*tmp1+0.04*vg0*vi0)

             QCLig= dt*cmi*iDE(i,k)*PI*6.*Eig*(NYM(i,k)*NGM(i,k))*tmp3*iGI31*iGG31*(0.5* &
                    iLAMg2*iLAMi3+2.*iLAMg*iLAMi4+5.*iLAMi5)
             NCLig= dt*PIov4*Eig*(NYM(i,k)*NGM(i,k))*GI31*GG31*tmp3*(GI33*GG31*iLAMi2+   &
                    2.*GI32*GG32*iLAMi*iLAMg+GI31*GG33*iLAMg2)

             QCLig= min(QCLig, (QIM(i,k)))
             NCLig= min(QCLig*(NYM(i,k)*iQIM), NCLig)
          else
             QCLig= 0.;   NCLig= 0.
          endif

       else
          QCLcg= 0.;   QCLrg= 0.;   QCLig= 0.
          NCLcg= 0.;   NCLrg= 0.;   NCLig= 0.
       endif

       ! Collection by HAIL:
       if (QHM(i,k)>epsQ) then

         ! cloud:
          if (QCM(i,k)>epsQ) then
             Ech  = exp(-8.68e-7*Dc**(-1.6)*Dh)    !Ziegler (1985) A24

             QCLch= dt*gam*afh*cmr*Ech*PIov4*iDE(i,k)*(NCM(i,k)*NHM(i,k))*iGC5*iGH31*    &
                    (GC13*GH13*iLAMc3*iLAMhB2+2.*GC14*GH12*iLAMc4*iLAMhB1+GC15*GH11*     &
                    iLAMc5*iLAMhB0)

             NCLch= dt*gam*afh*PIov4*Ech*(NCM(i,k)*NHM(i,k))*iGC5*iGH31*(GC5*GH13*       &
                    iLAMhB2+2.*GC11*GH12*iLAMc*iLAMhB1+GC12*GH11*iLAMc2*iLAMhB0)

             QCLch= min(QCLch, QCM(i,k))
             NCLch= min(NCLch, NCM(i,k))
          else
             QCLch= 0.;   NCLch= 0.
          endif

          ! rain:
          if (QRM(i,k)>epsQ) then
             tmp1= vh0-vr0
             tmp3= sqrt(tmp1*tmp1+0.04*vh0*vr0)
             QCLrh= dt*cmr*Erh*PIov4*iDE(i,k)*(NHM(i,k)*NRM(i,k))*iGR31*iGH31*tmp3*      &
                    (GR36*GH31*iLAMr5+2.*GR35*GH32*iLAMr4*iLAMh+GR34*GH33*iLAMr3*iLAMh2)

             NCLrh= dt*PIov4*Erh*(NHM(i,k)*NRM(i,k))*iGR31*iGH31*tmp3*(GR33*GH31*        &
                    iLAMr2+2.*GR32*GH32*iLAMr*iLAMh+GR31*GH33*iLAMh2)

             QCLrh= min(QCLrh, QRM(i,k))
             NCLrh= min(NCLrh, QCLrh*(NRM(i,k)*iQRM))
          else
             QCLrh= 0.;   NCLrh= 0.
          endif

          ! ice:
          if (QIM(i,k)>epsQ) then
             tmp1 = vh0-vi0
             tmp3 = sqrt(tmp1*tmp1+0.04*vh0*vi0)

             QCLih= dt*cmi*iDE(i,k)*PI*6.*Eih*(NYM(i,k)*NHM(i,k))*tmp3*iGI31*iGH31*(0.5* &
                    iLAMh2*iLAMi3+2.*iLAMh*iLAMi4+5.*iLAMi5)

             NCLih= dt*PIov4*Eih*(NYM(i,k)*NHM(i,k))*GI31*GH31*tmp3*(GI33*GH31*iLAMi2+   &
                    2.*GI32*GH32*iLAMi*iLAMh+GI31*GH33*iLAMh2)

             QCLih= min(QCLih, QIM(i,k))
             NCLih= min(QCLih*(NYM(i,k)*iQIM), NCLih)
          else
             QCLih= 0.;   NCLih= 0.
          endif

          ! snow:
          if (QNM(i,k)>epsQ) then
             tmp1 = vh0-vs0
             tmp3 = sqrt(tmp1*tmp1+0.04*vh0*vs0)
             tmp4 = iLAMs2*iLAMs2

             if (snowSpherical) then
               !hardcoded for dms=3:
                QCLsh= dt*cms*iDE(i,k)*PI*6.*Esh*(NNM(i,k)*NHM(i,k))*tmp3*iGS31*iGH31*  &
                       (0.5*iLAMh2*iLAMs2*iLAMs+2.*iLAMh*tmp4+5.*tmp4*iLAMs)
             else
               !hardcoded for dms=2:
                QCLsh= dt*cms*iDE(i,k)*PI*0.25*Esh*tmp3*NNM(i,k)*NHM(i,k)*iGS31*iGH31*  &
                       (GH33*GS33*iLAMh**2.*iLAMs**2. + 2.*GH32*GS34*iLAMh*iLAMs**3. +  &
                        GH31*GS35*iLAMs**4.)
             endif

             NCLsh= dt*PIov4*Esh*(NNM(i,k)*NHM(i,k))*GS31*GH31*tmp3*(GS33*GH31*iLAMs2+  &
                    2.*GS32*GH32*iLAMs*iLAMh+GS31*GH33*iLAMh2)

             QCLsh= min(QCLsh, (QNM(i,k)))
             NCLsh= min((NNM(i,k)*iQNM)*QCLsh, NCLsh, (NNM(i,k)))
          else
             QCLsh= 0.;   NCLsh= 0.
          endif

         !wet growth:
          VENTh= Avx*GH32*iLAMh**(2.+alpha_h) + Bvx*ScTHRD*sqrt(gam*afh*iMUkin)*GH09*    &
                 iLAMh**(2.5+0.5*bfh+alpha_h)
          QHwet= max(0., dt*PI2*(DE(i,k)*CHLC*Cdiff*DELqvs-Ka*Tc)*No_h*iDE(i,k)/(CHLF+   &
                 CPW*Tc)*VENTh+(QCLih*iEih+QCLsh*iEsh)*(1.-CPI_*Tc/(CHLF+CPW*Tc)) )

       else
          QCLch= 0.;   QCLrh= 0.;   QCLih= 0.;   QCLsh= 0.;   QHwet= 0.
          NCLch= 0.;   NCLrh= 0.;   NCLsh= 0.;   NCLih= 0.
       endif

       IF (TM(i,k)>TRPL .and. warmphase_ON) THEN
          !**********!
          !  T > To  !
          !**********!

          ! MELTING of frozen particles:
          !  ICE:
          QMLir   = QIM(i,k)  !all pristine ice melts in one time step
          QIM(i,k)= 0.
          NMLir   = NYM(i,k)

          !  SNOW:
          if (QNM(i,k)>epsQ) then
             QMLsr= dt*(PI2*iDE(i,k)*iCHLF*No_s*VENTs*(Ka*Tc-CHLC*Cdiff*DELqvs) + CPW*   &
                    iCHLF*Tc*(QCLcs+QCLrs)*idt)
             QMLsr= min(max(QMLsr,0.), QNM(i,k))
             NMLsr= NNM(i,k)*iQNM*QMLsr
          else
             QMLsr= 0.;   NMLsr= 0.
          endif

          !  GRAUPEL:
          if (QGM(i,k)>epsQ) then
             VENTg= Avx*GG32*iLAMg*iLAMg+Bvx*ScTHRD*sqrt(gam*afg*iMUkin)*GG09*iLAMg**    &
                    (2.5+0.5*bfg+alpha_g)
             QMLgr= dt*(PI2*iDE(i,k)*iCHLF*No_g*VENTg*(Ka*Tc-CHLC*Cdiff*DELqvs) + CPW*   &
                    iCHLF*Tc*(QCLcg+QCLrg)*idt)
             QMLgr= min(max(QMLgr,0.), QGM(i,k))
             NMLgr= NGM(i,k)*iQGM*QMLgr
          else
             QMLgr= 0.;   NMLgr= 0.
          endif

          !  HAIL:
          if (QHM(i,k)>epsQ.and.Tc>5.) then
             VENTh= Avx*GH32*iLAMh**(2.+alpha_h) + Bvx*ScTHRD*sqrt(gam*afh*iMUkin)*GH09* &
                    iLAMh**(2.5+0.5*bfh+alpha_h)
             QMLhr= dt*(PI2*iDE(i,k)*iCHLF*No_h*VENTh*(Ka*Tc-CHLC*Cdiff*DELqvs) + CPW/   &
                    CHLF*Tc*(QCLch+QCLrh)*idt)
             QMLhr= min(max(QMLhr,0.), QHM(i,k))
             NMLhr= NHM(i,k)*iQHM*QMLhr
             if(QCLrh>0.) NMLhr= NMLhr*0.1   !Prevents problems when hail is ML & CL
          else
             QMLhr= 0.;   NMLhr= 0.
          endif

         ! Cold (sub-zero) source/sink terms:
          QNUvi= 0.;   QFZci= 0.;   QVDvi= 0.;   QVDvs= 0.;   QVDvg= 0.
          QCLis= 0.;   QCNis1=0.;   QCNis2=0.;   QCNig= 0.
          QCNgh= 0.;   QIMsi= 0.;   QIMgi= 0.;   QCLir= 0.;   QCLri= 0.
          QCLrs= 0.;   QCLgr= 0.;   QCLrg= 0.;   QCNis= 0.;   QVDvh= 0.
          QCNsg= 0.;   QCLsr= 0.

          NNUvi= 0.;   NFZci= 0.;   NCLgr= 0.;   NCLrg= 0.;   NCNig= 0.
          NCLis= 0.;   NVDvi= 0.;   NVDvs= 0.;   NVDvg= 0.;   NVDvh= 0.
          NCNsg= 0.;   NCNgh= 0.;   NiCNis=0.;   NsCNis=0.;   NCLrs= 0.
          NIMsi= 0.;   NIMgi= 0.;   NCLir= 0.;   NCLri= 0.;   NCLsr= 0.

       ELSE
          !----------!
          !  T < To  !
          !----------!
          if (QSI(i,k)>1.e-20) then
             Si = QM(i,k)/QSI(i,k)
          else
             Si = 0.
          endif
          tmp2  = TM(i,k)*TM(i,k)
          iABi  = 1./( CHLS*CHLS/(Ka*RGASV*tmp2) + 1./(DE(i,k)*(QSI(i,k))*Cdiff) )

          ! Warm-air-only source/sink terms:
          QMLir= 0.;   QMLsr= 0.;   QMLgr= 0.;   QMLhr= 0.
          NMLir= 0.;   NMLsr= 0.;   NMLgr= 0.;   NMLhr= 0.

          !Probabilistic freezing (Bigg) of rain:
          if (Tc<Tc_FZrh .and. QRM(i,k)>epsQ .and. hail_ON) then
             !note: - (Tc<-10.C) condition is based on Pruppacher-Klett (1997) Fig. 9-41
             !      - Small raindrops will freeze to hail. However, if after all S/S terms
             !        are added Dh<Dh_min, then hail will be converted to graupel. Thus,
             !        probabilistic freezing of small rain is effectively a source of graupel.
             NrFZrh= -dt*Bbigg*(exp(Abigg*Tc)-1.)*DE(i,k)*QRM(i,k)*idew
             Rz= 1.  !N and Z (and Q) are conserved for FZrh with triple-moment
           ! The Rz factor serves to conserve reflectivity when a rain distribution
           !  converts to an distribution with a different shape parameter, alpha.
           !  (e.g. when rain freezes to hail)  The factor Rz non-conserves N while
           !  acting to conserve Z for double-moment.  See Ferrier, 1994 App. D)
           ! Rz= (gamma(7.d0+alpha_h)*GH31*GR34*GR34)/(GR36(i,k)*GR31*   &
           !      gamma(4.d0+alpha_h)*gamma(4.d0+alpha_h))
             NhFZrh= Rz*NrFZrh
             QFZrh = NrFZrh*(QRM(i,k)*iNRM)
          else
             QFZrh= 0.;   NrFZrh= 0.;  NhFZrh= 0.
          endif

          !--------!
          !  ICE:  !
          !--------!
          ! Homogeneous freezing of cloud to ice:
          if (dblMom_c) then
             if (QCM(i,k)>epsQ) then
                tmp2  = Tc*Tc; tmp3= tmp2*Tc; tmp4= tmp2*tmp2
                JJ    = (10.**max(-20.,(-606.3952-52.6611*Tc-1.7439*tmp2-0.0265*tmp3-    &
                         1.536e-4*tmp4)))
                tmp1  = 1.e6*(DE(i,k)*(QCM(i,k)*iNCM)*icmr) !i.e. Dc[cm]**3
                FRAC  = 1.-exp(-JJ*PIov6*tmp1*dt)
                if (Tc>-30.) FRAC= 0.
                if (Tc<-50.) FRAC= 1.
                QFZci= FRAC*QCM(i,k)
                NFZci= FRAC*NCM(i,k)
             else
                QFZci= 0.;   NFZci= 0.
             endif
          else
             !Homogeneous freezing of cloud to ice:  (simplified)
             if (QCM(i,k)>epsQ .and. Tc<-35.) then
                FRAC= 1.  !if T<-35
                QFZci= FRAC*QCM(i,k)
                NFZci= FRAC*N_c_SM
             else
                QFZci= 0.;   NFZci= 0.
             endif
          endif

          if (dblMom_i) then
            !Primary ice nucleation:
            NNUvi= 0.;   QNUvi= 0.
            if (primIceNucl==1) then

               NuDEPSOR= 0.;   NuCONT= 0.
               if (QSI(i,k)>1.e-20) then
                  Simax = min(Si, satw_peak*QSW(i,k)/QSI(i,k))
               else
                  Simax = 0.
               endif
               tmp1    = T(i,k)-7.66
               NNUmax  = max(0., DE(i,k)/mio*(Q(i,k)-QSS(i,k))/(1.+ck6*(QSS(i,k)/(tmp1*    &
                        tmp1))))
               !Deposition/sorption nucleation:
               if (Tc<-5. .and. Si>1.) then
                  NuDEPSOR= max(0., 1.e3*exp(12.96*(Simax-1.)-0.639)-(NYM(i,k))) !Meyers(1992)
               endif
               !Contact nucleation:
               if (QCM(i,k)>epsQ .and. Tc<-2. .and. WZ(i,k)>0.001) then
                  GG     =  1.*idew/(RGASV*(TM(i,k))/((QSW(i,k)*HPS(i)*S(i,k))/EPS1)/      &
                              Cdiff+CHLC/Ka/(TM(i,k))*(CHLC/RGASV/(TM(i,k))-1.))  !CP00a
                  Swmax  =  SxFNC(WZ(i,k),Tc,HPS(i)*S(i,k),QSW(i,k),QSI(i,k),CCNtype,1)
                  if (QSW(i,k)>1.e-20) then
                     ssat=  min((QM(i,k)/QSW(i,k)), Swmax) -1.
                  else
                     ssat= 0.
                  endif
                  Tcc    =  Tc + GG*ssat*CHLC/Kdiff                            !C86_eqn64
                  Na     =  exp(4.11-0.262*Tcc)                                !W95_eqn60/M92_2.6
                  Kn     =  LAMa0*(TM(i,k))*p0/(T0*(HPS(i)*S(i,k))*Ra)         !W95_eqn59
                  PSIa   =  -kBoltz*Tcc/(6.*pi*Ra*MUdyn)*(1.+Kn)               !W95_eqn58
                  ft     =  0.4*(1.+1.45*Kn+0.4*Kn*exp(-1./Kn))*(Ka+2.5*Kn*KAPa)/          &
                           (1.+3.*Kn)/(2.*Ka+5.*KAPa*Kn+KAPa)                  !W95_eqn57
                  Dc     =  (DE(i,k)*(QCM(i,k)*iNCM)*icmr)**thrd
                  F1     =  PI2*Dc*Na*(NCM(i,k))                               !W95_eqn55
                  F2     =  Ka/(HPS(i)*S(i,k))*(Tc-Tcc)                        !W95_eqn56
                  NuCONTA= -F1*F2*RGASV*(TM(i,k))/CHLC*iDE(i,k)                !diffusiophoresis
                  NuCONTB=  F1*F2*ft*iDE(i,k)                                  !thermeophoresis
                  NuCONTC=  F1*PSIa                                            !Brownian diffusion
                  NuCONT =  max(0.,(NuCONTA+NuCONTB+NuCONTC)*dt)
               endif
               !Total primary ice nucleation:
               if (icephase_ON) then
                  NNUvi= min(NNUmax, NuDEPSOR + NuCONT )
                  QNUvi= mio*iDE(i,k)*NNUvi
                  QNUvi= min(QNUvi,(Q(i,k)))
               endif

            elseif (primIceNucl==2) then
               if (Tc<-5. .and. Si>1.08) then !following Thompson etal (2006)
                  NNUvi= max(N_Cooper(TRPL,T(i,k))-NYM(i,k),0.)
                  QNUvi= min(mio*iDE(i,k)*NNUvi, Q(i,k))
               endif
           !elseif (primIceNucl==3) then
           !! (for alternative [future] ice nucleation parameterizations)
           !   NNUvi=...
           !   QNUvi=...
            endif !if (primIceNucl==1)

          else !dblMom_i
          !Ice initiation (single-moment):
             if (QIM(i,k)<=epsQ .and. Tc<-5. .and. Si>1.08) then !following Thompson etal (2006)
                NNUvi = N_Cooper(TRPL,T(i,k))
                QNUvi= mio*iDE(i,k)*NNUvi
                QNUvi= min(QNUvi,Q(i,k))
             endif
          endif !dblMom_i


          IF (QIM(i,k)>epsQ) THEN

             !Deposition/sublimation:
!            No_i  = NYM(i,k)*iGI31/iLAMi**(1.+alpha_i)
!            VENTi= Avx*GI32*iLAMi**(2.+alpha_i)+Bvx*ScTHRD*sqrt(gam*afi*iMUkin)*GI6*    &
!                     iLAMi**(2.5+0.5*bfi+alpha_i)
             No_i  = NYM(i,k)*iGI31/iLAMi    !optimized for alpha_i=0
             VENTi= Avx*GI32*iLAMi*iLAMi+Bvx*ScTHRD*sqrt(gam*afi*iMUkin)*GI6*iLAMi**     &
                    (2.5+0.5*bfi+alpha_i)
            !Note: ice crystal capacitance is implicitly C = 0.5*D*capFact_i
             QVDvi= dt*capFact_i*iABi*(PI2*(Si-1.)*No_i*VENTi)

             ! Prevent overdepletion of vapor:
             tmp1  = T(i,k)-7.66
             VDmax = (Q(i,k)-QSS(i,k))/(1.+ck6*(QSS(i,k))/(tmp1*tmp1))
             if(Si>=1.) then
                QVDvi= min(max(QVDvi,0.),VDmax)
             else
                if (VDmax<0.) QVDvi= max(QVDvi,VDmax)
               !IF prevents subl.(QVDvi<0 at t) changing to dep.(VDmax>0 at t*)  2005-06-28
             endif
             if (.not. iceDep_ON) QVDvi= 0. !suppresses depositional growth
             NVDvi= min(0., (NYM(i,k)*iQIM)*QVDvi) !dNi/dt=0 for deposition

             ! Conversion to snow:
             !   +depostion of ice:
             mi= DE(i,k)*(QIM(i,k)*iNYM)
             if (mi<=0.5*mso.and.abs(0.5*mso-mi)>1.e-20) then
                QCNis1= (mi/(mso-mi))*QVDvi
             else
                QCNis1= QVDvi + (1.-0.5*mso/mi)*QIM(i,k)
             endif
             QCNis1= max(0., QCNis1)
             !   +aggregation of ice:
             if(Di<0.5*Dso) then
                Ki    = PIov6*Di*Di*vi0*Eii*Xdisp
                tmp1  = log(Di/Dso)
                tmp2  = tmp1*tmp1*tmp1
                QCNis2= -dt*0.5*(QIM(i,k)*NYM(i,k))*Ki/tmp2
             else
                Ki= 0.;   QCNis2= 0.
             endif
             !   +total conversion rate:
             QCNis = QCNis1 + QCNis2
             NsCNis= DE(i,k)*imso*QCNis                               !source for snow (Ns)
             NiCNis= (DE(i,k)*imso*QCNis1 + 0.5*Ki*NYM(i,k)*NYM(i,k)) !sink for ice (Ni)
             NiCNis= min(NiCNis, NYM(i,k)*0.1) !Prevents overdepl. of NY when final QI>0

             if (.not.(snow_ON)) then
                QCNis= 0.; NiCNis= 0.; NsCNis= 0.  !Suppress SNOW initiation
             endif

             ! 3-component freezing (collisions with rain):
             if (QRM(i,k)>epsQ .and. QIM(i,k)>epsQ .and. Tc<Tc_3comp(1)) then
                tmp1 = vr0-vi0
                tmp3 = sqrt(tmp1*tmp1+0.04*vr0*vi0)

                QCLir= dt*cmi*Eri*PIov4*iDE(i,k)*(NRM(i,k)*NYM(i,k))*iGI31*iGR31*tmp3*   &
                       (GI36*GR31*iLAMi5+2.*GI35*GR32*iLAMi4*iLAMr+GI34*GR33*iLAMi3*     &
                       iLAMr2)

                NCLri= dt*PIov4*Eri*(NRM(i,k)*NYM(i,k))*iGI31*iGR31*tmp3*(GI33*GR31*     &
                       iLAMi2+2.*GI32*GR32*iLAMi*iLAMr+GI31*GR33*iLAMr2)

                QCLri= dt*cmr*Eri*PIov4*iDE(i,k)*(NYM(i,k)*NRM(i,k))*iGR31*iGI31*tmp3*   &
                       (GR36*GI31 *iLAMr5+2.*GR35*GI32*iLAMr4*iLAMi+GR34*GI33*iLAMr3*    &
                       iLAMi2)

               !note: For explicit eqns, both NCLri and NCLir are mathematically identical)
                NCLir= min(QCLir*(NYM(i,k)*iQIM), NCLri)
                QCLri= min(QCLri, (QRM(i,k)))
                QCLir= min(QCLir, (QIM(i,k)))
                NCLri= min(NCLri, (NRM(i,k)))
                NCLir= min(NCLir, (NYM(i,k)))

                !Determine destination of 3-comp.freezing:
                tmp1= max(Di,Dr)
                dey= (dei*Di*Di*Di+dew*Dr*Dr*Dr)/(tmp1*tmp1*tmp1)
                if (dey>0.5*(deg+deh) .and. Dr>Dr_3cmpThrs .and. hail_ON) then
                   Dirg= 0.;  Dirh= 1.
                else
                   Dirg= 1.;  Dirh= 0.
                endif
                if (.not. grpl_ON) Dirg= 0.

             else
                QCLir= 0.;  NCLir= 0.;  QCLri= 0.
                NCLri= 0.;  Dirh = 0.;  Dirg= 0.
             endif

             !  Rime-splintering (ice multiplication):
             ff= 0.
             if(Tc>=-8..and.Tc<=-5.) ff= 3.5e8*(Tc +8.)*thrd
             if(Tc> -5..and.Tc< -3.) ff= 3.5e8*(-3.-Tc)*0.5
             NIMsi= DE(i,k)*ff*QCLcs
             NIMgi= DE(i,k)*ff*QCLcg
             QIMsi= mio*iDE(i,k)*NIMsi
             QIMgi= mio*iDE(i,k)*NIMgi

          ELSE

             QVDvi= 0.;  QCNis= 0.;   QCNig= 0.
             QIMsi= 0.;  QIMgi= 0.;   QCLri= 0.;   QCLir= 0.
             NVDvi= 0.;  NCLir= 0.;   NIMsi= 0.;   NCNig= 0.
             NiCNis=0.;  NsCNis=0.;   NIMgi= 0.;   NCLri= 0.

          ENDIF
          !---------!
          !  SNOW:  !
          !---------!
          IF (QNM(i,k)>epsQ) THEN

            !Deposition/sublimation:
             !note: - snow crystal capacitance is implicitly C = 0.5*D*capFact_s
             !      - No_s and VENTs are computed above
             QVDvs = dt*capFact_s*iABi*(PI2*(Si-1.)*No_s*VENTs - CHLS*CHLF/(Ka*RGASV*    &
                     TM(i,k)*TM(i,k))*QCLcs*idt)

             ! Prevent overdepletion of vapor:
             tmp1  = T(i,k)-7.66
             VDmax = (Q(i,k)-QSS(i,k))/(1.+ck6*(QSS(i,k))/(tmp1*tmp1))  !KY97_A.33
             if(Si>=1.) then
                QVDvs= min(max(QVDvs,0.),VDmax)
             else
                if (VDmax<0.) QVDvs= max(QVDvs,VDmax)
                !IF prevents subl.(QVDvs<0 at t) changing to dep.(VDmax>0 at t*)
             endif
             NVDvs= -min(0.,(NNM(i,k)*iQNM)*QVDvs)  !pos. quantity

             ! Conversion to graupel:
             if (QCLcs>CNsgThres*QVDvs .and. 0.99*deg>des) then
               !note: The (deg>des) condition equates to (Ds>330microns) for m(D)=0.069D^2
               !      relation for snow, which implies a variable bulk density.  The physical
               !      assumption in the QCNsg equation is that snow converts to graupel due
               !      to densification from riming.
               !      The 0.99 is to prevent overflow if des~deg
                QCNsg= (deg/(deg-des))*QCLcs
             else
                QCNsg= 0.
             endif
             if (.not. grpl_ON) QCNsg= 0.
             NCNsg= DE(i,k)*imgo*QCNsg
             NCNsg= min(NCNsg, (0.5*NNM(i,k)*iQNM)*QCNsg) !Prevents incorrect Ns-depletion

             ! 3-component freezing (collisions with rain):
             if (QRM(i,k)>epsQ .and. QNM(i,k)>epsQ .and. Tc<Tc_3comp(2)) then
                tmp1 = vs0-vr0
                tmp2 = sqrt(tmp1*tmp1+0.04*vs0*vr0)
                tmp6 = iLAMs2*iLAMs2*iLAMs

                QCLrs= dt*cmr*Ers*PIov4*iDE(i,k)*NNM(i,k)*NRM(i,k)*iGR31*iGS31*tmp2*     &
                       (GR36*GS31*iLAMr5+2.*GR35*GS32*iLAMr4*iLAMs+GR34*GS33*iLAMr3*     &
                       iLAMs2)

                NCLrs= dt*0.25e0*PI*Ers*(NNM(i,k)*NRM(i,k))*iGR31*iGS31*tmp2*(GR33*      &
                       GS31*iLAMr2+2.*GR32*GS32*iLAMr*iLAMs+GR31*GS33*iLAMs2)

                if (snowSpherical) then
                  !hardcoded for dms=3:
                   QCLsr= dt*cms*Ers*PIov4*iDE(i,k)*(NRM(i,k)*NNM(i,k))*iGS31*iGR31*     &
                          tmp2*(GS36*GR31*tmp6+2.*GS35*GR32*iLAMs2*iLAMs2*iLAMr+GS34*    &
                          GR33*iLAMs2*iLAMs*iLAMr2)
                else
                  !hardcoded for dms=2:
                   QCLsr= dt*cms*iDE(i,k)*PI*0.25*ERS*tmp2*NNM(i,k)*NRM(i,k)*iGS31*      &
                          iGR31*(GR33*GS33*iLAMr**2.*iLAMs**2. + 2.*GR32*GS34*iLAMr*     &
                          iLAMs**3. +GR31*GS35*iLAMs**4.)
                endif

               !note: For explicit eqns, NCLsr = NCLrs
                NCLsr= min(QCLsr*(NNM(i,k)*iQNM), NCLrs)
                QCLrs= min(QCLrs, QRM(i,k));  QCLsr= min(QCLsr, QNM(i,k))
                NCLrs= min(NCLrs, NRM(i,k));  NCLsr= min(NCLsr, NYM(i,k))

                ! Determine destination of 3-comp.freezing:
                Dsrs= 0.;   Dsrg= 0.;    Dsrh= 0.
                tmp1= max(Ds,Dr)
                tmp2= tmp1*tmp1*tmp1
                dey = (des*Ds*Ds*Ds + dew*Dr*Dr*Dr)/tmp2
                if (dey<=0.5*(des+deg)                        ) Dsrs= 1.  !snow
                if (dey >0.5*(des+deg) .and. dey<0.5*(deg+deh)) Dsrg= 1.  !graupel
                if (dey>=0.5*(deg+deh)) then
                   Dsrh= 1.                                               !hail
                   if (.not.hail_ON .or. Dr<Dr_3cmpThrs) then
                      Dsrg= 1.;   Dsrh= 0.                                !graupel
                   endif
                endif
                if (.not. grpl_ON) Dsrg=0.
             else
                QCLrs= 0.;   QCLsr= 0.;   NCLrs= 0.;   NCLsr= 0.
             endif

          ELSE

             QVDvs= 0.;  QCLcs= 0.;  QCNsg= 0.;  QCLsr= 0.;  QCLrs= 0.
             NVDvs= 0.;  NCLcs= 0.;  NCLsr= 0.;  NCLrs= 0.;  NCNsg= 0.

          ENDIF
          !------------!
          !  GRAUPEL:  !
          !------------!
          IF (QGM(i,k)>epsQ) THEN

           !Conversion to hail:    (Dho given by S-L limit)
             if (WZ(i,k)>w_CNgh .and. ((QCM(i,k)+QRM(i,k))>epsQ) .and. hail_ON) then
                tmp1 = 1.1e+4*DE(i,k)*(QCM(i,k)+QRM(i,k)) - 1.3e+3*DE(i,k)*QIM(i,k) + 1.
                tmp1 = max(1.,tmp1)  !to prevent div-by-zero
                Dho  = 0.01*(exp(min(20.,-Tc/tmp1))-1.)
                Dho  = 2.0*Dho !correction factor [to account error in Ziegler (1985), as per Young (1993)]
                Dho  = min(1., max(0.0001,Dho))  !smallest Dho=0.1mm; largest=1m
                ratio= Dg/Dho
                if (ratio>r_CNgh) then
                   QCNgh= (0.5*ratio)*(QCLcg+QCLrg+QCLig)
                   QCNgh= min(QCNgh,QGM(i,k)+QCLcg+QCLrg+QCLig)
                   NCNgh= DE(i,k)*QCNgh*icmh/(Dho*Dho*Dho)
                else
                   QCNgh= 0.
                   NCNgh= 0.
                endif
             endif

          !3-component freezing (collisions with rain):
             if (QRM(i,k)>epsQ .and. QGM(i,k)>epsQ .and. Tc<Tc_3comp(3)) then
                tmp1 = vg0-vr0
                tmp2 = sqrt(tmp1*tmp1 + 0.04*vg0*vr0)
                tmp8 = iLAMg2*iLAMg      ! iLAMg**3
                tmp9 = tmp8*iLAMg        ! iLAMg**4
                tmp10= tmp9*iLAMg        ! iLAMg**5

                QCLrg= dt*cmr*Erg*PIov4*iDE(i,k)*(NGM(i,k)*NRM(i,k))*iGR31*iGG31*tmp2*   &
                       (GR36*GG31*iLAMr5+2.*GR35*GG32*iLAMr4*iLAMg+GR34*GG33*iLAMr3*     &
                       iLAMg2)

                NCLrg= dt*PIov4*Erg*(NGM(i,k)*NRM(i,k))*iGR31*iGG31*tmp2*(GR33*GG31*     &
                       iLAMr2+2.*GR32*GG32*iLAMr*iLAMg+GR31*GG33*iLAMg2)

                QCLgr= dt*cmg*Erg*PIov4*iDE(i,k)*(NRM(i,k)*NGM(i,k))*iGG31*iGR31*tmp2*   &
                       (GG36*GR31*tmp10+2.*GG35*GR32*tmp9*iLAMr+GG34*GR33*tmp8*iLAMr2)

               !(note: For explicit eqns, NCLgr= NCLrg)
                NCLgr= min(NCLrg, QCLgr*(NGM(i,k)*iQGM))
                QCLrg= min(QCLrg, QRM(i,k));  QCLgr= min(QCLgr, QGM(i,k))
                NCLrg= min(NCLrg, NRM(i,k));  NCLgr= min(NCLgr, NGM(i,k))

               ! Determine destination of 3-comp.freezing:
                tmp1= max(Dg,Dr)
                tmp2= tmp1*tmp1*tmp1
                dey = (deg*Dg*Dg*Dg + dew*Dr*Dr*Dr)/tmp2
                if (dey>0.5*(deg+deh) .and. Dr>Dr_3cmpThrs .and. hail_ON) then
                   Dgrg= 0.;  Dgrh= 1.
                else
                   Dgrg= 1.;  Dgrh= 0.
                endif
             else
                QCLgr= 0.;  QCLrg= 0.;  NCLgr= 0.;  NCLrg= 0.
             endif

          ELSE

             QVDvg= 0.;  QCNgh= 0.;  QCLgr= 0.;  QCLrg= 0.
             NVDvg= 0.;  NCNgh= 0.;  NCLgr= 0.;  NCLrg= 0.

          ENDIF
          !---------!
          !  HAIL:  !
          !---------!
          IF (QHM(i,k)>epsQ) THEN

            !Wet growth:
             if (QHwet<(QCLch+QCLrh+QCLih+QCLsh) .and. Tc>-40.) then
                QCLih= min(QCLih*iEih, QIM(i,k))  !change Eih to 1. in CLih
                NCLih= min(NCLih*iEih, NYM(i,k))  !  "    "
                QCLsh= min(QCLsh*iEsh, QNM(i,k))  !change Esh to 1. in CLsh
                NCLsh= min(NCLsh*iEsh, NNM(i,k))  !  "    "
                tmp3 = QCLrh
                QCLrh= QHwet-(QCLch+QCLih+QCLsh)  !actual QCLrh minus QSHhr
                QSHhr= tmp3-QCLrh                 !QSHhr used here only
                NSHhr= DE(i,k)*QSHhr/(cmr*Drshed*Drshed*Drshed)
             else
                NSHhr= 0.
             endif
          ELSE
             QVDvh= 0.;   NVDvh= 0.;   NSHhr= 0.
          ENDIF

       ENDIF  ! ( if Tc<0C Block )

     !------------  End of source/sink term calculation  -------------!


     !----- Prevent mass transfer from accretion during melting:  ---!
     ! (only if exepcted NX (X=s,g,h) < 0 after source/sinks added)
     ! The purpose is to prevent mass tranfer from cloud to X and then
     ! to vapor (during "prevent overdepletion" code) if all of X
     ! would otherwise be completely depleted (e.g. due to melting).

       if (dblMom_s) then
         !estimate NN after S/S:
          tmp1= NN(i,k) +NsCNis -NVDvs -NCNsg -NMLsr -NCLss -NCLsr -NCLsh +NCLsrs
          if (tmp1<epsN) then
             QCLcs= 0.
             NCLcs= 0.
             QCLrs= 0.
             NCLrs= 0.
             if (Dsrs==1.) then
                QCLrs= 0.
                QCLsr= 0.
             endif
          endif
       endif
       if (dblMom_g) then
         !estimate NG after S/S:
          tmp2= NG(i,k) +NCNig +NCNsg -NCLgr -NVDvg -NMLgr +NCLirg +NCLsrg +NCLgrg -NCNgh
          if (tmp2<epsN) then
             QCLcg= 0.
             NCLcg= 0.
             QCLrg= 0.
             NCLrg= 0.
             if (Dirg==1.) then
                QCLri= 0.
                QCLir= 0.
             endif
             if (Dsrg==1.) then
                QCLrs= 0.
                QCLsr= 0.
             endif
          endif
       endif
       if (dblMom_h) then
         !estimate NH after S/S:
          tmp3= NH(i,k) +NhFZrh +NCNgh -NMLhr -NVDvh +NCLirh +NCLsrh +NCLgrh
          if (tmp3<epsN) then
             QCLch= 0.
             NCLch= 0.
             QCLrh= 0.
             NCLrh= 0.
             if (Dirh==1.) then
                QCLri= 0.
                QCLir= 0.
             endif
             if (Dsrh==1.) then
                QCLrs= 0.
                QCLsr= 0.
             endif
          endif
       endif
     !=====

     !-- Adjustment of source/sink terms to prevent  overdepletion: --!
       do niter= 1,2

          ! (1) Vapor:
          source= Q(i,k) +dim(-QVDvi,0.)+dim(-QVDvs,0.)+dim(-QVDvg,0.)+dim(-QVDvh,0.)
          sink  = QNUvi+dim(QVDvi,0.)+dim(QVDvs,0.)
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             QNUvi= ratio*QNUvi;   NNUvi= ratio*NNUvi
             if(QVDvi>0.) then
               QVDvi= ratio*QVDvi; NVDvi= ratio*NVDvi
             endif
             if(QVDvs>0.) then
               QVDvs=ratio*QVDvs;  NVDvs=ratio*NVDvs
             endif
             QVDvg= ratio*QVDvg;   NVDvg= ratio*NVDvg
             QVDvh= ratio*QVDvh;   NVDvh= ratio*NVDvh
          endif

          ! (2) Cloud:
          source= QC(i,k)
          sink  = QCLcs+QCLcg+QCLch+QFZci
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             QFZci= ratio*QFZci;   NFZci= ratio*NFZci
             QCLcs= ratio*QCLcs;   NCLcs= ratio*NCLcs
             QCLcg= ratio*QCLcg;   NCLcg= ratio*NCLcg
             QCLch= ratio*QCLch;   NCLch= ratio*NCLch
          endif

          ! (3) Rain:
          source= QR(i,k)+QMLsr+QMLgr+QMLhr+QMLir
          sink  = QCLri+QCLrs+QCLrg+QCLrh+QFZrh
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             QCLrg= ratio*QCLrg;   QCLri= ratio*QCLri;   NCLri= ratio*NCLri
             QCLrs= ratio*QCLrs;   NCLrs= ratio*NCLrs
             NCLrg= ratio*NCLrg;   QCLrh= ratio*QCLrh;   NCLrh= ratio*NCLrh
             QFZrh= ratio*QFZrh;   NrFZrh=ratio*NrFZrh;  NhFZrh=ratio*NhFZrh
             if (ratio==0.) then
                Dirg= 0.; Dirh= 0.; Dgrg= 0.; Dgrh= 0.
                Dsrs= 0.; Dsrg= 0.; Dsrh= 0.
             endif
          endif

          ! (4) Ice:
          source= QI(i,k)+QNUvi+dim(QVDvi,0.)+QFZci
          sink  = QCNig+QCNis+QCLir+dim(-QVDvi,0.)+QCLis+QCLig+QCLih+QMLir
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             QMLir= ratio*QMLir;    NMLir= ratio*NMLir
             if (QVDvi<0.) then
                QVDvi= ratio*QVDvi; NVDvi= ratio*NVDvi
             endif
             QCNis=  ratio*QCNis;   NiCNis= ratio*NiCNis;   NsCNis= ratio*NsCNis
             QCLir=  ratio*QCLir;   NCLir=  ratio*NCLir;    QCLig=  ratio*QCLig
             QCLis=  ratio*QCLis;   NCLis=  ratio*NCLis;    QCNig=  ratio*QCNig
             QCLih=  ratio*QCLih;   NCLih=  ratio*NCLih;    NCNig=  ratio*NCNig
             if (ratio==0.) then
                Dirg= 0.; Dirh= 0.
             endif
          endif

          ! (5) Snow:
          source= QN(i,k)+QCNis+dim(QVDvs,0.)+QCLis+Dsrs*(QCLrs+QCLsr)+QCLcs
          sink  = dim(-QVDvs,0.)+QCNsg+QMLsr+QCLsr+QCLsh
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             if(QVDvs<=0.) then
                QVDvs= ratio*QVDvs;   NVDvs= ratio*NVDvs
             endif
             QCNsg= ratio*QCNsg;   NCNsg= ratio*NCNsg;   QMLsr= ratio*QMLsr
             NMLsr= ratio*NMLsr;   QCLsr= ratio*QCLsr;   NCLsr= ratio*NCLsr
             QCLsh= ratio*QCLsh;   NCLsh= ratio*NCLsh
             if (ratio==0.) then
                Dsrs= 0.; Dsrg= 0.; Dsrh= 0.
             endif
          endif

          !  (6) Graupel:
          source= QG(i,k)+QCNig+QCNsg+dim(QVDvg,0.)+Dirg*(QCLri+QCLir)+Dgrg*(QCLrg+      &
                  QCLgr)+QCLcg+Dsrg*(QCLrs+QCLsr)+QCLig
          sink  = dim(-QVDvg,0.)+QMLgr+QCNgh+QCLgr
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             QVDvg= ratio*QVDvg;   NVDvg= ratio*NVDvg;   QMLgr= ratio*QMLgr
             NMLgr= ratio*NMLgr;   QCNgh= ratio*QCNgh;   NCNgh= ratio*NCNgh
             QCLgr= ratio*QCLgr;   NCLgr= ratio*NCLgr
             if (ratio==0.) then
                Dgrg= 0.; Dgrh= 0.
             endif
          endif

          !  (7) Hail:
          source= QH(i,k)+dim(QVDvh,0.)+QCLch+QCLrh+Dirh*(QCLri+QCLir)+QCLih+QCLsh+    &
                  Dsrh*(QCLrs+QCLsr)+QCNgh+Dgrh*(QCLrg+QCLgr)+QFZrh
          sink  = dim(-QVDvh,0.)+QMLhr
          sour  = max(source,0.)
          if(sink>sour) then
             ratio= sour/sink
             QVDvh= ratio*QVDvh;   NVDvh= ratio*NVDvh
             QMLhr= ratio*QMLhr;   NMLhr= ratio*NMLhr
          endif

       enddo
       !---------------  End of source/sink term adjustment  ------------------!

      !Compute N-tendencies for destination categories of 3-comp.freezing:
       NCLirg= 0.;  NCLirh= 0.;  NCLsrs= 0.;  NCLsrg= 0.
       NCLsrh= 0.;  NCLgrg= 0.;  NCLgrh= 0.

       if (QCLir+QCLri>0.) then
          tmp1  = max(Dr,Di)
          tmp2  = tmp1*tmp1*tmp1*PIov6
          NCLirg= Dirg*DE(i,k)*(QCLir+QCLri)/(deg*tmp2)
          NCLirh= Dirh*DE(i,k)*(QCLir+QCLri)/(deh*tmp2)
       endif

       if (QCLsr+QCLrs>0.) then
          tmp1  = max(Dr,Ds)
          tmp2  = tmp1*tmp1*tmp1*PIov6
          NCLsrs= Dsrs*DE(i,k)*(QCLsr+QCLrs)/(des*tmp2)
          NCLsrg= Dsrg*DE(i,k)*(QCLsr+QCLrs)/(deg*tmp2)
          NCLsrh= Dsrh*DE(i,k)*(QCLsr+QCLrs)/(deh*tmp2)
       endif

       if (QCLgr+QCLrg>0.) then
          tmp1  = max(Dr,Dg)
          tmp2  = tmp1*tmp1*tmp1*PIov6
          NCLgrg= Dgrg*DE(i,k)*(QCLgr+QCLrg)/(deg*tmp2)
          NCLgrh= Dgrh*DE(i,k)*(QCLgr+QCLrg)/(deh*tmp2)
       endif

       !========================================================================!
       !           Add all source/sink terms to all predicted moments:          !
       !========================================================================!

      !Diagnostic S/S terms:  (to facilitate output of 3D variables for diagnostics)
      !SS01(i,k)= QVDvs*idt  (e.g., for depositional growth rate of snow, kg kg-1 s-1)

       ! Q-Source/Sink Terms:
       Q(i,k) = Q(i,k)  -QNUvi -QVDvi -QVDvs -QVDvg -QVDvh
       QC(i,k)= QC(i,k) -QCLcs -QCLcg -QCLch -QFZci
       QR(i,k)= QR(i,k) -QCLri +QMLsr -QCLrs -QCLrg +QMLgr -QCLrh +QMLhr -QFZrh +QMLir
       QI(i,k)= QI(i,k) +QNUvi +QVDvi +QFZci -QCNis -QCLir -QCLis -QCLig                 &
                        -QMLir -QCLih +QIMsi +QIMgi -QCNig
       QG(i,k)= QG(i,k) +QCNsg +QVDvg +QCLcg -QCLgr -QMLgr -QCNgh -QIMgi +QCLig +QCNig   &
                        +Dirg*(QCLri+QCLir) +Dgrg*(QCLrg+QCLgr) +Dsrg*(QCLrs+QCLsr)
       QN(i,k)= QN(i,k) +QCNis +QVDvs +QCLcs -QCNsg -QMLsr -QIMsi -QCLsr +QCLis -QCLsh   &
                        +Dsrs*(QCLrs+QCLsr)
       QH(i,k)= QH(i,k) +Dirh*(QCLri+QCLir) -QMLhr +QVDvh +QCLch +Dsrh*(QCLrs+QCLsr)     &
                        +QCLih +QCLsh +QFZrh +QCLrh +QCNgh +Dgrh*(QCLrg+QCLgr)

       ! N-Source/Sink Terms:
       if (dblMom_c) NC(i,k)= NC(i,k) -NCLcs -NCLcg -NCLch -NFZci
       if (dblMom_r) NR(i,k)= NR(i,k) -NCLri -NCLrs -NCLrg -NCLrh +NMLsr +NMLgr +NMLhr   &
                                      -NrFZrh +NMLir +NSHhr
       if (dblMom_i) NY(i,k)= NY(i,k) +NNUvi +NVDvi +NFZci -NCLir -NCLis -NCLig -NCLih   &
                                      -NMLir +NIMsi +NIMgi -NCNig -NiCNis
       if (dblMom_s) NN(i,k)= NN(i,k) +NsCNis -NVDvs -NCNsg -NMLsr -NCLss -NCLsr -NCLsh  &
                                      +NCLsrs
       if (dblMom_g) NG(i,k)= NG(i,k) +NCNig +NCNsg -NCLgr -NVDvg -NMLgr +NCLirg +NCLsrg &
                                      +NCLgrg -NCNgh
       if (dblMom_h) NH(i,k)= NH(i,k) +NhFZrh +NCNgh -NMLhr -NVDvh +NCLirh +NCLsrh       &
                                      +NCLgrh

       T(i,k)= T(i,k)   +LFP*(QCLri+QCLcs+QCLrs+QFZci-QMLsr+QCLcg+QCLrg-QMLir-QMLgr      &
                        -QMLhr+QCLch+QCLrh+QFZrh) +LSP*(QNUvi+QVDvi+QVDvs+QVDvg+QVDvh)

     !Prevent overdepletion:
       IF (dblMom_c) THEN
         if(QC(i,k)<epsQ .or. NC(i,k)<epsN)    then
            tmp1= min(QC(i,k), maxQlClip)
            Q(i,k) = Q(i,k) + tmp1
            T(i,k) = T(i,k) - LCP*tmp1
            QC(i,k)= 0.
            NC(i,k)= 0.
         endif
       ELSE
         if(QC(i,k)<epsQ)    then
            Q(i,k) = Q(i,k) + QC(i,k)
            T(i,k) = T(i,k) - LCP*QC(i,k)
            QC(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_r) THEN
         if (QR(i,k)<epsQ .or. NR(i,k)<epsN)   then
            tmp1= min(QR(i,k), maxQlClip)
            Q(i,k) = Q(i,k) + tmp1
            T(i,k) = T(i,k) - LCP*tmp1
            QR(i,k)= 0.
            NR(i,k)= 0.
         endif
       ELSE
         if (QR(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QR(i,k)
            T(i,k) = T(i,k) - LCP*QR(i,k)
           QR(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_i) THEN
         if (QI(i,k)<epsQ .or. NY(i,k)<epsN)   then
            tmp1= min(QI(i,k), maxQiClip)
            Q(i,k) = Q(i,k) + tmp1
            T(i,k) = T(i,k) - LSP*tmp1
            QI(i,k)= 0.
            NY(i,k)= 0.
         endif
       ELSE
         if (QI(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QI(i,k)
            T(i,k) = T(i,k) - LSP*QI(i,k)
            QI(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_s) THEN
         if (QN(i,k)<epsQ .or. NN(i,k)<epsN)   then
            tmp1= min(QN(i,k), maxQiClip)
            Q(i,k) = Q(i,k) + tmp1
            T(i,k) = T(i,k) - LSP*tmp1
            QN(i,k)= 0.
            NN(i,k)= 0.
         endif
       ELSE
         if (QN(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QN(i,k)
            T(i,k) = T(i,k) - LSP*QN(i,k)
            QN(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_g) THEN
         if (QG(i,k)<epsQ .or. NG(i,k)<epsN)   then
            tmp1= min(QG(i,k), maxQiClip)
            Q(i,k) = Q(i,k) + tmp1
            T(i,k) = T(i,k) - LSP*tmp1
            QG(i,k)= 0.
            NG(i,k)= 0.
         endif
       ELSE
         if (QG(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QG(i,k)
            T(i,k) = T(i,k) - LSP*QG(i,k)
            QG(i,k)= 0.
         endif
       ENDIF

       IF (dblMom_h) THEN
         if (QH(i,k)<epsQ .or. NH(i,k)<epsN)   then
            tmp1= min(QH(i,k), maxQiClip)
            Q(i,k) = Q(i,k) + tmp1
            T(i,k) = T(i,k) - LSP*tmp1
            QH(i,k)= 0.
            NH(i,k)= 0.
         else if (QH(i,k)>epsQ .and. NH(i,k)>epsN) then
          !Conversion to graupel of hail is small:
            Dh= (DE(i,k)*QH(i,k)/NH(i,k)*icmh)**thrd
            if (Dh<Dh_min) then
               QG(i,k)= QG(i,k) + QH(i,k)
               NG(i,k)= NG(i,k) + NH(i,k)
               QH(i,k)= 0.
               NH(i,k)= 0.
            endif
         endif
       ELSE
         if (QH(i,k)<epsQ)   then
            Q(i,k) = Q(i,k) + QH(i,k)
            T(i,k) = T(i,k) - LSP*QH(i,k)
            QH(i,k)= 0.
         endif
       ENDIF
       Q(i,k)=  max(Q(i,k),0.)
       NY(i,k)= min(NY(i,k), Ni_max)

      ENDIF  !if (activePoint)
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------------!
  !                    End of ice phase microphysics (Part 2)                        !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !                       PART 3: Warm Microphysics Processes                        !
  !                                                                                  !
  !  Equations for warm-rain coalescence based on Cohard and Pinty (2000a,b; QJRMS)  !
  !  Condensation/evaportaion equations based on Kong and Yau (1997; Atmos-Ocean)    !
  !  Equations for rain reflectivity (ZR) based on Milbrandt and Yau (2005b; JAS)    !
  !----------------------------------------------------------------------------------!

  ! Part 3a - Warm-rain Coallescence:

 IF (warmphase_ON) THEN

  DO k= 2,nk
     DO i= 1,ni

        RCAUTR= 0.;  CCACCR= 0.;  Dc= 0.;  iLAMc= 0.;  L  = 0.
        RCACCR= 0.;  CCSCOC= 0.;  Dr= 0.;  iLAMr= 0.;  TAU= 0.
        CCAUTR= 0.;  CRSCOR= 0.;  SIGc= 0.;  DrINIT= 0.
        iLAMc3= 0.;  iLAMc6= 0.;  iLAMr3= 0.;  iLAMr6= 0.

        if (dblMom_r) then
           rainPresent= (QRM(i,k)>epsQ .and. NRM(i,k)>epsN)
        else
           rainPresent= (QRM(i,k)>epsQ)
        endif

        if (.not. dblMom_c) NCM(i,k)= N_c_SM
        if (QCM(i,k)>epsQ .and. NCM(i,k)>epsN) then
           iLAMc = iLAMDA_x(DE(i,k),QCM(i,k),1./NCM(i,k),icexc9,thrd)
           iLAMc3= iLAMc*iLAMc*iLAMc
           iLAMc6= iLAMc3*iLAMc3
           Dc    = iLAMc*(GC2*iGC1)**thrd
           SIGc  = iLAMc*( GC3*iGC1- (GC2*iGC1)*(GC2*iGC1) )**sixth
           L     = 0.027*DE(i,k)*QCM(i,k)*(6.25e18*SIGc*SIGc*SIGc*Dc-0.4)
           if (SIGc>SIGcTHRS) TAU= 3.7/(DE(i,k)*(QCM(i,k))*(0.5e6*SIGc-7.5))
        endif

        if (rainPresent) then
           if (dblMom_r) then
              Dr = Dm_x(DE(i,k),QRM(i,k),1./NRM(i,k),icmr,thrd)
             !Drop-size limiter [prevents initially large drops from melted hail]
              if (Dr>3.e-3) then
                 tmp1    = (Dr-3.e-3);  tmp2= (Dr/DrMAX); tmp3= tmp2*tmp2*tmp2
                 NRM(i,k)= NRM(i,k)*max((1.+2.e4*tmp1*tmp1),tmp3)
                 tmp1    = DE(i,k)*QRM(i,k)*icmr
                 Dr      = (tmp1/NRM(i,k))**thrd
              endif
           else
              NRM(i,k)= GR50*sqrt(sqrt(GR31*iGR34*DE(i,k)*QRM(i,k)*icmr))
              Dr = Dm_x(DE(i,k),QRM(i,k),1./NRM(i,k),icmr,thrd)
           endif
           iLAMr = iLAMDA_x(DE(i,k),QRM(i,k),1./NRM(i,k),icexr9,thrd)
           iLAMr3= iLAMr*iLAMr*iLAMr
           iLAMr6= iLAMr3*iLAMr3
        endif

        !  Autoconversion:
        if (QCM(i,k)>epsQ .and. SIGc>SIGcTHRS .and. autoconv_ON) then
           RCAUTR= min( max(L/TAU,0.), QCM(i,k)*idt )
           DrINIT= max(83.e-6, 12.6e-4/(0.5e6*SIGc-3.5))  !initiation regime Dr
           DrAUT = max(DrINIT, Dr)                     !init. or feeding DrAUT
           CCAUTR= RCAUTR*DE(i,k)/(cmr*DrAUT*DrAUT*DrAUT)

           ! ---------------------------------------------------------------------------- !
           ! NOTE: The formulation for CCAUTR here (dNr/dt|initiation) does NOT follow
           !       eqn (18) in CP2000a, but rather it comes from the F90 code provided
           !       by J-P Pinty (subroutine: 'rain_c2r2.f90').
           !       (See notes: 2001-10-17; 2001-10-22)
           !
           !       Similarly, the condition for the activation of accretion and self-
           !       collection depends on whether or not autoconversion is in the feeding
           !       regime (see notes 2002-01-07).  This is apparent in the F90 code, but
           !       NOT in CP2000a.
           ! ---------------------------------------------------------------------------- !

           ! cloud self-collection: (dNc/dt_autoconversion)   {CP eqn(25)}
           if (dblMom_c) CCSCOC= min(KK2*NCM(i,k)*NCM(i,k)*GC3*iGC1*iLAMc6, NCM(i,k)*    &
                                 idt)  !{CP00a eqn(25)}
        endif

        ! Accretion, rain self-collection, and collisional breakup:
        if (((QRM(i,k))>1.2*max(L,0.)*iDE(i,k).or.Dr>max(5.e-6,DrINIT)).and.rainAccr_ON  &
             .and. rainPresent) then

           !  Accretion:                                                      !{CP00a eqn(22)}
           if (QCM(i,k)>epsQ.and.L>0.) then
              if (Dr.ge.100.e-6) then
                 CCACCR = KK1*(NCM(i,k)*NRM(i,k))*(GC2*iGC1*iLAMc3+GR34*iGR31*iLAMr3)
                 RCACCR = cmr*iDE(i,k)*KK1*(NCM(i,k)*NRM(i,k))*iLAMc3*(GC3*iGC1*iLAMc3+  &
                          GC2*iGC1*GR34*iGR31*iLAMr3)
              else
                 CCACCR = KK2*(NCM(i,k)*NRM(i,k))*(GC3*iGC1*iLAMc6+GR37*iGR31*iLAMr6)

!                  RCACCR= cmr*iDE(i,k)*KK2*(NCM(i,k)*NRM(i,k))*iLAMc3*                  &
!                          (GC4*iGR31*iLAMc6+GC2*iGC1*GR37*iGR31*iLAMr6)
!++  The following calculation of RCACCR avoids overflow:
                 tmp1   = cmr*iDE(i,k)
                 tmp2   = KK2*(NCM(i,k)*NRM(i,k))*iLAMc3
                 RCACCR = tmp1 * tmp2
                 tmp1   = GC4*iGR31
                 tmp1   = (tmp1)*iLAMc6
                 tmp2   = GC2*iGC1
                 tmp2   = tmp2*GR37*iGR31
                 tmp2   = (tmp2)*iLAMr6
                 RCACCR = RCACCR * (tmp1 + tmp2)
!++
              endif
              CCACCR = min(CCACCR,(NC(i,k))*idt)
              RCACCR = min(RCACCR,(QC(i,k))*idt)
            endif

           if (dblMom_r) then
            !Rain self-collection:
              tmp1= NRM(i,k)*NRM(i,k)
              if (Dr.ge.100.e-6) then
                 CRSCOR= KK1*tmp1*GR34*iGR31*iLAMr3                        !{CP00a eqn(24)}
              else
                 CRSCOR= KK2*tmp1*GR37*iGR31*iLAMr6                        !{CP00a eqn(25)}
              endif
            !Raindrop breakup:                                             !{CP00a eqn(26)}
              Ec= 1.
              if (Dr >=  600.e-6) Ec= exp(-2.5e3*(Dr-6.e-4))
              if (Dr >= 2000.e-6) Ec= 0.
              CRSCOR= min(Ec*CRSCOR,(0.5*NR(i,k))*idt) !0.5 prevents depletion of NR
           endif

        endif  !accretion/self-collection/breakup

        ! Prevent overdepletion of cloud:
        source= QC(i,k)
        sink  = (RCAUTR+RCACCR)*dt
        if (sink>source) then
           ratio = source/sink
           RCAUTR= ratio*RCAUTR
           RCACCR= ratio*RCACCR
           CCACCR= ratio*CCACCR
        endif

        ! Apply tendencies:
        QC(i,k)= max(0., QC(i,k)+(-RCAUTR-RCACCR)*dt )
        QR(i,k)= max(0., QR(i,k)+( RCAUTR+RCACCR)*dt )
        F12KY(i,k) = RCAUTR + RCACCR                  !units: kg kg-1 s-1

        if (dblMom_c) NC(i,k)= max(0., NC(i,k)+(-CCACCR-CCSCOC)*dt )
        if (dblMom_r) NR(i,k)= max(0., NR(i,k)+( CCAUTR-CRSCOR)*dt )

        if (dblMom_r) then
           if (QR(i,k)>epsQ .and. NR(i,k)>epsN) then
              Dr = Dm_x(DE(i,k),QR(i,k),1./NR(i,k),icmr,thrd)
              if (Dr>3.e-3) then
                 tmp1= (Dr-3.e-3);   tmp2= tmp1*tmp1
                 tmp3= (Dr/DrMAX);   tmp4= tmp3*tmp3*tmp3
                 NR(i,k)= NR(i,k)*(max((1.+2.e4*tmp2),tmp4))
              elseif (Dr<Dhh) then
              !Convert small raindrops to cloud:
                 QC(i,k)= QC(i,k) + QR(i,k)
                 NC(i,k)= NC(i,k) + NR(i,k)
                 QR(i,k)= 0.;   NR(i,k)= 0.
              endif
           else
              QR(i,k)= 0.;   NR(i,k)= 0.
           endif  !(Qr,Nr>eps)
        endif

     ENDDO
  ENDDO

  ! Part 3b - Condensation/Evaporation:

  DO k=1,nk
     DO i=1,ni

        DEo     = DE(i,nk)
        gam     = sqrt(DEo*iDE(i,k))
        QSS(i,k)= sngl(FOQSA(T(i,k), PS(i)*S(i,k)))  ! Re-calculates QS with new T (w.r.t. liquid)
        if (QSS(i,k)>1.e-20) then
           ssat = Q(i,k)/QSS(i,k)-1.
        else
           ssat = 0.
        endif
        Tc      = T(i,k)-TRPL
        Cdiff   = max(1.62e-5, (2.2157e-5 + 0.0155e-5*Tc)) *1.e5/(S(i,k)*PS(i))
        MUdyn   = max(1.51e-5, (1.7153e-5 + 0.0050e-5*Tc))
        MUkin   = MUdyn*iDE(i,k)
        iMUkin  = 1./MUkin
        Ka      = max(2.07e-2, (2.3971e-2 + 0.0078e-2*Tc))
        ScTHRD  = (MUkin/Cdiff)**thrd ! i.e. Sc^(1/3)

        !Condensation/evaporation:
        ! Capacity of evap/cond in one time step is determined by saturation
        ! adjustment technique [Kong and Yau, 1997 App.A].  Equation for rain evaporation rate
        ! comes from Cohard and Pinty, 2000a.  Explicit condensation rate is not considered
        ! (as it is in Ziegler, 1985), but rather complete removal of supersaturation is assumed.

        X= Q(i,k)-QSS(i,k)
        if (dblMom_r) then
           rainPresent= (QR(i,k)>epsQ .and. NR(i,k)>epsN)
        else
           rainPresent= (QR(i,k)>epsQ)
        endif
        IF(X>0. .or. QC(i,k)>epsQ .or. rainPresent) THEN
           tmp1 = T(i,k)-35.86
           X    = X/(1.+ck5*QSS(i,k)/(tmp1*tmp1))
           if (X<(-QC(i,k))) then
              D= 0.
              if(rainPresent) then
                 if(QM(i,k)<QSW(i,k)) then
                    MUkin = (1.715e-5+5.e-8*Tc)*iDE(i,k)
                    iMUkin= 1./MUkin
                    if (dblMom_r) then
                        Dr   = Dm_x(DE(i,k),QR(i,k),1./NR(i,k),icmr,thrd)
                        iLAMr= iLAMDA_x(DE(i,k),QR(i,k),1./NR(i,k),icexr9,thrd)
                        LAMr = 1./iLAMr
                        !note: The following coding of 'No_r=...' prevents overflow:
                       !No_r_DM= NR(i,k)*LAMr**(1.+alpha_r))*iGR31
                        No_r_DM= sngl(dble(NR(i,k))*dble(LAMr)**dble(1.+alpha_r))*iGR31
                        No_r   = No_r_DM
                    else
                       iLAMr = sqrt(sqrt( (QR(i,k)*DE(i,k))/(GR34*cmr*No_r) ))
                    !note: No_r= No_r_SM   is already done (in Part 1)
                    endif
                    !note: There is an error in MY05a_eq(8) for VENTx (corrected in code)
                    VENTr= Avx*GR32*iLAMr**cexr5 + Bvx*ScTHRD*sqrt(gam*afr*iMUkin)*GR17* &
                           iLAMr**cexr6
                    ABw  = CHLC*CHLC/(Ka*RGASV*T(i,k)*T(i,k))+1./(DE(i,k)*(QSS(i,k))*    &
                           Cdiff)
                    QREVP= -dt*(PI2*ssat*No_r*VENTr/ABw)
                !!  QREVP= 0.  !to suppress evaporation of rain
                    if ((QR(i,k))>QREVP) then             !Note: QREVP is [(dQ/dt)*dt]
                       DEL= -QREVP
                    else
                       DEL= -QR(i,k)
                    endif
                    D= max(X+QC(i,k), DEL)
                    FEVPKY(i,k)=-1.*D/QR(i,k)
                 endif  !QM< QSM
              endif   !QR<eps & NR<eps
              X= D - QC(i,k)

              QR(i,k)= QR(i,k) + D
              if (QR(i,k)>0. .and. dblMom_r)                                              &
                   NR(i,k)= max(0.,NR(i,k)+D*NR(i,k)/QR(i,k)) !(dNr/dt)|evap
              ! The above expression of (dNr/dt)|evap is from Ferrier, 1994.
              ! In CP2000a, Nr is not affected by evap. (except if Qr goes to zero).
              QC(i,k)= 0.;   NC(i,k)= 0.
              T(i,k) = T(i,k) + LCP*X
              Q(i,k) = Q(i,k) - X

           else  ![if(X >= -QC)]

              ! Nucleation of cloud droplets:
              if (ssat>0. .and. WZ(i,k)>0.001 .and. dblMom_c)                             &
                   NC(i,k)= max(NC(i,k),NccnFNC(WZ(i,k),TM(i,k),HPS(i)*S(i,k),CCNtype))

              ! All supersaturation is removed (condensed onto cloud field).
              T(i,k)  = T(i,k)  + LCP*X
              Q(i,k)  = Q(i,k)  - X
              QC(i,k) = QC(i,k) + X
              if (dblMom_c) then
                  if (X<0.) then
                     if (QC(i,k)>0.) then
                        NC(i,k)= max(0., NC(i,k) + X*NC(i,k)/QC(i,k) ) !(dNc/dt)|evap
                     else
                        NC(i,k)= 0.
                     endif
                  endif
                  if (QC(i,k)>0..and.NC(i,k)==0.) NC(i,k)= 1.e7 !prevents non-zero_Q & zero_N
              endif

           endif

        ENDIF

       !Protect against negative values due to overdepletion:
        if (dblMom_r) then
            if (QR(i,k)<epsQ.or.NR(i,k)<epsN)  then
               Q(i,k) = Q(i,k) + QR(i,k)
               T(i,k) = T(i,k) - QR(i,k)*LCP
               QR(i,k)= 0.;  NR(i,k)= 0.
            endif
        else
            if (QR(i,k)<epsQ)  then
               Q(i,k) = Q(i,k) + QR(i,k)
               T(i,k) = T(i,k) - QR(i,k)*LCP
               QR(i,k)= 0.
            endif
        endif

     ENDDO
  ENDDO    !cond/evap [k-loop]

 ENDIF  !if warmphase_ON

  !----------------------------------------------------------------------------------!
  !                    End of warm-phase microphysics (Part 3)                       !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !                            PART 4:  Sedimentation                                !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !  Sedimentation is computed using a modified version of the box-Lagrangian        !
  !  scheme.  Sedimentation is only computed for columns containing non-zero         !
  !  hydrometeor quantities (at at least one level).                                 !
  !----------------------------------------------------------------------------------!

 IF (sedi_ON) THEN

   fluxM_r= 0.;  fluxM_i= 0.;  fluxM_s= 0.;  fluxM_g= 0.;  fluxM_h= 0.
   RT_rn1 = 0.;  RT_rn2 = 0.;  RT_fr1 = 0.;  RT_fr2 = 0.;  RT_sn1 = 0.
   RT_sn2 = 0.;  RT_sn3 = 0.;  RT_pe1 = 0.;  RT_pe2 = 0.;  RT_peL = 0.

!--  RAIN sedimentation:
   if (DblMom_r) then
     call SEDI_main_2(QR,NR,1,Q,T,DE,iDE,gamfact_r,epsQr_sedi,epsN,afr,bfr,cmr,dmr,      &
                      ckQr1,ckQr2,icexr9,LCP,ni,nk,VrMax,DrMax,dt,DZ,fluxM_r,ktop_sedi,  &
                      GRAV,massFlux3D=massFlux3D_r)
   else  !if DblMom_r
      call SEDI_main_1b(QR,1,T,DE,iDE,gamfact_r,epsQr_sedi,afr,bfr,icmr,dmr,ckQr1,       &
                        icexr9,ni,nk,VrMax,DrMax,dt,DZ,fluxM_r,No_r_SM,ktop_sedi,GRAV,   &
                        massFlux3D=massFlux3D_r)
   endif  !if DblMom_r


!--  ICE sedimentation:
   if (DblMom_i) then
     call SEDI_main_2(QI,NY,2,Q,T,DE,iDE,gamfact,epsQi_sedi,epsN,afi,bfi,cmi,dmi,ckQi1,  &
                      ckQi2,ckQi4,LSP,ni,nk,ViMax,DiMax,dt,DZ,fluxM_i,ktop_sedi,GRAV,    &
                      massFlux3D=massFlux3D_i)
   else
     call SEDI_main_1b(QI,2,T,DE,iDE,gamfact,epsQi_sedi,afi,bfi,icmi,dmi,ckQi1,ckQi4,    &
                     ni,nk,ViMax,DiMax,dt,DZ,fluxM_i,-99.,ktop_sedi,GRAV,                &
                     massFlux3D=massFlux3D_i)
   endif


!--  SNOW sedimentation:
   if (DblMom_s) then
     call SEDI_main_2(QN,NN,3,Q,T,DE,iDE,gamfact,epsQs_sedi,epsN,afs,bfs,cms,dms,ckQs1,  &
                      ckQs2,iGS20,LSP,ni,nk,VsMax,DsMax,dt,DZ,fluxM_s,ktop_sedi,GRAV,    &
                      massFlux3D=massFlux3D_s)
   else
     call SEDI_main_1b(QN,3,T,DE,iDE,gamfact,epsQs_sedi,afs,bfs,icms,dms,ckQs1,iGS20,    &
                      ni,nk,VsMax,DsMax,dt,DZ,fluxM_s,-99.,ktop_sedi,GRAV,massFlux3D=    &
                      massFlux3D_s)
   endif


!--  GRAUPEL sedimentation:
   if (DblMom_g) then
     call SEDI_main_2(QG,NG,4,Q,T,DE,iDE,gamfact,epsQg_sedi,epsN,afg,bfg,cmg,dmg,ckQg1,  &
                      ckQg2,ckQg4,LSP,ni,nk,VgMax,DgMax,dt,DZ,fluxM_g,ktop_sedi,GRAV,    &
                      massFlux3D=massFlux3D_g)
   else
     call SEDI_main_1b(QG,4,T,DE,iDE,gamfact,epsQg_sedi,afg,bfg,icmg,dmg,ckQg1,ckQg4,    &
                      ni,nk,VgMax,DgMax,dt,DZ,fluxM_g,No_g_SM,ktop_sedi,GRAV,            &
                      massFlux3D=massFlux3D_g)
   endif


!--  HAIL sedimentation:
   if (DblMom_h) then
     call SEDI_main_2(QH,NH,5,Q,T,DE,iDE,gamfact,epsQh_sedi,epsN,afh,bfh,cmh,dmh,ckQh1,  &
                      ckQh2,ckQh4,LSP,ni,nk,VhMax,DhMax,dt,DZ,fluxM_h,ktop_sedi,GRAV,    &
                      massFlux3D=massFlux3D_h)
   else
     call SEDI_main_1b(QH,5,T,DE,iDE,gamfact,epsQh_sedi,afh,bfh,icmh,dmh,ckQh1,ckQh4,    &
                      ni,nk,VhMax,DhMax,dt,DZ,fluxM_h,No_h_SM,ktop_sedi,GRAV,            &
                      massFlux3D=massFlux3D_h)
   endif

  !-- Diagnostics for chemistry module:
   do k=1,nk
      do i=1,ni
         PRFLXKY(i,k) = massFlux3D_r(i,k)
         SWFLXKY(i,k) = massFlux3D_i(i,k) + massFlux3D_s(i,k) + massFlux3D_g(i,k) + massFlux3D_h(i,k)
      enddo
   enddo
  !==

!=======  End of sedimentation for each category ========!

!---  Impose constraints on size distribution parameters ---!
   do k= 1,nk
      do i= 1,ni

       !snow:
         if (QN(i,k)>epsQ .and. NN(i,k)>epsN) then

         !Impose No_s max for snow: (assumes alpha_s=0.)
            iLAMs  = max( iLAMmin2, iLAMDA_x(DE(i,k),QN(i,k), 1./NN(i,k),iGS20,idms) )
            tmp1   = min(NN(i,k)/iLAMs,No_s_max)                 !min. No_s
            NN(i,k)= tmp1**(dms/(1.+dms))*(iGS20*DE(i,k)*QN(i,k))**(1./(1.+dms)) !impose Nos_max

         !Impose LAMDAs_min (by increasing LAMDAs if it is below LAMDAs_min2 [2xLAMDAs_min])
            iLAMs  = max( iLAMmin2, iLAMDA_x(DE(i,k),QN(i,k),1./NN(i,k),iGS20,idms) )
            tmp2   = 1./iLAMs                                   !LAMs before adjustment
           !adjust value of lamdas_min to be applied:
           !  This adjusts for multiple corrections (each time step).  The factor 0.6 was obtained by
           !  trial-and-error to ultimately give reasonable LAMs profiles, smooth and with min LAMs~lamdas_min
            tmp4   = 0.6*lamdas_min
            tmp5   = 2.*tmp4
            tmp3   = tmp2 + tmp4*(max(0.,tmp5-tmp2)/tmp5)**2.   !LAMs after adjustment
            tmp3   = max(tmp3, lamdas_min)                      !final correction
            NN(i,k)= NN(i,k)*(tmp3*iLAMs)**dms                  !re-compute NN after LAMs adjustment
         endif

      enddo !i-loop
   enddo !k-loop
!===

  !Compute melted (liquid-equivalent) volume fluxes [m3 (liquid) m-2 (sfc area) s-1]:
  !  (note: For other precipitation schemes in RPN-CMC physics, this is computed in 'vkuocon6.ftn')
   RT_rn1 = fluxM_r *idew
   RT_sn1 = fluxM_i *idew
   RT_sn2 = fluxM_s *idew
   RT_sn3 = fluxM_g *idew
   RT_pe1 = fluxM_h *idew

!----
  !Compute sum of solid (unmelted) volume fluxes [m3 (bulk hydrometeor) m-2 (sfc area) s-1]:
  !(this is the precipitation rate for UNmelted total snow [i+s+g])

  !    Note:  In 'calcdiagn.ftn', the total solid precipitation (excluding hail) SN is computed
  !           from the sum of the liq-eq.vol fluxes, tss_sn1 + tss_sn2 + tss_sn3.  With the
  !           accumulation of SND (in 'calcdiag.ftn'), the solid-to-liquid ratio for the total
  !           accumulated "snow" (i+s+g) can be compute as SND/SN.  Likewise, the instantaneous
  !           solid-to-liquid ratio of solid precipitation is computed (in 'calcdiag.ftn') as
  !           RS2L = RSND/RSN.

   do i= 1,ni

      fluxV_i= fluxM_i(i)*idei
      fluxV_g= fluxM_g(i)*ideg
     !Compute unmelted volume flux for snow:
         ! note: This is based on the strict calculation of the volume flux, where vol=(pi/6)D^3,
         !       and remains in the integral calculation Fv = int[V(D)*vol(D)*N(D)*dD].
         !       For a constant density (ice and graupel), vol(D) = m(D)/dex, dex comes out of
         !       integral and Fv_x=Fm_x/dex
         !       Optimized for alpha_s = 0.
      if (QN(i,nk)>epsQ .and. NN(i,nk)>epsN .and. fluxM_s(i)>0.) then
         tmp1= 1./iLAMDA_x(DE(i,nk),QN(i,nk),1./NN(i,nk),iGS20,idms) !LAMDA_s
         fluxV_s= fluxM_s(i)*rfact_FvFm*tmp1**(dms-3.)
      else
         fluxV_s=0.
      endif

     !total solid unmelted volume flux, before accounting for partial melting:
      tmp1= fluxV_i + fluxV_g + fluxV_s

     !liquid-fraction of partially-melted solid precipitation:
     !  The physical premise is that if QR>0, QN+QI+QG>0, and T>0, then QR
     !  originates from melting and can be used to estimate the liquid portion
     !  of the partially-melted solid hydrometeor.
      tmp2= QR(i,nk) + QI(i,nk) + QN(i,nk) + QG(i,nk)
      if (T(i,nk)>TRPL .and. tmp2>epsQ) then
         fracLiq= QR(i,nk)/tmp2
      else
         fracLiq= 0.
      endif

     !Tend total volume flux towards the liquid-equivalent as the liquid-fraction increases to 1:
      tmp3= RT_sn1(i) + RT_sn2(i) + RT_sn3(i)      !total liquid-equivalent volume flux    (RSN,  Fv_sol)
      RT_snd(i)= (1.-fracLiq)*tmp1 + fracLiq*tmp3  !total volume flux with partial melting (RSND, Fvsle_sol)
      !Note:  Calculation of instantaneous solid-to-liquid ratio [RS2L = RSND/RSN]
      !       is based on the above quantities and is done on 'calcdiag.ftn'.

   enddo  !i-loop
!====

 !++++
 ! Diagnose surface precipitation types:
 !
 ! The following involves diagnostic conditions to determine surface precipitation rates
 ! for various precipitation elements defined in Canadian Meteorological Operational Internship
 ! Program module 3.1 (plus one for large hail) based on the sedimentation rates of the five
 ! sedimenting hydrometeor categories.
 !
 ! With the diagnostics shut off (precipDiag_ON=.false.), 5 rates are just the 5 category
 ! rates, with the other 6 rates just 0.  The model output variables will have:
 !
 !   total liquid = RT_rn1 [RAIN]
 !   total solid  = RT_sn1 [ICE] + RT_sn2 [SNOW] + RT_sn3 [GRAUPEL] + RT_pe1 [HAIL]
 !
 ! With the diagnostics on, the 5 sedimentation rates are partitioned into 9 rates,
 ! with the following model output variable:
 !
 !  total liquid = RT_rn1 [liquid rain] + RT_rn2 [liquid drizzle]
 !  total solid  = RT_fr1 [freezing rain] + RT_fr2 [freezing drizzle] + RT_sn1 [ice crystals] +
 !                 RT_sn2 [snow] + RT_sn3 [graupel] + RT_pe1 [ice pellets] + RT_pe2 [hail]
 !
 ! NOTE:  - The above total liquid/solid rates are computed in 'calcdiag.ftn' (as R2/R4).
 !        - RT_peL [large hail] is a sub-set of RT_pe2 [hail] and is thus not added to the total
 !          solid precipitation.

!   call tmg_start0(98,'mmCalcDIAG')

   IF (precipDiag_ON) THEN
      DO i= 1,ni

         DE(i,nk)= S(i,nk)*PS(i)/(RGASD*T(i,nk))

      !rain vs. drizzle:
         if (DblMom_r) then
            N_r= NR(i,nk)
         else
            N_r= (No_r*GR31)**(3./(4.+alpha_r))*(GR31*iGR34*DE(i,nk)*QR(i,nk)*icmr)**    &
                  ((1.+alpha_r)/(4.+alpha_r))             !i.e. NR = f(No_r,QR)
         endif
         if (QR(i,nk)>epsQ .and. N_r>epsN) then
            Dm_r(i,nk)= (DE(i,nk)*icmr*QR(i,nk)/N_r)**thrd
            if (Dm_r(i,nk)>Dr_large) then  !Dr_large is rain/drizzle size threshold
               RT_rn2(i)= RT_rn1(i);   RT_rn1(i)= 0.
            endif
         endif

      !liquid vs. freezing rain or drizzle:
         if (T(i,nk)<TRPL) then
            RT_fr1(i)= RT_rn1(i);   RT_rn1(i)= 0.
            RT_fr2(i)= RT_rn2(i);   RT_rn2(i)= 0.
         endif

      !ice pellets vs. hail:
         if (T(i,nk)>(TRPL+5.0)) then
         !note: The condition (T_sfc<5C) for ice pellets is a simply proxy for the presence
         !      of a warm layer aloft, though which falling snow or graupel will melt to rain,
         !      over a sub-freezinglayer, where the rain will freeze into the 'hail' category
            RT_pe2(i)= RT_pe1(i);   RT_pe1(i)= 0.
         endif

      !large hail:
         if (QH(i,nk)>epsQ) then
            if (DblMom_h) then
               N_h= NH(i,nk)
            else
               N_h= (No_h_SM*GH31)**(3./(4.+alpha_h))*(GH31*iGH34*DE(i,nk)*QH(i,nk)*     &
                  icmh)**((1.+alpha_h)/(4.+alpha_h))   !i.e. Nh = f(No_h,Qh)
            endif
            Dm_h(i,nk)= Dm_x(DE(i,nk),QH(i,nk),1./N_h,icmh,thrd)
            if (DM_h(i,nk)>Dh_large) RT_peL(i)= RT_pe2(i)
            !note: large hail (RT_peL) is a subset of the total hail (RT_pe2)
         endif

      ENDDO
   ENDIF  !if (precipDiag_ON)
 !
 !++++

 ELSE

    massFlux3D_r= 0.
    massFlux3D_s= 0.

 ENDIF  ! if (sedi_ON)

 where (Q<0.) Q= 0.

 !-----------------------------------------------------------------------------------!
 !                     End of sedimentation calculations (Part 4)                    !
 !-----------------------------------------------------------------------------------!


 !===================================================================================!
 !                             End of microphysics scheme                            !
 !===================================================================================!

 !-----------------------------------------------------------------------------------!
 !                    Calculation of diagnostic output variables:                    !

 ! Compute effective radii for cloud and ice (to be passed to radiation scheme):
 !  - based of definition, r_eff = M_r(3)/M_r(2) = 0.5*M_D(3)/M_D(2),
 !    where the pth moment w.r.t. diameter, M_D(p), is given by MY2005a, eqn (2)
  do k = 1,nk
     do i = 1,ni

      !cloud:
        if (QC(i,k)>epsQ .and. NC(i,k)>epsN) then
          !hardcoded for alpha_c = 1. and mu_c = 3.
           iNC   = 1./NC(i,k)
           iLAMc = iLAMDA_x(DE(i,k),QC(i,k),iNC,icexc9,thrd)
           reff_c(i,k) = 0.664639*iLAMc
        else
           reff_c(i,k) = 0.
        endif

      !ice:
        if (QI(i,k)>epsQ .and. NY(i,k)>epsN) then
          !hardcoded for alpha_i = 0. and mu_i = 1.
           iNY   = 1./NY(i,k)
           iLAMi = max( iLAMmin2, iLAMDA_x(DE(i,k),QI(i,k),iNY,icexi9,thrd) )
           reff_i1(i,k) = 1.5*iLAMi
        else
           reff_i1(i,k) = 0.
        endif

      !snow:
        if (QN(i,k)>epsQ .and. NN(i,k)>epsN) then
          !hardcoded for alpha_s = 0. and mu_s = 1.
           iNN   = 1./NN(i,k)
           iLAMs  = max( iLAMmin2, iLAMDA_x(DE(i,k),QN(i,k),iNN,iGS20,idms) )
           reff_i2(i,k) = 1.5*iLAMs
        else
           reff_i2(i,k) = 0.
        endif

      !graupel:
        if (QG(i,k)>epsQ .and. NG(i,k)>epsN) then
          !hardcoded for alpha_g = 0. and mu_g = 1.
           iNG   = 1./NG(i,k)
           iLAMg = max( iLAMmin1, iLAMDA_x(DE(i,k),QG(i,k),iNG,iGG99,thrd) )
           reff_i3(i,k) = 1.5*iLAMg
        else
           reff_i3(i,k) = 0.
        endif

      !hail:
        if (QH(i,k)>epsQ .and. NH(i,k)>epsN) then
          !hardcoded for alpha_h = 0. and mu_h = 1.
           iNH   = 1./NH(i,k)
           iLAMh = max( iLAMmin1, iLAMDA_x(DE(i,k),QH(i,k),iNH,iGH99,thrd) )
           reff_i4(i,k) = 1.5*iLAMh
        else
           reff_i4(i,k) = 0.
        endif

     ! calculate a 'binary' cloud fraction (0 or 1) based on supersaturation
        cloud_bin(i,k)=0.
        cloudliq=0.
        cloudice=0.
        cloudsnow=0.
        if (QC(i,k).ge.epsQ .and. (QM(i,k)/QSW(i,k)-1.).gt.1.e-6) cloudliq=1.
        if (QI(i,k).ge.epsQ .and. reff_i1(i,k).lt.100.e-6) cloudice=1.
        if (QN(i,k).ge.epsQ .and. reff_i2(i,k).lt.100.e-6) cloudsnow=1.
        cloud_bin(i,k)=cloudice+cloudsnow+cloudliq
        cloud_bin(i,k)=max(min(cloud_bin(i,k),1.),0.)

     enddo
  enddo

  IF (calcDiag) THEN

   !For reflectivity calculations:
     ZEC= minZET
     cxr=            icmr*icmr   !for rain (and cloud)
     cxi= 1./fdielec*icmr*icmr   !for all frozen categories
     Gzc= (6.+alpha_c)*(5.+alpha_c)*(4.+alpha_c)/((3.+alpha_c)*(2.+alpha_c)*(1.+alpha_c))
     Gzr= (6.+alpha_r)*(5.+alpha_r)*(4.+alpha_r)/((3.+alpha_r)*(2.+alpha_r)*(1.+alpha_r))
     Gzi= (6.+alpha_i)*(5.+alpha_i)*(4.+alpha_i)/((3.+alpha_i)*(2.+alpha_i)*(1.+alpha_i))
     if (snowSpherical) then  !dms=3
        Gzs= (6.+alpha_s)*(5.+alpha_s)*(4.+alpha_s)/((3.+alpha_s)*(2.+alpha_s)*          &
             (1.+alpha_s))
     else                     !dms=2
        Gzs= (4.+alpha_s)*(3.+alpha_s)/((2.+alpha_s)*(1.+alpha_s))
     endif
     Gzg= (6.+alpha_g)*(5.+alpha_g)*(4.+alpha_g)/((3.+alpha_g)*(2.+alpha_g)*(1.+alpha_g))
     Gzh= (6.+alpha_h)*(5.+alpha_h)*(4.+alpha_h)/((3.+alpha_h)*(2.+alpha_h)*(1.+alpha_h))

     do k= 1,nk
       do i= 1,ni
           DE(i,k)= S(i,k)*PS(i)/(RGASD*T(i,k))
           tmp9= DE(i,k)*DE(i,k)

        !Compute N_x for single-moment categories:
           if (DblMom_c) then
              N_c= NC(i,k)
           else
              N_c= N_c_SM
           endif
           if (DblMom_r) then
              N_r= NR(i,k)
           else
              N_r= (No_r_SM*GR31)**(3./(4.+alpha_r))*(GR31*iGR34*DE(i,k)*QR(i,k)*icmr)** &
                   ((1.+alpha_r)/(4.+alpha_r))             !i.e. NR = f(No_r,QR)
           endif
           if (DblMom_i) then
              N_i= NY(i,k)
           else
              N_i= N_Cooper(TRPL,T(i,k))
           endif
           if (DblMom_s) then
              N_s= NN(i,k)
           else
              No_s= Nos_Thompson(TRPL,T(i,k))
              N_s = (No_s*GS31)**(dms/(1.+dms+alpha_s))*(GS31*iGS34*DE(i,k)*QN(i,k)*     &
                    icms)**((1.+alpha_s)/(1.+dms+alpha_s))
           endif
           if (DblMom_g) then
              N_g= NG(i,k)
           else
              N_g= (No_g_SM*GG31)**(3./(4.+alpha_g))*(GG31*GG34*DE(i,k)*QG(i,k)*icmg)**  &
                   ((1.+alpha_g)/(4.+alpha_g))             !i.e. NX = f(No_x,QX)
           endif
           if (DblMom_h) then
              N_h= NH(i,k)
           else
              N_h= (No_h_SM*GH31)**(3./(4.+alpha_h))*(GH31*iGH34*DE(i,k)*QH(i,k)*icmh)** &
                   ((1.+alpha_h)/(4.+alpha_h))             !i.e. NX = f(No_x,QX)
           endif

        !Total equivalent reflectivity:     (units of [dBZ])
           tmp1= 0.;  tmp2= 0.;  tmp3= 0.;  tmp4= 0.;  tmp5= 0.;   tmp6= 0.
           if (QC(i,k)>epsQ .and. N_c>epsN) tmp1 = cxr*Gzr*tmp9*QC(i,k)*QC(i,k)/N_c
           if (QR(i,k)>epsQ .and. N_r>epsN) tmp2 = cxr*Gzr*tmp9*QR(i,k)*QR(i,k)/N_r
           if (QI(i,k)>epsQ .and. N_i>epsN) tmp3 = cxi*Gzi*tmp9*QI(i,k)*QI(i,k)/N_i
           if (QN(i,k)>epsQ .and. N_s>epsN) tmp4 = cxi*Gzs*tmp9*QN(i,k)*QN(i,k)/N_s
           if (QG(i,k)>epsQ .and. N_g>epsN) tmp5 = cxi*Gzg*tmp9*QG(i,k)*QG(i,k)/N_g
           if (QH(i,k)>epsQ .and. N_h>epsN) tmp6 = cxi*Gzh*tmp9*QH(i,k)*QH(i,k)/N_h
          !Modifiy dielectric constant for melting ice-phase categories:
           if ( T(i,k)>TRPL) then
             tmp3= tmp3*fdielec
             tmp4= tmp4*fdielec
             tmp5= tmp5*fdielec
             tmp6= tmp6*fdielec
           endif
           ZET(i,k) = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6   != Zc+Zr+Zi+Zs+Zg+Zh
           if (ZET(i,k)>0.) then
              ZET(i,k)= 10.*log10((ZET(i,k)*Zfact))      !convert to dBZ
           else
              ZET(i,k)= minZET
           endif
           ZET(i,k)= max(ZET(i,k),minZET)
           ZEC(i)= max(ZEC(i),ZET(i,k))  !composite (max in column) of ZET

         !Mean-mass diameters:  (units of [m])
           Dm_c(i,k)= 0.;   Dm_r(i,k)= 0.;   Dm_i(i,k)= 0.
           Dm_s(i,k)= 0.;   Dm_g(i,k)= 0.;   Dm_h(i,k)= 0.
           if(QC(i,k)>epsQ.and.N_c>epsN) Dm_c(i,k)=Dm_x(DE(i,k),QC(i,k),1./N_c,icmr,thrd)
           if(QR(i,k)>epsQ.and.N_r>epsN) Dm_r(i,k)=Dm_x(DE(i,k),QR(i,k),1./N_r,icmr,thrd)
           if(QI(i,k)>epsQ.and.N_i>epsN) Dm_i(i,k)=Dm_x(DE(i,k),QI(i,k),1./N_i,icmi,thrd)
           if(QN(i,k)>epsQ.and.N_s>epsN) Dm_s(i,k)=Dm_x(DE(i,k),QN(i,k),1./N_s,icms,idms)
           if(QG(i,k)>epsQ.and.N_g>epsN) Dm_g(i,k)=Dm_x(DE(i,k),QG(i,k),1./N_g,icmg,thrd)
           if(QH(i,k)>epsQ.and.N_h>epsN) Dm_h(i,k)=Dm_x(DE(i,k),QH(i,k),1./N_h,icmh,thrd)

         !Supercooled liquid water:
           SLW(i,k)= 0.
           if (T(i,k)<TRPL) SLW(i,k)= DE(i,k)*(QC(i,k)+QR(i,k))   !(units of [kg/m3])

         !Visibility:
          !VIS1:  component through liquid cloud (fog) [m]
          ! (following parameterization of Gultepe and Milbrandt, 2007)
           tmp1= QC(i,k)*DE(i,k)*1.e+3             !LWC [g m-3]
           tmp2= N_c*1.e-6                         !Nc  [cm-3]
           if (tmp1>0.005 .and. tmp2>1.) then
              VIS1(i,k)= max(epsVIS,1000.*(1.13*(tmp1*tmp2)**(-0.51))) !based on FRAM [GM2007, eqn (4)
             !VIS1(i,k)= max(epsVIS,min(maxVIS, (tmp1*tmp2)**(-0.65))) !based on RACE [GM2007, eqn (3)
           else
              VIS1(i,k)= 3.*maxVIS  !gets set to maxVIS after calc. of VIS
           endif

          !VIS2: component through rain  !based on Gultepe and Milbrandt, 2008, Table 2 eqn (1)
           tmp1= massFlux3D_r(i,k)*idew*3.6e+6                        !rain rate [mm h-1]
           if (tmp1>0.01) then
              VIS2(i,k)= max(epsVIS,1000.*(-4.12*tmp1**0.176+9.01))   ![m]
           else
              VIS2(i,k)= 3.*maxVIS
           endif

          !VIS3: component through snow  !based on Gultepe and Milbrandt, 2008, Table 2 eqn (6)
           tmp1= massFlux3D_s(i,k)*idew*3.6e+6                        !snow rate, liq-eq [mm h-1]
           if (tmp1>0.01) then
              VIS3(i,k)= max(epsVIS,1000.*(1.10*tmp1**(-0.701)))      ![m]
           else
              VIS3(i,k)= 3.*maxVIS
           endif

          !VIS:  visibility due to reduction from all components 1, 2, and 3
          !      (based on sum of extinction coefficients and Koschmieders's Law)
           VIS(i,k) = min(maxVIS, 1./(1./VIS1(i,k) + 1./VIS2(i,k) + 1./VIS3(i,k)))
           VIS1(i,k)= min(maxVIS, VIS1(i,k))
           VIS2(i,k)= min(maxVIS, VIS2(i,k))
           VIS3(i,k)= min(maxVIS, VIS3(i,k))

        enddo  !i-loop
     enddo     !k-loop

    !Diagnostic levels:
     h_CB = noVal_h_XX   !height (AGL) of cloud base
     h_SN = noVal_h_XX   !height (AGL) of snow level [conventional snow (not just QN>0.)]
     h_ML1= noVal_h_XX   !height (AGL) of melting level [first 0C isotherm from ground]
     h_ML2= noVal_h_XX   !height (AGL) of melting level [first 0C isotherm from top]
                         ! note: h_ML2 = h_ML1 implies only 1 melting level
     tmp1= 1./GRAV
     do i= 1,ni
        CB_found= .false.;   SN_found= .false.;   ML_found= .false.
        do k= nk,2,-1
          !cloud base:
           if ((QC(i,k)>epsQ2.or.QI(i,k)>epsQ2) .and. .not.CB_found)  then
              h_CB(i) = GZ(i,k)*tmp1
              CB_found= .true.
           endif
          !snow level:
           if ( ((QN(i,k)>epsQ2 .and. Dm_s(i,k)>minSnowSize)  .or.                       &
                 (QG(i,k)>epsQ2 .and. Dm_g(i,k)>minSnowSize)) .and. .not.SN_found) then
              h_SN(i) = GZ(i,k)*tmp1
              SN_found= .true.
           endif
          !melting level: (height of lowest 0C isotherm)
           if (T(i,k)>TRPL .and. T(i,k-1)<TRPL .and. .not.ML_found) then
              h_ML1(i) = GZ(i,k)*tmp1
              ML_found= .true.
           endif
        enddo
     enddo
     do i= 1,ni
        ML_found= .false.  !from top
        do k= 2,nk
          !melting level: (height of highest 0C isotherm)
           if (T(i,k)>TRPL .and. T(i,k-1)<TRPL .and. .not.ML_found) then
              h_ML2(i) = GZ(i,k)*tmp1
              ML_found= .true.
           endif
        enddo
     enddo

  ENDIF
 !                                                                                   !

!-------------
!Convert N from #/m3 to #/kg:
! note: - at this point, NX is updated NX (at t+1); NXTEND is NX before S/S (at t*)
!       - NXM is no longer used (it does not need a unit conversion)
  do k= 1,nk
     DE(:,k) = S(:,k)*PS(:)/(RGASD*T(:,k))     !air density at time (t)
     iDE(:,k)= 1./DE(:,k)
  enddo
  NC= NC*iDE;   NCTEND= NCTEND*iDE
  NR= NR*iDE;   NRTEND= NRTEND*iDE
  NY= NY*iDE;   NYTEND= NYTEND*iDE
  NN= NN*iDE;   NNTEND= NNTEND*iDE
  NG= NG*iDE;   NGTEND= NGTEND*iDE
  NH= NH*iDE;   NHTEND= NHTEND*iDE
!=============

 !-----------------------------------------------------------------------------------!
 !   Compute the tendencies of  T, Q, QC, etc. (to be passed back to model dynamics) !
 !   and reset the fields to their initial (saved) values at time {*}:               !

      do k= 1,nk
         do i= 1,ni

            tmp1=T_TEND(i,k);  T_TEND(i,k)=(T(i,k) -T_TEND(i,k))*iDT;  T(i,k) = tmp1
            tmp1=Q_TEND(i,k);  Q_TEND(i,k)=(Q(i,k) -Q_TEND(i,k))*iDT;  Q(i,k) = tmp1
            tmp1=QCTEND(i,k);  QCTEND(i,k)=(QC(i,k)-QCTEND(i,k))*iDT;  QC(i,k)= tmp1
            tmp1=QRTEND(i,k);  QRTEND(i,k)=(QR(i,k)-QRTEND(i,k))*iDT;  QR(i,k)= tmp1
            tmp1=QITEND(i,k);  QITEND(i,k)=(QI(i,k)-QITEND(i,k))*iDT;  QI(i,k)= tmp1
            tmp1=QNTEND(i,k);  QNTEND(i,k)=(QN(i,k)-QNTEND(i,k))*iDT;  QN(i,k)= tmp1
            tmp1=QGTEND(i,k);  QGTEND(i,k)=(QG(i,k)-QGTEND(i,k))*iDT;  QG(i,k)= tmp1
            tmp1=QHTEND(i,k);  QHTEND(i,k)=(QH(i,k)-QHTEND(i,k))*iDT;  QH(i,k)= tmp1

            if (DblMom_c) then
             tmp1=NCTEND(i,k);  NCTEND(i,k)=(NC(i,k)-NCTEND(i,k))*iDT; NC(i,k)= tmp1
            endif
            if (DblMom_r) then
             tmp1=NRTEND(i,k);  NRTEND(i,k)=(NR(i,k)-NRTEND(i,k))*iDT; NR(i,k)= tmp1
            endif
            if (DblMom_i) then
             tmp1=NYTEND(i,k);  NYTEND(i,k)=(NY(i,k)-NYTEND(i,k))*iDT; NY(i,k)= tmp1
            endif
            if (DblMom_s) then
             tmp1=NNTEND(i,k);  NNTEND(i,k)=(NN(i,k)-NNTEND(i,k))*iDT; NN(i,k)= tmp1
            endif
            if (DblMom_g) then
             tmp1=NGTEND(i,k);  NGTEND(i,k)=(NG(i,k)-NGTEND(i,k))*iDT; NG(i,k)= tmp1
            endif
            if (DblMom_h) then
             tmp1=NHTEND(i,k);  NHTEND(i,k)=(NH(i,k)-NHTEND(i,k))*iDT; NH(i,k)= tmp1
            endif

         enddo
      enddo
 !                                                                                   !
 !-----------------------------------------------------------------------------------!
 END SUBROUTINE MYDMOM_MAIN
 !___________________________________________________________________________________!

   real function des_OF_Ds(Ds_local,desMax_local,eds_local,fds_local)
   !Computes density of equivalent-volume snow particle based on (pi/6*des)*Ds^3 = cms*Ds^dms
      real :: Ds_local,desMax_local,eds_local,fds_local
!     des_OF_Ds= min(desMax_local, eds_local*Ds_local**fds_local)
      des_OF_Ds= min(desMax_local, eds_local*exp(fds_local*log(Ds_local)))   !IBM optimization
   end function des_OF_Ds


   real function Dm_x(DE_local,QX_local,iNX_local,icmx_local,idmx_local)
   !Computes mean-mass diameter
      real :: DE_local,QX_local,iNX_local,icmx_local,idmx_local
     !Dm_x = (DE_local*QX_local*iNX_local*icmx_local)**idmx_local
      Dm_x = exp(idmx_local*log(DE_local*QX_local*iNX_local*icmx_local))     !IBM optimization
   end function Dm_x


   real function iLAMDA_x(DE_local,QX_local,iNX_local,icex_local,idmx_local)
   !Computes 1/LAMDA ("slope" parameter):
      real :: DE_local,QX_local,iNX_local,icex_local,idmx_local
     !iLAMDA_x = (DE_local*QX_local*iNX_local*icex_local)**idmx_local
      iLAMDA_x = exp(idmx_local*log(DE_local*QX_local*iNX_local*icex_local)) !IBM optimization
   end function


   real function N_Cooper(TRPL_local,T_local)
   !Computes total number concentration of ice as a function of temperature
   !according to parameterization of Cooper (1986):
      real :: TRPL_local,T_local
      N_Cooper= 5.*exp(0.304*(TRPL_local-max(233.,T_local)))
   end function N_Cooper

   real function Nos_Thompson(TRPL_local,T_local)
   !Computes the snow intercept, No_s, as a function of temperature
   !according to the parameterization of Thompson et al. (2004):
      real :: TRPL_local,T_local
      Nos_Thompson= min(2.e+8, 2.e+6*exp(-0.12*min(-0.001,T_local-TRPL_local)))
   end function Nos_Thompson

!===================================================================================================!

END MODULE my_dmom_mod


