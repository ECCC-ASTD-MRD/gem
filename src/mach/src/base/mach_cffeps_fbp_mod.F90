!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_cffeps_fbp_mod.ftn90
! Creation         : Kerry Anderson - Jan 6, 2011
!                  : (C To Fortran)
!                    A. Akingunola, J. Chen, and P. Makar - Fall 2018
! Description      : Fire Behavior Prediction System for CFFEPS
!/* This subroutine represents equations from
!
!Wotton, B.M.; Alexander, M.E.; Taylor, S.W.  2009,.  Updates and revisions to the
!    1992 Canadian Forest Fire Behavior Prediction System .  Natural Resources
!    Canada, Canadian Forest Service, Great Lakes Forestry Centre, Sault Ste. Marie,
!    Ontario, Canada.  Infomration Report GLC-X-10, 45 p.
!
!Updates are indicated with a "- 2009" in the equation number comment on the right
!
!   Discrepancies between BMW's fbp.c and ST-X-3
!
!   eqn 32: BMW has 33.5 instead of 35.5 (see ROScalc and Slopecalc)
!   eqn 57: BMW had forced SFC to 2.0 for C6 (this is now gone)
!   b[O1b]: BMW has 0.0829 instead of 0.0310
!   c[O1a]: BMW has 1.41 instead of 1.4
!*/
!
!// Note: There is no attempt to account for greened-up deciduous (D2)
!============================================================================
!
 module mach_cffeps_fbp_mod

   use mach_cffeps_mod,  only: fbp_type, feps_type, cffeps_fueltypes,  &
                               fuel_size, max_fuels, cffeps_fc_method, &
                               cffeps_buieff
   
   implicit none
   private
   public  :: cffeps_fbpcalc, FMCcalc
   
   real(kind=4), parameter :: pi = 3.1415926, NoBUI = -1.0
   
   real(kind=4), dimension(max_fuels), parameter :: a = &
                    (/  90.,  110.,   110.,  110.,   30.,   30., &
                        45.,   30.,   30.,     0.,    0.,  120., &
                       100.,   75.,   40.,   55.,  190.,   250., &
                         0.,   0. /)
   real(kind=4), dimension(max_fuels), parameter :: b = &
                    (/.0649, .0282,  .0444, .0293, .0697, .0800, &
                      .0305, .0232, .0232,     0.,    0., .0572, &
                      .0404, .0297, .0438, .0829, .0310,  .0350, &
                         0.,   0.  /)
   real(kind=4), dimension(max_fuels), parameter :: c = &
                    (/  4.5,   1.5,    3.0,   1.5,   4.0,   3.0, &
                        2.0,   1.6,    1.6,    0.,    0.,   1.4, &
                       1.48,   1.3,    1.7,   3.2,   1.4,   1.7, &
                         0.,   0. /)
   real(kind=4), dimension(max_fuels), parameter :: q = &
                    (/ 0.90,  0.70,   0.75,  0.80,  0.80,  0.80, &
                       0.85,  0.90,   0.90,    .8,    .8,    .8, &
                         .8,  0.75,   0.75,  0.75,  1.00,  1.00, &
                         0.,   0.  /)
   real(kind=4), dimension(max_fuels), parameter :: BUIo = &
                    (/ 72.,   64.,   62.,    66.,   56.,   62.,  &
                      106.,   32.,   32.,    50.,   50.,   50.,  &
                       50.,   38.,   63.,    31.,    1.,    1.,  &
                        0.,    0./)
   real(kind=4), dimension(max_fuels), parameter :: CBHs = &
                    (/  2.,    3.,     8.,    4.,   18.,    7.,   &
                       10.,    0.,     0.,    6.,    6.,    6.,   &
                        6.,    0.,     0.,    0.,    0.,    0.,   &
                        0.,    0./)
   real(kind=4), dimension(max_fuels), parameter :: CFLs = &
                    (/ .75,   .80,   1.15,  1.20,  1.20,  1.80,   &
                       .50,    0.,     0.,   0.8,   0.8,   0.8,   &
                       0.8,    0.,     0.,    0.,    0.,    0.,   &
                        0.,    0./)
                        
!/* These values are pulled from the SFC calculations */
   real(kind=4), dimension(max_fuels), parameter :: SFLs = &
                    (/ 1.5,   5.0,    5.0,   5.0,   5.0,   5.0,   &
                       0.0,   1.5,    1.5,  3.25,  3.25,   5.0,   &
                       5.0,   0.0,    0.0,   0.0,   0.0,   0.0,   &
                       0.0,   0.0/)
   real(kind=4), dimension(max_fuels), parameter :: FFLs = &
                    (/ 0.0,   0.0,    0.0,   0.0,   0.0,   0.0,   &
                       2.0,   0.0,    0.0,   0.0,   0.0,   0.0,   &
                       0.0,   4.0,   10.0,  12.0,   0.0,   0.0,   &
                       0.0,   .0 /)
   real(kind=4), dimension(max_fuels), parameter :: WFLs = &
                    (/ 0.0,   0.0,    0.0,   0.0,   0.0,   0.0,   &
                       1.5,   0.0,    0.0,   0.0,   0.0,   0.0,   &
                       0.0,   4.0,    6.0,  20.0,   0.0,   0.0,   &
                       0.0,   0.0 /)

!/* These values taken from Anderson 2002 */
   real(kind=4), dimension(max_fuels), parameter :: depths = &
                    (/ 3.4,  10.0,    6.5,   4.6,   5.0,   5.0,   &
                       5.0,   2.4,    2.4,   5.0,   5.0,   7.5,   &
                       7.5,   7.4,    7.4,   7.4,   0.0,   0.0,   &
                       0.0,   0.0/) !  /* values corrected 2020-01-13 */
   real(kind=4), dimension(max_fuels), parameter :: rhoBs = &
                    (/0.045, 0.034,   0.026, 0.028, 0.093, 0.05,  &
                       0.02, 0.108,   0.108, 0.061, 0.061, 0.106, &
                      0.106, 0.078,   0.132,   0.1,   0.0,   0.0, &
                      0.0,   0.0/) !  /* values corrected 2020-01-13 */
   real(kind=4), dimension(max_fuels), parameter :: Ris = &
                    (/0.05,   0.0,   0.15,  0.15,  0.15,  0.15,   &
                      0.15,  0.59,   0.59,  0.25,  0.25,  0.15,   &
                      0.15,  0.15,   0.15,  0.15,   0.0,   0.0,   &
                       0.0,   0.0/)

   contains
  
!    double t,               /* Hours since ignition */
!           Depth,           /* Forest Floor Depth [cm] */
!           BD,              /* Bulk Density [gm/cm^3] */!           ISI,             /* ISI */
!           WS,              /* wind speed [kmh] */
!           WD,              /* wind direction [degrees] */
!           GS,              /* Slope [percent] */
!           Aspect,          /* Aspect [degrees] */
!           PC,              /* Percent Confier for M1/M2 */
!           PDF,             /* Percent Dead Fir for M3/M4 */
!           Cured,           /* Percent Cured for O1a/O1b (85% default) */
!           GFL,             /* Grass Fuel Load [kg/m^2] (0.3 default) */
!           CBH,             /* Crown to Base Height [m] (FBP defaults)*/
!           CFL,             /* Crown Fuel Load [kg/m^2] (FBP defaults) */
!           SFL,             /* Surface Fuel Load [kg/m^2] (FBP defaults) */
!           FFL,             /* Fine Fuel Load [kg/m^2] (FBP defaults) */
!           WFL,             /* Woody Fuel Load [kg/m^2] (FBP defaults) */
!           LFL,             /* Litter Fuel Load [kg/m^2] -- 2019-11-20 */
!           DFL,             /* Duff (Fermentation) Fuel Load [kg/m^2] -- 2019-11-20 */
!           HFL,             /* Humus Fuel Load [kg/m^2] -- 2019-11-20 */
!           FMC,             /* Foliar Moisture Content if known */
!           SH,              /* C6 Stand Height [m] - 2009 */
!           SD,              /* C6 Stand Density [stems/ha] - 2009 */
!           theta,           /* elliptical direction of calculation */
!           SFC,             /* Surface Fuel Consumption [kg/m^2] */
! 
! /* outputs */
!           ROS,             /* Rate of Spread [m/min] */
!           FROS,            /* Flank rate of Spread [m/min] */
!           BROS,            /* Back Rate of Spread [m/min] */
!           TROS,            /* Rate of Spread at angle theta [m/min] */
!           CFB,             /* Crown Fraction Burned */
!           HFI,             /* Head Fire Intensity [kW/m] */
!           TFC,             /* Total Fuel Consumption [kg/m^2]  */
!
!*******************************************************************************!
!* Main FBP routine
   subroutine cffeps_fbpcalc(ifuel, ws, wd, ffmc, feps, fbp, istep)
      implicit none
       
      integer(kind=4), intent(in)    :: ifuel
      integer(kind=4), intent(in)    :: istep
      real(kind=4),    intent(in)    :: ws, wd
      real(kind=4),    intent(inout) :: ffmc
      type(feps_type), intent(inout) :: feps
      type(fbp_type),  intent(out)   :: fbp
      
      real(kind=4) :: fmc, bui, dmc, dc, pc, pdf, cured, cbh
      real(kind=4) :: sfl, cfl, ffl, wfl, lfl, dfl, gs, aspect, sh, gfl, sd
      real(kind=4) :: cfc, sfc, tfc, cfb, ros, fros, bros, tros, theta, hfi
      
      real(kind=4) :: sfc_c2, sfc_d1, lfc, dfc, hfc, m, depth, bd, hfl, sfcmax
      real(kind=4) :: work, zero, cent, rsz, rsf, rsf_c2, rsf_d1, sf
      real(kind=4) :: cf, isz, isf, isf_c2, isf_d1, ecc
      integer(kind=4) :: ifuel_c2, ifuel_d1, method

      character(len=FUEL_SIZE) :: fueltype
      real(kind=4) :: isi, lb, wsv, saz, wse, ff
      real(kind=4) :: bfw, bisi, raz, wsx, wsy, waz
!      real(kind=4) :: alpha, lbt
      
      fueltype = cffeps_fueltypes(ifuel)
      
      fmc   = feps%fmc
      dmc   = feps%dmc
      dc    = feps%dc
      
      bui   = feps%bui
      
      pc    = feps%pc
      pdf   = feps%pdf
      cured = feps%curing
      cbh   = feps%cbh
      cfl   = feps%cfl
      gfl   = feps%gfl
      gs    = feps%gs
      sd    = feps%sd
      sh    = feps%sh
      aspect= feps%aspect * pi/ 180.0

      ffl   = feps%ffl
      sfl   = feps%sfl
      wfl   = feps%wfl
      lfl   = feps%lfl
      dfl   = feps%dfl
      
      method = cffeps_fc_method
      
      ifuel_c2 = 2   ! ifuel for "C2"
      ifuel_d1 = 8   ! ifuel for "D1"
         
      ! These are still set to a fixed value ATM
      depth  = 7.0
      bd     = 0.035
      hfl    = 0.0
      theta  = 0.0

      if (istep == 0 .and. lfl >= 0.0) then
         sfl = lfl + dfl + hfl
      end if

      ! FFMC consistency check
      if (ffmc < 0.0 .or. ffmc > 101.0) ffmc = 0.0
!
      if (dmc < 0.0 .or. dmc > 400.0) dmc = 0.0
      if (dc < 0.0 .or. dc > 1200.0) dc = 0.0
!
      if (sfl < 0.0) method = 1
      
      ! Set defaults, if input values are undefined
!      if (t < 0.0) t = -t
!/* Convert time from hours to minutes */
!      t = t * 60.0
      
      if (istep == 0) then
         if (lfl > 0.0 .or. dfl > 0.0) method = 4
         
         if (method == 1) then
            cured  = -1.0
            gs     = -1.0
            aspect = -pi / 180.0
            pc     = -1.0
            pdf    = -1.0 
            cbh    = -1.0
            gfl    = -1.0
            sfl    = -1.0
            cfl    = -1.0
            ffl    = -1.0
            wfl    = -1.0
            sh     = -1.0
            sd     = -1.0
         end if
         
         if (pc < 0.0) pc = 50.0
         if (pdf < 0.0) pdf = 35.0
         if (cured < 0.0) cured = 95.0
         if (gs < 0.0 .or. gs > 200.0 .or. &
             aspect < 2.0 * pi .or. aspect > 2.0 * pi) gs = 0.0
         if (gfl < 0.0) gfl = 0.35

!/* BUIcalc: Calculate the BUI value from DMC and DC values.  Determine whether
!       the BUI effect is being used. */  --- updated fix 2015-11-09 KRA
         if ((dmc * dc) == 0.0) then
            bui = 0.0
         else
            if (dmc <= (0.4 * dc)) then      ! /* 27a */
               bui = 0.8 * dmc * dc / (dmc + 0.4 * dc)
            else                             ! /* 27b */
               bui = dmc - (1.0 - 0.8 * dc / (dmc + 0.4 * dc)) *     &
                     (0.92 + ((0.0114 * dmc)**1.7))
            end if
            if (bui < 0.0 .or. bui > 1000.0) &
               bui = 0.0    !/* This turns off BUI effect */
         end if

!/* presently, we do not accept a zero CBH; use near zero if necessary */
         if (cbh <= 0.0 .or. cbh > 50.0) then
            if (ifuel == 6 .and. sd > 0.0 .and. sh > 0.0) then
               cbh = -11.2 + 1.06 * sh + 0.00170 * sd     !  /* 91 */
               cbh = max(0.0, cbh)
            else
               cbh = CBHs(ifuel)
            end if
         end if     
!
      ! Default FBP SFLs, WFLs, FFLs (e.g. SFL[C2, C3, C4, ...] = 5.0 kg/m^2
      ! OR; use the spatial SFL, WFL, FFL, GFL as provided by Dan Thompson
         if (cfl <= 0.0 .or. cfl > 2.0) then
            cfl = CFLs(ifuel)
         end if
      
!/* Ditto for FFL */
         if (ffl <= 0.0 .or. ffl > 20.0 .or. method == 1) then
            ffl = FFLs(ifuel)
         end if
      
!/* Ditto for SFL */
         if (sfl <= 0.0 .or. sfl > 200.0 .or. method == 1) then
            sfl = SFLs(ifuel)
         end if
      
!/* Ditto for WFL */
         if (wfl <= 0.0 .or. wfl > 50.0 .or. method == 1) then
            wfl = WFLs(ifuel)
         end if
!
! Update FBP parameters in the saved FEPS structure
         feps%bui    = bui
         feps%curing = cured
         feps%gfl    = gfl
         feps%gs     = gs
         feps%sd     = sd
         feps%sh     = sh      
         feps%aspect = aspect * 180 / pi
         feps%pc     = pc
         feps%pdf    = pdf
         feps%cbh    = cbh
         feps%cfl    = cfl
         feps%sfl    = sfl
         feps%ffl    = ffl
         feps%wfl    = wfl

      end if
      
!/* presently, we do not accept a zero BD,; use near zero if necessary */
      if (bd <= 0.0 .or. bd > 0.50) then 
         bd = 1000.0 * rhoBs(ifuel)
      else
         bd = 1000.0 * bd   !      // NOTE: BD used in FFFC is kg/m^3 (not g/cm^3)
      end if
  
!/* presently, we do not accept a zero Depth,; use near zero if necessary */
      if (depth <= 0.0 .or. depth > 50.0) then
         depth = depths(ifuel)
      end if

      if (ifuel == 17 .or. ifuel == 18) then
         sfcmax = 0.35  ! default grass fuel load
      else
         sfcmax = SFLs(ifuel) + WFLs(ifuel) + FFLs(ifuel)         
      end if
      
!****************************************************************************** 
!/* Surface Fuel Consumption (SFC) calculation */ - SFCcalc
!
      if (method < 4) then
         sfc = SFCcalc(ifuel, bui, ffmc, pc, gfl, sfl, ffl, wfl)
         ! set maximum SFC to the default FBP SFLs
         if (method == 3) sfc = min(sfc, sfcmax)
      !
      else if (method == 4) then ! // use de Groot et al (2007) with fuels
         if (ifuel <= 4 .or. ifuel == 8) then
         ! 'C1', 'C2', 'C3', 'C4', or 'D1'
            sfc = FFFCcalc(ifuel, depth, bd, sfl, ffmc, dmc, dc, bui)
         else if (ifuel == 10 .or. ifuel == 11) then
         ! 'M1' or 'M2'
            sfc_c2 = FFFCcalc(ifuel_c2, depth, bd, sfl, ffmc, dmc, dc, bui)
            sfc_d1 = FFFCcalc(ifuel_d1, depth, bd, sfl, ffmc, dmc, dc, bui)
            sfc = (pc * 1.0e-2 * sfc_c2) + ((100.0 - pc) * 1.0e-2 * sfc_d1)
                                                               ! /* 17 */
         else
            ! default calculations (method==2)
            sfc = SFCcalc(ifuel, bui, ffmc, pc, gfl, sfl, ffl, wfl)
         end if
      !
      else if (method == 5) then ! // use de Groot et al (2007) without fuels
         sfc = FFFCcalc(20, depth, bd, sfl, ffmc, dmc, dc, bui)
      !
      else ! // (method>=6) totally experimental
         lfl = max(lfl, 0.0)
         dfl = max(dfl, 0.0)
         hfl = max(hfl, 0.0)
         
         lfc = lfl
         ! /* substitute DMC for BUI */
         dfc = SFCcalc(ifuel, dmc, ffmc, pc, gfl, sfl, ffl, wfl)
         hfc = 0.0
         if (dfc > dfl) then
            dfc = dfl
            hfc = SFCcalc(ifuel, 0.4 * dc, ffmc, pc, gfl, sfl, ffl, wfl) - dfc
            hfc = max(hfc, 0.0)
            if (hfc > hfl) hfc = hfl
         end if
         sfc = lfc + dfc + hfc

         if (istep == 0) then
            feps%lfl = lfl
            feps%dfl = dfl
         end if
      end if
!****************************************************************************** 

      if (.not. cffeps_buieff) bui = -1.0 !/* This turns off BUI effect */

      m = 147.2 * (101.0 - ffmc) / (59.5 + ffmc)                ! /* 46 */
      ff = 91.9 * exp(-0.1386 * m) * (1.0 + (m**5.31) / 4.93e7) ! /* 45 */
!
!/* Corrections to reorient SAZ */
      waz = wd + pi
      if (waz > 2.0 * pi) waz = waz - 2.0 * pi

!/* nb: BMW's data set appears to have aspect not saz */
      saz = aspect + pi
      if (saz > 2.0 * pi) saz = saz - 2.0 * pi
!
!****************************************************************************** 
!/* Effect of Slope on Rate of Spread */ - slopecalc
!
      if (gs > 0.0 .and. ffmc > 0.) then
         if (gs >= 70.0) then
            sf = 10.0
         else
            sf = exp(3.533 * (gs * 1.0e-2)**1.2)          !  /* 39 */
         end if

         isz = ISIcalc(ff, zero)
         rsz = ROScalc(ifuel, isz, NoBUI, fmc, sfc, pc, pdf, cured, cbh)
         rsf = rsz * sf                                   !  /* 40 */
         
         zero = 0.0
         if (ifuel <= 8 .or. ifuel == 14 .or. ifuel == 15 .or. ifuel == 16) then
         ! case ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "D1", "S1", "S2", "S3")
            work = max(0.01, 1.0 - (rsf / a(ifuel))**(1.0 / c(ifuel)))
            isf = log(work) / (-b(ifuel))                 ! /* 41(a,b) - 2009 */
         else if(ifuel == 10 .or. ifuel == 11) then
         ! case ("M1", "M2")
            rsz = ROScalc(ifuel_c2, isz, NoBUI, fmc, sfc, pc, pdf, cured, cbh)
            rsf_c2 = rsz * sf                             ! /* 40 */ 
            
            rsz = ROScalc(ifuel_d1, isz, NoBUI, fmc, sfc, pc, pdf, cured, cbh)
            rsf_d1 = rsz * sf                             ! /* 40 */
              
            work = max(0.01, 1.0 - (rsf_c2 / a(ifuel_c2))**(1.0 / c(ifuel_c2)))
            isf_c2 = log(work) / (-b(ifuel_c2))           ! /* 41(a,b) - 2009 */
            work = max(0.01, 1.0 - (rsf_d1 / a(ifuel_d1))**(1.0 / c(ifuel_d1)))
            isf_d1 = log(work) / (-b(ifuel_d1))           ! /* 41(a,b) - 2009 */
              
            isf = pc * 1.0e-2 * isf_c2 + (1.0 - pc * 1.0e-2) * isf_d1 
                                                          ! /* 42a - 2009 */
         else if(ifuel == 12 .or. ifuel == 13) then
         ! case ("M3", "M4")
            cent = 100.0
            rsz = ROScalc(ifuel, isz, NoBUI, fmc, sfc, pc, cent, cured, cbh)
            rsf = rsz * sf                                ! /* 40 */ 
            rsz = ROScalc(ifuel_d1, isz, NoBUI, fmc, sfc, pc, pdf, cured, cbh)
            rsf_d1 = rsz * sf                             ! /* 40 */
              
            work = max(0.01, 1.0 - (rsf / a(ifuel))** (1.0 / c(ifuel)))
            isf = log(work) / (-b(ifuel))                 ! /* 41(a,b) - 2009 */
            work = max(0.01, 1.0 - (rsf_d1 / a(ifuel_d1))**(1.0 / c(ifuel_d1)))
            isf_d1 = log(work) / (-b(ifuel_d1))           ! /* 41(a,b) - 2009 */
              
            isf = pdf * 1.0e-2 * isf + (1.0 - pdf * 1.0e-2) * isf_d1 
                                                          ! /* 42a - 2009 */
         else if(ifuel == 17 .or. ifuel == 18) then
         !  case ("O1a", "O1b")
            if (cured < 58.8) then
               cf = 0.005 * (exp(0.061 * cured) - 1.0)    ! /* 35a - 2009 */
            else
               cf = 0.176 + 0.02 * (cured - 58.8)         ! /* 35b - 2009 */
            end if

            work = max(0.01, 1.0 - (rsf / (cf * a(ifuel)))**(1.0 / c(ifuel)))
            isf = log(work) / (-b(ifuel))                 ! /* 43(a,b) - 2009 */
         end if
         
!//       wse = log(isf/(0.208 * ff))/0.05039                       ! /* 44 */

         wse = 1.0 / 0.05039 * log(isf / (0.208 * ff))  ! /* 44a , 44d- 2009 */
         if (wse > 40.0) then                           ! /* 44e - 2009 */
            if (isf < (0.999 * 2.496 * ff)) then
               wse = 28.0 - (1.0 / 0.0818 * log(1.0 - isf / (2.496 * ff)))
                                                        ! /* 44b - 2009 */
            else
               wse = 112.45                             ! /* 44c - 2009 */
            end if
         end if
         
         wsx = ws * sin(waz) + wse * sin(saz)           ! /* 47 */
         wsy = ws * cos(waz) + wse * cos(saz)           ! /* 48 */
         
         wsv = sqrt(wsx * wsx + wsy * wsy)              ! /* 49 */
         raz = acos(wsy / wsv) !   /* in radians */       /* 50 */
         if (wsx < 0.0) raz = 2.0 * pi - raz            ! /* 51 */
!
      else
         wsv = ws
         raz = waz
      end if
!
!****************************************************************************** 
   
      isi = ISIcalc(ff, wsv)
!
!******************************************************************************
! Calculate CFB and ROS
      if (trim(fueltype) == "C6") then
         ! /* We use C6calc to calculate CFB */
         call C6calc(ifuel, isi, bui, fmc, sfc, cbh, ros, cfb)
      else
         ros = ROScalc(ifuel, isi, bui, fmc, sfc, pc, pdf, cured, cbh)
         if (cfl > 0.0 .and. sfc > 0.0) then
            cfb = CFBcalc(fmc, sfc, ros, cbh)
         else
            cfb = 0.0
         end if
      end if
!
!******************************************************************************
!** LBcalc
      if (fueltype(1:2) == "O1") then
         if (wsv >= 1.0) then
            lb = 1.1 * wsv**0.464 !/* corrected from "+" to "*" in the errata; 80 */
         else
            lb = 1.0                                               ! /* 81 */
         end if
      else
         lb = 1.0 + 8.729 * (1.0 - exp(-0.030 * wsv))**2.155       ! /* 79 */
      end if
!
!******************************************************************************
!** BROScalc      
      bfw  = exp(-0.05039 * wsv)                                   ! /* 75 */
      bisi = 0.208 * bfw * ff                                      ! /* 76 */
!/* Note the BUI effect is captured in ROScalc */
      bros = ROScalc(ifuel, bisi, bui, fmc, sfc, pc, pdf, cured, cbh) ! /* 77 */
      fros = (ros + bros) * 0.5 / lb                               ! /* 89 */
!/* TROS is the rate of spread towards angle theta */
      ecc  = sqrt(1.0 - 1.0 / (lb * lb))                !      /* eccentricity */
       ! /* note: this is the old method using the focus as the ignition point */
      tros = ros * (1.0 - ecc) / (1.0 - ecc * cos(theta - raz))

!******************************************************************************
! Calculate TFC
      cfc = cfl * cfb
      select case (trim(fueltype))
         case ("M1", "M2")
            cfc = pc * 1.0e-2 * cfc
         case ("M3", "M4")
            cfc = pdf * 1.0e-2 * cfc
      end select
      tfc = sfc + cfc

!******************************************************************************
! Calculate HFI
      hfi = 300.0 * tfc * ros
      
!*** LBtcalc      
!*      if (accel) then
!*         lbt = lb
!*      else if (t > 0.0) then
!*         select case (trim(fueltype))
!*           case ("C1", "O1a", "O1b", "S1", "S2", "S3")
!*              alpha = 0.115                                  !  /* page 41 */
!*           case default
!*              alpha = 0.115 - 18.8 * (cfb**2.5) * exp(-8.0 * cfb)  !  /* 72 */
!*         end select
!*         lbt = (lb - 1.0) * (1.0 - exp(-alpha * t)) + 1.0    ! /* 81 - 2009 */
!*      else
!*         lbt = 1.0      
!*      end if

      fbp%sfc  = sfc
      
      fbp%cfb  = cfb
      fbp%ros  = ros
      fbp%fros = fros
      fbp%bros = bros
      fbp%tfc  = tfc
      fbp%hfi  = hfi
      
      return
   end subroutine cffeps_fbpcalc
   
!****************************************************************************** 
!/* Foliar Moisture Content (FMC) calculation - FMCcalc
!   Note that 0.5 is added before the integer conversion in equations 2 and 4
!   Note that equations 1 and 3 use positive longitude values for Canada */
!
!    int    Dj,              /* Julian Day */
!           D0,              /* Julian day of minimum FMC */
!    double ELEV,            /* Elevation [m ASL] */
!           LAT,             /* Latitude [decimal degrees] */
!           LON,             /* Longitude [decimal degrees] */
   real(kind=4) function FMCcalc(dj, lat, lon, elev)
      implicit none
      integer(kind=4), intent(in) :: dj
      real(kind=4),    intent(in) :: lat, lon, elev
            
      integer(kind=4) :: d0, nd
      real(kind=4)    :: latn, lonn, fmc
!
      fmc = -1.0
         
!/* Make LON positive for the Western Hemisphere */
      lonn = -lon

!/* Calculate d0, date of min FMC (assuming it is not provided) */
      if (elev <= 0.0) then
         latn = 46.0 + 23.4 * exp(-0.0360 *(150. + lonn))      ! /* 1 */
         d0 = int(151.0 * lat / latn + 0.5)                   ! /* 2 (+0.5) */
      else
         latn = 43.0 + 33.7 * exp(-0.0351 * (150.0 + lonn))    ! /* 3 */
         d0 = int(142.1 * lat / latn + (0.0172 * elev) + 0.5) ! /* 4 (+0.5) */
      end if

      nd = abs(dj - d0)                                          ! /* 5 */

      if (nd < 30) then
         fmc = 85.0 + 0.0189 * nd * nd                           ! /* 6 */
      else if (nd < 50) then
         fmc = 32.9 + 3.17 * nd - 0.0288 * nd * nd               ! /* 7 */
      else
         fmc = 120.0                                             ! /* 8 */
      end if
!
      FMCcalc = fmc
      
      return
   end function FMCcalc
!****************************************************************************** 

!****************************************************************************** 
!/* Surface Fuel Consumption (SFC) calculation */ - SFCcalc
   real(kind=4) function SFCcalc(ifuel, bui, ffmc, pc, gfl, sfl, ffl, wfl)
      implicit none
            
      integer(kind=4), intent(in) :: ifuel
      real(kind=4),   intent(in)  :: bui, ffmc, pc, gfl, sfl, ffl, wfl

      real(kind=4) :: ffc, wfc, sfc

      sfc = -1.0
      if (ifuel == 1) then
!       case ("C1")
!//       sfc = 1.5 * (1. - exp(-0.230 * (ffmc - 81.0)))     ! /* 9 */
         if (ffmc > 84.0) then                              ! /* 9a - 2009 */
            sfc = sfl * 0.5 + sfl * 0.5 * sqrt(1.0 - exp(-0.23 * (ffmc - 84.0)))
         else                                               ! /* 9b - 2009 */
            sfc = sfl * 0.5 - sfl * 0.5 * sqrt(1.0 - exp(-0.23 * (84.0 - ffmc)))
         end if

      else if (ifuel == 2 .or. ifuel == 12 .or. ifuel == 13) then
!      case ("C2", "M3", "M4")
         sfc = sfl * (1.0 - exp(-0.0115 * bui))             ! /* 10 */

      else if (ifuel == 3 .or. ifuel == 4) then
!      case ("C3", "C4")
         sfc = sfl * (1.0 - exp(-0.0164 * bui))**2.24       ! /* 11 */

      else if (ifuel == 5 .or. ifuel == 6) then
!      case ("C5", "C6")
         sfc = sfl * (1.0 - exp(-0.0149 * bui))**2.48       ! /* 12 */

      else if (ifuel == 7) then
!      case ("C7")
         if (ffmc > 70.0) then
            ffc = ffl * (1.0 - exp(-0.104 * (ffmc - 70.0))) ! /* 13 */
         else
            ffc = 0.0
         end if
         wfc = wfl * (1.0 - exp(-0.0201 * bui))             ! /* 14 */
         sfc = ffc + wfc                                    ! /* 15 */

      else if (ifuel == 8) then
!       case ("D1")
         sfc = sfl * (1.0 - exp(-0.0183 * bui))             ! /* 16 */

      else if (ifuel == 10 .or. ifuel == 11) then
!       case ("M1", "M2")
         ffc = pc * 1.0e-2 * (SFLs(2) * (1.0 - exp(-0.0115 * bui)))
         wfc = (100.0 - pc) * 1.0e-2 * (SFLs(8) * (1.0 - exp(-0.0183 * bui)))
         sfc = ffc + wfc                                    ! /* 17 */

      else if (ifuel == 17 .or. ifuel == 18) then
!       case ("O1a", "O1b")
         sfc = gfl                                          ! /* 18 */
         if (gfl < 0.0) sfc = 0.35

      else if (ifuel == 14) then
!       case ("S1")
         ffc = ffl * (1.0 - exp(-0.025 * bui))              ! /* 19 */
         wfc = wfl * (1.0 - exp(-0.034 * bui))              ! /* 20 */
         sfc = ffc + wfc                                    ! /* 25 */

      else if (ifuel == 15) then
!       case ("S2")
         ffc = ffl * (1.0 - exp(-0.013 * bui))             ! /* 19 */
         wfc = wfl * (1.0 - exp(-0.060 * bui))              ! /* 20 */
         sfc = ffc + wfc                                    ! /* 25 */

      else if (ifuel == 16) then
!       case ("S3")
         ffc = ffl * (1.0 - exp(-0.0166 * bui))            ! /* 19 */
         wfc = wfl * (1.0 - exp(-0.0210 * bui))            ! /* 20 */
         sfc = ffc + wfc                                    ! /* 25 */
      else
         sfc = -1.0
      end if

      SFCcalc = sfc
      return
   end function SFCcalc
!****************************************************************************** 
!
!****************************************************************************** 
!/* Forest Floor Fuel Consumption (de Groot et al 2009) */
!/* this is a test of the FFFC as a substitute for SFC */
   real(kind=4) function FFFCcalc(ifuel, depth, bd, load, ffmc, dmc, dc, bui)
      implicit none
            
      integer(kind=4), intent(in) :: ifuel
      real(kind=4),   intent(in)  :: depth, bd, load, ffmc, dmc, dc, bui

      real(kind=4) :: fffc

      ! /* The general equations: n = 128 */
      ! /* r^2 */
      if (dc > 0.0 .and. depth > 0.0 .and. load > 0.0) then
         fffc = 1.1866 * exp(-7.388 + 0.754 * log(dc) + 0.691 * log(depth) + &
                0.608*log(bd))                           !         /* 0.796 */
      else if (dc > 0.0 .and. load > 0.0) then
         fffc = 1.1852 * exp(-4.252 + 0.710 * log(dc) + 0.671 * log(load))
                                                         !         /* 0.795 */
      else if (load > 0.0) then
         fffc = 1.2384 * exp(-0.882 + 1.082 * log(load)) !         /* 0.741 */
      else if (dc > 0.0) then
         fffc = 1.2821 * exp(-7.672 + 1.478 * log(dc))   !         /* 0.699 */
      end if


      !/* These are fuel specific equations.  Though some perform poorer than
      !   the general equations, they are used when applicable */
      if (ifuel == 1) then       ! /* n = 6 */
         if (ffmc > 0.0) fffc = -6.142 + 0.083 * ffmc    !        /* 0.856 */
      else if (ifuel == 2) then  ! /* n = 30 */
         if (load > 0.0 .and. dmc > 0.0) &
            fffc = 0.721 + 0.187 * load + 0.030 * dmc    !        /* 0.206 */
      else if (ifuel == 3 .or. ifuel == 4) then
         if (depth > 0.0 .and. bd > 0.0 .and. dc > 0.0) then
            fffc = -0.965 + 0.181 * depth + 0.012 * bd + 0.003 * dc
                                                         !        /* 0.429 */
         else if (dc > 0.0 .and. load > 0.0) then
            fffc = 1.0959 * exp(-3.486 + 0.612 * log(dc) + 0.484 * log(load))
                                                         !        /* 0.639 */
         end if
      
      else if (ifuel == 8) then
         if (dc > 0.0 .and. bui > 0.0) &
            fffc = 0.924 + 0.023 * dc - 0.082 * bui      !        /* 0.906 */
      end if

      ! Just use the default equation
      if (ifuel > 18 .or. fffc < 0.0) then
         if (dc > 0.0 .and. load > 0.0) &
            fffc = 1.1852 * exp(-4.252 + 0.710 * log(dc) + 0.671 * log(load))
                                                         !        /* 0.795 */
      end if

      FFFCcalc = max(fffc, 0.0)
      return
   end function FFFCcalc
!****************************************************************************** 
!
!****************************************************************************** 
   real(kind=4) function ISIcalc(ff, wsv)
      implicit none
            
      real(kind=4), intent(in)  :: ff, wsv

      real(kind=4) :: fw

         if (wsv < 40.0) then
            fw = exp(0.05039 * wsv)                                ! /* 53 */
         else
            fw = 12.0 * (1.0 - exp(-0.0818 * (wsv - 28.0)))        ! /* 53a */
         end if
         
         isicalc = 0.208 * fw * ff                                 ! /* 52 */
      
      return
   end function ISIcalc
!****************************************************************************** 

!****************************************************************************** 
   subroutine C6calc(ifuel, isi, bui, fmc, sfc, cbh, ros, cfb)
      implicit none
            
      integer(kind=4), intent(in)  :: ifuel
      real(kind=4),    intent(in)  :: isi, bui, fmc, sfc, cbh
      real(kind=4),    intent(out) :: ros, cfb

      real(kind=4) :: t, h, fme, rsc, rsi, rss, fme_avg, be

      fme_avg = 0.778                                ! /* page 37 */
      t = 1500.0 - 2.75 * fmc                        ! /* 59 */
      h = 460.0 + 25.9 * fmc                         ! /* 60 */
      fme = (1.5 - 0.00275 * fmc)**4 / (460.0 + 25.9 * fmc) * 1000.0 ! /* 61 */
      rsi = 30.0 * (1.0 - exp(-0.08 * isi))**3                       ! /* 62 */
      !BEcalc
      if (bui > 0.0 .and. BUIo(ifuel) > 0.0) then
         be = exp(50.0 * log(q(ifuel)) * (1.0 / bui - 1.0 / BUIo(ifuel)))
                                                                     ! /* 54 */
      else
         be = 1.0
      end if
      rss = rsi * be                                                 ! /* 63 */
      rsc = 60.0 * (1.0 - exp(-0.0497 * isi)) * fme / fme_avg        ! /* 64 */

      if (rsc > rss .and. sfc > 0.0) then
         cfb = CFBcalc(fmc, sfc, rss, cbh)
         ros = rss + (cfb) * (rsc - rss)                             ! /* 65 */
      else
         cfb = 0.0
         ros = rss
      end if
      
      return
   end subroutine C6calc
   
!****************************************************************************** 
!/* Crown Fraction Burned (CFB) calculation */
   real(kind=4) function CFBcalc(fmc, sfc, ros, cbh)
      implicit none
            
      real(kind=4),  intent(in)  :: fmc, sfc, ros, cbh
      
      real(kind=4) :: csi, rso
      
      CFBcalc = 0.0
      csi = 0.001 * (cbh**1.5) * (460.0 + 25.9 * fmc)**1.5    ! /* 56 */
      rso = csi / (300.0 * sfc)                               ! /* 57 */
      if (ros > rso) CFBcalc = 1.0 - exp(-0.23 * (ros - rso)) ! /* 58 */

      return
   end function CFBcalc
   
!****************************************************************************** 
!/* Rate of Spread calculations */
!
   recursive function ROScalc(ifuel, isi, bui, fmc, sfc, pc, pdf, cured, &
                              cbh) result(ros)
      implicit none
            
      integer(kind=4), intent(in)  :: ifuel
      real(kind=4),  intent(in)    :: isi, bui, fmc, sfc, pc, pdf, cured, cbh
      
      real(kind=4)    :: rsi, rsi_m34, cf, be, ros, cfb
      integer(kind=4) :: ifuel_c2, ifuel_d1
      
      rsi = -1.0
      
      ifuel_c2 = 2
      ifuel_d1 = 8
!/*Note that only preliminary RSS calculations are done for C6 in this routine*/
      if (ifuel <= 8 .or. ifuel == 14 .or. ifuel == 15 .or. ifuel == 16) then
        ! case ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "D1", "S1", "S2", "S3")
         rsi = a(ifuel) * (1.0 - exp(-b(ifuel) * isi))**c(ifuel)  !  /* 26 */
      else if(ifuel == 10) then
        ! case ("M1")
         rsi =       pc * 1.0e-2 * ROScalc(ifuel_c2, isi, NoBUI, fmc, sfc, &
                                           pc, pdf, cured, cbh) +          &
               (100.0 - pc) * 1.0e-2 * ROScalc(ifuel_d1, isi, NoBUI, fmc, sfc, &
                                               pc, pdf, cured, cbh)     ! /* 27 */
      else if(ifuel == 11) then
        ! case ("M2")
         rsi =       pc * 1.0e-2 * ROScalc(ifuel_c2, isi, NoBUI, fmc, sfc, &
                                           pc, pdf, cured, cbh) +          &
                0.2 * (100.0 - pc) * 1.0e-2 * ROScalc(ifuel_d1, isi, NoBUI, fmc,&
                                              sfc, pc, pdf, cured, cbh) ! /* 28 */
      else if(ifuel == 12) then
        ! case ("M3")
         rsi_m34 = a(ifuel) * (1.0 - exp(-b(ifuel) * isi))**c(ifuel) 
                                                           ! /* 30 - 2009 */
         rsi = pdf * 1.0e-2 * rsi_m34 + &
               (1.0 - pdf * 1.0e-2) * ROScalc(ifuel_d1, isi, NoBUI, fmc, sfc, &
                                              pc, pdf, cured, cbh) ! /* 29 - 2009 */
      else if(ifuel == 13) then
       ! case ("M4")
         rsi_m34 = a(ifuel) * (1.0 - exp(-b(ifuel) * isi))**c(ifuel) 
                                                           ! /* 32 - 2009 */
         rsi = pdf * 1.0e-2 * rsi_m34 + &
               0.2 * (1.0 - pdf * 1.0e-2) * ROScalc(ifuel_d1, isi, NoBUI, fmc,&
                                            sfc, pc, pdf, cured, cbh) ! /* 29 - 2009 */
      else if(ifuel == 17 .or. ifuel == 18) then
       ! case ("O1a", "O1b")
         if (cured < 58.8) then
            cf = 0.005 * (exp(0.061 * cured) - 1.0)            ! /* 35a - 2009 */
         else
            cf = 0.176 + 0.02 * (cured - 58.8)                 ! /* 35b - 2009 */
         end if
!/* RSI has been substituted for ROS in eqn 36 */
         rsi = a(ifuel) * ((1.0 - exp(-b(ifuel) * isi))**c(ifuel)) * cf 
                                                            ! /* 36 */
      end if
      
      if (ifuel == 6) then
       ! case ("C6")
         call C6calc(ifuel, isi, bui, fmc, sfc, cbh, ros, cfb)
                                !/* included here for completeness */
      else
         !BEcalc
         if (bui > 0.0 .and. BUIo(ifuel) > 0.0) then
            be = exp(50.0 * log(q(ifuel)) * (1.0 / bui - 1.0 / BUIo(ifuel)))
                                                                  ! /* 54 */
         else
            be = 1.0
         end if
         ros = be * rsi
      end if
      
      ros = max(0.000001, ros)
      
      return
   end function ROScalc
   
end module mach_cffeps_fbp_mod