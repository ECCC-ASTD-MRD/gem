
module linoz
   implicit none
   private
   public :: linoz3

contains

   !/@*
   subroutine linoz3(pvars, dt, kount, ni, nk, nkm1)
      use iso_c_binding
      use debug_mod, only: init2nan
      use tdpack_const, only: CAPPA, CONSOL2, GRAV, PI, STEFAN, RGASD
      use phy_options
      use phybusidx
      use phymem, only: phyvar
      use linoz_param
      use linoz_tend_mod, only: linoz_tend
      use tendency, only: apply_tendencies
      implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      !@object Stratospheric Ozone chemistry
      ! Produces the ozone tendency due to stratospheric 
      ! photochemistry (based on LINOZ scheme)
      !@arguments
      ! pvars    list of all phy vars (meta + slab data)
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! nkm1     vertical dimension nk-1
      ! dt       timestep
      ! kount    number of timesteps

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, nkm1, nk, kount
      real,    intent(in) :: dt

      !@author J. de Grandpre (ARQI): February 2013
      !@revisions
      ! * I. Ivanova (2015) - Major revision
      !*@/
#include <rmn/msg.h>

      include "phyinput.inc"
      include "ccc_tracegases.cdk"

      external ccc_tracedata

      ! declaration of local arguments
      !
      !          - input -
      ! o3_vmr   linoz tracer (t+)        3D  o3    vmr (mole /mole)
      ! ch4_vmr  linoz tracer (t+)        3D  ch4   vmr (mole /mole)
      ! n2o_vmr  linoz tracer (t+)        3D  n2o   vmr (mole /mole)
      ! f11_vmr  linoz tracer (t+)        3D  cfc11 vmr (mole /mole)
      ! f12_vmr  linoz tracer (t+)        3D  cfc12 vmr (mole /mole)
      ! o3c_vmr  climatology FK-HALOE blended with ERA5 zonavg  2D  o3    vmr (mole /mole)
      ! ztplus   temperature (K) mid-thermo layer
      ! zhuplus  specific humidity (kg /kg)
      ! zhuplus  surface pressure
      ! shtj     sigma levels on momentum/flux levels (nk)
      !
      !          - input/output -
      ! zo3lplus  linoz tracer (t*)        3D  o3    mmr (micro g /kg air)
      ! ziarplus age of air                3D  air   seconds

      real, dimension(ni, nkm1) :: ptop, pbot       ! nk-1 layers
      real, dimension(ni, nk)   :: shtj             ! nk "flux" levels

      real, dimension(ni, nkm1) :: o3_vmr, o3c_vmr, ch4_vmr, n2o_vmr, f11_vmr, f12_vmr ! nk-1 layers
      real, dimension(ni, nkm1) :: hu_vmr, o8_tend
      real, dimension(ni, nkm1) :: o3new, ch4new, f11new, f12new, n2onew          !GMV3

      integer :: i, k
      real :: delage, delrlx

      real, parameter :: T195 = 195. !DegK
      
      real, parameter :: radlat55 = 55. *pi/180. !Radians <- 55 Degrees latitude 

#define PHYPTRDCL
#include "linoz_ptr.hf"

      !----------------------------------------------------------------
      if (.not.llinoz) return

      call msg_toall(MSG_DEBUG, 'linoz [BEGIN]')
      if (timings_L) call timing_start_omp(416, 'linoz', 46)

#undef PHYPTRDCL
#include "linoz_ptr.hf"

      call init2nan(o3_vmr, o3c_vmr, ch4_vmr, n2o_vmr, f11_vmr, f12_vmr)
      call init2nan(hu_vmr, o8_tend) 
      call init2nan(o3new, ch4new, n2onew, f11new, f12new)  

      ! Calculate age of air
      !
      ! ---  change in age of air (age in units of years)
      !   --- assuming points with a water vapour mixing ratio
      !         greater than 20 ppmv are in the troposphere (ozone smaller than 150 ppbv)
      !   --- tropospheric points are relaxed back to zero age (here
      !         defined with an offset of 2 years) with a time constant
      !         of 0.1 days
      delage = dt/(365.0*86400.0)
      delrlx = exp(-1.0*dt/8640.0)

      IF_KOUNT0: if (kount == 0) then

         if (.not.ISPHYIN('o3ce')) zo3ce = -1.
         if (.not.ISPHYIN('ttce')) zttce = 0.
         if (.not.ISPHYIN('lin4')) zlin4 = 0.
         if (.not.ISPHYIN('lin5')) zlin5 = 0.
         if (.not.ISPHYIN('lin6')) zlin6 = 0.
         if (.not.ISPHYIN('lin7')) zlin7 = 0.

         ! Initialize ozone with climatology, if not found in analysis
         ! micro g /kg air <-- kg /kg air
         !   climato_phase2_v2 (Paul V)
         if (.not.(ISDYNIN('o3l') .or. ISPHYIN('tr/o3l:p'))) then
            if (minval(zo3ce) >= 0.) zo3lplus = zo3ce * 1E+9
         endif

         IF_LINGHG1: if (llingh) then

            if (.not.ISPHYIN('lin8' )) zlin8  = 0.
            if (.not.ISPHYIN('lin9' )) zlin9  = 0.
            if (.not.ISPHYIN('lin10')) zlin10 = 0.
            if (.not.ISPHYIN('lin11')) zlin11 = 0.

            ! Initialize GHG with climatology, if not found in analysis
            ! micro g /kg air <--      kg /kg air
            if (.not.(ISDYNIN('ch4l') .or. ISPHYIN('tr/ch4l:p'))) zch4lplus = zch4 * 1E+9
            if (.not.(ISDYNIN('n2ol') .or. ISPHYIN('tr/n2ol:p'))) zn2olplus = zn2o * 1E+9
            if (.not.(ISDYNIN('f11l') .or. ISPHYIN('tr/f11l:p'))) zf11lplus = zcf11 * 1E+9
            if (.not.(ISDYNIN('f12l') .or. ISPHYIN('tr/f12l:p'))) zf12lplus = zcf12 * 1E+9

         end if IF_LINGHG1

      endif IF_KOUNT0


      ! Calculate shtj, sigma at flux levels
      
      do i = 1, ni
         shtj(i,1) = zsigw(i,1) * sqrt(zsigw(i,1) / zsigw(i,2))
         shtj(i,nk) = 1.0
      enddo
      do k = 2, nkm1
         do i = 1, ni
            shtj(i,k) = sqrt(zsigw(i,k-1) * zsigw(i,k))
         enddo
      enddo

      ! --------------------
      !  Loop on longitudes
      ! -------------------
      DO_I1: do i = 1, ni

         ! ----------------------------
         ! Loop on all vertical levels  
         ! ----------------------------

         DO_K1: do k = 1, nkm1 !nk-1 layers, 1=TOA-layer, nk=SFC-layer

            ptop(i,k) = shtj(i,k)  *zpplus(i)       ! pressure (Pa) at the upper interface
            pbot(i,k) = shtj(i,k+1)*zpplus(i)       ! pressure (Pa) at the bottom interface

            IF_AGE: if (age_linoz) then
               ! Calculate age of air
               !
               ! ---  change in age of air (age in units of years)
               !   --- assuming points with a water vapour mixing ratio
               !         greater than 20 ppmv are in the troposphere (ozone smaller than 150 ppbv)
               !   --- tropospheric points are relaxed back to zero age (here
               !         defined with an offset of 2 years) with a time constant
               !         of 0.1 days

               ! Troposphere
               !       if (zo3lplus (i,k).lt.150.0E-9) then  ! ppb <-- mole /mole vmr
               ! PBL pressure 900 hPa 
               if (ptop(i,k) >= 0.900E+05) then !900 hPa 
                  zairplus(i,k) = 2.0 + (zairplus(i,k) - 2.0)*delrlx   !years
               else
                  ! Stratosphere
                  zairplus(i,k) = zairplus(i,k) + delage               !years
               endif
            endif IF_AGE

            ! Set lower limit on species as in BIRA (units mole /mole)
            !mole /mole vmr <-- micro g/kg air
            o3_vmr(i,k) = zo3lplus(i,k) *1E-9 * mwt_air/ mwt_o3
            o3_vmr(i,k) = max(o3_vmr(i,k), QepsO3)

            ! FK-HALOE ozone climatology from cccmarad blended with ERA5 below 1 hPa
            o3c_vmr(i,k) = zo3fk(i,k) * mwt_air/ mwt_o3
            ! 'climato_phase2_v2' (Paul V)
            if (zo3ce(i,k) >= 0. .and. ptop(i,k) > p_linoz_meso) &
                 o3c_vmr(i,k) = zo3ce(i,k) * mwt_air/ mwt_o3

            ! Ignore P-L term on RHS above p_linoz_c4=10hPa
            if (ptop(i,k) < p_linoz_c4) zlin4(i,k) = 0.

         enddo DO_K1
      enddo DO_I1

      IF_HET: if (linoz_het) then
         !           Hetchem begin : if hetchem...

         ! Computation of ozone tendency 
         ! Heterogeneous chemistry in [55-90N] and [55-90S] latitude bands
         ! (tau_hetchm=10 days)

         zo9chmtd(:,nk) = 0.
         DO_K2: do k = 1, nkm1
            DO_I2: do i = 1, ni

               o8_tend(i,k) = 0.
               zo9chmtd(i,k) = 0.
               IF_ZLAT: if (.not.(abs(zdlat(i)) < radlat55)) then 

                  if (zcosas(i) > 1.E-3) then
                     ! With halogens/cold tracer between 22-12km
                     if (zgztherm(i,k) <= 22000. .and. zgztherm(i,k) >= 12000.) &
                          o8_tend(i,k) = -o3_vmr(i,k)*zohalplus(i,k)/tau_hetchm   ! mole /mole sec-1
                  endif

                  ! Computation of cold tracer tendency (OHAL)
                  ! N.H. Halogens activation [+55N, +90N] (tau_activ=4 hours)
                  if (ztplus(i,k) <= T195 .and. zcosas(i) >  1.E-3) then  !195K threshold & sunlight
                     ! Activation between 22-12km
                     if (zgztherm(i,k) <= 22000. .and. zgztherm(i,k) >= 12000.) &
                          zo9chmtd(i,k) = (1. - zohalplus(i,k))/tau_activ     !sec-1
                  endif

               endif IF_ZLAT

               ! Halogen de-activation in both hemispheres (tau_deactn=5 days)
               if (zcosas(i) > 1.E-3) &
                    zo9chmtd(i,k) = zo9chmtd(i,k) - zohalplus(i,k)/tau_deactn         !sec-1
               
            enddo DO_I2
         enddo DO_K2
         
      endif IF_HET
      
      ! Total Column LINOZ ozone (D.U.)
      if (associated(zo3col)) then
         call linoz_xcol(o3_vmr, pbot, ptop, zo3col,ni,nkm1)
         if (associated(zo3tc)) zo3tc(:) = zo3col(:,nkm1)   !D.U.
      endif

      ! Tropospheric column LINOZ ozone (molec cm-2)
      if (associated(zo3xcol)) then
         hu_vmr(:,1:nkm1) = consth * zhuplus(:,1:nkm1)       !mole /mole vmr <-- kg/kg air ! Units ppmv
         call tropo_xcol(o3_vmr, hu_vmr, pbot, ptop, zo3xcol, ni, nkm1)
         if (associated(zo3xtc)) zo3xtc(:) = zo3xcol(:,nkm1)  ! (molec /cm2)
      endif

      ! Total column ERA5-FK-HALOE ozone climatology (D.U.)
      if (associated(zo3ccol)) then
         call linoz_xcol(o3c_vmr, pbot, ptop, zo3ccol,ni,nkm1)
         if (associated(zo3ctc)) zo3ctc(:) = zo3ccol(:,nkm1)  !D.U.
      endif

      ! GHG Linoz table of coefficients
      IF_LINGHG2: if (llingh) then

         !#TODO: see if we can avoid using all the _vmr 3d fields to avoid memory copy and reduce memory usage
         DO_K3: do k = 1, nkm1
            DO_I3: do i = 1, ni

               !mole /mole vmr <-- micro g /kg
               ch4_vmr(i,k) = zch4lplus(i,k) * 1E-9 * mwt_air / mwt_ch4
               n2o_vmr(i,k) = zn2olplus(i,k) * 1E-9 * mwt_air / mwt_n2o
               f11_vmr(i,k) = zf11lplus(i,k) * 1E-9 * mwt_air / mwt_f11
               f12_vmr(i,k) = zf12lplus(i,k) * 1E-9 * mwt_air / mwt_f12

               ! Set lower limit on species as in BIRA (units mole /mole)
               ch4_vmr(i,k) = max(ch4_vmr(i,k), Qeps)
               n2o_vmr(i,k) = max(n2o_vmr(i,k), Qeps)
               f11_vmr(i,k) = max(f11_vmr(i,k), Qeps)
               f12_vmr(i,k) = max(f12_vmr(i,k), Qeps)
               
            enddo DO_I3
         enddo DO_K3

      end if IF_LINGHG2

      !    Linoz tendencies 

      call linoz_tend( &
           o3_vmr                                       , & !input, mole/mole vmr
           zo3col                                       , & !input, D.U.
           ztplus, zpplus, shtj, zhuplus                , & !input
           o3c_vmr,zttce,zo3ccol,zlin4,zlin5,zlin6,zlin7, & !input mole/mole vmr ERA-3 ozone climato in troposphere
           o8_tend                                      , & !input mole/mole sec-1
           o3new                                  , & !output, mole /mole  !GMV3
           zo1chmtd, zo4chmtd, zo6chmtd, zo7chmtd , & !output, mole /mole /sec
           zo8chmtd                               , & !output, mole /mole /sec
           zo3chmtd                               , & !output, mole /mole /sec
           dt, ni, nkm1, nk)                          !input
      
      IF_LINGHG: if (llingh) then
         call linoz_tend_ghg( &
              ch4_vmr,n2o_vmr,f11_vmr,f12_vmr        , & !input, mole/mole vmr
              zpplus, shtj, zhuplus                  , & !input
              zlin8,zlin9,zlin10,zlin11              , & !input
              ch4new ,n2onew, f11new ,f12new         , & !output, mole /mole  !GMV3
              zch4chmtd,zn2ochmtd,zf11chmtd,zf12chmtd, & !output, mole /mole /sec
              dt, ni, nkm1, nk)                          !input
      endif IF_LINGHG
      
      ! Ozone tendency micro g /kg air sec-1
      ! micro g /kg air sec-1 <-- mole /mole sec-1
      zo1chmtd = zo1chmtd*1E+9 * mwt_o3 /mwt_air
      zo4chmtd = zo4chmtd*1E+9 * mwt_o3 /mwt_air
      zo6chmtd = zo6chmtd*1E+9 * mwt_o3 /mwt_air
      zo7chmtd = zo7chmtd*1E+9 * mwt_o3 /mwt_air
      if (linoz_het) zo8chmtd = zo8chmtd*1E+9 * mwt_o3 /mwt_air
      zo3chmtd = zo3chmtd*1E+9 * mwt_o3 /mwt_air

      ! Apply linoz tendencies micro g /kg air
      call apply_tendencies(zo3lplus,zo3chmtd, ztdmaskxdt, ni, nk, nkm1)

      ! Diagnostic level 1.5m: copy down the bottom hybrid level
      zo3lplus(:,nk) = zo3lplus(:,nkm1)
      
      ! Apply cold tracer tendencies
      if (linoz_het) &
           call apply_tendencies(zohalplus, zo9chmtd, ztdmaskxdt, ni, nk, nkm1)

      ! GHG Prognostic LINOZ

      IF_LINGHG3: if (llingh) then

         do i = 1,ni
            do k = 1, nkm1    !1=TOA, nk-1=SFC
               ptop(i,k) = shtj(i,k)  *zpplus(i)       ! pressure (Pa) at the upper interface
               pbot(i,k) = shtj(i,k+1)*zpplus(i)       ! pressure (Pa) at the bottom interface
            end do !k-loop
         end do !i-loop

         IF_OUT: if (out_linoz) then
            
            ! Total column prognostic ch4 (molec cm-2)
            if (associated(zch4col)) then
               call ghg_xcol(ch4_vmr, pbot, ptop, zch4col,ni,nkm1)
               if (associated(zch4tc)) zch4tc(:) = zch4col(:,nkm1)
            endif

            ! Total column prognostic n2o (molec cm-2)
            if (associated(zn2ocol)) then
               call ghg_xcol(n2o_vmr, pbot, ptop, zn2ocol,ni,nkm1)
               if (associated(zn2otc)) zn2otc(:) = zn2ocol(i,nkm1)
            endif

            ! Total column prognostic f11 (molec cm-2)
            if (associated(zf11col)) then
               call ghg_xcol(f11_vmr, pbot, ptop, zf11col,ni,nkm1)
               if (associated(zf11tc)) zf11tc(:) = zf11col(:,nkm1)
            endif

            ! Total column prognostic f12 (molec cm-2)
            if (associated(zf12col)) then
               call ghg_xcol(f12_vmr, pbot, ptop, zf12col,ni,nkm1)
               if (associated(zf12tc)) zf12tc(:) = zf12col(:,nkm1)
            endif
            
         endif IF_OUT

         ! Linoz GHG tendencies micro g /kg air sec-1
         ! micro g /kg air sec-1 <-- mole /mole sec-1
         zch4chmtd = zch4chmtd *1E+9 * mwt_ch4 /mwt_air 
         zn2ochmtd = zn2ochmtd *1E+9 * mwt_n2o /mwt_air
         zf11chmtd = zf11chmtd *1E+9 * mwt_f11 /mwt_air
         zf12chmtd = zf12chmtd *1E+9 * mwt_f12 /mwt_air

         ! Apply linoz tendencies micro g /kg air
         call apply_tendencies(zch4lplus, zn2olplus, zf11lplus, zf12lplus, &
              &                zch4chmtd, zn2ochmtd, zf11chmtd, zf12chmtd, &
              &                ztdmaskxdt, ni, nk, nkm1)

      end if IF_LINGHG3

      if (timings_L) call timing_stop_omp(416)
      call msg_toall(MSG_DEBUG, 'linoz [END]')
      !----------------------------------------------------------------
      return
   end subroutine linoz3

end module linoz
