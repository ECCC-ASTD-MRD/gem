!-------------------------------------- LICENCE BEGIN --------------------------
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

module linoz
   implicit none
   private
   public :: linoz3

contains

   !/@*
   subroutine linoz3(dbus, vbus, fbus, dt, kount, ni, nkm1, nk)
      use iso_c_binding
      use debug_mod, only: init2nan
      use tdpack_const, only: CAPPA, CONSOL2, GRAV, PI, STEFAN, RGASD
      use phy_options
      use phybus
      use linoz_param
      use tendency, only: apply_tendencies
      implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      !@object Stratospheric Ozone chemistry
      ! Produces the ozone tendency due to stratospheric 
      ! photochemistry (based on LINOZ scheme)
      !@arguments
      ! n        horizonal index
      ! nkm1     vertical  index nk-1
      ! nk       vertical  index
      ! dt       timestep
      ! trnch    index of vertical slice n*nk
      ! kount    number of timesteps
      ! vbus     volatile bus
      ! fbus     permanent bus
      ! dbus     dynamics bus

      integer, intent(in) :: ni, nkm1, nk, kount
      real,    intent(in) :: dt
      real, dimension(:), pointer, contiguous :: dbus, fbus, vbus

      !@author J. de Grandpre (ARQI): February 2013
      !@revisions
      ! * I. Ivanova (2015) - Complete remake
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
      real, dimension(ni, nkm1) :: o3new, ch4new, f11new, f12new, n2onew          !GMV3

      integer :: i, k, ni2
      real :: delage, delrlx

#define PHYPTRDCL
#include "linoz_ptr.hf"

      !----------------------------------------------------------------
      if (.not.llinoz) return

      call msg_toall(MSG_DEBUG, 'linoz [BEGIN]')
      if (timings_L) call timing_start_omp(416, 'linoz', 46)

#undef PHYPTRDCL
#include "linoz_ptr.hf"

      call init2nan(o3new, ch4new, n2onew, f11new, f12new)   !GMV3
      
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

      ni2 = int(float(ni)/2.)

      IF_KOUNT0: if (kount == 0) then

         if (.not.any(phyinread_list_S(1:phyinread_n) == 'o3ce')) zo3ce = -1.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'ttce')) zttce = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin4')) zlin4 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin5')) zlin5 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin6')) zlin6 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin7')) zlin7 = 0.

         ! Initialize ozone with climatology, if not found in analysis
         if (.not.any(dyninread_list_s == 'o3l') .or.   &
             .not.any(phyinread_list_s(1:phyinread_n) == 'tr/o3l:p')) then

            ! print*,'linoz: (0) Initialize O3L with climatology KOUNT0=', kount

            !    FK-HALOE kg /kg
            ! zo3lplus = zo3fk * 1E+9      !micro g /kg air <--      kg /kg air

            !   climato_phase2_v2 (Paul V)
            if (minval(zo3ce) >= 0.) zo3lplus = zo3ce * 1E+9  !micro g /kg air <-- kg /kg air
         endif

         IF_LINGHG1: if (llingh) then

            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin8' )) zlin8  = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin9' )) zlin9  = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin10')) zlin10 = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin11')) zlin11 = 0.

            ! Initialize GHG with climatology, if not found in analysis
            if (.not.any(dyninread_list_s == 'ch4l') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/ch4l:p')) zch4lplus = zch4 * 1E+9   !micro g /kg air <--      kg /kg air
            if (.not.any(dyninread_list_s == 'n2ol') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/n2ol:p')) zn2olplus = zn2o * 1E+9   !micro g /kg air <--      kg /kg air
            if (.not.any(dyninread_list_s == 'f11l') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/f11l:p')) zf11lplus = zcf11 * 1E+9  !micro g /kg air <--      kg /kg air
            if (.not.any(dyninread_list_s == 'f12l') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/f12l:p')) zf12lplus = zcf12 * 1E+9  !micro g /kg air <--      kg /kg air

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

            !        if(local_dbg .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
            !           write(lun_out,*) 'linoz: ', i,trnch,k, &
            !              'o3lplus, o3c_vmr=',                  &
            !               zo3lplus(i,k), o3c_vmr(i,k)   !,        &
            !              'pplus, shtj, tplus, huplus=',      &
            !               zpplus (i), shtj(i,k), ztplus(i,k), zhuplus(i,k)

            ptop(i,k) = shtj(i,k)  *zpplus(i)       ! pressure (Pa) at the upper interface
            pbot(i,k) = shtj(i,k+1)*zpplus(i)       ! pressure (Pa) at the bottom interface

            if (age_linoz) then
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
            end if 

            ! Set lower limit on species as in BIRA (units mole /mole)
            o3_vmr(i,k) = zo3lplus(i,k) *1E-9 * mwt_air/ mwt_o3      !mole /mole vmr <-- micro g/kg air
            o3_vmr(i,k) = max(o3_vmr(i,k), QepsO3)

            ! FK-HALOE ozone climatology from cccmarad blended with ERA5 below 1 hPa
            o3c_vmr(i,k) = zo3fk(i,k) * mwt_air/ mwt_o3           !mole /mole vmr  <-- kg /kg
            if (zo3ce(i,k) >= 0.) then
               if (ptop(i,k) > p_linoz_meso) &
                  ! 'climato_phase2_v2' (Paul V)
                  o3c_vmr(i,k) = zo3ce(i,k) * mwt_air/ mwt_o3           !mole /mole vmr <-- kg /kg air
            end if 

            ! Ignore P-L term on RHS above p_linoz_c4=10hPa
            if (ptop(i,k) < p_linoz_c4 ) zlin4(i,k) = 0.
            ! Reduce P-L by 15%
            ! zlin4(i,k) = 0.85 * zlin4(i,k)   !reduce P-L by 15%

         enddo DO_K1
      enddo DO_I1   

      ! Total Column LINOZ ozone (D.U.)
      call linoz_xcol(o3_vmr, pbot, ptop, zo3col,ni,nkm1)
      zo3tc (:) =zo3col (:,nkm1)   !D.U.

      ! Total column ERA5-FK-HALOE ozone climatology (D.U.)
      call linoz_xcol(o3c_vmr, pbot, ptop, zo3ccol,ni,nkm1)
      zo3ctc(:) = zo3ccol(:,nkm1)  !D.U.

      ! GHG Linoz table of coefficients
      IF_LINGHG2: if (llingh) then

        ch4_vmr(:,1:nkm1) = zch4lplus(:,1:nkm1) *1E-9 * mwt_air/ mwt_ch4             !mole /mole vmr <-- micro g /kg
        n2o_vmr(:,1:nkm1) = zn2olplus(:,1:nkm1) *1E-9 * mwt_air/ mwt_n2o             !mole /mole vmr <-- micro g /kg
        f11_vmr(:,1:nkm1) = zf11lplus(:,1:nkm1) *1E-9 * mwt_air/ mwt_f11             !mole /mole vmr <-- micro g /kg
        f12_vmr(:,1:nkm1) = zf12lplus(:,1:nkm1) *1E-9 * mwt_air/ mwt_f12             !mole /mole vmr <-- micro g /kg

        ! Set lower limit on species as in BIRA (units mole /mole)
        ch4_vmr = max(ch4_vmr ,Qeps)
        n2o_vmr = max(n2o_vmr ,Qeps)
        f11_vmr = max(f11_vmr ,Qeps)
        f12_vmr = max(f12_vmr ,Qeps)

      end if IF_LINGHG2

      !    Linoz tendencies 

      ! print*,'linoz: (1) o3chmtd (:,nkm1)',nkm1,zo3chmtd(:,nkm1)
      ! print*,'linoz: (1) o3chmtd (:,nk)',nk,zo3chmtd(:,nk)

      call linoz_tend( &
           o3_vmr                                       , & !input, mole/mole vmr
           zo3col                                       , & !input, D.U.
           ztplus, zpplus, shtj, zhuplus                , & !input
           o3c_vmr,zttce,zo3ccol,zlin4,zlin5,zlin6,zlin7, & !input mole/mole vmr ERA-3 ozone climato in troposphere
           o3new                                  , & !output, mole /mole  !GMV3
           zo1chmtd, zo4chmtd, zo6chmtd, zo7chmtd , & !output, mole /mole /sec
           zo3chmtd                               , & !output, mole /mole /sec
           dt, ni, nkm1, nk)                          !input
      
      IF_LINGHG: if (llingh) then
         call linoz_tend_ghg( &
              ch4_vmr,n2o_vmr,f11_vmr,f12_vmr        , & !input, mole/mole vmr
              zpplus, shtj, zhuplus                  , & !input
              zlin8,zlin9,zlin10,zlin11                      , & !input
              ch4new ,n2onew, f11new ,f12new         , & !output, mole /mole  !GMV3
              zch4chmtd,zn2ochmtd,zf11chmtd,zf12chmtd, & !output, mole /mole /sec
              dt, ni, nkm1, nk)                          !input
      endif IF_LINGHG
      
      ! Ozone tendency micro g /kg air sec-1
      zo1chmtd=zo1chmtd*1E+9 * mwt_o3 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
      zo4chmtd=zo4chmtd*1E+9 * mwt_o3 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
      zo6chmtd=zo6chmtd*1E+9 * mwt_o3 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
      zo7chmtd=zo7chmtd*1E+9 * mwt_o3 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
      zo3chmtd=zo3chmtd*1E+9 * mwt_o3 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1

      ! print*,'linoz: (2) o3chmtd (:,nkm1)',nkm1,zo3chmtd(:,nkm1)
      ! print*,'linoz: (2) o3chmtd (:,nk)',nk,zo3chmtd(:,nk)

      ! Apply linoz tendencies micro g /kg air
      call apply_tendencies(zo3lplus,zo3chmtd, ztdmask, ni, nk, nkm1)

      ! Diagnostic level 1.5m: copy down the bottom hybrid level
      zo3lplus(:,nk)=zo3lplus(:,nkm1)

      !do i =1,ni
      !   do k = 1, nkm1    !1=TOA, nk-1=SFC

      !        if(local_dbg  .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
      !           write(lun_out,*) 'linoz: ', i,trnch,k, &
      !              ' o3vmr,  o3chmtd=',  &
      !                o3_vmr(i,k), zo3chmtd(i,k)

      !   end do ! k-loop
      !end do !i-loop

      ! GHG Prognostic LINOZ

      IF_LINGHG3: if (llingh) then

         do i = 1,ni
            do k = 1, nkm1    !1=TOA, nk-1=SFC
               ptop(i,k) = shtj(i,k)  *zpplus(i)       ! pressure (Pa) at the upper interface
               pbot(i,k) = shtj(i,k+1)*zpplus(i)       ! pressure (Pa) at the bottom interface
            end do !k-loop
         end do !i-loop

         if (out_linoz) then
         ! Total column prognostic ch4 (molec cm-2)
         call ghg_xcol(ch4_vmr, pbot, ptop, zch4col,ni,nkm1)
         zch4tc(:) =zch4col(:,nkm1)

         ! Total column prognostic n2o (molec cm-2)
         call ghg_xcol(n2o_vmr, pbot, ptop, zn2ocol,ni,nkm1)
         zn2otc(:) =zn2ocol(i,nkm1)

         ! Total column prognostic f11 (molec cm-2)
         call ghg_xcol(f11_vmr, pbot, ptop, zf11col,ni,nkm1)
         zf11tc(:) =zf11col(:,nkm1)

         ! Total column prognostic f12 (molec cm-2)
         call ghg_xcol(f12_vmr, pbot, ptop, zf12col,ni,nkm1)
         zf12tc(:) =zf12col(:,nkm1)
         end if 

         ! Linoz GHG tendencies micro g /kg air sec-1
         zch4chmtd = zch4chmtd *1E+9 * mwt_ch4 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
         zn2ochmtd = zn2ochmtd *1E+9 * mwt_n2o /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
         zf11chmtd = zf11chmtd *1E+9 * mwt_f11 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1
         zf12chmtd = zf11chmtd *1E+9 * mwt_f12 /mwt_air ! micro g /kg air sec-1 <-- mole /mole sec-1

         ! Apply linoz tendencies micro g /kg air
         call apply_tendencies(zch4lplus,zch4chmtd, ztdmask, ni, nk, nkm1)
         call apply_tendencies(zn2olplus,zn2ochmtd, ztdmask, ni, nk, nkm1)
         call apply_tendencies(zf11lplus,zf11chmtd, ztdmask, ni, nk, nkm1)
         call apply_tendencies(zf12lplus,zf12chmtd, ztdmask, ni, nk, nkm1)


      end if IF_LINGHG3


      if (timings_L) call timing_stop_omp(416)
      call msg_toall(MSG_DEBUG, 'linoz [END]')
      !----------------------------------------------------------------
      return
   end subroutine linoz3

end module linoz
