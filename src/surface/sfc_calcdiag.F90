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
!-------------------------------------- LICENCE END ----------------------------

module sfc_calcdiag
   implicit none
   private
   public :: sfc_calcdiag3

contains

   !/@*
   subroutine sfc_calcdiag3(fbus,vbus,moyhr,acchr,dt,kount,step_driver,trnch,ni)
      use tdpack_const, only: CHLC, CHLF
      use sfc_options
      use sfcbus_mod
      use svs_configs
      
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object Calculates averages and accumulators of tendencies and diagnostics
      !@Arguments
      !          - Input/Output -
      ! fbus     permanent bus
      !          - input -
      ! vbus     volatile (output) bus
      !          - input -
      ! kount    timestep number
      ! dt       length of timestep
      ! trnch    slice number
      ! n        horizontal running length
      integer :: moyhr, acchr, kount, step_driver, trnch, ni
      real, dimension(:), pointer, contiguous :: fbus, vbus
      real :: dt
      !*@/
#include <rmn/msg.h>
      include "sfcinput.cdk"

      real, pointer, dimension(:), contiguous :: zaccevap, zdrain, zdrainaf, zfvapliqaf,&
           zisoil, zlatflaf, zleg, zlegaf, &
           zler, zleraf, zles, zlesaf, zletr, zletraf, zlev, zlevaf, zoverfl, &
           zoverflaf, zrofinlak, zrofinlakaf, zrootdp, zwflux, zwfluxaf, zwsoil, zinsmavg

      real, pointer, dimension(:,:), contiguous :: zlatflw, zrunofftot, &
           zrunofftotaf, zwatflow
      
      logical :: lacchr
      integer :: i, k, moyhr_steps
      real :: moyhri, tmp_r
      !-------------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'sfc_calcdiag [BEGIN]')
      if (timings_L) call timing_start_omp(480, 'sfc_calcdiag', 46)

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

      MKPTR1D(zaccevap, accevap)
      MKPTR1D(zdrain, drain)
      MKPTR1D(zdrainaf, drainaf)
      MKPTR1D(zfvapliqaf, fvapliqaf)
      MKPTR1D(zisoil, isoil)
      MKPTR1D(zlatflaf, latflaf)
      MKPTR1D(zleg, leg)
      MKPTR1D(zlegaf, legaf)
      MKPTR1D(zler, ler)
      MKPTR1D(zleraf, leraf)
      MKPTR1D(zles, les)
      MKPTR1D(zlesaf, lesaf)
      MKPTR1D(zletr, letr)      
      MKPTR1D(zletraf, letraf)
      MKPTR1D(zlev, lev)
      MKPTR1D(zlevaf, levaf)
      MKPTR1D(zoverfl, overfl)
      MKPTR1D(zoverflaf, overflaf)
      MKPTR1D(zrofinlak, rofinlak)
      MKPTR1D(zrofinlakaf, rofinlakaf)
      MKPTR1D(zrootdp, rootdp)
      MKPTR1D(zwflux, wflux)
      MKPTR1D(zwfluxaf, wfluxaf)
      MKPTR1D(zinsmavg, insmavg)

      MKPTR2D(zlatflw,latflw)
      MKPTR2D(zrunofftot,runofftot)
      MKPTR2D(zrunofftotaf,runofftotaf)
      MKPTR2D(zwatflow,watflow)

      
      lacchr = .false.
      if (acchr > 0) then
         lacchr = (mod(step_driver-1, acchr) == 0)
      elseif (acchr == 0) then
         lacchr = (step_driver-1 == 0)
      endif
      
      ! Common ISBA/SVS accumulators
      IF_ISBA_SVS: if (schmsol == 'ISBA' .or. schmsol == 'SVS') then

         IF_RESET: if (kount == 0 .or. lacchr) then

            ! Reset accumulators at t=T+00hr 
            
            ! Accum. of base drainage                      
            zdrainaf(:) = 0.
            
            ! Accumulation of surf. evaporation (HFLQ) (kg/m2 or mm)
            zfvapliqaf(:) = 0.
            
            !Accumulation of total surface runoff, per tile and total 
            zrunofftotaf(:,:) = 0.

         endif IF_RESET

         IF_KOUNT_NE_0: if (kount /= 0) then
            !       ACCUMULATE RUNNOFF FOR EACH SURFACE TYPE
            do k=1,nsurf+1
               do i=1,ni
                  zrunofftotaf(i,k) = zrunofftotaf(i,k) + zrunofftot(i,k)
               enddo
            enddo
         ENDIF IF_KOUNT_NE_0
      endif IF_ISBA_SVS


      ! ISBA only accumulators/calc.
      IF_ISBA: if (schmsol == 'ISBA') then

         nullify(zwsoil)
         if (wsoil > 0) zwsoil(1:ni) => fbus(wsoil+ni:)

         if (moyhr > 0) then
            if (kount == 0 .or. mod(step_driver-1, moyhr) == 0) then
               zinsmavg(:) = 0.
            endif
         endif

         IF_RESET_ISBA: if (kount == 0 .or. lacchr) then
            zlegaf(:) = 0.
            zleraf(:) = 0.
            zletraf(:) = 0.
            zlevaf(:) = 0.
            zlesaf(:) = 0.
            zoverflaf(:) = 0.
            zwfluxaf(:) = 0.
         endif IF_RESET_ISBA

         IF_KOUNT0: if (kount /= 0) then

            !#Take care of integrated soil moisture in kg/m**2
            IF_MOYHR: if (moyhr > 0) then
               do i = 1, ni
                  zinsmavg(i) = zinsmavg(i) + &
                       1000. * (zisoil(i) + zwsoil(i)) * zrootdp(i)
               end do

               if (mod(step_driver, moyhr) == 0) then
                  moyhr_steps = moyhr
                  if (step_driver == 0 .and. kount > 0) moyhr_steps = min(moyhr,kount)
                  moyhri = 1./float(moyhr_steps)
                  zinsmavg(:) = zinsmavg(:) * moyhri
               endif
            endif IF_MOYHR

            !VDIR NODEP
            do i = 1, ni

               !# Accumulation of evaporation (in kg/m2)
               tmp_r = dt / CHLC
               zlegaf(i)  = zlegaf(i)  + zleg(i)  * tmp_r
               zleraf(i)  = zleraf(i)  + zler(i)  * tmp_r
               zletraf(i) = zletraf(i) + zletr(i) * tmp_r
               zlevaf(i)  = zlevaf(i)  + zlev(i)  * tmp_r
               zlesaf(i)  = zlesaf(i)  + zles(i)  * dt / (CHLC+CHLF)

               !# Accumulation of drained water
               !#  (soil base water flux, in kg/m2 or mm);
               !#  factor 1000 is for density of water.
               zdrainaf(i) = zdrainaf(i) - 1000. * zdrain(i) * zrootdp(i)

               !# Accumulation of surface runoff (in kg/m2 or mm)
               zoverflaf(i) = zoverflaf(i) + zoverfl(i)

               !# Accumulation of upwards surface water flux (in kg/m2 or mm)
               zwfluxaf(i) = zwfluxaf(i) + zwflux(i) * dt

            enddo

         endif IF_KOUNT0

      endif IF_ISBA

      ! SVS only accumulators
      IF_SVS: if (schmsol == 'SVS') then
         
         ! Reset accumulators at t=T+00hr 
         IF_RESET_SVS: if (kount == 0 .or. lacchr) then
            
            ! Reset accumulators at t=T+00hr 
            
            !Accumulation of actual surf. evap. (kg/m2 or mm)    
            zaccevap(:) = 0.
            
            !Accumulation  of lateral flow at all levels (kg/m2 = mm) 
            zlatflaf(:) = 0.
            
         endif IF_RESET_SVS

         IF_KOUNT_NE_0_SVS:if (kount /= 0) then

            ! E. Gaborit: 
            ! Taken out from svs_update since latflw and watflow have to be 
            ! modified prior to the accumulation in case a lake model is run
            do i=1,ni
               !       ACCUMULATE LATERAL FLOW FOR EACH SOIL LAYER
               do k=1,khyd
                  zlatflaf(i)  = zlatflaf(i) + zlatflw(i,k)       
               enddo
               !       ACCUMULATION OF DRAINAGE (BASE FLOW)
               ! Drainage is the vertical flow across the bottom of the lowermost 
               ! active layer: level = # khyd + 1
               zdrainaf(i) = zdrainaf(i) + zwatflow(i,khyd+1)
            enddo           

         ENDIF IF_KOUNT_NE_0_SVS

      endif IF_SVS

      ! CSLM only accumulators
      IF_CSLM: if (schmlake == 'CSLM') then
         
         ! Reset accumulators at t=T+00hr 
         IF_RESET_CSLM: if (kount == 0 .or. lacchr) then
            zrofinlakaf(:)=0.0
         ENDIF IF_RESET_CSLM

         IF_KOUNT_NE_0_CSLM:if (kount /= 0) then
            do i=1,ni
               ! Accumulate ROFINLAK for water balance verification
               zrofinlakaf(i) = zrofinlakaf(i) + zrofinlak(i)
            enddo
         ENDIF IF_KOUNT_NE_0_CSLM
         
      ENDIF IF_CSLM

      
      if (timings_L) call timing_stop_omp(480)
      call msg_toall(MSG_DEBUG, 'sfc_calcdiag [END]')
      !-------------------------------------------------------------------
      return
   end subroutine sfc_calcdiag3

end module sfc_calcdiag
