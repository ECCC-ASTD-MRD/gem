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
   subroutine sfc_calcdiag3(f,v,fsiz,vsiz,moyhr,acchr,dt,trnch,kount,step_driver,ni,nk)
      use tdpack_const, only: CHLC, CHLF
      use sfc_options
      use sfcbus_mod
      implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object Calculates averages and accumulators of tendencies and diagnostics
      !@Arguments
      !          - Input/Output -
      ! f        permanent bus
      !          - input -
      ! v        volatile (output) bus
      !          - input -
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! trnch    slice number
      ! kount    timestep number
      ! dt       length of timestep
      ! n        horizontal running length
      ! nk       vertical dimension
      integer :: fsiz, vsiz, moyhr, acchr, trnch, kount, step_driver, ni, nk
      real, target :: f(fsiz), v(vsiz)
      real :: dt
      !*@/
#include <msg.h>
      include "sfcinput.cdk"

      real, pointer, dimension(:) :: zdrain, zdrainaf, zisoil, zleg, zlegaf, &
           zler, zleraf, zles, zlesaf, zletr, zletraf, zlev, zlevaf, zoverfl, &
           zoverflaf, zrootdp, zwflux, zwfluxaf, zwsoil, zinsmavg

      integer :: i, moyhr_steps
      real :: moyhri, tmp_r
      !-------------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'sfc_calcdiag [BEGIN]')
      if (timings_L) call timing_start_omp(480, 'sfc_alcdiag', 46)

      IF_ISBA: if (schmsol == 'ISBA') then

#define MKPTR1D(PTR1,NAME2,BUS) nullify(PTR1); if (NAME2 > 0) PTR1(1:ni) => BUS(NAME2:)

         MKPTR1D(zdrain, drain, f)
         MKPTR1D(zdrainaf, drainaf, f)
         MKPTR1D(zisoil, isoil, f)
         MKPTR1D(zleg, leg, v)
         MKPTR1D(zlegaf, legaf, f)
         MKPTR1D(zler, ler, v)
         MKPTR1D(zleraf, leraf, f)
         MKPTR1D(zles, les, v)
         MKPTR1D(zlesaf, lesaf, f)
         MKPTR1D(zletr, letr, v)
         MKPTR1D(zletraf, letraf, f)
         MKPTR1D(zlev, lev, v)
         MKPTR1D(zlevaf, levaf, f)
         MKPTR1D(zoverfl, overfl, v)
         MKPTR1D(zoverflaf, overflaf, f)
         MKPTR1D(zrootdp, rootdp, f)
         MKPTR1D(zwflux, wflux, v)
         MKPTR1D(zwfluxaf, wfluxaf, f)
         MKPTR1D(zinsmavg, insmavg, f)

         nullify(zwsoil)
         if (wsoil > 0) zwsoil(1:ni) => f(wsoil+ni:)

         if (moyhr > 0) then
            if (kount == 0 .or. mod(step_driver-1, moyhr) == 0) then
               zinsmavg(:) = 0.
            endif
         endif

         IF_RESET: if (kount == 0 .or. &
              (acchr > 0 .and. mod(step_driver-1, acchr) == 0) .or. &
              (acchr == 0 .and. step_driver-1 == 0)) then
            zlegaf(:) = 0.
            zleraf(:) = 0.
            zletraf(:) = 0.
            zlevaf(:) = 0.
            zlesaf(:) = 0.
            zdrainaf(:) = 0.
            zoverflaf(:) = 0.
            zwfluxaf(:) = 0.
         endif IF_RESET

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

      if (timings_L) call timing_stop_omp(480)
      call msg_toall(MSG_DEBUG, 'sfc_calcdiag [END]')
      !-------------------------------------------------------------------
      return
   end subroutine sfc_calcdiag3

end module sfc_calcdiag
