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
subroutine phyexe(d, f, v, dsiz, fsiz, vsiz, trnch, kount, ni, nk)
   use debug_mod, only: init2nan
   use apply_rad_tendencies, only: apply_rad_tendencies1
   use calcdiag, only: calcdiag1
   use diagnosurf, only: diagnosurf5
   use ens_perturb, only: ens_ptp_apply
   use extdiag, only: extdiag3
   use gwd, only: gwd9
   use linoz, only: linoz3
   use metox, only: metox3
   use phy_status, only: phy_error_L
   use phy_options
   use phybus
   use phystepinit, only: phystepinit3
   use precipitation, only: precipitation4
   use prep_cw, only: prep_cw3
   use radiation, only: radiation3
   use sfc_calcdiag, only: sfc_calcdiag3
   use surface, only: surface1
   use tendency, only: tendency5
   use turbulence, only: turbulence2
   use lhn_mod, only: lhn2
   implicit none
!!!#include <arch_specific.hf>
   !@object this is the main interface subroutine for the cmc/rpn unified physics
   !@arguments
   !          - input -
   ! d        dynamics input field
   !          - input/output -
   ! f        historic variables for the physics
   !          - output -
   ! v        physics tendencies and other output fields from the physics
   !          - input -
   ! dsiz     dimension of d
   ! fsiz     dimension of f
   ! vsiz     dimension of v
   ! trnch    slice number
   ! kount    timestep number
   ! n        horizontal running length
   ! nk       vertical dimension

   integer :: dsiz,fsiz,vsiz,trnch,kount,ni,nk
   real    :: d(dsiz), f(fsiz), v(vsiz)

   !@author L. Spacek (oct 2011)
   !@notes
   !          phy_exe is called by all the models that use the cmc/rpn
   !          common physics library. it returns tendencies to the
   !          dynamics.
   !*@/
#include <msg.h>
#include <rmnlib_basics.hf>
   include "tables.cdk"
   include "physteps.cdk"

   integer :: iverb, nkm1
   character(len=64) :: tmp_S

   real, dimension(ni,nk) :: uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0
   real, dimension(ni,nk) :: seloc, ficebl
   !----------------------------------------------------------------
   write(tmp_S, '(i6,i6,a)') kount, trnch, ' (phyexe)'
   call msg_verbosity_get(iverb)
   if (debug_trace_L) call msg_verbosity(MSG_DEBUG)
   call msg_toall(MSG_DEBUG, trim(tmp_S)//' [BEGIN]')

   call init2nan(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, seloc, ficebl)

   nkm1 = nk-1

   call inichamp4(kount, trnch, ni, nk)
   if (phy_error_L) return

   call phystepinit3(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, v, d, f,&
        seloc, delt, vsiz, dsiz, fsiz, kount, trnch, ni, nk)
   if (phy_error_L) return

   call radiation3(d, dsiz, f, fsiz, v, vsiz, ni, nk, kount, trnch)
   if (phy_error_L) return

   TURBULENT_FLUX_CONSISTENCY: if (pbl_flux_consistency) then

      call metox3(d, v, f, dsiz, vsiz, fsiz, ni, nk)
      if (phy_error_L) return

      call linoz3(d, v, f, dsiz, vsiz, fsiz, delt, kount, trnch, ni, nkm1, nk)
      if (phy_error_L) return

      call gwd9(d, f, v, dsiz, fsiz, vsiz, std_p_prof, delt, kount, trnch, ni, nk, nkm1)
      if (phy_error_L) return

      call apply_rad_tendencies1(d, dsiz, v, vsiz, f, fsiz, ni, nk, nkm1)
      if (phy_error_L) return

      call surface1(trnch, kount, delt, ni, nk)
      if (phy_error_L) return

   else

      call surface1(trnch, kount, delt, ni, nk)
      if (phy_error_L) return

      call metox3(d, v, f, dsiz, vsiz, fsiz, ni, nk)
      if (phy_error_L) return

      call linoz3(d, v, f, dsiz, vsiz, fsiz, delt, kount, trnch, ni, nkm1, nk)
      if (phy_error_L) return

      call gwd9(d, f, v, dsiz, fsiz, vsiz, std_p_prof, delt, kount, trnch, ni, nk, nkm1)
      if (phy_error_L) return

      call apply_rad_tendencies1(d, dsiz, v, vsiz, f, fsiz, ni, nk, nkm1)
      if (phy_error_L) return

   endif TURBULENT_FLUX_CONSISTENCY

   call turbulence2(d, f, v, dsiz, fsiz, vsiz, ficebl, seloc, delt, kount, trnch, ni, nk)
   if (phy_error_L) return

   call precipitation4(tplus0, huplus0, d, dsiz, f, fsiz, v, vsiz, delt, ni, nk, kount, trnch)
   if (phy_error_L) return

   call prep_cw3(f, fsiz, d, dsiz, v, vsiz, ficebl, ni, nk)
   if (phy_error_L) return

   call tendency5(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, v, d, &
        1./delt, vsiz, dsiz, kount, ni, nk)
   if (phy_error_L) return

   call lhn2(d, f, v, dsiz, fsiz, vsiz, delt, ni, nk, kount)
   if (phy_error_L) return

   call ens_ptp_apply(d, v, f, dsiz, fsiz, vsiz, ni, nk, kount)
   if (phy_error_L) return

   call calcdiag1(tplus0, huplus0, qcplus0, d, f, v, delt, kount, ni, nk)
   if (phy_error_L) return

   call sfc_calcdiag3(f, v, fsiz, vsiz, moyhr, acchr, delt, kount, step_driver, ni)
   if (phy_error_L) return

   call chm_exe2(d, f, v, dsiz, fsiz, vsiz, trnch, kount)
   if (phy_error_L) return

   call diagnosurf5(ni, trnch)
   if (phy_error_L) return

   call extdiag3(d, f, v, dsiz,fsiz, vsiz, trnch, ni, nk)

   call msg_toall(MSG_DEBUG, trim(tmp_S)//' [END]')
   call msg_verbosity(iverb)
   !----------------------------------------------------------------
   return
end subroutine phyexe
