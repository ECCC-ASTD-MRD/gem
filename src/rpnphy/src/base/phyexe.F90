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

module phyexe
   private
   public :: phyexe1
   
contains

!/@*
subroutine phyexe1(dbus, fbus, vbus, trnch, kount, ni, nk)
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
   use chm_headers_mod, only: chm_exe
   implicit none
!!!#include <arch_specific.hf>
   !@object this is the main interface subroutine for the cmc/rpn unified physics
   !@arguments
   !          - input -
   ! dbus     dynamics input field
   !          - input/output -
   ! fbus     historic variables for the physics
   !          - output -
   ! vbus     physics tendencies and other output fields from the physics
   !          - input -
   ! trnch    slice number
   ! kount    timestep number
   ! n        horizontal running length
   ! nk       vertical dimension

   integer, intent(in) :: trnch, kount, ni, nk
   real, dimension(:), pointer, contiguous :: dbus, fbus, vbus

   !@author L. Spacek (oct 2011)
   !@notes
   !          phy_exe is called by all the models that use the cmc/rpn
   !          common physics library. it returns tendencies to the
   !          dynamics.
   !*@/
#include <rmn/msg.h>
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

   call init2nan(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0)
   call init2nan(seloc, ficebl)

   nkm1 = nk-1

   call inichamp4(kount, trnch, ni, nk)
   if (phy_error_L) return

   call phystepinit3(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, &
        vbus, dbus, fbus, seloc, delt, kount, trnch, ni, nk)
   if (phy_error_L) return

   call radiation3(dbus, fbus, vbus, ni, nk, kount, trnch)
   if (phy_error_L) return

   call metox3(dbus, vbus, fbus, ni, nk)
   if (phy_error_L) return
   
   call linoz3(dbus, vbus, fbus, delt, kount, ni, nkm1, nk)
   if (phy_error_L) return
   
   call gwd9(dbus, fbus, vbus, std_p_prof, delt, kount, trnch, ni, nk, nkm1)
   if (phy_error_L) return
   
   call apply_rad_tendencies1(dbus, vbus, fbus, ni, nk, nkm1)
   if (phy_error_L) return
   
   call surface1(dbus, fbus, vbus, trnch, kount, delt, ni, nk)
   if (phy_error_L) return

   call turbulence2(dbus, fbus, vbus, ficebl, seloc, delt, kount, trnch, ni, nk)
   if (phy_error_L) return

   call precipitation4(dbus, fbus, vbus, delt, ni, nk, kount, trnch)
   if (phy_error_L) return

   call prep_cw3(fbus, dbus, vbus, ficebl, ni, nk)
   if (phy_error_L) return

   call tendency5(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, vbus, dbus, &
        1./delt, kount, ni, nk)
   if (phy_error_L) return

   call lhn2(dbus, fbus, vbus, delt, ni, nk, kount)
   if (phy_error_L) return

   call ens_ptp_apply(dbus, vbus, fbus, ni, nk, kount)
   if (phy_error_L) return

   call calcdiag1(dbus, fbus, vbus, delt, kount, ni, nk)
   if (phy_error_L) return

   call sfc_calcdiag3(fbus, vbus, moyhr, acchr, delt, kount, step_driver, trnch, ni)
   if (phy_error_L) return

   call chm_exe(dbus, fbus, vbus, trnch, kount)
   if (phy_error_L) return

   call diagnosurf5(ni, trnch)
   if (phy_error_L) return

   call extdiag3(dbus, fbus, vbus, trnch, ni, nk)

   call msg_toall(MSG_DEBUG, trim(tmp_S)//' [END]')
   call msg_verbosity(iverb)
   !----------------------------------------------------------------
   return
end subroutine phyexe1

end module phyexe
