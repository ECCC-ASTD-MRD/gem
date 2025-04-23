
module phyexe
   private
   public :: phyexe1

contains

!/@*
subroutine phyexe1(pvars, kount, ni, nk, trnch)
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
   use phymem, only: phyvar
   use phystepinit, only: phystepinit3
   use phystepend, only: phystepend1
   use precipitation, only: precipitation4
   use radiation, only: radiation3
   use sfc_calcdiag, only: sfc_calcdiag3
   use surface, only: surface1
   use inichamp, only: inichamp4
   use tendency, only: tendency5
   use turbulence, only: turbulence2
   use lhn_mod, only: lhn2
#ifdef HAVE_MACH
   use chm_headers_mod, only: chm_exe
#endif
   implicit none
!!!#include <arch_specific.hf>
   !@object this is the main interface subroutine for the cmc/rpn unified physics
   !@arguments
   !          - input/output -
   ! pvars    list of all phy vars (meta + slab data)
   !          - input -
   ! trnch    slice number
   ! kount    timestep number
   ! ni       horizontal running length
   ! nk       vertical dimension

   type(phyvar), pointer, contiguous :: pvars(:)
   integer, intent(in) :: trnch, kount, ni, nk

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

   !----------------------------------------------------------------
   write(tmp_S, '(i6,i6,a)') kount, trnch, ' (phyexe)'
   call msg_verbosity_get(iverb)
   if (debug_trace_L) call msg_verbosity(MSG_DEBUG)
   call msg_toall(MSG_DEBUG, trim(tmp_S)//' [BEGIN]')

   nkm1 = nk-1

   call inichamp4(pvars, kount, ni, nk)
   if (phy_error_L) return

   call phystepinit3(pvars, delt, kount, ni, nk, trnch)
   if (phy_error_L) return

   call radiation3(pvars, kount, ni, nk, trnch)
   if (phy_error_L) return

   call metox3(pvars, ni, nk)
   if (phy_error_L) return

   call linoz3(pvars, delt, kount, ni, nk, nkm1)
   if (phy_error_L) return

   call gwd9(pvars, std_p_prof, delt, kount, ni, nk, nkm1, trnch)
   if (phy_error_L) return

   call apply_rad_tendencies1(pvars, ni, nk, nkm1)
   if (phy_error_L) return

   call surface1(pvars, delt, kount, ni, nk, trnch)
   if (phy_error_L) return

   call turbulence2(pvars, delt, kount, ni, nk, trnch)
   if (phy_error_L) return

   call precipitation4(pvars, delt, kount, ni, nk)
   if (phy_error_L) return

   call phystepend1(pvars, ni, nk)
   if (phy_error_L) return

   call tendency5(pvars, delt, kount, ni, nk)
   if (phy_error_L) return

   call lhn2(pvars, delt, kount, ni, nk)
   if (phy_error_L) return

   call ens_ptp_apply(pvars, kount, ni, nk)
   if (phy_error_L) return

   call calcdiag1(pvars, delt, kount, ni, nk)
   if (phy_error_L) return

   call sfc_calcdiag3(pvars, moyhr, acchr, delt, kount, step_driver, ni)
   if (phy_error_L) return

#ifdef HAVE_MACH
   call chm_exe(pvars, trnch, kount)
   if (phy_error_L) return
#endif
   
   call diagnosurf5(pvars, ni, trnch)
   if (phy_error_L) return

   call extdiag3(pvars, ni, nk, trnch)

   call msg_toall(MSG_DEBUG, trim(tmp_S)//' [END]')
   call msg_verbosity(iverb)
   !----------------------------------------------------------------
   return
end subroutine phyexe1

end module phyexe
