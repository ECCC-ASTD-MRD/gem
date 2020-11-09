!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it 
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

module phy_step_mod
   private
   public :: phy_step

contains

  !/@*
  function phy_step (F_stepcount, F_stepdriver) result(F_istat)
    use wb_itf_mod, only: WB_OK, WB_IS_OK, wb_get
    use series_mod, only: series_stepinit, series_stepend
    use phy_status, only: phy_error_L, phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
    use phy_options, only: delt, cond_infilter, sgo_tdfilter, lhn_filter, sfcflx_filter_order
    use phygridmap, only: phydim_ni, phydim_nj, phydim_nk
    use physlb_mod, only: physlb1
    use cpl_itf   , only: cpl_step
    use phybus, only: perbus, dynbus, volbus
    use ens_perturb, only: ens_spp_stepinit, ENS_OK
    implicit none

    !@objective Apply the physical processes: CMC/RPN package
    integer, intent(in) :: F_stepcount     !Step kount incrementing from 0
    integer, intent(in) :: F_stepdriver    !Step number of the driving model
    integer :: F_istat                     !Return status (RMN_OK or RMN_ERR)

    !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014

    !@revision
    !*@/
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>
#include <msg.h>

    include "physteps.cdk"

    integer, external :: phyput_input_param, sfc_get_input_param

    integer, save :: pslic

    integer :: istat, sfcflxfilt
    real :: cond_sig, gwd_sig, lhn_sig
    !---------------------------------------------------------------
    F_istat = RMN_ERR
    if (phy_init_ctrl == PHY_NONE) then
       F_istat = PHY_NONE
       return
    else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
       call msg(MSG_ERROR,'(phy_step) Physics not properly initialized.')
       return
    endif
    if (F_stepcount == 0) then
      if (.not.WB_IS_OK(wb_get('dyn/cond_infilter',cond_sig))) cond_sig = -1.
      if (.not.WB_IS_OK(wb_get('dyn/sgo_tdfilter',gwd_sig))) gwd_sig = -1.
      if (cond_infilter /= cond_sig .or. sgo_tdfilter /= gwd_sig) then
         call msg(MSG_ERROR, '(phy_step) cond_infilter or sgo_tdfilter requested in but not supported by dyn')
         return
      endif
      if (.not.WB_IS_OK(wb_get('dyn/lhn_filter', lhn_sig))) lhn_sig = -1.
      if (lhn_filter /= lhn_sig) then
         call msg(MSG_ERROR, '(phy_step) lhn_filter requested in but not supported by dyn')
         return
      endif
      if (.not.WB_IS_OK(wb_get('dyn/sfcflx_filter_order', sfcflxfilt))) sfcflxfilt = -1
      if (sfcflx_filter_order /= sfcflxfilt) then
         call msg(MSG_ERROR, '(phy_step) sfcflx_filter_order requested in but not supported by dyn')
         return
      endif
    endif

    istat = series_stepinit(F_stepcount)
    if (.not.RMN_IS_OK(istat)) &
         call msg(MSG_ERROR,'(phy_step) Problem in series step init')

    if (ens_spp_stepinit(F_stepcount) /= ENS_OK) &
         call msg(MSG_ERROR,'(phy_step) Problem in ensemble step init')
    
    pslic = 0
    step_kount  = F_stepcount
    step_driver = F_stepdriver
    istat = WB_OK
    istat = min(phyput_input_param(),istat)
    istat = min(sfc_get_input_param(),istat)
    if (istat /= WB_OK) call msg(MSG_ERROR,'(phy_step)')

    call cpl_step(F_stepcount, F_stepdriver)

!$omp parallel
    call physlb1(dynbus, perbus, volbus, &
         size(dynbus,1), size(perbus,1), size(volbus,1), &
         F_stepcount, phydim_ni, phydim_nj, phydim_nk, pslic)
!$omp end parallel
    if (phy_error_L) return

    istat = series_stepend()
    if (.not.RMN_IS_OK(istat)) &
         call msg(MSG_ERROR,'(phy_step) Problem in series step end')

    if (phy_error_L) return

    call timing_start2(460, 'phystats', 46)
    call phystats(F_stepcount, delt)
    call timing_stop(460)
    if (phy_error_L) return

    F_istat = RMN_OK
    !---------------------------------------------------------------
    return
  end function phy_step

end module phy_step_mod

