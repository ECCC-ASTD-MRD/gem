
module phy_step_mod
   use wb_itf_mod, only: WB_OK, WB_IS_OK, wb_get
   use clib_itf_mod, only: clib_toupper
   use str_mod, only: str_concat
   use series_mod, only: series_stepinit, series_stepend
   use phy_status, only: phy_error_L, phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use phy_options, only: delt, sgo_tdfilter, lhn_filter, sfcflx_filter_order, debug_initonly_L, nphyoutlist, nphystepoutlist, phyoutlist_S, phystepoutlist_S
   use phygridmap, only: phydim_ni, phydim_nj, phydim_nk
   use physlb, only: physlb1
   use phymem, only: phymem_find_idxv, phymem_getmeta, phymeta

#ifdef HAVE_NEMO
   use cpl_itf   , only: cpl_step
#endif

   use ens_perturb, only: ens_spp_stepinit, ENS_OK
   
   private
   public :: phy_step
   
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>
#include <rmn/msg.h>

contains

  !/@*
  function phy_step (F_stepcount, F_stepdriver) result(F_istat)
    implicit none

    !@objective Apply the physical processes: CMC/RPN package
    integer, intent(in) :: F_stepcount     !Step kount incrementing from 0
    integer, intent(in) :: F_stepdriver    !Step number of the driving model
    integer :: F_istat                     !Return status (RMN_OK or RMN_ERR)

    !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014

    !@revision
    !*@/

    include "physteps.cdk"

    integer, external :: phyput_input_param, sfc_get_input_param

    integer, save :: pslic
    logical, save :: do_phyoutlist_L = .true.

    integer :: istat, sfcflxfilt, itmp
    real :: gwd_sig, lhn_sig
    !---------------------------------------------------------------
    F_istat = RMN_ERR
    if (phy_init_ctrl == PHY_NONE) then
       F_istat = PHY_NONE
       return
    else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
       call msg(MSG_ERROR,'(phy_step) Physics not properly initialized.')
       return
    endif
    
    if (debug_initonly_L) then
       call msg(MSG_WARNING,'(phy_step) debug_initonly_L - skipping')
       F_istat = RMN_OK
       return
    endif

    IF_KOUNT0: if (F_stepcount == 0) then
      if (.not.WB_IS_OK(wb_get('dyn/sgo_tdfilter',gwd_sig))) gwd_sig = -1.
      if (sgo_tdfilter /= gwd_sig) then
         call msg(MSG_ERROR, '(phy_step) sgo_tdfilter requested in but not supported by dyn')
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
    endif IF_KOUNT0

    istat = series_stepinit(F_stepcount)
    if (.not.RMN_IS_OK(istat)) &
         call msg(MSG_ERROR,'(phy_step) Problem in series step init')

    if (ens_spp_stepinit(F_stepcount) /= ENS_OK) &
         call msg(MSG_ERROR,'(phy_step) Problem in ensemble step init')

    if (do_phyoutlist_L) then
       do_phyoutlist_L = .false.
       call priv_vname2oname(F_stepdriver, nphyoutlist, phyoutlist_S)
    endif

    if (.not.associated(phystepoutlist_S)) &
         allocate(phystepoutlist_S(max(1,nphyoutlist)))
    phystepoutlist_S(:) = ' '
    nphystepoutlist = -1
    if (nphyoutlist > 0) then
       istat = wb_get('itf_phy/PHYSTEPOUT_N', nphystepoutlist)
       if (.not.WB_IS_OK(istat)) nphystepoutlist = -1
       if (nphystepoutlist > 0) then
          istat = wb_get('itf_phy/PHYSTEPOUT_V', phystepoutlist_S, itmp)
          if (.not.WB_IS_OK(istat)) nphystepoutlist = -1
       endif
    endif
    nphystepoutlist = min(nphystepoutlist, size(phystepoutlist_S))
    call priv_vname2oname(F_stepdriver, nphystepoutlist, phystepoutlist_S, phyoutlist_S(1:nphyoutlist))

    pslic = 0
    step_kount  = F_stepcount
    step_driver = F_stepdriver
    istat = WB_OK
    istat = min(phyput_input_param(),istat)
    istat = min(sfc_get_input_param(),istat)
    if (istat /= WB_OK) call msg(MSG_ERROR,'(phy_step)')

#ifdef HAVE_NEMO
    call cpl_step(F_stepcount, F_stepdriver)
#endif

!$omp parallel
    call physlb1(F_stepcount, phydim_ni, phydim_nj, phydim_nk, pslic)
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

  
  subroutine priv_vname2oname(F_stepdriver, F_n, F_olist_S, F_olistfull_S)
     integer, intent(IN) :: F_stepdriver
     integer, intent(INOUT) :: F_n
     character(len=*), intent(INOUT) :: F_olist_S(:)
     character(len=*), intent(IN), optional :: F_olistfull_S(:)

     integer :: n, nn, istat, ntmp, idxlist(1), msglvl
     character(len=1024) :: msg_S, list_S, tmp_S, prefix_S
     type(phymeta), pointer    :: vmeta
     !---------------------------------------------------------------
     prefix_S = '(PHYOUTREQ)'
     msglvl = MSG_INFO
     if (present(F_olistfull_S)) then
        prefix_S = '(PHYOUTSTEP)'
        msglvl = MSG_INFOPLUS
     endif
     
     ntmp = 0
     DO_OLIST: do n=1,F_n
        nn = phymem_find_idxv(idxlist, F_olist_S(n), F_npath='VO', F_bpath='EPVD', F_quiet=.true.)
        if (nn > 0) then
           istat = phymem_getmeta(vmeta, idxlist(1))
           if (.not.RMN_IS_OK(istat)) nn = 0
        endif
        if (nn <= 0) then
           call msg(msglvl, trim(prefix_S)//' Requested var not in bus: '//F_olist_S(n))
           cycle DO_OLIST
        endif
        tmp_S = vmeta%oname(1:4)
        istat = clib_toupper(tmp_S)
        if (ntmp > 0) then
           if (any(tmp_S == F_olist_S(1:ntmp))) then
              call msg(msglvl, trim(prefix_S)//' Dropping dup: '//tmp_S)
              cycle DO_OLIST
           endif
        endif
        if (present(F_olistfull_S)) then
           if (.not.any(tmp_S == F_olistfull_S(:))) then
              call msg(msglvl, trim(prefix_S)//' Not in full list, Dropping : '//tmp_S)
              cycle DO_OLIST
           endif
        endif
        ntmp = ntmp+1
        F_olist_S(ntmp) = tmp_S
     enddo DO_OLIST
     F_n = ntmp
     if (ntmp < size(F_olist_S)) F_olist_S(ntmp+1:) = ' '
     
     if (F_n == 0) then
        list_S = '[No phy output requested]'
     else if (F_n < 0) then
        list_S = '[phyoutlist Not Provided]'
        if (present(F_olistfull_S)) list_S = '[phystepoutlist Not Provided]'
     else
        call str_concat(list_S, F_olist_S(1:F_n), ', ')
     endif
     write(msg_S, '(a, i8,a,i4,a)') '[step=', F_stepdriver, '] (1:', F_n, ')'
     call msg(msglvl, trim(prefix_S)//' '//trim(msg_S)//': '//trim(list_S))
     !---------------------------------------------------------------
     return
  end subroutine priv_vname2oname

end module phy_step_mod

