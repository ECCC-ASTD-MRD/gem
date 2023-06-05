!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 

!/@*
function phy_output(F_dateo,F_dt,F_step,F_gid_l,F_gid_g,F_reduc_core,F_outcfg_S,F_basedir_S) result(F_istat)
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use ezgrid_mod
   use hgrid_wb
   use vgrid_wb
   use phy_output_mod, only: phy_output1_4
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   !@objective
   !@arguments
   integer,intent(in) :: F_dateo,F_dt,F_step,F_gid_l,F_gid_g,F_reduc_core(4)
   character(len=*),intent(in) :: F_outcfg_S,F_basedir_S
   !@return
   integer :: F_istat
   !*@/
   logical,save :: init_L = .false.
   integer :: istat,i0,j0,lni,lnj
   integer,target :: ip1(1)
!!$   integer,pointer :: ip1m(:),ip1t(:)
   integer,pointer :: p_ip1(:)
   !---------------------------------------------------------------------
   F_istat = RMN_ERR
   if (phy_init_ctrl == PHY_NONE) then
       F_istat = PHY_NONE
       return
    else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
       call msg(MSG_ERROR,'(phy_output) Physics not properly initialized.')
       return
   endif

   if (.not.init_L) then
      istat = ezgrid_find_ij0(F_gid_l,F_gid_g,i0,j0,lni,lnj)
      istat = hgrid_wb_put('local#',F_gid_g,i0,j0,lni,lnj)
      ip1(1) = 0
      p_ip1 => ip1
      istat = vgrid_wb_put('surf',VGRID_SURF_TYPE,p_ip1)
!!$      nullify(ip1m,ip1t)
!!$      istat = vgd_get(F_vgrid,'VIPM',ip1m)
!!$      istat = vgd_get(F_vgrid,'VIPT',ip1t)
!!$      istat = vgrid_wb_put('ref-m',F_vgrid,ip1m)
!!$      istat = vgrid_wb_put('ref-t',F_vgrid,ip1t)
      !TODO: gmm put p0? for ref-t/ref-m
      !TODO: output_set_postproc(F_id,F_step)
      init_L = .true.
   endif
   F_istat = phy_output1_4(F_dateo,F_dt,F_step,F_outcfg_S,F_basedir_S,F_reduc_core)
   !---------------------------------------------------------------------
   return
end function phy_output
