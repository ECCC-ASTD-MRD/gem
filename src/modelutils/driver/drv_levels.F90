!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 

!/@*
module drv_levels_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use wb_itf_mod
   use config_mod
   use vGrid_Descriptors
   use vgrid_wb
   implicit none
   private
   !@objective Manage Levels description
   !@author
   !  Stephane Chamberland, 2012-02
   !@public_functions
   public :: drv_levels_config, drv_levels_init
   !@public_vars
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   include "drv_dyn_itf.inc"

   character(len=*),parameter :: WB_LVL_SEC = 'levels_cfgs/'
   real,parameter :: PREF = 100000.
   integer,save   :: l_nk=0,surf_idx=0
   logical,save :: stag_L=.false.

contains


   !/@*
   function drv_levels_config(F_cfg_basename_S) result(F_istat)
      implicit none
      !@objective Read levels config from file to WB
      !@arguments
      character(len=*),intent(in) :: F_cfg_basename_S  !-
      !@returns
      integer :: F_istat
      !@author
      !  Stephane Chamberland, Feb 2008
   !*@/
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] levels_config')
      F_istat = config_read(F_cfg_basename_S,'levels_cfgs')
      call msg(MSG_DEBUG,'[END] levels_config')
      !---------------------------------------------------------------------
      return
   end function drv_levels_config


   !/@*
   function drv_levels_init() result(F_istat)
      implicit none
      !@objective Initialize drv grid
      !@returns
      integer :: F_istat
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
      !@revisions
      !  2012-02, Stephane Chamberland: RPNPhy offline
   !*@/
      real,pointer :: std_p_prof(:)
      integer,pointer :: ip1_m(:),ip1_t(:),p_ip1(:),p_ip1_m(:),p_ip1_t(:)
      integer,target :: ip1(1)
      integer :: istat,ip1mt(4),i0,in,i1,in1,nkm,nkt
      real :: r14(4)
      character(len=256) :: tmp_S
      type(vgrid_descriptor) :: vcoor
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] levels_init')
      F_istat = dyn_levels_init(l_nk,vcoor,stag_L,surf_idx)      

      nullify(std_p_prof,ip1_m,ip1_t,p_ip1,p_ip1_m,p_ip1_t)
      IF_LVL_OK: if (RMN_IS_OK(F_istat)) then
         F_istat = min(vgd_get(vcoor,key='VIPM',value=ip1_m,quiet=.true.),F_istat)
         !# If type 5005 skip level nk+1 (diag level)
         nkm = min(surf_idx, size(ip1_m))
         p_ip1_m(1:nkm) => ip1_m(1:nkm)
         F_istat = min(vgrid_wb_put('ref-m',vcoor,p_ip1_m),F_istat)
         if (stag_L) then
            F_istat = min(vgd_get(vcoor,key='VIPT',value=ip1_t,quiet=.true.),F_istat)
         else
            ip1_t => ip1_m
         endif
         nkt = min(surf_idx, size(ip1_t))
         p_ip1_t(1:nkt) => ip1_t(1:nkt)
         F_istat = min(vgrid_wb_put('ref-t',vcoor,p_ip1_t),F_istat)
         !TODO: write ip1_m(l_nk) or convip ip1(1) or 0
!!$         ipkind = RMN_CONV_ARBITRARY
!!$         zp1 = 0.
!!$         ipkind = RMN_CONV_SIGMA
!!$         zp1 = 0.
!!$         call convip(ip1(1),zp1,ipkind,RMN_CONV_P2IPNEW,' ',.not.RMN_CONV_USEFORMAT_L)
!!$         ip1(1) = ip1_m(l_nk)
         ip1(1) = 0
         p_ip1 => ip1
         F_istat = min(vgrid_wb_put('surf',VGRID_SURF_TYPE,p_ip1),F_istat)

         F_istat = min(vgd_levels(vcoor,p_ip1_m,std_p_prof,PREF,in_log=.false.),F_istat)
         F_istat = min(wb_put('std_p_prof_m',std_p_prof),F_istat)

         if (RMN_IS_OK(F_istat)) then
            i0 = lbound(p_ip1_m,1); in = ubound(p_ip1_m,1)
            i1 = min(i0+1,in) ; in1 = max(in-1,i0)
            ip1mt(1) = p_ip1_m(i0) ; ip1mt(2) = p_ip1_m(i1)
            ip1mt(3) = p_ip1_m(in1) ; ip1mt(4) = p_ip1_m(in)
            write(tmp_S,'(a,i3,a,i3,a,i9,a,i9,a,i9,a,i9)') &
                 '(drv_levels) ip1_m(',i0,':',in,')=', ip1mt(1),', ', &
                 ip1mt(2),', ...,', ip1mt(3),', ', ip1mt(4)
            call msg(MSG_INFO,tmp_S)
            i0 = lbound(p_ip1_t,1); in = ubound(p_ip1_t,1)
            i1 = min(i0+1,in) ; in1 = max(in-1,i0)
            ip1mt(1) = p_ip1_t(i0) ; ip1mt(2) = p_ip1_t(i1)
            ip1mt(3) = p_ip1_t(in1) ; ip1mt(4) = p_ip1_t(in)
            write(tmp_S,'(a,i3,a,i3,a,i9,a,i9,a,i9,a,i9)') &
                 '(drv_levels) ip1_t(',i0,':',in,')=', ip1mt(1),', ', &
                 ip1mt(2),', ...,', ip1mt(3),', ', ip1mt(4)
            call msg(MSG_INFO,tmp_S)
            i0 = lbound(std_p_prof,1); in = ubound(std_p_prof,1)
            i1 = min(i0+1,in) ; in1 = max(in-1,i0)
            r14(1) = std_p_prof(i0) ; r14(2) = std_p_prof(i1)
            r14(3) = std_p_prof(in1) ; r14(4) = std_p_prof(in)
            write(tmp_S,'(a,i3,a,i3,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') &
                 '(drv_levels) std_p_prof_m(',i0,':',in,')=', r14(1),', ', &
                 r14(2),', ...,', r14(3),', ', r14(4)
            call msg(MSG_INFO,tmp_S)
         endif

         if (stag_L) then
            if (associated(std_p_prof)) deallocate(std_p_prof,stat=istat)
            F_istat = min(vgd_levels(vcoor,p_ip1_t,std_p_prof,PREF,in_log=.false.),F_istat)
         endif
         F_istat = min(wb_put('std_p_prof_t',std_p_prof),F_istat)
         if (RMN_IS_OK(F_istat)) then
            i0 = lbound(std_p_prof,1); in = ubound(std_p_prof,1)
            i1 = min(i0+1,in) ; in1 = max(in-1,i0)
            r14(1) = std_p_prof(i0) ; r14(2) = std_p_prof(i1)
            r14(3) = std_p_prof(in1) ; r14(4) = std_p_prof(in)
            write(tmp_S,'(a,i3,a,i3,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') &
                 '(drv_levels) std_p_prof_t(',i0,':',in,')=', r14(1),', ', &
                 r14(2),', ...,', r14(3),', ', r14(4)
            call msg(MSG_INFO,tmp_S)
         endif
         if (associated(std_p_prof)) deallocate(std_p_prof,stat=istat)
         if (associated(ip1_m)) deallocate(ip1_m,stat=istat)
         if (associated(ip1_t).and.stag_L) deallocate(ip1_t,stat=istat)

         F_istat = min(wb_put(WB_LVL_SEC//'nk',l_nk),F_istat)
         F_istat = min(wb_put(WB_LVL_SEC//'stag_L',stag_L),F_istat)
         F_istat = min(wb_put(WB_LVL_SEC//'surf_idx',surf_idx),F_istat)
      endif IF_LVL_OK

      if (RMN_IS_OK(F_istat)) then
         write(tmp_S,'(a,i4,a)') '(drv_levels) Initialization OK: [nk=',l_nk,']'
         call msg(MSG_INFO,tmp_S)
      else
         call msg(MSG_ERROR,'(drv_levels) Problem in Initialisation')
      endif
      call msg(MSG_DEBUG,'[BEGIN] levels_init')
      !---------------------------------------------------------------------
      return
   end function drv_levels_init


end module drv_levels_mod
