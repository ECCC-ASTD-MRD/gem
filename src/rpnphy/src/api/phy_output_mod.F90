!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------

!/@*
module phy_output_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_getcwd
   use wb_itf_mod, only: wb_get
   use phy_itf
   use phygridmap, only: phydim_ni, phydim_nk
   use output_mod
   use config_mod
   use mu_jdate_mod
   use hgrid_wb
   use vgrid_wb
   use ptr_store
   use convert_units_mod
   use statfld_dm_mod
   implicit none
   private
   !@objective
   !@author Stephane Chamberland, April 2012
   !@revisions
   ! 001    M. Abrahamowicz    August 2013
   !        Test on 4 character instead of 2 for unit conversions
   !        in phy_output1
   !@public_functions
   public :: phy_output,phy_output0,phy_output1_4,phy_output1_8
   !@public_params
   !@public_vars
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface phy_output
      module procedure phy_output0
      module procedure phy_output1_4
      module procedure phy_output1_8
   end interface

contains

   !/@*
   function phy_output0(F_step) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_step
      !@return
      integer :: F_istat
   !*@/
      character(len=*),parameter :: OUTCFG_NAME = 'outcfg.out'
      integer,save :: reduc_core(4)
      integer(INT64),save :: jdateo = -1
      real(REAL64),save :: dt_8 = 0.D0
      character(len=1024),save :: outcfg_S,basedir_S
      character(len=1024) :: dateo_S,config_dir0_S,pwd_S
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      if (jdateo == -1) then
         F_istat = min(wb_get('time_run_start',dateo_S),F_istat)
         jdateo = jdate_from_print(dateo_S)
         F_istat = min(wb_get('time_dt',dt_8),F_istat)
         F_istat = min(wb_get('path/output',basedir_S),F_istat)
         F_istat = min(wb_get('path/config_dir0',config_dir0_S),F_istat)
         F_istat = min(config_cp2localdir(OUTCFG_NAME,config_dir0_S),F_istat)
         F_istat = min(clib_getcwd(pwd_S),F_istat)
         outcfg_S = trim(pwd_S)//'/'//OUTCFG_NAME
         reduc_core(1:4) = (/1,-1,1,-1/) !#model == core
      endif
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_output) Problem getting config info, cannot produce output')
         jdateo = -1
         return
      endif
      F_istat = phy_output1_8(jdateo,nint(dt_8),F_step,outcfg_S,basedir_S,reduc_core)
      !---------------------------------------------------------------------
      return
   end function phy_output0


   !/@*
   function phy_output1_4(F_dateo,F_dt,F_step,F_outcfg_S,F_basedir_S,F_reduc_core) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_dateo
      integer,intent(in) :: F_dt,F_step,F_reduc_core(4)
      character(len=*),intent(in) :: F_outcfg_S,F_basedir_S
      !@return
      integer :: F_istat
      !*@/
      integer(INT64) :: jdateo
      !---------------------------------------------------------------------
      jdateo = jdate_from_cmc(F_dateo)
      F_istat = phy_output1_8(jdateo,F_dt,F_step,F_outcfg_S,F_basedir_S,F_reduc_core)
      !---------------------------------------------------------------------
      return
   end function phy_output1_4


   !/@*
   function phy_output1_8(F_dateo,F_dt,F_step,F_outcfg_S,F_basedir_S,F_reduc_core) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer(INT64),intent(in) :: F_dateo
      integer,intent(in) :: F_dt,F_step,F_reduc_core(4)
      character(len=*),intent(in) :: F_outcfg_S,F_basedir_S
      !@return
      integer :: F_istat
   !*@/

      character(len=*),parameter :: PHYTAG = 'p'
      character(len=*),parameter :: HGRID_S = 'local#'
      integer, parameter :: MAX_ITEM = 999
      integer,parameter :: STAT_PERCISION = 4

      integer,save :: out_id = -1, p_ni = 0
      integer,save :: l_ijk(3) = (/0,0,0/)
      integer,save :: u_ijk(3) = (/0,0,0/)

      integer :: reduc_core(4),ivar,n_items,istat,k0,kn,grid_id,gi0,gj0,lni,lnj,hx,hy,p_nk,nk2
      character(len=4)   :: mylist_S(MAX_ITEM),outname_S,inname_S
      character(len=512) :: msg_S,vgrid_S,varname_S
      real,pointer :: data3dr4(:,:,:)
      integer,target :: ip1list(MAX_ITEM)
      integer,pointer :: p_ip1list(:)
      type(phymeta) :: mymeta
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      IF_INIT: if (out_id < 0) then
         call msg(MSG_INFO,'(Phy_output) Init Begin')
         reduc_core(1:4) = (/1,-1,1,-1/) !#model == core
         out_id = output_new(F_outcfg_S,PHYTAG,F_dateo,F_dt,F_reduc_core)
         F_istat = min(F_istat,out_id)
         if (.not.RMN_IS_OK(F_istat)) then
            call msg(MSG_ERROR,'(Phy_output) Problem initializing output module')
         endif
         p_ni = phydim_ni; p_nk = phydim_nk
         if (p_ni <= 0) F_istat = RMN_ERR
         istat = RMN_OK
         istat = min(output_set_diag_level(out_id,p_nk),istat)
         istat = min(hgrid_wb_get(HGRID_S,grid_id,gi0,gj0,lni,lnj,hx,hy),istat)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(Phy_output) Problem getting grid info')
         else
            l_ijk(1:2) = (/1-hx,1-hy/) ; l_ijk(3) = 1
            u_ijk(1:2) = (/lni+hx,lnj+hy/) ; u_ijk(3) = 1
         endif
         F_istat = min(F_istat,istat)
         call collect_error(F_istat)
         if (.not.RMN_IS_OK(F_istat)) then
            out_id  = -1
            call msg(MSG_ERROR,'(Phy_output) Problem in Init')
         else
            call msg(MSG_INFO,'(Phy_output) Init End')
         endif
      endif IF_INIT
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return

      F_istat = RMN_OK
      n_items = output_getlist(out_id,F_step,mylist_S)

      write(msg_S,'(a,I5.5)') '(phy_output) Step=',F_step
      if (n_items > 0) then
         call msg(MSG_INFO,trim(msg_S)//' [BEGIN]')
      else
         call msg(MSG_INFO,trim(msg_S)//' No Output')
         return
      endif

      F_istat = min(output_set_basedir(out_id,F_basedir_S),F_istat)

      VARLOOP: do ivar = 1,n_items
         nullify(data3dr4)

         istat = phy_getmeta(mymeta,mylist_S(ivar),F_npath='O',F_quiet=.true.)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING,'(phy_output) Skipping var (not in bus): '//trim(mylist_S(ivar)))
            cycle
         endif

         varname_S = mymeta%vname
         outname_S = mylist_S(ivar)
         inname_S = mymeta%iname
         nk2 = mymeta%nk * mymeta%fmul * max(1,mymeta%mosaic)

         if (nk2 == 1) then
            kn = 1
            vgrid_s = 'surf'
         else
            kn = nk2
            if (mymeta%nk == 1) then
               vgrid_s = 'surf'
               if (mymeta%fmul > 1) then
                  write(vgrid_s,'(a,i3.3)') 'soil',kn
                  if (.not.RMN_IS_OK(vgrid_wb_exists(vgrid_s))) then
                     call priv_soil_ip1list(ip1list,kn)
                     p_ip1list => ip1list(1:kn)
                     istat = vgrid_wb_put(vgrid_s,VGRID_GROUND_TYPE,p_ip1list)
                  endif
               endif
            else
               vgrid_s = 'ref-m'
               if (mymeta%stag > 0) vgrid_s = 'ref-t'
               !TODO-FR: vgrid_s=?? when nk>1.and.fmul>1
               !TODO-FR: vgrid_s=?? when mosaic>0
            endif
         endif
         u_ijk(3) = kn
         call ptr_store_get(data3dr4,l_ijk,u_ijk)
         if (.not.associated(data3dr4)) then
            call msg(MSG_ERROR,trim(msg_S)//' Problem getting workspace for: '//trim(outname_S))
            F_istat = RMN_ERR
            cycle VARLOOP
         endif
         k0 = 1
         l_ijk(3) = k0
!!$         istat = phy_get(data3dr4,mylist_S(ivar),F_npath='O',F_bpath=trim(mymeta%bus),F_start=l_ijk,F_end=u_ijk)
         istat = phy_get(data3dr4,mylist_S(ivar),F_npath='O',F_start=l_ijk,F_end=u_ijk)
         l_ijk(3) = 1
         if (RMN_IS_OK(istat)) then

            !TODO: automate unit convertion convert units, using vardict to get units
            select case(trim(mylist_S(ivar)(1:4)))
            case('sd')
               istat = convert_units(data3dr4,'m','cm')
            case('sndp')
               istat = convert_units(data3dr4,'m','cm')
            case('svdp')
               istat = convert_units(data3dr4,'m','cm')
            case('la')
               istat = convert_units(data3dr4,'rad','deg')
            case('lo')
               istat = convert_units(data3dr4,'rad','deg')
            end select
!!$            call statfld_dm(data3dr4(:,:,k0:kn),mylist_S(ivar),F_step,'phy_output',STAT_PERCISION) !TODO: only if msgLevel > INFOPLUS
            istat = output_writevar(out_id,F_step,mylist_S(ivar),data3dr4,HGRID_S,vgrid_S)
         else
            call msg(MSG_ERROR,trim(msg_S)//' Problem getting data for: '//trim(outname_S))
         endif
         call ptr_store_free(data3dr4)
         F_istat = min(istat,F_istat)
      enddo VARLOOP
      istat = output_close(out_id,F_step)

      call collect_error(F_istat)
      if (n_items>=0) then
         if (RMN_IS_OK(F_istat)) then
            F_istat = n_items
            write(msg_S,'(a,I5,a)') trim(msg_S)//' (nvar=',F_istat,')'
            call msg(MSG_INFO,trim(msg_S)//' [END] ok')
         else
            call msg(MSG_INFO,trim(msg_S)//' [END] with errors')
         endif
      endif
      !---------------------------------------------------------------------
      return
   end function phy_output1_8


   !/@*
   subroutine priv_soil_ip1list(my_ip1list,my_nk)
      implicit none
      !@arguments
      integer :: my_ip1list(:),my_nk
      !*@/
      integer :: k,mykind
      real :: zp1
      !---------------------------------------------------------------------
      mykind = RMN_CONV_ARBITRARY
      do k=1,my_nk
         zp1 = real(k)
         call convip_plus(my_ip1list(k), zp1, mykind, RMN_CONV_P2IPNEW, ' ', .not.RMN_CONV_USEFORMAT_L)
      enddo
      !---------------------------------------------------------------------
      return
   end subroutine priv_soil_ip1list

end module phy_output_mod
