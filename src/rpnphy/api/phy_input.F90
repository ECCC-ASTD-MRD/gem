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

module phy_input_mod
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_tolower
   use wb_itf_mod, only: WB_IS_OK, wb_get_meta, wb_get
   use input_mod, only: input_new, input_nbvar, input_set_basedir, &
        input_set_filename, input_setgridid, input_isvarstep, input_meta, input_get
   use inputio_mod, only: inputio_new, inputio_nbvar, inputio_set_filename, &
        inputio_nbvar, inputio_isvarstep, inputio_meta, inputio_get, &
        INPUT_FILES_GEOP, INPUT_FILES_CLIM, INPUTIO_T
   use mu_jdate_mod
   use ptopo_utils, only: PTOPO_BLOC, PTOPO_IODIST, ptopo_iotype
   use statfld_dm_mod, only: statfld_dm
   use str_mod, only: str_concat, str_encode_num
   use vGrid_Descriptors, only: vgrid_descriptor, vgd_free
   use vgrid_wb, only: vgrid_wb_get, vgrid_wb_put

   use phyfoldmeta_mod, only: phyfold
   use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
   use phygridmap, only: phydim_ni, phydim_nj , phydim_nk, phy_lcl_ni, &
        phy_lcl_nj, phy_lcl_i0, phy_lcl_in, phy_lcl_j0, phy_lcl_jn, phy_lcl_gid, &
        phy_lclcore_gid, drv_glb_gid, phy_glbcore_gid, phy_comm_io_id
   use phyinputdiag_mod, only: phyinputdiag
   use phy_options, only: jdateo, delt, dyninread_list_s, intozot, phystat_input_l, phystat_2d_l, phystat_dble_l, ninblocx, ninblocy, input_type, debug_trace_L, radia
   use physimple_transforms_mod, only: physimple_transforms3d
   use phy_status, only: PHY_NONE, PHY_CTRL_INI_OK, phy_init_ctrl, phy_error_l
   use phy_typedef, only: phymeta

   private
   public :: phy_input

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <mu_gmm.hf>
#include <msg.h>

   include "phyinput.inc"
   include "buses.cdk"

   integer, external :: phyent2per, phyfillbus

   logical, parameter :: IS_DIR = .true.
 
   character(len=32), parameter  :: VGRID_M_S = 'ref-m'
   character(len=32), parameter  :: VGRID_T_S = 'ref-t'
   character(len=32), parameter  :: PHY_VGRID_M_S = 'phy-m'
   character(len=32), parameter  :: PHY_VGRID_T_S = 'phy-t'
   character(len=32), parameter  :: PHY_REFP0_S = 'PHYREFP0:M'
   character(len=32), parameter  :: PHY_REFP0_LS_S = 'PHYREFP0LS:M'

   character(len=32), parameter  :: OZONEFILENAME_S = 'ozone_clim.fst'

contains

   !/@*
   function phy_input(pre_fold_opr_clbk, F_step, F_incfg_S, F_basedir_S, &
        F_geoname_S) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,external :: pre_fold_opr_clbk
      integer,intent(in) :: F_step
      character(len=*) :: F_incfg_S    !- physics_input_table path
      character(len=*) :: F_basedir_S  !- base path for input data file
      character(len=*) :: F_geoname_S  !- name of geophys file
      !@return
      integer :: F_istat
      !@author Michel Desgagne - Spring 2011
      !@revision
      !  2011-06 Stephane Chamberland
      !  2015-09 Stephane Chamberland: use phy_getmeta
      !  2017-09 Stephane Chamberland: use inputio_mod
      !*@/

      integer, save :: inputid = -1
      integer, save :: nbvar = 0
      type(INPUTIO_T), save :: inputobj

      integer :: ivar, istat, istat2, tmidx, nread, iverb, icat, icat0, icat1
      integer :: idt, ismandatory, readlist_nk(PHYINREAD_MAX)
      real, pointer, dimension(:,:)   :: refp0, pw_p0, refp0ls, pw_p0ls
      real, pointer, dimension(:,:,:) :: data, data2
      character(len=4) :: inname_S, inname2_S
      character(len=32) :: varname_S, varname2_S, readlist_S(PHYINREAD_MAX), horiz_interp_S, vgrid_S, str32, refp0_S, refp0ls_S
      character(len=512) :: str512, dummylist_S(10)
      type(vgrid_descriptor) :: vgridm
      type(phymeta) :: meta1, meta2
      type(phymetaplus) :: meta1plus, meta2plus
      real :: vmin, vmax
      ! ---------------------------------------------------------------------
      call msg_verbosity_get(iverb)
      if (debug_trace_L) call msg_verbosity(MSG_DEBUG)

      F_istat = RMN_ERR
      if (phy_init_ctrl == PHY_NONE) then
         F_istat = PHY_NONE
         return
      else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_input) Physics not properly initialized.')
         return
      endif

      istat = fstopc('MSGLVL','WARNIN',RMN_OPT_SET)

      !# Retrieve input from the model dynamics into the dynamics bus
      istat = phyfillbus(F_step)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR, '(phy_input) problem filling buses')
         return
      endif

      !# Read Ozone and chem related entries
      call priv_ozone(F_step)

      !# Init the input module
      call priv_dyninreadlist(F_step, dyninread_list_s)
      idt = nint(delt)
      F_istat = priv_init(inputid, inputobj, nbvar, jdateo, idt, F_incfg_S, &
           F_basedir_S, F_geoname_S, OZONEFILENAME_S)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR, '(phy_input) problem in init')
         return
      endif

      call chm_load_emissions2(F_basedir_S, F_step, inputobj, nbvar)
      if (phy_error_l) then
         F_istat = RMN_ERR
         call msg(MSG_ERROR, '(phy_input) problem in chm_load_emissions')
         return
      endif

      !# Re-set the phy read list
      phyinread_n = 0
      phyinread_jdateo = jdateo
      phyinread_dt = idt
      phyinread_step = F_step
      phyinread_list_nk(:) = 0
      phyinread_list_S(:)  = ''

      nread = 0
      readlist_nk(:) = 0
      readlist_S(:) = ' '

      !# Copy surface ref fields for vertical interpolation on phy grid
      istat = vgrid_wb_get(VGRID_M_S, vgridm, F_sfcfld_S=refp0_S, &
           F_sfcfld2_S=refp0ls_S)
      istat = vgd_free(vgridm)
      nullify(pw_p0, pw_p0ls, refp0, refp0ls)
      istat = gmm_get(refp0_S, pw_p0)
      istat = gmm_get(PHY_REFP0_S, refp0)
      if (refp0ls_S /= '') then
         istat = gmm_get(refp0ls_S, pw_p0ls)
         istat = gmm_get(PHY_REFP0_LS_S, refp0ls)
      endif
      if (associated(refp0) .and. associated(pw_p0)) then
         refp0(:,:) = pw_p0(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn)
      endif
      if (associated(refp0ls) .and. associated(pw_p0ls)) then
         refp0ls(:,:) = pw_p0ls(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn)
      endif

      !# Read-interpolate-fold vars
      tmidx = -1
      F_istat = RMN_OK
      if (input_type /= 'DIST') istat = rpn_comm_bloc(ninblocx, ninblocy)
      VARLOOP: do ivar=1,nbvar
         if (input_type == 'GEM_4.8') then
            istat = input_isvarstep(inputid, ivar, F_step)
         else
            istat = inputio_isvarstep(inputobj, ivar, F_step)
         endif
         if (.not.RMN_IS_OK(istat)) then
            cycle VARLOOP !var not requested at this step
         endif
         if (input_type == 'GEM_4.8') then
            istat = input_meta(inputid, ivar, inname_S, inname2_S, &
                 dummylist_S, horiz_interp_S, F_mandatory=ismandatory, &
                 F_vmin=vmin, F_vmax=vmax, F_cat0=icat0, F_cat1=icat1)
         else
            istat = inputio_meta(inputobj%cfg, ivar, inname_S, inname2_S,  &
                 dummylist_S, horiz_interp_S, F_mandatory=ismandatory, &
                 F_vmin=vmin, F_vmax=vmax, F_cat0=icat0, F_cat1=icat1)
         endif
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(phy_input) problem getting input varname')
            cycle VARLOOP
         endif

         istat = phygetmetaplus(meta1plus,inname_S, F_npath='IOV', &
              F_bpath='EDPV', F_quiet=.true., F_shortmatch=.false.)
         meta1 = meta1plus%meta
         if (RMN_IS_OK(istat) .and. inname2_S /= ' ') then
            istat = phygetmetaplus(meta2plus,inname2_S, F_npath='IOV', &
                 F_bpath='EDPV', F_quiet=.true., F_shortmatch=.false.)
            meta2 = meta2plus%meta
            varname2_S = meta2%vname
         else
            varname2_S = ' '
         endif
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_INFO,'(phy_input) ignoring var, not declared in bus: '//trim(inname_S)//' : '//trim(inname2_S))
            cycle VARLOOP !# var not needed
         endif

         varname_S  = meta1%vname
         if (ismandatory == -1) ismandatory = meta1%init

         vgrid_S = PHY_VGRID_M_S !#MOM/SLB
         if (meta1%stag == BUSPAR_STAG_THERMO .or. meta1%stag == BUSPAR_STAG_ENERGY) vgrid_S = PHY_VGRID_T_S !#THERMO/SLC, ENERGY/SLS

         if (meta1%nk > 1) then
            if (icat0 < 0) then
               icat0 = 1
               icat1 = meta1%fmul
            else
               icat0 = min(max(0, icat0), meta1%fmul)
               icat1 = min(max(icat0, icat1), meta1%fmul)
            endif
         else
            icat0 = 1
            icat1 = 1
         endif
         ICATLOOP: do icat = icat0, icat1
            inname_S = meta1%iname
            if (inname2_S /= ' ') inname2_S = meta2%iname
            if (icat1 > 1) then
               inname_S = trim(inname_S)//trim(str_encode_num(icat))
               if (inname2_S /= ' ') &
                    inname2_S = trim(inname2_S)//trim(str_encode_num(icat))
            endif
!!$            print *,'(phy_input)0 ',trim(varname_S),icat,'/',icat1,':',trim(inname_S),meta1%nk,meta1%fmul ; call flush(6)

            nullify(data, data2)
            if (input_type == 'GEM_4.8') then
               istat = input_get(inputid, ivar, F_step, phy_lcl_gid, vgrid_S, data, data2, F_ovname1_S=inname_S, F_ovname2_S=inname2_S) !#TODO: ovname
            else
               istat = inputio_get(inputobj, ivar, F_step, data, data2, &
                    F_vgrid_S=vgrid_S, F_ovname1_S=inname_S, F_ovname2_S=inname2_S)
            endif
!!$            print *,'(phy_input)1 ',trim(varname_S),icat,'/',icat1,':',trim(inname_S),meta1%nk,meta1%fmul ; call flush(6)
            if (.not.(RMN_IS_OK(istat) .and. &
                 associated(data) .and. &
                 (inname2_S == ' ' .or. associated(data2)))) then
               if (associated(data)) deallocate(data,stat=istat)
               if (associated(data2)) deallocate(data2,stat=istat)
               nullify(data, data2)
               if (ismandatory == 0) then
                  call msg(MSG_WARNING,'(phy_input) missing optional var: '//trim(inname_S)//' : '//trim(varname_S)//' ('//trim(inname2_S)//' : '//trim(varname2_S)//')')
                  cycle VARLOOP
               endif
               call msg(MSG_ERROR,'(phy_input) missing var: '//trim(inname_S)//' : '//trim(varname_S))
               F_istat = RMN_ERR
!!$            cycle VARLOOP
               call msg_verbosity(iverb)
               return
            endif

            if (inname_S == 'tm') tmidx = meta1plus%index

            vmin = max(vmin, meta1%vmin)
            vmax = min(vmax, meta1%vmax)
            istat = priv_fold(F_step, varname_S, inname_S, meta1%bus, data, icat, &
                 &           readlist_nk, readlist_S, nread, horiz_interp_S, &
                 &           meta1%n(3), vmin, vmax, pre_fold_opr_clbk)
            if (inname2_S /= ' ' .and. RMN_IS_OK(istat)) then
               istat = priv_fold(F_step, varname2_S, inname2_S, meta2%bus, data2, icat, &
                    &           readlist_nk, readlist_S, nread, horiz_interp_S, &
                    &           meta2%n(3), vmin, vmax, pre_fold_opr_clbk)
            endif

#ifndef __GFORTRAN__
            if (associated(data)) deallocate(data,stat=istat2)
            if (associated(data2)) deallocate(data2,stat=istat2)
#endif
            nullify(data, data2)

            call collect_error(istat)
            if (.not.RMN_IS_OK(istat)) then
               F_istat = RMN_ERR
               exit VARLOOP
            endif
         enddo ICATLOOP
      enddo VARLOOP

!!$      F_istat = min(priv_checklist(readlist_nk,readlist_S,nread,F_step),F_istat)
      F_istat = min(priv_checklist(readlist_S,nread,F_step),F_istat)

      if (RMN_IS_OK(F_istat) .and. nread > 0) then
         call msg(MSG_INFO,'(phy_input) All needed var were found')
         F_istat = phyent2per(readlist_S, nread, F_step)
      endif

      phyinread_n = nread
      if (nread > 0) then
         phyinread_list_nk(1:nread) = readlist_nk(1:nread)
         phyinread_list_S(1:nread)  = readlist_S(1:nread)
         str512 = ' '
         do ivar = 1, nread
            write(str32,'(a,i2,a)') trim(readlist_S(ivar))//'(',phyinread_list_nk(ivar),')'
            str512 = trim(str512)//', '//trim(str32)
         enddo
         call msg(MSG_INFO,'(phy_input) Read: '//trim(str512))
      endif
      call msg_verbosity(iverb)
      ! ---------------------------------------------------------------------
      return
   end function phy_input


   !/@*
   function priv_init(my_inputid, my_inputobj, my_nbvar, my_jdateo, my_idt, &
        my_incfg_S, my_basedir_S, my_geoname_S, my_ozonename_S) result(my_istat)
      implicit none
      integer, intent(inout) :: my_inputid
      type(INPUTIO_T), intent(inout) :: my_inputobj
      integer, intent(inout) :: my_nbvar
      integer(INT64), intent(in) :: my_jdateo
      integer, intent(in) :: my_idt
      character(len=*), intent(in) :: my_incfg_S, my_basedir_S, my_geoname_S, my_ozonename_S
      integer :: my_istat
      !*@/
      logical, save :: is_init_L = .false.
      integer, save :: istatus = RMN_ERR

      integer :: istat, iotype

      character(len=32) :: refp0ls_S, refp0ls2_S

      type(vgrid_descriptor) :: vgridm, vgridt
      type(gmm_metadata) :: mymeta

      integer, pointer :: ip1list(:), ip1list2(:)
      real, pointer, dimension(:,:)   :: refp0, refp0ls
      ! ---------------------------------------------------------------------
      my_istat = istatus
      if (input_type == 'GEM_4.8') ptopo_iotype = PTOPO_BLOC
      if (is_init_L) return
      is_init_L = .true.
 
      if (any((/phydim_ni, phydim_nj/) == 0)) then
         call msg(MSG_ERROR, '(phy_input) problem getting bus size')
         my_istat = RMN_ERR
         return
      endif
      my_istat = RMN_OK
      IF_GEM48: if (input_type == 'GEM_4.8') then
         my_inputid = input_new(my_jdateo, my_idt, my_incfg_S)
         istat = my_inputid
         if (RMN_IS_OK(istat)) then
            call phyinputdiag(my_inputid)
            my_nbvar = input_nbvar(my_inputid)
            istat = input_set_basedir(my_inputid, my_basedir_S)
            istat = min(input_setgridid(my_inputid, phy_lclcore_gid), istat)
         endif
         if (RMN_IS_OK(istat)) then
            istat = input_set_filename(my_inputid, 'geop', my_geoname_S, &
                 IS_DIR, INPUT_FILES_GEOP)
            if (radia /= 'NIL') then
               istat = min(istat, &
                    input_set_filename(my_inputid, 'ozon', my_ozonename_S, &
                    IS_DIR, INPUT_FILES_CLIM))
            endif
         endif
      else
         if (input_type == 'BLOC') then
            iotype = PTOPO_BLOC
            istat = inputio_new(my_inputobj, my_jdateo, my_idt, my_incfg_S, &
                 my_basedir_S, phy_lcl_gid, phy_lclcore_gid, phy_comm_io_id, &
                 F_li0=1, F_lin=phy_lcl_ni, F_lj0=1, F_ljn=phy_lcl_nj, &
                 F_iotype=iotype)
         else
            iotype = PTOPO_IODIST
            istat = inputio_new(my_inputobj, my_jdateo, my_idt, my_incfg_S, &
                 my_basedir_S, drv_glb_gid, phy_glbcore_gid, phy_comm_io_id, &
                 F_li0=phy_lcl_i0, F_lin=phy_lcl_in, F_lj0=phy_lcl_j0, F_ljn=phy_lcl_jn, &
                 F_iotype=iotype)
         endif
         if (RMN_IS_OK(istat)) then
            call phyinputdiag(my_inputobj)
            my_nbvar = inputio_nbvar(my_inputobj)
            istat = inputio_set_filename(my_inputobj, 'geop', my_geoname_S, &
                 IS_DIR, INPUT_FILES_GEOP)
            if (radia /= 'NIL') then
               istat = min(istat, &
                    inputio_set_filename(my_inputobj, 'ozon', my_ozonename_S, &
                    IS_DIR, INPUT_FILES_CLIM))
            endif
         endif
      endif IF_GEM48
      if (.not.RMN_IS_OK(istat)) &
           call msg(MSG_ERROR, '(phy_input) input module initialization problem')
      my_istat = min(istat, my_istat)

      nullify(ip1list)
      istat = vgrid_wb_get(VGRID_M_S, vgridm, ip1list, F_sfcfld2_S=refp0ls_S)
      refp0ls2_S = ''
      if (refp0ls_S /= '') refp0ls2_S = PHY_REFP0_LS_S
      ip1list2 => ip1list
      if (size(ip1list) > phydim_nk) then
         !#TODO: check: should we keep phydim_nk instead of phydim_nk+1 (diag level)
         ip1list2(phydim_nk) = ip1list2(phydim_nk+1)
         ip1list2 => ip1list(1:phydim_nk)
      endif
      istat = vgrid_wb_put(PHY_VGRID_M_S, vgridm, ip1list2, PHY_REFP0_S, &
           refp0ls2_S, F_overwrite_L=.true.)
      istat = vgd_free(vgridm)
      if (associated(ip1list)) deallocate(ip1list,stat=istat)

      nullify(ip1list, ip1list2)
      istat = vgrid_wb_get(VGRID_T_S, vgridt, ip1list, F_sfcfld2_S=refp0ls_S)
      refp0ls2_S = ''
      if (refp0ls_S /= '') refp0ls2_S = PHY_REFP0_LS_S
      ip1list2 => ip1list
      if (size(ip1list) > phydim_nk) then
         !#TODO: check: should we keep phydim_nk instead of phydim_nk+1 (diag level)
         ip1list2(phydim_nk) = ip1list2(phydim_nk+1)
         ip1list2 => ip1list(1:phydim_nk)
      endif
      istat = vgrid_wb_put(PHY_VGRID_T_S, vgridt, ip1list2, PHY_REFP0_S, &
           refp0ls2_S, F_overwrite_L=.true.)
      istat = vgd_free(vgridt)
      if (associated(ip1list)) deallocate(ip1list,stat=istat)

      mymeta = GMM_NULL_METADATA
      mymeta%l(1) = gmm_layout(1,phy_lcl_ni,0,0,phy_lcl_ni)
      mymeta%l(2) = gmm_layout(1,phy_lcl_nj,0,0,phy_lcl_nj)
      nullify(refp0, refp0ls)
      istat = gmm_create(PHY_REFP0_S, refp0, mymeta)
      if (refp0ls2_S /= '') &
           istat = gmm_create(PHY_REFP0_LS_S, refp0ls, mymeta)

      call collect_error(my_istat)
      istatus = my_istat
      ! ---------------------------------------------------------------------
      return
   end function priv_init

   
   !/@*
   subroutine priv_dyninreadlist(my_step, my_list_s)
      implicit none
      integer,intent(in) :: my_step
      character(len=*), pointer :: my_list_s(:)
      integer :: istat, ivar, type1, sizeof1, ntr, options1
      character(len=256) :: str256
      ! ---------------------------------------------------------------------
      !#TODO: rework when dyn set itf_phy/READ_TRACERS every step
      if (.not.associated(my_list_s) .or. my_step == 0) then
         ntr   = 0
         istat = wb_get_meta('itf_phy/READ_TRACERS', type1, sizeof1, ntr, options1)
         if (.not.WB_IS_OK(istat)) ntr = 0
         allocate(my_list_s(max(1,ntr)))
         my_list_s(:) = ''
         if (ntr > 0) then
            istat = wb_get('itf_phy/READ_TRACERS', my_list_s, ntr)
            if (.not.WB_IS_OK(istat)) then
               call msg (MSG_ERROR, &
                    '(phy_input) Unable to retrieve dyn read tracers WB entry')
               my_list_s(:) = ''
               return
            endif
            do ivar = 1, ntr
               istat = clib_tolower(my_list_s(ivar))
            end do
            call str_concat(str256, my_list_s, ', ')
            call msg(MSG_INFO, &
                 '(phy_input) List of dyn read tracers: '//trim(str256))
         else
            call msg(MSG_INFO,'(phy_input) No list of dyn read tracers found')
         endif
      endif
      if (associated(my_list_s) .and. my_step /= 0) my_list_s(:) = ''  
      ! ---------------------------------------------------------------------
      return
   end subroutine priv_dyninreadlist


   !/@*
   subroutine priv_ozone(my_step)
      implicit none
      !@objective check for daily update to climatological ozone
      integer,intent(in) :: my_step
      !*@/
      integer, save :: curdd0 = -1
      integer, save :: curmo0 = -1
      integer :: curdd, curmo
      integer(INT64) :: dt_8, istep_8, jdatev
      ! ---------------------------------------------------------------------
      if (radia == 'NIL') return
      
      if (intozot) then
         dt_8    = delt
         istep_8 = my_step
         jdatev  = jdateo + istep_8 * dt_8
         curdd   = jdate_day_of_month(jdatev)
         curmo   = jdate_month(jdatev)
         if (curdd /= curdd0 .or. curmo /= curmo0) then
            call intozon2(curdd, curmo)
            curdd0 = curdd
            curmo0 = curmo
         endif
      else
         if (curdd0 == -1 .or. curmo0 == -1) then
            curdd = 15  !#15th of the month for historical/legacy reasons
            !#TODO: should we do this instead? curdd = jdate_day_of_month(jdateo)
            curmo = jdate_month(jdateo)
            call intozon2(curdd, curmo)
            curdd0 = curdd
            curmo0 = curmo
         endif
      endif
      ! ---------------------------------------------------------------------
      return
   end subroutine priv_ozone


   !/@*
   function priv_fold(my_step, my_varname_S, my_inname_S, my_bus_S, my_data, &
        my_icat, my_readlist_nk, my_readlist_S, my_nread, my_horiz_interp_S, &
        my_nk, my_vmin, my_vmax, pre_fold_opr_clbk) result(my_istat)
      implicit none
      character(len=*) :: my_varname_S, my_inname_S, my_bus_S, &
           my_readlist_S(:), my_horiz_interp_S
      real, dimension(:,:,:), pointer :: my_data
      integer :: my_step, my_icat, my_readlist_nk(:), my_nread, my_nk, my_istat
      real    :: my_vmin, my_vmax
      integer,external :: pre_fold_opr_clbk
      !*@/
      character(len=64) :: msg_S, name_S
      integer :: minxyz(3), maxxyz(3), istat, k, stat_precision
      real, pointer :: data2(:,:,:)
      real, allocatable, target :: dataarr(:,:,:)
      ! ---------------------------------------------------------------------
      minxyz = lbound(my_data)
      maxxyz = ubound(my_data)
      call physimple_transforms3d(my_varname_S, my_inname_S, my_data)
      allocate(dataarr(minxyz(1):maxxyz(1),minxyz(2):maxxyz(2),minxyz(3):maxxyz(3)))
      dataarr(:,:,:) = my_data(:,:,:)
      my_istat = pre_fold_opr_clbk(dataarr, my_varname_S, my_horiz_interp_S, &
           minxyz(1), maxxyz(1), minxyz(2), maxxyz(2), minxyz(3), maxxyz(3))
      
      dataarr(:,:,:) = min(max(my_vmin, dataarr(:,:,:)), my_vmax)

      name_S = my_varname_S
      if (my_icat > 1) write(name_S,'(a,a,i2.2,a)') trim(my_varname_S),'[',my_icat,']'
      IF_STATS: if (phystat_input_l) then
         stat_precision = 4
         if (phystat_dble_l) stat_precision = 8
         if (phystat_2d_l) then
            do k = minxyz(3), maxxyz(3)
               write(msg_S,'(a,i4.4)') trim(my_inname_S)//' => '//trim(name_S)//' ',k
               data2 => dataarr(:,:,k:k)
               call statfld_dm(data2, msg_S, my_step, 'phy_input', stat_precision)
            enddo
         else
            write(msg_S,'(a)') trim(my_inname_S)//' => '//trim(name_S)
            call statfld_dm(dataarr, msg_S, my_step, 'phy_input', stat_precision)
         endif
      endif IF_STATS

      if (maxxyz(3) == minxyz(3) .and. my_nk > 1) then
         !# Put read data into diag level
         data2(1:,1:,my_nk:) => dataarr(:,:,minxyz(3):)
         minxyz = lbound(data2)
         maxxyz = ubound(data2)
         maxxyz(3) = my_nk
         minxyz(3) = maxxyz(3)
      else
         data2(1:,1:,1:) => dataarr(:,:,:)
         minxyz = lbound(data2)
         maxxyz = ubound(data2)
      endif
      !#TODO: re-use metadata from above, use phyfoldmeta
      if (my_icat == 1) then
         istat = phyfold(data2, trim(my_varname_S), trim(my_bus_S), minxyz, maxxyz)
      else
         istat = phyfold(data2, trim(my_varname_S), trim(my_bus_S), minxyz, maxxyz, my_icat)
      endif
      deallocate(dataarr)      
      !#TODO: WARNING: (phyfold) Horizontal sub domaine Not yet supported
      if (.not.RMN_IS_OK(istat)) then
         my_istat = RMN_ERR
         call msg_toall(MSG_WARNING, '(phy_input) Problem folding var: '//trim(name_S))
      endif
      if (RMN_IS_OK(my_istat)) then
         my_nread = min(my_nread + 1,size(my_readlist_S))
         my_readlist_nk(my_nread) = maxxyz(3) - minxyz(3) + 1
         my_readlist_S(my_nread)  = name_S
         istat = clib_tolower(my_readlist_S(my_nread))
      endif
      ! ---------------------------------------------------------------------
      return
   end function priv_fold


!!$   function priv_checklist(F_readlist_nk, F_readlist_S, F_nread, F_step) result(F_istat)
!!$      implicit none
!!$      !@objective Check if all needed var are read
!!$      integer,intent(in) :: F_nread,F_step,F_readlist_nk(:)
      !#TODO: check F_readlist_nk

   !/@*
   function priv_checklist(F_readlist_S, F_nread, F_step) result(F_istat)
      implicit none
      !@objective Check if all needed var are read
      integer,intent(in) :: F_nread,F_step
      character(len=*) :: F_readlist_S(:)
      integer :: F_istat
      !*@/
      logical, parameter :: NOSHORTMATCH_L = .false.
      integer, parameter :: NVARMAX = 512
      integer, parameter :: MUST_INIT = 1
      integer :: nvars,nvars2,ivar,istat
      type(phymetaplus), pointer :: metalist(:)
      character(len=512) :: str512
      ! ---------------------------------------------------------------------
      F_istat = RMN_OK
      if (F_step /= 0) return
      call msg(MSG_INFO,'(phy_input) Checking for mandatory variables.')

      nullify(metalist)
      nvars = phygetmetaplus(metalist, F_name=' ', F_npath='V', F_bpath='EPV', &
           F_maxmeta=NVARMAX, F_quiet=.true., F_shortmatch=NOSHORTMATCH_L)
      !#NOTE: For dynamic bus, the init bit has a different meaning and is checked in phyfillbus

      str512 = ''
      nvars2 = 0
      do ivar = 1,nvars
         if (metalist(ivar)%meta%init == MUST_INIT) then
            nvars2 = nvars2 + 1
            metalist(nvars2)%meta%vname = metalist(ivar)%meta%vname
            istat = clib_tolower(metalist(nvars2)%meta%vname)
            str512 = trim(str512)//', '//trim(metalist(nvars2)%meta%vname)
         end if
      end do
      call msg(MSG_INFO,'(phy_input) Required variables: '//str512(2:len_trim(str512)))

      do ivar=1,F_nread
         istat = clib_tolower(F_readlist_S(ivar))
      enddo

      do ivar = 1,nvars2
         if (F_nread == 0 .or. &
              .not.any(F_readlist_S(1:F_nread) == metalist(ivar)%meta%vname)) then
            F_istat = RMN_ERR
            call msg(MSG_ERROR,'(phy_input) Missing mandatory var (physics_input_table missing entry?): '//trim(metalist(ivar)%meta%vname))
         endif
      end do

      deallocate(metalist,stat=istat)
      ! ---------------------------------------------------------------------
      return
   end function priv_checklist


end module phy_input_mod
