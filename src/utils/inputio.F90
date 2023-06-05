!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

!/@*
module inputio_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   use rmn_gmm, only: gmm_metadata, gmm_get, gmm_create, GMM_NULL_FLAGS, GMM_FLAG_RSTR
   use vGrid_Descriptors, only: vgrid_descriptor, vgd_get, vgd_free, VGD_OK, operator(==)
   use incfg2_mod
   use inputio_files_mod
   use fstmpi_mod, only: fstmpi_rdhint_3d_r4, fstmpi_rdhint_3d_r4_vect
   use fstmpio_mod
   use hinterp4yy_mod, only: HINTERP4YY_NONE
   use mu_jdate_mod, only: jdate_to_cmc, jdate_from_cmc, jdate_month, jdate_year, jdate_midmonth, MU_JDATE_ANY
   use cmcdate_mod, only: cmcdate_toprint
   use ptopo_utils, only: PTOPO_BLOC
   use str_mod, only: str_concat_i
   use time_interp_mod
   use vinterp_mod, only: vinterp
   use vgrid_wb, only: vgrid_wb_get, vgrid_wb_put, vgrid_wb_is_press_kind
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland, 2017-05
   !@description
   ! Public functions
   public :: inputio_new, inputio_add, inputio_nbvar, inputio_meta
   public :: inputio_set_meta, inputio_isvarstep
   public :: inputio_set_filename, inputio_set_basedir, inputio_close_files
   public :: inputio_get, inputio_set_lcl_scope
   ! Public constants
   public :: INPUT_FILES_ANAL, INPUT_FILES_CLIM, INPUT_FILES_GEOP
   public :: INCFG_STRLEN
   public :: INCFG_KNOWN_H_INT, INCFG_KNOWN_V_INT, INCFG_KNOWN_T_INT
   public :: INCFG_KNOWN_LVL_TYPE, INCFG_LVL_ARBITRARY, INCFG_LVL_MOM
   public :: INCFG_LVL_THERMO, INCFG_LVL_TDIAG, INCFG_LVL_MDIAG
   public :: INPUT_ERR, INPUT_SKIP, INPUT_OK
   ! Public Types
   public :: INCFG_T, INPUTIO_FILES_T, INPUTIO_T
   !*@/

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

#define MK_ID2CHAR(ID) trim(achar(65+ID))
#define ADD_PREFIX(ID,VN) 'io/'//MK_ID2CHAR(ID)//'/'//VN
#define ADD_PREFIX2(ID,ID2,VN) 'io/'//MK_ID2CHAR(ID)//MK_ID2CHAR(ID2)//'/'//VN
#define ADD_PREFIX3(ID,ID2,DAT,VN) 'io/'//MK_ID2CHAR(ID)//MK_ID2CHAR(ID2)//'/'//DAT//'/'//VN
#define RM_PREFIX2(VN) VN(7:)
!!$#define RM_PREFIX3(VN) VN(7:)

   integer, parameter :: INPUT_ERR = -1  !# RMN_ERR
   integer, parameter :: INPUT_SKIP = 0
   integer, parameter :: INPUT_OK = 1

   integer, parameter :: NMAX_LEVELS = 1024
   
   character(len=*), parameter :: ALT_VGRID_H_FLDNAME = 'GZ'

   interface inputio_add
      module procedure incfg_add_string
      module procedure incfg_add_kv
   end interface inputio_add

   interface inputio_meta
      module procedure incfg_meta
   end interface inputio_meta

   interface inputio_set_meta
      module procedure incfg_set_meta
   end interface inputio_set_meta

   integer, save :: m_id = 0

   type :: INPUTIO_T
      integer :: id
      logical :: init_L = .false.
      type(INCFG_T) :: cfg
      type(INPUTIO_FILES_T) :: fid
      integer :: l_i0, l_in, l_j0, l_jn
   end type INPUTIO_T

   type :: INPUTIO_FLD_T
      integer :: id
      integer :: nkeys
      logical :: needonce_L
      character(len=64) :: vn1_S, vn2_S
      character(len=64) :: vgrid_S
      character(len=64) :: sfc_S(2) !# SFC + SLS
      character(len=64) :: alt_S    !# alt varname for vgrid def
      integer(INT64) :: jdatev
      logical :: salloc_L, dalloc_L, altalloc_L
      real, dimension(:,:,:), pointer :: sfc, psfc !# SFC + SLS
      real, dimension(:,:,:), pointer :: dalt, palt !# alt vgrid def data
      real, dimension(:,:,:), pointer :: d1, d2
      real, dimension(:,:,:), pointer :: p1, p2
      integer :: k1(NMAX_LEVELS)
      integer :: k2(NMAX_LEVELS)
      integer :: ks(2)
      integer :: kalt(NMAX_LEVELS)
      integer :: hstats(NMAX_LEVELS)
      integer :: hstat
      integer :: hgridid, hgridcoreid
      integer :: l_i0, l_in, l_j0, l_jn
   end type INPUTIO_FLD_T

contains

   !/@*
   function inputio_new(F_inputobj, F_jdateo, F_dt, F_filename_S, &
        F_basedir_S, F_hgridid, F_hgridcoreid, F_commgid, &
        F_ip1list, F_vgrid_m_S, F_vgrid_t_S, &
        F_li0, F_lin, F_lj0, F_ljn, F_iotype) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_T), intent(out) :: F_inputobj
      integer(INT64), intent(in) :: F_jdateo
      integer, intent(in) :: F_dt
      character(len=*), intent(in), optional :: F_filename_S !full/rel path of config file
      character(len=*), intent(in), optional :: F_basedir_S
      integer, intent(in), optional :: F_hgridid, F_hgridcoreid
      integer, intent(in), optional :: F_commgid
      integer, intent(in), optional :: F_ip1list(:)
      character(len=*), intent(in), optional :: F_vgrid_m_S
      character(len=*), intent(in), optional :: F_vgrid_t_S
      integer, intent(in), optional :: F_iotype
      integer, intent(in), optional :: F_li0, F_lin, F_lj0, F_ljn
      !@return
      integer :: F_istat
      !*@/
      integer :: hgridid, hgridcoreid, commgid
      character(len=1024) :: string_S, filename_S, basedir_S, vgrid_m_S, vgrid_t_S
      !----------------------------------------------------------------------
      F_inputobj%init_L = .false.
      F_inputobj%id = -1
      F_inputobj%l_i0 = -1*huge(1)
      F_inputobj%l_j0 = -1*huge(1)
      F_inputobj%l_in = huge(1)
      F_inputobj%l_jn = huge(1)
      if (present(F_li0)) F_inputobj%l_i0 = F_li0
      if (present(F_lj0)) F_inputobj%l_j0 = F_lj0
      if (present(F_lin)) F_inputobj%l_in = F_lin
      if (present(F_ljn)) F_inputobj%l_jn = F_ljn

      filename_S = ''
      basedir_S = ''
      hgridid = -1
      hgridcoreid = -1
      commgid = -1
      vgrid_m_S = ''
      vgrid_t_S = ''
      if (present(F_filename_S)) filename_S = F_filename_S
      if (present(F_basedir_S)) basedir_S = F_basedir_S
      if (present(F_hgridid)) hgridid = F_hgridid
      if (present(F_hgridcoreid)) hgridcoreid = F_hgridcoreid
      if (present(F_commgid)) commgid = F_commgid
      if (present(F_vgrid_m_S)) vgrid_m_S = F_vgrid_m_S
      if (present(F_vgrid_t_S)) vgrid_t_S = F_vgrid_t_S

      write(string_S, '(a,1x,i0,1x,i0,1x,a)') '(inputio) new [BEGIN]', F_jdateo, F_dt, trim(filename_S)
      call msg(MSG_DEBUG, string_S)
 
      if (present(F_ip1list)) then
         F_istat = incfg_new(F_inputobj%cfg, F_jdateo, F_dt, &
              filename_S, F_ip1list, hgridid, hgridcoreid, commgid, &
              vgrid_m_S, vgrid_t_S)
      else
         F_istat = incfg_new(F_inputobj%cfg, F_jdateo, F_dt, &
              filename_S, F_hgridid=hgridid, F_hgridcoreid=hgridcoreid, &
              F_commgid=commgid, F_vgrid_m_S=vgrid_m_S, F_vgrid_t_S=vgrid_t_S)
      endif

      if (RMN_IS_OK(F_istat)) then
         if (present(F_iotype)) then
            F_istat = inputio_files_new(F_inputobj%fid, basedir_S, F_iotype)
         else
            F_istat = inputio_files_new(F_inputobj%fid, basedir_S)
         endif
      endif

      if (RMN_IS_OK(F_istat)) then
         F_inputobj%init_L = .true.
         m_id = m_id + 1
         F_inputobj%id = m_id
      endif

      write(string_S, '(a,1x,i0)') '(inputio) new [END]', F_istat
      call msg(MSG_DEBUG, string_S)
      return
      !----------------------------------------------------------------------
   end function inputio_new


   !/@*
   function inputio_set_lcl_scope(F_inputobj, F_li0, F_lin, F_lj0, F_ljn) &
        result(F_istat)
      implicit none
      type(INPUTIO_T), intent(inout) :: F_inputobj
      integer, intent(in) :: F_li0, F_lin, F_lj0, F_ljn
      integer :: F_istat
      !----------------------------------------------------------------------
      F_istat = RMN_OK
      F_inputobj%l_i0 = F_li0
      F_inputobj%l_in = F_lin
      F_inputobj%l_j0 = F_lj0
      F_inputobj%l_jn = F_ljn
      !----------------------------------------------------------------------
      return
   end function inputio_set_lcl_scope


   !/@*
   function inputio_set_filename(F_inputobj, F_key_S, F_filename_S, &
        F_isdir_L, F_type) result(F_istat)
      implicit none
      type(INPUTIO_T), intent(inout) :: F_inputobj
      character(len=*),intent(in) :: F_key_S,F_filename_S
      logical,optional,intent(in) :: F_isdir_L
      integer,optional,intent(in) :: F_type
      integer :: F_istat
      !*@/
      logical :: isdir_L
      !------------------------------------------------------------------
      isdir_L = .false.
      if (present(F_isdir_L)) isdir_L = F_isdir_L
      if (present(F_type)) then
         F_istat = inputio_files_set_name(F_inputobj%fid, F_key_S, F_filename_S, &
              isdir_L, F_type)
      else
         F_istat = inputio_files_set_name(F_inputobj%fid, F_key_S, F_filename_S, &
              isdir_L)
      endif
      !------------------------------------------------------------------
      return
   end function inputio_set_filename


   !/@*
   function inputio_set_basedir(F_inputobj, F_basedir_S) result(F_istat)
      implicit none
      type(INPUTIO_T), intent(inout) :: F_inputobj
      character(len=*),intent(in) :: F_basedir_S
      integer :: F_istat
      !*@/
      !------------------------------------------------------------------
      F_istat = inputio_files_set_basedir(F_inputobj%fid, F_basedir_S)
      !------------------------------------------------------------------
      return
   end function inputio_set_basedir


   !/@*
   function inputio_nbvar(F_inputobj) result(F_nbvar)
      implicit none
      type(INPUTIO_T), intent(inout) :: F_inputobj
      integer :: F_nbvar
      !*@/
      !------------------------------------------------------------------
      F_nbvar = incfg_nbvar(F_inputobj%cfg)
      !------------------------------------------------------------------
      return
   end function inputio_nbvar


   !/@*
   function inputio_isvarstep(F_inputobj, F_index, F_istep) result(F_istat)
      implicit none
      type(INPUTIO_T), intent(inout) :: F_inputobj
      integer, intent(in) :: F_index, F_istep
      integer :: F_istat
      !*@/
      !------------------------------------------------------------------
      F_istat = incfg_isvarstep(F_inputobj%cfg, F_index, F_istep)
      !----------------------------------------------------------------------
      return
   end function inputio_isvarstep


   !/@*
   function inputio_close_files(F_inputobj) result(F_istat)
      implicit none
      type(INPUTIO_T), intent(inout) :: F_inputobj
      integer :: F_istat
      !*@/
      !----------------------------------------------------------------------
      F_istat = inputio_files_close(F_inputobj%fid)
      !----------------------------------------------------------------------
      return
   end function inputio_close_files


   !/@*
   function inputio_get(F_inputobj, F_ivar, F_step, F_data1, F_data2, &
        F_hgridid, F_hgridcoreid, F_vgrid_S, F_varname_S, F_varname2_S, &
        F_ovname1_S, F_ovname2_S) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_T), intent(inout), target :: F_inputobj
      integer, intent(in) :: F_ivar, F_step
      real, pointer :: F_data1(:,:,:), F_data2(:,:,:)
      integer, intent(in),optional :: F_hgridid, F_hgridcoreid
      character(len=*),intent(in),optional :: F_vgrid_S
      character(len=*),intent(out),optional :: F_varname_S, F_varname2_S
      character(len=*),intent(in),optional :: F_ovname1_S, F_ovname2_S  !#input name overrides
      !@return
      integer :: F_istat
      !*@/
      integer, parameter :: MAXITER = 4
      integer :: istat, niter, nn
      integer(INT64) :: jdatev
      type(INCFG_T), pointer :: cfg
      type(INCFG_VAR_T), pointer :: cfgvar
      type(INPUTIO_FILES_T), pointer :: cfgfile
      type(INPUTIO_FLD_T), target :: fld, fld1, fld2, fld1v, fld2v, fld12vt
      character(len=256) :: vgrid_S, msg_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) get [BEGIN]')
      F_istat = RMN_ERR
      cfg => F_inputobj%cfg
      cfgvar => F_inputobj%cfg%v(F_ivar)
      cfgfile => F_inputobj%fid

      call inputio_files_set_iotype(cfgfile)

      fld = priv_newfld(F_inputobj)
      istat = incfg_time(cfg, F_step, jdatev)
      fld%jdatev = jdatev
      istat = incfg_meta(cfg, F_ivar, fld%vn1_S, fld%vn2_S)
      if (present(F_ovname1_S)) then
         if (F_ovname1_S /= '') then
            fld%vn1_S = F_ovname1_S
            istat = clib_tolower(fld%vn1_S)
            if (present(F_ovname2_S)) then
               if (F_ovname2_S /= '') then
                  fld%vn2_S = F_ovname2_S
                  istat = clib_tolower(fld%vn2_S)
               endif
            endif
         endif
      endif
      if (present(F_varname_S)) F_varname_S = fld%vn1_S
      if (present(F_varname2_S)) F_varname2_S = fld%vn2_S
      fld%hgridid = cfg%hgridid
      fld%hgridcoreid = cfg%hgridcoreid
      if (present(F_hgridid)) then
         fld%hgridid = F_hgridid
         fld%hgridcoreid = F_hgridid
         if (present(F_hgridcoreid)) fld%hgridcoreid = F_hgridcoreid
      endif
      !#TODO: auto vgrid swap from meta?
      vgrid_S = cfg%vgrid_m_S
      if (present(F_vgrid_S)) vgrid_S = F_vgrid_S
 
      msg_S = cfgvar%files_S(1)
      do nn = 2, size(cfgvar%files_S)
         if (cfgvar%files_S(nn) == '') exit
         msg_S = trim(msg_S)//'+'//cfgvar%files_S(nn)(1:4)
      enddo
      call msg(MSG_INFO, '(inputio) Looking for '//trim(fld%vn1_S)//' '// &
           trim(fld%vn2_S)//' in file "'//trim(msg_S)// &
           '" (interpolation h/v/t: '//trim(cfgvar%hint_S)//'/'// &
           trim(cfgvar%vint_S)//'/'//trim(cfgvar%tint_S)//') (tv='// &
           trim(cfgvar%typvar_S)//')')

      istat = -1
      niter = 0
      do while (istat < 0 .and. niter < MAXITER)
         niter = niter + 1
         istat = priv_tint_status(cfgvar, jdatev, fld)
         if (RMN_IS_OK(istat) .or. istat == RMN_ERR) exit
         istat = priv_frh(cfg, F_ivar, cfgvar, cfgfile, fld, istat, niter, vgrid_S)
         if (RMN_IS_OK(istat)) istat = priv_tint_set(cfgvar, fld)
         istat = -1
      enddo
      if (RMN_IS_OK(istat) .and. &
           .not.any(cfgvar%tint_S(1:4) == (/'none', 'any '/))) then
         fld%dalloc_L = .true.
         fld%salloc_L = .true.
         istat = priv_freefld(fld)
      endif

      fld1 = priv_newfld(F_inputobj)
      fld2 = priv_newfld(F_inputobj)
      fld1v = priv_newfld(F_inputobj)
      fld2v = priv_newfld(F_inputobj)
      fld12vt = priv_newfld(F_inputobj)

      if (RMN_IS_OK(istat)) then
         fld1 = priv_newfld(F_inputobj, fld)
         istat = priv_tint_get(cfgvar, fld1, TIME_INTERP_PREV)
         if (RMN_IS_OK(istat)) then
            fld1v = priv_newfld(F_inputobj, fld1, F_reset_L=.true.)
            istat = priv_vint(cfgvar, fld1, fld1v, vgrid_S)
         endif
      endif

      if (RMN_IS_OK(istat)) then
         fld2 = priv_newfld(F_inputobj, fld1, F_reset_L=.true.)
         fld2v = priv_newfld(F_inputobj, fld1v, F_reset_L=.true.)
         if (.not.any(cfgvar%tint_S(1:4) == (/'none', 'any '/))) then
            istat = priv_tint_get(cfgvar, fld2, TIME_INTERP_NEXT)
            if (RMN_IS_OK(istat)) then
               fld2v = priv_newfld(F_inputobj, fld2, F_reset_L=.true.)
               istat = priv_vint(cfgvar, fld2, fld2v, vgrid_S)
            else
               istat = priv_tint_status(cfgvar, jdatev, fld2)
               fld2 = priv_newfld(F_inputobj, fld1, F_reset_L=.true.)
            endif
         endif
      endif

      if (RMN_IS_OK(istat)) then
         fld12vt = priv_newfld(F_inputobj, fld1v, F_reset_L=.true.)
         fld12vt%jdatev = jdatev
         if (RMN_IS_OK(istat)) &
              F_istat = priv_tint(cfgvar, fld1v, fld2v, fld12vt)
      endif

      if (.not.RMN_IS_OK(F_istat)) then
         msg_S = '(inputio) Could not get: '//trim(fld%vn1_S)//' '// &
              trim(fld%vn2_S)//' in file "'//trim(msg_S)// &
              '" (interpolation h/v/t: '//trim(cfgvar%hint_S)//'/'// &
              trim(cfgvar%vint_S)//'/'//trim(cfgvar%tint_S)//') (tv='// &
              trim(cfgvar%typvar_S)//')'
         if (cfgvar%mandatory /= 0) then
            call msg(MSG_WARNING, msg_S)
         else
            call msg(MSG_INFO, msg_S)
         endif
      endif

      call collect_error(F_istat)
      write(msg_S, '(a,1x,i0)') trim(fld%vn1_S), F_istat

      istat = priv_freefld(fld, fld12vt)
      istat = priv_freefld(fld1, fld12vt)
      istat = priv_freefld(fld2, fld12vt)
      istat = priv_freefld(fld1v, fld12vt)
      istat = priv_freefld(fld2v, fld12vt)
      if (RMN_IS_OK(F_istat)) then
         !#TODO?: F_istat = fst_checkalloc(F_data1, F_ni, F_nj, F_nk, F_realloc_L)
         if (associated(fld12vt%d1)) F_data1 => fld12vt%d1
         if (associated(fld12vt%d2) .and. fld%vn2_S /= ' ') F_data2 => fld12vt%d2
      else
         istat = priv_freefld(fld12vt)
      endif
      call msg(MSG_DEBUG, '(inputio) get [END] '//trim(msg_S))
      !----------------------------------------------------------------------
      return
   end function inputio_get


   !==== Private Functions =================================================

   !/@*
   function priv_newfld(F_inputobj, F_fldin, F_reset_L) result(F_fldout)
      implicit none
      type(INPUTIO_T), intent(in) :: F_inputobj
      type(inputio_fld_t), intent(in), optional:: F_fldin
      logical, intent(in), optional:: F_reset_L
      type(inputio_fld_t):: F_fldout
      !*@/
      logical :: reset_L
      !----------------------------------------------------------------------
      F_fldout%id = F_inputobj%id
      F_fldout%nkeys = -1
      F_fldout%needonce_L = .false.
      F_fldout%vn1_S = ''
      F_fldout%vn2_S = ''
      F_fldout%vgrid_S = ''
      F_fldout%sfc_S(:) = ''
      F_fldout%alt_S = ''
      F_fldout%jdatev = MU_JDATE_ANY
      F_fldout%salloc_L = .false.
      F_fldout%dalloc_L = .false.
      F_fldout%altalloc_L = .false.
      nullify(F_fldout%sfc)
      nullify(F_fldout%d1)
      nullify(F_fldout%d2)
      nullify(F_fldout%dalt)
      nullify(F_fldout%psfc)
      nullify(F_fldout%p1)
      nullify(F_fldout%p2)
      nullify(F_fldout%palt)
      F_fldout%k1(:) = -1
      F_fldout%k2(:) = -1
      F_fldout%ks(:) = -1
      F_fldout%kalt(:) = -1
      F_fldout%hstats = -1
      F_fldout%hstat = -1
      F_fldout%hgridid = -1
      F_fldout%hgridcoreid = -1
      F_fldout%l_i0 = F_inputobj%l_i0
      F_fldout%l_j0 = F_inputobj%l_j0
      F_fldout%l_in = F_inputobj%l_in
      F_fldout%l_jn = F_inputobj%l_jn
      if (present(F_fldin)) then
         F_fldout%nkeys = F_fldin%nkeys
         F_fldout%needonce_L = F_fldin%needonce_L
         F_fldout%vn1_S = F_fldin%vn1_S
         F_fldout%vn2_S = F_fldin%vn2_S
         F_fldout%vgrid_S = F_fldin%vgrid_S
         F_fldout%sfc_S(:) = F_fldin%sfc_S(:)
         F_fldout%alt_S = F_fldin%alt_S
         F_fldout%jdatev = F_fldin%jdatev
         F_fldout%k1(:) = F_fldin%k1(:)
         F_fldout%k2(:) = F_fldin%k2(:)
         F_fldout%ks(:) = F_fldin%ks(:)
         F_fldout%kalt(:) = F_fldin%kalt(:)
         F_fldout%hgridid = F_fldin%hgridid
         F_fldout%hgridcoreid = F_fldin%hgridcoreid
         reset_L = .false.
         if (present(F_reset_L)) reset_L = F_reset_L
         if (.not.reset_L) then
            if (associated(F_fldin%sfc)) F_fldout%sfc => F_fldin%sfc
            if (associated(F_fldin%d1)) F_fldout%d1 => F_fldin%d1
            if (associated(F_fldin%d2)) F_fldout%d2 => F_fldin%d2
            if (associated(F_fldin%dalt)) F_fldout%dalt => F_fldin%dalt
            if (associated(F_fldin%psfc)) F_fldout%psfc => F_fldin%psfc
            if (associated(F_fldin%p1)) F_fldout%p1 => F_fldin%p1
            if (associated(F_fldin%p2)) F_fldout%p2 => F_fldin%p2
            if (associated(F_fldin%palt)) F_fldout%palt => F_fldin%palt
            F_fldout%hstats = F_fldin%hstats
            F_fldout%hstat = F_fldin%hstat
         endif
      endif
      !----------------------------------------------------------------------
      return
   end function priv_newfld


   !/@*
   function priv_freefld(F_fld, F_fld2) result(F_istat)
      implicit none
      type(inputio_fld_t), intent(inout):: F_fld
      type(inputio_fld_t), intent(in), optional:: F_fld2
      integer :: F_istat
      !*@/
      integer :: istat
      !----------------------------------------------------------------------
      F_istat = RMN_OK

      if (F_fld%salloc_L) then
         F_fld%salloc_L = .false.
         if (associated(F_fld%psfc)) deallocate(F_fld%psfc, stat=istat)
      endif
      nullify(F_fld%psfc, F_fld%sfc)
      
      if (F_fld%altalloc_L) then
         F_fld%altalloc_L = .false.
         if (associated(F_fld%palt)) deallocate(F_fld%palt, stat=istat)
      endif
      nullify(F_fld%palt, F_fld%dalt)

      if (F_fld%dalloc_L) then
         if (present(F_fld2)) then
            if (associated(F_fld%p1)) then
               if (.not.associated(F_fld2%p1, F_fld%p1)) then
                  F_fld%dalloc_L = .false.
                  deallocate(F_fld%p1, stat=istat)
               endif
            endif
            if (associated(F_fld%p2)) then
               if (.not.associated(F_fld2%p2, F_fld%p2)) &
                    deallocate(F_fld%p2, stat=istat)
            endif
         else
            F_fld%dalloc_L = .false.
            if (associated(F_fld%p1)) deallocate(F_fld%p1, stat=istat)
            if (associated(F_fld%p2)) deallocate(F_fld%p2, stat=istat)
         endif
      endif
      nullify(F_fld%d1, F_fld%d2, F_fld%p1, F_fld%p2)
      !----------------------------------------------------------------------
      return
   end function priv_freefld


   !/@*
   function priv_frh(F_cfg, F_ivar, F_cfgvar, F_cfgfile, F_fld, F_status, &
        F_iter, F_vgridout_S) result(F_istat)
      implicit none
      type(INCFG_T), intent(in) :: F_cfg
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FILES_T), intent(inout) :: F_cfgfile
      type(INPUTIO_FLD_T), target, intent(inout) :: F_fld
      integer, intent(in)  :: F_status  !# time_interp_status
      integer, intent(in)  :: F_ivar, F_iter
      character(len=*), intent(in) :: F_vgridout_S
      integer :: F_istat
      !*@/
      integer :: ifile, fileidx
      character(len=256) :: msg_S
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) frh [BEGIN]')
      F_istat = RMN_ERR
      do ifile = 1, size(F_cfgvar%files_S)
         if (F_cfgvar%files_S(ifile) == ' ') cycle
         fileidx = inputio_files_get_idx(F_cfgfile, F_cfgvar%files_S(ifile))
         if (.not.RMN_IS_OK(fileidx)) then
            call msg(MSG_INFO, '(inputio) Ignoring unknown file: '// &
                 trim(F_cfgvar%files_S(ifile)))
            cycle
         endif
         F_istat = priv_find(F_cfg, F_ivar, F_cfgvar, F_cfgfile, F_fld, &
              fileidx, F_status, F_iter, F_vgridout_S)
         if (RMN_IS_OK(F_istat)) then
            exit
         else
            msg_S = '(inputio) '//trim(F_fld%vn1_S)//' '//trim(F_fld%vn2_S)// &
                 ' Not found in file "'//trim(F_cfgvar%files_S(ifile))//'"'
            if (ifile < size(F_cfgvar%files_S)) then
               if (F_cfgvar%files_S(ifile+1) /= '') msg_S = trim(msg_S)// &
                    '; looking in "'//trim(F_cfgvar%files_S(ifile+1))//'"'
               call msg(MSG_INFO, msg_S)
            endif
         endif
      enddo
!!$      F_istat = min(F_istat, minval(F_fld%k1(1:F_fld%nkeys)))
      if (RMN_IS_OK(F_istat)) then
         F_istat = priv_read(F_cfg, F_cfgvar, F_cfgfile, fileidx, F_fld)
      endif
      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) frh [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_frh


   !/@*
   function priv_find(F_cfg, F_ivar, F_cfgvar, F_cfgfile, F_fld, F_fileidx, &
        F_status, F_iter, F_vgridout_S) result(F_istat)
      implicit none
      type(INCFG_T), intent(in) :: F_cfg
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FILES_T), intent(inout) :: F_cfgfile
      type(INPUTIO_FLD_T), target, intent(inout) :: F_fld
      integer, intent(in)  :: F_fileidx
      integer, intent(in)  :: F_status  !# time_interp_status
      integer, intent(in)  :: F_ivar, F_iter
      character(len=*), intent(in) :: F_vgridout_S
      integer :: F_istat
      !*@/
      integer :: ftype, tint, ip2, ip2m1, ip2p1, datevfuzz, fuzztype, cmcdatev
      integer :: istat, nip1, funit, nn, ntypvar, itypvar, nkeys
      integer(INT64) :: jdatev, jdatev0, jdatevm1, jdatevp1
      integer,target :: ip1list(NMAX_LEVELS)
      type(vgrid_descriptor) :: vgrid, vgrid0
      integer, pointer :: pip1list(:), pk1(:), pk2(:)
      character(len=32) :: dummy_S, typvar_S, vn_S, msg_S, lvl_type_S, typvarlist_S(5), alt_S
      logical :: ispressin_L, ispressout_L, usealtin_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) find [BEGIN]')
      F_istat = RMN_ERR
      pk1 => F_fld%k1(:)
      pk2 => F_fld%k2(:)

      F_fld%nkeys = 0
      F_fld%k1(:) = RMN_ERR
      F_fld%k2(:) = RMN_ERR

      funit = inputio_files_open(F_cfgfile, F_fileidx)
      ftype = inputio_files_get_type(F_cfgfile, F_fileidx)
      call priv_clim_date(ftype, F_fld%jdatev, jdatev, &
           jdatevm1, jdatevp1, ip2, ip2m1, ip2p1)
      tint = time_interp_typecode(F_cfgvar%tint_S)
      if (ftype == INPUT_FILES_GEOP) tint = TIME_INTERP_ANY

      datevfuzz = 0
      fuzztype  = FST_FIND_NEAR
      select case(F_status)
      case(TIME_INTERP_NOT_FOUND)
         fuzztype = FST_FIND_LE
         if (tint == TIME_INTERP_NEXT) fuzztype = FST_FIND_GE
      case(TIME_INTERP_NEED_PREV)
         fuzztype = FST_FIND_LE
         if (tint == TIME_INTERP_NEXT) fuzztype = FST_FIND_GE
      case(TIME_INTERP_NEED_NEXT)
         fuzztype = FST_FIND_GE
!!$         if (any(tint == (/TIME_INTERP_NEXT, TIME_INTERP_STEP/))) fuzztype = FST_FIND_GT
         if (F_iter > 2) fuzztype = FST_FIND_GT
      end select
      if (tint /= TIME_INTERP_NONE) datevfuzz = huge(datevfuzz)

      if (fuzztype == FST_FIND_LE) then
         jdatev = jdatevm1
         ip2 = ip2m1
      else
         jdatev = jdatevp1
         ip2 = ip2p1
      endif

      if (tint == TIME_INTERP_ANY) then
         jdatev = MU_JDATE_ANY
         ip2 = RMN_ANY_I
      endif
 
      jdatev0 = jdatev
      if (ftype == INPUT_FILES_CLIM) then
         jdatev = MU_JDATE_ANY
      elseif (ftype == INPUT_FILES_GEOP) then
         jdatev = MU_JDATE_ANY
         ip2 = RMN_ANY_I
      endif

      cmcdatev = jdate_to_cmc(jdatev)
      ip1list = -1
      istat = incfg_meta(F_cfg, F_ivar, dummy_S, F_ip1list=ip1list, &
           F_lvl_type_S=lvl_type_S, F_nip1=nip1, F_typvar_S=typvar_S)
      if (lvl_type_S(1:5) == 'tdiag') then
         ip1list(1) = FST_FIND_DIAG_T
         nip1 = 1
      elseif (lvl_type_S(1:5) == 'mdiag') then
         ip1list(1) = FST_FIND_DIAG_M
         nip1 = 1
      endif
      if (ip1list(1) == -1) then
         nullify(pip1list)
      else
         pip1list => ip1list(1:nip1)
      endif

      IF_VINT: if (F_cfgvar%vint_S /= 'none' .and. F_vgridout_S /= '') then
         if (F_fld%vn2_S == '') then
            F_fld%nkeys = fstmpio_find_3d_0(pk1, funit, F_fld%vn1_S, &
                 cmcdatev, pip1list, ip2, RMN_ANY_I, &
                 datevfuzz, fuzztype, F_typvar_S=F_cfgvar%typvar_S, F_vgrid=vgrid)
         else
            F_fld%nkeys = fstmpio_find_3d_vect(pk1, pk2, funit, &
                 F_fld%vn1_S, F_fld%vn2_S, &
                 cmcdatev, pip1list, ip2, RMN_ANY_I, &
                 datevfuzz, fuzztype, F_typvar_S=F_cfgvar%typvar_S, F_vgrid=vgrid)
         endif
      else  !IF_VINT
         if (F_fld%vn2_S == '') then
            F_fld%nkeys = fstmpio_find_3d_0(pk1, funit, F_fld%vn1_S, &
                 cmcdatev, pip1list, ip2, RMN_ANY_I, &
                 datevfuzz, fuzztype, F_typvar_S=F_cfgvar%typvar_S)
         else
            F_fld%nkeys = fstmpio_find_3d_vect(pk1, pk2, funit, &
                 F_fld%vn1_S, F_fld%vn2_S, &
                 cmcdatev, pip1list, ip2, RMN_ANY_I, &
                 datevfuzz, fuzztype, F_typvar_S=F_cfgvar%typvar_S)
         endif
      endif IF_VINT

      if (ftype == INPUT_FILES_CLIM .or. ftype == INPUT_FILES_GEOP) then
         F_fld%jdatev = jdatev0
      else
         F_fld%jdatev = jdate_from_cmc(cmcdatev)
      endif
      if (F_fld%nkeys > 0) F_istat = RMN_OK

!       F_fld%nkeys = size(F_fld%ip1list)
!       F_istat = maxval(F_fld%k1(1:F_fld%nkeys))
!       if (RMN_IS_OK(F_istat)) then
!          if (F_fld%ftype == INPUT_FILES_CLIM) then
!             F_fld%jdatev = jdatev
!          elseif (itype == INPUT_FILES_GEOP) then
!             F_fld%jdatev = F_fld%jdatev
!          endif
!       else
!          F_fld%nkeys = 0
!          return
!       endif

      F_fld%ks(:) = RMN_ERR
      F_fld%kalt(:) = RMN_ERR
      F_fld%vgrid_S = ''
      F_fld%sfc_S = ''
      F_fld%alt_S = ''
      IF_VINT2: if (RMN_IS_OK(F_istat) .and. F_cfgvar%vint_S /= 'none' &
           .and. F_vgridout_S /= '') then
         ispressin_L = vgrid_wb_is_press_kind(vgrid)
         ispressout_L = vgrid_wb_is_press_kind(F_vgridout_S)
         usealtin_L = (ispressin_L .and. .not.ispressout_L)
         if (usealtin_L) then
            istat = vgrid_wb_get(F_vgridout_S, vgrid0, F_altfld_S=alt_S)
            usealtin_L = (alt_S == '')
         endif
         if (usealtin_L) then
            write(dummy_S, '(i0)') jdate_to_cmc(jdatev)
            F_fld%alt_S = ADD_PREFIX3(F_fld%id, F_fileidx, trim(adjustl(dummy_S)), ALT_VGRID_H_FLDNAME)
!!$         else if (ispressout_L .and. .not.ispressin_L) then  !#TODO
!!$            write(dummy_S, '(i0)') jdate_to_cmc(jdatev)
!!$            F_fld%alt_S = ADD_PREFIX3(F_fld%id, F_fileidx, trim(adjustl(dummy_S)), ALT_VGRID_P_FLDNAME)
         else
            istat = vgd_get(vgrid, 'RFLD', F_fld%sfc_S(1))
            F_fld%sfc_S(1) = ADD_PREFIX2(F_fld%id, F_fileidx, F_fld%sfc_S(1))
            if (istat /= VGD_OK) F_fld%sfc_S(1) = ' '
            istat = vgd_get(vgrid, 'RFLS', F_fld%sfc_S(2))
            F_fld%sfc_S(2) = ADD_PREFIX2(F_fld%id, F_fileidx, F_fld%sfc_S(2))
            if (istat /= VGD_OK) F_fld%sfc_S(2) = ' '
         endif
         F_fld%vgrid_S = ADD_PREFIX2(F_fld%id, F_fileidx, F_fld%vn1_S)
         istat = vgrid_wb_get(F_fld%vgrid_S, vgrid0)
         if (.not.RMN_IS_OK(istat)) then
            istat = vgrid_wb_put(F_fld%vgrid_S, vgrid, pip1list, &
                 F_fld%sfc_S(1), F_fld%sfc_S(2), F_altfld_S=F_fld%alt_S)
         elseif (.not.(vgrid0 == vgrid)) then
            call msg(MSG_WARNING,'(inputio) vgrid, ignoring an inconsistant vgrid for '//trim(F_fld%vn1_S))
            return
         endif
         istat = vgd_free(vgrid0)
         istat = vgd_free(vgrid)

         if (ftype == INPUT_FILES_CLIM) then
            jdatev = MU_JDATE_ANY
            cmcdatev = RMN_ANY_DATE
         elseif (ftype == INPUT_FILES_GEOP) then
            ip2 = RMN_ANY_I
            jdatev = MU_JDATE_ANY
            cmcdatev = RMN_ANY_DATE
         else
            jdatev = F_fld%jdatev
         endif
         fuzztype = FST_FIND_NEAR
         datevfuzz = 0
         typvarlist_S(1) = F_cfgvar%typvar_S
         ntypvar = 1
         if (any(typvar_S == (/'r', 'R'/))) ntypvar = 0
         if (any(ftype == (/INPUT_FILES_CLIM, INPUT_FILES_GEOP/))) then
            typvarlist_S(ntypvar+1:ntypvar+4) = (/'C', 'A', 'P', 'X'/)
         else
            typvarlist_S(ntypvar+1:ntypvar+4) = (/'A', 'P', 'C', 'X'/)
         endif
         ntypvar = ntypvar + 4

         !#TODO: try to reduce code dup with next (sfc) code bloc
         istat = priv_altvcoor_status(F_fld%alt_S)
         IF_ALT: if (.not.RMN_IS_OK(istat)) then
!!$            vn_S = F_fld%alt_S
!!$            vn_S = RM_PREFIX3(vn_S)
            vn_S = ALT_VGRID_H_FLDNAME
            pk1 => F_fld%kalt(:)
            DO_TYPVAR1: do itypvar = 1, ntypvar
               nkeys = fstmpio_find_3d_0(pk1, funit, vn_S, &
                    cmcdatev, pip1list, ip2, RMN_ANY_I, &
                    datevfuzz, fuzztype, F_typvar_S=typvarlist_S(itypvar))
               if (.not.RMN_IS_OK(nkeys)) then
                  !#TODO: should we do this?, what if the alt field is date independent
                  cmcdatev = RMN_ANY_DATE
                  nkeys = fstmpio_find_3d_0(pk1, funit, vn_S, &
                       cmcdatev, pip1list, ip2, RMN_ANY_I, &
                       datevfuzz, fuzztype, F_typvar_S=typvarlist_S(itypvar))
               endif               
               if (RMN_IS_OK(nkeys)) exit DO_TYPVAR1
            enddo DO_TYPVAR1
            if (nkeys /= F_fld%nkeys .or. .not.RMN_IS_OK(nkeys)) then
               F_istat = RMN_ERR
               F_fld%nkeys = 0
               call msg(MSG_WARNING, '(inputio) Problem finding alt vgd coor field (' &
                    //trim(vn_S)//') for: '//trim(F_fld%vn1_S)//' datev='//trim(cmcdate_toprint(cmcdatev))//' typv='//typvarlist_S(itypvar))
            endif
         endif IF_ALT

         DO_SFC: do nn = 1, size(F_fld%sfc_S)
            IF_SFC: if (F_fld%sfc_S(nn) /= ' ') then
               vn_S = F_fld%sfc_S(nn)
               vn_S = RM_PREFIX2(vn_S)
               ip1list = -1
               pip1list => ip1list(1:1)
               pk1 => F_fld%ks(nn:nn)
               DO_TYPVAR: do itypvar = 1, ntypvar
                  istat = fstmpio_find_3d_0(pk1, funit, vn_S, &
                       cmcdatev, pip1list, ip2, RMN_ANY_I, &
                       datevfuzz, fuzztype, F_typvar_S=typvarlist_S(itypvar))
                  if (.not.RMN_IS_OK(istat)) then
                     !#TODO: should we do this?, what if the sfc field is date independent like ME
                     cmcdatev = RMN_ANY_DATE
                     istat = fstmpio_find_3d_0(pk1, funit, vn_S, &
                          cmcdatev, pip1list, ip2, RMN_ANY_I, &
                          datevfuzz, fuzztype, F_typvar_S=typvarlist_S(itypvar))
                  endif
                  if (RMN_IS_OK(istat)) exit DO_TYPVAR
               enddo DO_TYPVAR
               if (.not.RMN_IS_OK(istat)) then
                  F_istat = RMN_ERR
                  F_fld%nkeys = 0
                  call msg(MSG_WARNING, '(inputio) Problem finding sfc ref field (' &
                       //trim(RM_PREFIX2(F_fld%sfc_S(nn)))//') for: '// &
                       trim(RM_PREFIX2(F_fld%vn1_S)))
               endif
            endif IF_SFC
         enddo DO_SFC

      endif IF_VINT2

      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) find [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_find


   !/@*
   function priv_read(F_cfg, F_cfgvar, F_cfgfile, F_fileidx, F_fld) result(F_istat)
      implicit none
      type(INCFG_T), intent(in) :: F_cfg
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FILES_T), intent(inout) :: F_cfgfile
      type(INPUTIO_FLD_T), target, intent(inout) :: F_fld
      integer, intent(inout) :: F_fileidx
      integer :: F_istat
      !*@/
      real,parameter :: MB2PA = 100.
      character(len=1024) :: msg_S, tmp_S, vn_S
      character(len=8), target :: hints_S(1)
      character(len=8), pointer :: phints_S(:)
      integer :: istat, ivar
      integer, target :: funit(1), hstats(2), hstats2(NMAX_LEVELS)
      integer, pointer :: pfunit(:), phstats(:), pk1(:), pk2(:)
      logical :: isassoc_L, isassoc2_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) read [BEGIN]')
      F_istat = RMN_ERR

      hints_S(1) = F_cfgvar%hint_S
      phints_S => hints_S(1:1)
      funit(1) = inputio_files_open(F_cfgfile, F_fileidx)
      pfunit => funit(1:1)
      phstats => F_fld%hstats(:)
      pk1 => F_fld%k1(1:F_fld%nkeys)
      pk2 => F_fld%k2(1:F_fld%nkeys)

      isassoc_L = associated(F_fld%p1)
      F_fld%hstats = RMN_ERR
      IF_BLOCIO: if (F_cfgfile%iotype == PTOPO_BLOC) then
         if (F_fld%vn2_S == ' ') then
            F_istat = fstmpi_rdhint_3d_r4( &
                 F_fld%p1, phstats, pk1, phints_S, pfunit, &
                 F_fld%hgridid, F_fld%hgridcoreid)
         else
            F_istat = fstmpi_rdhint_3d_r4_vect( &
                 F_fld%p1, F_fld%p2, phstats, pk1, pk2, phints_S, pfunit, &
                 F_fld%hgridid, F_fld%hgridcoreid)
         endif
      else !IF_BLOCIO
         if (F_fld%vn2_S == ' ') then
            F_istat = fstmpio_rdhint_3d_r4( &
                 F_fld%p1, phstats, pk1, phints_S, pfunit, &
                 F_cfg%commgid, F_fld%hgridid, F_fld%hgridcoreid)
         else
            F_istat = fstmpio_rdhint_3d_r4_vect( &
                 F_fld%p1, F_fld%p2, phstats, pk1, pk2, phints_S, pfunit, &
                 F_cfg%commgid, F_fld%hgridid, F_fld%hgridcoreid)
         endif
      endif IF_BLOCIO
      F_fld%hstat = maxval(F_fld%hstats(1:F_fld%nkeys))
      F_istat = min(F_istat, minval(F_fld%hstats(1:F_fld%nkeys)))

      !#TODO: level by level (hstats) warning...
      call str_concat_i(tmp_S, pk1, ', ')
      write(msg_S, '(i4,1x,a)') size(pk1), ' levels [ip1='//trim(tmp_S)//']'
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFO, '(inputio) Read ' &
              //trim(F_fld%vn1_S)//' '//trim(F_fld%vn2_S)//': '//msg_S)
         if (.not.isassoc_L) F_fld%dalloc_L = .true.
      else
         call msg(MSG_WARNING, '(inputio) Problem reading: ' &
              //trim(F_fld%vn1_S)//' '//trim(F_fld%vn2_S)//': '//msg_S)
         return
      endif

!!$         nullify(F_fld%palt)
      isassoc2_L = associated(F_fld%palt)
      READ_ALT: if (F_fld%alt_S /= '' .and. F_fld%kalt(1) >= 0) then
         pk1 => F_fld%kalt(1:F_fld%nkeys)
         phstats => hstats2(:)
         phints_S(1) = 'cubic'
         IF_BLOCIO1: if (F_cfgfile%iotype == PTOPO_BLOC) then
            istat = fstmpi_rdhint_3d_r4( &
                 F_fld%palt, phstats, pk1, phints_S, pfunit, &
                 F_fld%hgridid, F_fld%hgridcoreid)
         else
            istat = fstmpio_rdhint_3d_r4( &
                 F_fld%palt, phstats, pk1, phints_S, pfunit, &
                 F_cfg%commgid, F_fld%hgridid, F_fld%hgridcoreid)
         endif IF_BLOCIO1
         !#TODO: review hstat and error checking
!!$         hstat = maxval(hstats2(1:F_fld%nkeys))
         if (RMN_IS_OK(istat)) then
!!$            F_fld%hstat = max(F_fld%hstat, maxval(phstats))
            if (.not.isassoc2_L) F_fld%altalloc_L = .true.
            if (associated(F_fld%palt)) then
!!$               vn_S = RM_PREFIX3(F_fld%alt_S)
               vn_S = ALT_VGRID_H_FLDNAME
               if (any(vn_S == (/'P0','p0'/))) then
                  F_fld%palt(:,:,:) = F_fld%palt(:,:,:) * MB2PA
                  !#TODO: does GZ need unit conversion?
               endif
            endif
         else
            F_istat = RMN_ERR
            call msg(MSG_WARNING, '(inputio) Problem reading alt vgd coor fld: ' &
                 //trim(F_fld%alt_S))
            return
         endif
         
      endif READ_ALT

      isassoc2_L = associated(F_fld%psfc)
      READ_SFC: if (F_fld%sfc_S(1) /= '' .and. F_fld%ks(1) >= 0) then
         pk1 => F_fld%ks(1:2)
         phstats => hstats(:)
         phints_S(1) = 'cubic'
         IF_BLOCIO2: if (F_cfgfile%iotype == PTOPO_BLOC) then
            istat = fstmpi_rdhint_3d_r4( &
                 F_fld%psfc, phstats, pk1, phints_S, pfunit, &
                 F_fld%hgridid, F_fld%hgridcoreid)
         else
            istat = fstmpio_rdhint_3d_r4( &
                 F_fld%psfc, phstats, pk1, phints_S, pfunit, &
                 F_cfg%commgid, F_fld%hgridid, F_fld%hgridcoreid)
         endif IF_BLOCIO2
         if (RMN_IS_OK(istat)) then
            F_fld%hstat = max(F_fld%hstat, maxval(hstats))
            if (.not.isassoc2_L) F_fld%salloc_L = .true.
            if (associated(F_fld%psfc)) then
               do ivar = 1, size(F_fld%sfc_S)
                  vn_S = RM_PREFIX2(F_fld%sfc_S(ivar))
                  if (any(vn_S == (/'P0  ','p0  ','P0LS','p0ls'/))) then
                     F_fld%psfc(:,:,ivar) = F_fld%psfc(:,:,ivar) * MB2PA
                  endif
               enddo
            endif
         else
            F_istat = RMN_ERR
            call msg(MSG_WARNING, '(inputio) Problem reading sfc ref fld: ' &
                 //trim(F_fld%sfc_S(1))//' '//trim(F_fld%sfc_S(2)))
            return
         endif
      endif READ_SFC

      if (RMN_IS_OK(F_istat)) call priv_set_scope(F_fld)

      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) read [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_read


   !/@*
   subroutine priv_clim_date(F_filetype, F_jdatev0, F_jdatev, F_jdatevm1, &
        F_jdatevp1, F_ip2, F_ip2m1, F_ip2p1)
      implicit none
      !@objective
      !@arguments
      integer, intent(in) :: F_filetype
      integer(INT64), intent(in) :: F_jdatev0
      integer(INT64), intent(out) :: F_jdatev, F_jdatevm1, F_jdatevp1
      integer, intent(out) :: F_ip2, F_ip2m1, F_ip2p1
      !*@/
      integer :: year
      !------------------------------------------------------------------
      F_jdatev = F_jdatev0
      F_jdatevm1 = F_jdatev0
      F_jdatevp1 = F_jdatev0
      F_ip2 = RMN_ANY_I ; F_ip2m1 = RMN_ANY_I ; F_ip2p1 = RMN_ANY_I

      if (F_filetype == INPUT_FILES_ANAL) return

      F_jdatev = MU_JDATE_ANY
      F_jdatevm1 = MU_JDATE_ANY
      F_jdatevp1 = MU_JDATE_ANY

      if (F_filetype /= INPUT_FILES_CLIM) return

      F_ip2 = jdate_month(F_jdatev0)
      F_ip2m1 = F_ip2
      year = jdate_year(F_jdatev0)
      F_jdatevm1 = jdate_midmonth(year, F_ip2m1)
      if (F_jdatevm1 > F_jdatev0) then
         F_ip2m1 = F_ip2 - 1
         if (F_ip2m1 < 1) then
            F_ip2m1 = F_ip2m1 + 12
            year = year - 1
         endif
         F_jdatevm1 = jdate_midmonth(year, F_ip2m1)
      endif
      F_ip2p1 = F_ip2m1 + 1
      if (F_ip2p1 > 12) then
         F_ip2p1 = F_ip2p1 - 12
         year = year + 1
      endif
      F_jdatevp1 = jdate_midmonth(year, F_ip2p1)
      !------------------------------------------------------------------
      return
   end subroutine priv_clim_date


   !/@*
   function priv_tint_status(F_cfgvar, F_jdatev, F_fld) result(F_istat)
      implicit none
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      integer(INT64), intent(in) :: F_jdatev
      type(inputio_fld_t), intent(inout) :: F_fld
      integer :: F_istat
      !*@/
      character(len=64) :: vn_S, msg_S
      integer :: tint
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) tint_status [BEGIN]')
      F_istat = RMN_ERR

      IF_NONE: if (any(F_cfgvar%tint_S(1:4) == (/'none', 'any '/))) then
         !#TODO: what if (F_fld%nkeys == 0)
         if (F_fld%nkeys == 0) then
            print *,'WARNING: (inputio) priv_tint_status nkeys==0'
            call flush(6)
         endif
         if (F_fld%nkeys > 0) F_istat = RMN_OK
         if (F_fld%nkeys < 0) then
            F_istat = TIME_INTERP_NOT_FOUND
            F_fld%nkeys = 0
         endif
      else
         tint = time_interp_typecode(F_cfgvar%tint_S)
         vn_S = ADD_PREFIX(F_fld%id, F_fld%vn1_S)
         F_istat = time_interp_status(vn_S, F_jdatev, tint)
      endif IF_NONE
      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) tint_status [END] '//trim(msg_S))
       !------------------------------------------------------------------
      return
   end function priv_tint_status


   !/@*
   function priv_tint_set(F_cfgvar, F_fld) result(F_istat)
      implicit none
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FLD_T), intent(in) :: F_fld
      integer :: F_istat
      !*@/
      integer :: istat, tint
      logical :: inrestart_L
      character(len=64) :: vn_S, msg_S
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) tint_set [BEGIN]')
      F_istat = RMN_ERR

      IF_NONE: if (any(F_cfgvar%tint_S(1:4) == (/'none', 'any '/))) then
         if (F_fld%nkeys > 0) F_istat = RMN_OK
      else

         inrestart_L = .not.F_fld%needonce_L

         vn_S = ADD_PREFIX(F_fld%id, F_fld%vn1_S)
         F_istat = time_interp_set(F_fld%d1, vn_S, F_fld%jdatev, F_fld%vgrid_S, &
              F_fld%sfc_S(1), inrestart_L, F_flag=F_fld%hstat, &
              F_sfcfld2_S=F_fld%sfc_S(2))
         if (RMN_IS_OK(F_istat) .and. F_fld%vn2_S /= ' ') then
            vn_S = ADD_PREFIX(F_fld%id, F_fld%vn2_S)
            F_istat = time_interp_set(F_fld%d2, vn_S, F_fld%jdatev, F_fld%vgrid_S, &
                 F_fld%sfc_S(1), inrestart_L, F_flag=F_fld%hstat, &
                 F_sfcfld2_S=F_fld%sfc_S(2))
         endif

         if (RMN_IS_OK(F_istat) .and. F_fld%sfc_S(1) /= ' ') then
            vn_S = RM_PREFIX2(F_fld%sfc_S(2))
            vn_S = trim(F_fld%sfc_S(1))//trim(vn_S)
            tint = time_interp_typecode(F_cfgvar%tint_S)
            istat = time_interp_status(vn_S, F_fld%jdatev, tint)
            if (.not.RMN_IS_OK(istat)) then
               F_istat = time_interp_set(F_fld%sfc, vn_S, F_fld%jdatev, &
                    inrestart_L, F_flag=F_fld%hstat)
            endif
         endif

      endif IF_NONE

      if (RMN_IS_OK(F_istat) .and. F_fld%alt_S /= ' ') &
           F_istat = priv_altvcoor_set(F_fld%alt_S, F_fld%palt)
!!$      F_istat = priv_altvcoor_set(F_fld%alt_S, F_fld%dalt)

      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) tint_set [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_tint_set


   !/@*
   function priv_tint_get(F_cfgvar, F_fld, F_next_prev) result(F_istat)
      implicit none
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FLD_T), intent(inout) :: F_fld
      integer, intent(in) :: F_next_prev
      integer :: F_istat
      !*@/
      character(len=64) :: vn_S, msg_S
      integer(INT64) :: jdatev
      integer :: istat
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) tint_get [BEGIN]')
      F_istat = RMN_ERR

      IF_NONE: if (any(F_cfgvar%tint_S(1:4) == (/'none', 'any '/))) then
         if (F_fld%nkeys > 0 .and. associated(F_fld%d1) .and. &
              (F_fld%vn2_S == '' .or. associated(F_fld%d2))) F_istat = RMN_OK
      else

         F_fld%dalloc_L = .true.
         F_fld%salloc_L = .true.
         istat = priv_freefld(F_fld)

         vn_S = ADD_PREFIX(F_fld%id, F_fld%vn1_S)
         F_istat = time_interp_retrieve(vn_S, F_next_prev, F_fld%jdatev, &
              F_data=F_fld%p1, F_vgrid_S=F_fld%vgrid_S, F_sfcfld_S=F_fld%sfc_S(1), &
              F_sfcfld2_S=F_fld%sfc_S(2), F_flag=F_fld%hstat)
         if (associated(F_fld%p1)) F_fld%d1 => F_fld%p1
         if (RMN_IS_OK(F_istat) .and. F_fld%vn2_S /= ' ') then
            !#TODO: if (associated(F_fld%p2) .and. F_fld%dalloc_L) deallocate(F_fld%p2)
            vn_S = ADD_PREFIX(F_fld%id, F_fld%vn2_S)
            F_istat = time_interp_retrieve(vn_S, F_next_prev, jdatev, &
                 F_data=F_fld%p2)
            if (associated(F_fld%p2)) F_fld%d2 => F_fld%p2
        endif
         F_fld%dalloc_L = .false.

         if (RMN_IS_OK(F_istat) .and. F_fld%sfc_S(1) /= ' ') then
            vn_S = RM_PREFIX2(F_fld%sfc_S(2))
            vn_S = trim(F_fld%sfc_S(1))//trim(vn_S)
            F_istat = time_interp_retrieve(vn_S, F_next_prev, jdatev, &
                 F_data=F_fld%psfc)
            F_fld%salloc_L = .false.
            if (associated(F_fld%psfc)) F_fld%sfc => F_fld%psfc
         end if

      endif IF_NONE
      !#TODO:
!!$      if (.not.RMN_IS_OK(F_istat)) then
!!$         nullify(F_fld%d1, F_fld%d2)
!!$      endif

      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) tint_get [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_tint_get


   !/@*
   function priv_vint(F_cfgvar, F_fldin, F_fldout, F_vgrid_S) result(F_istat)
      implicit none
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FLD_T), intent(in) :: F_fldin
      type(INPUTIO_FLD_T), intent(inout) :: F_fldout
      character(len=*),intent(in) :: F_vgrid_S
      integer :: F_istat
      !*@/
      integer :: ivar, nlinbot, istat, lijk(3), uijk(3), lijk2(3), uijk2(3), nk
      logical :: same_sfc_L, reallocated_L, ispressin_L, ispressout_L
      character(len=256) :: vn_S, msg_S
      real, pointer, dimension(:,:,:) :: din, dout
      real, pointer, dimension(:,:) :: sfcin, sfcout, slsin, slsout
      integer,pointer :: ip1list(:)
      type(vgrid_descriptor) :: vgrid
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) vint [BEGIN]')
      F_istat = RMN_ERR
      if (.not.(associated(F_fldin%d1) .and. &
           (F_fldin%vn2_S == ' ' .or. associated(F_fldin%d2)))) then
         write(msg_S, '(i0)') F_istat
         call msg(MSG_DEBUG, '(inputio) vint [END] for "'//trim(F_fldin%vn1_S)// &
              '", No associated ptr '//trim(msg_S))
         return
      endif
      F_fldout%vgrid_S = F_vgrid_S

      IF_VINT: if (F_cfgvar%vint_S == 'none') then

         F_istat = RMN_OK

         !TODO: if associated(F_fldout%d1) check bounds and copy
         if (associated(F_fldin%d1)) F_fldout%d1 => F_fldin%d1
         if (associated(F_fldin%d2)) F_fldout%d2 => F_fldin%d2
         if (associated(F_fldin%p1)) F_fldout%p1 => F_fldin%p1
         if (associated(F_fldin%p2)) F_fldout%p2 => F_fldin%p2

!!$         do ivar = 1,2
!!$            vn_S = F_fldin%vn1_S
!!$            din => F_fldin%d1
!!$            dout => F_fldout%d1
!!$            if (ivar == 2) then
!!$               vn_S = F_fldin%vn2_S
!!$               din => F_fldin%d2
!!$               dout => F_fldout%d2
!!$            endif
!!$            if (vn_S /= ' ') then
!!$               if (all(shape(din) == shape(dout))) then
!!$                  dout(:,:,:) = din(:,:,:)
!!$               else
!!$                  F_istat = RMN_ERR
!!$               endif
!!$            endif
!!$         enddo

      else !IF_VINT

         lijk = lbound(F_fldin%d1)
         uijk = ubound(F_fldin%d1)

         !#TODO: why next 4 lines
!!$         nullify(ip1list)
!!$         istat = vgrid_wb_get(F_fldin%vgrid_S, vgrid, ip1list)
!!$         istat = vgd_free(vgrid)
!!$         if (associated(ip1list)) deallocate(ip1list, stat=istat)

         nullify(ip1list)
         istat = vgrid_wb_get(F_vgrid_S, vgrid, ip1list)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_DEBUG, '(inputio) vint [END] for "'//trim(F_fldin%vn1_S)// &
                 '", Problem getting vgrid: '//trim(F_vgrid_S))
            return
         endif
         istat = vgd_free(vgrid)
         lijk(3) = 1
         uijk(3) = size(ip1list)
         if (associated(ip1list)) deallocate(ip1list, stat=istat)

         reallocated_L = .false.
         if (.not.associated(F_fldout%p1)) then
            F_fldout%dalloc_L = .true.
            allocate(F_fldout%p1(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3)))
            F_fldout%p1 = 0.
            reallocated_L = .true.
            F_fldout%d1 => F_fldout%p1
         endif
         if (F_fldout%vn2_S /= ' ' .and. .not.associated(F_fldout%p2)) then
            allocate(F_fldout%p2(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3)))
            F_fldout%p2 = 0.
            reallocated_L = .true.
            F_fldout%d2 => F_fldout%p2
         endif

         nk = uijk(3)
         lijk = lbound(F_fldin%d1)
         uijk = ubound(F_fldin%d1)
         uijk(3) = nk
         lijk2 = lbound(F_fldout%d1)
         uijk2 = ubound(F_fldout%d1)
         if (.not.(all(lijk == lijk2) .and. all(uijk == uijk2))) then
            call msg(MSG_ERROR, '(inputio) vint, for "'//trim(F_fldin%vn1_S)// &
                 '", provided array 1 has wrong shape')
            print *,trim(F_fldin%vn1_S),lijk,":",uijk," != ",lijk2,":",uijk2
            return
         endif
         if (F_fldout%vn2_S /= ' ') then
            lijk2 = lbound(F_fldout%d2)
            uijk2 = ubound(F_fldout%d2)
            if (.not.(all(lijk == lijk2) .and. all(uijk == uijk2))) then
               call msg(MSG_ERROR, '(inputio) vint, for "'//trim(F_fldin%vn1_S)// &
                    '", provided array 2 has wrong shape')
               print *,trim(F_fldin%vn1_S),lijk,":",uijk," != ",lijk2,":",uijk2
               return
            endif
         endif


         nlinbot = 0
         same_sfc_L = .false.
         if (any(F_cfgvar%vint_S(1:4) == (/'l-co','line'/))) &
              nlinbot = size(F_fldout%d1,3)
         if (any(F_cfgvar%vint_S(1:4) == (/'l-co','c-co'/)) .and. &
              F_fldin%hstat == HINTERP4YY_NONE) then
            ispressin_L = vgrid_wb_is_press_kind(F_fldin%vgrid_S)
            ispressout_L = vgrid_wb_is_press_kind(F_fldout%vgrid_S)
            same_sfc_L = (ispressin_L .eqv. ispressout_L)
         endif
         nullify(sfcin, sfcout, slsin, slsout)
         if (associated(F_fldin%sfc)) then
            sfcin => F_fldin%sfc(:,:,1)
            if (F_fldin%sfc_S(2) /= '') slsin => F_fldin%sfc(:,:,2)
         endif

!!$         if (present(F_same_sfc_L)) then
!!$            if (F_same_sfc_L) then
!!$               if (associated(F_fldin%sfc)) then
!!$                  if (F_fldin%sfc_S(1) /= '') sfcout => F_fldin%sfc(:,:,1)
!!$                  if (F_fldin%sfc_S(2) /= '') slsout => F_fldin%sfc(:,:,2)
!!$               endif
!!$            endif
!!$         endif

         F_istat = RMN_OK
         do ivar = 1,2
            vn_S = F_fldin%vn1_S
            din => F_fldin%d1
            dout => F_fldout%d1
            if (ivar == 2) then
               vn_S = F_fldin%vn2_S
               din => F_fldin%d2
               dout => F_fldout%d2
            endif
            if (vn_S /= ' ') then
               if (associated(sfcout)) then
                  istat = vinterp( &
                       dout, F_fldout%vgrid_S,  &
                       din, F_fldin%vgrid_S, &
                       sfcout, sfcin,  &
                       nlinbot, vn_S, same_sfc_L, &
                       slsout, slsin)
                  F_istat = min(istat, F_istat)
               else
                  istat = vinterp( &
                       dout, F_fldout%vgrid_S,  &
                       din, F_fldin%vgrid_S, &
                       F_sfcfldin=sfcin,  &
                       F_nlinbot=nlinbot, F_msg_S=vn_S,  &
                       F_use_same_sfcfld_L=same_sfc_L, &
                       F_sfcfldin2=slsin)
                  F_istat = min(istat, F_istat)
               endif
            endif
         enddo

      endif IF_VINT
      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) vint [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_vint


   !/@*
   function priv_tint(F_cfgvar, F_fldin1, F_fldin2, F_fldout) result(F_istat)
      implicit none
      type(INCFG_VAR_T), intent(in) :: F_cfgvar
      type(INPUTIO_FLD_T), intent(in) :: F_fldin1
      type(INPUTIO_FLD_T), intent(in) :: F_fldin2
      type(INPUTIO_FLD_T), intent(inout) :: F_fldout
      integer :: F_istat
      !*@/
      logical :: was_assoc_L
      integer :: tint
      character(len=256) :: msg_S, tmp_S
      !------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio) tint [BEGIN]')
      F_istat = RMN_ERR

      IF_TINT: if (any(F_cfgvar%tint_S(1:4) == (/'none', 'any '/))) then

!!$      if (F_fldin1%nkeys > 0) then
         !TODO: if associated(F_fldout%d1) check bounds and copy
         if (associated(F_fldin1%d1)) F_fldout%d1 => F_fldin1%d1
         if (associated(F_fldin1%d2)) F_fldout%d2 => F_fldin1%d2
         if (associated(F_fldin1%p1)) F_fldout%p1 => F_fldin1%p1
         if (associated(F_fldin1%p2)) F_fldout%p2 => F_fldin1%p2
!!$      endif
         F_istat = RMN_ERR
         if (associated(F_fldin1%d1) .and. &
              (F_fldout%vn2_S == '' .or. associated(F_fldin1%d2))) F_istat = RMN_OK

      else

         was_assoc_L = associated(F_fldout%d1)

         tint = time_interp_typecode(F_cfgvar%tint_S)
         F_istat = time_interp_get(F_fldout%d1, F_fldin1%d1, F_fldin2%d1, &
              F_fldout%jdatev, F_fldin1%jdatev, F_fldin2%jdatev, &
              tint, F_varname_S=F_fldin1%vn1_S)
         if (F_fldin1%vn2_S /= ' ' .and. RMN_IS_OK(F_istat)) then
            F_istat = time_interp_get(F_fldout%d2, F_fldin1%d2, F_fldin2%d2, &
                 F_fldout%jdatev, F_fldin1%jdatev, F_fldin2%jdatev, &
                 tint, F_varname_S=F_fldin1%vn2_S)
         endif

         F_fldout%dalloc_L = (.not.was_assoc_L) .and. associated(F_fldout%d1)

         if (RMN_IS_OK(F_istat)) then
            if (F_istat > 0 .and. F_istat < nint(TIME_INTERP_WEIGHT_FACT)) then
               write(tmp_S,'(i6)') (F_istat*100)/nint(TIME_INTERP_WEIGHT_FACT)
               call msg(MSG_INFO, '(inputio) Got time interpolated value for '//&
                    trim(F_fldout%vn1_S)//' '//trim(F_fldout%vn2_S)// &
                    ' (interpolation h/v/t: '//trim(F_cfgvar%hint_S)//'/'// &
                    trim(F_cfgvar%vint_S)//'/'//trim(F_cfgvar%tint_S)//') (tv='// &
                    trim(F_cfgvar%typvar_S)//') (t_int_weight='//trim(tmp_S)//'%)')
            endif
         else
            call msg(MSG_WARNING, &
                 '(inputio) Problem getting time interpolated value for '//&
                 trim(F_fldout%vn1_S)//' '//trim(F_fldout%vn2_S)// &
                 ' (interpolation h/v/t: '//trim(F_cfgvar%hint_S)//'/'// &
                 trim(F_cfgvar%vint_S)//'/'//trim(F_cfgvar%tint_S)//') (tv='// &
                 trim(F_cfgvar%typvar_S)//')')
         endif

      endif IF_TINT

      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(inputio) tint [END] '//trim(msg_S))
      !------------------------------------------------------------------
      return
   end function priv_tint


   !/@*
   subroutine priv_set_scope(F_fld)
      implicit none
      type(INPUTIO_FLD_T), intent(inout) :: F_fld
      !*@/
      integer :: lijk(3), uijk(3)
      !------------------------------------------------------------------
      !#TODO? nullify(F_fld%d1, F_fld%d2, F_fld%sfc, F_fld%dalt)
      if (associated(F_fld%p1)) then
         lijk = lbound(F_fld%p1)
         uijk = ubound(F_fld%p1)
         lijk(1) = max(F_fld%l_i0, lijk(1))
         lijk(2) = max(F_fld%l_j0, lijk(2))
         uijk(1) = min(F_fld%l_in, uijk(1))
         uijk(2) = min(F_fld%l_jn, uijk(2))
         F_fld%d1 => F_fld%p1(lijk(1):uijk(1),lijk(2):uijk(2),:)
      endif
      if (associated(F_fld%p2)) then
         lijk = lbound(F_fld%p2)
         uijk = ubound(F_fld%p2)
         lijk(1) = max(F_fld%l_i0, lijk(1))
         lijk(2) = max(F_fld%l_j0, lijk(2))
         uijk(1) = min(F_fld%l_in, uijk(1))
         uijk(2) = min(F_fld%l_jn, uijk(2))
         F_fld%d2 => F_fld%p2(lijk(1):uijk(1),lijk(2):uijk(2),:)
      endif
      if (associated(F_fld%psfc)) then
         lijk = lbound(F_fld%psfc)
         uijk = ubound(F_fld%psfc)
         lijk(1) = max(F_fld%l_i0, lijk(1))
         lijk(2) = max(F_fld%l_j0, lijk(2))
         uijk(1) = min(F_fld%l_in, uijk(1))
         uijk(2) = min(F_fld%l_jn, uijk(2))
         F_fld%sfc => F_fld%psfc(lijk(1):uijk(1),lijk(2):uijk(2),:)
      endif
      if (associated(F_fld%palt)) then
         lijk = lbound(F_fld%palt)
         uijk = ubound(F_fld%palt)
         lijk(1) = max(F_fld%l_i0, lijk(1))
         lijk(2) = max(F_fld%l_j0, lijk(2))
         uijk(1) = min(F_fld%l_in, uijk(1))
         uijk(2) = min(F_fld%l_jn, uijk(2))
         F_fld%dalt => F_fld%palt(lijk(1):uijk(1),lijk(2):uijk(2),:)
      endif
      !------------------------------------------------------------------
      return
   end subroutine priv_set_scope

   
   !/@*
   function priv_altvcoor_status(F_name_S) result(F_istat)
      implicit none
      character(len=*), intent(in) :: F_name_S
      integer :: F_istat
      !*@/
      real, pointer :: data3d(:,:,:)
      type(gmm_metadata) :: meta3d
      !------------------------------------------------------------------
!!$         if (F_fld%alt_S /= ' ') istat = gmm_getmeta(F_name_S, meta3d)  !getmeta is too verbose
      F_istat = RMN_OK
      if (F_name_S /= ' ') F_istat = gmm_get(F_name_S, data3d, meta3d)
      !------------------------------------------------------------------
      return
   end function priv_altvcoor_status

   
   !/@*
   function priv_altvcoor_set(F_name_S, F_data3d) result(F_istat)
      implicit none
      character(len=*), intent(in) :: F_name_S
      real, pointer :: F_data3d(:,:,:)
      integer :: F_istat
      !*@/
      integer :: lijk(3), uijk(3), lni, lnj, lnk
      type(gmm_metadata) :: meta3d
      real, pointer :: data3d(:,:,:)
      !------------------------------------------------------------------
      F_istat = RMN_OK
      if (F_name_S == ' ' .or. .not.associated(F_data3d)) return
      F_istat = priv_altvcoor_status(F_name_S)
      if (RMN_IS_OK(F_istat)) return
      
      !#TODO: only keep the 2 time frames used for time interpolation
      lijk = lbound(F_data3d)
      uijk = ubound(F_data3d)
      lni = uijk(1) - lijk(1) + 1
      lnj = uijk(2) - lijk(2) + 1
      lnk = uijk(3) - lijk(3) + 1
      call gmm_build_meta3D(meta3d, &
           lijk(1), uijk(1), 0, 0, lni, &
           lijk(2), uijk(2), 0, 0, lnj, &
           lijk(3), uijk(3), 0, 0, lnk, &
           0, GMM_NULL_FLAGS)
      nullify(data3d)
      F_istat = gmm_create(F_name_S, data3d, meta3d, GMM_FLAG_RSTR)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING, '(inputio) Problem saving in gmm the alt vgd coor field: ' &
              //trim(F_name_S))
         return
      endif
      data3d = F_data3d
      !------------------------------------------------------------------
      return
   end function priv_altvcoor_set
   
end module inputio_mod

