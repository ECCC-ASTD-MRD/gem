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

#include <rmn/msg.h>

!/@*
module incfg2_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower, clib_isfile, clib_isreadok
   use str_mod, only: str_tab2space, str_toreal
   use mu_jdate_mod
   implicit none
   private
   !@objective
   !@author Stephane Chamberland, 2017-08
   !@description
   ! Public functions
   public :: incfg_new, incfg_add, incfg_add_string, &
        incfg_add_kv, incfg_nbvar, incfg_meta, incfg_set_meta, &
        incfg_isvarstep, incfg_varindex, &
        incfg_time, incfg_sethgridid, &
        incfg_gethgridid
   ! Public constants
   integer, parameter, public :: INCFG_STRLEN = 32
   character(len=INCFG_STRLEN), parameter, public :: INCFG_KNOWN_H_INT(3) = (/ &
        'cubic           ', 'linear          ', 'nearest         ' &
        /) !TODO-later: INCFG_KNOWN_H_INT = none or crop
   character(len=INCFG_STRLEN), parameter, public :: INCFG_KNOWN_V_INT(5) = (/ &
        'none            ', 'linear          ', 'cubic           ', &
        'l-cond          ', 'c-cond          ' &
        /) !TODO-later: INCFG_KNOWN_V_INT = nearest
   character(len=INCFG_STRLEN), parameter, public :: INCFG_KNOWN_LVL_TYPE(5) = (/ &
        'arbitrary       ', 'momentum        ', 'thermodynamic   ', &
        'tdiagnostic     ', 'mdiagnostic     ' &
        /)
   integer, parameter, public :: INCFG_LVL_ARBITRARY = 1
   integer, parameter, public :: INCFG_LVL_MOM = 2
   integer, parameter, public :: INCFG_LVL_THERMO = 3
   integer, parameter, public :: INCFG_LVL_TDIAG = 4
   integer, parameter, public :: INCFG_LVL_MDIAG = 5

   character(len=INCFG_STRLEN), parameter, public :: INCFG_KNOWN_T_INT(7) = (/ &
        'none            ', 'linear          ', 'nearest         ', &
        'any             ', 'step            ', 'next            ', &
        'increment       '/)
   ! Public Type
   public :: INCFG_T, INCFG_VAR_T

   ! File format
   !   one input var desc per line
   !   empty lines and lines starting with # are ignored
   !   a line is multiple key=var separated with ;
   !   spaces are ignored
   !   case is ignored
   !   recognized keys
   !   * in     : variable name as found in file (mandatory)
   !   * in2    : variable name as found in file, second component of vector fields
   !   * search : list of files to search var in, comma separated (mandatory)
   !              accepted values: GEOP, CLIM, ANAL, INREP,
   !                               MANAL, MINREP, MINPUT
   !              var are search in until found in files starting from 1st one
   !              Note: from GEM 4.8 onward,
   !                     MANAL  could be used instead of ANAL (TAKS_INPUT/MODEL_ANALYSIS)
   !                     MINREP should be used instead of INREP (TAKS_INPUT/MODEL_INREP)
   !                     MINPUT could be used to point to TAKS_INPUT/MODEL_INPUT
   !   * freq   : reading frequency (defaut=0)
   !              accepted format: TYPE,first,interval,last
   !              TYPE: STEP,MONTH,HOUR,MINUTE (optional, default=STEP)
   !              first:    start reading from "first", TYPE units (mandatory)
   !              interval: read every interval (optional, default=0, inital)
   !              last :    read up to "last", TYPE units (optional, default=-1 infinite)
   !   * hinterp: spatial horiz. interp. technique (defaut=cubic)
   !              accepted values: nearest, linear, cubic
   !   * vinterp: spatial vertical interp. technique (defaut=none)
   !              accepted values:
   !                none,
   !                linear,
   !                cubic,
   !                l-cond (linear conditional, interpolate only if different vgrid or hgrid)
   !                c-cond (cubic  conditional, interpolate only if different vgrid or hgrid)
   !   * tinterp: temporal interp. technique (defaut=none)
   !              accepted values: nearest, linear, step, any, none
   !   * levels : list of levels (defaul=0)
   !              accepted format: FIRST_LEVEL, LAST_LEVEL
   !              LAST_LEVEL is optional
   !              accepted values:
   !                0     (surface),
   !                1     (arbitrary level 1),
   !                1,26  (arbitrary level 1 to 26),
   !                -1    (caller provided list of levels),
   !                tdiag (thermo diag level),
   !                mdiag (momentum diag level)
   !   * cat    : list of categories (for 4D vars) (default=-1)
   !              accepted format: FIRST_CAT, LAST_CAT
   !              LAST_CAT is optional
   !              accepted values:
   !                1     (1st category),
   !                1,26  (categories 1 to 26),
   !                -1    (caller provided list of categories, all categories),
   !   * typvar : typvar param of RPN-std files (default=' ')
   !              accepted values: C, A, P, R, ...
   !   * mandatory : define if the var is mandatory or not (defaut=default)
   !                 accepted values:
   !                 default (up to the caller to decide), true, false
   !   * vmin   : min value of the field... f = max(vmin,f)
   !   * vmax   : max value of the field... f = min(f,vmax)
   !  key "interp" (now hinterp) is still accepted for backward compatibility
   !  key "timeint"(now tinterp) is still accepted for backward compatibility
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   real, parameter :: EPSILON = 1.e-5
   integer, parameter :: NFILE_MAX = 4
   integer, parameter :: NVAR_MAX = 512

!!$   type :: INCFG_T
!!$      logical :: init_L
!!$      integer :: id
!!$      integer :: maxbufsize
!!$      character(len=64) :: vgrid_m_S, vgrid_t_S
!!$   end type INCFG_T

   type :: INCFG_VAR_T
      character(len=INCFG_STRLEN) :: vn_S(2), files_S(NFILE_MAX),  &
           hint_S, vint_S, tint_S, freq_type_S, lvl_type_S, typvar_S
      integer :: freq_t0, freq_dt, freq_t1, lvl0, lvl1, mandatory, cat0, cat1
      integer, pointer :: ip1s(:)
      real :: vmin, vmax
      logical :: active_L
   end type INCFG_VAR_T

   type :: INCFG_T
      logical :: init_L
      character(len=64) :: vgrid_m_S, vgrid_t_S
      integer, pointer :: ip1s(:)  !list of ip1 for lvl=-1
      integer(INT64) :: jdateo  !Initial date [jdate sec]
      integer :: dt         !time step length [seconds]
      integer :: n          !number of elements
      integer :: hgridid     !dest grid id
      integer :: hgridcoreid !dest core grid id
      integer :: commgid     !rpn_comm_create_2dgrid id
      type(incfg_var_T), pointer :: v(:)
   end type INCFG_T

   logical, save :: m_is_init_L = .false.
   type(incfg_var_T), save :: INCFG_VAR_DEFAULT

   interface incfg_add
      module procedure incfg_add_string
      module procedure incfg_add_kv
   end interface incfg_add

   interface incfg_time
      module procedure incfg_time_4
      module procedure incfg_time_8
   end interface incfg_time


contains

   !/@*
   function incfg_new(F_incfgobj, F_jdateo, F_dt, F_filename_S, F_ip1list, &
        F_hgridid, F_hgridcoreid, F_commgid, F_vgrid_m_S, F_vgrid_t_S) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(out) :: F_incfgobj
      integer(INT64), intent(in) :: F_jdateo
      integer, intent(in) :: F_dt
      character(len=*), intent(in), optional :: F_filename_S !full/rel path of config file
      integer, intent(in), optional :: F_ip1list(:)
      integer, intent(in), optional :: F_hgridid
      integer, intent(in), optional :: F_hgridcoreid
      integer, intent(in), optional :: F_commgid
!!$      integer, intent(in), optional :: F_maxbufsize
      character(len=*), intent(in), optional :: F_vgrid_m_S
      character(len=*), intent(in), optional :: F_vgrid_t_S
      !@return
      integer :: F_istat
      !*@/
      integer :: istat, nvar, hgridid, hgridcoreid, commgid
      character(len=1024) :: string_S, filename_S, vgrid_m_S, vgrid_t_S
      !----------------------------------------------------------------------
      filename_S = ''
      hgridid = -1
      hgridcoreid = -1
      commgid = -1
      vgrid_m_S = ''
      vgrid_t_S = ''
      if (present(F_filename_S)) filename_S = F_filename_S
      if (present(F_hgridid)) hgridid = F_hgridid
      if (present(F_hgridcoreid)) hgridcoreid = F_hgridcoreid
      if (present(F_commgid)) commgid = F_commgid
      if (present(F_vgrid_m_S)) vgrid_m_S = F_vgrid_m_S
      if (present(F_vgrid_t_S)) vgrid_t_S = F_vgrid_t_S
!!$      if (present(F_maxbufsize)) F_incfgobj%maxbufsize = F_maxbufsize
      write(string_S, '(a,1x,i0,1x,i0,1x,a)') '(incfg) new [BEGIN]', F_jdateo, F_dt, trim(filename_S)
      call msg(MSG_DEBUG, string_S)

      F_istat = RMN_ERR

      if (.not.m_is_init_L) then
         INCFG_VAR_DEFAULT%vn_S(1:2)  = (/'         ', '         '/)
         INCFG_VAR_DEFAULT%files_S(1:4) = (/'                ', '                ', &
              &                           '                ', '                '/)
         INCFG_VAR_DEFAULT%hint_S = INCFG_KNOWN_H_INT(1)
         INCFG_VAR_DEFAULT%vint_S = INCFG_KNOWN_V_INT(1)
         INCFG_VAR_DEFAULT%tint_S = INCFG_KNOWN_T_INT(1)
         INCFG_VAR_DEFAULT%freq_type_S = 'step'
         INCFG_VAR_DEFAULT%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_ARBITRARY)
         INCFG_VAR_DEFAULT%freq_t0 = 0
         INCFG_VAR_DEFAULT%freq_dt = 0
         INCFG_VAR_DEFAULT%freq_t1 = 0
         INCFG_VAR_DEFAULT%lvl0 = 0
         INCFG_VAR_DEFAULT%lvl1 = 0
         INCFG_VAR_DEFAULT%cat0 = -1
         INCFG_VAR_DEFAULT%cat1 = -1
         INCFG_VAR_DEFAULT%typvar_S = ' '
         INCFG_VAR_DEFAULT%mandatory = -1
         INCFG_VAR_DEFAULT%ip1s => null()
         INCFG_VAR_DEFAULT%vmin = -1*huge(INCFG_VAR_DEFAULT%vmin)
         INCFG_VAR_DEFAULT%vmax = huge(INCFG_VAR_DEFAULT%vmax)
         INCFG_VAR_DEFAULT%active_L = .true.
         m_is_init_L = .true.
      endif

      F_incfgobj%jdateo = F_jdateo
      F_incfgobj%dt = F_dt
      F_incfgobj%n = 0
      F_incfgobj%hgridid = -1
      F_incfgobj%hgridcoreid = -1
      F_incfgobj%commgid = commgid
      nullify(F_incfgobj%ip1s)
      F_incfgobj%vgrid_m_S = vgrid_m_S
      F_incfgobj%vgrid_t_S = vgrid_t_S

      nullify(F_incfgobj%v)
      allocate(F_incfgobj%v(NVAR_MAX), stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(incfg) Problem allocating memory.')
         return
      endif

      F_istat = RMN_OK
      F_incfgobj%init_L = .true.
      if (hgridid >= 0) then
         if (hgridcoreid >= 0) then
            F_istat = min(F_istat, &
                 incfg_sethgridid(F_incfgobj, hgridid, hgridcoreid))
         else
            F_istat = min(F_istat, &
                 incfg_sethgridid(F_incfgobj, hgridid))
         endif
      endif

      nvar = 0
      if (present(F_filename_S) .and. filename_S /= '') then
         call msg(MSG_INFO, '(incfg) Parsing file: '//trim(filename_S))
         nvar = priv_parse_cfgfile(F_incfgobj, F_filename_S)
         if (.not.RMN_IS_OK(nvar)) then
            call msg(MSG_ERROR,'(incfg) Problem parsing file: '//trim(filename_S))
            deallocate(F_incfgobj%v, stat=istat)
            F_incfgobj%init_L = .false.
            F_istat = RMN_ERR
         endif
      endif

      if (RMN_IS_OK(F_istat)) then
         F_incfgobj%n = nvar
         if (present(F_ip1list)) istat = priv_set_ip1s(F_incfgobj, F_ip1list)
      else
         F_incfgobj%init_L = .false.
      endif

      write(string_S, '(a,1x,i0,1x,i0,1x,i0)') '(incfg) new [END]', F_istat, nvar, F_incfgobj%n
      call msg(MSG_DEBUG, string_S)
      !----------------------------------------------------------------------
      return
   end function incfg_new


   !TODO-later: incfg_free


   !/@*
   function incfg_add_string(F_incfgobj, F_string_S) result(F_istat)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      character(len=*), intent(in) :: F_string_S
      integer :: F_istat
      !*@/
      integer, external :: str_split2keyval
      integer, parameter :: NMAX = 32
      integer :: nkeys
      character(len=1024) :: string_S, kv_S(2, NMAX)
      !----------------------------------------------------------------------
      write(string_S, '(a)') '(incfg) add_string [BEGIN] '//trim(F_string_S)
      call msg(MSG_DEBUG, string_S)

      string_S = F_string_S
      call str_tab2space(string_S)
      string_S = adjustl(string_S)
      if (string_S(1:1) /= '#' .and. string_S /= ' ') then
         kv_S = ' '
         nkeys = str_split2keyval(kv_S, string_S, NMAX)
         F_istat = incfg_add(F_incfgobj, kv_S)
      endif
      write(string_S, '(a,1x,i0)') '(incfg) add_string [END]', F_istat
      call msg(MSG_DEBUG, string_S)
      !----------------------------------------------------------------------
      return
   end function incfg_add_string


   !/@*
   function incfg_add_kv(F_incfgobj, F_kv_S) result(F_istat)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      character(len=*), intent(in) :: F_kv_S(:,:)
      integer :: F_istat
      !*@/
      integer :: index, istat, k, ii
      character(len=1024) :: key_S, val_S, string_S
      !----------------------------------------------------------------------
      write(string_S, '(5a)') '(incfg) add_kv [BEGIN] ', trim(F_kv_S(1, 1)), &
           ':=:', trim(F_kv_S(2, 1)), '; ...'
      call msg(MSG_DEBUG, string_S)

      F_istat = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(F_istat)) return

      F_istat = RMN_ERR
      index = F_incfgobj%n + 1
      F_incfgobj%v(index) = INCFG_VAR_DEFAULT
      string_S = ''
!!$      print *, 'incfg_add_kv - nkeys=', size(F_kv_S, 2)
      DOKEYS: do k = 1, size(F_kv_S, 2)
!!$         print *, 'incfg_add_kv - k=', k, trim(F_kv_S(1, k))//'='//trim(F_kv_S(2, k))//';'
         if (F_kv_S(1, k) == ' ') exit DOKEYS
         key_S = F_kv_S(1, k)
         val_S = F_kv_S(2, k)
         istat = clib_tolower(key_S)
         istat = clib_tolower(val_S)
         if (key_S(1:4) == 'inte') key_S = 'hinterp' !#Backward compatibility
         if (key_S(1:4) == 'time') key_S = 'tinterp' !#Backward compatibility
         string_S = trim(string_S)//trim(key_S)//'='//trim(val_S)//';'
         call msg(MSG_DEBUG, 'incfg_add_kv key:'//trim(key_S)//'='//&
              trim(val_S) //';')
         select case(key_S(1:4))
         case('in  ')
            F_incfgobj%v(index)%vn_S(1) = val_S
         case('in2 ')
            F_incfgobj%v(index)%vn_S(2) = val_S
         case('sear') !search
            call priv_parse_file(F_incfgobj, index, val_S)
         case('freq') !freq
            call priv_parse_freq(F_incfgobj, index, val_S)
         case('leve') !levels
            call priv_parse_level(F_incfgobj, index, val_S)
         case('cat ')
            call priv_parse_cat(F_incfgobj, index, val_S)
         case('hint') !hinterp
            do ii=1, size(INCFG_KNOWN_H_INT)
               if (val_S(1:4) == INCFG_KNOWN_H_INT(ii)(1:4)) then
                  val_S = INCFG_KNOWN_H_INT(ii)
                  F_incfgobj%v(index)%hint_S = INCFG_KNOWN_H_INT(ii)
               endif
            enddo
            if (F_incfgobj%v(index)%hint_S == ' ') then
               call msg(MSG_WARNING, '(incfg) using default interp - ignoring unknown interp= ' &
                    //trim(val_S))
            endif
         case('vint') !vinterp
            do ii=1, size(INCFG_KNOWN_V_INT)
               if (val_S(1:4) == INCFG_KNOWN_V_INT(ii)(1:4)) then
                  F_incfgobj%v(index)%vint_S = INCFG_KNOWN_V_INT(ii)
                  val_S = INCFG_KNOWN_V_INT(ii)
               endif
            enddo
            if (F_incfgobj%v(index)%vint_S == ' ') then
               call msg(MSG_WARNING, '(incfg) using default vinterp - ignoring unknown vinterp= ' &
                    //trim(val_S))
            endif
         case('tint') !timeint
            F_incfgobj%v(index)%tint_S = ' '
            do ii=1, size(INCFG_KNOWN_T_INT)
               if (val_S(1:4) == INCFG_KNOWN_T_INT(ii)(1:4)) then
                  F_incfgobj%v(index)%tint_S = INCFG_KNOWN_T_INT(ii)
                  val_S = INCFG_KNOWN_T_INT(ii)
               endif
            enddo
            if (F_incfgobj%v(index)%tint_S == ' ') then
               F_incfgobj%v(index)%tint_S = INCFG_KNOWN_T_INT(1)
               call msg(MSG_WARNING, '(incfg) using default interp - ignoring unknown timeint= ' &
                    //trim(val_S))
            endif
         case('typv') !typvar
            F_incfgobj%v(index)%typvar_S = val_S(1:1)
         case('mand') !mandatory
            F_incfgobj%v(index)%mandatory = -1
            if (val_S(1:4) == 'fals') F_incfgobj%v(index)%mandatory = 0
            if (val_S(1:4) == 'true') F_incfgobj%v(index)%mandatory = 1
         case('vmin')
            istat = str_toreal( F_incfgobj%v(index)%vmin, val_S)
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_WARNING, '(incfg) invalid VMIN= '//trim(val_S))
               return
            endif
         case('vmax')
            istat = str_toreal( F_incfgobj%v(index)%vmax, val_S)
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_WARNING, '(incfg) invalid VMAX= '//trim(val_S))
               return
            endif
         case default ! skip
            call msg(MSG_DEBUG, '(incfg) ignoring param: '// &
                 trim(key_S)//'='//trim(val_S))
         end select
      enddo DOKEYS
      if (F_incfgobj%v(index)%vn_S(1) == ' ' .or. &
           F_incfgobj%v(index)%files_S(1) == ' ' .or. &
           F_incfgobj%v(index)%freq_type_S == 'error') then
         call msg(MSG_WARNING, '(incfg) ignoring var: '//trim(string_S))
         return
      else
         call msg(MSG_INFO, '(incfg) Add var: '//trim(string_S))
      endif
      F_incfgobj%n = index
      F_istat = index
      write(string_S, '(a,1x,i0)') '(incfg) add_kv [END]', F_istat
      call msg(MSG_DEBUG, string_S)
      !----------------------------------------------------------------------
      return
   end function incfg_add_kv


   !/@*
   function incfg_nbvar(F_incfgobj) result(F_nbvar)
      implicit none
      type(INCFG_T), intent(in) :: F_incfgobj
      integer :: F_nbvar
      !*@/
      !------------------------------------------------------------------
      F_nbvar = F_incfgobj%n
      !------------------------------------------------------------------
      return
   end function incfg_nbvar


   !/@*
   function incfg_sethgridid(F_incfgobj, F_hgridid, F_hgridcoreid) result(F_istat)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_hgridid
      integer, intent(in), optional :: F_hgridcoreid
      integer :: F_istat
      !*@/
      !------------------------------------------------------------------
      F_istat = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(F_istat)) return
      if (F_hgridid < 0) then
         call msg(MSG_WARNING, '(incfg) invalid hgridid')
         F_istat = RMN_ERR
         return
      endif
      F_incfgobj%hgridid = F_hgridid
      if (present(F_hgridcoreid)) then
         if (F_hgridcoreid < 0) then
            call msg(MSG_WARNING, '(incfg) invalid hgridcoreid')
            F_istat = RMN_ERR
            return
         endif
         F_incfgobj%hgridcoreid = F_hgridcoreid
      else
         if (F_incfgobj%hgridcoreid < 0) &
              F_incfgobj%hgridcoreid = F_hgridid
      endif
      !------------------------------------------------------------------
      return
   end function incfg_sethgridid


   !/@*
   function incfg_gethgridid(F_incfgobj, F_hgridcoreid) result(F_hgridid)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(out), optional :: F_hgridcoreid
      integer :: F_hgridid
      !*@/
      integer :: istat
      !------------------------------------------------------------------
      F_hgridid = -1
      istat = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(istat)) return
      F_hgridid = F_incfgobj%hgridid
      if (present(F_hgridcoreid)) then
         F_hgridcoreid = F_incfgobj%hgridcoreid
      endif
      !------------------------------------------------------------------
      return
   end function incfg_gethgridid


   !/@*
   function incfg_varindex(F_incfgobj, F_varname_S, F_varname2_S) result(F_index)
      implicit none
      type(INCFG_T), intent(in) :: F_incfgobj
      character(len=*), intent(in) :: F_varname_S
      character(len=*), intent(in), optional :: F_varname2_S
      integer :: F_index
      !*@/
      character(len=INCFG_STRLEN) :: varname_S, varname2_S
      integer :: ii
      !------------------------------------------------------------------
      F_index = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(F_index)) return

      F_index = RMN_ERR
      varname_S = F_varname_S
      varname2_S = ' '
      if (present(F_varname2_S)) varname2_S = F_varname2_S
      ii = clib_tolower(varname_S)
      ii = clib_tolower(varname2_S)
      do ii=1, F_incfgobj%n
         if (varname_S == F_incfgobj%v(ii)%vn_S(1) .and. &
              (varname2_S == F_incfgobj%v(ii)%vn_S(2) .or. varname_S == ' ')) then
            F_index = ii
            exit
         endif
      enddo
      !------------------------------------------------------------------
      return
   end function incfg_varindex


   !/@*
   function incfg_isvarstep(F_incfgobj, F_index, F_istep) result(F_istat)
      implicit none
      type(INCFG_T), intent(in) :: F_incfgobj
      integer, intent(in) :: F_index, F_istep
      integer :: F_istat
      !*@/
      integer :: delta, mo0, mo1, yy0, yy1
      integer(INT64) :: jdatev, step_8, dt_8
      !------------------------------------------------------------------
      F_istat = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(F_istat)) return

      F_istat = RMN_ERR
      if (.not.F_incfgobj%v(F_index)%active_L) return

      select case(F_incfgobj%v(F_index)%freq_type_S(1:2))
      case('st') !STEP
         delta = F_istep
      case('mo') !MONTH
         yy0 = jdate_year(F_incfgobj%jdateo)
         mo0 = jdate_month(F_incfgobj%jdateo)
         step_8 = F_istep
         dt_8   = F_incfgobj%dt
         jdatev = F_incfgobj%jdateo + step_8 * dt_8
         yy1 = jdate_year(jdatev)
         mo1 = jdate_month(jdatev)
         delta = (mo1-mo0) + (yy1-yy0)*12
      end select

      if (delta == F_incfgobj%v(F_index)%freq_t0) then
         F_istat = RMN_OK
      else if (F_incfgobj%v(F_index)%freq_dt > 0 .and. &
           delta >= F_incfgobj%v(F_index)%freq_t0 .and. &
           delta <= F_incfgobj%v(F_index)%freq_t1) then
         delta = delta - F_incfgobj%v(F_index)%freq_t0
         if (mod(delta, F_incfgobj%v(F_index)%freq_dt) == 0) F_istat = RMN_OK
      endif
      !------------------------------------------------------------------
      return
   end function incfg_isvarstep


!!$   !/@*
!!$   function incfg_listvarstep(F_incfgobj, F_istep, F_indexlist) result(F_istat)
!!$      implicit none
!!$      type(INCFG_T), intent(in) :: F_incfgobj
!!$      integer,intent(in) :: F_istep
!!$      integer,intent(out) :: F_indexlist(:)
!!$      integer :: F_istat
!!$      !*@/
!!$      integer :: istat, nbvar, ivar
!!$      !------------------------------------------------------------------
!!$      F_indexlist = -1
!!$      F_istat = priv_check_obj(F_incfgobj)
!!$
!!$      F_istat = 0
!!$      nbvar = incfg_nbvar(F_incfgobj)
!!$      do ivar = 1, nbvar
!!$         istat = incfg_isvarstep(F_incfgobj, ivar, F_istep)
!!$         if (RMN_IS_OK(istat)) then
!!$            if (F_istat >= size(F_indexlist)) then
!!$               F_istat = -1 * F_istat
!!$               call msg(MSG_WARNING, '(incfg) listvarstep, output list is too small')
!!$               return
!!$            endif
!!$            F_istat = F_istat + 1
!!$            F_indexlist(F_istat) = F_id
!!$         end if
!!$      enddo
!!$      !------------------------------------------------------------------
!!$      return
!!$   end function incfg_listvarstep


   !/@*
   function incfg_meta(F_incfgobj, F_index, F_varname_S, F_varname2_S, &
        F_files_S, F_h_int_S, F_v_int_S, F_t_int_S, F_lvl_type_S, F_lvl0, F_lvl1, &
        F_needonce_L, F_ip1list, F_nip1, F_typvar_S, F_mandatory, &
        F_vmin, F_vmax, F_active_L, F_cat0, F_cat1) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(in) :: F_incfgobj
      integer, intent(in) :: F_index
      character(len=*), intent(inout) :: F_varname_S
      character(len=*), intent(inout), optional :: F_varname2_S
      character(len=*), intent(inout), optional :: F_files_S(:)
      character(len=*), intent(inout), optional :: F_h_int_S, F_v_int_S, &
           F_t_int_S, F_lvl_type_S, F_typvar_S
      integer, intent(out), optional :: F_lvl0, F_lvl1, F_ip1list(:), &
           F_nip1, F_mandatory, F_cat0, F_cat1
      real, intent(out), optional :: F_vmin, F_vmax
      logical, intent(out), optional :: F_needonce_L, F_active_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer :: nfiles, mylen, ii
      !----------------------------------------------------------------------
      F_varname_S = ' '
      F_istat = priv_check_obj(F_incfgobj, F_index)
      if (.not.RMN_IS_OK(F_istat)) return

      mylen = min(len_trim(F_incfgobj%v(F_index)%vn_S(1)), len(F_varname_S))
      F_varname_S(1:mylen) = trim(F_incfgobj%v(F_index)%vn_S(1))

      if (present(F_varname2_S)) then
         F_varname2_S = ' '
         mylen = min(len_trim(F_incfgobj%v(F_index)%vn_S(2)), len(F_varname2_S))
         F_varname2_S(1:mylen) = trim(F_incfgobj%v(F_index)%vn_S(2))
      endif

      if (present(F_files_S)) then
         F_files_S(:) = ' '
         nfiles = min(NFILE_MAX, size(F_files_S))
         do ii=1, nfiles
            mylen = min(len_trim(F_incfgobj%v(F_index)%files_S(ii)), &
                 len(F_files_S(ii)))
            F_files_S(ii)(1:mylen) = trim(F_incfgobj%v(F_index)%files_S(ii))
         enddo
      endif

      if (present(F_typvar_S)) then
         F_typvar_S = F_incfgobj%v(F_index)%typvar_S
      endif

      if (present(F_h_int_S)) then
         F_h_int_S = ' '
         mylen = min(len_trim(F_incfgobj%v(F_index)%hint_S), len(F_h_int_S))
         F_h_int_S(1:mylen) = trim(F_incfgobj%v(F_index)%hint_S)
      endif

      if (present(F_v_int_S)) then
         F_v_int_S = ' '
         mylen = min(len_trim(F_incfgobj%v(F_index)%vint_S), len(F_v_int_S))
         F_v_int_S(1:mylen) = trim(F_incfgobj%v(F_index)%vint_S)
      endif

      if (present(F_t_int_S)) then
         F_t_int_S = ' '
         mylen = min(len_trim(F_incfgobj%v(F_index)%tint_S), len(F_t_int_S))
         F_t_int_S(1:mylen) = trim(F_incfgobj%v(F_index)%tint_S)
      endif

      if (present(F_lvl_type_S)) then
         F_lvl_type_S = ' '
         mylen = min(len_trim(F_incfgobj%v(F_index)%lvl_type_S), len(F_lvl_type_S))
         F_lvl_type_S(1:mylen) = trim(F_incfgobj%v(F_index)%lvl_type_S)
      endif

      if (present(F_lvl0)) F_lvl0 = F_incfgobj%v(F_index)%lvl0
      if (present(F_lvl1)) F_lvl1 = F_incfgobj%v(F_index)%lvl1

      if (present(F_cat0)) F_cat0 = F_incfgobj%v(F_index)%cat0
      if (present(F_cat1)) F_cat1 = F_incfgobj%v(F_index)%cat1

      if (present(F_nip1) .and. present(F_ip1list)) then
         F_nip1 = priv_get_ip1s(F_incfgobj, F_index, F_ip1list)
      else if (present(F_nip1)) then
         F_nip1 = priv_get_ip1s(F_incfgobj, F_index)
      else if (present(F_ip1list)) then
         F_ip1list = -1
         ii = priv_get_ip1s(F_incfgobj, F_index, F_ip1list)
      endif

      if (present(F_vmin)) F_vmin = F_incfgobj%v(F_index)%vmin
      if (present(F_vmax)) F_vmax = F_incfgobj%v(F_index)%vmax

      if (present(F_needonce_L)) F_needonce_L = &
           (F_incfgobj%v(F_index)%freq_dt == 0 .or. &
           F_incfgobj%v(F_index)%freq_t1 <= F_incfgobj%v(F_index)%freq_t0)

      if (present(F_mandatory)) F_mandatory = F_incfgobj%v(F_index)%mandatory

      if (present(F_active_L)) F_active_L = F_incfgobj%v(F_index)%active_L
      !----------------------------------------------------------------------
      return
   end function incfg_meta


   !/@*
   function incfg_set_meta(F_incfgobj, F_index, F_varname_S, F_varname2_S, &
        F_files_S, F_h_int_S, F_v_int_S, F_t_int_S, F_lvl_type_S, &
        F_lvl0, F_lvl1, F_ip1list, F_nip1, F_typvar_S,  &
        F_mandatory, F_vmin, F_vmax, F_active_L, F_cat0, F_cat1) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_index
      character(len=*), intent(in), optional :: F_varname_S, F_varname2_S
      character(len=*), intent(in), optional :: F_files_S(:)
      character(len=*), intent(in), optional :: F_h_int_S, F_v_int_S, &
           F_t_int_S, F_lvl_type_S, F_typvar_S
      integer, intent(in), optional :: F_lvl0, F_lvl1, F_ip1list(:), &
           F_nip1, F_mandatory, F_cat0, F_cat1
      real, intent(in), optional :: F_vmin, F_vmax
      logical, intent(in), optional :: F_active_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer :: istat, ii, nfiles
      character(len=INCFG_STRLEN) :: val_S
      !----------------------------------------------------------------------
      F_istat = priv_check_obj(F_incfgobj, F_index)
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_varname_S)) then
         F_incfgobj%v(F_index)%vn_S(1) = trim(F_varname_S)
         istat = clib_tolower(F_incfgobj%v(F_index)%vn_S(1))
      endif

      if (present(F_varname2_S)) then
         F_incfgobj%v(F_index)%vn_S(2) = trim(F_varname2_S)
         istat = clib_tolower(F_incfgobj%v(F_index)%vn_S(2))
      endif

      if (present(F_files_S)) then
         F_incfgobj%v(F_index)%files_S(:) = ' '
         nfiles = min(NFILE_MAX, size(F_files_S))
         do ii=1, nfiles
            F_incfgobj%v(F_index)%files_S(ii) = trim(F_files_S(ii))
            istat = clib_tolower(F_incfgobj%v(F_index)%files_S(ii))
         enddo
      endif

      if (present(F_typvar_S)) then
         F_incfgobj%v(F_index)%typvar_S = trim(F_typvar_S)
         istat = clib_tolower(F_incfgobj%v(F_index)%typvar_S)
      endif

      if (present(F_h_int_S)) then
         val_S = trim(F_h_int_S)
         istat = clib_tolower(val_S)
         istat = RMN_ERR
         do ii=1, size(INCFG_KNOWN_H_INT)
            if (val_S(1:4) == INCFG_KNOWN_H_INT(ii)(1:4)) then
               F_incfgobj%v(F_index)%hint_S = INCFG_KNOWN_H_INT(ii)
               istat = RMN_OK
               exit
            endif
         enddo
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING, '(incfg) ignoring unknown hinterp= '//trim(val_S))
         endif
      endif

      if (present(F_v_int_S)) then
         val_S = trim(F_v_int_S)
         istat = clib_tolower(val_S)
         istat = RMN_ERR
         do ii=1, size(INCFG_KNOWN_V_INT)
            if (val_S(1:4) == INCFG_KNOWN_V_INT(ii)(1:4)) then
               F_incfgobj%v(F_index)%vint_S = INCFG_KNOWN_V_INT(ii)
               istat = RMN_OK
               exit
            endif
         enddo
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING, '(incfg) ignoring unknown vinterp= '//trim(val_S))
         endif
      endif

      if (present(F_t_int_S)) then
         val_S = trim(F_t_int_S)
         istat = clib_tolower(val_S)
         istat = RMN_ERR
         do ii=1, size(INCFG_KNOWN_T_INT)
            if (val_S(1:4) == INCFG_KNOWN_T_INT(ii)(1:4)) then
               F_incfgobj%v(F_index)%tint_S = INCFG_KNOWN_T_INT(ii)
               istat = RMN_OK
               exit
            endif
         enddo
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING, '(incfg) ignoring unknown tinterp= '//trim(val_S))
         endif
      endif

      if (present(F_lvl_type_S)) then
         val_S = trim(F_lvl_type_S)
         istat = clib_tolower(val_S)
         istat = RMN_ERR
         do ii=1, size(INCFG_KNOWN_LVL_TYPE)
            if (val_S(1:3) == INCFG_KNOWN_LVL_TYPE(ii)(1:3)) then
               F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(ii)
               istat = RMN_OK
               exit
            endif
         enddo
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING, '(incfg) ignoring unknown lvl_type= '//trim(val_S))
         endif
      endif

      !TODO: check if valid
      if (present(F_lvl0)) F_incfgobj%v(F_index)%lvl0 = F_lvl0

      !TODO: check if valid
      if (present(F_lvl1)) F_incfgobj%v(F_index)%lvl1 = F_lvl1

      !TODO: check if valid
      if (present(F_cat0)) F_incfgobj%v(F_index)%cat0 = F_cat0

      !TODO: check if valid
      if (present(F_cat1)) F_incfgobj%v(F_index)%cat1 = F_cat1

      !TODO: allow update of nip1 and ip1s
!!$      if (present(F_nip1) .and. present(F_ip1list)) then
!!$         F_nip1 = priv_get_ip1s(F_incfgobj, F_index, F_ip1list)
!!$      else if (present(F_nip1)) then
!!$         F_nip1 = priv_get_ip1s(F_incfgobj, F_index)
!!$      else if (present(F_ip1list)) then
!!$         F_ip1list = -1
!!$         ii = priv_get_ip1s(F_incfgobj, F_index, F_ip1list)
!!$      endif
      if (present(F_nip1) .or. present(F_ip1list)) then
         call msg(MSG_WARNING, '(incfg) Updating nip1, ip1list not yet implemented -- Ignored')
      endif

      if (present(F_vmin)) F_incfgobj%v(F_index)%vmin = F_vmin
      if (present(F_vmax)) F_incfgobj%v(F_index)%vmax = F_vmax

      if (present(F_mandatory)) F_incfgobj%v(F_index)%mandatory = F_mandatory

      if (present(F_active_L)) F_incfgobj%v(F_index)%active_L = F_active_L
      !----------------------------------------------------------------------
      return
   end function incfg_set_meta


   !/@*
   function incfg_time_4(F_incfgobj, F_step, F_datev, F_dateo, F_dt) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_step
      integer, intent(out) :: F_datev, F_dateo, F_dt
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer(INT64) :: jdatev, jdateo
      !----------------------------------------------------------------------
      F_datev = RMN_ERR
      F_dateo = RMN_ERR
      F_istat =  incfg_time_8(F_incfgobj, F_step, jdatev, jdateo, F_dt)
      if (RMN_IS_OK(F_istat)) then
         F_datev = jdate_to_cmc(jdatev)
         F_dateo = jdate_to_cmc(jdateo)
      endif
      !----------------------------------------------------------------------
      return
   end function incfg_time_4


   !/@*
   function incfg_time_8(F_incfgobj, F_step, F_datev, F_dateo, F_dt) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_step
      integer(INT64), intent(out) :: F_datev
      integer(INT64), intent(out), optional :: F_dateo
      integer, intent(out), optional :: F_dt
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer(INT64) :: step_8, dt_8
      !----------------------------------------------------------------------
      F_istat = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(F_istat)) return

      step_8 = F_step
      dt_8   = F_incfgobj%dt
      F_datev = F_incfgobj%jdateo + step_8 * dt_8
      if (present(F_dateo)) F_dateo = F_incfgobj%jdateo
      if (present(F_dt)) F_dt = F_incfgobj%dt
      !----------------------------------------------------------------------
      return
   end function incfg_time_8

   !==== Private Functions =================================================


   !/@*
   function priv_check_obj(F_incfgobj, F_index)  result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(in) :: F_incfgobj
      integer, intent(in), optional :: F_index
      !@return
      integer :: F_istat
      !*@/
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (.not.(m_is_init_L .and. F_incfgobj%init_L .and. &
           associated(F_incfgobj%v))) then
         call msg(MSG_WARNING, '(incfg) incfg obj not properly init')
         return
      endif
       if (present(F_index)) then
         if (F_index < 1 .or. F_index > F_incfgobj%n .or. &
              F_index > size(F_incfgobj%v)) then
            call msg(MSG_WARNING, '(incfg) index out of range')
            return
         endif
      endif
      F_istat = RMN_OK
      !----------------------------------------------------------------------
      return
   end function priv_check_obj


   !/@*
   function priv_parse_cfgfile(F_incfgobj, F_filename_S) result(F_nvar)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(inout) :: F_incfgobj
      character(len=*), intent(in) :: F_filename_S !full/rel path of config file
      !@return
      integer :: F_nvar
      !@author Stephane Chamberland, 2011-07
      !*@/
      integer :: istat, fileid
      character(len=1024) :: string_S
      !---------------------------------------------------------------------
      write(string_S, '(a)') 'priv_parse_cfgfile: '//trim(F_filename_S)
      call msg(MSG_DEBUG, string_S)

      F_nvar = RMN_ERR

      istat = clib_isfile(trim(F_filename_S))
      istat = min(clib_isreadok(trim(F_filename_S)), istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, &
              '(incfg) File not found or not readable: '//trim(F_filename_S))
         return
      endif

      fileid = 0
      istat = fnom(fileid, F_filename_S, 'SEQ/FMT+R/O+OLD', 0)
      if (.not.RMN_IS_OK(istat) .or. fileid <= 0) then
         call msg(MSG_WARNING, '(incfg) Problem opening file: '//trim(F_filename_S))
         return
      endif

      F_nvar = 0
      DOFILE: do
         read(fileid, '(a)', iostat=istat) string_S
         if (istat /=0 .or. F_nvar+1 > NVAR_MAX) exit DOFILE
         call str_tab2space(string_S)
         string_S = adjustl(string_S)
         if (string_S(1:1)/='#' .and. string_S/=' ') then
            istat = incfg_add_string(F_incfgobj, string_S)
            if (.not.RMN_IS_OK(istat)) then
!!$               F_nvar = RMN_ERR
!!$               exit DOFILE
               cycle DOFILE
            endif
            F_nvar = F_nvar + 1
         endif
      enddo DOFILE
      istat = fclos(fileid)
      !----------------------------------------------------------------------
      return
   end function priv_parse_cfgfile


   !/@*
   subroutine priv_parse_file(F_incfgobj, F_index, F_string_S)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_index
      character(len=*), intent(in) :: F_string_S
      !*@/
      integer :: ii
      character(len=1024) :: parts_S(NFILE_MAX+1)
      !------------------------------------------------------------------
      call str_split2list(parts_S, F_string_S, ',', NFILE_MAX+1)
      ii = 1
      ITEMLOOP: do
         if (ii == NFILE_MAX+1) exit ITEMLOOP
         if (parts_S(ii) == ' ') exit ITEMLOOP
         F_incfgobj%v(F_index)%files_S(ii) = parts_S(ii)
         ii = ii + 1
      enddo ITEMLOOP
      if (parts_S(NFILE_MAX+1) /= ' ') &
           call msg(MSG_WARNING, '(incfg) ignoring extra files (max of 4): '// &
           &                     trim(parts_S(NFILE_MAX+1)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_file


   !/@*
   subroutine priv_parse_freq(F_incfgobj, F_index, F_string_S)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_index
      character(len=*), intent(in) :: F_string_S
      !*@/
      integer, parameter :: NMAX = 5
      character(len=1024) :: parts_S(NMAX)
      integer :: ii, val, istat
      real :: rfact, rt0, rdt, rt1
      !------------------------------------------------------------------
      call str_split2list(parts_S, F_string_S, ',', NMAX)
      rfact = 1.
      F_incfgobj%v(F_index)%freq_type_S = 'error'
      ii = 2
      select case(parts_S(1)(1:2))
      case('st') !STEP
         F_incfgobj%v(F_index)%freq_type_S = 'step'//'       '
      case('mo') !MONTH
         F_incfgobj%v(F_index)%freq_type_S = 'month'//'       '
      case('we') !WEEK
         F_incfgobj%v(F_index)%freq_type_S = 'step'//'       '
         rfact = real(60*60*24*7)/real(F_incfgobj%dt)
      case('da') !DAY
         F_incfgobj%v(F_index)%freq_type_S = 'step'//'       '
         rfact = real(60*60*24)/real(F_incfgobj%dt)
      case('ho') !HOURS
         F_incfgobj%v(F_index)%freq_type_S = 'step'//'       '
         rfact = real(60*60)/real(F_incfgobj%dt)
      case('mi') !MINUTES
         F_incfgobj%v(F_index)%freq_type_S = 'step'//'       '
         rfact = real(60)/real(F_incfgobj%dt)
      case('  ')
         return
      case default
         F_incfgobj%v(F_index)%freq_type_S = 'step'//'       '
         ii = 1
      end select

      read(parts_S(ii), fmt=*, iostat=istat) val
      if (istat /= 0 .or. parts_S(ii) == ' ') then
         F_incfgobj%v(F_index)%freq_type_S = 'error'//'       '
         call msg(MSG_WARNING, '(incfg) ignoring invalid freq: '//trim(F_string_S))
         return
      endif
      F_incfgobj%v(F_index)%freq_t0 = max(0, val)

      ii = ii + 1
      F_incfgobj%v(F_index)%freq_dt = 0
      if (parts_S(ii) /= ' ') then
         read(parts_S(ii), fmt=*, iostat=istat) val
         if (istat /= 0) then
            F_incfgobj%v(F_index)%freq_type_S = 'error'//'       '
            call msg(MSG_WARNING, '(incfg) ignoring invalid freq: '// &
                 trim(F_string_S))
            return
         endif
         F_incfgobj%v(F_index)%freq_dt = max(0, val)
      endif

      ii = ii + 1
      F_incfgobj%v(F_index)%freq_t1 = F_incfgobj%v(F_index)%freq_t0
      if (F_incfgobj%v(F_index)%freq_dt > 0) then
         F_incfgobj%v(F_index)%freq_t1 = -1
         if (parts_S(ii) /= ' ') then
            read(parts_S(ii), fmt=*, iostat=istat) val
            if (istat /= 0) then
               F_incfgobj%v(F_index)%freq_type_S = 'error'//'       '
               call msg(MSG_WARNING, '(incfg) ignoring invalid freq: '// &
                    trim(F_string_S))
               return
            endif
            F_incfgobj%v(F_index)%freq_t1 = max(0, val)
         endif
      endif

      if (abs(1. - rfact) > EPSILON) then
         rdt = real(F_incfgobj%v(F_index)%freq_dt) * rfact
         rt0 = real(F_incfgobj%v(F_index)%freq_t0) * rfact
         if (abs(rt0 - real(nint(rt0))) > EPSILON .or. &
              abs(rdt - real(nint(rdt))) > EPSILON) then
            F_incfgobj%v(F_index)%freq_type_S = 'error'//'       '
            call msg(MSG_WARNING, '(incfg) cannot satisfy frequency: '// &
                 trim(F_string_S))
            return
         endif
         F_incfgobj%v(F_index)%freq_t0 = nint(rt0)
         F_incfgobj%v(F_index)%freq_dt = nint(rdt)
         if (F_incfgobj%v(F_index)%freq_t1 > 0) then
            rt1 = floor(real(F_incfgobj%v(F_index)%freq_t1) * rfact)
            F_incfgobj%v(F_index)%freq_t1 = nint(max(rt0, rt1))
         endif
      endif
      if (F_incfgobj%v(F_index)%freq_t1 < 0) then
         F_incfgobj%v(F_index)%freq_t1 = huge(val)
      end if

      if (parts_S(NMAX) /= ' ') &
           call msg(MSG_WARNING, '(incfg) ignoring extra freq params: '// &
           &                      trim(parts_S(NMAX)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_freq


   !/@*
   subroutine priv_parse_level(F_incfgobj, F_index, F_string_S)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_index
      character(len=*), intent(in) :: F_string_S
      !*@/
      integer, parameter :: NMAX = 3
      character(len=1024) :: parts_S(NMAX), tmp_S
      integer :: val, istat, n
      !------------------------------------------------------------------
      call str_split2list(parts_S, F_string_S, ',', NMAX)

      F_incfgobj%v(F_index)%lvl0 = -1
      F_incfgobj%v(F_index)%lvl1 = -1

      n = 1
      tmp_S = parts_S(n)
      select case(tmp_S(1:3))
      case('arb') !arbitrary
         F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_ARBITRARY)
         return
      case('tdi') !diagnostic themo
         F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_TDIAG)
         F_incfgobj%v(F_index)%lvl0 = 0
         F_incfgobj%v(F_index)%lvl1 = 0
         return
      case('mdi') !diagnostic momentum
         F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_MDIAG)
         F_incfgobj%v(F_index)%lvl0 = 0
         F_incfgobj%v(F_index)%lvl1 = 0
         return
      case('mom') !momentum
         F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_MOM)
         return
      case('the') !thermodynamic
         F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_THERMO)
         return
      case default ! skip
!!$         call msg(MSG_WARNING, '(incfg) ignoring invalid level: '//trim(F_string_S))
         F_incfgobj%v(F_index)%lvl_type_S = INCFG_KNOWN_LVL_TYPE(INCFG_LVL_ARBITRARY)
      end select

      read(parts_S(n), fmt=*, iostat=istat) val
      if (istat /= 0 .or. parts_S(n) == ' ') then
         call msg(MSG_WARNING, '(incfg) ignoring invalid level: '//trim(F_string_S))
         return
      endif
      F_incfgobj%v(F_index)%lvl0 = max(-1, val) !#max(0, val)
      n = n+1

      F_incfgobj%v(F_index)%lvl1 = F_incfgobj%v(F_index)%lvl0
      if (val >= 0) then
         read(parts_S(n), fmt=*, iostat=istat) val
         if (istat == 0 .and. parts_S(n) /= ' ') then
            F_incfgobj%v(F_index)%lvl1 = max(0, val)
         endif
         n = n+1
      endif

      if (parts_S(n) /= ' ') call msg(MSG_WARNING, &
           '(incfg) ignoring extra level params: '//trim(parts_S(n)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_level


   !/@*
   subroutine priv_parse_cat(F_incfgobj, F_index, F_string_S)
      implicit none
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_index
      character(len=*), intent(in) :: F_string_S
      !*@/
      integer, parameter :: NMAX = 3
      character(len=1024) :: parts_S(NMAX)
      integer :: val, istat, n
      !------------------------------------------------------------------
      call str_split2list(parts_S, F_string_S, ',', NMAX)

      F_incfgobj%v(F_index)%cat0 = -1
      F_incfgobj%v(F_index)%cat1 = -1

      n = 1
      read(parts_S(n), fmt=*, iostat=istat) val
      if (istat /= 0 .or. parts_S(n) == ' ') then
         call msg(MSG_WARNING, '(incfg) ignoring invalid cat: '//trim(F_string_S))
         return
      endif
      F_incfgobj%v(F_index)%cat0 = max(-1, val) !#max(0, val)
      n = n+1
      F_incfgobj%v(F_index)%cat1 = F_incfgobj%v(F_index)%cat0
      if (val >= 0) then
         read(parts_S(n), fmt=*, iostat=istat) val
         if (istat == 0 .and. parts_S(n) /= ' ') then
            F_incfgobj%v(F_index)%cat1 = max(0, val)
         endif
         n = n+1
      endif

      if (parts_S(n) /= ' ') call msg(MSG_WARNING, &
           '(incfg) ignoring extra cat params: '//trim(parts_S(n)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_cat


   !/@*
   function priv_set_ip1s(F_incfgobj, F_ip1list) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(inout) :: F_incfgobj
      integer, intent(in) :: F_ip1list(:)
      !@return
      integer :: F_istat
      !*@/
      integer :: istat
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(incfg) set_ip1 [BEGIN]')
      F_istat = priv_check_obj(F_incfgobj)
      if (.not.RMN_IS_OK(F_istat)) return
      allocate(F_incfgobj%ip1s(size(F_ip1list)), stat=istat)
      F_incfgobj%ip1s(:) = F_ip1list(:)
      F_istat = RMN_OK
      call msg(MSG_DEBUG, '(incfg) set_ip1 [END]')
      !----------------------------------------------------------------------
      return
   end function priv_set_ip1s


   !/@*
   function priv_get_ip1s(F_incfgobj, F_index, F_ip1list) result(F_nip1)
      implicit none
      !@objective
      !@arguments
      type(INCFG_T), intent(in) :: F_incfgobj
      integer, intent(in) :: F_index
      integer, intent(out), optional :: F_ip1list(:)
      !@return
      integer :: F_nip1
      !*@/
      integer :: lvl0, lvl1, k, nmax, ipkind
      real :: zp1
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(incfg) get_ip1 [BEGIN]')
      F_nip1 = RMN_ERR

      lvl0 = -1
      lvl1 = -1
      if (F_index > 0) then
         lvl0 = F_incfgobj%v(F_index)%lvl0
         lvl1 = F_incfgobj%v(F_index)%lvl1
      endif

      nmax = 999999
      F_nip1 = 1
      if (present(F_ip1list)) then
         nmax = size(F_ip1list)
         F_ip1list(1) = 0
      endif
      if (lvl0 < 0) then
         if (associated(F_incfgobj%ip1s)) then
            F_nip1 = min(nmax, size(F_incfgobj%ip1s))
            if (present(F_ip1list)) &
                 F_ip1list(1:F_nip1) = F_incfgobj%ip1s(1:F_nip1)
         else
            F_nip1 = 1
            F_ip1list(1) = -1
         endif
      else if (lvl0 > 0) then
         ipkind = RMN_CONV_ARBITRARY
         F_nip1 = 0
         do k=lvl0, lvl1
            F_nip1 = F_nip1 + 1
            if (F_nip1 > nmax) cycle
            if (present(F_ip1list)) then
               zp1 = real(k)
               call convip_plus(F_ip1list(F_nip1), zp1, ipkind, RMN_CONV_P2IPNEW, &
                    ' ', .not.RMN_CONV_USEFORMAT_L)
            endif
         enddo
      endif

      !# Backward compatibility when levels is not specified in the input table
      if (present(F_ip1list) .and. F_nip1 == 1 .and. &
           F_incfgobj%v(F_index)%vint_S /= INCFG_KNOWN_V_INT(1)) then
         if (F_ip1list(1) == 0) F_ip1list(1) = -1
      endif

      call msg(MSG_DEBUG, '(incfg) get_ip1 [END]')
      !----------------------------------------------------------------------
      return
   end function priv_get_ip1s

end module incfg2_mod
