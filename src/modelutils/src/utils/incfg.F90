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
module incfg_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower, clib_isfile, clib_isreadok
   use str_mod, only: str_tab2space, str_toreal
   use mu_jdate_mod
   implicit none
   private
   !@objective 
   !@author Michel Desgagne, 2011-04
   !@revision
   !  2011-04 Stephane Chamberland - module, options and automation
   !  2013-11 Stephane Chamberland - typvar
   !@description
   ! Public functions
   public :: incfg_new, incfg_new_4, incfg_new_8, incfg_add, incfg_nbvar, &
        incfg_meta, incfg_isvarstep, incfg_varindex,incfg_time,incfg_add_string, &
        incfg_add_kv,incfg_setgridid,incfg_getgridid
   ! Public constants
   integer,parameter,public :: INCFG_STRLEN = 32
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

   real,parameter :: EPSILON = 1.e-5
   integer,parameter :: NLIST_MAX = 16
   integer,parameter :: NFILE_MAX = 4
   integer,parameter :: NVAR_MAX = 512

   type :: incfg_var_T
      character(len=INCFG_STRLEN) :: vn(2),file(NFILE_MAX),h_int,v_int, &
           t_int,freq_type,lvl_type,typvar
      integer :: freq_t0,freq_dt,freq_t1,lvl0,lvl1,mandatory,cat0,cat1
      real :: vmin, vmax
   end type incfg_var_T

   type :: incfg_T
      integer,pointer :: ip1list(:)  !list of ip1 for lvl=-1
      integer(INT64) :: dateo  !Initial date [jdate sec]
      integer :: dt         !time step length [seconds]
      integer :: n          !number of elements
      integer :: gridid     !dest grid id
      type(incfg_var_T),pointer :: v(:)
   end type incfg_T

   character(len=INCFG_STRLEN),parameter :: KNOWN_H_INT(3) = (/ &
        'cubic           ','linear          ','nearest         ' &
        /) !TODO-later: KNOWN_H_INT = none or crop
   character(len=INCFG_STRLEN),parameter :: KNOWN_V_INT(5) = (/ &
        'none            ','linear          ','cubic           ', &
        'l-cond          ','c-cond          ' &
        /) !TODO-later: KNOWN_V_INT = nearest
   character(len=INCFG_STRLEN),parameter :: KNOWN_LEVEL_TYPE(5) = (/ &
        'arbitrary       ','momentum        ','thermodynamic   ', &
        'tdiagnostic     ','mdiagnostic     ' &
        /)
   integer,parameter :: LEVEL_ARBITRARY = 1
   integer,parameter :: LEVEL_TDIAG = 4
   integer,parameter :: LEVEL_MDIAG = 5

   character(len=INCFG_STRLEN),parameter :: KNOWN_T_INT(7) = (/ &
        'none            ','linear          ','nearest         ', &
        'any             ','step            ','next            ', &
        'increment       '/)

   type(incfg_var_T),save :: INCFG_VAR_DEFAULT

   type(incfg_T),pointer,save :: m_list(:) => null()
   integer,save :: m_nlist = 0

   interface incfg_new
      module procedure incfg_new_4
      module procedure incfg_new_8
   end interface incfg_new

   interface incfg_add
      module procedure incfg_add_string
      module procedure incfg_add_kv
   end interface

   interface incfg_time
      module procedure incfg_time_4
      module procedure incfg_time_8
   end interface incfg_time

contains

   !/@*
   function incfg_new_4(F_dateo,F_dt,F_filename_S,F_ip1list) result(F_id)
      implicit none
      !@arguments
      integer,intent(in) :: F_dateo,F_dt
      character(len=*),intent(in),optional :: F_filename_S !full/rel path of config file
      integer,intent(in),optional :: F_ip1list(:)
      !@return
      integer :: F_id !config id
      !*@/
      integer(INT64) :: jdateo
      character(len=1024) :: filename_S
      !----------------------------------------------------------------------
      filename_S = ''
      if (present(F_filename_S)) filename_S = F_filename_S
      jdateo = jdate_from_cmc(F_dateo)
      if (present(F_ip1list)) then
         F_id = incfg_new_8(jdateo,F_dt,filename_S,F_ip1list)
      else
         F_id = incfg_new_8(jdateo,F_dt,filename_S)
      endif
      !----------------------------------------------------------------------
      return
   end function incfg_new_4

   !/@*
   function incfg_new_8(F_dateo,F_dt,F_filename_S,F_ip1list) result(F_id)
      implicit none
      !@objective
      !@arguments
      integer(INT64),intent(in) :: F_dateo
      integer,intent(in) :: F_dt
      character(len=*),intent(in),optional :: F_filename_S !full/rel path of config file
      integer,intent(in),optional :: F_ip1list(:)
      !@return
      integer :: F_id !config id
      !@author Stephane Chamberland, 2011-04
      !@revision
      !  2012-04, S.Chamberland: add ip1list
      !*@/
      logical,save :: is_init_L = .false.
      integer :: istat,nvar
      character(len=1024) :: string_S,filename_S
      !----------------------------------------------------------------------
      filename_S = ''
      if (present(F_filename_S)) filename_S = F_filename_S
      write(string_S, '(a,1x,i0,1x,i0,1x,a)') '(incfg) new [BEGIN]',F_dateo,F_dt,trim(filename_S)
      call msg(MSG_DEBUG,string_S)

      F_id = RMN_ERR

      if (.not.is_init_L) then
         allocate(m_list(NLIST_MAX),stat=istat)
         if (istat /= 0) return
         INCFG_VAR_DEFAULT%vn(1:2) = (/'         ','         '/)
         INCFG_VAR_DEFAULT%file(1:4) = (/'                ','                ', &
              '                ','                '/)
         INCFG_VAR_DEFAULT%h_int = KNOWN_H_INT(1)
         INCFG_VAR_DEFAULT%v_int = KNOWN_V_INT(1)
         INCFG_VAR_DEFAULT%t_int = KNOWN_T_INT(1)
         INCFG_VAR_DEFAULT%freq_type = 'step'
         INCFG_VAR_DEFAULT%lvl_type = KNOWN_LEVEL_TYPE(LEVEL_ARBITRARY)
         INCFG_VAR_DEFAULT%freq_t0 = 0
         INCFG_VAR_DEFAULT%freq_dt = 0
         INCFG_VAR_DEFAULT%freq_t1 = 0
         INCFG_VAR_DEFAULT%lvl0 = 0
         INCFG_VAR_DEFAULT%lvl1 = 0
         INCFG_VAR_DEFAULT%cat0 = -1
         INCFG_VAR_DEFAULT%cat1 = -1
         INCFG_VAR_DEFAULT%typvar = ' '
         INCFG_VAR_DEFAULT%mandatory = -1
         INCFG_VAR_DEFAULT%vmin = -1*huge(INCFG_VAR_DEFAULT%vmin)
         INCFG_VAR_DEFAULT%vmax = huge(INCFG_VAR_DEFAULT%vmax)
         is_init_L = .true.
      endif

      if (m_nlist == NLIST_MAX) then
         call msg(MSG_WARNING,'(incfg) Too many opened incfg.')
         return
      endif

      m_nlist = m_nlist + 1

      nullify(m_list(m_nlist)%v)
      allocate(m_list(m_nlist)%v(NVAR_MAX),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(incfg) Problem allocating memory.')
         return
      endif

      m_list(m_nlist)%dateo = F_dateo
      m_list(m_nlist)%dt = F_dt
      m_list(m_nlist)%n = 0
      m_list(m_nlist)%gridid = -1
      nullify(m_list(m_nlist)%ip1list)

      nvar = 0
      if (present(F_filename_S) .and. filename_S /= '') then
         nvar = priv_parse_cfgfile(m_nlist,F_filename_S)
         if (.not.RMN_IS_OK(nvar)) then
            deallocate(m_list(m_nlist)%v,stat=istat)
            m_nlist = m_nlist - 1
            return
         endif
      endif

      if (present(F_ip1list)) then
         istat = priv_set_ip1list(m_nlist,F_ip1list)
      endif

      F_id = m_nlist
      m_list(m_nlist)%n = nvar

      write(string_S, '(a,1x,i0,1x,i0,1x,i0)') '(incfg) new [END]',m_nlist,nvar,m_list(m_nlist)%n
      call msg(MSG_DEBUG,string_S)
      !----------------------------------------------------------------------
      return
   end function incfg_new_8


   !TODO-later: incfg_free


   !/@*
   function incfg_add_string(F_id,F_string_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_string_S
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-07
      !*@/
      integer,external :: str_split2keyval
      integer,parameter :: NMAX = 32
      integer :: nkeys
      character(len=1024) :: string_S,kv_S(2,NMAX)
      !----------------------------------------------------------------------
      write(string_S, '(a,1x,i0,1x,a)') '(incfg) add_string [BEGIN]',F_id,trim(F_string_S)
      call msg(MSG_DEBUG,string_S)

      F_istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(F_istat)) return

      string_S = F_string_S
      call str_tab2space(string_S)
      string_S = adjustl(string_S)
      if (string_S(1:1)/='#' .and. string_S/=' ') then
         kv_S = ' '
         nkeys = str_split2keyval(kv_S,string_S,NMAX)
         F_istat = incfg_add_kv(F_id,kv_S)
      endif
      write(string_S, '(a,1x,i0)') '(incfg) add_string [END]',F_istat
      call msg(MSG_DEBUG,string_S)
      !----------------------------------------------------------------------
      return
   end function incfg_add_string


   !/@*
   function incfg_add_kv(F_id,F_kv_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_kv_S(:,:)
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-07
      !*@/
      integer :: index, istat,k,ii
      character(len=1024) :: key_S,val_S,string_S
      !----------------------------------------------------------------------
      write(string_S, '(a,1x,i0,1x,a)') '(incfg) add_kv [BEGIN]',F_id,trim(F_kv_S(1,1))//':=:'// &
           trim(F_kv_S(2,1))//'; ...'
      call msg(MSG_DEBUG,string_S)

      F_istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(F_istat)) return

      F_istat = RMN_ERR
      if (.not.associated(m_list(F_id)%v)) return

      index = m_list(F_id)%n + 1
      m_list(F_id)%v(index) = INCFG_VAR_DEFAULT
      string_S = ''
!!$      print *,'incfg_add_kv - nkeys=',size(F_kv_S,2)
      DOKEYS: do k = 1, size(F_kv_S,2)
!!$         print *,'incfg_add_kv - k=',k,trim(F_kv_S(1,k))//'='//trim(F_kv_S(2,k))//';'
         if (F_kv_S(1,k) == ' ') exit DOKEYS
         key_S = F_kv_S(1,k)
         val_S = F_kv_S(2,k)
         istat = clib_tolower(key_S)
         istat = clib_tolower(val_S)
         string_S = trim(string_S)//trim(key_S)//'='//trim(val_S)//';'
         call msg(MSG_DEBUG,'incfg_add_kv key:'//trim(key_S)//'='//trim(val_S)//';')
         if (key_S(1:4) == 'inte') key_S = 'hinterp' !#Backward compatibility
         if (key_S(1:4) == 'time') key_S = 'tinterp' !#Backward compatibility
         select case(key_S(1:4))
         case('in  ')
            m_list(F_id)%v(index)%vn(1) = val_S
         case('in2 ')
            m_list(F_id)%v(index)%vn(2) = val_S
         case('sear') !search
            call priv_parse_file(F_id,index,val_S)
         case('freq') !freq
            call priv_parse_freq(F_id,index,val_S)
         case('leve') !levels
            call priv_parse_level(F_id,index,val_S)
         case('cat ')
            call priv_parse_cat(F_id,index,val_S)
         case('hint') !hinterp
            do ii=1,size(KNOWN_H_INT)
               if (val_S(1:4) == KNOWN_H_INT(ii)(1:4)) &
                    m_list(F_id)%v(index)%h_int = KNOWN_H_INT(ii)
            enddo
            if (m_list(F_id)%v(index)%h_int == ' ') then
               call msg(MSG_WARNING,'(incfg) using default interp - ignoring unknown interp= ' &
                    //trim(val_S))
            endif
         case('vint') !vinterp
            do ii=1,size(KNOWN_V_INT)
               if (val_S(1:4) == KNOWN_V_INT(ii)(1:4)) &
                    m_list(F_id)%v(index)%v_int = KNOWN_V_INT(ii)
            enddo
            if (m_list(F_id)%v(index)%v_int == ' ') then
               call msg(MSG_WARNING,'(incfg) using default vinterp - ignoring unknown vinterp= ' &
                    //trim(val_S))
            endif
         case('tint') !timeint
            m_list(F_id)%v(index)%t_int = ' '
            do ii=1,size(KNOWN_T_INT)
               if (val_S(1:4) == KNOWN_T_INT(ii)(1:4)) &
                    m_list(F_id)%v(index)%t_int = KNOWN_T_INT(ii)
            enddo
            if (m_list(F_id)%v(index)%t_int == ' ') then
               m_list(F_id)%v(index)%t_int = KNOWN_T_INT(1)
               call msg(MSG_WARNING,'(incfg) using default interp - ignoring unknown timeint= ' &
                    //trim(val_S))
            endif
         case('typv') !typvar
            m_list(F_id)%v(index)%typvar = val_S(1:1)
         case('mand') !mandatory
            m_list(F_id)%v(index)%mandatory = -1
            if (val_S(1:4) == 'fals') m_list(F_id)%v(index)%mandatory = 0
            if (val_S(1:4) == 'true') m_list(F_id)%v(index)%mandatory = 1
         case('vmin')
            istat = str_toreal( m_list(F_id)%v(index)%vmin, val_S)
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_WARNING,'(incfg) invalid VMIN= '//trim(val_S))
               return
            endif
         case('vmax')
            istat = str_toreal( m_list(F_id)%v(index)%vmax, val_S)
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_WARNING,'(incfg) invalid VMAX= '//trim(val_S))
               return
            endif
         case default ! skip
            call msg(MSG_DEBUG,'(incfg) ignoring param: '//trim(key_S)//'='//trim(val_S))
         end select
      enddo DOKEYS
      if (m_list(F_id)%v(index)%vn(1) == ' ' .or. &
           m_list(F_id)%v(index)%file(1) == ' ' .or. &
           m_list(F_id)%v(index)%freq_type == 'error') then
         call msg(MSG_WARNING,'(incfg) ignoring var: '//trim(string_S))
         return
      endif
      m_list(F_id)%n = index
      F_istat = index
      write(string_S, '(a,1x,i0)') '(incfg) add_kv [END]',F_istat
      call msg(MSG_DEBUG,string_S)
      !----------------------------------------------------------------------
      return
   end function incfg_add_kv


   !/@*
   function incfg_nbvar(F_id) result(F_nbvar)
      implicit none
      integer,intent(in) :: F_id
      integer :: F_nbvar
      !*@/
      integer :: istat
      !------------------------------------------------------------------
      F_nbvar = 0
      istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(istat)) return
      F_nbvar = m_list(F_id)%n
      !------------------------------------------------------------------
      return
   end function incfg_nbvar


   !/@*
   function incfg_setgridid(F_id,F_gridid) result(F_istat)
      implicit none
      integer,intent(in) :: F_id,F_gridid
      integer :: F_istat
      !*@/
      !------------------------------------------------------------------
      F_istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(F_istat)) return
      m_list(F_id)%gridid = F_gridid
      !------------------------------------------------------------------
      return
   end function incfg_setgridid


   !/@*
   function incfg_getgridid(F_id) result(F_gridid)
      implicit none
      integer,intent(in) :: F_id
      integer :: F_gridid
      !*@/
      integer :: istat
      !------------------------------------------------------------------
      F_gridid = -1
      istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(istat)) return
      F_gridid = m_list(F_id)%gridid
      !------------------------------------------------------------------
      return
   end function incfg_getgridid


   !/@*
   function incfg_varindex(F_id,F_varname_S,F_varname2_S) result(F_index)
      implicit none
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_varname_S
      character(len=*),intent(in),optional :: F_varname2_S
      integer :: F_index
      !*@/
      character(len=INCFG_STRLEN) :: varname_S,varname2_S
      integer :: ii
      !------------------------------------------------------------------
      F_index = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(F_index)) return

      F_index = RMN_ERR
      varname_S = F_varname_S
      varname2_S = ' '
      if (present(F_varname2_S)) varname2_S = F_varname2_S
      ii = clib_tolower(varname_S)
      ii = clib_tolower(varname2_S)
      do ii=1,m_list(F_id)%n
         if (varname_S == m_list(F_id)%v(ii)%vn(1) .and. &
              (varname2_S == m_list(F_id)%v(ii)%vn(2) .or. varname_S == ' ')) then
            F_index = ii
            exit
         endif
      enddo
      !------------------------------------------------------------------
      return
   end function incfg_varindex


   !/@*
   function incfg_isvarstep(F_id,F_index,F_istep) result(F_istat)
      implicit none
      integer,intent(in) :: F_id,F_index,F_istep
      integer :: F_istat
      !*@/
      integer :: delta,mo0,mo1,yy0,yy1
      integer(INT64) :: jdatev, step_8, dt_8
      !------------------------------------------------------------------
      F_istat = priv_check_idx(F_id,F_index)
      if (.not.RMN_IS_OK(F_istat)) return

      F_istat = RMN_ERR
      select case(m_list(F_id)%v(F_index)%freq_type(1:2))
      case('st') !STEP
         delta = F_istep
      case('mo') !MONTH
         yy0 = jdate_year(m_list(F_id)%dateo)
         mo0 = jdate_month(m_list(F_id)%dateo)
         step_8 = F_istep
         dt_8   = m_list(F_id)%dt
         jdatev = m_list(F_id)%dateo + step_8 * dt_8
         yy1 = jdate_year(jdatev)
         mo1 = jdate_month(jdatev)
         delta = (mo1-mo0) + (yy1-yy0)*12
      end select

      if (delta == m_list(F_id)%v(F_index)%freq_t0) then
         F_istat = RMN_OK
      else if (m_list(F_id)%v(F_index)%freq_dt > 0 .and. &
           delta >= m_list(F_id)%v(F_index)%freq_t0 .and. &
           delta <= m_list(F_id)%v(F_index)%freq_t1) then
         delta = delta - m_list(F_id)%v(F_index)%freq_t0
         if (mod(delta,m_list(F_id)%v(F_index)%freq_dt) == 0) F_istat = RMN_OK
      endif
      !------------------------------------------------------------------
      return
   end function incfg_isvarstep


   !/@*
   function incfg_meta(F_id,F_index,F_varname_S,F_varname2_S,F_files_S, &
        F_h_int_S,F_v_int_S,F_t_int_S,F_lvl_type_S,F_lvl0,F_lvl1,F_needonce_L, &
        F_ip1list,F_nip1,F_typvar_S,F_mandatory,F_vmin,F_vmax,F_cat0,F_cat1) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments 
      integer,intent(in) :: F_id,F_index
      character(len=*),intent(inout) :: F_varname_S
      character(len=*),intent(inout),optional :: F_varname2_S
      character(len=*),intent(inout),optional :: F_files_S(:)
      character(len=*),intent(inout),optional :: F_h_int_S,F_v_int_S,F_t_int_S, &
           F_lvl_type_S,F_typvar_S
      integer,intent(out),optional :: F_lvl0,F_lvl1,F_ip1list(:),F_nip1,F_mandatory,F_cat0,F_cat1
      real,intent(out),optional :: F_vmin, F_vmax
      logical,intent(out),optional :: F_needonce_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer :: nfiles,mylen,ii
      !----------------------------------------------------------------------
      F_varname_S = ' '
      F_istat = priv_check_idx(F_id,F_index)
      if (.not.RMN_IS_OK(F_istat)) return

      mylen = min(len_trim(m_list(F_id)%v(F_index)%vn(1)),len(F_varname_S))
      F_varname_S(1:mylen) = trim(m_list(F_id)%v(F_index)%vn(1))

      if (present(F_varname2_S)) then
         F_varname2_S = ' '
         mylen = min(len_trim(m_list(F_id)%v(F_index)%vn(2)),len(F_varname2_S))
         F_varname2_S(1:mylen) = trim(m_list(F_id)%v(F_index)%vn(2))
      endif

      if (present(F_files_S)) then
         F_files_S(:) = ' '
         nfiles = min(NFILE_MAX,size(F_files_S))
         do ii=1,nfiles
            mylen = min(len_trim(m_list(F_id)%v(F_index)%file(ii)), &
                 len(F_files_S(ii)))
            F_files_S(ii)(1:mylen) = trim(m_list(F_id)%v(F_index)%file(ii))
         enddo
      endif

      if (present(F_typvar_S)) then
         F_typvar_S = m_list(F_id)%v(F_index)%typvar
      endif

      if (present(F_h_int_S)) then
         F_h_int_S = ' '
         mylen = min(len_trim(m_list(F_id)%v(F_index)%h_int),len(F_h_int_S))
         F_h_int_S(1:mylen) = trim(m_list(F_id)%v(F_index)%h_int)
      endif

      if (present(F_v_int_S)) then
         F_v_int_S = ' '
         mylen = min(len_trim(m_list(F_id)%v(F_index)%v_int),len(F_v_int_S))
         F_v_int_S(1:mylen) = trim(m_list(F_id)%v(F_index)%v_int)
      endif

      if (present(F_t_int_S)) then
         F_t_int_S = ' '
         mylen = min(len_trim(m_list(F_id)%v(F_index)%t_int),len(F_t_int_S))
         F_t_int_S(1:mylen) = trim(m_list(F_id)%v(F_index)%t_int)
      endif

      if (present(F_lvl_type_S)) then
         F_lvl_type_S = ' '
         mylen = min(len_trim(m_list(F_id)%v(F_index)%lvl_type),len(F_lvl_type_S))
         F_lvl_type_S(1:mylen) = trim(m_list(F_id)%v(F_index)%lvl_type)
      endif

      if (present(F_lvl0)) F_lvl0 = m_list(F_id)%v(F_index)%lvl0
      if (present(F_lvl1)) F_lvl1 = m_list(F_id)%v(F_index)%lvl1

      if (present(F_cat0)) F_cat0 = m_list(F_id)%v(F_index)%cat0
      if (present(F_cat1)) F_cat1 = m_list(F_id)%v(F_index)%cat1

      if (present(F_nip1) .and. present(F_ip1list)) then
         F_nip1 = priv_get_ip1list(F_id,F_index,F_ip1list)
      else if (present(F_nip1)) then
         F_nip1 = priv_get_ip1list(F_id,F_index)
      else if (present(F_ip1list)) then
         F_ip1list = -1
         ii = priv_get_ip1list(F_id,F_index,F_ip1list)
      endif

      if (present(F_vmin)) F_vmin = m_list(F_id)%v(F_index)%vmin
      if (present(F_vmax)) F_vmax = m_list(F_id)%v(F_index)%vmax

      if (present(F_needonce_L)) F_needonce_L = &
           (m_list(F_id)%v(F_index)%freq_dt == 0 .or. &
           m_list(F_id)%v(F_index)%freq_t1 <= m_list(F_id)%v(F_index)%freq_t0)

      if (present(F_mandatory)) F_mandatory = m_list(F_id)%v(F_index)%mandatory
      !----------------------------------------------------------------------
      return
   end function incfg_meta


   !/@*
   function incfg_time_4(F_id,F_step,F_datev,F_dateo,F_dt) result(F_istat)
      implicit none
      !@objective
      !@arguments 
      integer,intent(in) :: F_id,F_step
      integer,intent(out) :: F_datev,F_dateo,F_dt
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer(INT64) :: jdatev,jdateo
      !----------------------------------------------------------------------
      F_datev = RMN_ERR
      F_dateo = RMN_ERR
      F_istat =  incfg_time_8(F_id,F_step,jdatev,jdateo,F_dt)
      if (RMN_IS_OK(F_istat)) then
         F_datev = jdate_to_cmc(jdatev)
         F_dateo = jdate_to_cmc(jdateo)
      endif
      !----------------------------------------------------------------------
      return
   end function incfg_time_4


   !/@*
   function incfg_time_8(F_id,F_step,F_datev,F_dateo,F_dt) result(F_istat)
      implicit none
      !@objective
      !@arguments 
      integer,intent(in) :: F_id,F_step
      integer(INT64),intent(out) :: F_datev,F_dateo
      integer,intent(out) :: F_dt
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer(INT64) :: step_8, dt_8
      !----------------------------------------------------------------------
      F_istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(F_istat)) return
 
      step_8 = F_step
      dt_8   = m_list(F_id)%dt
      F_datev = m_list(F_id)%dateo + step_8 * dt_8
      F_dateo = m_list(F_id)%dateo
      F_dt = m_list(F_id)%dt
      !----------------------------------------------------------------------
      return
   end function incfg_time_8

   !==== Private Functions =================================================


   !/@*
   function priv_check_idx(F_id,F_index)  result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      integer,intent(in),optional :: F_index
      !@return
      integer :: F_istat
      !*@/
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_id < 1 .or. F_id > m_nlist) then
         call msg(MSG_WARNING,'(incfg) invalid incfg id')
         return
      endif
      if (present(F_index)) then
         if (F_index < 1 .or. F_index > m_list(F_id)%n .or. &
              F_index > size(m_list(F_id)%v)) then
            call msg(MSG_WARNING,'(incfg) index out of range')
            return
         endif
      endif
      F_istat = RMN_OK
      !----------------------------------------------------------------------
      return
   end function priv_check_idx


   !/@*
   function priv_parse_cfgfile(F_id,F_filename_S) result(F_nvar)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_filename_S !full/rel path of config file
      !@return
      integer :: F_nvar
      !@author Stephane Chamberland, 2011-07
      !*@/
      integer :: istat,fileid
      character(len=1024) :: string_S
      !---------------------------------------------------------------------
      write(string_S, '(a,1x,i0,1x,a)') 'priv_parse_cfgfile:',F_id,trim(F_filename_S)
      call msg(MSG_DEBUG,string_S)

      F_nvar = RMN_ERR

      istat = clib_isfile(trim(F_filename_S))
      istat = min(clib_isreadok(trim(F_filename_S)),istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(incfg) File not found or not readable: '//trim(F_filename_S))
         return
      endif
 
      fileid = 0
      istat = fnom(fileid,F_filename_S,'SEQ/FMT+R/O+OLD',0)
      if (.not.RMN_IS_OK(istat) .or. fileid <= 0) then
         call msg(MSG_WARNING,'(incfg) Problem opening file: '//trim(F_filename_S))
         return
      endif

      F_nvar = 0
      DOFILE: do
         read(fileid,'(a)',iostat=istat) string_S
         if (istat /=0 .or. F_nvar+1 > NVAR_MAX) exit DOFILE
         call str_tab2space(string_S)
         string_S = adjustl(string_S)
         if (string_S(1:1)/='#' .and. string_S/=' ') then
            istat = incfg_add_string(F_id,string_S)
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
   subroutine priv_parse_file(F_id,F_index,F_string_S)
      implicit none
      character(len=*),intent(in) :: F_string_S
      integer,intent(in) :: F_id,F_index
      !*@/
      integer :: ii
      character(len=1024) :: parts_S(NFILE_MAX+1)
      !------------------------------------------------------------------
      call str_split2list(parts_S,F_string_S,',',NFILE_MAX+1)
      ii = 1
      ITEMLOOP: do
         if (ii == NFILE_MAX+1) exit ITEMLOOP
         if (parts_S(ii) == ' ') exit ITEMLOOP
         m_list(F_id)%v(F_index)%file(ii) = parts_S(ii)
         ii = ii + 1
      enddo ITEMLOOP
      if (parts_S(NFILE_MAX+1) /= ' ') &
           call msg(MSG_WARNING,'(incfg) ignoring extra files (max of 4): '// &
           trim(parts_S(NFILE_MAX+1)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_file


   !/@*
   subroutine priv_parse_freq(F_id,F_index,F_string_S)
      implicit none
      character(len=*),intent(in) :: F_string_S
      integer,intent(in) :: F_id,F_index
      !*@/
      integer,parameter :: NMAX = 5
      character(len=1024) :: parts_S(NMAX)
      integer :: ii,val,istat
      real :: rfact,rt0,rdt,rt1
      !------------------------------------------------------------------
      call str_split2list(parts_S,F_string_S,',',NMAX)
      rfact = 1.
      m_list(F_id)%v(F_index)%freq_type = 'error'
      ii = 2
      select case(parts_S(1)(1:2))
      case('st') !STEP
         m_list(F_id)%v(F_index)%freq_type = 'step'//'       '
      case('mo') !MONTH
         m_list(F_id)%v(F_index)%freq_type = 'month'//'       '
      case('we') !WEEK
         m_list(F_id)%v(F_index)%freq_type = 'step'//'       '
         rfact = real(60*60*24*7)/real(m_list(F_id)%dt)
      case('da') !DAY
         m_list(F_id)%v(F_index)%freq_type = 'step'//'       '
         rfact = real(60*60*24)/real(m_list(F_id)%dt)
      case('ho') !HOURS
         m_list(F_id)%v(F_index)%freq_type = 'step'//'       '
         rfact = real(60*60)/real(m_list(F_id)%dt)
      case('mi') !MINUTES
         m_list(F_id)%v(F_index)%freq_type = 'step'//'       '
         rfact = real(60)/real(m_list(F_id)%dt)
      case('  ')
         return
      case default
         m_list(F_id)%v(F_index)%freq_type = 'step'//'       '
         ii = 1
      end select

      read(parts_S(ii),fmt=*,iostat=istat) val
      if (istat /= 0 .or. parts_S(ii) == ' ') then
         m_list(F_id)%v(F_index)%freq_type = 'error'//'       '
         call msg(MSG_WARNING,'(incfg) ignoring invalid freq: '//trim(F_string_S))
         return
      endif
      m_list(F_id)%v(F_index)%freq_t0 = max(0,val)

      ii = ii + 1
      m_list(F_id)%v(F_index)%freq_dt = 0
      if (parts_S(ii) /= ' ') then
         read(parts_S(ii),fmt=*,iostat=istat) val
         if (istat /= 0) then
            m_list(F_id)%v(F_index)%freq_type = 'error'//'       '
            call msg(MSG_WARNING,'(incfg) ignoring invalid freq: '//trim(F_string_S))
            return
         endif
         m_list(F_id)%v(F_index)%freq_dt = max(0,val)
      endif

      ii = ii + 1
      m_list(F_id)%v(F_index)%freq_t1 = m_list(F_id)%v(F_index)%freq_t0
      if (m_list(F_id)%v(F_index)%freq_dt > 0) then
         m_list(F_id)%v(F_index)%freq_t1 = -1
         if (parts_S(ii) /= ' ') then
            read(parts_S(ii),fmt=*,iostat=istat) val
            if (istat /= 0) then
               m_list(F_id)%v(F_index)%freq_type = 'error'//'       '
               call msg(MSG_WARNING,'(incfg) ignoring invalid freq: '//trim(F_string_S))
               return
            endif
            m_list(F_id)%v(F_index)%freq_t1 = max(0,val)
         endif
      endif

      if (abs(1. - rfact) > EPSILON) then
         rdt = real(m_list(F_id)%v(F_index)%freq_dt) * rfact
         rt0 = real(m_list(F_id)%v(F_index)%freq_t0) * rfact
         if (abs(rt0 - real(nint(rt0))) > EPSILON .or. &
              abs(rdt - real(nint(rdt))) > EPSILON) then
            m_list(F_id)%v(F_index)%freq_type = 'error'//'       '
            call msg(MSG_WARNING,'(incfg) cannot satisfy frequency: '//trim(F_string_S))
            return
         endif
         m_list(F_id)%v(F_index)%freq_t0 = nint(rt0)
         m_list(F_id)%v(F_index)%freq_dt = nint(rdt)
         if (m_list(F_id)%v(F_index)%freq_t1 > 0) then
            rt1 = floor(real(m_list(F_id)%v(F_index)%freq_t1) * rfact)
            m_list(F_id)%v(F_index)%freq_t1 = nint(max(rt0,rt1))
         endif
      endif
      if (m_list(F_id)%v(F_index)%freq_t1 < 0) then
         m_list(F_id)%v(F_index)%freq_t1 = huge(val)
      end if

      if (parts_S(NMAX) /= ' ') &
           call msg(MSG_WARNING,'(incfg) ignoring extra freq params: '//trim(parts_S(NMAX)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_freq


   !/@*
   subroutine priv_parse_level(F_id,F_index,F_string_S)
      implicit none
      character(len=*),intent(in) :: F_string_S
      integer,intent(in) :: F_id,F_index
      !*@/
      integer,parameter :: NMAX = 3
      character(len=1024) :: parts_S(NMAX),tmp_S
      integer :: val,istat,n
      !------------------------------------------------------------------
      call str_split2list(parts_S,F_string_S,',',NMAX)
 
      m_list(F_id)%v(F_index)%lvl0 = -1
      m_list(F_id)%v(F_index)%lvl1 = -1

      n = 1
      tmp_S = parts_S(n)
      select case(tmp_S(1:3))
      case('arb') !arbitrary
         m_list(F_id)%v(F_index)%lvl_type = KNOWN_LEVEL_TYPE(LEVEL_ARBITRARY)
         return
      case('tdi') !diagnostic themo
         m_list(F_id)%v(F_index)%lvl_type = KNOWN_LEVEL_TYPE(LEVEL_TDIAG)
         m_list(F_id)%v(F_index)%lvl0 = 0
         m_list(F_id)%v(F_index)%lvl1 = 0
         return
      case('mdi') !diagnostic momentum
         m_list(F_id)%v(F_index)%lvl_type = KNOWN_LEVEL_TYPE(LEVEL_MDIAG)
         m_list(F_id)%v(F_index)%lvl0 = 0
         m_list(F_id)%v(F_index)%lvl1 = 0
         return
!!$      case('mom') !momentum
!!$         m_list(F_id)%v(F_index)%lvl_type = 'momentum'
!!$         return
!!$      case('the') !thermodynamic
!!$         m_list(F_id)%v(F_index)%lvl_type = 'thermodynamic'
!!$         return
!!$      case default ! skip
!!$         call msg(MSG_WARNING,'(incfg) ignoring invalid level: '//trim(F_string_S))
      end select

      read(parts_S(n),fmt=*,iostat=istat) val
      if (istat /= 0 .or. parts_S(n) == ' ') then
         call msg(MSG_WARNING,'(incfg) ignoring invalid level: '//trim(F_string_S))
         return
      endif
      m_list(F_id)%v(F_index)%lvl0 = max(-1,val) !#max(0,val)
      n = n+1

      m_list(F_id)%v(F_index)%lvl1 = m_list(F_id)%v(F_index)%lvl0
      if (val >= 0) then
         read(parts_S(n),fmt=*,iostat=istat) val
         if (istat == 0 .and. parts_S(n) /= ' ') then
            m_list(F_id)%v(F_index)%lvl1 = max(0,val)
         endif
         n = n+1
      endif

      if (parts_S(n) /= ' ') &
           call msg(MSG_WARNING,'(incfg) ignoring extra level params: '//trim(parts_S(n)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_level


   !/@*
   subroutine priv_parse_cat(F_id,F_index,F_string_S)
      implicit none
      character(len=*),intent(in) :: F_string_S
      integer,intent(in) :: F_id,F_index
      !*@/
      integer,parameter :: NMAX = 3
      character(len=1024) :: parts_S(NMAX)
      integer :: val,istat,n
      !------------------------------------------------------------------
      call str_split2list(parts_S,F_string_S,',',NMAX)
 
      m_list(F_id)%v(F_index)%cat0 = -1
      m_list(F_id)%v(F_index)%cat1 = -1

      n = 1
      read(parts_S(n),fmt=*,iostat=istat) val
      if (istat /= 0 .or. parts_S(n) == ' ') then
         call msg(MSG_WARNING,'(incfg) ignoring invalid cat: '//trim(F_string_S))
         return
      endif
      m_list(F_id)%v(F_index)%cat0 = max(-1,val) !#max(0,val)
      n = n+1
      m_list(F_id)%v(F_index)%cat1 = m_list(F_id)%v(F_index)%cat0
      if (val >= 0) then
         read(parts_S(n),fmt=*,iostat=istat) val
         if (istat == 0 .and. parts_S(n) /= ' ') then
            m_list(F_id)%v(F_index)%cat1 = max(0,val)
         endif
         n = n+1
      endif

      if (parts_S(n) /= ' ') &
           call msg(MSG_WARNING,'(incfg) ignoring extra cat params: '//trim(parts_S(n)))
      !------------------------------------------------------------------
      return
   end subroutine priv_parse_cat


   !/@*
   function priv_set_ip1list(F_id,F_ip1list) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_ip1list(:)
      !@return
      integer :: F_istat
      !*@/
      integer :: istat
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(incfg) set_ip1 [BEGIN]')
      F_istat = priv_check_idx(F_id)
      if (.not.RMN_IS_OK(F_istat)) return
      allocate(m_list(F_id)%ip1list(size(F_ip1list)),stat=istat)
      m_list(F_id)%ip1list(:) = F_ip1list(:)
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(incfg) set_ip1 [END]')
      !----------------------------------------------------------------------
      return
   end function priv_set_ip1list


   !/@*
   function priv_get_ip1list(F_id,F_index,F_ip1list) result(F_nip1)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_index
      integer,intent(out),optional :: F_ip1list(:)
      !@return
      integer :: F_nip1
      !*@/
      integer :: lvl0,lvl1,k,nmax,ipkind
      real :: zp1
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(incfg) get_ip1 [BEGIN]')
      F_nip1 = RMN_ERR
 
      lvl0 = -1
      lvl1 = -1
      if (F_index > 0) then
         lvl0 = m_list(F_id)%v(F_index)%lvl0
         lvl1 = m_list(F_id)%v(F_index)%lvl1
      endif

      nmax = 999999
      F_nip1 = 1
      if (present(F_ip1list)) then
         nmax = size(F_ip1list)
         F_ip1list(1) = 0
      endif
      if (lvl0 < 0) then
         if (associated(m_list(F_id)%ip1list)) then
            F_nip1 = min(nmax,size(m_list(F_id)%ip1list))
            if (present(F_ip1list)) &
                 F_ip1list(1:F_nip1) = m_list(F_id)%ip1list(1:F_nip1)
         else
            F_nip1 = 1
            F_ip1list(1) = -1
         endif
      else if (lvl0 > 0) then
         ipkind = RMN_CONV_ARBITRARY
         F_nip1 = 0
         do k=lvl0,lvl1
            F_nip1 = F_nip1 + 1
            if (F_nip1 > nmax) cycle
            if (present(F_ip1list)) then
               zp1 = real(k)
               call convip_plus(F_ip1list(F_nip1),zp1,ipkind,RMN_CONV_P2IPNEW,' ', &
                    .not.RMN_CONV_USEFORMAT_L)
            endif
         enddo
      endif
      call msg(MSG_DEBUG,'(incfg) get_ip1 [END]')
      !----------------------------------------------------------------------
      return
   end function priv_get_ip1list

end module incfg_mod
