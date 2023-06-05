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
module outcfg_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_isfile, clib_tolower, clib_isreadok
   use str_mod, only: str_tab2space,str_rm_quotes,str_toint,str_tobool,str_toreal
   use sort_mod
   use mu_jdate_mod
   implicit none
   private
   !@objective
   !@author Stephane Chamberland, 2011-09
   !@description
   ! Public functions
   public :: outcfg_new, outcfg_new_4, outcfg_new_8, outcfg_getlist, &
        outcfg_time, outcfg_meta,outcfg_get_idxlist,outcfg_var_meta, &
        outcfg_level_meta,outcfg_grid_meta
   ! Public constants
   integer,parameter,public :: OUTCFG_NCFG_MAX = 8
   integer,parameter,public :: OUTCFG_STRLEN = 16
   character(len=1),parameter,public :: OUTCFG_LVLTYPE(2) = (/'p','m'/)
   integer,parameter,public :: OUTCFG_UNITS_STEP = 1
   integer,parameter,public :: OUTCFG_UNITS_1STOF_MONTH = 2
   ! File format
   !   one output command per line
   !   empty lines and lines starting with # are ignored
   !   spaces are ignored
   !   case is ignored
   !   recognized keys
   !   * closeopen : time before changing file (default=step,1)
   !                 accepted format:
   !                   closeopen = UNITS, INTERVAL
   !                     UNITS   : SECONDS, MINUTES, HOURS, DAYS, STEPS, MONTH_1STOF (optional, default=steps)
   !                     INTERVAL: UNITS multiplier [int] (ignored is UNITS = MONTH_1STOF)
   !                 closeopen is converted to steps and thus will be rounded
   !   * debug_L   : print extra debug info [T, F] (default=false)
   !   * datyp     : RPN STD default data type for compression [int]
   !                 (override with xnbits) (default=0)
   !   * etk_S     : etiket [string] (default='')
   !   * filt_pass : number of pass for horizontal filtering [int] (default=0)
   !                 (override with filtre)
   !   * filt_coef : coefficient for horizontal filtering [float] (default=0.)
   !                 (override with filtre)
   !                 accepted format: closeopen = UNITS, INTERVAL
   !   * filtre    : specific var filtering override
   !                 accepted formats:
   !                   filtre([VARLIST],pass,NPASS,coef,MYCOEF)
   !                   filtre([VARLIST],coef,MYCOEF)
   !                   filtre([VARLIST],pass,NPASS)
   !                 VARLIST : comma separated list of var
   !                 NPASS   : override filt_pass [int]
   !                 MYCOEF  : override filt_coef [float]
   !   * flip_L    : output vertical planes instead of horizontal [T, F]
   !                 (default=false)
   !   * fullplane_L : output on the full horizontal domain
   !                   instead of block sub-domain [T, F] (default=false)
   !   * grid      : define a named list of output sub-domain
   !                 accepted format:
   !                   grid = NAME, model
   !                   grid = NAME, core
   !                   grid = NAME, reduc, ????
   !                 TODO
   !   * ip3       : ip3 value [int] (default=0)
   !   * levels    : define a named list of output levels
   !                 accepted format:
   !                   levels = NAME, TODO??? press, cubic, linear, all, -1, surface, diag, 0, [LIST], <lvl0,lvlend,interval>
   !                 TODO
   !   * linbot    : number of bottom levels to vertically interpolate linearly [int]
   !   * nbits     : RPN STD default nb of bits for compression [int <= 32]
   !                 (override with xnbits) (default=16)
   !   * ndigits   : TODO [int] (default=0)
   !   * rewrite_L : overwrite existing field in file [T, F] (default=F)
   !   * units_S   : Fxc time encoding unit for filename extentions (default=hour)
   !                 accepted format:
   !                   units_S = UNITS
   !                     UNITS   : SECONDS, MINUTES, HOURS, DAYS, STEPS, MONTH
   !                 if UNIT == month, this trigger closeopen=MONTH_1STOF
   !   * vinterp_S : type of default vertical interpolation [cubic, linear] (deafult=cubic)
   !   * sortie_?  : request output for specific var list at specified sub-grid, levels and steps
   !                 accepted format:
   !                   sortie_TAG([VARLIST],grid,MYGRID,levels,MYLEVELS,steps,MYSTEPS)
   !                     TAG     : model component's tag (e.g. "p", "d", ...)
   !                     VARLIST : comma separated list of var to output/write
   !                     MYGRID  : previously defined output grid name/tag
   !                     MYLEVELS: previously defined output level set name/tag
   !                     MYSTEPS : previously defined output steps set name/tag
   !   * steps     : define a named list of output interval/steps
   !                 accepted formats:
   !                   steps = NAME,UNITS,LIST
   !                   steps = NAME,UNITS,[LIST]
   !                   steps = NAME,UNITS,<STEP0,STEPEND,INTERVAL>
   !                   steps = NAME,init,UNITS,LIST
   !                     NAME  : tag/name of the list
   !                     UNITS : units of givent list (optional, default=steps)
   !                             accepted values: HOURS, STEPS
   !                     LIST  : list steps/hours
   !                     STEP0,STEPEND,INTERVAL: loop from step0 to STEPEND with strides of INTERVAL
   !                     init  : ?TODO?
   !                 Note that steps are converted to steps and
   !                 might thus be rounded if HOURS is used.
   !   * xnbits    : specific var compression nbits and datyp override
   !                 accepted formats:
   !                   xnbits([VARLIST],bits,MYNBITS,datyp,MYDATYP)
   !                   xnbits([VARLIST],bits,MYNBITS)
   !                   xnbits([VARLIST],datyp,MYDATYP)
   !                 VARLIST : comma separated list of var
   !                 MYNBITS : override nbits [int]
   !                 MYDATYP : override datyp [int]
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   real,parameter :: EPSILON = 1.e-5
   integer,parameter :: NSETS_MAX = 32
   integer,parameter :: NLIST_MAX = 128
   integer,parameter :: NLEVELS_MAX = 32768
   integer,parameter :: LVL_TYPE_PRES = 1
   integer,parameter :: LVL_TYPE_MODEL = 2
   integer,parameter :: NBITS_MIN = 8
   integer,parameter :: NBITS_MAX = 32
   real,parameter :: NO_VAL = -2.
   real,parameter :: ALL_VAL = -1.
   real,parameter :: SURF_VAL = 0.

   type :: outcfg_levels_T
      character(len=OUTCFG_STRLEN) :: id_S
      character(len=OUTCFG_STRLEN) :: vint_S !vertical interpolation type
      integer :: type,linbot !type of levels, nb.lvls to force lin-vert-int
      integer :: n           !number of v elements
      real,pointer :: v(:)   !level list (pres or 1..nk number)
   end type outcfg_levels_T

   type :: outcfg_steps_T
      character(len=OUTCFG_STRLEN) :: id_S
      logical :: init_L
      integer :: n            !number of v elements
      integer,pointer :: v(:) !when to output in seconds
      integer :: l(3)         !loop i0,in,di (when to output in seconds)
   end type outcfg_steps_T

   type :: outcfg_grid_T
      character(len=OUTCFG_STRLEN) :: id_S
      character(len=OUTCFG_STRLEN) :: tag_S
      integer :: ij0n(5) !i0,in,j0,nj,dd
   end type outcfg_grid_T

   type :: outcfg_sortie_T
      integer :: n            !number of ON_S elements
      character(len=OUTCFG_STRLEN),pointer :: on_S(:)
      integer :: l_idx,s_idx,g_idx
   end type outcfg_sortie_T

   type :: outcfg_xnbits_T
      integer :: n            !number of ON_S elements
      character(len=OUTCFG_STRLEN),pointer :: n_S(:)
      integer :: nbits
      integer :: datyp
   end type outcfg_xnbits_T

   type :: outcfg_filtre_T
      integer :: n            !number of ON_S elements
      character(len=OUTCFG_STRLEN),pointer :: n_S(:)
      integer :: npass
      real :: coef
   end type outcfg_filtre_T

   type :: outcfg_T
      character(len=OUTCFG_STRLEN) :: tag_S,etk_S,vinterp_s,units_s
      integer(INT64) :: dateo  !Initial date [jdate sec]
      integer :: dt         !time step length [seconds]
      integer :: ip3, datyp, nbits,ndigits,linbot,filt_pass,closeopen,closeopen_units
      integer :: reduc_core(5)
      real    :: filt_coef
      logical :: debug_L, flip_L, rewrite_L, fullplane_L
      integer :: nl,ns,ng,no,nb,nf
      type(outcfg_levels_T),pointer :: l(:)
      type(outcfg_steps_T),pointer :: s(:)
      type(outcfg_grid_T),pointer :: g(:)
      type(outcfg_sortie_T),pointer :: o(:)
      type(outcfg_xnbits_T),pointer :: b(:)
      type(outcfg_filtre_T),pointer :: f(:)
   end type outcfg_T

   character(len=OUTCFG_STRLEN),parameter :: KNOWN_V_INT(3) = (/ &
        'cubic           ','linear          ','nearest         ' &
        /)
   character(len=OUTCFG_STRLEN),parameter :: KNOWN_T_UNITS(6) = (/ &
        'steps           ','seconds         ','minutes         ', &
        'hours           ','days            ','months          '/)

   character(len=OUTCFG_STRLEN),parameter :: &
        TAG_START = '<outcfg>', &
        TAG_END   = '</outcfg>'

   type(outcfg_T),save :: OUTCFG_VAR_DEFAULT

   type(outcfg_T),pointer,save :: m_cfgs(:) => null()
   integer,save :: m_ncfgs = 0

   interface outcfg_new
      module procedure outcfg_new_4
      module procedure outcfg_new_8
   end interface outcfg_new

   interface outcfg_time
      module procedure outcfg_time_4
      module procedure outcfg_time_8
   end interface outcfg_time

contains

   !/@*
   function outcfg_new_4(F_filename_S,F_tag_S,F_dateo,F_dt,F_reduc_core,F_intags_L) result(F_id)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_filename_S !full/rel path of config file
      character(len=*),intent(in) :: F_tag_S
      integer,intent(in) :: F_dateo,F_dt
      integer,intent(in),optional :: F_reduc_core(:)
      logical,intent(in),optional ::  F_intags_L
      !@return
      integer :: F_id !config id
      !@author Stephane Chamberland, 2011-09
      !*@/
      integer(INT64) :: jdateo
      logical ::  intags_L
      !----------------------------------------------------------------------
      intags_L = .false.
      if (present(F_intags_L)) intags_L = F_intags_L
      jdateo = jdate_from_cmc(F_dateo)
      if (present(F_reduc_core)) then
         F_id = outcfg_new_8(F_filename_S,F_tag_S,jdateo,F_dt,F_reduc_core,intags_L)
      else
         F_id = outcfg_new_8(F_filename_S,F_tag_S,jdateo,F_dt,F_intags_L=intags_L)
      endif
      !----------------------------------------------------------------------
      return
   end function outcfg_new_4


   !/@*
   function outcfg_new_8(F_filename_S,F_tag_S,F_dateo,F_dt,F_reduc_core,F_intags_L) result(F_id)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_filename_S !full/rel path of config file
      character(len=*),intent(in) :: F_tag_S
      integer(INT64),intent(in) :: F_dateo
      integer,intent(in) :: F_dt
      integer,intent(in),optional :: F_reduc_core(:)
      logical,intent(in),optional ::  F_intags_L
      !@return
      integer :: F_id !config id
      !@author Stephane Chamberland, 2011-09
      !*@/
      logical,save :: is_init_L = .false.
      logical :: intags_L
      integer :: istat,nvar,n
      character(len=1024) :: string_S
      !----------------------------------------------------------------------
      write(string_S, *) '(outcfg_new) ',trim(F_tag_S),F_dateo,F_dt,trim(F_filename_S)
      call msg(MSG_DEBUG,string_S)

      F_id = RMN_ERR

      if (.not.is_init_L) then
         allocate(m_cfgs(OUTCFG_NCFG_MAX),stat=istat)
         if (istat /= 0) return
         OUTCFG_VAR_DEFAULT%tag_S = ' '
         OUTCFG_VAR_DEFAULT%etk_S = ' '
         OUTCFG_VAR_DEFAULT%vinterp_s = 'cubic'
         OUTCFG_VAR_DEFAULT%units_s = ' '
         OUTCFG_VAR_DEFAULT%dateo = 0
         OUTCFG_VAR_DEFAULT%dt = 0
         OUTCFG_VAR_DEFAULT%ip3 = 0
         OUTCFG_VAR_DEFAULT%datyp = 0
         OUTCFG_VAR_DEFAULT%nbits = 16
         OUTCFG_VAR_DEFAULT%ndigits = 0
         OUTCFG_VAR_DEFAULT%linbot = 0
         OUTCFG_VAR_DEFAULT%closeopen = 1
         OUTCFG_VAR_DEFAULT%closeopen_units = OUTCFG_UNITS_STEP
         OUTCFG_VAR_DEFAULT%filt_pass = 0
         OUTCFG_VAR_DEFAULT%filt_coef = 0.
         OUTCFG_VAR_DEFAULT%debug_L = .false.
         OUTCFG_VAR_DEFAULT%flip_L = .false.
         OUTCFG_VAR_DEFAULT%rewrite_L = .false.
         OUTCFG_VAR_DEFAULT%fullplane_L = .false.
         OUTCFG_VAR_DEFAULT%reduc_core = (/1,-1,1,-1,1/)
         OUTCFG_VAR_DEFAULT%nl = 0
         OUTCFG_VAR_DEFAULT%ns = 0
         OUTCFG_VAR_DEFAULT%ng = 0
         OUTCFG_VAR_DEFAULT%no = 0
         OUTCFG_VAR_DEFAULT%nb = 0
         OUTCFG_VAR_DEFAULT%nf = 0
         nullify(OUTCFG_VAR_DEFAULT%l)
         nullify(OUTCFG_VAR_DEFAULT%s)
         nullify(OUTCFG_VAR_DEFAULT%g)
         nullify(OUTCFG_VAR_DEFAULT%o)
         nullify(OUTCFG_VAR_DEFAULT%b)
         nullify(OUTCFG_VAR_DEFAULT%f)
         is_init_L = .true.
      endif

      if (m_ncfgs == OUTCFG_NCFG_MAX) then
         call msg(MSG_WARNING,'(outcfg) Too many opened outcfg.')
         return
      endif

      m_ncfgs = m_ncfgs + 1
      m_cfgs(m_ncfgs) = OUTCFG_VAR_DEFAULT

      allocate( &
           m_cfgs(m_ncfgs)%l(NSETS_MAX), &
           m_cfgs(m_ncfgs)%s(NSETS_MAX), &
           m_cfgs(m_ncfgs)%g(NSETS_MAX), &
           m_cfgs(m_ncfgs)%o(NSETS_MAX), &
           m_cfgs(m_ncfgs)%b(NSETS_MAX), &
           m_cfgs(m_ncfgs)%f(NSETS_MAX), &
           stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(outcfg) Problem allocating memory.')
         return
      endif

      if (F_tag_S /= ' ') m_cfgs(m_ncfgs)%tag_S = '_'//F_tag_S
      m_cfgs(m_ncfgs)%dateo = F_dateo
      m_cfgs(m_ncfgs)%dt = F_dt
      if (present(F_reduc_core)) then
         n = min(size(m_cfgs(m_ncfgs)%reduc_core),size(F_reduc_core))
         m_cfgs(m_ncfgs)%reduc_core(1:n) = F_reduc_core(1:n)
      endif

      intags_L = .false.
      if (present(F_intags_L)) intags_L = F_intags_L

      nvar = RMN_ERR
!!$      if (present(F_filename_S)) then
      istat = priv_parse_cfgfile(m_ncfgs,F_filename_S,intags_L)
         if (.not.RMN_IS_OK(istat)) &
              call msg(MSG_ERROR,'(outcfg) Problem parsing config file: '//trim(F_filename_S))

         if (.not.RMN_IS_OK(istat)) then
            deallocate( &
                 m_cfgs(m_ncfgs)%l, &
                 m_cfgs(m_ncfgs)%s, &
                 m_cfgs(m_ncfgs)%g, &
                 m_cfgs(m_ncfgs)%o, &
                 m_cfgs(m_ncfgs)%b, &
                 m_cfgs(m_ncfgs)%f, &
                 stat=istat)
            m_ncfgs = m_ncfgs - 1
            return
         endif
!!$      endif

      F_id = m_ncfgs

      write(string_S, '(a,1x,i0,1x,i0,1x,a)') '(outcfg_new) ',F_id,m_ncfgs,trim(F_tag_S)//' : '//trim(F_filename_S)
      call msg(MSG_DEBUG,string_S)
      !----------------------------------------------------------------------
      return
   end function outcfg_new_8


   !TODO-later: outcfg_free

   !/@*
   function outcfg_parse_string(F_id, F_string_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_string_S
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-07
      !*@/
      integer, external :: str_split2keyval
      integer, parameter :: NMAX = 32
      integer :: istat
      character(len=2048) :: key_S, val_S, string_S, msg_S
      !----------------------------------------------------------------------
      write(string_S, '(a,i3,a)') '(outcfg_parse_string)', F_id, trim(F_string_S)
      call msg(MSG_DEBUG, string_S)

      F_istat = priv_str2kv(F_string_S, key_S, val_S)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(outcfg) invalid token: '//trim(F_string_S))
         F_istat = RMN_ERR
         return
      endif

      if (key_S == ' ') then
         call msg(MSG_DEBUG,'(outcfg) ignoring comment/empty:'//trim(F_string_S))
         return !ignore comment and empty lines
      endif

      istat = RMN_ERR
      IF_VAL: if (val_S /= ' ') then
         F_istat = RMN_OK
         istat = RMN_OK
         select case(trim(key_S))
         case('etk_s')
            m_cfgs(F_id)%etk_S = val_S
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(m_cfgs(F_id)%etk_S))
         case('ip3')
            istat = str_toint(m_cfgs(F_id)%ip3,val_S)
            write(msg_S,'(i7)') m_cfgs(F_id)%ip3
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('datyp')
            istat = str_toint(m_cfgs(F_id)%datyp,val_S)
            write(msg_S,'(i7)') m_cfgs(F_id)%datyp
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('nbits')
            istat = str_toint(m_cfgs(F_id)%nbits,val_S)
            m_cfgs(F_id)%nbits = max(NBITS_MIN,min(abs(m_cfgs(F_id)%nbits),NBITS_MAX))
            write(msg_S,'(i7)') m_cfgs(F_id)%nbits
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('ndigits')
            istat = str_toint(m_cfgs(F_id)%ndigits,val_S)
            write(msg_S,'(i7)') m_cfgs(F_id)%ndigits
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('linbot')
            istat = str_toint(m_cfgs(F_id)%linbot,val_S)
            write(msg_S,'(i7)') m_cfgs(F_id)%linbot
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('debug_l')
            istat = str_tobool(m_cfgs(F_id)%debug_L,val_S)
            if (m_cfgs(F_id)%debug_L) then
               istat = fstopc('MSGLVL','WARNIN',RMN_OPT_SET)
            else
               istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
            endif
            write(msg_S,'(l)') m_cfgs(F_id)%debug_L
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('flip_l')
            istat = str_tobool(m_cfgs(F_id)%flip_L,val_S)
            write(msg_S,'(l)') m_cfgs(F_id)%flip_L
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('rewrite_l')
            istat = str_tobool(m_cfgs(F_id)%rewrite_L,val_S)
            write(msg_S,'(l)') m_cfgs(F_id)%rewrite_L
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('fullplane_l')
            istat = str_tobool(m_cfgs(F_id)%fullplane_L,val_S)
            write(msg_S,'(l)') m_cfgs(F_id)%fullplane_L
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('filt_pass')
            istat = str_toint(m_cfgs(F_id)%filt_pass,val_S)
            m_cfgs(F_id)%filt_pass = max(0,min(abs(m_cfgs(F_id)%filt_pass),1024))
            write(msg_S,'(i7)') m_cfgs(F_id)%filt_pass
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('filt_coef')
            istat = str_toreal(m_cfgs(F_id)%filt_coef,val_S)
            write(msg_S, '(f0.3)') m_cfgs(F_id)%filt_coef
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('vinterp_s')
            istat = priv_str2str(m_cfgs(F_id)%vinterp_S,val_S,KNOWN_V_INT,1)
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(m_cfgs(F_id)%vinterp_S))
         case('units_s')
            istat = priv_str2str(m_cfgs(F_id)%units_S,val_S,KNOWN_T_UNITS,2)
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(m_cfgs(F_id)%units_S))
         case('closeopen')
            istat = priv_parse_closeopen(m_cfgs(F_id)%closeopen,m_cfgs(F_id)%closeopen_units,val_S,m_cfgs(F_id)%dt)
            if (m_cfgs(F_id)%closeopen_units == OUTCFG_UNITS_STEP) then
               write(msg_S,'(i7,a,i7)') m_cfgs(F_id)%closeopen,' * dt: units=',m_cfgs(F_id)%closeopen_units
            else
               write(msg_S,'(i7,a,i7)') m_cfgs(F_id)%closeopen,': units 1ST OF MONTH ',m_cfgs(F_id)%closeopen_units
            endif
            call msg(MSG_INFO,'(outcfg) '//trim(key_S)//': '//trim(msg_S))
         case('levels')
            istat = priv_parse_levels(F_id,val_S)
         case('grid')
            istat = priv_parse_grid(F_id,val_S)
         case('steps')
            istat = priv_parse_steps(F_id,val_S)

         case('xnbit')
            istat = priv_parse_xnbits(F_id,val_S)
         case('filtre')
            istat = priv_parse_filtre(F_id,val_S)
         case default
            if (key_S == 'sortie'//trim(m_cfgs(F_id)%tag_S)) then
               istat = priv_parse_sortie(F_id,val_S)
            else if (key_S(1:6) == 'sortie') then
               !TODO: use a caller's defined default (like sortie_d)
               call msg(MSG_DEBUG,'(outcfg) ignoring token: '//trim(key_S)//' [accepting only: sortie'//trim(m_cfgs(F_id)%tag_S)//']')
            else
               call msg(MSG_WARNING,'(outcfg) unrecognize token (skipping): '//trim(F_string_S))
            endif
         end select
      endif IF_VAL

      if (.not.RMN_IS_OK(istat)) then
         F_istat = RMN_ERR
         call msg(MSG_WARNING,'(outcfg) invalid token parameter: '//trim(F_string_S))
      else
         call msg(MSG_DEBUG,'(outcfg) '//trim(key_S)//' :=: '//trim(val_S))
      endif
      !----------------------------------------------------------------------
      return
   end function outcfg_parse_string


   !/@*
   function outcfg_getlist(F_id,F_istep,F_varlist_S) result(F_nvar)
      implicit none
      integer,intent(in) :: F_id,F_istep
      character(len=*),intent(inout) :: F_varlist_S(:)
      integer :: F_nvar
      !*@/
      integer :: n,nmax,istat,ivar
      !------------------------------------------------------------------
      F_nvar = RMN_ERR
      if (F_id < 1 .or. F_id > m_ncfgs) return

      F_nvar = 0
      F_varlist_S(:) = ' '
      nmax = size(F_varlist_S)

      do n=1,m_cfgs(F_id)%no
         istat = priv_is_stepset(F_id,m_cfgs(F_id)%o(n)%s_idx,F_istep)
         if (RMN_IS_OK(istat)) then
            do ivar=1,m_cfgs(F_id)%o(n)%n
               if (.not.any(m_cfgs(F_id)%o(n)%on_S(ivar) == F_varlist_S(1:F_nvar))) then
                  F_nvar  = min(F_nvar+1,nmax)
                  F_varlist_S(F_nvar) = m_cfgs(F_id)%o(n)%on_S(ivar)
               endif
            enddo
         endif
      enddo
      !------------------------------------------------------------------
      return
   end function outcfg_getlist


   !/@*
   function outcfg_get_idxlist(F_id,F_istep,F_varname_S,F_idxlist) result(F_nidx)
      implicit none
      integer,intent(in) :: F_id,F_istep
      character(len=*),intent(in) :: F_varname_S
      integer,intent(out) :: F_idxlist(:)
      integer :: F_nidx
      !*@/
      integer :: n,nmax,istat
      !------------------------------------------------------------------
      F_nidx = RMN_ERR
      if (F_id < 1 .or. F_id > m_ncfgs) return

      F_nidx = 0
      F_idxlist = 0
      nmax = size(F_idxlist)

      do n=1,m_cfgs(F_id)%no
         istat = priv_is_stepset(F_id,m_cfgs(F_id)%o(n)%s_idx,F_istep)
         if (RMN_IS_OK(istat)) then
            if (any(m_cfgs(F_id)%o(n)%on_S(:) == F_varname_S)) then
               F_nidx = min(F_nidx+1,nmax)
               F_idxlist(F_nidx) = n
            endif
         endif
      enddo
      !------------------------------------------------------------------
      return
   end function outcfg_get_idxlist


   !/@*
   function outcfg_meta(F_id,F_tag_S,F_units_S,F_closeopen,F_ndigits, &
        F_debug_L,F_rewrite_L,F_fullplane_L,F_closeopen_units) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(out),optional :: F_tag_S,F_units_S
      integer,intent(out),optional :: F_closeopen,F_closeopen_units,F_ndigits
      logical,intent(out),optional :: F_debug_L,F_rewrite_L,F_fullplane_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-09
      !*@/
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_id < 1 .or. F_id > m_ncfgs) return
      F_istat = RMN_OK

      if (present(F_tag_S)) F_tag_S = m_cfgs(F_id)%tag_S(2:)
      if (present(F_units_s)) F_units_S = m_cfgs(F_id)%units_S
      if (present(F_closeopen)) F_closeopen = m_cfgs(F_id)%closeopen
      if (present(F_closeopen_units)) F_closeopen_units = m_cfgs(F_id)%closeopen_units
      if (present(F_ndigits)) F_ndigits = m_cfgs(F_id)%ndigits
      if (present(F_debug_L)) F_debug_L = m_cfgs(F_id)%debug_L
      if (present(F_rewrite_L)) F_rewrite_L = m_cfgs(F_id)%rewrite_L
      if (present(F_fullplane_L)) F_fullplane_L = m_cfgs(F_id)%fullplane_L
      !----------------------------------------------------------------------
      return
   end function outcfg_meta


   !/@*
   function outcfg_var_meta(F_id,F_varidx,F_varname_S,&
        F_datyp,F_nbits,F_filt_pass,F_filt_coef, &
        F_rewrite_L,F_flip_L,F_etk_S,F_ip3,F_fullplane_L) result(F_istat)
      !TODO: F_varidx not used
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_varidx
      character(len=*),intent(in) :: F_varname_S
      integer,intent(out),optional :: F_datyp,F_nbits,F_filt_pass,F_ip3
      real,intent(out),optional :: F_filt_coef
      logical,intent(out),optional :: F_rewrite_L,F_flip_L,F_fullplane_L
      character(len=*),intent(out),optional :: F_etk_S
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2012-01
      !*@/
      character(len=64) :: varname_S
      integer :: ii,nvar,istat
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_id < 1 .or. F_id > m_ncfgs) return
      F_istat = RMN_OK

      varname_S = F_varname_S
      istat = clib_tolower(varname_S)

      if (F_varidx < 1) then
         call msg(MSG_DEBUG,'(outcfg) outcfg_var_meta: negative value for F_varidx -- Ignored/Unused')
      endif

      if (present(F_ip3)) F_ip3 = m_cfgs(F_id)%ip3
      if (present(F_etk_S)) F_etk_S = m_cfgs(F_id)%etk_S
      if (present(F_rewrite_L)) F_rewrite_L = m_cfgs(F_id)%rewrite_L
      if (present(F_fullplane_L)) F_fullplane_L = m_cfgs(F_id)%fullplane_L
      if (present(F_flip_L)) F_flip_L = m_cfgs(F_id)%flip_L
      if (present(F_filt_pass).and.present(F_filt_coef)) then
         F_filt_pass = m_cfgs(F_id)%filt_pass
         F_filt_coef = m_cfgs(F_id)%filt_coef
         do ii=1,m_cfgs(F_id)%nf
            nvar = m_cfgs(F_id)%f(ii)%n
            if (any(varname_S == m_cfgs(F_id)%f(ii)%n_S(1:nvar))) then
               F_filt_pass = m_cfgs(F_id)%f(ii)%npass
               F_filt_coef = m_cfgs(F_id)%f(ii)%coef
            endif
         enddo
      endif
      if (present(F_nbits)) then
         F_nbits = m_cfgs(F_id)%nbits
        do ii=1,m_cfgs(F_id)%nb
            nvar = m_cfgs(F_id)%b(ii)%n
            if (any(varname_S == m_cfgs(F_id)%b(ii)%n_S(1:nvar))) then
               F_nbits = m_cfgs(F_id)%b(ii)%nbits
            endif
         enddo
      endif
      if (present(F_datyp)) then
         F_datyp = m_cfgs(F_id)%datyp
         do ii=1,m_cfgs(F_id)%nb
            nvar = m_cfgs(F_id)%b(ii)%n
            if (any(varname_S == m_cfgs(F_id)%b(ii)%n_S(1:nvar))) then
               F_datyp = m_cfgs(F_id)%b(ii)%datyp
            endif
         enddo
      endif
      !----------------------------------------------------------------------
      return
   end function outcfg_var_meta


   !/@*
   function outcfg_level_meta(F_id,F_varidx,F_lvltype_S,F_v_int_S,F_levels,F_nlevels,F_nlinbot) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_varidx
      character(len=*),intent(out) :: F_lvltype_S,F_v_int_S
      real,intent(out) :: F_levels(:)
      integer,intent(out) :: F_nlevels,F_nlinbot
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2012-01
      !*@/
      integer :: l_idx
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      F_nlevels = 0
      if (F_id < 1 .or. F_id > m_ncfgs) return
      if (F_varidx < 1 .or. F_varidx > m_cfgs(F_id)%no) return

      l_idx = m_cfgs(F_id)%o(F_varidx)%l_idx

      F_lvltype_S = OUTCFG_LVLTYPE(m_cfgs(F_id)%l(l_idx)%type)
      F_v_int_S = m_cfgs(F_id)%l(l_idx)%vint_S
      F_nlinbot = m_cfgs(F_id)%l(l_idx)%linbot

      F_nlevels = m_cfgs(F_id)%l(l_idx)%n
      if (F_nlevels > size(F_levels)) return
      F_levels(1:F_nlevels) = m_cfgs(F_id)%l(l_idx)%v(1:F_nlevels)

      F_istat = F_nlevels
      !----------------------------------------------------------------------
      return
   end function outcfg_level_meta


   !/@*
   function outcfg_grid_meta(F_id,F_varidx,F_i0,F_j0,F_in,F_jn,F_dij,F_tag_S,F_fullplane_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_varidx
      integer,intent(out) :: F_i0,F_j0,F_in,F_jn,F_dij
      character(len=*),intent(out) :: F_tag_S
      logical,intent(out),optional :: F_fullplane_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2012-01
      !*@/
      integer :: g_idx
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_id < 1 .or. F_id > m_ncfgs) return
      if (F_varidx < 1 .or. F_varidx > m_cfgs(F_id)%no) return

      g_idx = m_cfgs(F_id)%o(F_varidx)%g_idx
      F_i0 = m_cfgs(F_id)%g(g_idx)%ij0n(1)
      F_in = m_cfgs(F_id)%g(g_idx)%ij0n(2)
      F_j0 = m_cfgs(F_id)%g(g_idx)%ij0n(3)
      F_jn = m_cfgs(F_id)%g(g_idx)%ij0n(4)
      F_dij = m_cfgs(F_id)%g(g_idx)%ij0n(5)
      F_tag_S = m_cfgs(F_id)%g(g_idx)%tag_S
      if (present(F_fullplane_L)) F_fullplane_L = m_cfgs(F_id)%fullplane_L
      F_istat = RMN_OK
      !----------------------------------------------------------------------
      return
   end function outcfg_grid_meta


   !/@*
   function outcfg_time_4(F_id,F_step,F_datev,F_dateo,F_dt,F_iscloseopen_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      integer,intent(out) :: F_datev,F_dateo,F_dt
      logical,intent(out),optional :: F_iscloseopen_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer(INT64) :: jdateo, jdatev
      !----------------------------------------------------------------------
      F_datev = RMN_ERR
      F_dateo = RMN_ERR
      if (present(F_iscloseopen_L)) then
         F_istat = outcfg_time_8(F_id,F_step,jdatev,jdateo,F_dt,F_iscloseopen_L)
      else
         F_istat = outcfg_time_8(F_id,F_step,jdatev,jdateo,F_dt)
      endif
      if (RMN_IS_OK(F_istat)) then
         F_datev = jdate_to_cmc(jdatev)
         F_dateo = jdate_to_cmc(jdateo)
      endif
      !----------------------------------------------------------------------
      return
   end function outcfg_time_4


   !/@*
   function outcfg_time_8(F_id,F_step,F_datev,F_dateo,F_dt,F_iscloseopen_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      integer(INT64),intent(out) :: F_datev,F_dateo
      integer,intent(out) :: F_dt
      logical,intent(out),optional :: F_iscloseopen_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-04
      !*@/
      integer :: yy1, mo1
      integer(INT64) :: n_sec_8,step_8, dt_8, jdate1stof
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_id < 1 .or. F_id > m_ncfgs) return
      F_istat = RMN_OK

      step_8   = F_step
      dt_8     = m_cfgs(F_id)%dt
      n_sec_8  = step_8 * dt_8
      F_datev  = m_cfgs(F_id)%dateo + n_sec_8
      F_dateo  = m_cfgs(F_id)%dateo
      F_dt     = m_cfgs(F_id)%dt

      if (present(F_iscloseopen_L)) then
         if (m_cfgs(F_id)%closeopen_units == OUTCFG_UNITS_1STOF_MONTH) then
            F_iscloseopen_L = .true.
            if (F_step > 0) then
               yy1 = jdate_year(F_datev)
               mo1 = jdate_month(F_datev)
               call mu_ymdhms2js(jdate1stof,  yy1,  mo1,  1, 0, 0, 0)
               n_sec_8 = jdate1stof - F_datev
               F_iscloseopen_L = .not.( &
                    n_sec_8<0. .or. n_sec_8>=m_cfgs(F_id)%dt)
            endif
         else !OUTCFG_UNITS_STEP
            F_iscloseopen_L = (mod(F_step,m_cfgs(F_id)%closeopen) == 0)
         endif
      endif
      !----------------------------------------------------------------------
      return
   end function outcfg_time_8

   !==== Private Functions =================================================

   !/@*
   function priv_parse_cfgfile(F_id,F_filename_S,F_intags_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_filename_S !full/rel path of config file
      logical,intent(in) :: F_intags_L
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-09
      !*@/
      logical :: insection_L
      integer :: istat,istat2,fileid,n
      character(len=2048) :: string_S,string1_S
      !---------------------------------------------------------------------
      write(string_S, '(a,1x,i0.1x,a)') '(priv_parse_cfgfile) ',F_id,trim(F_filename_S)
      call msg(MSG_DEBUG,string_S)

      F_istat = RMN_ERR

      istat = clib_isfile(trim(F_filename_S))
      istat = min(clib_isreadok(trim(F_filename_S)),istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(outcfg) File not found or not readable: '//trim(F_filename_S))
         return
      endif

      fileid = 0
      istat = fnom(fileid,F_filename_S,'SEQ/FMT+R/O+OLD',0)
      if (.not.RMN_IS_OK(istat) .or. fileid <= 0) then
         call msg(MSG_WARNING,'(outcfg) Problem opening file: '//trim(F_filename_S))
         return
      endif

      insection_L = .not.F_intags_L
      string1_S = ' '
      istat = 0
      DOFILE: do
         read(fileid,'(a)',iostat=istat2) string_S
         if (istat2 /= 0) exit DOFILE
         call str_tab2space(string_S)
         string_S = adjustl(string_S)
         if (string_S(1:1)/='#' .and. string_S/=' ') then
            if (F_intags_L) then
               if (any(string_S == (/TAG_START,TAG_END/))) then
                  insection_L = .true.
                  if (string_S == TAG_END)   insection_L = .false.
                  cycle DOFILE
               endif
            endif
            if (.not.insection_L) print *,'Ignoring: ',trim(string_S)
            if (.not.insection_L) cycle DOFILE
            string1_S = trim(string1_S)//trim(string_S)
            n = len_trim(string_S)
            if (string_S(n:n) == ';' .or. string_S(n:n) == ')') then
               istat = min(outcfg_parse_string(F_id,string1_S),istat)
               string1_S = ' '
            else
               cycle DOFILE
            endif
            if (.not.RMN_IS_OK(istat)) then
               !exit DOFILE
               cycle DOFILE
            endif
         endif
      enddo DOFILE
      if (string1_S /= ' ') then
         istat = RMN_ERR
         call msg(MSG_ERROR,'(outcfg) Token not recognized: '//trim(string1_S))
      endif
!!$      if (RMN_IS_OK(istat)) F_nvar = outcfg_nbvar(F_id)
      F_istat = istat
      istat = fclos(fileid)
      !----------------------------------------------------------------------
      return
   end function priv_parse_cfgfile


   !/@*
   function priv_parse_closeopen(F_out,F_closeopen_units,F_str_S,F_dt) result(F_istat)
      implicit none
      integer,intent(out) :: F_out,F_closeopen_units
      integer,intent(in)  :: F_dt
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      character(len=64) :: str1_S,str2_S
      integer :: itmp
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      F_out = 0

      call str_split(str1_S,str2_S,F_str_S,',')

      if (str1_S(1:5) == 'month') then
         F_closeopen_units = OUTCFG_UNITS_1STOF_MONTH
         F_out = 1 ; F_istat = RMN_OK
         return
      endif

      F_closeopen_units = OUTCFG_UNITS_STEP
      if (str2_S == ' ') then
         F_istat = str_toint(F_out,F_str_S)
         return
      endif

      F_istat = str_toint(F_out,str2_S)
      if (.not.RMN_IS_OK(F_istat)) return

      itmp = priv_time2steps(abs(F_out),str1_S,abs(F_dt))
      if (.not.RMN_IS_OK(itmp)) F_istat = RMN_ERR
      F_out = itmp * sign(1,F_out) * sign(1,F_dt)

      !TODO-later: adjust closeopen to dir change frequency (common denominator)
      !---------------------------------------------------------------------
      return
   end function priv_parse_closeopen


   !/@*
   function priv_time2steps(F_in,F_units_S,F_dt) result(F_out)
      implicit none
      integer,intent(in)  :: F_in,F_dt
      character(len=*),intent(in) :: F_units_S
      integer :: F_out
      !*@/
      !---------------------------------------------------------------------
      F_out = RMN_ERR

      SEL_CLOSE: select case(F_units_S(1:2))
      case('se') !seconds
         F_out = nint(max(1., real(F_in)/real(F_dt)))
      case('mi') !minutes
         F_out = nint(max(1., 60.*real(F_in)/real(F_dt)))
      case('ho') !hours
         F_out = nint(max(1., 3660.*real(F_in)/real(F_dt)))
      case('da') !days
         F_out = nint(max(1., 86400.*real(F_in)/real(F_dt)))
      case('st') !steps
         F_out = F_in
      case('  ') !steps
         F_out = F_in
      case default
         F_out = -1
      end select SEL_CLOSE
     !---------------------------------------------------------------------
      return
   end function priv_time2steps


   !/@*
   function priv_parse_levels(F_id,F_str_S) result(F_istat)
      implicit none
      integer,intent(in)  :: F_id
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer,parameter :: NMIN_ARGS = 2
      integer :: nl,nitems,istat,ii,nrlist,itmp
      real :: rlist(NLEVELS_MAX),loop(3),v,v2
      character(len=16) :: items_S(128)
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
!!$      print *,'levels='//trim(F_str_S)
      F_istat = RMN_ERR
      if (m_cfgs(F_id)%nl >= NSETS_MAX) then
         call msg(MSG_ERROR,'(outcfg) too many levels sets')
         return
      endif
      nl = m_cfgs(F_id)%nl + 1

      nitems = priv_split_args(items_S,F_str_S)
!!$      print *,'levels#',nitems,NMIN_ARGS,':',items_S(1:nitems)
      if (nitems < NMIN_ARGS .or. items_S(1) == ' ') return

      m_cfgs(F_id)%l(nl)%id_S = items_S(1)
      m_cfgs(F_id)%l(nl)%vint_S = m_cfgs(F_id)%vinterp_s
      m_cfgs(F_id)%l(nl)%type = RMN_ERR
      m_cfgs(F_id)%l(nl)%linbot = m_cfgs(F_id)%linbot
      m_cfgs(F_id)%l(nl)%n = 0
      allocate(m_cfgs(F_id)%l(nl)%v(NLEVELS_MAX),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(outcfg) probleme allocating level list mem')
         return
      endif
      m_cfgs(F_id)%l(nl)%v(:) = NO_VAL

      nrlist = 0
      rlist(:) = NO_VAL
      loop(:) = NO_VAL
      ii = 1
      do
         ii = ii + 1
         if (ii > nitems) exit
         select case(items_S(ii)(1:4))
         case('pres') !pres
            if (all(m_cfgs(F_id)%l(nl)%type /= (/RMN_ERR,LVL_TYPE_PRES/))) return
            m_cfgs(F_id)%l(nl)%type = LVL_TYPE_PRES
         case('cubi') !cubic
            m_cfgs(F_id)%l(nl)%vint_S = 'cubic'
         case('line') !linear
            m_cfgs(F_id)%l(nl)%vint_S = 'linear'
         case('linb') !linbot
            if (ii == nitems) return
            ii = ii + 1
            istat = str_toint(itmp,items_S(ii))
            if (.not.RMN_IS_OK(istat) .or. itmp<0) return
            m_cfgs(F_id)%l(nl)%linbot = itmp
         case('-1  ')
            nrlist = 1
            rlist(1) = ALL_VAL
            if (all(m_cfgs(F_id)%l(nl)%type /= (/RMN_ERR,LVL_TYPE_MODEL/))) return
            m_cfgs(F_id)%l(nl)%type = LVL_TYPE_MODEL
         case('all ')
            nrlist = 1
            rlist(1) = ALL_VAL
            if (all(m_cfgs(F_id)%l(nl)%type /= (/RMN_ERR,LVL_TYPE_MODEL/))) return
            m_cfgs(F_id)%l(nl)%type = LVL_TYPE_MODEL
         case('0   ')
            if (nrlist >= NLEVELS_MAX) return
            nrlist = nrlist + 1
            rlist(nrlist) = SURF_VAL
            if (all(m_cfgs(F_id)%l(nl)%type /= (/RMN_ERR,LVL_TYPE_MODEL/))) return
            m_cfgs(F_id)%l(nl)%type = LVL_TYPE_MODEL
         case('surf')
            if (nrlist >= NLEVELS_MAX) return
            nrlist = nrlist + 1
            rlist(nrlist) = SURF_VAL
            if (all(m_cfgs(F_id)%l(nl)%type /= (/RMN_ERR,LVL_TYPE_MODEL/))) return
            m_cfgs(F_id)%l(nl)%type = LVL_TYPE_MODEL
         case('diag')
            if (nrlist >= NLEVELS_MAX) return
            nrlist = nrlist + 1
            rlist(nrlist) = SURF_VAL
            if (all(m_cfgs(F_id)%l(nl)%type /= (/RMN_ERR,LVL_TYPE_MODEL/))) return
            m_cfgs(F_id)%l(nl)%type = LVL_TYPE_MODEL
         case('[   ')
            ii = priv_parse_rlist(items_S,ii,nitems,rlist,nrlist,NLEVELS_MAX)
            if (.not.RMN_IS_OK(ii)) return
         case('<   ')
            ii = priv_parse_loop(items_S,ii,nitems,loop)
            if (.not.RMN_IS_OK(ii)) return
         case default !eta/model or numbers
            if (any(items_S(ii)(1:2) == (/'et','mo'/))) then
               if (m_cfgs(F_id)%l(nl)%type /= RMN_ERR) return
               m_cfgs(F_id)%l(nl)%type = LVL_TYPE_MODEL
            else
               if (m_cfgs(F_id)%l(nl)%type /= RMN_ERR) m_cfgs(F_id)%l(nl)%type  = LVL_TYPE_MODEL
               if (nrlist >= NLEVELS_MAX) return
               nrlist = nrlist + 1
               istat = str_toreal(rlist(nrlist),items_S(ii))
               if (.not.RMN_IS_OK(istat)) return
            endif
         end select
      enddo

      if (all(loop >= NO_VAL) .and. loop(3) > 0.) then
         v2 = maxval(loop(1:2))
         v = minval(loop(1:2))
         do while (v <= v2)
            if (nrlist+1 > NLEVELS_MAX) exit
            nrlist = nrlist + 1
            rlist(nrlist) = v
            v = v + loop(3)
         enddo
      endif
      if (m_cfgs(F_id)%l(nl)%type == LVL_TYPE_PRES) then
         nrlist = sort(rlist(1:nrlist),SORT_DOWN,SORT_UNIQUE)
      else
         nrlist = sort(rlist(1:nrlist),SORT_UP,SORT_UNIQUE)
      endif

      if (m_cfgs(F_id)%l(nl)%type == LVL_TYPE_PRES) then
         if (nrlist > 0 .and. any(rlist(1:nrlist)<=0.)) then
            call msg(MSG_WARNING,'(outcfg) invalid list of pres levels - ignored')
            return
         endif
      endif

      m_cfgs(F_id)%nl = nl
      m_cfgs(F_id)%l(nl)%n = nrlist
      m_cfgs(F_id)%l(nl)%v = rlist
      F_istat = RMN_OK

      msg_S = 'model'
      if (m_cfgs(F_id)%l(nl)%type == LVL_TYPE_PRES) then
         if (m_cfgs(F_id)%l(nl)%vint_S(1:1) == 'c') then
            write(msg_S,'(a,i4,a)') 'pres ('//m_cfgs(F_id)%l(nl)%vint_S(1:6)//' linbot=',m_cfgs(F_id)%l(nl)%linbot,')'
         else
            write(msg_S,'(a,i4,a)') 'pres ('//m_cfgs(F_id)%l(nl)%vint_S(1:6)//')'
         endif
      endif
      msg_S = '(outcfg) levelset: '//m_cfgs(F_id)%l(nl)%id_S(1:4)//' '//trim(msg_S)
      if (nrlist > 0) then
         if (nrlist == 1) then
            write(msg_S,'(a,f5.0)') trim(msg_S),rlist(1)
         else if (nrlist == 2) then
            write(msg_S,'(a,2f5.0)') trim(msg_S),rlist(1),rlist(2)
         else
            write(msg_S,'(a,f5.0,a,f5.0)') trim(msg_S),rlist(1),', ...,',rlist(nrlist)
         endif
      endif
      call msg(MSG_INFO,msg_S)
      !---------------------------------------------------------------------
      return
   end function priv_parse_levels


   !/@*
   function priv_parse_grid(F_id,F_str_S) result(F_istat)
      implicit none
      integer,intent(in)  :: F_id
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer,parameter :: NMIN_ARGS = 2
      integer :: ng,nitems,istat,ii,reduc
      character(len=16) :: items_S(128),type_S
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
!!$      print *,'grid='//trim(F_str_S)
      F_istat = RMN_ERR
      if (m_cfgs(F_id)%ng >= NSETS_MAX) then
         call msg(MSG_ERROR,'(outcfg) too many grid sets')
         return
      endif
      ng = m_cfgs(F_id)%ng + 1

      nitems = priv_split_args(items_S,F_str_S)
!!$      print *,'grid#',nitems,NMIN_ARGS,':',items_S(1:nitems)
      if (nitems < NMIN_ARGS .or. items_S(1) == ' ') return

      m_cfgs(F_id)%g(ng)%id_S = items_S(1)
      m_cfgs(F_id)%g(ng)%tag_S = ' '
      m_cfgs(F_id)%g(ng)%ij0n = (/0,0,0,0,1/)

      reduc = 0
      ii = 1
      do
         ii = ii + 1
         if (ii > nitems) exit
         select case(items_S(ii)(1:4))
         case('mode') !model
            type_S = items_S(ii)
            m_cfgs(F_id)%g(ng)%ij0n = (/1,-1,1,-1,1/)
            reduc = size(m_cfgs(F_id)%g(ng)%ij0n)+1
         case('core')
            type_S = items_S(ii)
            m_cfgs(F_id)%g(ng)%ij0n = m_cfgs(F_id)%reduc_core
            reduc = size(m_cfgs(F_id)%g(ng)%ij0n)+1
!!$         case('free')
!!$            type_S = items_S(ii)
!!$            m_cfgs(F_id)%g(ng)%ij0n = (/3,-3,3,-3,1/)
!!$            reduc = size(m_cfgs(F_id)%g(ng)%ij0n)+1
         case('redu') !reduc
            type_S = items_S(ii)
            reduc = 1
         case default
            if (any(items_S(ii)(1:1) == (/'"',"'"/))) then
               m_cfgs(F_id)%g(ng)%tag_S = items_S(ii)
               call str_rm_quotes(m_cfgs(F_id)%g(ng)%tag_S)
            else if (reduc > 0 .and. reduc <= size(m_cfgs(F_id)%g(ng)%ij0n)) then
               istat = str_toint(m_cfgs(F_id)%g(ng)%ij0n(reduc),items_S(ii))
               if (.not.RMN_IS_OK(istat)) return
               reduc = reduc + 1
            else
               return
            endif
        end select
      enddo
      m_cfgs(F_id)%ng = ng
      F_istat = RMN_OK

      write(msg_S,'(a,5i6)') '(outcfg) grid set: '//m_cfgs(F_id)%g(ng)%id_S(1:4)//' '//type_S(1:5)//' : ',m_cfgs(F_id)%g(ng)%ij0n
      call msg(MSG_INFO,msg_S)
     !---------------------------------------------------------------------
      return
   end function priv_parse_grid


   !/@*
   function priv_parse_steps(F_id,F_str_S) result(F_istat)
      implicit none
      integer,intent(in)  :: F_id
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer,parameter :: NMIN_ARGS = 2
      integer :: ns,nitems,ii,istat,nrlist
      real :: rlist(NLIST_MAX),loop(3)
      character(len=16) :: items_S(128),timetype_S
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
!!$      print *,'steps='//trim(F_str_S)
      F_istat = RMN_ERR
      if (m_cfgs(F_id)%ns >= NSETS_MAX) then
         call msg(MSG_ERROR,'(outcfg) too many steps sets')
         return
      endif
      ns = m_cfgs(F_id)%ns + 1

      nitems = priv_split_args(items_S,F_str_S)
!!$      print *,'steps#',nitems,NMIN_ARGS,':',items_S(1:nitems)
      if (nitems < NMIN_ARGS .or. items_S(1) == ' ') return

      m_cfgs(F_id)%s(ns)%id_S = items_S(1)
      m_cfgs(F_id)%s(ns)%init_L = .false.
      m_cfgs(F_id)%s(ns)%n = 0
      allocate(m_cfgs(F_id)%s(ns)%v(NLIST_MAX),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(outcfg) probleme allocating steps list mem')
         return
      endif
      m_cfgs(F_id)%s(ns)%v(:) = int(NO_VAL)
      m_cfgs(F_id)%s(ns)%l(:) = int(NO_VAL)

      nrlist = 0
      rlist(:) = NO_VAL
      loop(:) = NO_VAL
      timetype_S = 'step'
      ii = 1
      do
         ii = ii + 1
         if (ii > nitems) exit
         select case(items_S(ii)(1:4))
         case('init')
            m_cfgs(F_id)%s(ns)%init_L = .true.
         case('hour')
            timetype_S = 'hour'
         case('step')
            timetype_S = 'step'
         case('[   ')
            ii = priv_parse_rlist(items_S,ii,nitems,rlist,nrlist,NLIST_MAX)
            if (.not.RMN_IS_OK(ii)) return
         case('<   ')
            ii = priv_parse_loop(items_S,ii,nitems,loop)
            if (.not.RMN_IS_OK(ii)) return
         case('-1  ')
            nrlist = 1
            rlist(1) = ALL_VAL
         case('all ')
            nrlist = 1
            rlist(1) = ALL_VAL
         case default !
            if (nrlist >= NLIST_MAX) return
            nrlist = nrlist + 1
            istat = str_toreal(rlist(nrlist),items_S(ii))
            if (.not.RMN_IS_OK(istat)) return
        end select
      enddo
      if (timetype_S == 'hour') then
!!$         where(rlist /= NO_VAL) rlist = rlist * 3600. / m_cfgs(F_id)%dt
!!$         where(loop /= NO_VAL) loop = loop * 3600. / m_cfgs(F_id)%dt
         where(rlist > 0.) rlist = rlist * 3600.
         where(loop > 0.) loop = loop * 3600.
      else
         where(rlist > 0.) rlist = rlist * real(m_cfgs(F_id)%dt)
         where(loop > 0.) loop = loop * real(m_cfgs(F_id)%dt)
      endif
      m_cfgs(F_id)%s(ns)%n = nrlist
      m_cfgs(F_id)%s(ns)%v = int(rlist)
      m_cfgs(F_id)%s(ns)%l = int(loop)
      m_cfgs(F_id)%ns = ns
      F_istat = RMN_OK

      where(rlist > 0.) rlist = rlist / real(m_cfgs(F_id)%dt)
      where(loop > 0.) loop = loop / real(m_cfgs(F_id)%dt)
      write(msg_S,'(a,5i6)') '(outcfg) stepsset: '//m_cfgs(F_id)%s(ns)%id_S(1:4)//' steps : '
      if (nrlist > 0) then
         if (nrlist == 1) then
            write(msg_S,'(a,i5)') trim(msg_S),int(rlist(1))
         else if (nrlist == 2) then
            write(msg_S,'(a,2i5)') trim(msg_S),int(rlist(1)),int(rlist(2))
         else
            write(msg_S,'(a,i5,a,i5)') trim(msg_S),int(rlist(1)),', ...,',int(rlist(nrlist))
         endif
      else
         write(msg_S,'(a,3i5,a)') trim(msg_S)//' loop(',int(loop),')'
      endif
      call msg(MSG_INFO,msg_S)
      !---------------------------------------------------------------------
      return
   end function priv_parse_steps


   !/@*
   function priv_parse_sortie(F_id,F_str_S) result(F_istat)
      implicit none
      integer,intent(in)  :: F_id
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer,parameter :: NVAR_MAX = 256
      integer,parameter :: NMIN_ARGS = 7
      integer :: no,nitems,ii,nvar,idx,istat
      character(len=16) :: items_S(128)
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
!!$      print *,F_id,'sortie='//trim(F_str_S) ; call flush(6)
      F_istat = RMN_ERR
      if (m_cfgs(F_id)%no >= NSETS_MAX) then
         call msg(MSG_ERROR,'(outcfg) too many sortie'//trim(m_cfgs(F_id)%tag_S)//' commands')
         return
      endif
      no = m_cfgs(F_id)%no + 1

      nitems = priv_split_args(items_S,F_str_S)
!!$      print *,'sortie#',nitems,NMIN_ARGS,':',items_S(1:nitems) ; call flush(6)
      if (nitems < NMIN_ARGS) return

      m_cfgs(F_id)%o(no)%s_idx = RMN_ERR
      m_cfgs(F_id)%o(no)%g_idx = RMN_ERR
      m_cfgs(F_id)%o(no)%l_idx = RMN_ERR
      m_cfgs(F_id)%o(no)%n = 0
      allocate(m_cfgs(F_id)%o(no)%on_S(NLIST_MAX),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(outcfg) probleme allocating sortie name list mem')
         return
      endif
      m_cfgs(F_id)%o(no)%on_S(:) = ' '

      nvar = 0
      ii = 0
      do
         ii = ii + 1
         if (ii > nitems) exit
         select case(items_S(ii)(1:4))
         case('grid')
            if (ii == nitems) return
            ii = ii + 1
            do idx=1,m_cfgs(F_id)%ng
               if (m_cfgs(F_id)%g(idx)%id_S == items_S(ii)) then
                  m_cfgs(F_id)%o(no)%g_idx = idx
                  exit
               endif
            enddo
            if (.not.RMN_IS_OK(m_cfgs(F_id)%o(no)%g_idx)) return
         case('step')
            if (ii == nitems) return
            ii = ii + 1
            do idx=1,m_cfgs(F_id)%ns
               if (m_cfgs(F_id)%s(idx)%id_S == items_S(ii)) then
                  m_cfgs(F_id)%o(no)%s_idx = idx
                  exit
               endif
            enddo
            if (.not.RMN_IS_OK(m_cfgs(F_id)%o(no)%s_idx)) return
         case('leve')
            if (ii == nitems) return
            ii = ii + 1
            do idx=1,m_cfgs(F_id)%nl
               if (m_cfgs(F_id)%l(idx)%id_S == items_S(ii)) then
                  m_cfgs(F_id)%o(no)%l_idx = idx
                  exit
               endif
            enddo
            if (.not.RMN_IS_OK(m_cfgs(F_id)%o(no)%l_idx)) return
         case default !
            if (.not.any(items_S(ii) == (/'[',']'/))) then
               if (nvar == NVAR_MAX) cycle
               nvar = nvar + 1
               m_cfgs(F_id)%o(no)%on_S(nvar) = items_S(ii)
               call str_rm_quotes(m_cfgs(F_id)%o(no)%on_S(nvar))
            endif
         end select
      enddo
      if (.not.(RMN_IS_OK(m_cfgs(F_id)%o(no)%s_idx) &
           .and.RMN_IS_OK(m_cfgs(F_id)%o(no)%g_idx) &
           .and.RMN_IS_OK(m_cfgs(F_id)%o(no)%l_idx))) then
!!$         call msg(MSG_ERROR,'(outcfg) Invalid sortie directive: '//trim(F_str_S))
         return
      endif
      m_cfgs(F_id)%no = no
      m_cfgs(F_id)%o(no)%n = nvar
      F_istat = RMN_OK

      msg_S = '(outcfg) sortie'//trim(m_cfgs(F_id)%tag_S)//' : (step,grid,levl)=('// &
           m_cfgs(F_id)%s(m_cfgs(F_id)%o(no)%s_idx)%id_S(1:4)// &
           m_cfgs(F_id)%g(m_cfgs(F_id)%o(no)%g_idx)%id_S(1:4)// &
           m_cfgs(F_id)%l(m_cfgs(F_id)%o(no)%l_idx)%id_S(1:4)//')'
      if (nvar == 1) then
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%o(no)%on_S(1)(1:4)
      else if (nvar == 2) then
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%o(no)%on_S(1)(1:4)//', '//m_cfgs(F_id)%o(no)%on_S(2)(1:4)
      else
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%o(no)%on_S(1)(1:4)//', ..., '//m_cfgs(F_id)%o(no)%on_S(nvar)(1:4)
      endif
      call msg(MSG_INFO,msg_S)
      !---------------------------------------------------------------------
      return
   end function priv_parse_sortie


   !/@*
   function priv_parse_xnbits(F_id,F_str_S) result(F_istat)
      implicit none
      integer,intent(in)  :: F_id
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer,parameter :: NVAR_MAX = 256
      integer,parameter :: NMIN_ARGS = 3
      integer :: nb,nitems,ii,nvar,istat
      character(len=16) :: items_S(128)
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (m_cfgs(F_id)%nb >= NSETS_MAX) then
         call msg(MSG_ERROR,'(outcfg) too many xnbits '//trim(m_cfgs(F_id)%tag_S)//' commands')
         return
      endif
      nb = m_cfgs(F_id)%nb + 1

      nitems = priv_split_args(items_S,F_str_S)
      if (nitems < NMIN_ARGS) return

      m_cfgs(F_id)%b(nb)%nbits = m_cfgs(F_id)%nbits
      m_cfgs(F_id)%b(nb)%datyp = m_cfgs(F_id)%datyp

      allocate(m_cfgs(F_id)%b(nb)%n_S(NLIST_MAX),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(outcfg) probleme allocating xnbits name list mem')
         return
      endif
      m_cfgs(F_id)%b(nb)%n_S(:) = ' '

      nvar = 0
      ii = 0
      do
         ii = ii + 1
         if (ii > nitems) exit
         select case(items_S(ii)(1:4))
         case('bits')
            if (ii == nitems) return
            ii = ii + 1
            istat = str_toint(m_cfgs(F_id)%b(nb)%nbits,items_S(ii))
            if (.not.RMN_IS_OK(istat)) return
            m_cfgs(F_id)%b(nb)%nbits = max(NBITS_MIN,min(abs(m_cfgs(F_id)%b(nb)%nbits),NBITS_MAX))
         case('daty')
            if (ii == nitems) return
            ii = ii + 1
            istat = str_toint(m_cfgs(F_id)%b(nb)%datyp,items_S(ii))
            if (.not.RMN_IS_OK(istat)) return
            !TODO: check if datyp is part of allowed dtypes
         case default
            if (.not.any(items_S(ii) == (/'[',']'/))) then
               if (nvar == NVAR_MAX) cycle
               nvar = nvar + 1
               m_cfgs(F_id)%b(nb)%n_S(nvar) = items_S(ii)
               call str_rm_quotes(m_cfgs(F_id)%b(nb)%n_S(nvar))
            endif
         end select
      enddo
      m_cfgs(F_id)%nb = nb
      m_cfgs(F_id)%b(nb)%n = nvar
      F_istat = RMN_OK

      write(msg_S,'(a,i3,a)') '(outcfg) xnbits (',m_cfgs(F_id)%b(nb)%nbits,') : '
      if (nvar == 1) then
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%b(nb)%n_S(1)(1:4)
      else if (nvar == 2) then
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%b(nb)%n_S(1)(1:4)//', '//m_cfgs(F_id)%b(nb)%n_S(2)(1:4)
      else
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%b(nb)%n_S(1)(1:4)//', ..., '//m_cfgs(F_id)%b(nb)%n_S(nvar)(1:4)
      endif
      call msg(MSG_INFO,msg_S)
      !---------------------------------------------------------------------
      return
   end function priv_parse_xnbits


   !/@*
   function priv_parse_filtre(F_id,F_str_S) result(F_istat)
      implicit none
      integer,intent(in)  :: F_id
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer,parameter :: NVAR_MAX = 256
      integer,parameter :: NMIN_ARGS = 5
      integer :: nf,nitems,ii,nvar,istat
      character(len=16) :: items_S(128)
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (m_cfgs(F_id)%nf >= NSETS_MAX) then
         call msg(MSG_ERROR,'(outcfg) too many filtre '//trim(m_cfgs(F_id)%tag_S)//' commands')
         return
      endif
      nf = m_cfgs(F_id)%nf + 1

      nitems = priv_split_args(items_S,F_str_S)
      if (nitems < NMIN_ARGS) return

      m_cfgs(F_id)%f(nf)%npass = m_cfgs(F_id)%filt_pass
      m_cfgs(F_id)%f(nf)%coef = m_cfgs(F_id)%filt_coef

      allocate(m_cfgs(F_id)%f(nf)%n_S(NLIST_MAX),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(outcfg) probleme allocating filtre name list mem')
         return
      endif
      m_cfgs(F_id)%f(nf)%n_S(:) = ' '

      nvar = 0
      ii = 0
      do
         ii = ii + 1
         if (ii > nitems) exit
         select case(items_S(ii)(1:4))
         case('pass')
            if (ii == nitems) return
            ii = ii + 1
            istat = str_toint(m_cfgs(F_id)%f(nf)%npass,items_S(ii))
            if (.not.RMN_IS_OK(istat)) return
            m_cfgs(F_id)%f(nf)%npass = max(0,min(abs(m_cfgs(F_id)%f(nf)%npass),1024))
         case('coef')
            if (ii == nitems) return
            ii = ii + 1
            istat = str_toreal(m_cfgs(F_id)%f(nf)%coef,items_S(ii))
            if (.not.RMN_IS_OK(istat)) return
         case default
            if (.not.any(items_S(ii) == (/'[',']'/))) then
               if (nvar == NVAR_MAX) cycle
               nvar = nvar + 1
               m_cfgs(F_id)%f(nf)%n_S(nvar) = items_S(ii)
               call str_rm_quotes(m_cfgs(F_id)%f(nf)%n_S(nvar))
            endif
         end select
      enddo
      m_cfgs(F_id)%nf = nf
      m_cfgs(F_id)%f(nf)%n = nvar
      F_istat = RMN_OK

      write(msg_S,'(a,i4,f6.3,a)') '(outcfg) filtre (',m_cfgs(F_id)%f(nf)%npass,m_cfgs(F_id)%f(nf)%coef,') : '
      if (nvar == 1) then
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%f(nf)%n_S(1)(1:4)
      else if (nvar == 2) then
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%f(nf)%n_S(1)(1:4)//', '//m_cfgs(F_id)%f(nf)%n_S(2)(1:4)
      else
         msg_S = trim(msg_S)//' '//m_cfgs(F_id)%f(nf)%n_S(1)(1:4)//', ..., '//m_cfgs(F_id)%f(nf)%n_S(nvar)(1:4)
      endif
      call msg(MSG_INFO,msg_S)
      !---------------------------------------------------------------------
      return
   end function priv_parse_filtre


   !/@*
   function priv_str2kv(F_string_S,F_key_S,F_val_S) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_string_S
      character(len=*),intent(out) :: F_key_S,F_val_S
      integer :: F_istat
      !*@/
      integer :: istat,n
      character(len=2048) :: string_S
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      F_key_S = ' '
      F_val_S = ' '
      string_S = F_string_S
      call str_tab2space(string_S)
      string_S = adjustl(string_S)
      IF_TOKEN: if (string_S(1:1)/='#' .and. string_S/=' ') then
         F_istat = RMN_ERR
         istat = clib_tolower(string_S)
         call str_split(F_key_S,F_val_S,string_S,'=')
         if (F_val_S == ' ') call str_split(F_key_S,F_val_S,string_S,'(')
         IF_VAL: if (F_val_S /= ' ') then
            F_istat = RMN_OK
            F_key_S = adjustl(F_key_S)
            F_val_S = adjustl(F_val_S)
            n = len_trim(F_val_S)
            if (F_val_S(n:n) == ';' .or. F_val_S(n:n) == ')') F_val_S(n:n) = ' '
         endif IF_VAL
      endif IF_TOKEN
      !---------------------------------------------------------------------
      return
   end function priv_str2kv


   !/@*
   function priv_str2str(F_out_S,F_str_S,F_list_S,F_minlen) result(F_istat)
      implicit none
      character(len=*),intent(out) :: F_out_S
      character(len=*),intent(in) :: F_str_S,F_list_S(:)
      integer,intent(in) :: F_minlen
      integer :: F_istat
      !*@/
      integer :: n,n2,i
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      F_out_S = ' '
      n2 = max(1,min(len_trim(F_str_S),F_minlen))
      do i=1,size(F_list_S)
         n = max(1,min(len_trim(F_list_S(i)),F_minlen))
         if (F_str_S(1:n2) == F_list_S(i)(1:n)) then
            F_istat = RMN_OK
            F_out_S = F_list_S(i)
            exit
        endif
      enddo
      !---------------------------------------------------------------------
      return
   end function priv_str2str


   !/@*
   function priv_split_args(F_items_S,F_str_S) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_str_S
      character(len=*),intent(out) :: F_items_S(:)
      integer :: F_istat
      !*@/
      integer :: ii,nmax,n
      character(len=1024) :: s1_S,s2_S,s0_S
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      F_items_S(:) = ' '
      nmax = size(F_items_S)
      s0_S = adjustl(F_str_S)
      ii = 0
      do
         call str_split(s1_S,s2_S,s0_S,',')
         ii = ii + 1
         if (ii >= nmax) then
            ii = RMN_ERR
            exit
         endif
         s1_S = adjustl(s1_S)
         n = len_trim(s1_S)
         if (any(s1_S(1:1) == (/'[','<'/))) then
            F_items_S(ii) = s1_S(1:1)
            ii = ii + 1
            s1_S = adjustl(s1_S(2:n))
            n = len_trim(s1_S)
         endif
         if (any(s1_S(n:n) == (/']','>'/))) then
            F_items_S(ii) = s1_S(1:max(1,n-1))
            ii = ii + 1
            s1_S = s1_S(n:n)
            n = len_trim(s1_S)
         endif
         F_items_S(ii) = s1_S
         if (s2_S == ' ') exit
         s0_S = adjustl(s2_S)
      enddo
      F_istat = ii
      !---------------------------------------------------------------------
      return
   end function priv_split_args


   !/@*
   function priv_parse_rlist(my_items_S,my_ii,my_nitems,my_rlist,my_nlist,my_nmax) result(my_iout)
      implicit none
      integer,intent(in) :: my_ii,my_nitems,my_nmax
      character(len=*),intent(in) :: my_items_S(my_nitems)
      real,intent(inout) :: my_rlist(:)
      integer,intent(inout) :: my_nlist
      integer :: my_iout
      !*@/
      integer :: ii,istat
      !---------------------------------------------------------------------
      my_iout = my_ii
      if (my_items_S(my_ii) /= '[') return
      my_iout = RMN_ERR
      if (my_ii+2 > my_nitems) return !needs [ and ]
      ii = my_ii
      do
         ii = ii + 1
         if (ii > my_nitems) return  !needs ]
         if (my_items_S(ii) == ']') exit
         if (my_nlist >= my_nmax) return
         istat = str_toreal(my_rlist(my_nlist+1),my_items_S(ii))
         if (.not.RMN_IS_OK(istat)) return
         my_nlist = my_nlist + 1
      enddo
      my_iout = ii
      !---------------------------------------------------------------------
      return
   end function priv_parse_rlist


   !/@*
   function priv_parse_loop(my_items_S,my_ii,my_nitems,my_loop) result(my_iout)
      implicit none
      integer,intent(in) :: my_ii,my_nitems
      character(len=*),intent(in) :: my_items_S(my_nitems)
      real,intent(inout) :: my_loop(3)
      integer :: my_iout
      !*@/
      integer :: istat
      !---------------------------------------------------------------------
      my_iout = my_ii
      if (my_items_S(my_ii) /= '<') return
      my_iout = RMN_ERR
      if (my_ii+4 > my_nitems) return
      if (my_items_S(my_ii+4) /= '>') return

      istat = str_toreal(my_loop(1),my_items_S(my_ii+1))
      istat = min(str_toreal(my_loop(2),my_items_S(my_ii+2)),istat)
      istat = min(str_toreal(my_loop(3),my_items_S(my_ii+3)),istat)
      if (.not.RMN_IS_OK(istat)) return
      my_iout = my_ii + 5
      !---------------------------------------------------------------------
      return
   end function priv_parse_loop


   !/@*
   function priv_is_stepset(F_id,F_s_idx,F_istep) result(F_istat)
      implicit none
      integer,intent(in) :: F_id,F_s_idx,F_istep
      integer :: F_istat
      !*@/
      integer :: now,n,delta
      !------------------------------------------------------------------
      F_istat = RMN_ERR
      now = F_istep * m_cfgs(F_id)%dt
      n = m_cfgs(F_id)%s(F_s_idx)%n
      if ((n>0 .and. m_cfgs(F_id)%s(F_s_idx)%v(1) == int(ALL_VAL)) .or. &
           any(m_cfgs(F_id)%s(F_s_idx)%v(1:n) == now)) then
         F_istat = RMN_OK
      else if (m_cfgs(F_id)%s(F_s_idx)%l(1) <= now .and. &
           now <= m_cfgs(F_id)%s(F_s_idx)%l(2)) then
         delta = now - m_cfgs(F_id)%s(F_s_idx)%l(1)
         if (mod(delta,m_cfgs(F_id)%s(F_s_idx)%l(3)) == 0) F_istat = RMN_OK
!!$            print *,now,m_cfgs(F_id)%s(F_s_idx)%l(1),delta,F_istat
      endif
      !------------------------------------------------------------------
      return
   end function priv_is_stepset

end module outcfg_mod
