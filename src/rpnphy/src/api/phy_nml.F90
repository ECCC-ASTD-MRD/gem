
module phy_nml_mod
   use clib_itf_mod, only: clib_isreadok, clib_toupper
   use debug_mod, only: init2nan_L
   use str_mod, only: str_concat, str_toreal
   use series_mod, only: series_nml
   use phy_status, only: PHY_ERROR, PHY_NONE, PHY_OK, PHY_CTRL_NML_OK, phy_init_ctrl
   use phy_options
   use cnv_options
#ifdef HAVE_NEMO
      use cpl_itf, only: cpl_nml
#endif
   use mixing_length, only: ML_CLOSURES
   private
   public :: phy_nml

!!!#include <arch_specific.hf>
#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   include "rpnphy_version.inc"

contains

   !/@*
   function phy_nml(F_namelist) result(F_status)
      implicit none
      !@object Set defaults values and read physics namelists to initialize physics configuration
      !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014
      !@revision

      !@Arguments
      ! F_namelist    File name containing the namelists to read
      character(len=*), intent(in) :: F_namelist

      !@return
      integer :: F_status
      !*@/

#ifdef HAVE_MACH
      integer, external :: chm_nml
#endif

      integer, external :: sfc_nml2, cnv_nml2, check_options2, msg_getUnit

      integer :: err, unout
      character(len=1024) :: msg_S
      !-------------------------------------------------------------------
      F_status      = PHY_ERROR
      phy_init_ctrl = PHY_ERROR
      unout         = msg_getUnit(MSG_INFO)

      !# Print Banner
      call msg(MSG_INFO, '   ************************************************************')
      write(msg_S, "(3x,'Package: ',a,5x,'version: ',a)") &
           trim(RPNPHY_NAME_S), trim(RPNPHY_VERSION_S)
      call msg(MSG_INFO, msg_S)
      write(msg_S, "(3x,'Release of: ',a,5x,'COMPILER: ',a)") &
           trim(RPNPHY_DSTP_S), trim(RPNPHY_EC_ARCH_S)
      call msg(MSG_INFO, msg_S)
      call msg(MSG_INFO, '   ************************************************************')

      !# Reading namelists physics_cfgs from file F_namelist
      !# Return code F_status will be set to 0 or 1 if no error occurs.

      err = phy_options_init()
      if (.not.RMN_IS_OK(err)) return
      err = phy_nml_read(F_namelist)
      if (.not.RMN_IS_OK(err)) return
      if (err == RMN_OK) then
         F_status = PHY_NONE
         return
      endif
      err = phy_nml_check()
      if (.not.RMN_IS_OK(err)) return

      !# Read surface namelist
      err = sfc_nml2(F_namelist)
      if (.not.RMN_IS_OK(err)) return

#ifdef HAVE_NEMO
      !# Read coupling namelist and initialize coupling configuration
      err   = cpl_nml(F_namelist, unout)
      if (.not.RMN_IS_OK(err)) then
         call msg(MSG_ERROR, '(phy_nml) Problem reading COUPLING namelist')
         return
      endif
#endif

      !# Read convection namelist
      err = cnv_nml2(F_namelist)
      if (.not.RMN_IS_OK(err)) return

      !# Read series namelists
      err = series_nml(F_namelist)
      if (.not.RMN_IS_OK(err)) return

#ifdef HAVE_MACH
      !# Read chemistry namelists and initialize chemistry configuration
      err = chm_nml(F_namelist, unout)
      if (err <= 0) then
         call msg(MSG_ERROR, '(phy_nml) Problem reading chemestry namelist')
         return
      endif
#endif

      !# 
      err = phy_nml_post_init()
      if (.not.RMN_IS_OK(err)) return
      call phy_nml_print()

      err = check_options2()
      if (.not.RMN_IS_OK(err)) then
         call msg(MSG_ERROR, '(phy_nml) invalid options values in check_options')
         return
      endif

      phy_init_ctrl = PHY_CTRL_NML_OK
      F_status = PHY_OK
      !-------------------------------------------------------------------
      return
   end function phy_nml


   function phy_nml_read(m_namelist) result(m_istat)
      implicit none
      character(len=*), intent(in) :: m_namelist
      integer :: m_istat
      character(len=1024) :: msg_S, namelist_S, name_S
      integer :: istat, unf
      !----------------------------------------------------------------
      m_istat = RMN_ERR
      name_S = '(phy_nml)'
      namelist_S = '&physics_cfgs'

      unf = 0
      istat = clib_isreadok(m_namelist)
      if (RMN_IS_OK(istat)) istat = fnom(unf, m_namelist, 'SEQ+OLD', 0)
      if (.not.RMN_IS_OK(istat)) then
!!$         write(msg_S,'(a,a,a,a)') trim(name_S), &
!!$              ' Using default config, ', &
!!$              'Namelist file Not found/readable: ', trim(m_namelist)
         write(msg_S,'(a,a,a,a)') trim(name_S), &
              ' Running without Physics, ', &
              'Namelist file Not found/readable: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
         return
      endif

      read(unf, nml=physics_cfgs, iostat=istat)
      if (istat == 0) then     !# Read ok
         m_istat = RMN_OK + 1
      else if (istat < 0) then !# EOF, nml not found
!!$         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
!!$              ' Using default config, Namelist ',&
!!$              trim(namelist_S), ' not available in file: ', trim(m_namelist)
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Running without Physics, Namelist ',&
              trim(namelist_S), ' not available in file: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
      else !# Error
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Namelist', trim(namelist_S), &
              ' invalid in file: ', trim(m_namelist)
         call msg(MSG_ERROR, msg_S)
      endif
      istat = fclos(unf)
      !----------------------------------------------------------------
      return
   end function phy_nml_read


   subroutine phy_nml_print()
      implicit none
      integer, external :: msg_getUnit
      integer :: unout, nrec
      !----------------------------------------------------------------
      unout = msg_getUnit(MSG_INFO)
      if (unout > 0) then
         write(unout, nml=physics_cfgs_p)
         do nrec = 1, size(indiag_list_s)
            if (indiag_list_s(nrec) == ' ') exit
         enddo
         write(unout, *) 'indiag_list_s = ',indiag_list_s(1:max(1,nrec-1))
      endif
      !----------------------------------------------------------------
      return
   end subroutine phy_nml_print


   function phy_nml_check() result(m_istat)
      implicit none
      integer :: m_istat
            
      integer :: istat
      character(len=1024) :: msg_S
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      if (priv_stropt(cond_conserve,   'cond_conserve',   COND_CONSERVE_OPT) == RMN_ERR) return
      if (priv_stropt(cond_sgspdf,     'cond_sgspdf',     COND_SGSPDF_OPT) == RMN_ERR) return
      if (priv_stropt(fluvert,         'fluvert',         FLUVERT_OPT) == RMN_ERR) return
      if (priv_stropt(gwdrag,          'gwdrag',          GWDRAG_OPT) == RMN_ERR) return
      
      !#TODO: indiag_list_s?

      if (priv_stropt(input_type,      'input_type',      INPUT_TYPE_OPT) == RMN_ERR) return
      istat = clib_toupper(kntrad_S)
      istat = clib_toupper(kntraduv_S)
      if (priv_stropt(longmel,         'longmel',         ML_CLOSURES(:)%name) == RMN_ERR) return
      if (priv_stropt(pbl_diss,        'pbl_diss',        PBL_DISS_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_dissheat,    'pbl_dissheat',    PBL_DISSHEAT_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_conserve,    'pbl_conserve',    PBL_CONSERVE_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_func_stab,   'pbl_func_stab',   PBL_FUNC_STAB_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_func_unstab, 'pbl_func_unstab', PBL_FUNC_UNSTAB_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_mlblac_max,  'pbl_mlblac_max',  PBL_MLBLAC_MAX_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_nonloc,      'pbl_nonloc',      PBL_NONLOC_OPT) == RMN_ERR) return
      if (priv_stropt(pbl_shal,        'pbl_shal',        PBL_SHAL_OPT) == RMN_ERR) return
      if (priv_stropt(pcptype,         'pcptype',         PCPTYPE_OPT) == RMN_ERR) return

      istat = clib_toupper(phystat_freq_S)
      !#TODO: phystat_list_s?

      if (priv_stropt(rad_atmpath,     'rad_atmpath',     RAD_ATMPATH_OPT) == RMN_ERR) return
      if (priv_stropt(rad_part_nomp,   'rad_part_nomp',   RAD_PART_NOMP_OPT) == RMN_ERR) return
      if (priv_stropt(rad_cond_rei,    'rad_cond_rei',    RAD_COND_REI_OPT, rei_const) == RMN_ERR) return
      if (priv_stropt(rad_cond_rew,    'rad_cond_rew',    RAD_COND_REW_OPT, rew_const) == RMN_ERR) return
      if (priv_stropt(rad_conserve,    'rad_conserve',    RAD_CONSERVE_OPT) == RMN_ERR) return
      if (priv_stropt(linoz_chm,       'linoz_chm',       LINOZ_CHM_OPT) == RMN_ERR) return
      if (priv_stropt(iuv_method,      'iuv_method',      IUV_METHOD_OPT) == RMN_ERR) return
      if (priv_stropt(radia,           'radia',           RADIA_OPT) == RMN_ERR) return
      if (priv_stropt(stcond,          'stcond',          STCOND_OPT) == RMN_ERR) return
      if (priv_stropt(tofd,            'tofd',            TOFD_OPT) == RMN_ERR) return
      if (priv_stropt(lhn,             'lhn',             LHN_OPT) == RMN_ERR) return
      istat = clib_toupper(lhn_start_S)
      istat = clib_toupper(lhn_stop_S)
      istat = clib_toupper(lhn_ramp_S)

      if (gwdrag == 'NIL') then
         call msg(MSG_INFO, '(phy_nml_check) gwdrag = NIL : setting sgo_tdfilter = -1')
         sgo_tdfilter = -1
      endif
      if (lhn == 'NIL') then
         call msg(MSG_INFO, '(phy_nml_check) lhn = NIL : setting lhn_filter = -1')
         lhn_filter = -1
      endif

      if (fluvert /= "SURFACE" .and. any(pcptype == (/'SPS_W19', 'SPS_FRC', 'SPS_H13'/))) then
         call msg(MSG_ERROR, '(phy_nml_check) pcptype = '//trim(pcptype)//' can only be used with fluvert = "SURFACE" ')
         return
      end if

      if (stcond(1:3) /= "MP_") mpdiag_for_sfc = .false.

      if (.not.any(sfcflx_filter_order == (/ -1, 2, 4 /))) then
         call msg(MSG_ERROR, '(phy_nml_check) sfcflx_filter_order must be -1, 2 or 4')
         return
      endif

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function phy_nml_check


   function phy_nml_post_init() result(m_istat)
      use mixing_length, only: ml_init,ml_put
      implicit none
      integer :: m_istat

      integer :: istat
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      !# Copy cnv_options var to phy_options var
      convec = deep
      conv_shal = shal
      conv_mid = mid

      !# Set flags for memory debugging
      if (debug_mem_L) init2nan_L = .true.

      !# Mixing length module init
      if (ml_init(longmel, ilongmel) /= PHY_OK) then
         call msg_toall(MSG_ERROR,'(phy_nml) Problem with ml_init')
         return
      endif
      istat = ml_put('mlblac_max',pbl_mlblac_max)
      if (istat /= PHY_OK) then
         call msg(MSG_ERROR,'(phy_nml) cannot configure mixing length module')
         return
      endif

      !# Set Linoz flags
      if (linoz_chm /= 'NIL' .and. .not.any(radia == (/'CCCMARAD ', 'CCCMARAD2'/)) ) then
         call msg(MSG_ERROR,'(phy_nml) Cannot use Linoz_chm without CCCMARAD or CCCMARAD2')
         return
      endif
      if (rad_linoz_l .and. radia /= 'CCCMARAD2') then
         call msg(MSG_ERROR,'(phy_nml) Cannot use Linoz_rad_L=T without CCCMARAD2')
         return
      endif
      llinoz  = ( &
           any(radia /= (/ &
           'CCCMARAD ', &
           'CCCMARAD2'/)) .and. &
           any(linoz_chm == (/ &
           'OZONE   ', &
           'OZONEGHG' /)))
      llingh  = ( &
           any(radia /= (/ &
           'CCCMARAD ', &
           'CCCMARAD2'/)) .and. &
           any(linoz_chm == (/ &
           'GHG     ', &
           'OZONEGHG' /)))

      if (lhn_filter >= 0. .and. lhn == 'NIL')  then
         call msg(MSG_WARNING,'(phy_nml) lhn_filter ignored since lhn="NIL"')
         lhn_filter = -1.
      endif

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function phy_nml_post_init

   
   function priv_stropt(F_var, F_name, F_optlist, F_realval) result(m_istat)
      implicit none
      character(len=*), intent(INOUT) :: F_var
      character(len=*), intent(IN)    :: F_name
      character(len=*), intent(IN)    :: F_optlist(:)
      real, intent(INOUT), optional   :: F_realval
      integer :: m_istat

      integer :: istat
      character(len=1024) :: msg_S, msg1_S
      !----------------------------------------------------------------
      m_istat = RMN_ERR
      istat = clib_toupper(F_var)

      if (present(F_realval)) then
         msg1_S = 'Should be a number or one of'
         istat = str_toreal(F_realval, F_var)
      else
         msg1_S = 'Should be one of'
         istat = RMN_ERR
      endif
      if (.not.(any(F_var == F_optlist) .or. RMN_IS_OK(istat))) then
        call str_concat(msg_S, F_optlist, ', ')
        call msg(MSG_ERROR,'(phy_nml_check) '//trim(F_name)//' = '//trim(F_var)//' : '//trim(msg1_S)//': '//trim(msg_S))
        return
      endif
      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function priv_stropt
   
end module phy_nml_mod
