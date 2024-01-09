!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module phy_nml_mod
   use clib_itf_mod, only: clib_isreadok, clib_toupper
   use debug_mod, only: init2nan_L
   use str_mod, only: str_concat, str_toreal
   use series_mod, only: series_nml
   use phy_status, only: PHY_ERROR, PHY_NONE, PHY_OK, PHY_CTRL_NML_OK, phy_init_ctrl
   use phy_options
   use cnv_options
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

      integer, external :: chm_nml, sfc_nml2, cnv_nml2, check_options2, msg_getUnit

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

      !# Read surface namelist
      err = sfc_nml2(F_namelist)
      if (.not.RMN_IS_OK(err)) return

      !# Read convection namelist
      err = cnv_nml2(F_namelist)
      if (.not.RMN_IS_OK(err)) return

      !# Read series namelists
      err = series_nml(F_namelist)
      if (.not.RMN_IS_OK(err)) return

      !# Read chemistry namelists and initialize chemistry configuration
      err = chm_nml(F_namelist, unout)
      if (err <= 0) then
         call msg(MSG_ERROR, '(phy_nml) Problem reading chemestry namelist')
         return
      endif
      
      !# Check namelist options validity and consistency
      err = phy_nml_check()
      if (.not.RMN_IS_OK(err)) return

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

      istat = clib_toupper(cond_conserve)
      istat = clib_toupper(fluvert)
      istat = clib_toupper(gwdrag)
      !#TODO: indiag_list_s?
      istat = clib_toupper(input_type)
      istat = clib_toupper(iuv_method)
      istat = clib_toupper(kntrad_S)
      istat = clib_toupper(kntraduv_S)
      istat = clib_toupper(linoz_chm)
      istat = clib_toupper(lhn_ramp_S)
      istat = clib_toupper(lhn_start_S)
      istat = clib_toupper(lhn_stop_S)
      istat = clib_toupper(longmel)
      istat = clib_toupper(pbl_diss)
      istat = clib_toupper(pbl_dissheat)
      istat = clib_toupper(pbl_conserve)
      istat = clib_toupper(pbl_func_stab)
      istat = clib_toupper(pbl_func_unstab)
      istat = clib_toupper(pbl_mlblac_max)
      istat = clib_toupper(pbl_nonloc)
      istat = clib_toupper(pbl_shal)
      istat = clib_toupper(pcptype)
      !#TODO: phystat_list_s?
      istat = clib_toupper(rad_atmpath)
      istat = clib_toupper(rad_cond_rei)
      istat = clib_toupper(rad_cond_rew)
      istat = clib_toupper(rad_conserve)
      istat = clib_toupper(rad_part_nomp)
      istat = clib_toupper(radia)
      istat = clib_toupper(stcond)
      istat = clib_toupper(tofd)

      if (.not.any(cond_conserve == COND_CONSERVE_OPT)) then
         call str_concat(msg_S,COND_CONSERVE_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) cond_conserve = '//trim(cond_conserve)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(fluvert == FLUVERT_OPT)) then
         call str_concat(msg_S,FLUVERT_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) fluvert = '//trim(fluvert)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(gwdrag == GWDRAG_OPT)) then
         call str_concat(msg_S,GWDRAG_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) gwdrag = '//trim(gwdrag)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(input_type == INPUT_TYPE_OPT)) then
         call str_concat(msg_S, INPUT_TYPE_OPT, ', ')
         call msg(MSG_ERROR,'(phy_nml_check) input_type = '//trim(input_type)//' : Should be one of: '//trim(msg_S))
         return
      endif
      
      if (.not.any(iuv_method == IUV_METHOD_OPT)) then
         call str_concat(msg_S, IUV_METHOD_OPT, ', ')
         call msg(MSG_ERROR,'(phy_nml_check) iuv_method = '//trim(iuv_method)//' : Should be one of: '//trim(msg_S))
         return
      endif
      
      if (.not.any(linoz_chm == LINOZ_CHM_OPT)) then
         call str_concat(msg_S, LINOZ_CHM_OPT, ', ')
         call msg(MSG_ERROR,'(phy_nml_check) linoz_chm = '//trim(linoz_chm)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(lhn == LHN_OPT)) then
         call str_concat(msg_S,LHN_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) lhn = '//trim(lhn)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(longmel == ML_CLOSURES(:)%name)) then
         call str_concat(msg_S,ML_CLOSURES(:)%name,', ')
         call msg(MSG_ERROR,'(phy_nml_check) longmel = '//trim(longmel)//' : Should be one of: '//trim(msg_S))
      endif

      if (.not.any(pbl_conserve == PBL_CONSERVE_OPT)) then
         call str_concat(msg_S,PBL_CONSERVE_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_conserve = '//trim(pbl_conserve)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(pbl_diss == PBL_DISS_OPT)) then
         call str_concat(msg_S,PBL_DISS_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_diss = '//trim(pbl_diss)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pbl_dissheat == PBL_DISSHEAT_OPT)) then
         call str_concat(msg_S,PBL_DISSHEAT_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_dissheat = '//trim(pbl_dissheat)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pbl_func_stab == PBL_FUNC_STAB_OPT)) then
         call str_concat(msg_S,PBL_FUNC_STAB_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_func_stab = '//trim(pbl_func_stab)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pbl_func_unstab == PBL_FUNC_UNSTAB_OPT)) then
         call str_concat(msg_S,PBL_FUNC_UNSTAB_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_func_unstab = '//trim(pbl_func_unstab)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pbl_mlblac_max == PBL_MLBLAC_MAX_OPT)) then
         call str_concat(msg_S,PBL_MLBLAC_MAX_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_mlblac_max = '//trim(pbl_mlblac_max)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pbl_nonloc == PBL_NONLOC_OPT)) then
         call str_concat(msg_S,PBL_NONLOC_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_nonloc = '//trim(pbl_nonloc)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pbl_shal == PBL_SHAL_OPT)) then
         call str_concat(msg_S,PBL_SHAL_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pbl_shal = '//trim(pbl_shal)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(pcptype == PCPTYPE_OPT)) then
         call str_concat(msg_S,PCPTYPE_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) pcptype = '//trim(pcptype)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (fluvert /= "SURFACE" .and. any(pcptype == (/'SPS_W19', 'SPS_FRC', 'SPS_H13'/))) then
         call msg(MSG_ERROR,'(phy_nml_check) pcptype = '//trim(pcptype)//' can only be used with fluvert = "SURFACE" ')
         return
      end if

      if (stcond(1:3) /= "MP_") mpdiag_for_sfc = .false.

      if (.not.any(rad_atmpath == RAD_ATMPATH_OPT)) then
         call str_concat(msg_S,RAD_ATMPATH_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) rad_atmpath = '//trim(rad_atmpath)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(rad_part_nomp == RAD_PART_NOMP_OPT)) then
         call str_concat(msg_S,RAD_PART_NOMP_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) rad_part_nomp = '//trim(rad_part_nomp)//' : Should be one of: '//trim(msg_S))
         return
      endif
 

      if (.not.RMN_IS_OK(str_toreal(rei_const,rad_cond_rei))) then
         if (.not.any(rad_cond_rei == RAD_COND_REI_OPT)) then
            call str_concat(msg_S,RAD_COND_REI_OPT,', ')
            call msg(MSG_ERROR,'(phy_nml_check) rad_cond_rei = '//trim(rad_cond_rei)//' : Should be a number (in microns) or one of: '//trim(msg_S))
            return
         endif
      endif
 
      if (.not.RMN_IS_OK(str_toreal(rew_const,rad_cond_rew))) then
         if (.not.any(rad_cond_rew == RAD_COND_REW_OPT)) then
            call str_concat(msg_S,RAD_COND_REW_OPT,', ')
            call msg(MSG_ERROR,'(phy_nml_check) rad_cond_rew = '//trim(rad_cond_rew)//' : Should be a number (in microns) or one of: '//trim(msg_S))
            return
         endif
      endif

      if (.not.any(rad_conserve == RAD_CONSERVE_OPT)) then
         call str_concat(msg_S,RAD_CONSERVE_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) rad_conserve = '//trim(rad_conserve)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(radia == RADIA_OPT)) then
         call str_concat(msg_S,RADIA_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) radia = '//trim(radia)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(stcond == STCOND_OPT)) then
         call str_concat(msg_S,STCOND_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) stcond = '//trim(stcond)//' : Should be one of: '//trim(msg_S))
         return
      endif

     if (.not.any(tofd == TOFD_OPT)) then
         call str_concat(msg_S,TOFD_OPT,', ')
         call msg(MSG_ERROR,'(phy_nml_check) tofd = '//trim(tofd)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(sfcflx_filter_order == (/ -1, 2, 4 /))) then
         call msg(MSG_ERROR,'(phy_nml_check) sfcflx_filter_order must be -1, 2 or 4')
         return
      endif
      
      if (advectke .and. deep == 'KFC2') then
         call msg(MSG_ERROR,'(phy_nml_check) Cannot use ADVECTKE with KFC2 (fix pending)')
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

      if (STCOND(1:2) /= 'MP') then
         IOPTIX = OPT_OPTIX_OLD
      else
         IOPTIX = OPT_OPTIX_NEW
      endif

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

end module phy_nml_mod
