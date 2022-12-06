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

!/@*
function cnv_nml2(F_namelist) result(F_istat)
   use clib_itf_mod, only: clib_isreadok, clib_toupper
   use tdpack_const
   use str_mod, only: str_concat,str_toreal
   use cnv_options
   implicit none
!!!#include <arch_specific.hf>
   !@Object Set defaults values and read convection namelist
   !@Arguments
   ! F_namelist    File name containing the namelists to read
   ! F_istat       1 if no error occurs, -1 otherwise
   character(len=*), intent(in) :: F_namelist
   integer :: F_istat
   !*@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   integer, parameter :: CNV_NML_ERR = RMN_ERR
   integer, parameter :: CNV_NML_OK  = RMN_OK + 1

   integer :: err
   !-------------------------------------------------------------------
   F_istat = CNV_NML_ERR

   err = cnv_options_init()
   if (.not.RMN_IS_OK(err)) return
   err = cnv_nml_read(F_namelist)
   if (.not.RMN_IS_OK(err)) return
   err = cnv_nml_check()
   if (.not.RMN_IS_OK(err)) return
   call cnv_nml_print()
   err = cnv_nml_post_init()
   if (.not.RMN_IS_OK(err)) return

   F_istat = CNV_NML_OK
  !-------------------------------------------------------------------
  return


contains


   function cnv_nml_read(m_namelist) result(m_istat)
      implicit none
      character(len=*), intent(in) :: m_namelist
      integer :: m_istat
      character(len=1024) :: msg_S, namelist_S, name_S
      integer :: istat, unf
      !----------------------------------------------------------------
      m_istat = RMN_ERR
      name_S = '(cnv_nml)'
      namelist_S = '&convection_cfgs'

      unf = 0
      istat = clib_isreadok(m_namelist)
      if (RMN_IS_OK(istat)) istat = fnom(unf, m_namelist, 'SEQ+OLD', 0)
      if (.not.RMN_IS_OK(istat)) then
         write(msg_S,'(a,a,a,a)') trim(name_S), &
              ' Using default config, ', &
              'Namelist file Not found/readable: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
         return
      endif

      read(unf, nml=convection_cfgs, iostat=istat)
      if (istat == 0) then     !# Read ok
         m_istat = RMN_OK + 1
      else if (istat < 0) then !# EOF, nml not found
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Using default config, Namelist ',&
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
   end function cnv_nml_read


   subroutine cnv_nml_print()
      implicit none
      integer, external :: msg_getUnit
      integer :: unout
      !----------------------------------------------------------------
      unout = msg_getUnit(MSG_INFO)
      if (unout > 0) then
         write(unout, nml=convection_cfgs)
      endif
      !----------------------------------------------------------------
      return
   end subroutine cnv_nml_print


   function cnv_nml_check() result(m_istat)
      use timestr_mod, only: timestr_check
      implicit none
      integer :: m_istat

      integer :: istat
      real :: mid_peff_tmp
      character(len=1024) :: msg_S
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      istat = clib_toupper(cmt_type)
      istat = clib_toupper(deep)
      istat = clib_toupper(deep_conserve)
      istat = clib_toupper(deep_timeent)
      istat = clib_toupper(deep_timeconv)
      istat = clib_toupper(deep_timerefresh)
      istat = clib_toupper(mid)
      istat = clib_toupper(mid_conserve)
      istat = clib_toupper(mid_emffrac)
      istat = clib_toupper(mid_emfmod)
      istat = clib_toupper(mid_peff)
      istat = clib_toupper(shal)
      istat = clib_toupper(shal_conserve)
      istat = clib_toupper(shal_timeconv)
      istat = clib_toupper(bkf_entrains)
      istat = clib_toupper(bkf_detrains)
      istat = clib_toupper(bkf_closures)

      if (.not.any(cmt_type == CMT_TYPE_OPT)) then
         call str_concat(msg_S,CMT_TYPE_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) cmt_type = '//trim(cmt_type)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(deep == DEEP_OPT)) then
         call str_concat(msg_S,DEEP_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) deep = '//trim(deep)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(deep_conserve == DEEP_CONSERVE_OPT)) then
         call str_concat(msg_S,DEEP_CONSERVE_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) deep_conserve = '//trim(deep_conserve)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(deep_timeent == DEEP_TIMEENT_OPT(min(2,size(DEEP_TIMEENT_OPT)):))) then
         if (RMN_IS_OK(timestr_check(deep_timeent))) then
            deep_timeent_sec = 1. !final update will be done in phy_init when timestep information is available
         else
            call str_concat(msg_S,DEEP_TIMEENT_OPT,', ')
            call msg(MSG_ERROR,'(cnv_nml_check) deep_timeent = '//trim(deep_timeent)//' : Should be one of: '//trim(msg_S))
            return
         end if
      end if

      if (.not.any(deep_timeconv == DEEP_TIMECONV_OPT(min(2,size(DEEP_TIMECONV_OPT)):))) then
         if (RMN_IS_OK(timestr_check(deep_timeconv))) then
            deep_timeconv_sec = 1. !final update will be done in phy_init when timestep information is available
         else
            call str_concat(msg_S,DEEP_TIMECONV_OPT,', ')
            call msg(MSG_ERROR,'(cnv_nml_check) deep_timeconv = '//trim(deep_timeconv)//' : Should be one of: '//trim(msg_S))
            return
         end if
      end if

      if (.not.any(deep_timerefresh == DEEP_TIMEREFRESH_OPT(min(2,size(DEEP_TIMEREFRESH_OPT)):))) then
         if (RMN_IS_OK(timestr_check(deep_timerefresh))) then
            deep_timerefresh_sec = 1. !final update will be done in phy_init when timestep information is available
         else
            call str_concat(msg_S,DEEP_TIMEREFRESH_OPT,', ')
            call msg(MSG_ERROR,'(cnv_nml_check) deep_timerefresh = '//trim(deep_timerefresh)//' : Should be one of: '//trim(msg_S))
            return
         end if
      end if

      if (.not. any(MID == MID_OPT)) then
         call str_concat(msg_S,MID_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) mid = '//trim(mid)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(mid_conserve == MID_CONSERVE_OPT)) then
         call str_concat(msg_S,MID_CONSERVE_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) mid_conserve = '//trim(mid_conserve)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(mid_emffrac == MID_EMFFRAC_OPT)) then
         call str_concat(msg_S,MID_EMFFRAC_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) mid_emffrac = '//trim(mid_emffrac)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(mid_emfmod == MID_EMFMOD_OPT)) then
         call str_concat(msg_S,MID_EMFMOD_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) mid_emfmod = '//trim(mid_emfmod)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (RMN_IS_OK(str_toreal(mid_peff_tmp,mid_peff))) then
         mid_peff_const = mid_peff_tmp
      else
         if (.not.any(mid_peff == MID_PEFF_OPT)) then
            call str_concat(msg_S,MID_PEFF_OPT,', ')
            call msg(MSG_ERROR,'(cnv_nml_check) mid_peff = '//trim(mid_peff)//' : Should be one of: '//trim(msg_S))
            return
         end if
      endif

      if (.not. any(SHAL == SHAL_OPT)) then
         call str_concat(msg_S,SHAL_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) shal = '//trim(shal)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(shal_conserve == SHAL_CONSERVE_OPT)) then
         call str_concat(msg_S,SHAL_CONSERVE_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) shal_conserve = '//trim(shal_conserve)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(shal_timeconv == SHAL_TIMECONV_OPT(min(2,size(SHAL_TIMECONV_OPT)):))) then
         if (RMN_IS_OK(timestr_check(shal_timeconv))) then
            shal_timeconv_sec = 1. !final update will be done in phy_init when timestep information is available
         else
            call str_concat(msg_S,SHAL_TIMECONV_OPT,', ')
            call msg(MSG_ERROR,'(cnv_nml_check) shal_timeconv = '//trim(shal_timeconv)//' : Should be one of: '//trim(msg_S))
            return
         end if
      end if

      if (.not.any(bkf_entrains == BKF_ENTRAINS_OPT)) then
         call str_concat(msg_S,BKF_ENTRAINS_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) bkf_entrains = '//trim(bkf_entrains)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(bkf_detrains == BKF_DETRAINS_OPT)) then
         call str_concat(msg_S,BKF_DETRAINS_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) bkf_detrains = '//trim(bkf_detrains)//' : Should be one of: '//trim(msg_S))
         return
      end if

      if (.not.any(bkf_closures == BKF_CLOSURES_OPT)) then
         call str_concat(msg_S,BKF_CLOSURES_OPT,', ')
         call msg(MSG_ERROR,'(cnv_nml_check) bkf_closures = '//trim(bkf_closures)//' : Should be one of: '//trim(msg_S))
         return
      end if

      !# check if triglat is valid
      if ( (min(triglat(1),triglat(2)) <  0.0) .or. &
           (max(triglat(1),triglat(2)) > 90.0) .or. &
           (triglat(1) > triglat(2)) ) then
         call msg(MSG_ERROR,'(cnv_nml_check) Wrong values for triglat')
         return
      end if

      !# check if bkf_rads is valid
      if ( (min(bkf_rads(1),bkf_rads(2)) <  10.0) .or. &
           (bkf_rads(2)-bkf_rads(1) < 0.0)  .or. &
           (bkf_rads(3) < 0.0) ) then
         call msg(MSG_ERROR,'(cnv_nml_check) Wrong values for bkf_rads; (1-2) shoud be .GE. 10 m')
         return
      endif

      !# check if bkf_tperts is valid
      if ( (min(bkf_tperts(1),bkf_tperts(2),bkf_tperts(3)) < 0.0) .or. &
           (bkf_tperts(2)-bkf_tperts(1) < 0.0) ) then
         call msg(MSG_ERROR,'(cnv_nml_check) Wrong values for bkf_tperts')
         return
      endif

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function cnv_nml_check


   function cnv_nml_post_init() result(m_istat)
      implicit none
      integer :: m_istat
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      !# conversion from degrees to radians
      triglat(1) = triglat(1) * PI/180.
      triglat(2) = triglat(2) * PI/180.

      !# set convective radius over water if not specified by user
      if (kfcradw < 0.) kfcradw = kfcrad

      !# cmt options
      select case(cmt_type)
      case ('ECMWF_PH2')
         cmt_type_i = CMT_ECMWF_PH2
         cmt_gki_cup = 0.
         cmt_gki_cdown = 0.
      case ('ECMWF')
         cmt_type_i = CMT_ECMWF
         cmt_gki_cup = 0.
         cmt_gki_cdown = 0.
      case ('GKI')
         cmt_type_i = CMT_GKI
         cmt_ecmwf_lambda = 0.
      case default
         cmt_type_i = CMT_NONE
         cmt_ecmwf_lambda = 0.
         cmt_gki_cup = 0.
         cmt_gki_cdown = 0.
      end select
     
      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function cnv_nml_post_init

end function cnv_nml2
