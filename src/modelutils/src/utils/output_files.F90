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

!/@
module output_files_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_isdir, clib_mkdir, clib_realpath
   use wb_itf_mod
   use outcfg_mod
   use fstmpi_mod
   use ptopo_utils
   use cmcdate_mod
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland, 2014-12
   !@description
   ! Public functionsw
   public :: output_files_set_basedir, output_files_get_basedir, &
        output_files_close,output_files_open,output_files_set_postproc, &
        output_files_get_closestep,output_files_filename
   ! Public constants
   !
   character(len=*),parameter,public :: OUTPUT_FILES_WBFMT='("O#/",i6.6,"/",a)'
!@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   integer,save :: m_postprocstep(OUTCFG_NCFG_MAX) = -1

contains

   !/@*
   function output_files_set_basedir(F_id,F_basedir_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_basedir_S
      !@return
      integer :: F_istat
      !*@/
      character(len=WB_MAXSTRINGLENGTH) :: basedir_S,name_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(output) set_basedir [BEGIN]')
      F_istat = clib_isdir(trim(F_basedir_S))
      if (.not.RMN_IS_OK(F_istat)) then
         F_istat = clib_mkdir(trim(F_basedir_S))
      endif
      F_istat = clib_realpath(trim(F_basedir_S),basedir_S)
      F_istat = clib_isdir(trim(basedir_S))
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(output) set_basedir, not such dir: '//trim(basedir_S))
         return
      endif
      write(name_S,OUTPUT_FILES_WBFMT) F_id,'basedir'
      F_istat = wb_put(name_S,basedir_S,WB_REWRITE_MANY)
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFO,'(output) set_basedir: '//trim(basedir_S))
      else
         call msg(MSG_WARNING,'(output) Probleme saving basedir: '//trim(basedir_S))
      endif
      call msg(MSG_DEBUG,'(output) set_basedir [END]')
      !----------------------------------------------------------------------
      return
   end function output_files_set_basedir


   !/@*
   function output_files_get_basedir(F_id,F_basedir_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(out) :: F_basedir_S
      !@return
      integer :: F_istat
      !*@/
      character(len=WB_MAXSTRINGLENGTH) :: basedir_S,name_S
      !----------------------------------------------------------------------
      write(name_S,OUTPUT_FILES_WBFMT) F_id,'basedir'
      F_istat = wb_get(name_S,basedir_S)
      if (basedir_S == ' ' .or. .not.RMN_IS_OK(F_istat)) then
         F_istat = clib_realpath(".",basedir_S)
      endif
      F_basedir_S = basedir_S
      !----------------------------------------------------------------------
      return
   end function output_files_get_basedir


   !/@*
   function output_files_close(F_id,F_step) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-12
      !*@/
      integer :: fileid,istat,mystep,itype
      character(len=1024) :: msg_S
      !----------------------------------------------------------------------
      F_istat = RMN_OK
      mystep = output_files_get_closestep(F_id,F_step)

      do itype=1,size(OUTCFG_LVLTYPE)
         fileid = priv_get_fileid(F_id,mystep,OUTCFG_LVLTYPE(itype))
         if (RMN_IS_OK(fileid)) then
            write(msg_S,'(i6,1x,a)') fileid
            call msg(MSG_DEBUG,'(Output) Close file:'//trim(msg_S))
            F_istat = fstmpi_close(fileid)
            istat = priv_set_fileid(F_id,mystep,OUTCFG_LVLTYPE(itype),RMN_ERR)
         endif
      enddo
      !----------------------------------------------------------------------
      return
   end function output_files_close


   !/@*
   function output_files_open(F_id,F_step,F_lvltyp_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      character(len=*),intent(in),optional :: F_lvltyp_S
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-12
      !*@/
      integer :: istat,mystep,itype,fileid,nlvltyp
      character(len=1024) :: filename_S,msg_S
      character(len=1) :: lvltyp_S(2)
      !----------------------------------------------------------------------
      mystep = output_files_get_closestep(F_id,F_step)
      F_istat = RMN_OK
      if (present(F_lvltyp_S)) then
         !TODO-later: if (any(F_lvltyp_S(1:1) == OUTCFG_LVLTYPE(:))) then
         nlvltyp = 1
         lvltyp_S(1) = F_lvltyp_S(1:1)
      else
         nlvltyp = size(OUTCFG_LVLTYPE)
         lvltyp_S(:) = OUTCFG_LVLTYPE(:)
      endif
      fileid = RMN_ERR
      do itype=1,nlvltyp
         fileid = priv_get_fileid(F_id,mystep,lvltyp_S(itype))
         if (RMN_IS_OK(fileid)) cycle
         !TODO-later: try to close previous file?

         istat = output_files_filename(F_id,mystep,lvltyp_S(itype),filename_S)
         if (RMN_IS_OK(istat)) then
            fileid = fstmpi_open(filename_S)
            if (RMN_IS_OK(fileid)) then
               write(msg_S,'(i6,1x,a)') fileid,trim(filename_S)
               call msg(MSG_DEBUG,'(Output) Open file:'//trim(msg_S))
               istat = priv_set_fileid(F_id,mystep,lvltyp_S(itype),fileid)
            endif
         endif
         F_istat = min(istat,fileid,F_istat)
      enddo
      if (nlvltyp == 1 .and. RMN_IS_OK(F_istat) .and. RMN_IS_OK(fileid)) F_istat = fileid
      !----------------------------------------------------------------------
      return
   end function output_files_open


   !/@*
   function output_files_set_postproc(F_id,F_step) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2012-01
      !*@/
      !----------------------------------------------------------------------
      F_istat = RMN_OK
      if (F_id > 0 .and. F_id <= size(m_postprocstep)) &
           m_postprocstep(F_id) = F_step
      !----------------------------------------------------------------------
      return
   end function output_files_set_postproc


   !/@*
   function output_files_get_closestep(F_id,F_step) result(F_closestep)
      implicit none
      !@arguments
      integer,intent(in) :: F_id,F_step
      !@return
      integer :: F_closestep
      !*@/
      character(len=64) :: tag_S,units_S
      integer :: closeopen,ndigits,istat,closeopen_units
      !----------------------------------------------------------------------
      istat = outcfg_meta(F_id,tag_S,units_S,closeopen,ndigits,F_closeopen_units=closeopen_units)
      if (closeopen_units == OUTCFG_UNITS_1STOF_MONTH) units_S = 'month'
      if (units_S(1:2) == 'mo') then
         F_closestep = F_step
         return
      endif

      closeopen = max(1,closeopen)
      if (.not.RMN_IS_OK(istat) .or. F_step == 0) then
         F_closestep = F_step
      else if (closeopen <= 1 .and. m_postprocstep(F_id) <= 1) then
         F_closestep = F_step
      else
         F_closestep = floor(float(F_step)/float(closeopen)) * closeopen
         if (m_postprocstep(F_id) > 0) F_closestep = max(F_closestep,&
              floor(float(F_step)/float(m_postprocstep(F_id))) * m_postprocstep(F_id))
      endif
      !----------------------------------------------------------------------
      return
   end function output_files_get_closestep


   !/@*
   function output_files_filename(F_id,F_closestep,F_lvltyp_S,F_filename_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_closestep
      character(len=*),intent(in) :: F_lvltyp_S
      character(len=*),intent(out) :: F_filename_S
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2011-12
      !*@/
      integer,external :: prog_filename
      real(REAL64),parameter :: EPSILON_8 = 1.D-12
      integer :: istat,datev,dateo,dt,closeopen,ndigits,prognum,pdate,ptime,hour,minutes,seconds,closeopen_units,delta_year
      character(len=64) :: tag_S,units_S,ext_s
      character(len=1024) :: dirname_S
      !----------------------------------------------------------------------
      call ptopo_init_var()
      F_istat = RMN_ERR
      F_filename_S = ' '

      istat = outcfg_time(F_id,F_closestep,datev,dateo,dt)
      istat = outcfg_meta(F_id,tag_S,units_S,closeopen,ndigits,F_closeopen_units=closeopen_units)

      if (closeopen_units == OUTCFG_UNITS_1STOF_MONTH) units_S = 'month'

      select case(units_S(1:2))
      case('st')
         prognum = F_closestep
         units_S = 'pas'
      case('da')
         prognum = ceiling(dble(F_closestep)*dble(dt)/86400.D0 - EPSILON_8)
      case('mi')
         prognum = ceiling(dble(F_closestep)*dble(dt)/60.D0 - EPSILON_8)
      case('se')
         prognum = ceiling(dble(F_closestep)*dble(dt) - EPSILON_8)
      case('mo')
         delta_year = cmcdate_year(datev) - cmcdate_year(dateo)
         prognum = cmcdate_month(datev) - cmcdate_month(dateo) + 12*delta_year
         units_S = 'Month'
      case default !hour
         prognum = ceiling(dble(F_closestep)*dble(dt)/3600.D0 - EPSILON_8)
      end select

      istat = newdate(dateo,pdate,ptime,RMN_DATE_STAMP2PRINT)
      hour = ptime/1000000
      minutes = (ptime - hour*1000000)/10000
      seconds = (ptime - hour*1000000 - minutes*10000)/100

      !#TODO: use timestr_prognum
      tag_s = trim(tag_S)//trim(F_lvltyp_S)
      F_filename_S = ' '
      istat = prog_filename(F_filename_S,tag_S,pdate,hour,minutes,seconds, &
           ptopo_grid_iblocx,ptopo_grid_iblocy,prognum,ndigits,units_S(1:1))
      if (RMN_IS_OK(istat)) then
         ext_S = ' ' !TODO-later: if init ext_S='_dgf'
         F_istat = output_files_get_basedir(F_id,dirname_S)
         F_filename_S = trim(dirname_S)//'/'//trim(F_filename_S)//trim(ext_S)
         F_istat = RMN_OK
      endif
      !----------------------------------------------------------------------
      return
   end function output_files_filename


   !==== Private Functions =================================================


   !/@*
   function priv_set_fileid(F_id,F_step,F_lvltyp_S,F_fileid) result(F_istat)
      implicit none
      !@arguments
      integer,intent(in) :: F_id,F_step,F_fileid
      character(len=*),intent(in) :: F_lvltyp_S
      !@return
      integer :: F_istat
      !*@/
      character(len=WB_MAXNAMELENGTH) :: keyname_S
      !----------------------------------------------------------------------
      write(keyname_S,'(a,i4.4,a,i6.6,a)') 'OUTFID/',F_id,'-',F_step,trim(F_lvltyp_S)
      F_istat = wb_put(keyname_S,F_fileid,WB_REWRITE_MANY)
      !----------------------------------------------------------------------
      return
   end function priv_set_fileid


   !/@*
   function priv_get_fileid(F_id,F_step,F_lvltyp_S) result(F_fileid)
      implicit none
      !@arguments
      integer,intent(in) :: F_id,F_step
      character(len=*),intent(in) :: F_lvltyp_S
      !@return
      integer :: F_fileid
      !*@/
      character(len=WB_MAXNAMELENGTH) :: keyname_S
      integer :: istat,iverb
      !----------------------------------------------------------------------
      write(keyname_S,'(a,i4.4,a,i6.6,a)') 'OUTFID/',F_id,'-',F_step,trim(F_lvltyp_S)
      iverb = wb_verbosity(WB_MSG_FATAL)
      istat = wb_get(keyname_S,F_fileid)
      iverb = wb_verbosity(iverb)
      if (.not.RMN_IS_OK(istat)) F_fileid = RMN_ERR
      !----------------------------------------------------------------------
      return
   end function priv_get_fileid

end module output_files_mod
