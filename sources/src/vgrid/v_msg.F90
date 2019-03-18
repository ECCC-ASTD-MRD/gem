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
#include "v_msg.h"

!/@*
module mod_v_msg
   implicit none
   public
   !@author  Stephane Chamberland, 2009-11
!@*/
   logical,save :: isInit_L = .false.
   logical,save :: canWrite_L = .true.
   integer,save :: v_msgLevelMin = V_MSG_INFO
   integer,save :: v_msgUnit(V_MSG_NBLEVELS0:V_MSG_NBLEVELS)
   character(len=V_MSG_MAXLEN),save :: v_msgFormat(V_MSG_NBLEVELS0:V_MSG_NBLEVELS) 
end module mod_v_msg


!/@*
subroutine v_msg_init()
   use mod_v_msg
   implicit none
   !@author  Stephane Chamberland, 2009-11
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) then
      isInit_L = .true.
      v_msgUnit(V_MSG_QUIET)   = V_MSG_STDOUT
      v_msgUnit(V_MSG_DEBUG)   = V_MSG_STDERR
      v_msgUnit(V_MSG_INFOPLUS)= V_MSG_STDOUT
      v_msgUnit(V_MSG_INFO)    = V_MSG_STDOUT
      v_msgUnit(V_MSG_IMPORTANT)= V_MSG_STDOUT
      v_msgUnit(V_MSG_WARNING) = V_MSG_STDOUT
      v_msgUnit(V_MSG_ERROR)   = V_MSG_STDERR
      v_msgUnit(V_MSG_CRITICAL)= V_MSG_STDERR
      v_msgUnit(V_MSG_VERBATIM)= V_MSG_STDOUT

      v_msgFormat(V_MSG_QUIET)   = '(a)'
      v_msgFormat(V_MSG_DEBUG)   = '("DEBUG: ",a)'
      v_msgFormat(V_MSG_INFOPLUS)= '(a)'
      v_msgFormat(V_MSG_INFO)    = '(a)'
      v_msgFormat(V_MSG_IMPORTANT)= '(a)'
      v_msgFormat(V_MSG_WARNING) = '("WARNING: ",a)'
      v_msgFormat(V_MSG_ERROR)   = '("ERROR: ",a)'
      v_msgFormat(V_MSG_CRITICAL)= '("CRITICAL ERROR: ",a)'
      v_msgFormat(V_MSG_VERBATIM)= '(a)'
   endif
   !---------------------------------------------------------------------
   return
end subroutine v_msg_init


!/@*
subroutine v_msg_set_p0only(F_myPE)
   use mod_v_msg, only: isInit_L,canWrite_L
   implicit none
   !@arguments
   integer,intent(in) :: F_myPE
   !@author  Stephane Chamberland, 2009-11
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   canWrite_L = .false.
   if (F_myPe == 0) canWrite_L = .true.
   !---------------------------------------------------------------------
   return
end subroutine v_msg_set_p0only


!/@*
subroutine v_msg_set_can_write(F_canWrite_L)
   use mod_v_msg, only: isInit_L,canWrite_L
   implicit none
   logical,intent(in) :: F_canWrite_L
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   canWrite_L = F_canWrite_L
   !---------------------------------------------------------------------
   return
end subroutine v_msg_set_can_write


!/@*
subroutine v_msg_verbosity(F_v_msgLevel)
   implicit none
   !@arguments
   integer,intent(in) :: F_v_msgLevel
   !@*/
   !---------------------------------------------------------------------
   call v_msg_set_minMessageLevel(F_v_msgLevel)
   !---------------------------------------------------------------------
   return
end subroutine v_msg_verbosity


!/@*
subroutine v_msg_verbosity_get(F_v_msgLevel)
   use mod_v_msg
   implicit none
   !@arguments
   integer,intent(out) :: F_v_msgLevel
   !@*/
   !---------------------------------------------------------------------
   F_v_msgLevel = v_msgLevelMin
   !---------------------------------------------------------------------
   return
end subroutine v_msg_verbosity_get


!/@*
subroutine v_msg_set_minMessageLevel(F_v_msgLevel)
   use mod_v_msg, only: isInit_L,v_msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_v_msgLevel
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   v_msgLevelMin = max(V_MSG_NBLEVELS0,min(F_v_msgLevel,V_MSG_NBLEVELS))
   !---------------------------------------------------------------------
   return
end subroutine v_msg_set_minMessageLevel


!/@*
subroutine v_msg_set_redirect2fileUnit(F_out_or_err,F_fileUnit)
   use mod_v_msg, only: isInit_L,v_msgUnit
   implicit none
   !@arguments
   integer,intent(in) :: F_out_or_err
   integer,intent(in) :: F_fileUnit
   !@*/
   integer :: i, i0=1, i1=V_MSG_NBLEVELS
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   if (F_out_or_err == V_MSG_STDERR) then
      i0=V_MSG_ERROR
   else if (F_out_or_err == V_MSG_STDOUT) then
      i1=V_MSG_WARNING
   else
      i0 = F_out_or_err
      i1 = F_out_or_err
   endif
   do i=i0,i1
      v_msgUnit(i) = F_fileUnit
   enddo
   !---------------------------------------------------------------------
   return
end subroutine v_msg_set_redirect2fileUnit


!/@*
subroutine v_msg_toall(F_v_msgLevel,F_Message)
   use mod_v_msg, only: isInit_L,canWrite_L,v_msgUnit,v_msgFormat,v_msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_v_msgLevel
   character(len=*),intent(in) :: F_Message
   !@author  Stephane Chamberland, 2009-11
   !@*/
   integer :: v_msgLevel
   logical :: canWrite_L_bk
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   canWrite_L_bk = canWrite_L
   canWrite_L = .true.
   call v_msg(F_v_msgLevel,F_Message)
   canWrite_L = canWrite_L_bk
   !---------------------------------------------------------------------
   return
end subroutine v_msg_toall


!/@*
subroutine v_msg(F_v_msgLevel,F_Message)
   use mod_v_msg, only: isInit_L,canWrite_L,v_msgUnit,v_msgFormat,v_msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_v_msgLevel
   character(len=*),intent(in) :: F_Message
   !@author  Stephane Chamberland, 2009-11
   !@*/
   integer :: v_msgLevel
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   if (canWrite_L .and. F_v_msgLevel>=v_msgLevelMin) then
      v_msgLevel = max(V_MSG_NBLEVELS0,min(F_v_msgLevel,V_MSG_NBLEVELS))
      !TODO: use C's fprintf so that it can be called from C
      write(v_msgUnit(v_msgLevel),trim(v_msgFormat(v_msgLevel))) trim(F_Message)
      call flush(v_msgUnit(v_msgLevel))
   endif
   !---------------------------------------------------------------------
   return
end subroutine v_msg


!/@*
function v_msg_getUnit(F_v_msgLevel) result(F_v_msgUnit)
   use mod_v_msg, only: isInit_L,canWrite_L,v_msgUnit,v_msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_v_msgLevel
   !@returns
   integer :: F_v_msgUnit
   !@author  Stephane Chamberland, 2009-11
   !@*/
   integer :: v_msgLevel
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   F_v_msgUnit = -1
   if (canWrite_L .and. F_v_msgLevel>=v_msgLevelMin) then
      v_msgLevel = max(V_MSG_NBLEVELS0,min(F_v_msgLevel,V_MSG_NBLEVELS))
      F_v_msgUnit = v_msgUnit(v_msgLevel)
   endif
   !---------------------------------------------------------------------
   return
end function v_msg_getUnit


!/@*
subroutine v_msg_getInfo(F_canWrite_L,F_v_msgLevelMin,F_v_msgUnit,F_v_msgFormat_S)
   use mod_v_msg
   implicit none
   !@arguments
   logical,intent(out) :: F_canWrite_L
   integer,intent(out) :: F_v_msgLevelMin,F_v_msgUnit
   character(len=*),intent(out) :: F_v_msgFormat_S
   !@author  Stephane Chamberland, 2012-02
   !@*/
   integer,external :: v_msg_getUnit
   !---------------------------------------------------------------------
   if (.not.isInit_L) call v_msg_init()
   F_canWrite_L = canWrite_L
   F_v_msgLevelMin = v_msgLevelMin
   F_v_msgLevelMin = max(V_MSG_NBLEVELS0,min(v_msgLevelMin,V_MSG_NBLEVELS))
   F_v_msgUnit = v_msgUnit(v_msgLevelMin)
   F_v_msgFormat_S = v_msgFormat(v_msgLevelMin)
   !---------------------------------------------------------------------
   return
end subroutine v_msg_getInfo
