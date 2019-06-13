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
#include <msg.h>

!/@*
module mod_msg
   implicit none
   public
   !@author  Stephane Chamberland, 2009-11
!@*/
   logical,save :: isInit_L = .false.
   logical,save :: canWrite_L = .true.
   integer,save :: msgLevelMin = MSG_INFO
   integer,save :: msgUnit(MSG_NBLEVELS0:MSG_NBLEVELS)
   character(len=MSG_MAXLEN),save :: msgFormat(MSG_NBLEVELS0:MSG_NBLEVELS)
end module mod_msg


!/@*
module mod_msg_buffer
   !@author  Stephane Chamberland, 2018-06
   use mod_msg, only: msgUnit, msgFormat
   implicit none
   private
   public :: buffer_verbosity, buffer_add, buffer_reset, buffer_flush
!@*/
   integer, parameter :: MSG_NMAXBUFFER = 256
   integer,save :: buffer_min_level = MSG_ERROR
   integer,save :: nBuffer = 0
   integer,save :: levelBuffer(MSG_NMAXBUFFER)
   character(len=MSG_MAXLEN),save :: msgBuffer(MSG_NMAXBUFFER)
   !#TODO: offer option to turn off buffering if it's too slooow

contains

   !/@*
   function buffer_verbosity(F_msgLevel) result(F_msgLevel0)
      implicit none
      !@arguments
      integer,intent(in),optional :: F_msgLevel
      !@returns
      integer :: F_msgLevel0
      !@*/
      F_msgLevel0 = buffer_min_level
      if (present(F_msgLevel)) then
!$omp single
         buffer_min_level = max(MSG_NBLEVELS0,min(F_msgLevel,MSG_NBLEVELS))
!$omp end single
      endif
      return
   end function buffer_verbosity


   !/@*
   subroutine buffer_add(F_msgLevel,F_Message)
      implicit none
      !@arguments
      integer,intent(in) :: F_msgLevel
      character(len=*),intent(in) :: F_Message
      !@*/
      integer :: n
      if (F_msgLevel < buffer_min_level) return
!$omp critical
      if (nBuffer == 0 .or. .not.( &
           F_msgLevel == levelBuffer(max(1,nBuffer)) .and. &
           F_Message == msgBuffer(max(1,nBuffer))) &
           ) then
         nBuffer = nBuffer + 1
         if (nBuffer > MSG_NMAXBUFFER) then
            do n = 2,MSG_NMAXBUFFER
               levelBuffer(n-1) = levelBuffer(n)
               msgBuffer(n-1) = msgBuffer(n)
            enddo
            nBuffer = MSG_NMAXBUFFER
         endif
         levelBuffer(nBuffer) = F_msgLevel
         msgBuffer(nBuffer) = F_Message
      endif
!$omp end critical
      return
   end subroutine buffer_add


   !/@*
   subroutine buffer_reset()
      implicit none
      !@*/
      !#TODO: does single sync threads or should we use something else
!$omp single
      nBuffer = 0
!$omp end single
      return
   end subroutine buffer_reset


   !/@*
   subroutine buffer_flush()
      implicit none
      !@*/
      integer :: n, nunits, m
      integer :: units(MSG_NBLEVELS-MSG_NBLEVELS0+1)
!$omp single
      nunits = MSG_NBLEVELS - buffer_min_level + 1
      units(1:nunits) = msgUnit(buffer_min_level:MSG_NBLEVELS)
      call isort(units,nunits)
      ! unique values only
      m=1
      do n=2,nunits
         if (units(n) /= units(n-1)) m = m + 1
         units(m) = units(n)
      enddo
      nunits = m
      do n = 1, nunits
         write(units(n),*) '==== Msg Buffer traceback [BEGIN] =============================='
      enddo
      do n = nBuffer,1,-1
         write(msgUnit(levelBuffer(n)),trim(msgFormat(levelBuffer(n)))) trim(msgBuffer(n))
      enddo
      do n = 1, nunits
         write(units(n),*) '==== Msg Buffer traceback [END] ================================'
      enddo
      nBuffer = 0
!$omp end single
      return
   end subroutine buffer_flush

end module mod_msg_buffer


!/@*
subroutine msg_init()
   use mod_msg
   implicit none
   !@author  Stephane Chamberland, 2009-11
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) then
!$omp single
      isInit_L = .true.
      msgUnit(MSG_QUIET)   = MSG_STDOUT
      msgUnit(MSG_DEBUG)   = MSG_STDERR
      msgUnit(MSG_INFOPLUS)= MSG_STDOUT
      msgUnit(MSG_INFO)    = MSG_STDOUT
      msgUnit(MSG_IMPORTANT)= MSG_STDOUT
      msgUnit(MSG_WARNING) = MSG_STDOUT
      msgUnit(MSG_ERROR)   = MSG_STDERR
      msgUnit(MSG_CRITICAL)= MSG_STDERR
      msgUnit(MSG_VERBATIM)= MSG_STDOUT

      msgFormat(MSG_QUIET)   = '(a)'
      msgFormat(MSG_DEBUG)   = '("DEBUG: ",a)'
      msgFormat(MSG_INFOPLUS)= '(a)'
      msgFormat(MSG_INFO)    = '(a)'
      msgFormat(MSG_IMPORTANT)= '(a)'
      msgFormat(MSG_WARNING) = '("WARNING: ",a)'
      msgFormat(MSG_ERROR)   = '("ERROR: ",a)'
      msgFormat(MSG_CRITICAL)= '("CRITICAL ERROR: ",a)'
      msgFormat(MSG_VERBATIM)= '(a)'
!$omp end single
   endif
   !---------------------------------------------------------------------
   return
end subroutine msg_init


!#TODO: review OMP thread safety of code below

!/@*
subroutine msg_set_p0only(F_myPE)
   use mod_msg, only: isInit_L,canWrite_L
   implicit none
   !@arguments
   integer,intent(in) :: F_myPE
   !@author  Stephane Chamberland, 2009-11
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   canWrite_L = .false.
   if (F_myPe == 0) canWrite_L = .true.
   !---------------------------------------------------------------------
   return
end subroutine msg_set_p0only


!/@*
subroutine msg_set_can_write(F_canWrite_L)
   use mod_msg, only: isInit_L,canWrite_L
   implicit none
   logical,intent(in) :: F_canWrite_L
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   canWrite_L = F_canWrite_L
   !---------------------------------------------------------------------
   return
end subroutine msg_set_can_write


!/@*
subroutine msg_verbosity(F_msgLevel)
   implicit none
   !@arguments
   integer,intent(in) :: F_msgLevel
   !@*/
   !---------------------------------------------------------------------
   call msg_set_minMessageLevel(F_msgLevel)
   !---------------------------------------------------------------------
   return
end subroutine msg_verbosity


!/@*
subroutine msg_verbosity_get(F_msgLevel)
   use mod_msg
   implicit none
   !@arguments
   integer,intent(out) :: F_msgLevel
   !@*/
   !---------------------------------------------------------------------
   F_msgLevel = msgLevelMin
   !---------------------------------------------------------------------
   return
end subroutine msg_verbosity_get


!/@*
subroutine msg_set_minMessageLevel(F_msgLevel)
   use mod_msg, only: isInit_L,msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_msgLevel
   !@*/
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   msgLevelMin = max(MSG_NBLEVELS0,min(F_msgLevel,MSG_NBLEVELS))
   !---------------------------------------------------------------------
   return
end subroutine msg_set_minMessageLevel


!/@*
subroutine msg_set_redirect2fileUnit(F_out_or_err,F_fileUnit)
   use mod_msg, only: isInit_L,msgUnit
   implicit none
   !@arguments
   integer,intent(in) :: F_out_or_err
   integer,intent(in) :: F_fileUnit
   !@*/
   integer :: i, i0=1, i1=MSG_NBLEVELS
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   if (F_out_or_err == MSG_STDERR) then
      i0=MSG_ERROR
   else if (F_out_or_err == MSG_STDOUT) then
      i1=MSG_WARNING
   else
      i0 = F_out_or_err
      i1 = F_out_or_err
   endif
   do i=i0,i1
      msgUnit(i) = F_fileUnit
   enddo
   !---------------------------------------------------------------------
   return
end subroutine msg_set_redirect2fileUnit


!/@*
subroutine msg_toall(F_msgLevel,F_Message)
   use mod_msg, only: isInit_L,canWrite_L,msgUnit,msgFormat,msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_msgLevel
   character(len=*),intent(in) :: F_Message
   !@author  Stephane Chamberland, 2009-11
   !@*/
   integer :: msgLevel
   logical :: canWrite_L_bk
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   canWrite_L_bk = canWrite_L
   canWrite_L = .true.
   call msg(F_msgLevel,F_Message)
   canWrite_L = canWrite_L_bk
   !---------------------------------------------------------------------
   return
end subroutine msg_toall


!/@*
subroutine msg(F_msgLevel,F_Message)
   use mod_msg, only: isInit_L,canWrite_L,msgUnit,msgFormat,msgLevelMin
   use mod_msg_buffer, only: buffer_add
   implicit none
   !@arguments
   integer,intent(in) :: F_msgLevel
   character(len=*),intent(in) :: F_Message
   !@author  Stephane Chamberland, 2009-11
   !@*/
   integer :: msgLevel
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   call buffer_add(F_msgLevel,F_Message)
   if (canWrite_L .and. F_msgLevel>=msgLevelMin) then
      msgLevel = max(MSG_NBLEVELS0,min(F_msgLevel,MSG_NBLEVELS))
      !TODO: use C's fprintf so that it can be called from C
      write(msgUnit(msgLevel),trim(msgFormat(msgLevel))) trim(F_Message)
      call flush(msgUnit(msgLevel))
   endif
   !---------------------------------------------------------------------
   return
end subroutine msg


!/@*
function msg_getUnit(F_msgLevel) result(F_msgUnit)
   use mod_msg, only: isInit_L,canWrite_L,msgUnit,msgLevelMin
   implicit none
   !@arguments
   integer,intent(in) :: F_msgLevel
   !@returns
   integer :: F_msgUnit
   !@author  Stephane Chamberland, 2009-11
   !@*/
   integer :: msgLevel
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   F_msgUnit = -1
   if (canWrite_L .and. F_msgLevel>=msgLevelMin) then
      msgLevel = max(MSG_NBLEVELS0,min(F_msgLevel,MSG_NBLEVELS))
      F_msgUnit = msgUnit(msgLevel)
   endif
   !---------------------------------------------------------------------
   return
end function msg_getUnit


!/@*
subroutine msg_getInfo(F_canWrite_L,F_msgLevelMin,F_msgUnit,F_msgFormat_S)
   use mod_msg
   implicit none
   !@arguments
   logical,intent(out) :: F_canWrite_L
   integer,intent(out) :: F_msgLevelMin,F_msgUnit
   character(len=*),intent(out) :: F_msgFormat_S
   !@author  Stephane Chamberland, 2012-02
   !@*/
   integer,external :: msg_getUnit
   !---------------------------------------------------------------------
   if (.not.isInit_L) call msg_init()
   F_canWrite_L = canWrite_L
   F_msgLevelMin = msgLevelMin
   F_msgLevelMin = max(MSG_NBLEVELS0,min(msgLevelMin,MSG_NBLEVELS))
   F_msgUnit = msgUnit(msgLevelMin)
   F_msgFormat_S = msgFormat(msgLevelMin)
   !---------------------------------------------------------------------
   return
end subroutine msg_getInfo
 

!/@*
subroutine msg_buffer_verbosity(F_msgLevel)
   use mod_msg_buffer, only: buffer_verbosity
   implicit none
   !@objective Set minimal level of messages to be buffered
   !@arguments
   integer,intent(in) :: F_msgLevel
   !@*/
   integer :: itmp
   !---------------------------------------------------------------------
   itmp = buffer_verbosity(F_msgLevel)
   !---------------------------------------------------------------------
   return
end subroutine msg_buffer_verbosity


!/@*
subroutine msg_buffer_verbosity_get(F_msgLevel)
   use mod_msg_buffer, only: buffer_verbosity
   implicit none
   !@objective Get minimal level of messages to be buffered
   !@arguments
   integer,intent(out) :: F_msgLevel
   !@*/
   !---------------------------------------------------------------------
   F_msgLevel = buffer_verbosity()
   !---------------------------------------------------------------------
   return
end subroutine msg_buffer_verbosity_get


!/@*
subroutine msg_buffer_reset()
   use mod_msg_buffer, only: buffer_reset
   implicit none
   !@objective Reset msg buffer (empty stack)
   !@*/
   !---------------------------------------------------------------------
   call buffer_reset()
   !---------------------------------------------------------------------
   return
end subroutine msg_buffer_reset


!/@*
subroutine msg_buffer_flush()
   use mod_msg_buffer, only: buffer_flush
   implicit none
   !@objective Flush (print and reset) msg buffer
   !@*/
   !---------------------------------------------------------------------
   call buffer_flush()
   !---------------------------------------------------------------------
   return
end subroutine msg_buffer_flush
