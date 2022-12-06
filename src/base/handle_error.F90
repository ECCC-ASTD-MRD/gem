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
#include "call_back_utils.h"
#include <rmn/msg.h>
#include "stop_mpi.h"

!/@*
module mod_handle_error
   implicit none
   !@author  Stephane Chamberland, 2009-11
!*@/
   logical, public, save :: debug_L = .false.
   logical, public, save :: hasCallBackFn = .false.
   integer, public, save :: callBackFnNb
end module mod_handle_error


!/@*
subroutine collect_error(F_errorCode) !TODO: optional F_comm_S
   use iso_c_binding
   use rpn_comm_itf_mod
   use mod_handle_error, only: debug_L
   implicit none
   integer,intent(inout) :: F_errorCode
!!!#include <arch_specific.hf>
!*@/
   integer :: errcode,err
   !---------------------------------------------------------------------
   if (debug_L .and. F_errorCode < STOP_OK) then
      call msg_toall(MSG_INFOPLUS,'(collect_error) error code detected')
   endif
   call rpn_comm_allreduce(F_errorCode,errcode,1,RPN_COMM_INTEGER,"MPI_MIN",RPN_COMM_MULTIGRID,err)
   F_errorCode = errcode
   !---------------------------------------------------------------------
   return
end subroutine collect_error


!/@*
subroutine handle_error_setdebug(F_debug_L)
   use mod_handle_error, only: debug_L
   implicit none
   !@objective 
   !@arguments
   logical,intent(in) :: F_debug_L
!*@/
   !---------------------------------------------------------------------
   debug_L = F_debug_L
   !---------------------------------------------------------------------
   return
end subroutine handle_error_setdebug


!/@*
subroutine handle_error_L(F_isOK_L,F_FromSubName,F_Message)
!!$   use mod_handle_error, only: hasCallBackFn,callBackFnNb
   implicit none
!!!#include <arch_specific.hf>
   !@objective 
   !@arguments
   logical :: F_isOK_L
   character(len=*) :: F_FromSubName
   character(len=*) :: F_Message
   !@author  Michel Desgagne, Ron McTaggart-Cowan
   !@revisions
   !  2009-12,  Stephane Chamberland
!*@/
   integer :: errcode
   !---------------------------------------------------------------------
   errcode = STOP_ERROR
   if (F_isOK_L) errcode = STOP_OK
   call handle_error(errcode,F_FromSubName,F_Message)
   !---------------------------------------------------------------------
   return
end subroutine handle_error_L


!/@*
subroutine handle_error(F_errorCode,F_FromSubName,F_Message)
   use iso_c_binding
   use rpn_comm_itf_mod
!!$   use mod_handle_error, only: hasCallBackFn,callBackFnNb,debug_L
   use mod_handle_error, only: debug_L
   implicit none
!!!#include <arch_specific.hf>
#include "stop_mpi.h"
   !@objective 
   !@arguments
   integer :: F_errorCode
   character(len=*) :: F_FromSubName
   character(len=*) :: F_Message
   !@author  Michel Desgagne, Ron McTaggart-Cowan
   !@revisions
   !  2009-11,  Stephane Chamberland
!*@/
   character(len=MSG_MAXLEN) :: message
   integer :: errcode,errcode2
!!$   integer, external :: callback_call_ftn
   logical :: debug_bk_L
   !---------------------------------------------------------------------
   if (debug_L .and. F_errorCode < STOP_OK) then
      write(message,'(2a)') &
           trim(F_Message)//"\n",&
           "= ABORTING - PREMATURED STOP IN: "//trim(F_FromSubName)//"\n"
      call msg_toall(MSG_CRITICAL,message)
   endif

   debug_bk_L = debug_L
   debug_L = .false.
   errcode = F_errorCode
   call collect_error(errcode)
   debug_L = debug_bk_L

   if (errcode<STOP_OK) then
      call msg_buffer_flush()
      errcode2 = errcode
!!$      if (hasCallBackFn) then
!!$         errcode2 = callback_call_ftn(callBackFnNb,errcode)
!!$      endif
!!$      call collect_error(errcode2)
      if (errcode2<STOP_OK) then
         call stop_mpi(STOP_ERROR,F_FromSubName,F_Message)
      else
         write(message,'(3a)') 'Recovered from an error in: ',trim(F_FromSubName)
         call msg(MSG_WARNING,message)
       endif
   endif
   !---------------------------------------------------------------------
   return
end subroutine handle_error


!/@*
subroutine handle_error_set_callBackFn(F_callBackFn)
   use mod_handle_error, only: hasCallBackFn,callBackFnNb
   implicit none
!!!#include <arch_specific.hf>
   !@objective Set a fn to be called before stop in case of error
   !@arguments
   integer,external :: F_callBackFn
   !@author  Stephane Chamberland, Nov 2009
   !@description
   !  the call back function should have the following interface 
   !  (with a differnt name!):
   !
   !  function myCallBack(F_errorCode) result(newErrorCode)
   !     integer :: F_errorCode
   !     integer :: newErrorCode
   !
   !  It should return the control to the caller 
   !  so the program can end gracefully
!*@/
   integer,external :: callback_register_ftn
   !---------------------------------------------------------------------
!   callBackFnNb = callback_register_ftn(F_callBackFn)
   if (callBackFnNb == UTILS_CALL_BACK_OK) &
        hasCallBackFn = .true.
   !---------------------------------------------------------------------
   return
end subroutine handle_error_set_callBackFn
