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
!/@*
subroutine stop_mpi(F_EndType,F_FromSubName,F_Message)
   use iso_c_binding
   use rpn_comm_itf_mod
   implicit none
!!!#include <arch_specific.hf>
#include "stop_mpi.h"
#include <rmn/msg.h>

   !@objective 
   !@arguments
   integer :: F_EndType !STOP_OK, STOP_ERROR
   character(len=*) :: F_FromSubName
   character(len=*) :: F_Message
   !@author  Michel Desgagne - July 2001
   !@revisions
   !  2009-11,  Stephane Chamberland
!*@/
   integer :: err,msgLevel
   character(len=MSG_MAXLEN) :: message,msg_part_S
   character(len=1024) :: msg_format = '(8a)'
   character(len=1024) :: msg_line_S = '==============================================\n'
   character(len=1024) :: msg_ok_S = '= STOPPING - NORMAL ENDING OF: '
   character(len=1024) :: msg_err_S = '= ABORTING - PREMATURED STOP IN: '
   !---------------------------------------------------------------------
   msgLevel = MSG_CRITICAL
   msg_part_S = msg_err_S
   if (min(STOP_OK,F_EndType) == STOP_OK) then
      msgLevel = MSG_INFO
      msg_part_S = msg_ok_S
   endif
   write(message,msg_format) &
        trim(F_Message),"\n",&
        trim(msg_line_S), &
        trim(msg_part_S), " ",trim(F_FromSubName), "\n", &
        trim(msg_line_S)
   call msg(msgLevel,message)

   call rpn_comm_Barrier(RPN_COMM_GRID, err)
   call rpn_comm_FINALIZE(err)
   !---------------------------------------------------------------------
   stop
end subroutine stop_mpi

