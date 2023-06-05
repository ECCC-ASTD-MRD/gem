!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

subroutine timing_init2(myproc, msg)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_getenv, CLIB_IS_OK
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   character(len=*), intent(in) :: msg
   integer, intent(in) :: myproc
   !@author M. Desgagne   -- Winter 2012 --
   !@revision
   ! v4_40 - Desgagne - initial version
   ! v4_80 - Desgagne - introduce timer_level and timer_cnt

#include <rmnlib_basics.hf>
   include "timing.cdk"

   character(len=16) :: dumc_S
   !-------------------------------------------------------------------
   if (.not.CLIB_IS_OK(clib_getenv('TMG_ON', Timing_S))) Timing_S = 'NO'
   call low2up(Timing_S, dumc_S)
   Timing_S = dumc_S

   if (Timing_S == 'YES') call tmg_init(myproc, msg)

   sum_tb = 0.D0; timer_cnt = 0 ; timer_level = 0 ; nam_subr_S = ''
   !-------------------------------------------------------------------
   return
end subroutine timing_init2
