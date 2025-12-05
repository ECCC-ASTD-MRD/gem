!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

module chm_utils_mod
! TODO faire l'entete
! This file contains general definitions that can be used in ANY modules
! All system related "magic numbers" should go here
use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64, logical_kinds
#include <rmn/msg.h>

   save

   integer,  parameter :: sp = real32
   integer,  parameter :: dp = real64
   integer,  parameter :: i4 = int32
   integer,  parameter :: i8 = int64
   integer,  parameter :: lk4 = logical_kinds(3)

   real(sp), parameter :: zero_sp = 0.0_sp
   real(sp), parameter :: tiny_sp = TINY(zero_sp)
   real(sp), parameter :: huge_sp = HUGE(zero_sp)
   real(dp), parameter :: zero_dp = 0.0_dp
   real(dp), parameter :: tiny_dp = TINY(zero_dp)
   real(dp), parameter :: huge_dp = HUGE(zero_dp)

   integer(i4), parameter      :: NOMV_LEN       = 4      ! length of variable name in a RPN std file
   integer(i4), parameter      :: LONG_VARNAME   = 16     ! Length of longer variable name
   integer(i4), parameter      :: NMLKEY_LEN     = 16     ! max length of a string namelist key
   integer(i4), parameter      :: MAX_DEBUG_VAR  = 18     ! max number of debug var for each kind (2D and 3D)
   integer(i4), parameter      :: DICTSTRING_LEN = 120    ! Length of dict string needed by chm_gestdict
   character(len=4), parameter :: UNDEFINED      = '*NA*' ! Generic undefined sequence for character variables
#if defined(MACH_TENDENCIES)
   integer(i4), parameter      :: MAX_TRACERS    = 160    ! max number of tracers to have their tendency computed
#endif
   integer(i4), parameter      :: chm_msg_debug  = MSG_DEBUG

   integer(i4)    :: chm_lun_out    ! This is the file unit used for output.
                                    ! It is assigned to STDOUT or /dev/null

   real(sp)       :: chm_timestep   ! Chemistry timestep (seconds)

   integer(i4)    :: dbg_id = -1    ! Debug horizontal grid point
   integer(i4)    :: dbg_jd = -1    ! Debug slice id

   logical(lk4)    :: chm_error_l = .false.

   contains
!============================================================================
! Name           : post_increment
!
! Description    : Return the value passed, and then increment it by 1
!
! Extra info : Equivalent to the i++ unaryoperator in C
!
! Arguments:   None
!
!============================================================================
      integer(i4) function post_increment(i)
         implicit none
         integer(i4), intent(inout) :: i
         post_increment = i
         i = i + 1
      end function post_increment

!============================================================================
! Name           : pre_increment
!
! Description    : Increment the value passed by 1 and return it
!
! Extra info : Equivalent to the ++i unary operator in C
!
! Arguments:   None
!
!============================================================================
      integer(i4) function pre_increment(i)
         implicit none
         integer(i4), intent(inout) :: i
         i = i + 1
         pre_increment = i
      end function pre_increment


      integer(i4) function ik(ix, kz, ni)
         implicit none
         integer(i4), intent(in) :: ix, kz, ni
         ik = (kz - 1) * ni + (ix - 1)
      end function ik

!============================================================================!
! Name           : chm_stop
! Description    : Premature stop in chemistry and clean exit in a MPI world.
! Extra info     : Similar to gem_stop on the GEM side.
!                  It should not use chm_lun_out, since chm_stop can be called
!                  when chm_lun_out is not set yet.
! Arguments:   IN
!                 F_name_S --> Name of the calling routine
!                 F_err    --> Error code
!===========================================================================
   subroutine chm_stop (F_name_S, F_err)
      implicit none
      character(len=*), intent(in) :: F_name_S
      integer(i4),      intent(in) :: F_err
      external                     :: handle_error

      call handle_error (F_err,F_name_S,'ABORT IN CHEMISTRY')

      return
   end subroutine chm_stop

end module chm_utils_mod
