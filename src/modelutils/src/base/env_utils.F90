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
module env_utils
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_getenv
   use str_mod
   implicit none
   private
   public :: env_get
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   interface env_get
      module procedure env_get_str
      module procedure env_get_int
      module procedure env_get_real
      module procedure env_get_bool
   end interface

   !TODO: env_put

contains

   !/@*
   function env_get_str(F_name_S,F_sval,F_default,F_normalize_L,F_okvalues) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_name_S
      character(len=*),intent(out) :: F_sval
      character(len=*),intent(in),optional :: F_default
      logical,intent(in),optional :: F_normalize_L
      character(len=*),intent(in),optional :: F_okvalues(:)
      !@returns
      integer :: F_istat
      !@author  Stephane Chamberland, 2020-10
      !*@/
      !---------------------------------------------------------------------
      F_istat = clib_getenv(trim(F_name_S), F_sval)
      if (.not.RMN_IS_OK(F_istat) .and. present(F_default)) then
         F_sval = F_default
         F_istat = RMN_OK
      endif
      if (.not.RMN_IS_OK(F_istat)) return
      if (present(F_normalize_L)) then
         if (F_normalize_L) call str_normalize(F_sval)
      endif
      if (present(F_okvalues)) then
         if (.not.any(F_sval == F_okvalues)) then
            F_sval = ''
            F_istat = RMN_ERR
         endif
      endif
      !---------------------------------------------------------------------
      return
   end function env_get_str
   
   !/@*
   function env_get_int(F_name_S,F_ival,F_default,F_min,F_max) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_name_S
      integer,intent(out) :: F_ival
      integer,intent(in),optional :: F_default,F_min,F_max
      !@returns
      integer :: F_istat
      !@author  Stephane Chamberland, 2012-02
      !*@/
      character(len=1024) :: tmp_S
      !---------------------------------------------------------------------
      F_istat = clib_getenv(trim(F_name_S),tmp_S)
      F_istat = min(str_toint(F_ival,tmp_S),F_istat)
      if (.not.RMN_IS_OK(F_istat) .and. present(F_default)) then
         F_ival = F_default
         F_istat = RMN_OK
      endif
      if (.not.RMN_IS_OK(F_istat)) return
      if (present(F_min)) F_ival = max(F_min,F_ival)
      if (present(F_max)) F_ival = min(F_ival,F_max)
      !---------------------------------------------------------------------
      return
   end function env_get_int


   !/@*
   function env_get_real(F_name_S,F_rval,F_default,F_min,F_max) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_name_S
      real,intent(out) :: F_rval
      real,intent(in),optional :: F_default,F_min,F_max
      !@returns
      integer :: F_istat
      !@author  Stephane Chamberland, 2012-02
      !*@/
      character(len=1024) :: tmp_S
      !---------------------------------------------------------------------
      F_istat = clib_getenv(trim(F_name_S),tmp_S)
      F_istat = min(str_toreal(F_rval,tmp_S),F_istat)
      if (.not.RMN_IS_OK(F_istat) .and. present(F_default)) then
         F_rval = F_default
         F_istat = RMN_OK
      endif
      if (.not.RMN_IS_OK(F_istat)) return
      if (present(F_min)) F_rval = max(F_min,F_rval)
      if (present(F_max)) F_rval = min(F_rval,F_max)
      !---------------------------------------------------------------------
      return
   end function env_get_real


   !/@*
   function env_get_bool(F_name_S,F_lval,F_default) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_name_S
      logical,intent(out) :: F_lval
      logical,intent(in),optional :: F_default
      !@returns
      integer :: F_istat
      !@author  Stephane Chamberland, 2012-02
      !*@/
      character(len=1024) :: tmp_S
      !---------------------------------------------------------------------
      F_istat = clib_getenv(trim(F_name_S),tmp_S)
      F_istat = min(str_tobool(F_lval,tmp_S),F_istat)
      if (.not.RMN_IS_OK(F_istat) .and. present(F_default)) then
         F_lval = F_default
         F_istat = RMN_OK
      endif
      !---------------------------------------------------------------------
      return
   end function env_get_bool

end module env_utils

