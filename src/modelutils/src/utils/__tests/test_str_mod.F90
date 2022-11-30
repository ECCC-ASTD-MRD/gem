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

#include <rmn/msg.h>

!/@
subroutine test_str_mod()
   use, intrinsic :: iso_fortran_env, only: INT64
   use str_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-01
!@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>

   character(len=512) :: string_S
   logical :: ok_L, tmp_L
   integer :: tmp_i,istat
   real :: tmp_r
   ! ---------------------------------------------------------------------
   call msg_set_minMessageLevel(MSG_DEBUG)
   !call msg_set_minMessageLevel(MSG_ERROR)

   string_S = ' "A quoted string" '
   call str_rm_quotes(string_S)
   ok_L = (string_S == "A quoted string")
   call testutils_assert_ok(ok_L,'test_str_mod:str_rm_quotes','double')

   string_S = "'I'm a quoted string'"
   call str_rm_quotes(string_S)
   ok_L = (string_S == "I'm a quoted string")
   call testutils_assert_ok(ok_L,'test_str_mod:str_rm_quotes','single')


   string_S = "\ta\ttabbed\tstring with some more"
   call str_tab2space(string_S)
   ok_L = (string_S == " a tabbed string with some more")
   call testutils_assert_ok(ok_L,'test_str_mod:str_tab2space','')


   string_S = "more"
   istat = str_toint(tmp_i,string_S)
   ok_L = (.not.RMN_IS_OK(istat))
   call testutils_assert_ok(ok_L,'test_str_mod:str_toint','bad value str')

!!$   string_S = "1.2"
!!$   istat = str_toint(tmp_i,string_S)
!!$   ok_L = (.not.RMN_IS_OK(istat))
!!$   call testutils_assert_ok(ok_L,'test_str_mod:str_toint','bad value real')
!!$   if (.not.ok_L) print *,istat,'test_str_mod:str_toint got:',tmp_i,' for ',trim(string_S)
!!$
!!$   string_S = "2.7"
!!$   istat = str_toint(tmp_i,string_S)
!!$   ok_L = (.not.RMN_IS_OK(istat))
!!$   call testutils_assert_ok(ok_L,'test_str_mod:str_toint','bad value real 2')
!!$   if (.not.ok_L) print *,istat,'test_str_mod:str_toint got:',tmp_i,' for ',trim(string_S)

   string_S = "223"
   istat = str_toint(tmp_i,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_i==223)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toint','ok value 1')

   string_S = "-99"
   istat = str_toint(tmp_i,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_i==-99)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toint','ok value 2')



   string_S = "bad value"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (.not.RMN_IS_OK(istat))
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','bad value')

   string_S = "223"
   istat =str_toreal (tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==223.)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 1')

   string_S = "-99"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==-99.)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 2')

   string_S = "123.0"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==123.)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 3')

   string_S = ".5"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==.5)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 4')

   string_S = "5.e2"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==5.e2)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 5')

   string_S = "7.e-3"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==7.e-3)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 6')

   string_S = "0.678E-6"
   istat = str_toreal(tmp_r,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_r==0.678E-6)
   call testutils_assert_ok(ok_L,'test_str_mod:str_toreal','ok value 6')


   string_S = "bad value"
   istat = str_tobool(tmp_L,string_S)
   ok_L = (.not.RMN_IS_OK(istat))
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','bad value')

   string_S = "true"
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 1')

   string_S = ".t."
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 2')

   string_S = ".T."
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 3')

   string_S = ".true."
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 4')

   string_S = "False"
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. .not.tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 5')

   string_S = ".f."
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. .not.tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 6')

   string_S = ".F."
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. .not.tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 7')

   string_S = ".false."
   istat = str_tobool(tmp_L,string_S)
   ok_L = (RMN_IS_OK(istat) .and. .not.tmp_L)
   call testutils_assert_ok(ok_L,'test_str_mod:str_tobool','ok value 8')
   ! ---------------------------------------------------------------------
   return
end subroutine test_str_mod
