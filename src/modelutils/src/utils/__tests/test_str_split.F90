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
subroutine test_str_split()
   use, intrinsic :: iso_fortran_env, only: INT64
   use str_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-01
!@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   integer,parameter :: NMAX = 32
   integer,parameter :: KEY = 1
   integer,parameter :: VAL = 2

   character(len=512) :: string_S,sep_S,part1_S,part2_S,parts_S(NMAX),kv_S(2,NMAX)
   integer :: nkeys,nn
   logical :: ok_L

   integer,external :: str_split2keyval,str_split2keyval0
   ! ---------------------------------------------------------------------
   call msg_set_minMessageLevel(MSG_DEBUG)
   !call msg_set_minMessageLevel(MSG_ERROR)

   string_S = ' "A quoted string" '
   call str_rm_quotes(string_S)
   ok_L = (string_S == "A quoted string")
   call testutils_assert_ok(ok_L,'test_str_split:str_rm_quotes','double')

   string_S = "'I'm a quoted string'"
   call str_rm_quotes(string_S)
   ok_L = (string_S == "I'm a quoted string")
   call testutils_assert_ok(ok_L,'test_str_split:str_rm_quotes','single')

   string_S = "\ta\ttabbed\tstring with some more"
   call str_tab2space(string_S)
   ok_L = (string_S == " a tabbed string with some more")
   call testutils_assert_ok(ok_L,'test_str_split:str_tab2space','')

   string_S = "CVS,is,a,comma,separated,string"
   call str_split(part1_S,part2_S,string_S,',')
   ok_L = (part1_S=='CVS' .and. part2_S=='is,a,comma,separated,string')
   call testutils_assert_ok(ok_L,'test_str_split:str_split','')

   string_S = "CVS,is,a,comma,separated,string"
   parts_S(:) = ' '
   call str_split2list(parts_S,string_S,',',NMAX)
   ok_L = (parts_S(1)=='CVS' .and. parts_S(2)=='is' .and.parts_S(3)=='a' .and. &
        parts_S(4)=='comma' .and. parts_S(5)=='separated' .and. &
        parts_S(6)=='string' .and. all(parts_S(7:NMAX)==' '))
   call testutils_assert_ok(ok_L,'test_str_split:str_split2list','')
   if (.not.ok_L) then
      do nn=1,7
         print *,'parts_S(',nn,') = ',trim(parts_S(nn))
      enddo
   endif

   string_S = 'VN=V0; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
   kv_S(:,:) = ' '
   nkeys = str_split2keyval0(kv_S,string_S,NMAX)
   ok_L = (nkeys == 4 .and. &
        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == 'V0' .and. &
        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == 'ANAL' .and. &
        kv_S(KEY,3) == 'interp' .and. kv_S(VAL,3) == 'NEAREST' .and. &
        kv_S(KEY,4) == 'FREQ' .and. kv_S(VAL,4) == '0' .and. &
        all(kv_S(:,5:NMAX) == ' ') )
   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval0','')
   if (.not.ok_L) then
      print *,'nkeys=',nkeys
      do nn=1,nkeys+1
         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
      enddo
   endif

   string_S = 'VN=V0; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
   kv_S(:,:) = ' '
   nkeys = str_split2keyval(kv_S,string_S,NMAX)
   ok_L = (nkeys == 4 .and. &
        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == 'V0' .and. &
        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == 'ANAL' .and. &
        kv_S(KEY,3) == 'interp' .and. kv_S(VAL,3) == 'NEAREST' .and. &
        kv_S(KEY,4) == 'FREQ' .and. kv_S(VAL,4) == '0' .and. &
        all(kv_S(:,5:NMAX) == ' ') )
   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval','')
   if (.not.ok_L) then
      print *,'nkeys=',nkeys
      do nn=1,nkeys+1
         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
      enddo
   endif

   string_S = 'VN=V0; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
   kv_S(:,:) = ' '
   kv_S(KEY,1) = 'FREQ'
   kv_S(KEY,2) = 'VN'
   kv_S(KEY,5) = 'TOTO' ; kv_S(VAL,5) = 'tata'
   nkeys = str_split2keyval(kv_S,string_S,NMAX)
   ok_L = (nkeys == 4 .and. &
        kv_S(KEY,1) == 'FREQ' .and. kv_S(VAL,1) == '0' .and. &
        kv_S(KEY,2) == 'VN' .and. kv_S(VAL,2) == 'V0' .and. &
        kv_S(KEY,3) == 'SEARCH' .and. kv_S(VAL,3) == 'ANAL' .and. &
        kv_S(KEY,4) == 'interp' .and. kv_S(VAL,4) == 'NEAREST' .and. &
        kv_S(KEY,5) == 'TOTO' .and. kv_S(VAL,5) == 'tata' .and. &
        all(kv_S(:,6:NMAX) == ' ') )
   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval',' with predef keys')
   if (.not.ok_L) then
      print *,'nkeys=',nkeys
      do nn=1,nkeys+1
         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
      enddo
   endif

   kv_S(:,:) = ' '
   string_S = 'VN="test = test"; SEARCH="my is"'
   nkeys = str_split2keyval(kv_S,string_S,NMAX)
   ok_L = (nkeys == 2 .and. &
        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == 'test = test' .and. &
        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == 'my is' .and. &
        all(kv_S(:,3:NMAX) == ' ') )
   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval',' with with = in value')
   if (.not.ok_L) then
      print *,'nkeys=',nkeys
      do nn=1,nkeys+1
         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
      enddo
   endif

!!$   kv_S(:,:) = ' '
!!$   string_S = 'VN="test"; SEARCH="my ; is"'
!!$   nkeys = str_split2keyval(kv_S,string_S,NMAX)
!!$   ok_L = (nkeys == 2 .and. &
!!$        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == 'test' .and. &
!!$        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == 'my ; is' .and. &
!!$        all(kv_S(:,3:NMAX) == ' ') )
!!$   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval',' with with ; in value')
!!$   if (.not.ok_L) then
!!$      print *,'nkeys=',nkeys
!!$      do nn=1,nkeys+1
!!$         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
!!$      enddo
!!$   endif

!!$   kv_S(:,:) = ' '
!!$   string_S = 'VN=" test"; SEARCH="my is"'
!!$   nkeys = str_split2keyval(kv_S,string_S,NMAX)
!!$   ok_L = (nkeys == 2 .and. &
!!$        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == ' test' .and. &
!!$        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == '"my is"' .and. &
!!$        all(kv_S(:,3:NMAX) == ' ') )
!!$   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval0',' quoted lead blank preseved')
!!$   if (.not.ok_L) then
!!$      print *,'nkeys=',nkeys
!!$      do nn=1,nkeys+1
!!$         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
!!$      enddo
!!$   endif


   kv_S(:,:) = ' '
   string_S = 'VN="test = test"; SEARCH="my is"'
   nkeys = str_split2keyval0(kv_S,string_S,NMAX)
   ok_L = (nkeys == 2 .and. &
        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == '"test = test"' .and. &
        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == '"my is"' .and. &
        all(kv_S(:,3:NMAX) == ' ') )
   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval0',' with with = in value')
   if (.not.ok_L) then
      print *,'nkeys=',nkeys
      do nn=1,nkeys+1
         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
      enddo
   endif

!!$   kv_S(:,:) = ' '
!!$   string_S = 'VN="test"; SEARCH="my ; is"'
!!$   nkeys = str_split2keyval0(kv_S,string_S,NMAX)
!!$   ok_L = (nkeys == 2 .and. &
!!$        kv_S(KEY,1) == 'VN' .and. kv_S(VAL,1) == 'test' .and. &
!!$        kv_S(KEY,2) == 'SEARCH' .and. kv_S(VAL,2) == 'my ; is' .and. &
!!$        all(kv_S(:,3:NMAX) == ' ') )
!!$   call testutils_assert_ok(ok_L,'test_str_split:str_split2keyval0',' with with ; in value')
!!$   if (.not.ok_L) then
!!$      print *,'nkeys=',nkeys
!!$      do nn=1,nkeys+1
!!$         print *,nn,trim(kv_S(1,nn)),' : ',trim(kv_S(2,nn))
!!$      enddo
!!$   endif


   ! ---------------------------------------------------------------------
   return
end subroutine test_str_split
