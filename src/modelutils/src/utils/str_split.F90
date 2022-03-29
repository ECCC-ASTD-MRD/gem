!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

!/@*
subroutine str_split(F_part1,F_part2,F_str,F_sep)
   implicit none
!!!#include <arch_specific.hf>
   !@objective Split a string at first occurence of a separator
   !@arguments
   character(len=*),intent(out) :: F_part1 !part of str before 1st occurence of sep
   character(len=*),intent(out) :: F_part2 !part of str after  1st occurence of sep
   character(len=*),intent(in) :: F_str    !string to split
   character(len=*),intent(in) :: F_sep    !separator
   !@author  S. Chamberland, 2010-03
   !*@/
   integer :: ii
   !---------------------------------------------------------------------
   ii = index(F_str,F_sep)
   if (ii > 0) then
     F_part1  = trim(adjustl(F_str(1:ii-1)))
     F_part2  = trim(adjustl(F_str(ii+len(F_sep):len(F_str))))
   else
     F_part1  = trim(adjustl(F_str))
     F_part2  = ''
   endif
   !---------------------------------------------------------------------
   return
end subroutine str_split


!/@*
subroutine str_split2list(F_parts,F_str,F_sep,F_nmax)
   implicit none
!!!#include <arch_specific.hf>
   !@objective Split a string at first occurence of a separator
   !@arguments
   integer,intent(in) :: F_nmax
   character(len=*),intent(in) :: F_str    !string to split
   character(len=*),intent(in) :: F_sep    !separator
   character(len=*),intent(out) :: F_parts(F_nmax)
   !@author  S. Chamberland, 2011-04
   !*@/
   integer :: ii
   character(len=1024) :: s1_S,s2_S,s0_S
   !---------------------------------------------------------------------
   if (F_nmax < 1) return
   F_parts(:) = ' '
   s0_S = F_str
   ii = 1
   do
      call str_split(s1_S,s2_S,s0_S,F_sep)
      F_parts(ii) = s1_S
      if (s2_S == ' ') exit
      if (ii == F_nmax) then
         F_parts(ii) = s0_S
         exit
      endif
      s0_S = s2_S
      ii = ii + 1
   enddo
   !---------------------------------------------------------------------
   return
end subroutine str_split2list

!TODO: inverse of str_split2keyval (str_keyval2string)

!/@*
function str_split2keyval(F_kv_S,F_string_S,F_nmax) result(F_nkeys)
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   use str_mod, only: str_tab2space, str_rm_quotes
   implicit none
   !@objective split a config like string (key1=val;key2=val) in a set of key/val
   !@arguments
   integer,intent(in) :: F_nmax
   character(len=*),intent(in) :: F_string_S
   character(len=*),intent(inout) :: F_kv_S(2,F_nmax)
   !@return
   integer :: F_nkeys
   !@description
   !  F_kv_S(1,:) has keys
   !  F_kv_S(2,:) has correponding values
   !
   !  this function has 2 modes:
   !  1) F_kv_S == ' '       : all found keys are parsed and kept
   !  2) F_kv_S(1,:) != ' '  : desired keynames are provided, possible blank ones left at the end will be filled with unknown keys
   !
   ! F_string format: key1=value1 ; key2="value 2" ; key3 = "value's"
   !
   ! F_string_S's keys should not contain space, tab, = or ';'
   ! F_string_S's values should not contain ';'
   ! spaces around '=', ';' are ignored
   ! tabs are converted to one space
   !*@/
   integer,parameter :: KEY = 1
   integer,parameter :: VAL = 2
   character(len=1024) :: s0_S,key_S,key2_S,key3_S,val_S,list_S(F_nmax)
   integer :: n,n2,idx,istat
   !------------------------------------------------------------------
   s0_S = F_string_S
   call str_tab2space(s0_S)
   s0_S = adjustl(s0_S)
 
   F_nkeys = 0
   do n=1,F_nmax
      if (F_kv_S(KEY,n) == '') exit
      F_nkeys = F_nkeys + 1
   enddo

   list_S = ' '
   call str_split2list(list_S,s0_S,';',F_nmax)
   LOOP1: do n=1,F_nmax
      if (list_S(n) == ' ') cycle
      call str_split(key_S,val_S,list_S(n),'=')
      call str_rm_quotes(key_S)
      if (key_S == ' ') cycle LOOP1
      call str_rm_quotes(val_S)
      key2_S = key_S
      istat = clib_tolower(key2_S)
      idx = F_nkeys + 1
      LOOP2: do n2=1,F_nkeys
         key3_S = F_kv_S(KEY,n2)
         istat = clib_tolower(key3_S)
         if (key2_S == key3_S) then
            idx = n2
            exit LOOP2
         endif
      enddo LOOP2
      if (idx > F_nmax) cycle
      if (idx > F_nkeys) F_nkeys = idx
      if (F_kv_S(KEY,idx) == ' ') F_kv_S(KEY,idx) = key_S
      F_kv_S(VAL,idx) = val_S
   enddo LOOP1
   !------------------------------------------------------------------
   return
end function str_split2keyval


!/@*
function str_split2keyval0(F_kv_S,F_string_S,F_nmax) result(F_nkeys)
   use str_mod, only: str_tab2space
   implicit none
   !@objective split a config like sting (key1=val;key2=val) in a set of key/val
   !@arguments
   integer,intent(in) :: F_nmax
   character(len=*),intent(in) :: F_string_S
   character(len=*),intent(inout) :: F_kv_S(2,F_nmax)
   !@description
   !  F_kv_S(1,:) has keys
   !  F_kv_S(2,:) has correponding values
   !@return
   integer,parameter :: KEY = 1
   integer,parameter :: VAL = 2
   integer :: F_nkeys
   !*@/
   character(len=1024) :: s1_S,s2_S,s0_S,s0b_S
   !------------------------------------------------------------------
   F_kv_S(:,:) = ' '
   s0_S = F_string_S
   call str_tab2space(s0_S)
   F_nkeys = 0
   DOITEMS: do
      call str_split(s1_S,s2_S,s0_S,';')
      s0_S = s2_S
      if (s1_S /= ' ') then
         s0b_S = s1_S
         call str_split(s1_S,s2_S,s0b_S,'=')
         if (s1_S /= ' ' .and. s2_S /= ' ') then
            F_nkeys = F_nkeys + 1
            if (F_nkeys > F_nmax) return
            F_kv_S(KEY,F_nkeys) = adjustl(s1_S)
            F_kv_S(VAL,F_nkeys) = adjustl(s2_S)
         endif
      endif
      if (s0_S == ' ') exit DOITEMS
   enddo DOITEMS
   !------------------------------------------------------------------
   return
end function str_split2keyval0
