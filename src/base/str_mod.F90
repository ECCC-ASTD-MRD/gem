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


!/@
module str_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   implicit none
   private
   !@objective Provide string manipulation functions
   !@author Stephane Chamberland, 2011-09
   !@description
   ! Public functions
   public :: str_normalize,str_rm_quotes,str_tab2space,str_toint,str_toreal,str_tobool, str_concat, str_concat_i, str_encode_num
   ! Public constants
   !
!@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

contains

   !TODO: include str_split fn
 
   !/@*
   subroutine str_normalize(F_str)
      implicit none
      !@objective 
      !@arguments
      character(len=*),intent(inout) :: F_str
      !@author  S. Chamberland, 2012-05
      !*@/
      integer :: istat
      !---------------------------------------------------------------------
      call str_tab2space(F_str)
      F_str = adjustl(F_str)
      istat = clib_tolower(F_str)
      !---------------------------------------------------------------------
      return
   end subroutine str_normalize


   !/@*
   subroutine str_rm_quotes(F_str)
      implicit none
      !@objective Remove leading and trailing quotes/blanks if any
      !@arguments
      character(len=*),intent(inout) :: F_str
      !@author  S. Chamberland, 2010-03
      !*@/
      integer :: ii
      !---------------------------------------------------------------------
      F_str = adjustl(F_str)
      if (F_str(1:1) == "'" .or. F_str(1:1) == '"') then
         ii = len_trim(F_str)
         if (F_str(ii:ii) == F_str(1:1)) then
            F_str(1:1)   = ' '
            F_str(ii:ii) = ' '
            F_str = adjustl(F_str)
         endif
      endif
      !---------------------------------------------------------------------
      return
   end subroutine str_rm_quotes


   !/@*
   subroutine str_tab2space(F_str) !,F_nspaces
      implicit none
      !@objective Replace tab char by space char (so that trim/adjustl work!)
      !@arguments
      character(len=*),intent(inout) :: F_str
      !integer,intent(in) :: F_nspaces !replace tab chars by F_nspaces spaces (not yet implemented - default to 1)
      !@author  S. Chamberland, 2011-04
      !*@/
      integer,parameter :: ASCII_TAB = 9
      integer :: ii
      !---------------------------------------------------------------------
      do ii=1,len_trim(F_str)
         if (iachar(F_str(ii:ii)) == ASCII_TAB) F_str(ii:ii) = ' '
      enddo
      !---------------------------------------------------------------------
      return
   end subroutine str_tab2space


   !/@*
   function str_toint(F_int,F_str_S) result(F_istat)
      implicit none
      integer,intent(out) :: F_int
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer :: istat
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      F_int = 0
      read(F_str_S,*,iostat=istat) F_int
      !TODO: prevent from casting a real
      if (istat == 0) F_istat = RMN_OK
      !---------------------------------------------------------------------
      return
   end function str_toint


   !/@*
   function str_toreal(F_real,F_str_S) result(F_istat)
      implicit none
      real,intent(out) :: F_real
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      integer :: istat
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      F_real = 0.
      read(F_str_S,*,iostat=istat) F_real
      if (istat == 0) F_istat = RMN_OK
      !---------------------------------------------------------------------
      return
   end function str_toreal


   !/@*
   function str_tobool(F_bool_L,F_str_S) result(F_istat)
      implicit none
      logical,intent(out) :: F_bool_L
      character(len=*),intent(in) :: F_str_S
      integer :: F_istat
      !*@/
      character(len=32) :: tmp_S
      integer :: istat
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      F_bool_L = .false.
      tmp_S = F_str_S
      call str_tab2space(tmp_S)
      tmp_S = adjustl(tmp_S)
      istat = clib_tolower(tmp_S)
      if (tmp_S(1:1) == 't' .or. tmp_S(1:2) == '.t') then
         F_bool_L = .true.
         F_istat = RMN_OK
      else if (tmp_S(1:1) == 'f' .or. tmp_S(1:2) == '.f') then
         F_istat = RMN_OK
      endif
      !---------------------------------------------------------------------
      return
   end function str_tobool

 
   !/@*
   subroutine str_concat(F_str_out_S,F_str_array_S,F_sep_S)
      implicit none
      character(len=*),intent(out) :: F_str_out_S
      character(len=*),intent(in)  :: F_str_array_S(:),F_sep_S
      !*@/
      integer :: i
      character(len=512) :: str
      !--------------------------------------------------------------------
      str = F_str_array_S(1)
      call str_tab2space(str)
      F_str_out_S = adjustl(str)
      do i = 2, size(F_str_array_S)
         str = F_str_array_S(i)
         call str_tab2space(str)
         F_str_out_S = trim(F_str_out_S)//F_sep_S//trim(adjustl(str))
      enddo
      !---------------------------------------------------------------------
      return
   end subroutine str_concat

   
   !/@*
   subroutine str_concat_i(F_str_out_S,F_i_array,F_sep_S)
      implicit none
      character(len=*),intent(out) :: F_str_out_S
      integer,         intent(in)  :: F_i_array(:)
      character(len=*),intent(in)  :: F_sep_S
      !*@/
      integer :: i
      character(len=512) :: str
      !--------------------------------------------------------------------
      write(str, '(i0)') F_i_array(1)
      call str_tab2space(str)
      F_str_out_S = adjustl(str)
      do i = 2, size(F_i_array)
         write(str, '(i0)') F_i_array(i)
         call str_tab2space(str)
         F_str_out_S = trim(F_str_out_S)//F_sep_S//trim(adjustl(str))
      enddo
      !---------------------------------------------------------------------
      return
   end subroutine str_concat_i
   

   !/@*
   function str_encode_num(F_val,F_codes_S) result(F_str_S)
      implicit none
      integer, intent(in) :: F_val
      character(len=*),intent(in), optional  :: F_codes_S
      character(len=512) :: F_str_S
      !*@/
      character(len=*), parameter :: DEFAULT_S = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer :: ibase, iexp, iexpval, idiv, val, x
      character(len=512) :: codes_S
      !--------------------------------------------------------------------
      F_str_S = ''
      codes_s = DEFAULT_S
      if (present(F_codes_S)) codes_S = F_codes_S
      ibase = len_trim(codes_S)
      if (ibase == 0) return
      if (F_val < 0) return
      if (F_val == 0) then
         F_str_S = codes_S(1:1)
         return
      endif
      F_str_S = ''
      iexp = int(log(real(F_val)) / log(real(ibase)))
      val = F_val
      do x = iexp,0,-1
         iexpval = ibase**x
         idiv = val / iexpval
         F_str_S = trim(F_str_S)//codes_S(idiv+1:idiv+1)
         val = val - (idiv*iexpval)
      enddo
      !---------------------------------------------------------------------
      return
   end function str_encode_num

end module str_mod
