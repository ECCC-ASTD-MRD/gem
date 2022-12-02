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
!-------------------------------------- LICENCE END --------------------------

!/@*
subroutine tracers_attributes2(F_attributes_S, F_wload, F_hzd, &
     F_monot, F_massc, F_min)
   implicit none
!!!#include <arch_specific.hf>
   !@objective
   !@author
   !@arguments
   character(len=*), intent(in) :: F_attributes_S
   integer, intent(out) :: F_wload, F_hzd, F_monot, F_massc
   real, intent(out) :: F_min
   !*@/
   character(len=2048) :: attributes
   logical :: change_default_L
   integer :: ind,deb,eqi
   integer, save :: default_wload = 0
   integer, save :: default_hzd   = 0
   integer, save :: default_monot = 1
   integer, save :: default_massc = 0
   real, save :: default_vmin = -1*huge(1.)
   !---------------------------------------------------------------
   F_wload = default_wload
   F_hzd   = default_hzd
   F_monot = default_monot
   F_massc = default_massc
   F_min   = default_vmin

   if (trim(F_attributes_S) == '') return
   ind = 0
   call low2up(F_attributes_S, attributes)

   change_default_L = .false.
44 if (ind > len(attributes)) return

   deb = ind+1
   ind = index(attributes(deb:), ",")
   if (ind == 0) then
      ind = len(attributes) + 1
   else
      ind = ind + deb - 1
   endif

   if (trim(attributes(deb:ind-1)) == 'DEFAULT') then
      change_default_L = .true.
      goto 44
   endif

   eqi = index(attributes(deb:ind-1), "=")

   if (change_default_L) then

      if (eqi > 0) then
         eqi = eqi + deb - 1

         if (trim(attributes(deb:eqi-1)) == 'WLOAD') then
            read(attributes(eqi+1:ind-1),*) default_wload
            F_wload = default_wload
         endif
         if (trim(attributes(deb:eqi-1)) == 'HZD') then
            read(attributes(eqi+1:ind-1),*) default_hzd
            F_hzd = default_hzd
         endif
         if (trim(attributes(deb:eqi-1)) == 'MONO') then
            read(attributes(eqi+1:ind-1),*) default_monot
            F_monot = default_monot
         endif
         if (trim(attributes(deb:eqi-1)) == 'MASS') then
            read(attributes(eqi+1:ind-1),*) default_massc
            F_massc = default_massc
         endif
         if (trim(attributes(deb:eqi-1)) == 'MIN') then
            read(attributes(eqi+1:ind-1),*) default_vmin
            F_min = default_vmin
         endif
      endif

   else

      if (eqi == 0) then
         if (trim(attributes(deb:ind-1)) == 'WLOAD') &
              F_wload = default_wload
         if (trim(attributes(deb:ind-1)) == 'HZD') &
              F_hzd = default_hzd
         if (trim(attributes(deb:ind-1)) == 'MONO') &
              F_monot = default_monot
         if (trim(attributes(deb:ind-1)) == 'MASS') &
              F_massc = default_massc
         if (trim(attributes(deb:ind-1)) == 'MIN') &
              F_min = default_vmin
      else
         eqi = eqi + deb - 1

         if (trim(attributes(deb:eqi-1)) == 'WLOAD') &
              read(attributes(eqi+1:ind-1),*) F_wload
         if (trim(attributes(deb:eqi-1)) == 'HZD') &
              read(attributes(eqi+1:ind-1),*) F_hzd
         if (trim(attributes(deb:eqi-1)) == 'MONO') &
              read(attributes(eqi+1:ind-1),*) F_monot
         if (trim(attributes(deb:eqi-1)) == 'MASS') &
              read(attributes(eqi+1:ind-1),*) F_massc
         if (trim(attributes(deb:eqi-1)) == 'MIN') &
              read(attributes(eqi+1:ind-1),*) F_min
      endif

   endif

   goto 44
   !---------------------------------------------------------------
   return
end subroutine tracers_attributes2


module tracers_attributes_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use str_mod, only: str_tab2space, str_toint, str_toreal
   private
   public :: tracers_attributes

contains

   !/@*
   function tracers_attributes(F_attributes_S, F_wload, F_hzd, &
        F_monot, F_massc, F_vmin, F_vmax, F_intp_S, F_ignore_L) result(F_istat)
      use clib_itf_mod, only: clib_toupper
      implicit none
#include <rmn/msg.h>
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@objective
      !@author
      !@arguments
      character(len=*), intent(in) :: F_attributes_S
      integer, intent(out), optional :: F_wload, F_hzd, F_monot, F_massc
      real, intent(out), optional :: F_vmin, F_vmax
      character(len=*), intent(out), optional :: F_intp_S
      logical, intent(in), optional :: F_ignore_L  !# Ignore unknown args
      !@return
      integer :: F_istat
      !*@/
      integer, parameter :: NMAX = 32

      integer, save :: default_wload = 0
      integer, save :: default_hzd   = 0
      integer, save :: default_monot = 1
      integer, save :: default_massc = 0
      real, save :: default_vmin  = -1*huge(1.)
      real, save :: default_vmax  = huge(1.)
      character(len=32), save :: default_intp_S  = 'TRICUB'

      character(len=2048) :: attributes_S, list_S(NMAX), &
           key_S, val_S, intp_S
      integer :: istat, n, n0
      integer :: wload, hzd, monot, massc
      real    :: vmin, vmax
      logical :: change_default_L, ignore_L
      !---------------------------------------------------------------
      F_istat = RMN_OK

      ignore_L = .false.
      if (present(F_ignore_L)) ignore_L = F_ignore_L

      attributes_S = F_attributes_S
      call str_tab2space(attributes_S)
      attributes_S = adjustl(attributes_S)
      istat = clib_toupper(attributes_S)

      list_S = ' '
      call str_split2list(list_S, attributes_S, ',', NMAX)

      change_default_L = (list_S(1) == 'DEFAULT')
      n0 = 1
      if (change_default_L) n0 = 2

      wload = default_wload
      hzd   = default_hzd
      monot = default_monot
      massc = default_massc
      vmin  = default_vmin
      vmax  = default_vmax
      intp_S  = default_intp_S

      LOOP1: do n = n0, NMAX
         if (list_S(n) == ' ') cycle
         call str_split(key_S, val_S, list_S(n), '=')
         select case(key_S)
         case('WLOAD')
            F_istat = min(F_istat, str_toint(wload, val_S))
         case('HZD')
            F_istat = min(F_istat, str_toint(hzd, val_S))
         case('MONO')
            F_istat = min(F_istat, str_toint(monot, val_S))
         case('MONOT')
            F_istat = min(F_istat, str_toint(monot, val_S))
         case('MASS')
            F_istat = min(F_istat, str_toint(massc, val_S))
         case('MASSC')
            F_istat = min(F_istat, str_toint(massc, val_S))
         case('MIN')
            F_istat = min(F_istat, str_toreal(vmin, val_S))
         case('VMIN')
            F_istat = min(F_istat, str_toreal(vmin, val_S))
         case('MAX')
            F_istat = min(F_istat, str_toreal(vmax, val_S))
         case('VMAX')
            F_istat = min(F_istat, str_toreal(vmax, val_S))
         case('INTP')
            intp_S  = val_S
         case default
            if (.not.ignore_L) then
               call msg(MSG_WARNING, '(tracers_attributes) unknown key:'//trim(key_S))
               F_istat = RMN_ERR
            endif
         end select
      end do LOOP1

      if (change_default_L) then
         default_wload = wload
         default_hzd   = hzd
         default_monot = monot
         default_massc = massc
         default_vmin  = vmin
         default_vmax  = vmax
         default_intp_S  = intp_S
      endif

      if (present(F_wload)) F_wload = wload
      if (present(F_hzd))   F_hzd   = hzd
      if (present(F_monot)) F_monot = monot
      if (present(F_massc)) F_massc = massc
      if (present(F_vmin))  F_vmin  = vmin
      if (present(F_vmax))  F_vmax  = vmax
      if (present(F_intp_S))  F_intp_S  = intp_S
      !---------------------------------------------------------------
      return
   end function tracers_attributes

end module tracers_attributes_mod
