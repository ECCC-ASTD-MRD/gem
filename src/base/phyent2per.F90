!-------------------------------------- LICENCE BEGIN -----------------------------
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
!-------------------------------------- LICENCE END -------------------------------

!/@*
function phyent2per(F_readlist_S, F_nread, F_step) result(F_istat)
   use statfld_dm_mod
   use phy_get_mod, only: phy_get
   use phy_options
   use phy_put_mod, only: phy_put
   use phy_typedef, only: phymeta
   implicit none
   !@objective Transfer read vars from busent to busper if need be
   integer,intent(in) :: F_nread, F_step
   character(len=*),intent(in) :: F_readlist_S(F_nread)
   integer :: F_istat
   !*@/
!!!#include <arch_specific.hf>
#include <msg.h>
#include <rmnlib_basics.hf>
   integer, parameter :: STAT_PRECISION = 8
   integer, parameter :: NVARS_MAX = 256
   integer, parameter :: NAME_ENT = 1
   integer, parameter :: NAME_PER = 2
   integer :: istat, ivar, nvar
   character(len=32) :: var_duets_S(2,NVARS_MAX)
   type(phymeta) :: meta_e, meta_p
   real, pointer :: ptr_e(:,:,:), ptr_p(:,:,:)
   ! ---------------------------------------------------------------------
   F_istat = RMN_OK
   if (F_nread == 0 .or. (F_step == 0 .and. &
        any(fluvert == (/ &
        'MOISTKE', &
        'CLEF   ', &
        'SIMPLE '  &
        /)))) return
 
   nvar = 0
   var_duets_S = ' '

!!$   nvar = nvar + 1
!!$   var_duets_S(NAME_PER,nvar) = ''
!!$   var_duets_S(NAME_ENT,nvar) = ''

   do ivar = 1, nvar
      if (any(F_readlist_S(1:F_nread) == var_duets_S(NAME_ENT,ivar))) then

         nullify(ptr_e ,ptr_p)
         istat = phy_get(ptr_e, var_duets_S(NAME_ENT,ivar), &
              F_npath='V', F_meta=meta_e, F_quiet=.true.)!#, F_folded=.true.)
         istat = phy_get(ptr_p, var_duets_S(NAME_PER,ivar), &
              F_npath='V', F_meta=meta_p, F_quiet=.true.)!#, F_folded=.true.)

         if (.not.(associated(ptr_e).and.associated(ptr_p))) then
            F_istat = RMN_ERR
            call msg(MSG_WARNING,'(phy_ent2per) Problem getting '//trim(var_duets_S(NAME_PER,ivar)))
            return
         endif

         call msg(MSG_INFO,'(phy_ent2per) updating '//trim(var_duets_S(NAME_PER,ivar))//' from '//trim(var_duets_S(NAME_ENT,ivar)))
         ptr_p = ptr_e

         istat = phy_put(ptr_p, var_duets_S(NAME_PER,ivar), F_npath='V')

      endif
   enddo
   ! ---------------------------------------------------------------------
   return
end function phyent2per
