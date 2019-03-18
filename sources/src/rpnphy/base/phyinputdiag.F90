 
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

module phyinputdiag_mod
   use phy_options
   use phy_typedef, only: phymeta
   use phy_getmeta_mod, only: phy_getmeta
   private
   public :: phyinputdiag, phyinputdiag_id, phyinputdiag_obj

   interface phyinputdiag
      module procedure phyinputdiag_id
      module procedure phyinputdiag_obj
   end interface

contains
 
   !/@*
   subroutine phyinputdiag_id(F_inputid)
      use input_mod
      implicit none
      !@objective Add diag level variables to input list
      !@arguments
      integer,intent(in) :: F_inputid
      !@author Stephane Chamberland,2014-11
      !*@/
#include <rmnlib_basics.hf>
#include <msg.h>
      logical,parameter:: SHORTMATCH_L = .true.
      integer,parameter:: MYMAX = 1024
      character(len=256) :: incfg_S
      character(len=32) :: inname_S,prefix_S,basename_S,time_S,ext_S
      integer :: istat,ivar,nvars
      type(phymeta), pointer :: metalist(:)
      ! ---------------------------------------------------------------------
      istat = RMN_OK

      IF_ISLIST: if (indiag_list_s(1) == 'DEFAULT LIST')  then

         call msg(MSG_INFO,'(phyinputdiag) Adding default in-diag list')

         incfg_S =  'freq=0; search=ANAL ; hinterp=cubic; mandatory=false; levels=mdiag'
         istat = min(input_add(F_inputid,'in=UU; IN2=VV; '//trim(incfg_S)),istat)

         call msg(MSG_INFO,'(phyinputdiag) Adding: UU/VV')
         incfg_S =  'freq=0; search=ANAL ; hinterp=cubic; mandatory=false; levels=tdiag'
         istat = min(input_add(F_inputid,'in=TT; '//trim(incfg_S)),istat)
         call msg(MSG_INFO,'(phyinputdiag) Adding: TT')
         istat = min(input_add(F_inputid,'in=HU; '//trim(incfg_S)),istat)
         call msg(MSG_INFO,'(phyinputdiag) Adding: HU')

         nullify(metalist)
         nvars = phy_getmeta(metalist, 'tr/', F_npath='V', F_bpath='D', &
              F_maxmeta=MYMAX, F_shortmatch=SHORTMATCH_L)
         do ivar = 1,nvars
            call gmmx_name_parts(metalist(ivar)%vname,prefix_S,basename_S,time_S,ext_S)
            if  (metalist(ivar)%vname /= 'tr/hu:m' .and. &
                 metalist(ivar)%vname /= 'tr/hu:p' .and. &
                 all(time_S /= (/':M',':m'/))) then
               inname_S = metalist(ivar)%iname
               if (inname_S == ' ') inname_S = metalist(ivar)%oname
               if (inname_S == ' ') inname_S = metalist(ivar)%vname
               if (any(dyninread_list_s(:) == inname_S)) then
                  istat = min(input_add(F_inputid,'in='//trim(inname_S)//'; '//trim(incfg_S)),istat)
                  call msg(MSG_INFO,'(phyinputdiag) Adding: '//trim(inname_S)//' for '//trim(metalist(ivar)%vname))
               else
                  call msg(MSG_INFO,'(phyinputdiag) Skipping (not read by dyn): '//trim(inname_S)//' for '//trim(metalist(ivar)%vname))
               endif
            endif
         enddo

      else if (indiag_list_s(1) /= ' ') then

         call msg(MSG_INFO,'(phyinputdiag) Adding user provided in-diag list')

         do ivar = 1,size(indiag_list_s)
            if (indiag_list_s(ivar) == ' ') exit
            incfg_S =  '; freq=0; search=ANAL ; hinterp=cubic; mandatory=false; levels='
            if (any(indiag_list_s(ivar) == (/'UU','VV'/))) then
               incfg_S = trim(incfg_S)//'mdiag'
            else
               incfg_S = trim(incfg_S)//'tdiag'
            endif
            istat = min(input_add(F_inputid,'in='//trim(indiag_list_s(ivar))//trim(incfg_S)),istat)
            call msg(MSG_INFO,'(phyinputdiag) Adding: '//trim(indiag_list_s(ivar)))
         enddo

      endif IF_ISLIST

      if (.not.RMN_IS_OK(istat)) &
           call msg(MSG_WARNING,'(phyinputdiag) Problem registering all resquested diag vars')
      ! ---------------------------------------------------------------------
      return
   end subroutine phyinputdiag_id

   !/@*
   subroutine phyinputdiag_obj(F_inputobj)
      use inputio_mod
      implicit none
      !@objective Add diag level variables to input list
      !@arguments
      type(INPUTIO_T), intent(inout) :: F_inputobj
      !@author Stephane Chamberland, 2017-09
      !*@/
#include <rmnlib_basics.hf>
#include <msg.h>
      logical, parameter:: SHORTMATCH_L = .true.
      integer, parameter:: MYMAX = 1024
      character(len=256) :: incfg_S
      character(len=32) :: inname_S, prefix_S, basename_S, time_S, ext_S
      integer :: istat, ivar, nvars
      type(phymeta), pointer :: metalist(:)
      ! ---------------------------------------------------------------------
      istat = RMN_OK

      IF_ISLIST: if (indiag_list_s(1) == 'DEFAULT LIST') then

         call msg(MSG_INFO, '(phyinputdiag) Adding default in-diag list')

         incfg_S = 'freq=0; search=ANAL ; hinterp=cubic; mandatory=false; levels=mdiag'
         istat = min(istat, &
              inputio_add(F_inputobj%cfg, 'in=UU; IN2=VV; '//trim(incfg_S)))

         call msg(MSG_INFO, '(phyinputdiag) Adding: UU/VV')
         incfg_S =  'freq=0; search=ANAL ; hinterp=cubic; mandatory=false; levels=tdiag'
         istat = min(inputio_add(F_inputobj%cfg, 'in=TT; '//trim(incfg_S)), istat)
         call msg(MSG_INFO, '(phyinputdiag) Adding: TT')
         istat = min(inputio_add(F_inputobj%cfg, 'in=HU; '//trim(incfg_S)), istat)
         call msg(MSG_INFO, '(phyinputdiag) Adding: HU')

         nullify(metalist)
         nvars = phy_getmeta(metalist, 'tr/', F_npath='V', F_bpath='D', &
              F_maxmeta=MYMAX, F_shortmatch=SHORTMATCH_L)
         do ivar = 1, nvars
            call gmmx_name_parts(metalist(ivar)%vname, prefix_S, basename_S, &
                 time_S, ext_S)
            if  (metalist(ivar)%vname /= 'tr/hu:m' .and. &
                 metalist(ivar)%vname /= 'tr/hu:p' .and. &
                 all(time_S /= (/':M', ':m'/))) then
               inname_S = metalist(ivar)%iname
               if (inname_S == ' ') inname_S = metalist(ivar)%oname
               if (inname_S == ' ') inname_S = metalist(ivar)%vname
               if (any(dyninread_list_s(:) == inname_S)) then
                  istat = min(inputio_add(F_inputobj%cfg, 'in='//trim(inname_S)// &
                       '; '//trim(incfg_S)), istat)
                  call msg(MSG_INFO, '(phyinputdiag) Adding: '//trim(inname_S)// &
                       ' for '//trim(metalist(ivar)%vname))
               else
                  call msg(MSG_INFO, &
                       '(phyinputdiag) Skipping (not read by dyn): '// &
                       trim(inname_S)//' for '//trim(metalist(ivar)%vname))
               endif
            endif
         enddo

      else if (indiag_list_s(1) /= ' ') then

         call msg(MSG_INFO, '(phyinputdiag) Adding user provided in-diag list')

         do ivar = 1, size(indiag_list_s)
            if (indiag_list_s(ivar) == ' ') exit
            incfg_S = '; freq=0; search=ANAL ; hinterp=cubic; mandatory=false; levels='
            if (any(indiag_list_s(ivar) == (/'UU', 'VV'/))) then
               incfg_S = trim(incfg_S)//'mdiag'
            else
               incfg_S = trim(incfg_S)//'tdiag'
            endif
            istat = min(istat, inputio_add(F_inputobj%cfg, &
                 'in='//trim(indiag_list_s(ivar))//trim(incfg_S)))
            call msg(MSG_INFO, '(phyinputdiag) Adding: '//trim(indiag_list_s(ivar)))
         enddo

      endif IF_ISLIST

      if (.not.RMN_IS_OK(istat)) &
           call msg(MSG_WARNING, '(phyinputdiag) Problem registering all resquested diag vars')
      ! ---------------------------------------------------------------------
      return
   end subroutine phyinputdiag_obj

end module phyinputdiag_mod
