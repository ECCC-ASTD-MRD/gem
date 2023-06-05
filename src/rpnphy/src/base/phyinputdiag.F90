 
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

module phyinputdiag
   use phy_options
   use phymem, only: phymeta, pvarlist, phymem_find, PHY_MAXVARS
   private
   public :: phyinputdiag1, phyinputdiag_id, phyinputdiag_obj

   logical,parameter:: SHORTMATCH_L = .true.
   logical,parameter:: QUIET_L = .true.
  
   interface phyinputdiag1
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
#include <rmn/msg.h>
      logical,parameter:: SHORTMATCH_L = .true.
      character(len=256) :: incfg_S
      character(len=32) :: inname_S,prefix_S,basename_S,time_S,ext_S
      integer :: istat,ivar,nvars, ivalist(PHY_MAXVARS)
      type(phymeta), pointer :: vmeta
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

         nvars = phymem_find(ivalist, 'tr/', 'V', 'D', QUIET_L, SHORTMATCH_L)
         do ivar = 1,nvars
            vmeta => pvarlist(ivalist(ivar))%meta
            call gmmx_name_parts(vmeta%vname,prefix_S,basename_S,time_S,ext_S)
            if  (vmeta%vname /= 'tr/hu:m' .and. &
                 vmeta%vname /= 'tr/hu:p' .and. &
                 all(time_S /= (/':M',':m'/))) then
               inname_S = vmeta%iname
               if (inname_S == ' ') inname_S = vmeta%oname
               if (inname_S == ' ') inname_S = vmeta%vname
               if (any(dyninread_list_s(:) == inname_S)) then
                  istat = min(input_add(F_inputid,'in='//trim(inname_S)//'; '//trim(incfg_S)),istat)
                  call msg(MSG_INFO,'(phyinputdiag) Adding: '//trim(inname_S)//' for '//trim(vmeta%vname))
               else
                  call msg(MSG_INFO,'(phyinputdiag) Skipping (not read by dyn): '//trim(inname_S)//' for '//trim(vmeta%vname))
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
#include <rmn/msg.h>
      logical, parameter:: SHORTMATCH_L = .true.
      character(len=256) :: incfg_S
      character(len=32) :: inname_S, prefix_S, basename_S, time_S, ext_S
      integer :: istat, ivar, nvars, ivalist(PHY_MAXVARS)
      type(phymeta), pointer :: vmeta
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

         nvars = phymem_find(ivalist, 'tr/', 'V', 'D', QUIET_L, SHORTMATCH_L)
         do ivar = 1, nvars
            vmeta => pvarlist(ivalist(ivar))%meta
            call gmmx_name_parts(vmeta%vname, prefix_S, basename_S, &
                 time_S, ext_S)
            if  (vmeta%vname /= 'tr/hu:m' .and. &
                 vmeta%vname /= 'tr/hu:p' .and. &
                 all(time_S /= (/':M', ':m'/))) then
               inname_S = vmeta%iname
               if (inname_S == ' ') inname_S = vmeta%oname
               if (inname_S == ' ') inname_S = vmeta%vname
               if (any(dyninread_list_s(:) == inname_S)) then
                  istat = min(inputio_add(F_inputobj%cfg, 'in='//trim(inname_S)// &
                       '; '//trim(incfg_S)), istat)
                  call msg(MSG_INFO, '(phyinputdiag) Adding: '//trim(inname_S)// &
                       ' for '//trim(vmeta%vname))
               else
                  call msg(MSG_INFO, &
                       '(phyinputdiag) Skipping (not read by dyn): '// &
                       trim(inname_S)//' for '//trim(vmeta%vname))
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

end module phyinputdiag
