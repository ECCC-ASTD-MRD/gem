!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

!/@*
module vinterp_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use vGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use vgrid_wb
   use samevgrid_mod
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2011-04
   !@description
   ! Public functions
   public :: vinterp,vinterp0,vinterp1,vinterp10,vinterp01
   ! Public constants
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
#include <mu_gmm.hf>

   interface vinterp
      module procedure vinterp0
      module procedure vinterp0ls
      module procedure vinterp1
      module procedure vinterp10
      module procedure vinterp01
   end interface

   real, parameter :: EPSILON_R4 = 1.e-6
   real, parameter :: P0MIN = 1.
   real, parameter :: P0MAX = 1.e6


contains

   !/@*
   function vinterp0( &
        F_dataout, F_vgridout, F_ip1listout, &
        F_datain, F_vgridin, F_ip1listin, &
        F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgridout,F_vgridin
      integer,intent(in) :: F_ip1listout(:),F_ip1listin(:)
      real,pointer :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S, tmp_S
      integer :: nlinbot
      real,pointer :: sfcfldout2(:,:), sfcfldin2(:,:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp0 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      if (.not.associated(F_datain)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, null pointer: '//trim(msg_S))
         return
      endif
      nlinbot = 0
      if (present(F_nlinbot)) nlinbot = F_nlinbot
      nullify(sfcfldout2, sfcfldin2)
      F_istat = vinterp0ls( &
           F_dataout, F_vgridout, F_ip1listout, F_sfcfldout, sfcfldout2, &
           F_datain, F_vgridin, F_ip1listin, F_sfcfldin, sfcfldin2, &
           F_nlinbot, F_msg_S)
      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp0 [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp0


   !/@*
   function vinterp0ls( &
        F_dataout, F_vgridout, F_ip1listout, F_sfcfldout, F_sfcfldout2,&
        F_datain, F_vgridin, F_ip1listin, F_sfcfldin, F_sfcfldin2, &
        F_nlinbot, F_msg_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgridout,F_vgridin
      integer,intent(in) :: F_ip1listout(:),F_ip1listin(:)
      real,pointer :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S, tmp_S
      integer,pointer :: samevert_list(:)
      integer :: istat, nlinbot, nijkout(3), nijkin(3), lijk(3), uijk(3), k
      integer, dimension(size(F_datain,dim=3),1) :: slist
      real,target :: sfcfld0(1,1)
      real,pointer :: levelsin(:,:,:), levelsout(:,:,:), sfcfldin(:,:), sfcfldin2(:,:)
      real,dimension(size(F_datain,dim=3)) :: scol
      real,dimension(size(F_datain,dim=1),size(F_datain,dim=2),size(F_datain,dim=3)) :: sdatain, slevelsin
!!$      real,dimension(size(F_datain,dim=1),size(F_datain,dim=2),size(F_datain,dim=3)),target :: sdatain, slevelsin
!!$      real,dimension(size(F_datain,dim=1),size(F_datain,dim=2),1),target :: ssfcsin
!!$      real,pointer :: datatmp(:,:,:)

      logical :: use_same_sfcfld_L, samevert_L, ok_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp0ls [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      if (.not.associated(F_datain)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, null pointer: '//trim(msg_S))
         return
      endif
      nlinbot = 0
      if (present(F_nlinbot)) nlinbot = F_nlinbot
      nullify(sfcfldin2)
      if (associated(F_sfcfldin) .and. .not.associated(F_sfcfldin,F_sfcfldout)) then
         use_same_sfcfld_L = .false.
         sfcfldin => F_sfcfldin
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
!!$         print *,'(vinterp) sfcfldin => F_sfcfldin',use_same_sfcfld_L
      else
         use_same_sfcfld_L = .true.
         sfcfldin => F_sfcfldout
         if (associated(F_sfcfldout2)) sfcfldin2 => F_sfcfldout2
!!$         print *,'(vinterp) sfcfldin => F_sfcfldout',use_same_sfcfld_L
      endif

      nijkin  = shape(F_datain)
      if (.not.associated(F_dataout)) then
         lijk = lbound(F_datain) ; uijk = ubound(F_datain)
         allocate(F_dataout(lijk(1):uijk(1),lijk(2):uijk(2),size(F_ip1listout))) !#,stat=istat
      endif
      nijkout = shape(F_dataout)
      ok_L = .true.
      if (associated(F_sfcfldout2)) ok_L = all(nijkout(1:2) == shape(F_sfcfldout2))
      if (ok_L .and. associated(sfcfldin2)) ok_L = all(nijkout(1:2) == shape(sfcfldin2))
      if (.not.( &
           ok_L .and. &
           all(nijkout(1:2) == nijkin(1:2)) .and. &
           all(nijkout(1:2) == shape(F_sfcfldout)) .and. &
           all(nijkout(1:2) == shape(sfcfldin)) .and. &
           nijkout(3) == size(F_ip1listout) .and. &
           nijkin(3)  == size(F_ip1listin) &
           )) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, not same shape: '//trim(msg_S))
         print *,'(vinterp) ok_L=', ok_L
         print *,'(vinterp) nijkout=', nijkout
         if (associated(F_sfcfldout2)) &
              print *,'(vinterp) shape(F_sfcfldout2)=', shape(F_sfcfldout2)
         if (ok_L .and. associated(sfcfldin2)) &
              print *,'(vinterp) shape(sfcfldin2)=', shape(sfcfldin2)
         print *,'(vinterp) nijkin(1:2)=', nijkin(1:2)
         print *,'(vinterp) shape(F_sfcfldout)=', shape(F_sfcfldout)
         print *,'(vinterp) shape(sfcfldin)=', shape(sfcfldin)
         print *,'(vinterp) size(F_ip1listout)=', size(F_ip1listout)
         print *,'(vinterp) nijkin(3)=', nijkin(3)
         print *,'(vinterp) size(F_ip1listin)=', size(F_ip1listin)
         call flush(6)
         return
      endif

      allocate(samevert_list(nijkout(3))) !#,stat=istat
      samevert_L = samevgrid(F_vgridin,F_ip1listin,F_vgridout,F_ip1listout,samevert_list)
      if (samevert_L .and. (.not.use_same_sfcfld_L) .and. &
           associated(F_sfcfldout) .and. associated(sfcfldin)) then
         if (any(abs(F_sfcfldout-sfcfldin) > EPSILON_R4)) samevert_L = .false.
         !#TODO: if (any(abs(F_sfcfldout2-sfcfldin2) > EPSILON_R4)) samevert_L = .false.
      endif

      tmp_S = 'Interpolating'
      nullify(levelsout,levelsin)
      if (samevert_L) then
         tmp_S = 'No Interpolation for'
         sfcfld0 = 100000. !TODO: make sure F_vgridin is pressure based
         sfcfldin => sfcfld0
         sfcfldin2 => sfcfld0
         F_istat = RMN_OK
      else
         !TODO: if vgrid needs sfcfield and sfcfld => null then error
         if (minval(F_sfcfldout) < P0MIN .or. maxval(F_sfcfldout) > P0MAX) then
            call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, sfcfldout provided '// &
                 'sfcfldout has wrong values: '//trim(msg_S))
            F_istat = RMN_ERR
            return
         endif
         if (associated(F_sfcfldout2)) then
            F_istat = vgd_levels(F_vgridout, F_ip1listout, levelsout, &
                 F_sfcfldout, in_log=.false., sfc_field_ls=F_sfcfldout2)
         else
            F_istat = vgd_levels(F_vgridout, F_ip1listout, levelsout, &
                 F_sfcfldout, in_log=.false.)
         endif
      endif
      if (minval(sfcfldin) < P0MIN .or. maxval(sfcfldin) > P0MAX) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, provided sfcfldin has wrong values: '//trim(msg_S))
         F_istat = RMN_ERR
         return
      endif
      if (associated(sfcfldin2)) then
         F_istat = min(vgd_levels(F_vgridin, F_ip1listin, levelsin, &
              sfcfldin, in_log=.false., sfc_field_ls=sfcfldin2), F_istat)
      else
         F_istat = min(vgd_levels(F_vgridin, F_ip1listin, levelsin, &
              sfcfldin, in_log=.false.), F_istat)
      endif

      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, problem getting levels: '//trim(msg_S))
         return
      endif
      call msg(MSG_INFO,'(vinterp) '//trim(tmp_S)//': '//trim(msg_S))


      if (samevert_L) then
!!$         !TODO: check that output levels are increasing as well
         do k=1,size(samevert_list)
            F_dataout(:,:,k) = F_datain(:,:,samevert_list(k))
         enddo
      else

         ! Sort input levels into increasing pressures
         !#TODO: sorting should be done in vgrid from file, could be ascending or
         ! descending depending on model (ascending for GEM), only check monotonicity in here
         scol = levelsin(1,1,:)
         do k=1,size(scol)
            slist(k,:) = minloc(scol)
            scol(slist(k,1)) = maxval(scol)+1.
         enddo

!!$         ssfcsin(:,:,1) = sfcfldin(:,:)
!!$         datatmp => ssfcsin(:,:,:)
!!$         call statfld_dm(datatmp, 'sfcfldin', 0, 'vinterp', 8)

         do k=1,size(scol)
            sdatain(:,:,k) = F_datain(:,:,slist(k,1))
            slevelsin(:,:,k) = levelsin(:,:,slist(k,1))
         enddo

!!$         do k=1,size(scol)
!!$             write(tmp_S, '(i3)') k
!!$             print *,'(vinterp) k,ip1,nlinbot=',k,F_ip1listin(k),nlinbot
!!$             datatmp => slevelsin(:,:,k:k)
!!$             call statfld_dm(datatmp, trim(tmp_S)//' slevelsin ', 0, 'vinterp', 8)
!!$          enddo
!!$         do k=1,size(scol)
!!$            write(tmp_S, '(i3)') k
!!$            datatmp => sdatain(:,:,k:k)
!!$            call statfld_dm(datatmp, trim(tmp_S)//' sdatain ', 0, 'vinterp', 8)
!!$         enddo

!!$         print *,'(vinterp) '//trim(msg_S)//' nlinbot=',nlinbot,nijkout
!!$         print *,'(vinterp) '//trim(msg_S)//' sfcfldin=',minval(sfcfldin),maxval(sfcfldin), &
!!$              sum(dble(sfcfldin))/dble(float(size(sfcfldin)))
!!$         print *,'(vinterp) '//trim(msg_S)//' slevelsin=',minval(slevelsin),maxval(slevelsin), &
!!$              sum(dble(slevelsin))/dble(float(size(slevelsin)))
!!$         print *,'(vinterp) '//trim(msg_S)//' F_sfcfldout=',minval(F_sfcfldout), &
!!$              maxval(F_sfcfldout),sum(dble(F_sfcfldout))/dble(float(size(F_sfcfldout)))
!!$         print *,'(vinterp) '//trim(msg_S)//' levelsout=',minval(levelsout),maxval(levelsout), &
!!$              sum(dble(levelsout))/dble(float(size(levelsout)))
!!$         print *,'(vinterp) '//trim(msg_S)//' sdatain=',minval(sdatain),maxval(sdatain), &
!!$              sum(dble(sdatain))/dble(float(size(sdatain)))

         call vte_intvertx4(F_dataout,sdatain,slevelsin,levelsout, &
              nijkout(1)*nijkout(2),nijkin(3),nijkout(3),&
              msg_S,nlinbot)

!!$         print *,'(vinterp) '//trim(msg_S)//' F_dataout=',minval(F_dataout), &
!!$              maxval(F_dataout),sum(F_dataout)/float(size(F_dataout))

     endif

      if (associated(samevert_list)) deallocate(samevert_list,stat=istat)
      if (associated(levelsout)) deallocate(levelsout,stat=istat)
      if (associated(levelsin)) deallocate(levelsin,stat=istat)
      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp0ls [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp0ls


   !/@*
   function vinterp1(F_dataout, F_vgridout_S, F_datain, F_vgridin_S, &
        F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S, F_use_same_sfcfld_L, &
        F_sfcfldout2, F_sfcfldin2) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      character(len=*),intent(in) :: F_vgridout_S,F_vgridin_S
      real,pointer,optional :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer,optional :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      logical,intent(in),optional :: F_use_same_sfcfld_L
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S,tmp_S
      integer :: nlinbot, istat
      type(vgrid_descriptor) :: vgridout,vgridin
      real,pointer :: sfcfldout(:,:),sfcfldin(:,:)
      real,pointer :: sfcfldout2(:,:),sfcfldin2(:,:)
      integer,pointer :: ip1listout(:),ip1listin(:)
      logical :: use_same_sfcfld_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp1 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nlinbot = 0 
      if (present(F_nlinbot)) nlinbot = F_nlinbot
      use_same_sfcfld_L = .false.
      if (present(F_use_same_sfcfld_L)) use_same_sfcfld_L = F_use_same_sfcfld_L

      nullify(ip1listout, ip1listin, sfcfldout, sfcfldin, sfcfldout2, sfcfldin2)
      if (present(F_sfcfldout)) then
         if (associated(F_sfcfldout)) sfcfldout => F_sfcfldout
      endif
      if (present(F_sfcfldin)) then
         if (associated(F_sfcfldin)) sfcfldin => F_sfcfldin
      endif
      if (present(F_sfcfldout2)) then
         if (associated(F_sfcfldout2)) sfcfldout2 => F_sfcfldout2
      endif
      if (present(F_sfcfldin2)) then
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
      endif

!!$      print *,'(vinterp) use_same_sfcfld_L=',use_same_sfcfld_L, associated(sfcfldout), associated(sfcfldin)

      F_istat = priv_vgrid_details(F_vgridout_S, vgridout, ip1listout, &
           sfcfldout, sfcfldout2)
      if (use_same_sfcfld_L) sfcfldin => sfcfldout
      if (use_same_sfcfld_L) sfcfldin2 => sfcfldout2
      if (RMN_IS_OK(F_istat)) &
           F_istat = priv_vgrid_details(F_vgridin_S, vgridin, ip1listin, &
                                        sfcfldin, sfcfldin2)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, missing vgrid info: '//trim(msg_S))
         return
      endif

      F_istat = vinterp0ls( &
           F_dataout, vgridout, ip1listout, sfcfldout, sfcfldout2, &
           F_datain, vgridin, ip1listin, sfcfldin, sfcfldin2, nlinbot, msg_S)

      istat = vgd_free(vgridout)
      if (associated(ip1listout)) deallocate(ip1listout,stat=istat)

      istat = vgd_free(vgridin)
      if (associated(ip1listin)) deallocate(ip1listin,stat=istat)

      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp1 [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp1


   !/@*
   function vinterp01(F_dataout, F_vgridout, F_ip1listout, &
        F_datain, F_vgridin_S, F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S, &
        F_sfcfldout2, F_sfcfldin2) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgridout
      integer,intent(in) :: F_ip1listout(:)
      character(len=*),intent(in) :: F_vgridin_S
      real,pointer,optional :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer,optional :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S,tmp_S
      integer :: nlinbot, istat
      type(vgrid_descriptor) :: vgridin
      real,pointer :: sfcfldout(:,:),sfcfldin(:,:)
      real,pointer :: sfcfldout2(:,:),sfcfldin2(:,:)
      integer,pointer :: ip1listin(:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp01 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nlinbot = 0
      if (present(F_nlinbot)) nlinbot = F_nlinbot

      nullify(ip1listin, sfcfldout, sfcfldin, sfcfldout2, sfcfldin2)
      if (present(F_sfcfldout)) then
         if (associated(F_sfcfldout)) sfcfldout => F_sfcfldout
      endif
      if (present(F_sfcfldin)) then
         if (associated(F_sfcfldin)) sfcfldin => F_sfcfldin
      endif
      if (present(F_sfcfldout2)) then
         if (associated(F_sfcfldout2)) sfcfldout2 => F_sfcfldout2
      endif
      if (present(F_sfcfldin2)) then
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
      endif
      F_istat = priv_vgrid_details(F_vgridin_S, vgridin, ip1listin, &
           sfcfldin, sfcfldin2)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, vgrid info: '//trim(msg_S))
         return
      endif

      F_istat = vinterp0ls( &
           F_dataout, F_vgridout, F_ip1listout, sfcfldout, sfcfldout2, &
           F_datain, vgridin, ip1listin, sfcfldin, sfcfldin2, nlinbot, msg_S)

      istat = vgd_free(vgridin)
      if (associated(ip1listin)) deallocate(ip1listin,stat=istat)

      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp01 [END] '//trim(tmp_S))
     !------------------------------------------------------------------
      return
   end function vinterp01


   !/@*
   function vinterp10(F_dataout, F_vgridout_S, F_datain, F_vgridin, &
        F_ip1listin, F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S, &
        F_sfcfldout2, F_sfcfldin2) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      character(len=*),intent(in) :: F_vgridout_S
      type(vgrid_descriptor),intent(in) :: F_vgridin
      integer,intent(in) :: F_ip1listin(:)
      real,pointer,optional :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer,optional :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S,tmp_S
      integer :: nlinbot, istat
      type(vgrid_descriptor) :: vgridout
      real,pointer :: sfcfldout(:,:),sfcfldin(:,:)
      real,pointer :: sfcfldout2(:,:),sfcfldin2(:,:)
      integer,pointer :: ip1listout(:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp10 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nlinbot = 0 
      if (present(F_nlinbot)) nlinbot = F_nlinbot

      nullify(ip1listout, sfcfldout, sfcfldin, sfcfldout2, sfcfldin2)
      if (present(F_sfcfldout)) then
         if (associated(F_sfcfldout)) sfcfldout => F_sfcfldout
      endif
      if (present(F_sfcfldin)) then
         if (associated(F_sfcfldin)) sfcfldin => F_sfcfldin
      endif
      if (present(F_sfcfldout2)) then
         if (associated(F_sfcfldout2)) sfcfldout2 => F_sfcfldout2
      endif
      if (present(F_sfcfldin2)) then
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
      endif

      F_istat = priv_vgrid_details(F_vgridout_S, vgridout, ip1listout, &
           sfcfldout, sfcfldout2)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, vgrid info: '//trim(msg_S))
         return
      endif

      F_istat = vinterp0ls( &
           F_dataout, vgridout, ip1listout, sfcfldout, sfcfldout2, &
           F_datain, F_vgridin, F_ip1listin, sfcfldin, sfcfldin2, &
           nlinbot, msg_S)

      istat = vgd_free(vgridout)
      if (associated(ip1listout)) deallocate(ip1listout,stat=istat)

      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp10 [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp10


   !==== Private Functions =================================================


   function priv_vgrid_details(F_vgrid_S, F_vgrid, F_ip1list, &
        F_sfcfld, F_sfcfld2) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_vgrid_S
      type(vgrid_descriptor),intent(out) :: F_vgrid
      integer,pointer :: F_ip1list(:)
      real,pointer :: F_sfcfld(:,:)
      real,pointer :: F_sfcfld2(:,:)
      integer :: F_istat, istat, vtype
      character(len=32) :: sfcfld_S, sfcfld2_S
      !------------------------------------------------------------------
      sfcfld2_S = ' '
      call vgrid_nullify(F_vgrid)
      F_istat = vgrid_wb_get(F_vgrid_S, F_vgrid, F_ip1list, vtype, &
           sfcfld_S, sfcfld2_S)
      if (RMN_IS_OK(F_istat)) then
         if (sfcfld_S /= ' ' .and. .not.associated(F_sfcfld)) &
              istat = gmm_get(sfcfld_S, F_sfcfld)
         if (sfcfld2_S /= ' ' .and. .not.associated(F_sfcfld2)) &
              istat = gmm_get(sfcfld2_S, F_sfcfld2)
      endif
      if (.not.(RMN_IS_OK(F_istat) .and. &
           (associated(F_sfcfld) .or. sfcfld_S == ' '))) then
         istat = vgd_free(F_vgrid)
         if (associated(F_ip1list)) deallocate(F_ip1list, stat=istat)
         F_istat = RMN_ERR
         return
      endif
      if (sfcfld2_S /= ' ' .and. .not.associated(F_sfcfld2)) then
         istat = vgd_free(F_vgrid)
         if (associated(F_ip1list)) deallocate(F_ip1list, stat=istat)
         F_istat = RMN_ERR
      endif
      !------------------------------------------------------------------
      return
   end function priv_vgrid_details

end module vinterp_mod
