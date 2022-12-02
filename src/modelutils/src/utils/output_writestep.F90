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
function output_writestep(F_id,F_step,F_callback) result(F_istat)
   use gmmx_mod
   use ptr_store
   use output_mod,only: output_getlist,output_writevar,output_close
   implicit none
   !@objective
   !@arguments
   integer,intent(in) :: F_id,F_step
   integer,external :: F_callback
   !@return
   integer :: F_istat
   !@author Stephane Chamberland, 2012-01
   !@description
   !  integer function callback(F_step,name_s,outname_s,data3d,lijk,uijk) result(istat)
   !     integer,intent(in) :: F_step,lijk(3),uijk(3)
   !     character(len=*),intent(inout) :: name_s    !VN in GMMX
   !     character(len=*),intent(inout) :: outname_s !ON in GMMX
   !     real,intent(inout) :: data3d(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3))
   !  end function callback
   !
   !  The callback function is called twice
   !  1) input:   F_step, name_s=' ', outname_s/=' '
   !     output:  name_s, outname_s
   !     return:  istat = RMN_ERR, 0 ,1
   !              return 0 if next call to F_callback will NOT modify the field
   !              return 1 if next call to F_callback will modify the field
   !     ignored: data3d, lijk, uijk
   !     action:  callback can fill a gmmx var with values 
   !              then return the name_s/outname_s to find it
   !  2) input:   F_step, name_s/=' ', outname_s/=' ', data3d
   !     output:  data3d
   !     return:  istat = RMN_ERR, 0, 1
   !              return 0 if ok to convert units for output with vardict
   !              return 1 to skip unit conversion (var ready to write as is)
   !     ignored: 
   !     action:  callback can modify the data (change units...)
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   logical,parameter :: NOPRINT = .false.
   integer,parameter :: MAX_ITEM = 999
   character(len=4)  :: mylist_S(MAX_ITEM)
   character(len=64)  :: name_S,outname_S,busname_S,hgrid_S,vgrid_S
   character(len=512) :: msg_S
   integer :: istat,ivar,n_items,idx,lijk(3),uijk(3)
   logical :: take_copy_L
   real,pointer :: data2dr4(:,:),data3dr4(:,:,:),data3dr4b(:,:,:)
   real :: dummy3d(1,1,1)
   !----------------------------------------------------------------------
   nullify(data2dr4,data3dr4,data3dr4b)
   F_istat = RMN_OK
   n_items = output_getlist(F_id,F_step,mylist_S)

   write(msg_S,'(a,I5.5)') '(Output) Step=',F_step
   if (n_items > 0) then
      call msg(MSG_INFO,trim(msg_S)//' [BEGIN]')
   endif

   VARLOOP: do ivar = 1,n_items
      call msg(MSG_INFOPLUS,trim(msg_S)//' var='//trim(mylist_S(ivar)))
      outname_S = mylist_S(ivar)
      name_S = ''
      lijk = 1; uijk = 1
      istat = F_callback(F_step,name_S,outname_S,dummy3d,lijk,uijk)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(Output) Call back step 1 error, ignoring: (ON='//trim(mylist_S(ivar))//')')
         F_istat = RMN_ERR
         cycle
      endif
      take_copy_L = (istat == 1)

      busname_S = ''
      hgrid_S = ''
      vgrid_S = ''
      idx = gmmx_meta(name_S,  &
           F_nameout_S=outname_S, &
           F_vb_S=busname_S,&
           F_hgrid_S=hgrid_S,&
           F_vgrid_S=vgrid_S)
      nullify(data2dr4,data3dr4)
      if (name_S == ' ' .or..not.RMN_IS_OK(idx)) then
         call msg(MSG_WARNING,'(Output) Could find any matching var, ignoring: '//trim(name_S)//' (ON='//trim(outname_S)//')')
         cycle
      endif
      outname_S = mylist_S(ivar)

      IF_GMM: if (any(busname_S(1:1)==(/'g','G'/))) then
         istat = gmmx_data(name_S,data2dr4,NOPRINT)
         if (.not.RMN_IS_OK(istat)) &
              istat = gmmx_data(name_S,data3dr4,NOPRINT)
         !TODO-later: support other type/shape
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING,'(Output) Could not get GMM var, ignoring: '//trim(name_S)//' (ON='//trim(outname_S)//')')
            F_istat = RMN_ERR
            cycle
         endif
      else
         call msg(MSG_WARNING,'(Output) Bus var not yet supported, ignoring: ' &
              //trim(name_S)//' (ON='//trim(outname_S)//' ; VB='// &
              trim(busname_S)//')')
         !TODO-later: implement bus unfolding
         F_istat = RMN_ERR
         cycle
      endif IF_GMM

      nullify(data3dr4b)
      IF_2D: if (associated(data2dr4)) then
         lijk(1:2) = lbound(data2dr4) ; lijk(3) = 1
         uijk(1:2) = ubound(data2dr4) ; uijk(3) = 1
         take_copy_L = .true. !TODO: find a way to associate ptr3d => ptr2d
         if (take_copy_L) then
            call ptr_store_get(data3dr4b,lijk,uijk)
            data3dr4b(:,:,1) = data2dr4(:,:)
!!$         else
!!$            data3dr4b(:,:,1:1) => data2dr4(:,:)
         endif
      else if (associated(data3dr4)) then
         lijk = lbound(data3dr4)
         uijk = ubound(data3dr4)
         if (take_copy_L) then
            call ptr_store_get(data3dr4b,lijk,uijk)
            data3dr4b = data3dr4
         else
            data3dr4b => data3dr4
         endif
      else
         call msg(MSG_WARNING,'(Output) Pointer not associated, ignoring: '//trim(name_S)//' (ON='//trim(outname_S)//')')
         F_istat = RMN_ERR
         cycle
      endif IF_2D

      if (.not.associated(data3dr4b)) then
         call msg(MSG_WARNING,'(Output) Could not get a working ptr, ignoring: '//trim(name_S)//' (ON='//trim(outname_S)//')')
         F_istat = RMN_ERR
         cycle
      endif

      istat = F_callback(F_step,name_S,outname_S,data3dr4b,lijk,uijk)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(Output) Call back step 2 error, ignoring: '//trim(name_S)//' (ON='//trim(outname_S)//')')
         F_istat = RMN_ERR
         cycle
!!$         else if (istat == 0) then
!!$            !TODO: Try automatic conversion with vardict... on a copy of the field
      endif
      F_istat = min(output_writevar(F_id,F_step,outname_S,data3dr4b,hgrid_S,vgrid_S),F_istat)
      if (take_copy_L) call ptr_store_free(data3dr4b)

   enddo VARLOOP

   istat = output_close(F_id,F_step)

   if (n_items>=0) then
      if (RMN_IS_OK(F_istat)) then
         F_istat = n_items
         call msg(MSG_INFO,trim(msg_S)//' [END] ok')
      else
         call msg(MSG_INFO,trim(msg_S)//' [END] with errors')
      endif
   endif
   !----------------------------------------------------------------------
   return
end function output_writestep

