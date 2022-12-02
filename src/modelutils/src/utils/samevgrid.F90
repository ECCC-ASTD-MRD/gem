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
module samevgrid_mod
   use vGrid_Descriptors
   use vgrid_wb
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2014-11
   !@description
   ! Public functions
   public :: samevgrid,samevgrid_vg,samevgrid_ip,samevgrid_vgip,samevgrid_str
   ! Public constants
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface samevgrid
      module procedure samevgrid_vg
      module procedure samevgrid_ip
      module procedure samevgrid_vgip
      module procedure samevgrid_str
   end interface samevgrid

contains

   !/@*
   function samevgrid_vg(F_vgridin,F_vgridout) result(F_same_L)
      !@arguments
      type(vgrid_descriptor),intent(in) :: F_vgridout,F_vgridin
      logical :: F_same_L
      !*@/
      !------------------------------------------------------------------
      F_same_L = (F_vgridout == F_vgridin)
      !------------------------------------------------------------------
      return
   end function samevgrid_vg

   !/@*
   function samevgrid_ip(F_ip1listin,F_ip1listout,F_sublistin) result(F_same_L)
      !@arguments
      integer,intent(in) :: F_ip1listin(:),F_ip1listout(:)
      integer,intent(out),optional :: F_sublistin(:)
      logical :: F_same_L
      !*@/
      integer :: k,ki
      !------------------------------------------------------------------
      !#TODO: if .not.present(F_sublistin) should we return F_same_L = (F_ip1listout(k) == F_ip1listin)

      F_same_L = .true.
      if (present(F_sublistin)) F_sublistin = 1
      !# TODO: convert to ip1 newstyle before comparing
      do k=1,size(F_ip1listout)
         if (.not.any(F_ip1listout(k) == F_ip1listin(:))) then
            F_same_L = .false.
            exit
         endif
      enddo
      if (F_same_L .and.present(F_sublistin)) then
         if (size(F_sublistin) < size(F_ip1listout)) then
            call msg(MSG_WARNING,'(samevgrid) size(F_sublistin) < size(F_ip1listout)')
         endif
         do k=1,min(size(F_sublistin),size(F_ip1listout))
            do ki=1,size(F_ip1listin)
               if (F_ip1listout(k) == F_ip1listin(ki)) then
                  F_sublistin(k) = ki
                  cycle
               endif
            enddo
         enddo
      endif
      !------------------------------------------------------------------
      return
   end function samevgrid_ip

   !/@*
   function samevgrid_vgip(F_vgridin,F_ip1listin,F_vgridout,F_ip1listout,F_sublistin) result(F_same_L)
      !@arguments
      type(vgrid_descriptor),intent(in) :: F_vgridin,F_vgridout
      integer,intent(in) :: F_ip1listin(:),F_ip1listout(:)
      integer,intent(out),optional :: F_sublistin(:)
      logical :: F_same_L
      !*@/
!!$      integer :: k
      !------------------------------------------------------------------
      !#TODO: should F_vgridout, F_vgridin matter?
      !#TODO: maybe we should check ip list consistency with provided vgrid
      F_same_L = (F_vgridout == F_vgridin)
!!$      if (present(F_sublistin)) F_sublistin = 1
!!$      if (F_same_L) then
!!$         if (present(F_sublistin)) then
!!$            do k=1,min(size(F_sublistin),size(F_ip1listout))
!!$               F_sublistin(k) = k
!!$            enddo
!!$         endif
!!$      else
         if (present(F_sublistin)) then
            F_same_L = samevgrid_ip(F_ip1listin,F_ip1listout,F_sublistin)
         else
            F_same_L = samevgrid_ip(F_ip1listin,F_ip1listout)
         endif
!!$      endif
      !------------------------------------------------------------------
      return
   end function samevgrid_vgip

   !/@*
   function samevgrid_str(F_vgridin_S,F_vgridout_S,F_sublistin) result(F_same_L)
      !@arguments
      character(len=*),intent(in) :: F_vgridin_S,F_vgridout_S
      integer,intent(out),optional :: F_sublistin(:)
      logical :: F_same_L
      !*@/
      type(vgrid_descriptor) :: vgridin,vgridout
      integer,pointer :: ip1listin(:),ip1listout(:)
      character(len=32) :: sfcfldin_S,sfcfldout_S
      integer :: istat,vtypein,vtypeout
      !------------------------------------------------------------------
      F_same_L = .false.
      if (present(F_sublistin)) F_sublistin = 1
      nullify(ip1listin,ip1listout)
      istat = vgrid_wb_get(F_vgridin_S,vgridin,ip1listin,vtypein,sfcfldin_S)
      istat = min(vgrid_wb_get(F_vgridout_S,vgridout,ip1listout,vtypeout,sfcfldout_S),istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(samevgrid) Problem getting vgrid details')
         return
      endif
      if (present(F_sublistin)) then
         F_same_L = samevgrid_vgip(vgridin,ip1listin,vgridout,ip1listout,F_sublistin)
      else
         F_same_L = samevgrid_vgip(vgridin,ip1listin,vgridout,ip1listout)
      endif
      if (associated(ip1listin)) deallocate(ip1listin,stat=istat)
      if (associated(ip1listout)) deallocate(ip1listout,stat=istat)
      istat = vgd_free(vgridin)
      istat = vgd_free(vgridout)
      !------------------------------------------------------------------
      return
   end function samevgrid_str

end module samevgrid_mod
