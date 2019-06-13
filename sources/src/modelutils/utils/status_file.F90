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
subroutine open_status_file3(fn)
   implicit none
!!!#include <arch_specific.hf>
   character(len=*),intent(in) :: fn
   !*@/
   integer,external :: fnom
   logical,save :: init_L = .false.
   integer :: ier,iun
   common /stfile/ iun
   if (.not.init_L) iun = 0
   init_L = .true.
   if (iun > 0) call close_status_file3()
   iun = 0
   ier = fnom(iun, fn, 'FMT',0)
   if (iun < 1 .or. ier < 0) iun = 0
   return
end subroutine open_status_file3


!/@*
subroutine write_status_file3(mesg)
   implicit none
!!!#include <arch_specific.hf>
   character(len=*),intent(in) :: mesg
   !*@/
   integer :: iun
   common /stfile/ iun
   if (iun > 0) then
      write(iun,'(a,";")') trim(mesg)
      call flush(iun)
   endif
   return
end subroutine write_status_file3


!/@*
subroutine close_status_file3()
   implicit none
!!!#include <arch_specific.hf>
   !*@/
   integer :: iun
   common /stfile/ iun
   if (iun > 0) call fclos(iun)
   return
end subroutine close_status_file3
