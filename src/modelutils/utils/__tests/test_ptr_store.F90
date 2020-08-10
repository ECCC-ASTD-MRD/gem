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

!/@
subroutine test_ptr_store()
   use clib_itf_mod
   use testutils
   use ptr_store
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
   !@/
#include <rmnlib_basics.hf>
   real,pointer,dimension(:,:,:) :: p1,p2,p3,p1b,p1c
   integer :: l_ijk(3),u_ijk(3),l_ijk2(3),u_ijk2(3)
   logical :: ok_L
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_ptr_store')

   l_ijk(:) = (/3,2,1/)
   u_ijk(:) = l_ijk(:) + (/2,3,4/)

   call ptr_store_get(p1,l_ijk,u_ijk)
   ok_l = (associated(p1) .and. &
        all(lbound(p1) == l_ijk) .and. &
        all(ubound(p1) == u_ijk))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p1')
   p1 = 1.

   call ptr_store_get(p2,l_ijk,u_ijk)
   ok_l = (associated(p2) .and. &
        all(lbound(p2) == l_ijk) .and. &
        all(ubound(p2) == u_ijk))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p2')
   p2 = 2.
   ok_l = (.not.associated(p2,p1))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p2 != p1')
   ok_l = (all(p2 /= p1))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p2 != p1')

   l_ijk2 = 1
   u_ijk2 = 9
   call ptr_store_get(p3,l_ijk2,u_ijk2)
   ok_l = (associated(p3) .and. &
        all(lbound(p3) == l_ijk2) .and. &
        all(ubound(p3) == u_ijk2))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p3')
   p3 = 3.

   p1b => p1
   call ptr_store_free(p1)
   ok_l = (.not.associated(p1))
   call testutils_assert_ok(ok_L,'test_ptr_store','free p1')

   call ptr_store_get(p1c,l_ijk,u_ijk)
   ok_l = (associated(p1c,p1b))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p2 != p1')
   ok_l = (all(p1c == 1.))
   call testutils_assert_ok(ok_L,'test_ptr_store','get p2 != p1')
   ! ---------------------------------------------------------------------
   return
end subroutine test_ptr_store
