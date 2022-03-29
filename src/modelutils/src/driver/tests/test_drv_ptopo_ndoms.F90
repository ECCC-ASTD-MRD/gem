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

#include <msg.h>

!/@
subroutine test_drv_ptopo_ndoms()
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_putenv
   use testutils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
!@/
#include <rmnlib_basics.hf>
   integer :: ndomains,dom_deb,istat
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_drv_ptopo_ndoms')

   call drv_ptopo_ndoms(ndomains,dom_deb,istat)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_drv_ptopo_ndoms','UM_EXEC_NDOMAINS not def')

   istat = clib_putenv('UM_EXEC_NDOMAINS=')
   call drv_ptopo_ndoms(ndomains,dom_deb,istat)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_drv_ptopo_ndoms','UM_EXEC_NDOMAINS empty')

   istat = clib_putenv('UM_EXEC_NDOMAINS=1:')
   call drv_ptopo_ndoms(ndomains,dom_deb,istat)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_drv_ptopo_ndoms','UM_EXEC_NDOMAINS invalid')

   istat = clib_putenv('UM_EXEC_NDOMAINS=2')
   call drv_ptopo_ndoms(ndomains,dom_deb,istat)
   call testutils_assert_eq((/ndomains,dom_deb,istat/),(/1,2,0/),'UM_EXEC_NDOMAINS=2')

   istat = clib_putenv('UM_EXEC_NDOMAINS=5:8')
   call drv_ptopo_ndoms(ndomains,dom_deb,istat)
   call testutils_assert_eq((/ndomains,dom_deb,istat/),(/4,5,0/),'UM_EXEC_NDOMAINS=5:8')
   ! ---------------------------------------------------------------------
   return
end subroutine test_drv_ptopo_ndoms
