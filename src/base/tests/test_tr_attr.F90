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
subroutine test_tr_attr()
   use testutils
   use tracers_attributes_mod, only: tracers_attributes
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04, 2017-09
!@/
#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   integer, parameter :: default_wload= 0
   integer, parameter :: default_hzd  = 0
   integer, parameter :: default_monot= 1
   integer, parameter :: default_massc= 0
   real   , parameter :: default_vmin = -1*huge(1.)
   character(len=*), parameter :: default_intp_S = 'CUBIC'

   character(len=2048) :: attributes_S, intp_S
   integer :: wload, hzd, monot, massc, istat
   real    :: vmin
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_fst')

   attributes_S = ''
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin)
   call testutils_assert_eq(wload, default_wload, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, default_monot, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  default_vmin, 'default vmin')

   attributes_S = 'wload=1,mono=2,min=2.1'
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')

   attributes_S = 'default,wload=1,mono=2,min=2.1'
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')

   attributes_S = ''
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')

   attributes_S = ''
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin, F_intp_S=intp_S)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')
   call testutils_assert_eq(intp_S, default_intp_S, 'default intp')

   attributes_S = 'default,intp=none'
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin, F_intp_S=intp_S)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')
   call testutils_assert_eq(intp_S, 'NONE', 'default intp')

   attributes_S = ''
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin, F_intp_S=intp_S)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')
   call testutils_assert_eq(intp_S, 'NONE', 'default intp')

   attributes_S = 'toto=1'
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin, F_intp_S=intp_S)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')
   call testutils_assert_eq(intp_S, 'NONE', 'default intp')
   call testutils_assert_eq(istat, RMN_ERR, 'istat not ignore')

   attributes_S = 'toto=1'
   istat = tracers_attributes(attributes_S, wload, hzd, monot, massc, vmin, F_intp_S=intp_S, F_ignore_L=.true.)
   call testutils_assert_eq(wload, 1, 'default wload')
   call testutils_assert_eq(hzd,   default_hzd, 'default hzd')
   call testutils_assert_eq(monot, 2, 'default monot')
   call testutils_assert_eq(massc, default_massc, 'default massc')
   call testutils_assert_eq(vmin,  2.1, 'default vmin')
   call testutils_assert_eq(intp_S, 'NONE', 'default intp')
   call testutils_assert_eq(istat, RMN_OK, 'istat ignore')

   ! ---------------------------------------------------------------------
   return
end subroutine test_tr_attr
