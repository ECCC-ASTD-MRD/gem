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

#include <rmn/msg.h>

!/@
subroutine test_fstmpio()
   use iso_c_binding
   use testutils
   use ptopo_utils
   use fstmpio_mod, only: fstmpio_set_iotype
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
!@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   integer :: istat, myproc, ndomains, idomain, ngrids, igrid
   character(len=512) :: dfiles_S, bcmk_S
   ! ---------------------------------------------------------------------
!!$   print *,'---------------- testutils_initmpi'
   myproc = testutils_initmpi()
!!$   print *,testutils_myproc,'---------------- ptopo_init_var'
   ndomains = 1
 
   call ptopo_init_var(ndomains, idomain, ngrids, igrid)

!!$   print *,testutils_myproc,ptopo_grid_ipe,'---------------- ptopo_io_set',testutils_npeio
   istat = ptopo_io_set(testutils_npeio)
!!$   print *,testutils_myproc,ptopo_grid_ipe,'---------------- '

   istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
   call msg_set_p0only(0)

   call fstmpio_set_iotype(PTOPO_IO)

   istat = clib_getenv('ATM_MODEL_DFILES', dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      print *, 'ERROR: ATM_MODEL_DFILES not defined'
      return
   endif
   bcmk_S = trim(dfiles_S)//'/bcmk/'

   call testutils_set_name('test_fstmpio_open_notfound')
   call test_fstmpio_open_notfound(bcmk_S)

   call testutils_set_name('test_fstmpio_find_0')
   call test_fstmpio_find_0(bcmk_S)

   call testutils_set_name('test_fstmpio_find_3d_0')
   call test_fstmpio_find_3d_0(bcmk_S)

   call testutils_set_name('test_fstmpio_find_3d_0b')
   call test_fstmpio_find_3d_0b(bcmk_S)

   call testutils_set_name('test_fstmpio_find_3d_vect')
   call test_fstmpio_find_3d_vect(bcmk_S)

   call testutils_set_name('test_fstmpio_find_3d_vectb')
   call test_fstmpio_find_3d_vectb(bcmk_S)

   call testutils_set_name('test_fstmpio_rdhint')
   call test_fstmpio_rdhint(bcmk_S)

   call testutils_set_name('test_fstmpio_rdhint_vect')
   call test_fstmpio_rdhint_vect(bcmk_S)

!!$   call testutils_set_name('test_fstmpio_meta')
!!$   call test_fstmpio_meta(bcmk_S)
!!$   call testutils_set_name('test_fstmpio_vgrid')
!!$   call test_fstmpio_vgrid(bcmk_S)
!!$   call testutils_set_name('test_fstmpio_read')
!!$   call test_fstmpio_read(bcmk_S)

  !!$   call testutils_set_name('test_fstmpio_find_fuzz')
  !!$   call test_fstmpio_find_fuzz(bcmk_S)

   call testutils_set_name('test_fstmpio_write')
   !# call test_fstmpio_write()

   call testutils_stats_print()
   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio


!/@
subroutine test_fstmpio_open_notfound(F_bcmk_S)
   use testutils
   use fstmpio_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, istat
   ! ---------------------------------------------------------------------
   funit = fstmpio_open(trim(F_bcmk_S)//'__does_not_exists__', FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .false., '')
   if (funit > 0) istat = fstmpio_close(funit)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_open_notfound


!/@
subroutine test_fstmpio_find_0(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat, datev, datev2, datev3, key, datevfuzz, kind, ip1
   character(len=512) :: filename_S, datev_S, dummy_S
   real :: zp1
   real(REAL64) :: nhours_8
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   istat = fstmpio_find_0(key, funit, 'NONE', datev, RMN_ANY_I, 0, 0)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'find not found')

   istat = fstmpio_find_0(key, funit, 'TT', datev, RMN_ANY_I, 0, 0)
   call testutils_assert_eq(RMN_IS_OK(key), .true., 'find')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz)
   call testutils_assert_ok(.not.RMN_IS_OK(key), 'fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_near')
   call testutils_assert_ok(datev2==datev, 'fuzz_near value')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LE)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_le')
   call testutils_assert_ok(datev2==datev, 'fuzz_le value')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GE)
   call testutils_assert_ok(.not.RMN_IS_OK(key), 'fuzz_ge not found')

   datev_S = '20090426.220000'
   call datp2f(datev2, datev_S)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GE)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_ge')
   call testutils_assert_ok(datev2==datev, 'fuzz_ge value')


   datev_S = '20090427.000000'
   call datp2f(datev2, datev_S)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GT)
   call testutils_assert_ok(.not.RMN_IS_OK(key), 'fuzz_gt not found')

   datev_S = '20090427.000000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600*6
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LT)
   call testutils_assert_ok(.not.RMN_IS_OK(key), 'fuzz_lt not found')

   datev_S = '20090427.000000'
   call datp2f(datev3, datev_S)
   nhours_8 = -40.D0/SEC_PER_HR
   call incdatr(datev2, datev3, nhours_8)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GT)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_gt as eq')
   call testutils_assert_ok(datev2==datev, 'fuzz_gt as eq value')
   call datf2p(datev_S, datev2)


   datev_S = '20090427.000000'
   call datp2f(datev3, datev_S)
   nhours_8 = -1.D0
   call incdatr(datev2, datev3, nhours_8)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GT)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_gt')
   call testutils_assert_ok(datev2==datev, 'fuzz_gt value')


   datev_S = '20090427.000000' 
   call datp2f(datev3, datev_S)
   nhours_8 = 40.D0/SEC_PER_HR
   call incdatr(datev2, datev3, nhours_8)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LT)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_lt as eq')
   call testutils_assert_ok(datev2==datev, 'fuzz_lt as eq value')


   datev_S = '20090427.000000' 
   call datp2f(datev3, datev_S)
   nhours_8 = 1.D0
   call incdatr(datev2, datev3, nhours_8)
   istat = fstmpio_find_0(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LT)
   call testutils_assert_ok(RMN_IS_OK(key), 'fuzz_lt')
   call testutils_assert_ok(datev2==datev, 'fuzz_lt value')
   call datf2p(datev_S, datev2)

   istat = fstmpio_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat), 'close')


   filename_S = trim(F_bcmk_S)//'geophy/Gem_geophy.fst'
   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_ok(RMN_IS_OK(funit), 'open')

   datev = RMN_ANY_DATE
   zp1 = 1.
   kind = RMN_CONV_ARBITRARY
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPOLD, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   istat = fstmpio_find_0(key, funit, 'J1', datev, ip1, RMN_ANY_I, RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key), 'find ip1>0 old')

   datev = RMN_ANY_DATE
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPNEW, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   istat = fstmpio_find_0(key, funit, 'J1', datev, ip1, RMN_ANY_I, RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key), 'find ip1>0 new')

   datev = RMN_ANY_DATE
   ip1 = 1200
   istat = fstmpio_find_0(key, funit, 'ME', datev, ip1, RMN_ANY_I, RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key), 'find ip1=1200')

   datev = RMN_ANY_DATE
   ip1 = 0
   istat = fstmpio_find_0(key, funit, 'ME', datev, ip1, RMN_ANY_I, RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key), 'find ip1=0 for 1200')

   datev = RMN_ANY_DATE
   ip1 = 0
   istat = fstmpio_find_0(key, funit, 'MG', datev, ip1, RMN_ANY_I, RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key), 'find ip1=0')

   datev = RMN_ANY_DATE
   ip1 = 1200
   istat = fstmpio_find_0(key, funit, 'MG', datev, ip1, RMN_ANY_I, RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key), 'find ip1=1200 for 0')

   istat = fstmpio_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat), 'close')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_find_0

!/@
subroutine test_fstmpio_find_3d_0(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
    use vGrid_Descriptors
  implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, datev, datev2, datevfuzz, nkeys
   character(len=512) :: filename_S, datev_S
   integer, pointer :: keys1(:), ip1s(:), ip1s2(:)
   type(vgrid_descriptor) :: vgrid1
   character(len=12) :: lvltyp_S
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'NONE', datev, ip1s, 0, 0)
   call testutils_assert_ok(nkeys <= 0, 'find not found')

   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find nip1s')
   call testutils_assert_eq(keys1(1), 1025, 'find keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 81921, 'find keys1(nkeys)')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')

   nullify(keys1)
   ip1s2 => ip1s(1:5)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find 5 nkeys1')
   call testutils_assert_eq(size(ip1s2), nkeys, 'find 5 nip1s')
   call testutils_assert_eq(keys1(1), 1025, 'find 5 keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 5121, 'find 5 keys1(nkeys)')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0, &
        F_vgrid=vgrid1, F_lvltyp_S=lvltyp_S)
   call testutils_assert_eq(nkeys, 80, 'find v nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find v associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find v nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find v nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find v ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find v ip1s(nkeys)')
!!$   print *,lvltyp_S
!!$   print *,vgrid1

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_find_3d_0

!/@
subroutine test_fstmpio_find_3d_0b(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   use vGrid_Descriptors
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer, parameter :: NMAXKEYS = 1000
   integer :: funit, datev, datev2, datevfuzz, nkeys
   character(len=512) :: filename_S, datev_S
   integer, dimension(NMAXKEYS), target :: tkeys1, tkeys2, tip1s, tip1s2
   integer, pointer :: keys1(:), keys2(:), ip1s(:), ip1s2(:)
   type(vgrid_descriptor) :: vgrid1
   character(len=12) :: lvltyp_S
   ! ---------------------------------------------------------------------
   keys1 => tkeys1
   keys2 => tkeys2
   ip1s => tip1s
   ip1s2 => tip1s2

   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

!!$   nullify(keys1, ip1s)
   keys1 = -1
   ip1s = -1
   nkeys = fstmpio_find_3d_0(keys1, funit, 'NONE', datev, ip1s, 0, 0)
   call testutils_assert_ok(nkeys <= 0, 'find not found')

!!$   nullify(keys1, ip1s)
   keys1 = -1
   ip1s = -1
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find associated keys1')
   call testutils_assert_eq(associated(ip1s, tip1s), .true., 'find associated ip1s')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), NMAXKEYS, 'find nip1s')
   call testutils_assert_eq(keys1(1), 1026, 'find keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 81922, 'find keys1(nkeys)')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')

!!$   nullify(keys1)
   keys1 = -1
   ip1s2 => ip1s(1:5)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find 5 nkeys1')
   call testutils_assert_eq(size(ip1s2), nkeys, 'find 5 nip1s')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find 5 associated keys1')
!!$   call testutils_assert_eq(associated(ip1s2, ip1s), .true., 'find 5 associated ip1s')
   call testutils_assert_eq(keys1(1), 1026, 'find 5 keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 5122, 'find 5 keys1(nkeys)')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

!!$   nullify(keys1, ip1s)
   keys1 = -1
   ip1s = -1
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0, &
        F_vgrid=vgrid1, F_lvltyp_S=lvltyp_S)
   call testutils_assert_eq(nkeys, 80, 'find v nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find v associated')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find v associated keys')
   call testutils_assert_eq(associated(ip1s, tip1s), .true., 'find v associated ip1s')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find v nkeys1')
   call testutils_assert_eq(size(ip1s), NMAXKEYS, 'find v nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find v ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find v ip1s(nkeys)')
!!$   print *,lvltyp_S
!!$   print *,vgrid1

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_find_3d_0b


subroutine test_fstmpio_find_3d_vect(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   use vGrid_Descriptors
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, datev, datev2, datevfuzz, nkeys, k
   character(len=512) :: filename_S, datev_S
   integer, pointer :: keys1(:), keys2(:), ip1s(:), ip1s2(:)
   type(vgrid_descriptor) :: vgrid1
   character(len=12) :: lvltyp_S
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   nullify(keys1, keys2, ip1s)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(keys2), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')
   do k=1,nkeys
      call testutils_assert_neq(keys1(k), keys2(k), 'diff keys')
   enddo

   nullify(keys1, keys2)
   ip1s2 => ip1s(1:5)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find 5 nkeys1')
   call testutils_assert_eq(size(keys2), nkeys, 'find 5 nkeys1')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

   nullify(keys1, keys2, ip1s)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, &
        ip1s, 0, 0, F_vgrid=vgrid1, F_lvltyp_S=lvltyp_S)
   call testutils_assert_eq(nkeys, 80, 'find v nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find v associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find v nkeys1')
   call testutils_assert_eq(size(keys2), nkeys, 'find v nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find v nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find v ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find v ip1s(nkeys)')
   do k=1,nkeys
      call testutils_assert_neq(keys1(k), keys2(k), 'diff v keys')
   enddo
!!$   print *,lvltyp_S
!!$   print *,vgrid1


   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, keys2, ip1s)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'find fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   nullify(keys1, keys2, ip1s)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'find fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'find fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_find_3d_vect


subroutine test_fstmpio_find_3d_vectb(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   use vGrid_Descriptors
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer, parameter :: NMAXKEYS = 1000
   integer :: funit, datev, nkeys, k
   character(len=512) :: filename_S, datev_S
   integer, dimension(NMAXKEYS), target :: tkeys1, tkeys2, tip1s, tip1s2
   integer, pointer :: keys1(:), keys2(:), ip1s(:), ip1s2(:)
   type(vgrid_descriptor) :: vgrid1
   character(len=12) :: lvltyp_S
   ! ---------------------------------------------------------------------
   keys1 => tkeys1
   keys2 => tkeys2
   ip1s => tip1s
   ip1s2 => tip1s2

   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

!!$   nullify(keys1, keys2, ip1s)
   keys1 = -1
   keys2 = -1
   ip1s = -1
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find associated keys1')
   call testutils_assert_eq(associated(keys2, tkeys2), .true., 'find associated keys2')
   call testutils_assert_eq(associated(ip1s, tip1s), .true., 'find associated ip1s')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find nkeys1')
   call testutils_assert_eq(size(keys2), NMAXKEYS, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), NMAXKEYS, 'find nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')
   do k=1,nkeys
      call testutils_assert_neq(keys1(k), keys2(k), 'diff keys')
   enddo

!!$   nullify(keys1, keys2)
   keys1 = -1
   keys2 = -1
   ip1s2 => ip1s(1:5)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find 5 associated keys1')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find 5 nkeys1')
   call testutils_assert_eq(size(keys2), NMAXKEYS, 'find 5 nkeys1')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

!!$   nullify(keys1, keys2, ip1s)
   keys1 = -1
   keys2 = -1
   ip1s = -1
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, &
        ip1s, 0, 0, F_vgrid=vgrid1, F_lvltyp_S=lvltyp_S)
   call testutils_assert_eq(nkeys, 80, 'find v nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find v associated')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find v nkeys1')
   call testutils_assert_eq(size(keys2), NMAXKEYS, 'find v nkeys1')
   call testutils_assert_eq(size(ip1s), NMAXKEYS, 'find v nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find v ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find v ip1s(nkeys)')
   do k=1,nkeys
      call testutils_assert_neq(keys1(k), keys2(k), 'diff v keys')
   enddo
!!$   print *,lvltyp_S
!!$   print *,vgrid1


   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_find_3d_vectb


!/@
subroutine test_fstmpio_meta(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, istat, datev, key, &
        dateo, deet, npas, ni, nj, nk, ip1, ip2, ip3, ig1, ig2, ig3, ig4, &
        nbits, datyp
   character(len=512) :: filename_S, datev_S, nomvar_S, etiket_S, typvar_S
   character(len=2) :: grtyp_S
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   istat = fstmpio_find(key, funit, 'TT', datev, RMN_ANY_I, 0, 0)
   istat = fstmpio_getmeta(key, nomvar_S, dateo, deet, npas, &
        ip1, ip2, ip3, etiket_S, typvar_S, ni, nj, nk, nbits, datyp, &
        grtyp_S, ig1, ig2, ig3, ig4)
!!$   print *,trim(nomvar_S), dateo, deet, npas, &
!!$        ip1, ip2, ip3, trim(etiket_S), trim(typvar_S), ni, nj, nk, nbits, datyp, &
!!$        trim(grtyp_S), ig1, ig2, ig3, ig4
   call testutils_assert_eq(nomvar_S, 'TT', 'meta nomvar_S')
   call testutils_assert_eq(etiket_S, 'G133K80P', 'meta etiket_S')
   call testutils_assert_eq(typvar_S, 'P', 'meta typvar_S')
   call testutils_assert_eq(grtyp_S, 'G', 'meta grtyp_S')
   call testutils_assert_eq((/ni, nj, nk/), (/200, 100, 1/), 'meta nijk')
   call testutils_assert_eq((/ip1, ip2, ip3/), (/97642568, 0, 0/), 'meta ip123')
   call testutils_assert_eq((/nbits, datyp/), (/12, 1/), 'meta nbits, datyp')
   call testutils_assert_eq((/ig1, ig2, ig3, ig4/), (/0, 0, 0, 0/), 'meta ig1234')
   call testutils_assert_eq((/dateo, deet, npas/), (/354514400, 900, 0/), 'meta dateo, deet, npas')

   istat = fstmpio_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat), 'close')
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_meta

!/@
subroutine test_fstmpio_vgrid(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   use vGrid_Descriptors
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, istat, datev, key, &
        dateo, deet, npas, ni, nj, nk, ip1, ip2, ip3, ig1, ig2, ig3, ig4, &
        nbits, datyp
   character(len=512) :: filename_S, datev_S, nomvar_S, etiket_S, typvar_S, lvltyp_S,sfcfld_S
   character(len=2) :: grtyp_S
   type(vgrid_descriptor) :: vgrid
   integer, pointer :: ip1list(:)
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   istat = fstmpio_find(key, funit, 'TT', datev, RMN_ANY_I, 0, 0)
   istat = fstmpio_getmeta(key, nomvar_S, dateo, deet, npas, &
        ip1, ip2, ip3, etiket_S, typvar_S, ni, nj, nk, nbits, datyp, &
        grtyp_S, ig1, ig2, ig3, ig4)

   nullify(ip1list)
   lvltyp_S = ' '
   istat = fstmpio_get_vgrid(funit, key, vgrid, ip1list, lvltyp_S)
   call testutils_assert_eq(size(ip1list), 158, 'meta size(ip1list)')
   call testutils_assert_eq(ip1list(1:3), (/97642568, 97690568, 97738568/), 'meta ip1list')
   call testutils_assert_eq('M', lvltyp_S, 'meta lvltyp_S')
   istat = vgd_get(vgrid, key='RFLD', value=sfcfld_S)
   call testutils_assert_eq('P0', sfcfld_S, 'meta sfcfld_S')

   istat = fstmpio_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat), 'close')
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_vgrid

!/@
subroutine test_fstmpio_read(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpio_mod
   use vGrid_Descriptors
   use hinterp4yy_mod
   use ptopo_utils
   use iso_c_binding
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"

   integer, parameter :: GNI = 80
   integer, parameter :: GNJ = 40
   integer, parameter :: HALO = 0
   logical, parameter :: ALONGX = .true.
   logical, parameter :: FILL = .false. !.true.
   integer :: funit(1), istat, datev, key, &
        dateo, deet, npas, ni, nj, nk, ip1, ip2, ip3, ig1, ig2, ig3, ig4, &
        nbits, datyp, outgridid, coregridid, rpncomm_gridid, k, &
        mini,maxi,minj,maxj,lni, lnimax, li0, lnj, lnjmax, lj0
   real, pointer :: data(:, :, :)
   character(len=512) :: filename_S, datev_S, nomvar_S, etiket_S, typvar_S, lvltyp_S
   character(len=2) :: grtyp_S
   type(vgrid_descriptor) :: vgrid
   integer, pointer :: ip1list(:), keylist(:), pkeylist(:), istatlist(:)
   character(len=8), pointer :: hintlist_S(:), phintlist_S(:)
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit(1) = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit(1)), .true., 'open')

   istat = fstmpio_find(key, funit(1), 'TT', datev, RMN_ANY_I, 0, 0)
   istat = fstmpio_getmeta(key, nomvar_S, dateo, deet, npas, &
        ip1, ip2, ip3, etiket_S, typvar_S, ni, nj, nk, nbits, datyp, &
        grtyp_S, ig1, ig2, ig3, ig4)

   nullify(ip1list)
   lvltyp_S = ' '
   istat = fstmpio_get_vgrid(funit(1), key, vgrid, ip1list, lvltyp_S)

   nk = size(ip1list)
   allocate(keylist(nk), hintlist_S(nk))

   nk = 0
   do k = 1, size(ip1list)
      istat = fstmpio_find(key, funit(1), 'TT', datev, ip1list(k), ip2, 0)
      if (RMN_IS_OK(key)) then
         nk = nk+1
         keylist(nk) = key
         hintlist_S(nk) = HINTERP4YY_CUBIC
      endif
   enddo
   istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
        ALONGX, FILL)
   istat = rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
        .not.ALONGX, FILL)
   rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)

   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   nullify(data)
   pkeylist => keylist(1:nk)
   phintlist_S => hintlist_S(1:nk)
   allocate(istatlist(1:nk))
!!$   print *,'dims0:',GNI, GNJ, mini, maxi, minj, maxj
   istat = fstmpio_rdhint(data, istatlist, pkeylist, phintlist_S, &
        funit, rpncomm_gridid, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'fstmpio_rdhint - status')
   call testutils_assert_eq(associated(data), .true., 'fstmpio_rdhint - alloc')
   call testutils_assert_eq(ubound(data,1), maxi, 'fstmpio_rdhint - ni')
   call testutils_assert_eq(ubound(data,2), maxj, 'fstmpio_rdhint - nj')
   call testutils_assert_eq(ubound(data,3), nk, 'fstmpio_rdhint - nk')
!!$   print *,'asso=',associated(data),': nk0=',size(ip1list), nk
!!$   if (associated(data)) then
!!$      print *,'lijk=',lbound(data), ': uijk=', ubound(data)
!!$      print *,'min=',minval(data), ': max=', maxval(data)
!!$   endif
!!$   call flush(6)
   istat = fstmpio_close(funit(1))
   call testutils_assert_ok(RMN_IS_OK(istat), 'close')
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_read

!!$!/@
!!$subroutine test_fstmpio_find_read(F_bcmk_S)
!!$   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
!!$   use testutils
!!$   use fstmpio_mod
!!$   implicit none
!!$   !@objective 
!!$   !@author Stephane Chamberland, 2017-04
!!$   !@argument
!!$   character(len=*), intent(in) :: F_bcmk_S
!!$!@/
!!$#include <rmnlib_basics.hf>
!!$   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
!!$   integer :: funit, istat, datev, datev2, datev3, key, gridid, datevfuzz, &
!!$        ni, nj, ig1, ig2, ig3, ig4, ip1, kind
!!$   real, pointer :: data(:, :, :)
!!$   character(len=512) :: filename_S, datev_S, dummy_S
!!$   character(len=2) :: grtyp_S
!!$   real :: zp1
!!$   real(REAL64) :: nhours_8
!!$   ! ---------------------------------------------------------------------
!!$   datev_S = '20090427.000000'
!!$   call datp2f(datev, datev_S)
!!$   filename_S = trim(F_bcmk_S)//'2009042700_000'
!!$
!!$   funit = fstmpio_open(filename_S, FST_READONLY)
!!$   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')
!!$
!!$   istat = fstmpio_find(key, funit, 'NONE', datev, RMN_ANY_I, 0, 0)
!!$   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'find not found')
!!$
!!$   istat = fstmpio_find(key, funit, 'TT', datev, RMN_ANY_I, 0, 0)
!!$   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'find')
!!$
!!$   istat = fstmpio_read(key, data, funit)
!!$   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'read status')
!!$   istat = RMN_ERR
!!$   if (associated(data) .and. &
!!$        size(data, 1) == 200 .and. &
!!$        size(data, 2) == 100 .and. &
!!$        size(data, 3) == 1  .and. &
!!$        nint(minval(data)*100.) == -5089 .and. & 
!!$        nint(maxval(data)*100.) == -2618) istat = RMN_OK
!!$   call testutils_assert_eq(istat, RMN_OK, 'read values')
!!$   if (associated(data)) deallocate(data, stat=istat)
!!$
!!$   gridid = fstmpio_get_gridid(funit, key)
!!$   call testutils_assert_eq(RMN_IS_OK(gridid), .true., 'get_gridid status')
!!$   istat = ezgprm(gridid, grtyp_S, ni, nj, ig1, ig2, ig3, ig4)
!!$   istat = RMN_ERR
!!$   call testutils_assert_eq(grtyp_S(1:1), 'G', 'get_gridid grtyp')
!!$   call testutils_assert_eq((/ni, nj/), (/200, 100/), 'get_gridid nij')
!!$   call testutils_assert_eq((/ig1, ig2, ig3, ig4/), (/0, 0, 0, 0/), 'get_gridid ig14')
!!$ !!$   ier = ezgxprm(gdid, ni, nj, grtyp, ig1, ig2, ig3, ig4, grref, ig1ref, ig2ref, ig3ref, ig4ref)
!!$ !!$   ier = ezgfstp(gdid, nomvarx, typvarx, etikx, nomvary, typvary, etiky, ip1, ip2, ip3, dateo, deet, npas, nbits)
!!$ !!$   ier = gdgaxes(gdid, ax, ay)
!!$   istat = gdrls(gridid)
!!$
!!$   datev_S = '20090427.020000'
!!$   call datp2f(datev2, datev_S)
!!$   datevfuzz = 3600
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz)
!!$   call testutils_assert_ok(.not.RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_near not found')
!!$
!!$   datevfuzz = 3600*6
!!$
!!$   datev_S = '20090427.020000'
!!$   call datp2f(datev2, datev_S)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_near')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_near value')
!!$
!!$   datev_S = '20090427.020000'
!!$   call datp2f(datev2, datev_S)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LE)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_le')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_le value')
!!$
!!$   datev_S = '20090427.020000'
!!$   call datp2f(datev2, datev_S)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GE)
!!$   call testutils_assert_ok(.not.RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_ge not found')
!!$
!!$   datev_S = '20090426.220000'
!!$   call datp2f(datev2, datev_S)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GE)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_ge')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_ge value')
!!$
!!$
!!$   datev_S = '20090427.000000'
!!$   call datp2f(datev2, datev_S)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GT)
!!$   call testutils_assert_ok(.not.RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_gt not found')
!!$
!!$   datev_S = '20090427.000000'
!!$   call datp2f(datev2, datev_S)
!!$   datevfuzz = 3600*6
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LT)
!!$   call testutils_assert_ok(.not.RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_lt not found')
!!$
!!$   datev_S = '20090427.000000'
!!$   call datp2f(datev3, datev_S)
!!$   nhours_8 = -40.D0/SEC_PER_HR
!!$   call incdatr(datev2, datev3, nhours_8)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GT)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_gt as eq')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_gt as eq value')
!!$   call datf2p(datev_S, datev2)
!!$
!!$
!!$   datev_S = '20090427.000000'
!!$   call datp2f(datev3, datev_S)
!!$   nhours_8 = -1.D0
!!$   call incdatr(datev2, datev3, nhours_8)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_GT)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_gt')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_gt value')
!!$
!!$
!!$   datev_S = '20090427.000000' 
!!$   call datp2f(datev3, datev_S)
!!$   nhours_8 = 40.D0/SEC_PER_HR
!!$   call incdatr(datev2, datev3, nhours_8)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LT)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_lt as eq')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_lt as eq value')
!!$
!!$
!!$   datev_S = '20090427.000000' 
!!$   call datp2f(datev3, datev_S)
!!$   nhours_8 = 1.D0
!!$   call incdatr(datev2, datev3, nhours_8)
!!$   istat = fstmpio_find(key, funit, 'TT', datev2, RMN_ANY_I, 0, 0, datevfuzz, FST_FIND_LT)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fuzz_lt')
!!$   call testutils_assert_ok(datev2==datev, 'test_fstmpio_find_read', 'fuzz_lt value')
!!$   call datf2p(datev_S, datev2)
!!$
!!$   istat = fstmpio_close(funit)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_find_read', 'fstmpio_close')
!!$
!!$
!!$   filename_S = trim(F_bcmk_S)//'geophy/Gem_geophy.fst'
!!$   funit = fstmpio_open(filename_S, FST_READONLY)
!!$   call testutils_assert_ok(RMN_IS_OK(funit), 'test_fstmpio_find_read2', 'fstmpio_open')
!!$
!!$   datev = RMN_ANY_DATE
!!$   zp1 = 1.
!!$   kind = RMN_CONV_ARBITRARY
!!$   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPOLD, dummy_S, .not.RMN_CONV_USEFORMAT_L)
!!$   istat = fstmpio_find(key, funit, 'J1', datev, ip1, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fstmpio_find ip1>0 old')
!!$
!!$   datev = RMN_ANY_DATE
!!$   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPNEW, dummy_S, .not.RMN_CONV_USEFORMAT_L)
!!$   istat = fstmpio_find(key, funit, 'J1', datev, ip1, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fstmpio_find ip1>0 new')
!!$
!!$   datev = RMN_ANY_DATE
!!$   ip1 = 1200
!!$   istat = fstmpio_find(key, funit, 'ME', datev, ip1, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fstmpio_find ip1=1200')
!!$
!!$   datev = RMN_ANY_DATE
!!$   ip1 = 0
!!$   istat = fstmpio_find(key, funit, 'ME', datev, ip1, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fstmpio_find ip1=0 for 1200')
!!$
!!$   datev = RMN_ANY_DATE
!!$   ip1 = 0
!!$   istat = fstmpio_find(key, funit, 'MG', datev, ip1, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fstmpio_find ip1=0')
!!$
!!$   datev = RMN_ANY_DATE
!!$   ip1 = 1200
!!$   istat = fstmpio_find(key, funit, 'MG', datev, ip1, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_find_read', 'fstmpio_find ip1=1200 for 0')
!!$
!!$   istat = fstmpio_close(funit)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_find_read2', 'fstmpio_close')
!!$
!!$   ! ---------------------------------------------------------------------
!!$   return
!!$end subroutine test_fstmpio_find_read
!!$
!!$
!!$!/@
!!$subroutine test_fstmpio_write()
!!$   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
!!$   use testutils
!!$   use fstmpio_mod
!!$   use ezgrid_mod
!!$   use ptopo_utils
!!$   implicit none
!!$   !@objective 
!!$   !@author Stephane Chamberland, 2012-01
!!$   !@argument
!!$!@/
!!$#include <rmnlib_basics.hf>
!!$#include <clib_interface_mu.hf>
!!$   include "rpn_comm.inc"
!!$   real, parameter :: MYVALUE = 3.3
!!$   integer, parameter :: NI0=50, NJ0=30, NK0=3, HALO=2
!!$   character(len=256) :: nomvar_S, filename_S
!!$   logical :: ok_L, ok2_L
!!$   integer :: funit, istat, gridid, grididh, grididfull, gridid2, lvlid, dateo, deet, npas, key, ig1, ig2, ig3, ig4, ip1, i, j, k, datev, ip1list(NK0)
!!$   real, pointer :: data2d(:, :), data2dh(:, :), data3d(:, :, :), data3dh(:, :, :)
!!$   real :: ax(1-HALO:NI0+HALO, 1), ay(1, 1-HALO:NJ0+HALO)
!!$   ! ---------------------------------------------------------------------
!!$   call ptopo_init_var()
!!$   write(filename_S, '(a, I3.3)') '__test_fstmpio_to-rm__.fst-', ptopo_grid_ipe
!!$   istat = clib_unlink(trim(filename_S))
!!$   funit = fstmpio_open(filename_S)
!!$   call testutils_assert_ok(RMN_IS_OK(funit), 'test_fstmpio_write:open', '')
!!$
!!$   nomvar_S = 'ZX'
!!$
!!$   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
!!$   do i=1-HALO, NI0+HALO
!!$      ax(i, 1) = 10.+float(ptopo_bloc_ipex*NI0+i)*0.25
!!$   enddo
!!$   do j=1-HALO, NJ0+HALO
!!$      ay(1, j) = float(ptopo_bloc_ipey*NJ0+j)*0.25
!!$   enddo
!!$   grididh = ezgdef_fmem(NI0+2*HALO, NJ0+2*HALO, 'Z', 'E', ig1, ig2, ig3, ig4, ax, ay)
!!$   gridid = ezgdef_fmem(NI0, NJ0, 'Z', 'E', ig1, ig2, ig3, ig4, ax(1:NI0, 1), ay(1, 1:NJ0))
!!$   grididfull = ezgrid_merge(gridid, RPN_COMM_BLOC_COMM, .true.)
!!$
!!$   allocate(data2d(NI0, NJ0), &
!!$        data2dh(1-HALO:NI0+HALO, 1-HALO:NJ0+HALO), &
!!$        data3d(NI0, NJ0, NK0), &
!!$        data3dh(1-HALO:NI0+HALO, 1-HALO:NJ0+HALO, NK0), &
!!$        stat=istat)
!!$   data2d = MYVALUE
!!$   data2dh = MYVALUE
!!$   data3d = MYVALUE
!!$   data3dh = MYVALUE
!!$   ip1 = 0
!!$ !!$   lvlid = 
!!$   do k=1, NK0
!!$      ip1list(k) = (k+2)*3
!!$   enddo
!!$
!!$   nomvar_S = 'ZX2d'
!!$   istat = fstmpio_write(funit, nomvar_S, data2d, gridid, ip1, F_npak=FST_NPAK_FULL32)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:write_2d_r4', '')
!!$
!!$   nomvar_S = 'ZH2d'
!!$   istat = fstmpio_write(funit, nomvar_S, data2dh, grididh, ip1, F_npak=FST_NPAK_FULL32, F_lni=NI0, F_lnj=NJ0)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:write_2d_r4', '')
!!$
!!$   nomvar_S = 'ZX3d'
!!$   istat = fstmpio_write(funit, nomvar_S, data3d, gridid, ip1list, F_npak=FST_NPAK_FULL32)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:write_3d_r4', '')
!!$
!!$   nomvar_S = 'ZH3d'
!!$   istat = fstmpio_write(funit, nomvar_S, data3dh, grididh, ip1list, F_npak=FST_NPAK_FULL32, F_lni=NI0, F_lnj=NJ0)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:write_3d_r4', '')
!!$   !TODO: test with a vgrid instead of ip1list
!!$
!!$   if (funit > 0) istat = fstmpio_close(funit)
!!$   deallocate(data2d, data2dh, data3d, stat=istat)
!!$
!!$   !- Checking
!!$
!!$   funit = fstmpio_open(filename_S, FST_READONLY)
!!$   call testutils_assert_ok(RMN_IS_OK(funit), 'test_fstmpio_write:open', '')
!!$
!!$   nomvar_S = 'ZX2d'
!!$   datev = RMN_ANY_DATE
!!$   istat = fstmpio_find(key, funit, nomvar_S, datev, RMN_ANY_I, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_write:find', '2d')
!!$
!!$   istat = fstmpio_read(key, data3d, funit, gridid2)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:read', '2d')
!!$   ok_L = .false. ; ok2_L = .false.
!!$   if (associated(data3d)) then
!!$      ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0, ptopo_bloc_npey*NJ0, 1/))
!!$      ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
!!$   endif
!!$   call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data2d shape')
!!$   call testutils_assert_ok(ok2_L, 'test_fstmpio_write:read', 'data2d')
!!$   if (.not.ok2_L) then
!!$      print *, 'data2d min, max:', minval(data3d), maxval(data3d)
!!$   endif
!!$   ok_L = ezgrid_samegrid(grididfull, gridid2)
!!$   call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data2d grid')
!!$
!!$
!!$   if (associated(data3d)) deallocate(data3d, stat=istat)
!!$   nomvar_S = 'ZH2d'
!!$   datev = RMN_ANY_DATE
!!$   istat = fstmpio_find(key, funit, nomvar_S, datev, RMN_ANY_I, RMN_ANY_I, RMN_ANY_I)
!!$   call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_write:find', '2dh')
!!$
!!$   istat = fstmpio_read(key, data3d, funit, gridid2)
!!$   call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:read', '2dh')
!!$   ok_L = .false. ; ok2_L = .false.
!!$   if (associated(data3d)) then
!!$      ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0, ptopo_bloc_npey*NJ0, 1/))
!!$      ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
!!$   endif
!!$   call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data2dh shape')
!!$   call testutils_assert_ok(ok2_L, 'test_fstmpio_write:read', 'data2dh')
!!$   if (.not.ok2_L) then
!!$      print *, 'data2d min, max:', minval(data3d), maxval(data3d)
!!$   endif
!!$ !!$   ok_L = ezgrid_samegrid(grididfull, gridid2)
!!$ !!$   call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data2dh grid')
!!$
!!$
!!$   nomvar_S = 'ZX3d'
!!$   do k=1, NK0
!!$      if (associated(data3d)) deallocate(data3d, stat=istat)
!!$      datev = RMN_ANY_DATE
!!$      istat = fstmpio_find(key, funit, nomvar_S, datev, ip1list(k), RMN_ANY_I, RMN_ANY_I)
!!$      call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_write:find', '3d')
!!$
!!$      istat = fstmpio_read(key, data3d, funit, gridid2)
!!$      call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:read', '3d')
!!$      ok_L = .false. ; ok2_L = .false.
!!$      if (associated(data3d)) then
!!$         ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0, ptopo_bloc_npey*NJ0, 1/))
!!$         ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
!!$      endif
!!$      call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data3d shape')
!!$      call testutils_assert_ok(ok2_L, 'test_fstmpio_write:read', 'data3d')
!!$      if (.not.ok2_L) then
!!$         print *, 'data3d min, max:', minval(data3d), maxval(data3d)
!!$      endif
!!$      ok_L = ezgrid_samegrid(grididfull, gridid2)
!!$      call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data3d grid')
!!$   enddo
!!$
!!$   nomvar_S = 'ZH3d'
!!$   do k=1, NK0
!!$      if (associated(data3d)) deallocate(data3d, stat=istat)
!!$      datev = RMN_ANY_DATE
!!$      istat = fstmpio_find(key, funit, nomvar_S, datev, ip1list(k), RMN_ANY_I, RMN_ANY_I)
!!$      call testutils_assert_ok(RMN_IS_OK(key), 'test_fstmpio_write:find', '3dh')
!!$
!!$      istat = fstmpio_read(key, data3d, funit, gridid2)
!!$      call testutils_assert_ok(RMN_IS_OK(istat), 'test_fstmpio_write:read', '3dh')
!!$      ok_L = .false. ; ok2_L = .false.
!!$      if (associated(data3d)) then
!!$         ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0, ptopo_bloc_npey*NJ0, 1/))
!!$         ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
!!$      endif
!!$      call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data3dh shape')
!!$      call testutils_assert_ok(ok2_L, 'test_fstmpio_write:read', 'data3dh')
!!$      if (.not.ok2_L) then
!!$         print *, 'data3d min, max:', minval(data3d), maxval(data3d)
!!$      endif
!!$      ok_L = ezgrid_samegrid(grididfull, gridid2)
!!$      call testutils_assert_ok(ok_L, 'test_fstmpio_write:read', 'data3d grid')
!!$   enddo
!!$
!!$   if (funit > 0) istat = fstmpio_close(funit)
!!$
!!$   if (associated(data2d)) deallocate(data2d, stat=istat)
!!$   if (associated(data3d)) deallocate(data3d, stat=istat)
!!$
!!$   istat = clib_unlink(trim(filename_S))
!!$   ! ---------------------------------------------------------------------
!!$   return
!!$end subroutine test_fstmpio_write


!/@
subroutine test_fstmpio_rdhint(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use iso_c_binding
   use testutils
   use ptopo_utils
   use fstmpio_mod
   use hinterp4yy_mod
   use statfld_dm_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"

   integer, parameter :: GNI = 80
   integer, parameter :: GNJ = 40
   integer, parameter :: HALO = 0
   logical, parameter :: ALONGX = .true.
   logical, parameter :: FILL = .false. !.true.
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat, datev, nkeys, outgridid, coregridid, &
        k, rpncomm_gridid, &
        mini,maxi,minj,maxj,lni, lnimax, li0, lnj, lnjmax, lj0
   real, pointer :: data1(:, :, :), data2(:, :, :)
   character(len=512) :: filename_S, datev_S, dummy_S, varname_S
   integer, target :: tkeys1(1), tip1s(1), fileid(1), tistats1(1)
   integer, pointer :: keys1(:), ip1s(:), istatlist(:), fidlist(:)
   character(len=8), target :: hintlist_S(1)
   character(len=8), pointer :: phintlist_S(:)

   integer :: ijkmin(3), ijkmax(3), ii, ilvl, vals2(10), allvals2(10,2,80)
   real(REAL64) :: mean, var, rmin, rmax
   logical :: ok_L
   ! ---------------------------------------------------------------------
   allvals2 = 0
   allvals2(:,   1,   1) = (/     967063,      93001,     550814,    1039260,  20,  28,   1,  49,  32,   1/) !P0
   allvals2(:,   2,   1) = (/     -36252,       5162,     -50106,     -26181,   7,   2,   1,  48,  38,   1/) !TT
   allvals2(:,   2,   2) = (/     -29434,       7552,     -51758,      -7742,  65,   9,   1,  51,  38,   1/) !TT
   allvals2(:,   2,   3) = (/     -19400,       9093,     -44699,       3769,  65,   9,   1,  55,  37,   1/) !TT
   allvals2(:,   2,   4) = (/     -11541,       9066,     -38341,       8620,  67,   9,   1,   5,  40,   1/) !TT
   allvals2(:,   2,   5) = (/     -10035,      11368,     -40158,       7943,  64,   7,   1,  14,  38,   1/) !TT
   allvals2(:,   2,   6) = (/     -14894,      14332,     -50057,       8663,   5,   2,   1,  62,  31,   1/) !TT
   allvals2(:,   2,   7) = (/     -22538,      15158,     -58211,       3006,   5,   4,   1,  64,  31,   1/) !TT
   allvals2(:,   2,   8) = (/     -29529,      14051,     -61202,      -6080,   6,   4,   1,  63,  30,   1/) !TT
   allvals2(:,   2,   9) = (/     -34839,      13000,     -64083,     -14432,   8,   4,   1,  62,  29,   1/) !TT
   allvals2(:,   2,  10) = (/     -39446,      12227,     -68657,     -20666,  12,   3,   1,  50,  21,   1/) !TT
   allvals2(:,   2,  11) = (/     -43496,      11389,     -70855,     -26918,   5,   4,   1,  48,  20,   1/) !TT
   allvals2(:,   2,  12) = (/     -46781,      10688,     -73347,     -31783,   2,   4,   1,  71,  30,   1/) !TT
   allvals2(:,   2,  13) = (/     -49298,      10207,     -75399,     -35558,   2,   3,   1,  52,  29,   1/) !TT
   allvals2(:,   2,  14) = (/     -51307,       9695,     -76205,     -39037,   1,   2,   1,   6,  27,   1/) !TT
   allvals2(:,   2,  15) = (/     -52942,       9074,     -76268,     -40794,  77,   2,   1,  34,  33,   1/) !TT
   allvals2(:,   2,  16) = (/     -54197,       8434,     -76630,     -41836,   2,   1,   1,  35,  33,   1/) !TT
   allvals2(:,   2,  17) = (/     -55306,       7769,     -76850,     -43272,  78,   1,   1,  36,  33,   1/) !TT
   allvals2(:,   2,  18) = (/     -56298,       7133,     -77354,     -44534,  79,   2,   1,  36,  33,   1/) !TT
   allvals2(:,   2,  19) = (/     -57087,       6575,     -77176,     -45338,  78,   2,   1,  37,  33,   1/) !TT
   allvals2(:,   2,  20) = (/     -57856,       6132,     -76745,     -46268,  73,   1,   1,  37,  33,   1/) !TT
   allvals2(:,   2,  21) = (/     -58715,       5785,     -76322,     -47493,  75,   1,   1,  37,  33,   1/) !TT
   allvals2(:,   2,  22) = (/     -59456,       5577,     -75673,     -47884,  74,   1,   1,  38,  33,   1/) !TT
   allvals2(:,   2,  23) = (/     -60177,       5603,     -75015,     -48397,  66,   1,   1,  19,  33,   1/) !TT
   allvals2(:,   2,  24) = (/     -60888,       5913,     -74496,     -48276,  67,   1,   1,  64,  35,   1/) !TT
   allvals2(:,   2,  25) = (/     -61603,       6377,     -73779,     -47714,  68,   1,   1,  65,  35,   1/) !TT
   allvals2(:,   2,  26) = (/     -62323,       6974,     -76740,     -47112,  23,  24,   1,  65,  35,   1/) !TT
   allvals2(:,   2,  27) = (/     -62973,       7657,     -78197,     -46899,  25,  23,   1,  66,  35,   1/) !TT
   allvals2(:,   2,  28) = (/     -63569,       8396,     -82581,     -46758,  24,  24,   1,  66,  35,   1/) !TT
   allvals2(:,   2,  29) = (/     -64100,       9091,     -85557,     -46786,  24,  24,   1,  67,  35,   1/) !TT
   allvals2(:,   2,  30) = (/     -64556,       9700,     -85684,     -47099,  24,  24,   1,  67,  35,   1/) !TT
   allvals2(:,   2,  31) = (/     -64883,      10200,     -85898,     -47138,  22,  23,   1,  68,  35,   1/) !TT
   allvals2(:,   2,  32) = (/     -65043,      10552,     -86236,     -47294,  67,  19,   1,  68,  35,   1/) !TT
   allvals2(:,   2,  33) = (/     -65072,      10767,     -86972,     -47524,  67,  19,   1,  68,  35,   1/) !TT
   allvals2(:,   2,  34) = (/     -64955,      10874,     -86955,     -47738,   6,  20,   1,  70,  36,   1/) !TT
   allvals2(:,   2,  35) = (/     -64707,      10846,     -88016,     -47636,   7,  19,   1,  71,  36,   1/) !TT
   allvals2(:,   2,  36) = (/     -64291,      10634,     -87785,     -47436,   9,  20,   1,  71,  36,   1/) !TT
   allvals2(:,   2,  37) = (/     -63743,      10213,     -85549,     -47371,   9,  20,   1,  67,  35,   1/) !TT
   allvals2(:,   2,  38) = (/     -63028,       9609,     -82699,     -47203,   9,  20,   1,  68,  34,   1/) !TT
   allvals2(:,   2,  39) = (/     -62245,       8957,     -80126,     -46827,   5,  19,   1,  79,  34,   1/) !TT
   allvals2(:,   2,  40) = (/     -61227,       8116,     -77089,     -45983,   5,  19,   1,  79,  34,   1/) !TT
   allvals2(:,   2,  41) = (/     -60130,       7202,     -73423,     -45525,   5,  19,   1,  79,  34,   1/) !TT
   allvals2(:,   2,  42) = (/     -59007,       6285,     -72362,     -45211,  63,  10,   1,  79,  34,   1/) !TT
   allvals2(:,   2,  43) = (/     -57808,       5382,     -71279,     -44236,  64,  10,   1,  68,  34,   1/) !TT
   allvals2(:,   2,  44) = (/     -56355,       4675,     -69341,     -43202,  65,  11,   1,  68,  34,   1/) !TT
   allvals2(:,   2,  45) = (/     -54876,       4484,     -66934,     -42362,  76,  10,   1,  68,  34,   1/) !TT
   allvals2(:,   2,  46) = (/     -53279,       4816,     -67583,     -41212,  28,   6,   1,  34,  29,   1/) !TT
   allvals2(:,   2,  47) = (/     -51585,       5665,     -69070,     -40254,  79,   5,   1,  34,  29,   1/) !TT
   allvals2(:,   2,  48) = (/     -49870,       6848,     -69216,     -37888,  79,   5,   1,   3,  21,   1/) !TT
   allvals2(:,   2,  49) = (/     -47909,       8338,     -68313,     -33952,  26,   5,   1,  54,  21,   1/) !TT
   allvals2(:,   2,  50) = (/     -45773,       9806,     -67495,     -29902,  79,   4,   1,  53,  20,   1/) !TT
   allvals2(:,   2,  51) = (/     -43209,      11157,     -66810,     -25851,  22,   4,   1,  53,  20,   1/) !TT
   allvals2(:,   2,  52) = (/     -40467,      12189,     -66882,     -21978,  21,   3,   1,  53,  19,   1/) !TT
   allvals2(:,   2,  53) = (/     -37656,      12943,     -67204,     -17993,  21,   3,   1,  54,  19,   1/) !TT
   allvals2(:,   2,  54) = (/     -34614,      13533,     -66352,     -15346,  22,   3,   1,  63,  21,   1/) !TT
   allvals2(:,   2,  55) = (/     -31659,      13954,     -64772,     -12016,  22,   3,   1,  26,  25,   1/) !TT
   allvals2(:,   2,  56) = (/     -28603,      14274,     -63439,      -8719,  17,   2,   1,  26,  25,   1/) !TT
   allvals2(:,   2,  57) = (/     -25801,      14494,     -62006,      -5844,  17,   2,   1,  26,  26,   1/) !TT
   allvals2(:,   2,  58) = (/     -22948,      14669,     -60340,      -2195,  17,   2,   1,  23,  17,   1/) !TT
   allvals2(:,   2,  59) = (/     -20261,      14811,     -58646,        975,  16,   2,   1,  23,  16,   1/) !TT
   allvals2(:,   2,  60) = (/     -17662,      14924,     -56998,       4046,  16,   2,   1,  22,  17,   1/) !TT
   allvals2(:,   2,  61) = (/     -15244,      15001,     -55438,       5858,  16,   2,   1,  22,  17,   1/) !TT
   allvals2(:,   2,  62) = (/     -12896,      15074,     -53863,       8068,  16,   2,   1,  40,  16,   1/) !TT
   allvals2(:,   2,  63) = (/     -10610,      15160,     -52256,      10503,  16,   2,   1,  40,  16,   1/) !TT
   allvals2(:,   2,  64) = (/      -8490,      15262,     -50730,      12964,  16,   2,   1,  41,  16,   1/) !TT
   allvals2(:,   2,  65) = (/      -6522,      15388,     -49252,      14922,  16,   2,   1,   1,  17,   1/) !TT
   allvals2(:,   2,  66) = (/      -4790,      15528,     -48006,      16427,  16,   2,   1,  64,  17,   1/) !TT
   allvals2(:,   2,  67) = (/      -3223,      15682,     -46929,      19135,  17,   2,   1,  14,  26,   1/) !TT
   allvals2(:,   2,  68) = (/      -1891,      15853,     -46843,      21628,  11,   3,   1,  14,  26,   1/) !TT
   allvals2(:,   2,  69) = (/       -682,      16051,     -46728,      24130,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  70) = (/        412,      16312,     -46376,      26677,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  71) = (/       1283,      16608,     -46151,      28856,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  72) = (/       2039,      16953,     -46168,      30894,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  73) = (/       2681,      17329,     -46250,      32808,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  74) = (/       3226,      17756,     -46694,      34476,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  75) = (/       3766,      18306,     -47714,      36035,  11,   3,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  76) = (/       4339,      18981,     -53070,      37558,  20,   1,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  77) = (/       4826,      19832,     -59190,      39103,  21,   1,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  78) = (/       5109,      20768,     -64109,      39871,  21,   1,   1,  17,  27,   1/) !TT
   allvals2(:,   2,  79) = (/       5263,      21477,     -67262,      40426,  22,   1,   1,  17,  26,   1/) !TT
   allvals2(:,   2,  80) = (/       4781,      22093,     -70946,      32562,  10,   3,   1,  59,  26,   1/) !TT

   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
        ALONGX, FILL)
   istat = min(rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
        .not.ALONGX, FILL), istat)
   rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'rpn_comm_topo')
   call testutils_assert_eq(RMN_IS_OK(rpncomm_gridid), .true., 'rpncomm_gridid')

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   tkeys1(1) = -1
   tip1s(1) = 0
   keys1 => tkeys1
   ip1s => tip1s
   nkeys = fstmpio_find_3d_0(keys1, funit, 'P0', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 1, 'find P0 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find P0 associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find P0 nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find P0 nip1s')
   call testutils_assert_eq(ip1s(1), 0, 'find P0 ip1s(1)')
   call testutils_assert_eq(minval(keys1) > 0, .true., 'keys1 P0')

   nullify(data1)
   istatlist => tistats1
   hintlist_S = HINTERP4YY_CUBIC
   phintlist_S => hintlist_S(1:1)
   fileid(1) = funit
   fidlist => fileid(1:1)
   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   istat = fstmpio_rdhint(data1, istatlist, keys1, &
        phintlist_S, fidlist, rpncomm_gridid, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status P0')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .true., 'status2 P0')
   call testutils_assert_eq(associated(data1), .true., 'alloc P0')
   call testutils_assert_eq(size(istatlist), nkeys, 'nstat P0')
   call testutils_assert_eq(ubound(data1,1), maxi, 'ni P0')
   call testutils_assert_eq(ubound(data1,2), maxj, 'nj P0')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk P0')
   call testutils_assert_not_naninf(data1, 'P0 validity')

   ijkmin = lbound(data1)
   ijkmax = ubound(data1)

   ii = 1
   varname_S = 'P0'
   do ilvl=1,ubound(data1,3)
      ok_L = .true.
      call statfld_dm(data1(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
!!$         print *,'stat1',ilvl, mean, var, rmin, rmax, ijkmin, ijkmax
!!$         call flush(6)
         vals2 = (/nint(1000.*mean), nint(1000.*var), nint(1000.*rmin), nint(1000.*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
         if (all(vals2(:) == allvals2(:,ii,ilvl))) then
            ok_L = .true.
         else
            ok_L = .false.
            write(dummy_S, '(i11,",",i11,",",i11,",",i11,",",i4,",",i4,",",i4,",",i4,",",i4,",",i4)') vals2
            print '("allvals2(:,", i4, ",", i4, ") = (/", a, "/) !", a)', &
                 ii, ilvl, trim(dummy_S), varname_S(1:2)
         endif
      endif
      write(dummy_S,*) ilvl
      call testutils_assert_ok(ok_L, 'stats '//varname_S(1:2)//' '//trim(dummy_S))
   enddo


   nullify(keys1, ip1s)
   nkeys = fstmpio_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find 3d 0 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')
   call testutils_assert_eq(minval(keys1) > 0, .true., 'keys1')

   nullify(data1)
   allocate(istatlist(1:nkeys))
   hintlist_S = HINTERP4YY_CUBIC
   phintlist_S => hintlist_S(1:1)
   fileid(1) = funit
   fidlist => fileid(1:1)
   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   istat = fstmpio_rdhint(data1, istatlist, keys1, &
        phintlist_S, fidlist, rpncomm_gridid, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .true., 'status2')
   call testutils_assert_eq(associated(data1), .true., 'alloc')
   call testutils_assert_eq(size(istatlist), nkeys, 'nstat')
   call testutils_assert_eq(ubound(data1,1), maxi, 'ni')
   call testutils_assert_eq(ubound(data1,2), maxj, 'nj')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk')
   call testutils_assert_not_naninf(data1, 'TT validity')

!!$   do k=1,nkeys
!!$      print *,k,minval(data1(:,:,k)),maxval(data1(:,:,k))
!!$   enddo

   ii = 2
   varname_S = 'TT'
   do ilvl=1,ubound(data1,3)
      ok_L = .true.
      call statfld_dm(data1(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
!!$         print *,'stat1',ilvl, mean, var, rmin, rmax, ijkmin, ijkmax
!!$         call flush(6)
         vals2 = (/nint(1000.*mean), nint(1000.*var), nint(1000.*rmin), nint(1000.*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
         if (all(vals2(:) == allvals2(:,ii,ilvl))) then
            ok_L = .true.
         else
            ok_L = .false.
            write(dummy_S, '(i11,",",i11,",",i11,",",i11,",",i4,",",i4,",",i4,",",i4,",",i4,",",i4)') vals2
            print '("allvals2(:,", i4, ",", i4, ") = (/", a, "/) !", a)', &
                 ii, ilvl, trim(dummy_S), varname_S(1:2)
         endif
      endif
      write(dummy_S,*) ilvl
      call testutils_assert_ok(ok_L, 'stats '//varname_S(1:2)//' '//trim(dummy_S))
   enddo


   nullify(data2)
   keys1(nkeys/2) = -1
   istat = fstmpio_rdhint(data2, istatlist, keys1, &
        phintlist_S, fidlist, rpncomm_gridid, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'not found status')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status2')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(1:nkeys/2-1))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(nkeys/2+1:nkeys))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status4')
   do k=1,nkeys
      if (k /= nkeys/2) then
         write(dummy_S,*) k
         call testutils_assert_eq(minval(data1(:,:,k)), minval(data2(:,:,k)), 'minval '//dummy_S)
         call testutils_assert_eq(maxval(data1(:,:,k)), maxval(data2(:,:,k)), 'maxval '//dummy_S)
      endif
!!$      print *,k,minval(data1(:,:,k)),maxval(data1(:,:,k)),':',minval(data2(:,:,k)),maxval(data2(:,:,k))
   enddo
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_rdhint


!/@
subroutine test_fstmpio_rdhint_vect(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use iso_c_binding
   use ptopo_utils
   use testutils
   use fstmpio_mod
   use hinterp4yy_mod
   use statfld_dm_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"

   integer, parameter :: GNI = 80
   integer, parameter :: GNJ = 40
   integer, parameter :: HALO = 0
   logical, parameter :: ALONGX = .true.
   logical, parameter :: FILL = .false. !.true.
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat, datev, nkeys, outgridid, coregridid, &
        k, rpncomm_gridid, &
        mini,maxi,minj,maxj,lni, lnimax, li0, lnj, lnjmax, lj0
   real, pointer :: data1(:, :, :),  data2(:, :, :)
   character(len=512) :: filename_S, datev_S, dummy_S, varname_S
   integer, target :: fileid(1)
   integer, pointer :: keys1(:), keys2(:), ip1s(:), istatlist(:), fidlist(:)
   character(len=8), target :: hintlist_S(1)
   character(len=8), pointer :: phintlist_S(:)

   integer :: ijkmin(3), ijkmax(3), ii, ilvl, vals2(10), allvals2(10,2,80)
   real(REAL64) :: mean, var, rmin, rmax
   logical :: ok_L
   ! ---------------------------------------------------------------------
   allvals2 = 0
   allvals2(:,   1,   1) = (/      46811,      66306,     -76966,     235444,  18,  27,   1,  72,  15,   1/) !UU                                                
   allvals2(:,   2,   1) = (/      -9617,      29628,    -140913,      70989,   2,  10,   1,  64,  10,   1/) !VV                                                
   allvals2(:,   1,   2) = (/      50592,      64273,     -73014,     206524,  12,  27,   1,  73,  15,   1/) !UU                                                
   allvals2(:,   2,   2) = (/        153,      22142,     -73931,      73521,  76,  33,   1,  62,  11,   1/) !VV                                                
   allvals2(:,   1,   3) = (/      53146,      60273,     -77207,     198210,  66,  27,   1,  72,  16,   1/) !UU                                                
   allvals2(:,   2,   3) = (/      -1761,      21707,     -67303,      68453,  76,  30,   1,  65,  10,   1/) !VV                                                
   allvals2(:,   1,   4) = (/      55806,      56564,     -79718,     182923,  66,  27,   1,  73,  16,   1/) !UU                                                
   allvals2(:,   2,   4) = (/      -2141,      21696,     -70261,      65128,  78,  29,   1,  65,  10,   1/) !VV                                                
   allvals2(:,   1,   5) = (/      55392,      52427,     -67234,     186839,  69,  26,   1,  62,   7,   1/) !UU                                                
   allvals2(:,   2,   5) = (/       -755,      19687,     -67306,      52249,  78,  29,   1,  67,  10,   1/) !VV                                                
   allvals2(:,   1,   6) = (/      51587,      48728,     -54948,     184545,  59,  26,   1,  62,   7,   1/) !UU                                                
   allvals2(:,   2,   6) = (/        492,      17411,     -61464,      50957,  79,  29,   1,  68,  10,   1/) !VV                                                
   allvals2(:,   1,   7) = (/      45724,      44769,     -64007,     171235,  48,  25,   1,  62,   7,   1/) !UU                                                
   allvals2(:,   2,   7) = (/        893,      14909,     -56177,      52395,  80,  29,   1,  66,   9,   1/) !VV                                                
   allvals2(:,   1,   8) = (/      39034,      41229,     -61640,     153072,  49,  25,   1,  62,   7,   1/) !UU                                                
   allvals2(:,   2,   8) = (/        796,      12863,     -43675,      51841,  80,  29,   1,  66,   9,   1/) !VV                                                
   allvals2(:,   1,   9) = (/      30868,      40973,     -60252,     137167,  50,  25,   1,  60,   7,   1/) !UU                                                
   allvals2(:,   2,   9) = (/        376,      11400,     -34705,      42650,  79,  30,   1,  68,   9,   1/) !VV                                                
   allvals2(:,   1,  10) = (/      22046,      44685,     -65506,     127794,  56,  24,   1,  60,   6,   1/) !UU                                                
   allvals2(:,   2,  10) = (/       -277,      10387,     -30964,      35206,  78,  30,   1,  71,   8,   1/) !VV                                                
   allvals2(:,   1,  11) = (/      14880,      49394,     -94907,     126467,  60,  21,   1,  63,   7,   1/) !UU                                                
   allvals2(:,   2,  11) = (/       -355,       9903,     -30708,      32381,  17,   9,   1,  71,   8,   1/) !VV                                                
   allvals2(:,   1,  12) = (/      10586,      50998,    -116093,     125884,  60,  20,   1,  63,   7,   1/) !UU                                                
   allvals2(:,   2,  12) = (/       -240,       9581,     -27900,      30188,  20,   8,   1,  44,  28,   1/) !VV                                                
   allvals2(:,   1,  13) = (/       8873,      48795,    -104931,     124034,  59,  20,   1,  64,   7,   1/) !UU                                                
   allvals2(:,   2,  13) = (/       -269,       9359,     -33438,      31997,  19,   9,   1,  29,  10,   1/) !VV                                                
   allvals2(:,   1,  14) = (/       8226,      45205,     -89114,     117824,  58,  20,   1,  64,   7,   1/) !UU                                                
   allvals2(:,   2,  14) = (/       -335,       9164,     -33798,      34265,  20,   9,   1,  29,  10,   1/) !VV                                                
   allvals2(:,   1,  15) = (/       7896,      41686,     -77869,     110488,  59,  19,   1,  64,   7,   1/) !UU                                                
   allvals2(:,   2,  15) = (/       -379,       9016,     -32717,      30478,  20,   9,   1,   4,  10,   1/) !VV                                                
   allvals2(:,   1,  16) = (/       8012,      38181,     -69713,     106002,  23,  21,   1,  65,   7,   1/) !UU                                                
   allvals2(:,   2,  16) = (/       -336,       8955,     -33193,      30845,  19,  10,   1,  52,   4,   1/) !VV                                                
   allvals2(:,   1,  17) = (/       8848,      34007,     -59876,     102193,  24,  21,   1,  65,   7,   1/) !UU                                                
   allvals2(:,   2,  17) = (/       -254,       8866,     -31031,      33565,  20,  10,   1,  50,   5,   1/) !VV                                                
   allvals2(:,   1,  18) = (/      10398,      29596,     -40371,      98545,  26,  22,   1,  65,   7,   1/) !UU                                                
   allvals2(:,   2,  18) = (/       -155,       8778,     -32738,      33662,  27,  33,   1,  50,   5,   1/) !VV                                                
   allvals2(:,   1,  19) = (/      12192,      26150,     -34403,      94992,   6,  17,   1,  63,   7,   1/) !UU                                                
   allvals2(:,   2,  19) = (/        -59,       8690,     -32721,      32737,  27,  33,   1,  50,   5,   1/) !VV                                                
   allvals2(:,   1,  20) = (/      13817,      24136,     -28851,      91257,  55,  25,   1,  63,   7,   1/) !UU                                                
   allvals2(:,   2,  20) = (/         73,       8671,     -33336,      31656,  28,  33,   1,  50,   5,   1/) !VV                                                
   allvals2(:,   1,  21) = (/      15097,      22884,     -26968,      86972,  59,  24,   1,  63,   7,   1/) !UU                                                
   allvals2(:,   2,  21) = (/        248,       8720,     -32376,      32079,  28,  33,   1,  50,   6,   1/) !VV                                                
   allvals2(:,   1,  22) = (/      15962,      21719,     -27086,      84631,  40,  24,   1,  61,   7,   1/) !UU                                                
   allvals2(:,   2,  22) = (/        367,       8882,     -32166,      31837,  24,   8,   1,  41,  36,   1/) !VV                                                
   allvals2(:,   1,  23) = (/      16464,      20442,     -22698,      83398,  42,  25,   1,  61,   7,   1/) !UU                                                
   allvals2(:,   2,  23) = (/        313,       9155,     -32900,      34237,  29,  33,   1,  40,  36,   1/) !VV                                                
   allvals2(:,   1,  24) = (/      16792,      19339,     -19602,      81908,  61,  26,   1,  61,   7,   1/) !UU                                                
   allvals2(:,   2,  24) = (/        123,       9410,     -33394,      35741,  25,   8,   1,  40,  36,   1/) !VV                                                
   allvals2(:,   1,  25) = (/      17076,      18442,     -17026,      80688,  63,  26,   1,  63,   7,   1/) !UU                                                
   allvals2(:,   2,  25) = (/        -41,       9659,     -33324,      36604,  25,   8,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  26) = (/      17311,      17858,     -17971,      80234,  40,  17,   1,  60,   8,   1/) !UU
   allvals2(:,   2,  26) = (/       -169,       9924,     -32934,      38468,  25,   8,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  27) = (/      17631,      17531,     -21469,      79993,  19,  23,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  27) = (/       -256,      10219,     -32713,      40358,  25,   8,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  28) = (/      18089,      17371,     -21542,      80828,  12,  23,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  28) = (/       -205,      10522,     -32417,      41956,  25,   8,   1,  41,  35,   1/) !VV
   allvals2(:,   1,  29) = (/      18526,      17460,     -25530,      81946,  11,  23,   1,  62,   8,   1/) !UU
   allvals2(:,   2,  29) = (/       -136,      10855,     -33895,      43149,  55,  11,   1,  41,  35,   1/) !VV
   allvals2(:,   1,  30) = (/      18989,      17836,     -29187,      83096,  20,  23,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  30) = (/       -133,      11209,     -33952,      44125,  53,  35,   1,  41,  35,   1/) !VV
   allvals2(:,   1,  31) = (/      19595,      18317,     -28982,      84679,  10,  22,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  31) = (/       -168,      11572,     -35819,      44988,  54,  11,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  32) = (/      20356,      18646,     -34799,      86111,   3,  23,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  32) = (/       -194,      12033,     -37337,      46235,  43,  26,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  33) = (/      21184,      18892,     -36708,      87307,   3,  22,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  33) = (/       -168,      12501,     -38417,      47079,  52,  34,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  34) = (/      22320,      19209,     -40761,      88877,   3,  22,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  34) = (/        -70,      12963,     -40699,      47568,  44,  26,   1,  41,  36,   1/) !VV
   allvals2(:,   1,  35) = (/      23503,      19519,     -37349,      90811,   3,  22,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  35) = (/         34,      13388,     -45429,      47755,  43,  25,   1,  44,  35,   1/) !VV
   allvals2(:,   1,  36) = (/      24885,      19924,     -32486,      92976,   4,  22,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  36) = (/         94,      13841,     -48100,      47530,  42,  29,   1,  44,  36,   1/) !VV
   allvals2(:,   1,  37) = (/      26335,      20579,     -29942,      95084,   4,  22,   1,  62,   8,   1/) !UU
   allvals2(:,   2,  37) = (/         66,      14290,     -54334,      50848,  43,  28,   1,  44,  36,   1/) !VV
   allvals2(:,   1,  38) = (/      27757,      21607,     -36087,      97285,  80,  21,   1,  62,   8,   1/) !UU
   allvals2(:,   2,  38) = (/        -56,      14856,     -63018,      52431,  43,  28,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  39) = (/      29011,      22726,     -41288,     101714,  80,  21,   1,  64,   9,   1/) !UU
   allvals2(:,   2,  39) = (/       -173,      15501,     -70027,      55299,  44,  27,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  40) = (/      30286,      24022,     -38935,     113990,  80,  21,   1,  38,  30,   1/) !UU
   allvals2(:,   2,  40) = (/       -279,      16253,     -69051,      57354,  44,  27,   1,  49,   8,   1/) !VV
   allvals2(:,   1,  41) = (/      31248,      25369,     -36542,     124357,  80,  21,   1,  38,  30,   1/) !UU
   allvals2(:,   2,  41) = (/       -282,      17065,     -61831,      59005,  43,  26,   1,  49,   9,   1/) !VV
   allvals2(:,   1,  42) = (/      31893,      26816,     -37040,     129091,  27,  22,   1,  38,  30,   1/) !UU
   allvals2(:,   2,  42) = (/       -264,      18151,     -62145,      60659,  41,  29,   1,  49,   9,   1/) !VV
   allvals2(:,   1,  43) = (/      32240,      28405,     -50277,     137245,  27,  21,   1,  39,  30,   1/) !UU
   allvals2(:,   2,  43) = (/       -219,      19601,     -69304,      68056,  41,  29,   1,  45,  33,   1/) !VV
   allvals2(:,   1,  44) = (/      32292,      29934,     -47938,     151494,  26,  21,   1,  39,  30,   1/) !UU
   allvals2(:,   2,  44) = (/       -222,      21478,     -78492,      84050,  55,  11,   1,  59,  30,   1/) !VV
   allvals2(:,   1,  45) = (/      32105,      31062,     -52278,     155941,  26,  21,   1,  39,  30,   1/) !UU
   allvals2(:,   2,  45) = (/       -276,      23280,     -88929,     106795,  55,  11,   1,  59,  30,   1/) !VV
   allvals2(:,   1,  46) = (/      31710,      31825,     -48779,     150138,  60,  12,   1,  59,   9,   1/) !UU
   allvals2(:,   2,  46) = (/       -343,      24985,     -99398,     111613,  42,  29,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  47) = (/      30952,      32143,     -49238,     152176,   8,  31,   1,  58,   9,   1/) !UU
   allvals2(:,   2,  47) = (/       -315,      26337,    -102878,     124960,  55,  11,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  48) = (/      29974,      32124,     -55268,     151108,   8,  31,   1,  58,   9,   1/) !UU
   allvals2(:,   2,  48) = (/       -245,      27161,    -104945,     127635,  55,  11,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  49) = (/      28749,      31813,     -56622,     149175,   8,  31,   1,  58,   9,   1/) !UU
   allvals2(:,   2,  49) = (/       -131,      27685,    -103763,     127046,  55,  11,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  50) = (/      27416,      31252,     -55297,     148767,   7,  31,   1,  61,   8,   1/) !UU
   allvals2(:,   2,  50) = (/         13,      27734,     -99300,     130150,  55,  11,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  51) = (/      25915,      30458,     -55207,     146592,   7,  31,   1,  32,  28,   1/) !UU
   allvals2(:,   2,  51) = (/        118,      27291,     -94716,     132152,  77,  33,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  52) = (/      24420,      29497,     -53094,     143726,   7,  31,   1,  32,  28,   1/) !UU
   allvals2(:,   2,  52) = (/        199,      26582,     -96389,     131709,  77,  33,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  53) = (/      22939,      28436,     -50545,     139686,   7,  31,   1,  65,   9,   1/) !UU
   allvals2(:,   2,  53) = (/        278,      25692,     -93352,     125702,  35,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  54) = (/      21358,      27225,     -47715,     137072,   7,  31,   1,  65,   9,   1/) !UU
   allvals2(:,   2,  54) = (/        311,      24547,    -100141,     118434,  35,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  55) = (/      19863,      26001,     -45027,     133439,   7,  31,   1,  65,   9,   1/) !UU
   allvals2(:,   2,  55) = (/        278,      23314,    -100926,     113373,  35,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  56) = (/      18383,      24741,     -41042,     126769,   7,  31,   1,  65,   9,   1/) !UU
   allvals2(:,   2,  56) = (/        207,      22018,     -98359,     109146,  35,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  57) = (/      17096,      23621,     -35728,     122885,  10,  22,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  57) = (/        156,      20866,     -93899,     104916,  35,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  58) = (/      15822,      22520,     -36220,     118599,  12,  23,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  58) = (/        158,      19723,     -86768,      99068,  35,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  59) = (/      14631,      21453,     -35498,     111851,  12,  23,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  59) = (/        153,      18678,     -83067,      90875,  34,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  60) = (/      13454,      20432,     -32403,     104627,  10,  22,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  60) = (/        110,      17700,     -79868,      81171,  34,  35,   1,  45,  35,   1/) !VV
   allvals2(:,   1,  61) = (/      12308,      19539,     -31001,      98720,   8,  21,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  61) = (/         53,      16839,     -75279,      75616,  34,  35,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  62) = (/      11176,      18717,     -29930,      92427,  20,   5,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  62) = (/         -4,      16040,     -69559,      70358,  34,  35,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  63) = (/      10050,      17962,     -30743,      84692,   2,  22,   1,  63,   9,   1/) !UU
   allvals2(:,   2,  63) = (/        -61,      15300,     -63077,      63196,  34,  35,   1,  45,  34,   1/) !VV
   allvals2(:,   1,  64) = (/       9008,      17245,     -31769,      84362,   2,  22,   1,  65,   9,   1/) !UU
   allvals2(:,   2,  64) = (/        -53,      14674,     -56347,      56844,  34,  35,   1,  45,  36,   1/) !VV
   allvals2(:,   1,  65) = (/       7996,      16586,     -35072,      83914,  75,  22,   1,  65,   9,   1/) !UU
   allvals2(:,   2,  65) = (/        -55,      14126,     -52841,      55763,  76,   5,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  66) = (/       7036,      16024,     -35697,      82811,  22,   5,   1,  65,  10,   1/) !UU
   allvals2(:,   2,  66) = (/        -65,      13677,     -54235,      57449,  76,   5,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  67) = (/       6080,      15600,     -37502,      80661,  22,   5,   1,  65,  10,   1/) !UU
   allvals2(:,   2,  67) = (/        -73,      13386,     -55628,      60062,  76,   5,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  68) = (/       5239,      15300,     -38150,      74803,  22,   5,   1,  65,  10,   1/) !UU
   allvals2(:,   2,  68) = (/        -47,      13224,     -55102,      61376,  76,   5,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  69) = (/       4460,      15101,     -38145,      69426,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  69) = (/         -9,      13149,     -53092,      60663,  77,   6,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  70) = (/       3702,      15037,     -38431,      69587,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  70) = (/         19,      13180,     -52356,      58669,  77,   5,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  71) = (/       2982,      15132,     -41873,      68506,  45,  15,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  71) = (/         22,      13296,     -55017,      54256,  27,   8,   1,  46,  34,   1/) !VV
   allvals2(:,   1,  72) = (/       2263,      15276,     -46470,      67252,  45,  15,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  72) = (/          2,      13494,     -60181,      50038,  27,   8,   1,  46,  33,   1/) !VV
   allvals2(:,   1,  73) = (/       1632,      15321,     -52677,      65783,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  73) = (/         -9,      13762,     -60685,      51256,  77,   6,   1,  46,  31,   1/) !VV
   allvals2(:,   1,  74) = (/       1199,      15262,     -56335,      64376,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  74) = (/         12,      14007,     -61021,      53819,  77,   6,   1,  46,  31,   1/) !VV
   allvals2(:,   1,  75) = (/        839,      15137,     -58823,      63029,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  75) = (/         28,      14201,     -60135,      53377,  77,   6,   1,  46,  31,   1/) !VV
   allvals2(:,   1,  76) = (/        473,      15026,     -60030,      61758,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  76) = (/         61,      14328,     -58565,      54911,  77,   6,   1,  33,   5,   1/) !VV
   allvals2(:,   1,  77) = (/        157,      14819,     -57462,      60165,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  77) = (/        164,      14315,     -56658,      62802,  77,   6,   1,  33,   5,   1/) !VV
   allvals2(:,   1,  78) = (/          9,      14164,     -52113,      57665,  22,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  78) = (/        274,      13785,     -53844,      56901,  77,   6,   1,  33,   5,   1/) !VV
   allvals2(:,   1,  79) = (/        -80,      12579,     -44406,      52783,   8,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  79) = (/        276,      12165,     -47879,      48935,  77,   6,   1,  33,   5,   1/) !VV
   allvals2(:,   1,  80) = (/       -180,      10777,     -39180,      44132,  21,   5,   1,  17,  11,   1/) !UU
   allvals2(:,   2,  80) = (/        257,      10266,     -41853,      42878,  77,   6,   1,  33,   5,   1/) !VV

   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fstmpio_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   nullify(keys1, keys2, ip1s)
   nkeys = fstmpio_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(keys2), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find nip1s')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')
   call testutils_assert_eq(minval(keys1) > 0, .true., 'keys1')
   call testutils_assert_eq(minval(keys2) > 0, .true., 'keys2')
   do k=1,nkeys
      call testutils_assert_neq(keys1(k), keys2(k), 'diff keys')
   enddo

   istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
        ALONGX, FILL)
   istat = min(rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
        .not.ALONGX, FILL), istat)
   rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'rpn_comm_topo')
   call testutils_assert_eq(RMN_IS_OK(rpncomm_gridid), .true., 'rpncomm_gridid')

   nullify(data1, data2)
   allocate(istatlist(1:nkeys))
   hintlist_S = HINTERP4YY_CUBIC
   phintlist_S => hintlist_S(1:1)
   fileid(1) = funit
   fidlist => fileid(1:1)
   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   istat = fstmpio_rdhint(data1, data2, istatlist, keys1, keys2, &
        phintlist_S, fidlist, rpncomm_gridid, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .true., 'status2')
   call testutils_assert_eq(associated(data1), .true., 'alloc1')
   call testutils_assert_eq(associated(data2), .true., 'alloc2')
   call testutils_assert_eq(size(istatlist), nkeys, 'nstat')
   call testutils_assert_eq(ubound(data1,1), maxi, 'ni')
   call testutils_assert_eq(ubound(data1,2), maxj, 'nj')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk')
   call testutils_assert_eq(ubound(data2,1), maxi, 'ni')
   call testutils_assert_eq(ubound(data2,2), maxj, 'nj')
   call testutils_assert_eq(ubound(data2,3), nkeys, 'nk')
   call testutils_assert_not_naninf(data1, 'UU validity')
   call testutils_assert_not_naninf(data2, 'VV validity')
   do k=1,nkeys
      call testutils_assert_neq(minval(data1(:,:,k)), minval(data2(:,:,k)), 'minval')
      call testutils_assert_neq(maxval(data1(:,:,k)), maxval(data2(:,:,k)), 'maxval')
      !! print *,k,keys1(k), keys2(k),':',minval(data1(:,:,k)),maxval(data1(:,:,k)),':',minval(data2(:,:,k)),maxval(data2(:,:,k))
   enddo

   do ilvl=1,ubound(data1,3)
      ii = 1
      varname_S = 'UU'
      call statfld_dm(data1(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
      ok_L = .true.
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         vals2 = (/nint(1000.*mean), nint(1000.*var), nint(1000.*rmin), nint(1000.*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
         if (all(vals2(:) == allvals2(:,ii,ilvl))) then
            ok_L = .true.
         else
            ok_L = .false.
            write(dummy_S, '(i11,",",i11,",",i11,",",i11,",",i4,",",i4,",",i4,",",i4,",",i4,",",i4)') vals2
            print '("allvals2(:,", i4, ",", i4, ") = (/", a, "/) !", a)', &
                 ii, ilvl, trim(dummy_S), varname_S(1:2)
         endif
      endif
      write(dummy_S,*) ilvl
      call testutils_assert_ok(ok_L, 'stats '//varname_S(1:2)//' '//trim(dummy_S))
      ii = 2
      varname_S = 'VV'
      call statfld_dm(data2(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
      ok_L = .true.
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         vals2 = (/nint(1000.*mean), nint(1000.*var), nint(1000.*rmin), nint(1000.*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
         if (all(vals2(:) == allvals2(:,ii,ilvl))) then
            ok_L = .true.
         else
            ok_L = .false.
            write(dummy_S, '(i11,",",i11,",",i11,",",i11,",",i4,",",i4,",",i4,",",i4,",",i4,",",i4)') vals2
            print '("allvals2(:,", i4, ",", i4, ") = (/", a, "/) !", a)', &
                 ii, ilvl, trim(dummy_S), varname_S(1:2)
         endif
      endif
      write(dummy_S,*) ilvl
      call testutils_assert_ok(ok_L, 'stats '//varname_S(1:2)//' '//trim(dummy_S))
   enddo

   nullify(data1, data2)
   keys1(nkeys/2) = -1
   keys2(nkeys/2) = -1
   istat = fstmpio_rdhint(data1, data2, istatlist, keys1, keys2, &
        phintlist_S, fidlist, rpncomm_gridid, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'not found status')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status2')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(1:nkeys/2-1))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(nkeys/2+1:nkeys))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status4')
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpio_rdhint_vect
