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
subroutine test_fst()
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use ptopo_utils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   integer :: istat, myproc, ndomains, idomain, ngrids, igrid
   character(len=512) :: dfiles_S, bcmk_S
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   ndomains = 1
   call ptopo_init_var(ndomains, idomain, ngrids, igrid)

   call testutils_set_name('test_fst')

   istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
   istat = clib_getenv('ATM_MODEL_DFILES',dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      print *,'ERROR: ATM_MODEL_DFILES not defined'
      return
   endif
   bcmk_S = trim(dfiles_S)//'/bcmk/'

   call testutils_set_name('test_fst_open_notfound')
   call test_fst_open_notfound(bcmk_S)

   call testutils_set_name('test_fst_find_read_0')
   call test_fst_find_read(bcmk_S)
!!$   call test_fst_find_fuzz(bcmk_S)

   call testutils_set_name('test_fst_get_vgrid')
   call test_fst_get_vgrid(bcmk_S)

   call testutils_set_name('test_fst_find_3d_0')
   call test_fst_find_3d_0(bcmk_S)

   call testutils_set_name('test_fst_find_3d_0b')
   call test_fst_find_3d_0b(bcmk_S)

   call testutils_set_name('test_fst_find_3d_vect')
   call test_fst_find_3d_vect(bcmk_S)

   call testutils_set_name('test_fst_find_3d_vectb')
   call test_fst_find_3d_vectb(bcmk_S)

   call testutils_set_name('test_fst_rdhint')
   call test_fst_rdhint(bcmk_S)

   call testutils_set_name('test_fst_rdhint_vect')
   call test_fst_rdhint_vect(bcmk_S)

! !!$   call test_fst_write()

   call testutils_stats_print()
   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fst


!/@
subroutine test_fst_open_notfound(F_bcmk_S)
   use testutils
   use fst_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, istat
   ! ---------------------------------------------------------------------
   funit = fst_open(trim(F_bcmk_S)//'__does_not_exists__',FST_READONLY)
   call testutils_assert_ok(.not.RMN_IS_OK(funit),'test_fst_open_notfound','')
   if (funit > 0) istat = fst_close(funit)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_open_notfound


!/@
subroutine test_fst_find_read(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   logical,parameter :: DIR_IS_OK_L = .true.
   real(REAL64),parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat,datev,datev2,datev3,key,gridid,datevfuzz, &
        ni, nj, ig1, ig2, ig3, ig4,ip1,kind
   real, pointer :: data(:,:,:)
   character(len=512) :: filename_S,datev_S,dummy_S
   character(len=2) :: grtyp_S
   real :: zp1
   real(REAL64) :: nhours_8
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev,datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S,FST_READONLY)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fst_find_read','fst_open')

   key = fst_find(funit,'NONE',datev,RMN_ANY_I,0,0)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fst_find_read','fst_find not found')

   key = fst_find(funit,'TT',datev,RMN_ANY_I,0,0)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find')

   istat = fst_read(key,data,funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_read')
   istat = RMN_ERR
   if (associated(data) .and. &
        size(data,1) == 200 .and. &
        size(data,2) == 100 .and. &
        size(data,3) == 1  .and. &
        nint(minval(data)*100.) == -5089 .and. & 
        nint(maxval(data)*100.) == -2618) istat = RMN_OK
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_read valuse')
   if (associated(data)) deallocate(data,stat=istat)

   gridid = fst_get_gridid(funit,key)
   call testutils_assert_ok(RMN_IS_OK(gridid),'test_fst_find_read','fst_get_gridid')
   istat = ezgprm(gridid, grtyp_S(1:1), ni, nj, ig1, ig2, ig3, ig4)
   istat = RMN_ERR
   if (grtyp_S(1:1) == 'G' .and. ni == 200 .and. nj == 100 .and. &
        ig1+ig2+ig3+ig4 == 0) istat = RMN_OK
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_get_gridid values')
!!$   ier = ezgxprm(gdid, ni, nj, grtyp, ig1, ig2, ig3, ig4, grref, ig1ref, ig2ref, ig3ref, ig4ref)
!!$   ier = ezgfstp(gdid, nomvarx, typvarx, etikx, nomvary, typvary, etiky, ip1, ip2, ip3, dateo, deet, npas, nbits)
!!$   ier = gdgaxes(gdid, ax, ay)
   istat = gdrls(gridid)

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   datevfuzz = 3600
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_near')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_near value')

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LE)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_le')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_le value')

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GE)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_ge not found')

   datev_S = '20090426.220000'
   call datp2f(datev2,datev_S)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GE)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_ge')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_ge value')


   datev_S = '20090427.000000'
   call datp2f(datev2,datev_S)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GT)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_gt not found')

   datev_S = '20090427.000000'
   call datp2f(datev2,datev_S)
   datevfuzz = 3600*6
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LT)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_lt not found')

   datev_S = '20090427.000000'
   call datp2f(datev3,datev_S)
   nhours_8 = -40.D0/SEC_PER_HR
   call incdatr(datev2,datev3,nhours_8)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_gt as eq')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_gt as eq value')
   call datf2p(datev_S,datev2)


   datev_S = '20090427.000000'
   call datp2f(datev3,datev_S)
   nhours_8 = -1.D0
   call incdatr(datev2,datev3,nhours_8)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_gt')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_gt value')


   datev_S = '20090427.000000' 
   call datp2f(datev3,datev_S)
   nhours_8 = 40.D0/SEC_PER_HR
   call incdatr(datev2,datev3,nhours_8)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_lt as eq')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_lt as eq value')


   datev_S = '20090427.000000' 
   call datp2f(datev3,datev_S)
   nhours_8 = 1.D0
   call incdatr(datev2,datev3,nhours_8)
   key = fst_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find_fuzz_lt')
   call testutils_assert_ok(datev2==datev,'test_fst_find_read','fst_find_fuzz_lt value')
   call datf2p(datev_S,datev2)

   istat = fst_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_close')

   !----

   filename_S = trim(F_bcmk_S)//'geophy/Gem_geophy.fst'
   funit = fst_open(filename_S,FST_READONLY)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fst_find_read2','fst_open')

   datev = RMN_ANY_DATE
   zp1 = 1.
   kind = RMN_CONV_ARBITRARY
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPOLD, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   key = fst_find(funit,'J1',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find ip1>0 old')

   datev = RMN_ANY_DATE
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPNEW, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   key = fst_find(funit,'J1',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find ip1>0 new')

   datev = RMN_ANY_DATE
   ip1 = 1200
   key = fst_find(funit,'ME',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find ip1=1200')

   datev = RMN_ANY_DATE
   ip1 = 0
   key = fst_find(funit,'ME',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find ip1=0 for 1200')

   datev = RMN_ANY_DATE
   ip1 = 0
   key = fst_find(funit,'MG',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find ip1=0')

   datev = RMN_ANY_DATE
   ip1 = 1200
   key = fst_find(funit,'MG',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find ip1=1200 for 0')

   istat = fst_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read2','fst_close')

   !----

   filename_S = trim(F_bcmk_S)
   funit = fst_open(filename_S,FST_READONLY,DIR_IS_OK_L)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fst_find_read','fst_open dir')
   key = fst_find(funit,'TT',datev,RMN_ANY_I,0,0)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find dir')

   istat = fst_read(key,data,funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_read dir')
   istat = RMN_ERR
   if (associated(data) .and. &
        size(data,1) == 200 .and. &
        size(data,2) == 100 .and. &
        size(data,3) == 1  .and. &
        nint(minval(data)*100.) == -5089 .and. & 
        nint(maxval(data)*100.) == -2618) istat = RMN_OK
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_read dir valuse')
   if (associated(data)) deallocate(data,stat=istat)

   gridid = fst_get_gridid(funit,key)
   call testutils_assert_ok(RMN_IS_OK(gridid),'test_fst_find_read','fst_get_gridid dir')
   istat = ezgprm(gridid, grtyp_S(1:1), ni, nj, ig1, ig2, ig3, ig4)
   istat = RMN_ERR
   if (grtyp_S(1:1) == 'G' .and. ni == 200 .and. nj == 100 .and. &
        ig1+ig2+ig3+ig4 == 0) istat = RMN_OK
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read','fst_get_gridid dir values')
!!$   ier = ezgxprm(gdid, ni, nj, grtyp, ig1, ig2, ig3, ig4, grref, ig1ref, ig2ref, ig3ref, ig4ref)
!!$   ier = ezgfstp(gdid, nomvarx, typvarx, etikx, nomvary, typvary, etiky, ip1, ip2, ip3, dateo, deet, npas, nbits)
!!$   ier = gdgaxes(gdid, ax, ay)
   istat = gdrls(gridid)

   datev = RMN_ANY_DATE
   zp1 = 1.
   kind = RMN_CONV_ARBITRARY
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPOLD, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   key = fst_find(funit,'J1',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_find_read','fst_find dir J1 ip1>0 old')

   istat = fst_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_find_read2','fst_close dir')
   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_find_read

!/@
subroutine test_fst_get_vgrid(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use vGrid_Descriptors
   use fst_mod
   use vgrid_from_file_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat, datev
   character(len=512) :: filename_S, datev_S, lvltyp_S
   integer, pointer :: keys1(:), ip1s(:)
   type(vgrid_descriptor) :: vgrid1
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

!!$   key = fst_find(funit,'TT',datev,RMN_ANY_I,0,0)
!!$   call testutils_assert_ok(RMN_IS_OK(key),'fst_find')
!!$
!!$   istat = fst_get_vgrid(funit, key, vgrid1, ip1s, lvltyp_S)
!!$   call testutils_assert_ok(RMN_IS_OK(istat),'fst_get_vgd status')
!!$   call testutils_assert_eq(size(ip1s), 80, 'fst_get_vgd nip1s')
!!$   call testutils_assert_eq(lvltyp_S, 'T', 'fst_get_vgd F_lvltyp_S')
!!$   call testutils_assert_eq(size(ip1s), 80, 'fst_get_vgd nip1s')

   istat = vgrid_from_file(funit, 'TT', datev, vgrid1, &
        ip1s, lvltyp_S, keys1)
   call testutils_assert_ok(RMN_IS_OK(istat),'fst_get_vgd status')
   call testutils_assert_eq(size(ip1s), 80, 'fst_get_vgd nip1s')
!!$   call testutils_assert_eq(lvltyp_S, 'T', 'fst_get_vgd F_lvltyp_S')
   call testutils_assert_eq(size(ip1s), 80, 'fst_get_vgd nip1s')
!!$   istat = fst_get_vgrid(F_fileid, key, F_vgrid, F_ip1s, F_lvltyp_S)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_get_vgrid


!/@
subroutine test_fst_rdhint(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   use hinterp4yy_mod
   use statfld_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer, parameter :: GNI = 80
   integer, parameter :: GNJ = 40
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat, datev, nkeys, outgridid, coregridid, &
        k, mini,maxi,minj,maxj
   real, pointer :: data1(:, :, :), data2(:, :, :)
   character(len=512) :: filename_S, datev_S, dummy_S, varname_S
   integer, target :: tkeys1(1), tip1s(1), fileid(1), tistats1(1)
   integer, pointer :: keys1(:), ip1s(:), istatlist(:), fidlist(:)
   character(len=8), target :: hintlist_S(1)
   character(len=8), pointer :: phintlist_S(:)

   integer :: ijkmin(3), ijkmax(3), ii, ilvl, vals2(10), allvals2(10,2,80)
   real(REAL64) :: mean, var, rmin, rmax
   ! ---------------------------------------------------------------------

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

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')


   mini = 1
   minj = 1
   maxi = GNI
   maxj = GNJ

   tkeys1(1) = -1
   tip1s(1) = 0
   keys1 => tkeys1
   ip1s => tip1s
   nkeys = fst_find_3d_0(keys1, funit, 'P0', datev, ip1s, 0, 0)
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
   istat = fst_rdhint(data1, istatlist, keys1, phintlist_S, fidlist, &
        outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status P0')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .true., 'status2 P0')
   call testutils_assert_eq(associated(data1), .true., 'alloc P0')
   call testutils_assert_eq(size(istatlist), nkeys, 'nstat P0')
   call testutils_assert_eq(ubound(data1,1), maxi, 'ni P0')
   call testutils_assert_eq(ubound(data1,2), maxj, 'nj P0')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk P0')

   ii = 1
   varname_S = 'P0'
   do ilvl=1,ubound(data1,3)
      call statfld(data1(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
      print *,ilvl, mean, var, rmin, rmax, ijkmin, ijkmax
      call flush(6)
      vals2 = (/nint(1000.*mean), nint(1000.*var), nint(1000.*rmin), nint(1000.*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
      if (all(vals2(:) == allvals2(:,ii,ilvl))) then
         write(dummy_S,*) ilvl
         call testutils_assert_ok(.true., 'stats '//varname_S(1:2)//' '//trim(dummy_S))
      else
         call testutils_assert_ok(.false., 'stats '//varname_S(1:2)//' '//trim(dummy_S))
         write(dummy_S, '(i11,",",i11,",",i11,",",i11,",",i4,",",i4,",",i4,",",i4,",",i4,",",i4)') vals2
         print '("allvals2(:,", i4, ",", i4, ") = (/", a, "/) !", a)', &
              ii, ilvl, trim(dummy_S), varname_S(1:2)
      endif
   enddo


   nullify(keys1, ip1s)
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find 3d 0 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find associated')

   nullify(data1)
   allocate(istatlist(1:nkeys))
   hintlist_S = HINTERP4YY_CUBIC
   phintlist_S => hintlist_S(1:1)
   fileid(1) = funit
   fidlist => fileid(1:1)
   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   istat = fst_rdhint(data1, istatlist, keys1, phintlist_S, fidlist, &
        outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .true., 'status2')
   call testutils_assert_eq(associated(data1), .true., 'alloc')
   call testutils_assert_eq(ubound(data1,1), GNI, 'ni')
   call testutils_assert_eq(ubound(data1,2), GNJ, 'nj')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk')
!!$   do k=1,nkeys
!!$      print *,k,minval(data1(:,:,k)),maxval(data1(:,:,k))
!!$   enddo

   nullify(data2)
   keys1(nkeys/2) = -1
   istat = fst_rdhint(data2, istatlist, keys1, phintlist_S, fidlist, &
        outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'not found status')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status2')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(1:nkeys/2-1))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(nkeys/2+1:nkeys))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status4')
   do k=1,nkeys
      if (k /= nkeys/2) then
         call testutils_assert_eq(minval(data1(:,:,k)), minval(data2(:,:,k)), 'minval')
         call testutils_assert_eq(maxval(data1(:,:,k)), maxval(data2(:,:,k)), 'maxval')
      endif
!!$      print *,k,minval(data1(:,:,k)),maxval(data1(:,:,k)),':',minval(data2(:,:,k)),maxval(data2(:,:,k))
   enddo

   nullify(data1)
   allocate(istatlist(1:nkeys))
   hintlist_S = HINTERP4YY_CUBIC
   phintlist_S => hintlist_S(1:1)
   fileid(1) = funit
   fidlist => fileid(1:1)
   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   keys1(nkeys/2) = -1
   istat = fst_rdhint(data1, istatlist, keys1, phintlist_S, fidlist, &
        outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status 1miss')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .false., 'status2 1miss')
   call testutils_assert_eq(associated(data1), .true., 'alloc 1miss')
   call testutils_assert_eq(ubound(data1,1), GNI, 'ni 1miss')
   call testutils_assert_eq(ubound(data1,2), GNJ, 'nj 1miss')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk 1miss')
   do k=1,nkeys
      if (k == nkeys/2) then
         call testutils_assert_eq(RMN_IS_OK(istatlist(k)), .false., 'hstatus nk/2 1miss')
      else
         call testutils_assert_eq(RMN_IS_OK(istatlist(k)), .true., 'hstatus 1miss')
      endif
   enddo

   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_rdhint

!/@
subroutine test_fst_rdhint_vect(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   use hinterp4yy_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer, parameter :: GNI = 80
   integer, parameter :: GNJ = 40
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, istat, datev, nkeys, outgridid, coregridid, k
   real, pointer :: data1(:, :, :),  data2(:, :, :)
   character(len=512) :: filename_S, datev_S
   integer, target :: fileid(1)
   integer, pointer :: keys1(:), keys2(:), ip1s(:), istatlist(:), fidlist(:)
   character(len=8), target :: hintlist_S(1)
   character(len=8), pointer :: phintlist_S(:)
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   nullify(keys1, keys2, ip1s)
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find associated')

   nullify(data1, data2)
   allocate(istatlist(1:nkeys))
   hintlist_S = HINTERP4YY_CUBIC
   phintlist_S => hintlist_S(1:1)
   fileid(1) = funit
   fidlist => fileid(1:1)
   outgridid = ezqkdef(GNI, GNJ, 'G', 0,0,0,0,0)
   coregridid = outgridid
   istat = fst_rdhint(data1, data2, istatlist, keys1, keys2, &
        phintlist_S, fidlist, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'status')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist)), .true., 'status2')
   call testutils_assert_eq(associated(data1), .true., 'alloc')
   call testutils_assert_eq(ubound(data1,1), GNI, 'ni1')
   call testutils_assert_eq(ubound(data1,2), GNJ, 'nj1')
   call testutils_assert_eq(ubound(data1,3), nkeys, 'nk1')
   call testutils_assert_eq(ubound(data2,1), GNI, 'ni2')
   call testutils_assert_eq(ubound(data2,2), GNJ, 'nj2')
   call testutils_assert_eq(ubound(data2,3), nkeys, 'nk2')
   do k=1,nkeys
      call testutils_assert_neq(minval(data1(:,:,k)), minval(data2(:,:,k)), 'minval')
      call testutils_assert_neq(maxval(data1(:,:,k)), maxval(data2(:,:,k)), 'maxval')
!!$      print *,k,keys1(k), keys2(k),':',minval(data1(:,:,k)),maxval(data1(:,:,k)),':',minval(data2(:,:,k)),maxval(data2(:,:,k))
   enddo

   keys1(nkeys/2) = -1
   istat = fst_rdhint(data1, data2, istatlist, keys1, keys2, &
        phintlist_S, fidlist, outgridid, coregridid)
   call testutils_assert_eq(RMN_IS_OK(istat), .true., 'not found status')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status2')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(1:nkeys/2-1))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(minval(istatlist(nkeys/2+1:nkeys))), .true., 'not found status3')
   call testutils_assert_eq(RMN_IS_OK(istatlist(nkeys/2)), .false., 'not found status4')
   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_rdhint_vect


!/@
subroutine test_fst_find_3d_0(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, datev, datev2, datevfuzz, nkeys
   character(len=512) :: filename_S, datev_S
   integer, pointer :: keys1(:), ip1s(:), ip1s2(:)
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   nullify(keys1, ip1s)
   nkeys = fst_find_3d_0(keys1, funit, 'NONE', datev, ip1s, 0, 0)
   call testutils_assert_ok(nkeys <= 0, 'find not found')

   nullify(keys1, ip1s)
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), nkeys, 'find nip1s')
   call testutils_assert_eq(keys1(1), 1036, 'find keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 81932, 'find keys1(nkeys)')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')

   nullify(keys1)
   ip1s2 => ip1s(1:5)
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find 5 nkeys1')
   call testutils_assert_eq(size(ip1s2), nkeys, 'find 5 nip1s')
   call testutils_assert_eq(keys1(1), 1036, 'find 5 keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 5132, 'find 5 keys1(nkeys)')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, ip1s)
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   nullify(keys1, ip1s)
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_find_3d_0


!/@
subroutine test_fst_find_3d_0b(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
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
   ! ---------------------------------------------------------------------
   keys1 => tkeys1
   keys2 => tkeys2
   ip1s => tip1s
   ip1s2 => tip1s2

   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

!!$   nullify(keys1, ip1s)
   keys1 = -1
   ip1s = -1
   nkeys = fst_find_3d_0(keys1, funit, 'NONE', datev, ip1s, 0, 0)
   call testutils_assert_ok(nkeys <= 0, 'find not found')

!!$   nullify(keys1, ip1s)
   keys1 = -1
   ip1s = -1
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find associated keys1')
   call testutils_assert_eq(associated(ip1s, tip1s), .true., 'find associated ip1s')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), NMAXKEYS, 'find nip1s')
   call testutils_assert_eq(keys1(1), 1037, 'find keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 81933, 'find keys1(nkeys)')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')

!!$   nullify(keys1)
   keys1 = -1
   ip1s2 => ip1s(1:5)
!!$   call testutils_assert_eq(associated(ip1s2, ip1s), .true., 'find 5 associated ip1s 0')
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find 5 nkeys1')
   call testutils_assert_eq(size(ip1s2), nkeys, 'find 5 nip1s')
   call testutils_assert_eq(associated(keys1, tkeys1), .true., 'find 5 associated keys1')
!!$   call testutils_assert_eq(associated(ip1s2, ip1s), .true., 'find 5 associated ip1s')
!!$   call testutils_assert_eq(associated(ip1s2, tip1s), .true., 'find 5 associated tip1s')
   call testutils_assert_eq(keys1(1), 1037, 'find 5 keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 5133, 'find 5 keys1(nkeys)')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, ip1s)
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
!!$   nullify(keys1, ip1s)
   keys1 = -1
   ip1s = -1
   nkeys = fst_find_3d_0(keys1, funit, 'TT', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_find_3d_0b


subroutine test_fst_find_3d_vect(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer :: funit, datev, datev2, datevfuzz, nkeys, k
   character(len=512) :: filename_S, datev_S
   integer, pointer :: keys1(:), keys2(:), ip1s(:), ip1s2(:)
   ! ---------------------------------------------------------------------
   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

   nullify(keys1, keys2, ip1s)
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s, 0, 0)
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
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), nkeys, 'find 5 nkeys1')
   call testutils_assert_eq(size(keys2), nkeys, 'find 5 nkeys1')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, keys2, ip1s)
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'find fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   nullify(keys1, keys2, ip1s)
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'find fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'find fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_find_3d_vect

subroutine test_fst_find_3d_vectb(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2017-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer, parameter :: NMAXKEYS = 1000
   integer :: funit, datev, datev2, datevfuzz, nkeys, k
   character(len=512) :: filename_S, datev_S
   integer, dimension(NMAXKEYS), target :: tkeys1, tkeys2, tip1s, tip1s2
   integer, pointer :: keys1(:), keys2(:), ip1s(:), ip1s2(:)
   ! ---------------------------------------------------------------------
   keys1 => tkeys1
   keys2 => tkeys2
   ip1s => tip1s
   ip1s2 => tip1s2

   datev_S = '20090427.000000'
   call datp2f(datev, datev_S)
   filename_S = trim(F_bcmk_S)//'2009042700_000'

   funit = fst_open(filename_S, FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit), .true., 'open')

!!$   nullify(keys1, keys2, ip1s)
   keys1 = -1
   keys2 = -1
   ip1s = -1
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s, 0, 0)
   call testutils_assert_eq(nkeys, 80, 'find nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(keys2) .and. associated(ip1s), .true., 'find associated')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find nkeys1')
   call testutils_assert_eq(size(keys2), NMAXKEYS, 'find nkeys1')
   call testutils_assert_eq(size(ip1s), NMAXKEYS, 'find nip1s')
   call testutils_assert_eq(keys1(1), 1730575, 'find keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 2154511, 'find keys1(nkeys)')
   call testutils_assert_eq(keys2(1), 1731599, 'find keys2(1)')
   call testutils_assert_eq(keys2(nkeys), 2155535, 'find keys2(nkeys)')
   call testutils_assert_eq(ip1s(1), 97642568, 'find ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find ip1s(nkeys)')
   do k=1,nkeys
      call testutils_assert_neq(keys1(k), keys2(k), 'diff keys')
   enddo

!!$   nullify(keys1, keys2)
   keys1 = -1
   keys2 = -1
   ip1s2 => ip1s(1:5)
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev, ip1s2, 0, 0)
   call testutils_assert_eq(nkeys, 5, 'find 5 nkeys')
   call testutils_assert_eq(associated(keys1) .and. associated(ip1s2), .true., 'find associated')
   call testutils_assert_eq(size(keys1), NMAXKEYS, 'find 5 nkeys1')
   call testutils_assert_eq(size(keys2), NMAXKEYS, 'find 5 nkeys1')
   call testutils_assert_eq(size(ip1s2), nkeys, 'find 5 nip1s')
   call testutils_assert_eq(keys1(1), 1730575, 'find 5 keys1(1)')
   call testutils_assert_eq(keys1(nkeys), 1738767, 'find 5 keys1(nkeys)')
   call testutils_assert_eq(keys2(1), 1731599, 'find 5 keys2(1)')
   call testutils_assert_eq(keys2(nkeys), 1739791, 'find 5 keys2(nkeys)')
   call testutils_assert_eq(ip1s2(1), 97642568, 'find 5 ip1s(1)')
   call testutils_assert_eq(ip1s2(nkeys), 96569992, 'find 5 ip1s(nkeys)')

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
   datevfuzz = 3600
   nullify(keys1, keys2, ip1s)
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_ok(nkeys <= 0, 'find fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2, datev_S)
!!$   nullify(keys1, keys2, ip1s)
   keys1 = -1
   keys2 = -1
   ip1s = -1
   nkeys = fst_find_3d_vect(keys1, keys2, funit, 'UU', 'VV', datev2, ip1s, 0, 0, datevfuzz)
   call testutils_assert_eq(nkeys, 80, 'find fuzz_near nkeys')
   call testutils_assert_ok(datev2==datev, 'find fuzz_near value')
   call testutils_assert_eq(ip1s(1), 97642568, 'find fuzz_near ip1s(1)')
   call testutils_assert_eq(ip1s(nkeys), 93423264, 'find fuzz_near ip1s(nkeys)')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_find_3d_vectb



!/@
subroutine test_fst_write()
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fst_mod
   use ezgrid_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-01
   !@argument
!@/
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>
   real,parameter :: MYVALUE = 3.3
   integer,parameter :: NI0=50,NJ0=30,NK0=3
   character(len=256) :: nomvar_S,filename_S
   character(len=2) :: grtyp_S,grtyp2_S,grref_S,grref2_S
   logical :: ok_L,ok2_L
   integer :: funit,istat,gridid,gridid2,key,ig1,ig2,ig3,ig4,ip1,i,j,k,datev,ip1list(NK0),ni,nj,ni2, nj2, ig12, ig22, ig32, ig42,nij(2),nij2(2),ij0(2),ij02(2),ig14(4),ig142(4)
   real,pointer :: data2d(:,:),data3d(:,:,:)
   real :: ax(NI0,1),ay(1,NJ0)
   real,pointer :: ax1(:,:),ax2(:,:),ay1(:,:),ay2(:,:)
   ! ---------------------------------------------------------------------
   filename_S = '__test_fst_to-be-removed__.fst'
   istat = clib_unlink(trim(filename_S))
   funit = fst_open(filename_S)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fst_write:open','')

   nomvar_S = 'ZX'

   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
   do i=1,NI0
      ax(i,1) = 10.+float(i)*0.25
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)*0.25
   enddo
   gridid = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)

   allocate(data2d(NI0,NJ0),data3d(NI0,NJ0,NK0),stat=istat)
   data2d = MYVALUE
   data3d = MYVALUE
   ip1 = 0
!!$   lvlid = 
   do k=1,NK0
      ip1list(k) = (k+2)*3
   enddo

   nomvar_S = 'ZX2d'
   istat = fst_write(funit,nomvar_S,data2d,gridid,ip1,F_npak=FST_NPAK_FULL32)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_write:write_2d_r4','')

   nomvar_S = 'ZX3d'
   istat = fst_write(funit,nomvar_S,data3d,gridid,ip1list,F_npak=FST_NPAK_FULL32)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_write:write_3d_r4','')

   !TODO: test with a vgrid instead of ip1list

   if (funit > 0) istat = fst_close(funit)
   deallocate(data2d,data3d,stat=istat)

   !- Checking

   funit = fst_open(filename_S,FST_READONLY)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fst_write:open','')

   nomvar_S = 'ZX2d'
   datev = RMN_ANY_DATE
   key = fst_find(funit,nomvar_S,datev,RMN_ANY_I,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fst_write:find','2d')

   istat = fst_read(key,data3d,funit,gridid2)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_write:read','2d')
   ok_L = .false. ; ok2_L = .false.
   if (associated(data3d)) then
      ok_L = all(shape(data3d)==(/NI0,NJ0,1/))
      ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
   endif
   call testutils_assert_ok(ok_L,'test_fst_write:read','data2d shape')
   call testutils_assert_ok(ok2_L,'test_fst_write:read','data2d')
   if (.not.ok2_L) then
      print *,'data2d min,max:',minval(data3d),maxval(data3d)
   endif
   ok_L = ezgrid_samegrid(gridid,gridid2)
   call testutils_assert_ok(ok_L,'test_fst_write:read','data2d grid')
   if (.not.ok_L) then
      print *,'gridid1,2=',gridid,gridid2
      istat = ezgprm(gridid, grtyp_S(1:1), ni, nj, ig1, ig2, ig3, ig4)
      istat = ezgprm(gridid2, grtyp2_S(1:1), ni2, nj2, ig12, ig22, ig32, ig42)
      call testutils_assert_eq(grtyp2_S(1:1),grtyp_S(1:1),'grtyp_S')
      call testutils_assert_eq((/ni2, nj2/),(/ni, nj/),'nij')
      call testutils_assert_eq((/ig12, ig22, ig32, ig42/),(/ig1, ig2, ig3, ig4/),'ig1234')
      istat  = ezgrid_params(gridid,nij,grtyp_S,grref_S,ig14,ij0,ax1,ay1)
      istat  = ezgrid_params(gridid2,nij2,grtyp2_S,grref2_S,ig142,ij02,ax2,ay2)
      call testutils_assert_eq(grtyp2_S(1:1)//grref2_S(1:1),grtyp_S(1:1)//grref_S(1:1),'grtyp_S//grref_S')
      call testutils_assert_eq(grtyp2_S(1:1)//grref2_S(1:1),'ZE','grtyp_S//grref_S ZE')
      call testutils_assert_eq(nij2,nij,'nij')
      call testutils_assert_eq(nij2,(/NI0,NJ0/),'NIJ0')
      call testutils_assert_eq(ig142,ig14,'ig14')
      call testutils_assert_eq(ij02,ij0,'ij0')
      call testutils_assert_eq(ax2(:,1),ax1(:,1),'ax')
      call testutils_assert_eq(ay2(1,:),ay1(1,:),'ay')
      call testutils_assert_eq(ax2(:,1),ax(:,1),'ax0')
      call testutils_assert_eq(ay2(1,:),ay(1,:),'ay0')
   endif


   nomvar_S = 'ZX3d'

   do k=1,NK0
      if (associated(data3d)) deallocate(data3d,stat=istat)
      datev = RMN_ANY_DATE
      key = fst_find(funit,nomvar_S,datev,ip1list(k),RMN_ANY_I,RMN_ANY_I)
      call testutils_assert_ok(RMN_IS_OK(key),'test_fst_write:find','3d')

      istat = fst_read(key,data3d,funit,gridid2)
      call testutils_assert_ok(RMN_IS_OK(istat),'test_fst_write:read','3d')
      ok_L = .false. ; ok2_L = .false.
      if (associated(data3d)) then
         ok_L = all(shape(data3d)==(/NI0,NJ0,1/))
         ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
      endif
      call testutils_assert_ok(ok_L,'test_fst_write:read','data3d shape')
      call testutils_assert_ok(ok2_L,'test_fst_write:read','data3d')
      if (.not.ok2_L) then
         print *,'data3d min,max:',minval(data3d),maxval(data3d)
      endif
      ok_L = ezgrid_samegrid(gridid,gridid2)
      call testutils_assert_ok(ok_L,'test_fst_write:read','data3d grid')
   enddo

   if (funit > 0) istat = fst_close(funit)

   if (associated(data2d)) deallocate(data2d,stat=istat)
   if (associated(data3d)) deallocate(data3d,stat=istat)

   istat = clib_unlink(trim(filename_S))
   ! ---------------------------------------------------------------------
   return
end subroutine test_fst_write
