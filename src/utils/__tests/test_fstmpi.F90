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
subroutine test_fstmpi()
   use testutils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   integer :: istat,myproc
   character(len=512) :: dfiles_S,bcmk_S
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()

   istat = clib_getenv('ATM_MODEL_DFILES',dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      print *,'ERROR: ATM_MODEL_DFILES not defined'
      return
   endif
   bcmk_S = trim(dfiles_S)//'/bcmk/'
   call testutils_set_name('test_fstmpi_open_notfound')
   call test_fstmpi_open_notfound(bcmk_S)
   call testutils_set_name('test_fstmpi_find_read')
   call test_fstmpi_find_read(bcmk_S)
!!$   call testutils_set_name('test_fstmpi_find_fuzz')
!!$   call test_fstmpi_find_fuzz(bcmk_S)
   call testutils_set_name('test_fstmpi_write')
   call test_fstmpi_write()

   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpi


!/@
subroutine test_fstmpi_open_notfound(F_bcmk_S)
   use testutils
   use fstmpi_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
   integer :: funit, istat
   ! ---------------------------------------------------------------------
   funit = fstmpi_open(trim(F_bcmk_S)//'__does_not_exists__',FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit),.false.,'')
   if (funit > 0) istat = fstmpi_close(funit)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpi_open_notfound


!/@
subroutine test_fstmpi_find_read(F_bcmk_S)
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use testutils
   use fstmpi_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@argument
   character(len=*), intent(in) :: F_bcmk_S
!@/
#include <rmnlib_basics.hf>
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

   funit = fstmpi_open(filename_S,FST_READONLY)
   call testutils_assert_eq(RMN_IS_OK(funit),.true.,'open')

   key = fstmpi_find(funit,'NONE',datev,RMN_ANY_I,0,0)
   call testutils_assert_eq(RMN_IS_OK(funit),.true.,'find not found')

   key = fstmpi_find(funit,'TT',datev,RMN_ANY_I,0,0)
   call testutils_assert_eq(RMN_IS_OK(funit),.true.,'find')

   istat = fstmpi_read(key,data,funit)
   call testutils_assert_eq(RMN_IS_OK(funit),.true.,'read status')
   istat = RMN_ERR
   if (associated(data) .and. &
        size(data,1) == 200 .and. &
        size(data,2) == 100 .and. &
        size(data,3) == 1  .and. &
        nint(minval(data)*100.) == -5089 .and. & 
        nint(maxval(data)*100.) == -2618) istat = RMN_OK
   call testutils_assert_eq(istat,RMN_OK,'read values')
   if (associated(data)) deallocate(data,stat=istat)

   gridid = fstmpi_get_gridid(funit,key)
   call testutils_assert_eq(RMN_IS_OK(gridid),.true.,'get_gridid status')
   istat = ezgprm(gridid, grtyp_S, ni, nj, ig1, ig2, ig3, ig4)
   istat = RMN_ERR
   call testutils_assert_eq(grtyp_S(1:1),'G','get_gridid grtyp')
   call testutils_assert_eq((/ni,nj/),(/200,100/),'get_gridid nij')
   call testutils_assert_eq((/ig1,ig2,ig3,ig4/),(/0,0,0,0/),'get_gridid ig14')
!!$   ier = ezgxprm(gdid, ni, nj, grtyp, ig1, ig2, ig3, ig4, grref, ig1ref, ig2ref, ig3ref, ig4ref)
!!$   ier = ezgfstp(gdid, nomvarx, typvarx, etikx, nomvary, typvary, etiky, ip1, ip2, ip3, dateo, deet, npas, nbits)
!!$   ier = gdgaxes(gdid, ax, ay)
   istat = gdrls(gridid)

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   datevfuzz = 3600
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_near not found')

   datevfuzz = 3600*6

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_near')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_near value')

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LE)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_le')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_le value')

   datev_S = '20090427.020000'
   call datp2f(datev2,datev_S)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GE)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_ge not found')

   datev_S = '20090426.220000'
   call datp2f(datev2,datev_S)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GE)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_ge')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_ge value')


   datev_S = '20090427.000000'
   call datp2f(datev2,datev_S)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GT)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_gt not found')

   datev_S = '20090427.000000'
   call datp2f(datev2,datev_S)
   datevfuzz = 3600*6
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LT)
   call testutils_assert_ok(.not.RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_lt not found')

   datev_S = '20090427.000000'
   call datp2f(datev3,datev_S)
   nhours_8 = -40.D0/SEC_PER_HR
   call incdatr(datev2,datev3,nhours_8)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_gt as eq')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_gt as eq value')
   call datf2p(datev_S,datev2)


   datev_S = '20090427.000000'
   call datp2f(datev3,datev_S)
   nhours_8 = -1.D0
   call incdatr(datev2,datev3,nhours_8)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_GT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_gt')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_gt value')


   datev_S = '20090427.000000' 
   call datp2f(datev3,datev_S)
   nhours_8 = 40.D0/SEC_PER_HR
   call incdatr(datev2,datev3,nhours_8)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_lt as eq')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_lt as eq value')


   datev_S = '20090427.000000' 
   call datp2f(datev3,datev_S)
   nhours_8 = 1.D0
   call incdatr(datev2,datev3,nhours_8)
   key = fstmpi_find(funit,'TT',datev2,RMN_ANY_I,0,0,datevfuzz,FST_FIND_LT)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find_fuzz_lt')
   call testutils_assert_ok(datev2==datev,'test_fstmpi_find_read','fstmpi_find_fuzz_lt value')
   call datf2p(datev_S,datev2)

   istat = fstmpi_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_find_read','fstmpi_close')


   filename_S = trim(F_bcmk_S)//'geophy/Gem_geophy.fst'
   funit = fstmpi_open(filename_S,FST_READONLY)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fstmpi_find_read2','fstmpi_open')

   datev = RMN_ANY_DATE
   zp1 = 1.
   kind = RMN_CONV_ARBITRARY
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPOLD, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   key = fstmpi_find(funit,'J1',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find ip1>0 old')

   datev = RMN_ANY_DATE
   call convip_plus(ip1, zp1, kind, RMN_CONV_P2IPNEW, dummy_S, .not.RMN_CONV_USEFORMAT_L)
   key = fstmpi_find(funit,'J1',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find ip1>0 new')

   datev = RMN_ANY_DATE
   ip1 = 1200
   key = fstmpi_find(funit,'ME',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find ip1=1200')

   datev = RMN_ANY_DATE
   ip1 = 0
   key = fstmpi_find(funit,'ME',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find ip1=0 for 1200')

   datev = RMN_ANY_DATE
   ip1 = 0
   key = fstmpi_find(funit,'MG',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find ip1=0')

   datev = RMN_ANY_DATE
   ip1 = 1200
   key = fstmpi_find(funit,'MG',datev,ip1,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_find_read','fstmpi_find ip1=1200 for 0')

   istat = fstmpi_close(funit)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_find_read2','fstmpi_close')

   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpi_find_read


!/@
subroutine test_fstmpi_write()
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use iso_c_binding
   use testutils
   use fstmpi_mod
   use ezgrid_mod
   use ptopo_utils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-01
   !@argument
!@/
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>
   include "rpn_comm.inc"
   real,parameter :: MYVALUE = 3.3
   integer,parameter :: NI0=50,NJ0=30,NK0=3,HALO=2
   character(len=256) :: nomvar_S,filename_S
   logical :: ok_L,ok2_L
   integer :: funit,istat,gridid,grididh,grididfull,gridid2,key,ig1,ig2,ig3,ig4,ip1,i,j,k,datev,ip1list(NK0)
   real,pointer :: data2d(:,:),data2dh(:,:),data3d(:,:,:),data3dh(:,:,:)
   real :: ax(1-HALO:NI0+HALO,1),ay(1,1-HALO:NJ0+HALO)
   ! ---------------------------------------------------------------------
   call ptopo_init_var()
   write(filename_S,'(a,I3.3)') '__test_fstmpi_to-rm__.fst-',ptopo_grid_ipe
   istat = clib_unlink(trim(filename_S))
   funit = fstmpi_open(filename_S)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fstmpi_write:open','')

   nomvar_S = 'ZX'

   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
   do i=1-HALO,NI0+HALO
      ax(i,1) = 10.+float(ptopo_bloc_ipex*NI0+i)*0.25
   enddo
   do j=1-HALO,NJ0+HALO
      ay(1,j) = float(ptopo_bloc_ipey*NJ0+j)*0.25
   enddo
   grididh = ezgdef_fmem(NI0+2*HALO,NJ0+2*HALO, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)
   gridid = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax(1:NI0,1), ay(1,1:NJ0))
   grididfull = ezgrid_merge(gridid,RPN_COMM_BLOC_COMM,.true.)

   allocate(data2d(NI0,NJ0),&
        data2dh(1-HALO:NI0+HALO,1-HALO:NJ0+HALO), &
        data3d(NI0,NJ0,NK0),&
        data3dh(1-HALO:NI0+HALO,1-HALO:NJ0+HALO,NK0), &
        stat=istat)
   data2d = MYVALUE
   data2dh = MYVALUE
   data3d = MYVALUE
   data3dh = MYVALUE
   ip1 = 0
!!$   lvlid = 
   do k=1,NK0
      ip1list(k) = (k+2)*3
   enddo

   nomvar_S = 'ZX2d'
   istat = fstmpi_write(funit,nomvar_S,data2d,gridid,ip1,F_npak=FST_NPAK_FULL32)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:write_2d_r4','')

   nomvar_S = 'ZH2d'
   istat = fstmpi_write(funit,nomvar_S,data2dh,grididh,ip1,F_npak=FST_NPAK_FULL32,F_lni=NI0,F_lnj=NJ0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:write_2d_r4','')

   nomvar_S = 'ZX3d'
   istat = fstmpi_write(funit,nomvar_S,data3d,gridid,ip1list,F_npak=FST_NPAK_FULL32)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:write_3d_r4','')

   nomvar_S = 'ZH3d'
   istat = fstmpi_write(funit,nomvar_S,data3dh,grididh,ip1list,F_npak=FST_NPAK_FULL32,F_lni=NI0,F_lnj=NJ0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:write_3d_r4','')
   !TODO: test with a vgrid instead of ip1list

   if (funit > 0) istat = fstmpi_close(funit)
   deallocate(data2d,data2dh,data3d,stat=istat)

   !- Checking

   funit = fstmpi_open(filename_S,FST_READONLY)
   call testutils_assert_ok(RMN_IS_OK(funit),'test_fstmpi_write:open','')

   nomvar_S = 'ZX2d'
   datev = RMN_ANY_DATE
   key = fstmpi_find(funit,nomvar_S,datev,RMN_ANY_I,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_write:find','2d')

   istat = fstmpi_read(key,data3d,funit,gridid2)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:read','2d')
   ok_L = .false. ; ok2_L = .false.
   if (associated(data3d)) then
      ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0,ptopo_bloc_npey*NJ0,1/))
      ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
   endif
   call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data2d shape')
   call testutils_assert_ok(ok2_L,'test_fstmpi_write:read','data2d')
   if (.not.ok2_L) then
      print *,'data2d min,max:',minval(data3d),maxval(data3d)
   endif
   ok_L = ezgrid_samegrid(grididfull,gridid2)
   call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data2d grid')


   if (associated(data3d)) deallocate(data3d,stat=istat)
   nomvar_S = 'ZH2d'
   datev = RMN_ANY_DATE
   key = fstmpi_find(funit,nomvar_S,datev,RMN_ANY_I,RMN_ANY_I,RMN_ANY_I)
   call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_write:find','2dh')

   istat = fstmpi_read(key,data3d,funit,gridid2)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:read','2dh')
   ok_L = .false. ; ok2_L = .false.
   if (associated(data3d)) then
      ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0,ptopo_bloc_npey*NJ0,1/))
      ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
   endif
   call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data2dh shape')
   call testutils_assert_ok(ok2_L,'test_fstmpi_write:read','data2dh')
   if (.not.ok2_L) then
      print *,'data2d min,max:',minval(data3d),maxval(data3d)
   endif
!!$   ok_L = ezgrid_samegrid(grididfull,gridid2)
!!$   call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data2dh grid')


   nomvar_S = 'ZX3d'
   do k=1,NK0
      if (associated(data3d)) deallocate(data3d,stat=istat)
      datev = RMN_ANY_DATE
      key = fstmpi_find(funit,nomvar_S,datev,ip1list(k),RMN_ANY_I,RMN_ANY_I)
      call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_write:find','3d')

      istat = fstmpi_read(key,data3d,funit,gridid2)
      call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:read','3d')
      ok_L = .false. ; ok2_L = .false.
      if (associated(data3d)) then
         ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0,ptopo_bloc_npey*NJ0,1/))
         ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
      endif
      call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data3d shape')
      call testutils_assert_ok(ok2_L,'test_fstmpi_write:read','data3d')
      if (.not.ok2_L) then
         print *,'data3d min,max:',minval(data3d),maxval(data3d)
      endif
      ok_L = ezgrid_samegrid(grididfull,gridid2)
      call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data3d grid')
   enddo

   nomvar_S = 'ZH3d'
   do k=1,NK0
      if (associated(data3d)) deallocate(data3d,stat=istat)
      datev = RMN_ANY_DATE
      key = fstmpi_find(funit,nomvar_S,datev,ip1list(k),RMN_ANY_I,RMN_ANY_I)
      call testutils_assert_ok(RMN_IS_OK(key),'test_fstmpi_write:find','3dh')

      istat = fstmpi_read(key,data3d,funit,gridid2)
      call testutils_assert_ok(RMN_IS_OK(istat),'test_fstmpi_write:read','3dh')
      ok_L = .false. ; ok2_L = .false.
      if (associated(data3d)) then
         ok_L = all(shape(data3d)==(/ptopo_bloc_npex*NI0,ptopo_bloc_npey*NJ0,1/))
         ok2_L = all(abs(data3d-MYVALUE)<1.e-5)
      endif
      call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data3dh shape')
      call testutils_assert_ok(ok2_L,'test_fstmpi_write:read','data3dh')
      if (.not.ok2_L) then
         print *,'data3d min,max:',minval(data3d),maxval(data3d)
      endif
      ok_L = ezgrid_samegrid(grididfull,gridid2)
      call testutils_assert_ok(ok_L,'test_fstmpi_write:read','data3d grid')
   enddo

   if (funit > 0) istat = fstmpi_close(funit)

   if (associated(data2d)) deallocate(data2d,stat=istat)
   if (associated(data3d)) deallocate(data3d,stat=istat)

   istat = clib_unlink(trim(filename_S))
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstmpi_write
