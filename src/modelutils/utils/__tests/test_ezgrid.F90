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
subroutine test_ezgrid()
   use, intrinsic :: iso_fortran_env, only: INT64
   use testutils
   use ezgrid_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   integer,parameter :: NI0=11,NJ0=13
   integer :: istat,ni,nj,ig1,ig2,ig3,ig4,ig1ref,ig2ref,ig3ref,ig4ref, &
        i0,in,i,j,gid_l1,gid_l2,nij(2),ig14(4),igp14(4),ig14_l2_orig(4),ij0(2),gid_z1,gid_z2,gid_z1b
   real :: ax(NI0,1),ay(1,NJ0),lat0,lon0,dlat,dlon
   real,pointer :: ax2(:,:),ay2(:,:),lat(:,:),lon(:,:)
   character(len=2) :: grtyp_S,grref_S
   integer,pointer :: griddata(:)
   logical :: ok_L,ok2_L
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_ezgrid')

   do i=1,NI0
      ax(i,1) = float(i)
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)
   enddo

   lat0 = 45. ; lon0=270.; dlat=1. ; dlon=1.
   call cxgaig('L', ig1,ig2,ig3,ig4, lat0,lon0,dlat,dlon)
   gid_l1 = ezgdef_fmem(NI0,NJ0, 'L', ' ', ig1,ig2,ig3,ig4, ax, ay)
   call cxgaig('L', ig1,ig2,ig3,ig4, 30., 180., 1.5, 1.5)
   ig14_l2_orig = (/ig1,ig2,ig3,ig4/)
   gid_l2 = ezgdef_fmem(NI0,NJ0, 'L', ' ', ig1,ig2,ig3,ig4, ax, ay)

   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
   do i=1,NI0
      ax(i,1) = 10.+float(i)*0.25
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)*0.25
   enddo
   gid_z1 = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)
   ax = ax + 20.
   ay = ay + 5.
   gid_z1b = ezgdef_fmem(NI0,NJ0, 'Z', 'E', ig1,ig2,ig3,ig4, ax, ay)
   ig1=0 ; ig2=10 ; ig3=43200 ; ig4=43150
   gid_z2 = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)

   !- ezgrid_params: L
   istat  = ezgrid_params(gid_l2,nij,grtyp_S,grref_S,ig14,ij0)
   ok_L  = (grtyp_S(1:1)=='L' .and. all(ig14 == ig14_l2_orig) .and. all(ij0 == (/1,1/)) .and. all(nij == (/NI0,NJ0/)))
   call testutils_assert_ok(ok_L,'test_ezgrid','params L')

   !- ezgrid_params: Z/E
   nullify(ax2,ay2)
   istat  = ezgrid_params(gid_z2,nij,grtyp_S,grref_S,ig14,ij0,ax2,ay2,igp14)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'params Z/E status')
   call testutils_assert_eq(grtyp_S(1:1)//grref_S(1:1),'ZE','params Z/E grtyp_S//grref_S')
   call testutils_assert_eq(nij,(/NI0,NJ0/),'params Z/E nij')
   call testutils_assert_eq(ig14,(/0,10,43200,43150/),'params Z/E ig14')
   call testutils_assert_eq(ij0,(/1,1/),'params Z/E ij0')
   call testutils_assert_eq(ax2(:,1),ax(:,1),'params Z/E ax')
   call testutils_assert_eq(ay2(1,:),ay(1,:),'params Z/E ay')
!!$   print *,'igp14=',igp14

   !- ezgrid_latlon
   nullify(lat,lon)
   istat = ezgrid_latlon(gid_l1,lat,lon)
!!$   ok_L = (RMN_IS_OK(istat) .and. &
!!$        associated(lat) .and. associated(lon) .and. &
!!$        all((/minval(lat),maxval(lat),minval(lon),maxval(lon)/) == (/lat0,lat0+dlat*float(NJ0-1),lon0,lon0+dlon*float(NI0-1)/)))
   ok_L = (RMN_IS_OK(istat) .and. associated(lat) .and. associated(lon))
   call testutils_assert_ok(ok_L,'test_ezgrid','latlon')
!!$   if (.not.ok_L) then
!!$      print *,'Expected:',lat0,lat0+dlat*float(NJ0-1),lon0,lon0+dlon*float(NI0-1)
!!$      print *,'Found   :',minval(lat),maxval(lat),minval(lon),maxval(lon)
!!$   endif
   do i=1,NI0
      ax2(i,1) = lon0+float(i-1)*dlon
   enddo
   do j=1,NJ0
      ay2(1,j) = lat0+float(j-1)*dlat
   enddo
   call testutils_assert_eq(lat(1,:),ay2(1,:),'lat')
   call testutils_assert_eq(lon(:,1),ax2(:,1),'lon')
   deallocate(lat,lon,stat=istat)

   !- ezgrid_sameproj: L - L ok
   ok_L = ezgrid_sameproj(gid_l1,gid_l1)
   call testutils_assert_ok(ok_L,'test_ezgrid','sameproj(gid_l1,gid_l1)')
   ok_L = ezgrid_samegrid(gid_l1,gid_l1)
   call testutils_assert_ok(ok_L,'test_ezgrid','samegrid(gid_l1,gid_l1)')
   ok_L = ezgrid_sameproj(gid_l2,gid_l2)
   call testutils_assert_ok(ok_L,'test_ezgrid','sameproj(gid_l2,gid_l2)')
   ok_L = ezgrid_samegrid(gid_l2,gid_l2)
   call testutils_assert_ok(ok_L,'test_ezgrid','samegrid(gid_l2,gid_l2)')

   !- ezgrid_sameproj: L - L not
   ok_L = ezgrid_sameproj(gid_l1,gid_l2)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','sameproj(gid_l1,gid_l2)')
   ok_L = ezgrid_samegrid(gid_l1,gid_l2)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','samegrid(gid_l1,gid_l2)')

   !- ezgrid_sameproj: Z/E - Z/E ok
   ok_L = ezgrid_sameproj(gid_z1,gid_z1)
   call testutils_assert_ok(ok_L,'test_ezgrid','sameproj(gid_z1,gid_z1)')
   ok_L = ezgrid_samegrid(gid_z1,gid_z1)
   call testutils_assert_ok(ok_L,'test_ezgrid','samegrid(gid_z1,gid_z1)')

   ok_L = ezgrid_sameproj(gid_z1,gid_z1b)
   call testutils_assert_ok(ok_L,'test_ezgrid','sameproj(gid_z1,gid_z1b)')
   ok_L = ezgrid_samegrid(gid_z1,gid_z1b)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','samegrid(gid_z1,gid_z1b)')

   !- ezgrid_sameproj: Z/E - Z/E not
   ok_L = ezgrid_sameproj(gid_z1,gid_z2)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','sameproj(gid_z1,gid_z2)')
   ok_L = ezgrid_samegrid(gid_z1,gid_z2)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','samegrid(gid_z1,gid_z2)')

!!$   !- ezgrid_sameproj: L - Z/E ok
!!$   ok_L = ezgrid_sameproj(gid_l1,gid_z1)

   !- ezgrid_sameproj: L - Z/E not
   ok_L = ezgrid_sameproj(gid_l1,gid_z1)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','sameproj(gid_l1,gid_z1)')
   ok_L = ezgrid_samegrid(gid_l1,gid_z1)
   call testutils_assert_ok(.not.ok_L,'test_ezgrid','samegrid(gid_l1,gid_z1)')

!!$   !- ezgrid_samegrid: L - Z/E ok but not supported yet
!!$   ok_L = ezgrid_samegrid(gid_l1,gid_z1)
!!$
!!$   !- ezgrid_samegrid: L - Z/E not
!!$   ok_L = ezgrid_samegrid(gid_l1,gid_z2)
!!$   ok_L = ezgrid_samegrid(gid_l1,gid_z1b)


   !- ezgrid_serialize: L
   nullify(griddata)
   istat  = ezgrid_serialize(gid_l2,griddata)
   gid_l1 = ezgrid_unserialize(griddata)
   ok_L = ezgrid_samegrid(gid_l1,gid_l2)
   call testutils_assert_ok(ok_L,'test_ezgrid','serialize/unserialize L grid')

   !- ezgrid_serialize: Z/E
   nullify(griddata)
!!$   print *,'test_ezgrid 1' ; call flush(6)
   istat  = ezgrid_serialize(gid_z1,griddata)
!!$   print *,'test_ezgrid 2' ; call flush(6)
   gid_z2 = ezgrid_unserialize(griddata)
!!$   print *,'test_ezgrid 3' ; call flush(6)
   ok_L = ezgrid_samegrid(gid_z1,gid_z2)
!!$   print *,'test_ezgrid 4' ; call flush(6)
   call testutils_assert_ok(ok_L,'test_ezgrid','serialize/unserialize Z/E grid')

   !- ezgrid_sub
   gid_l1 = ezgrid_sub(gid_l2,-1,-1,999,999)
   call testutils_assert_ok(.not.RMN_IS_OK(gid_l1),'test_ezgrid','sub L grid error')

   gid_z2 = ezgrid_sub(gid_z1,-1,-1,999,999)
   ok_L = ezgrid_samegrid(gid_z1,gid_z2)
   call testutils_assert_ok(ok_L,'test_ezgrid','sub Z grid 0')

   gid_z2 = ezgrid_sub(gid_z1b,3,2,NI0-4,NJ0-5)
   nullify(ax2,ay2)
   istat  = ezgrid_params(gid_z2,nij,grtyp_S,grref_S,ig14,ij0,ax2,ay2)
   ok_L  = (grtyp_S(1:1)=='Z' .and. grref_S(1:1)=='E' .and. all(ig14 == (/900,0,43200,43100/)) .and. all(ij0 == (/1,1/)) .and. all(nij == (/NI0-6,NJ0-6/)))
   ok2_L  = (all(ax2 == ax(3:NI0-4,:)).and. all(ay2 == ay(:,2:NJ0-5)))
   call testutils_assert_ok(ok_L.and.ok2_L,'test_ezgrid','sub Z grid 1')
   if (.not.ok_L) then
      print *,grtyp_S(1:1),grref_S(1:1),ig14,ij0,nij
      print *,'ZE',(/900,0,43200,43100/),(/1,1/),(/NI0-6,NJ0-6/)
   endif
   if (.not.ok2_L) then
      print *,ax2
      print *,ax(3:NI0-4,:)
      print *,ay2
      print *,ay(:,2:NJ0-5)
   endif
   ! ---------------------------------------------------------------------
   return
end subroutine test_ezgrid
