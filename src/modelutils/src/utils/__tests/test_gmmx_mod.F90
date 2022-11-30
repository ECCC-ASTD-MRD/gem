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

#include <rmn/msg.h>

!/@
subroutine test_gmmx_mod()
   use testutils
   use gmmx_mod
   use hgrid_wb
   use vgrid_wb
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-01
!@/
#include <rmnlib_basics.hf>

   integer,parameter :: NI0=4
   integer,parameter :: NJ0=5
   integer,parameter :: NK0=3
   character(len=512) :: name_S,meta_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S
   logical :: ok_L
   integer :: istat,ig1,ig2,ig3,ig4,gid_l1,i,j
   integer,target:: ip1list(NK0)
   integer,pointer:: ip1list_p(:)
   real :: ax(NI0,1),ay(1,NJ0)
   real,pointer :: ptr3d(:,:,:),ptr2d(:,:)
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_gmmx_mod')

   !- Fail Creation

   istat = gmmx_new(' ')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, no meta')

   istat = gmmx_new('hgrid=a ; vgrid=b')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, missing VN meta')

   istat = gmmx_new('vn=my ; hgrid=a')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, missing Vgrid meta')

   istat = gmmx_new(' ',' ',' ')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, no name')
   
   istat = gmmx_new('MYVAR',' ','vgrid')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, no hgrid')

   istat = gmmx_new('MYVAR','hgrid',' ')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, no vgrid')

   istat = gmmx_new('MYVAR','hgrid','vgrid')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, non existant h/vgrid')

   !- Creation

   do i=1,NI0
      ax(i,1) = float(i)
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)
   enddo
   ig1=0 ; ig2=0 ; ig3=100 ; ig4=100 
   gid_l1 = ezgdef_fmem(NI0,NJ0, 'L', ' ', ig1,ig2,ig3,ig4, ax, ay)
   istat = hgrid_wb_put('localz',gid_l1)

   ip1list = (/123,234,345/)
   ip1list_p => ip1list
   istat = vgrid_wb_put('myvgrid',VGRID_GROUND_TYPE,ip1list_p)

   istat = gmmx_new('vn=MYVARr43d0 ; hgrid=localz ; vgrid=myvgrid')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, r4 3d from string')

   istat = gmmx_new('MYVARr43d','localz','myvgrid')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, r4 3d')

   istat = gmmx_new('MYVARr43d2','localz','myvgrid',F_meta_S='in=i2 ; on=o2')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, r4 3d2')

   istat = gmmx_new('MYVARr43d3','localz','myvgrid',F_meta_S='in=i3 ; on=o3')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, r4 3d3')

   istat = gmmx_new('MYVARr43d4','localz','myvgrid','g',F_meta_S='in=i4 ; on=o4')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, r4 3d4')

   istat = gmmx_new('MYVARr43d5','localz','myvgrid','g',F_meta_S='in=i5 ; on=o5 ; mozaic=4')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, r4 3d extra param')

   istat = gmmx_new('MYVARbus1','localz','myvgrid','P')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, bus1')

   istat = gmmx_new('MYVARbus2','localz','myvgrid','P')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'new, bus2')

   !- Fail Creation

   istat = gmmx_new('MYVARr43d5','localz','myvgrid','g')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, duplicate name gmm')

   istat = gmmx_new('MYVARbus2','localz','myvgrid','P')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'new, duplicate name bus')

   !TODO: test duplicate in/out name

   !- Fail Meta Queries

   name_S = '__NO_SUCH_VAR__'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'meta, no such name')

   name_S = ' '
   in_name_S = 'zz'
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'meta, no such IN')

   name_S = 'MYVARr43d'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = 'P'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'meta wrong VB')

   name_S = 'MYVARr43d2'
   in_name_S = 'i3'
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'meta wrong IN')

   name_S = ' '
   in_name_S = ' '
   out_name_S = 'o3'
   busname_S = 'P'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'meta unmatched IN/VB')

  !- Meta Queries

   name_S = 'MYVARr43d'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta w/ name')
   call testutils_assert_eq(busname_S,'g','meta w/ name VB')
   call testutils_assert_eq(hgrid_S,'localz','meta w/ name HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','meta w/ name VG')

   name_S = 'MYVARr43d3'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta w/ name 2')
   call testutils_assert_eq(busname_S,'g','meta w/ name 2 VB')
   call testutils_assert_eq(hgrid_S,'localz','meta w/ name 2 HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','meta w/ name 2 VG')
   call testutils_assert_eq(in_name_S,'i3','meta w/ name 2 IN')
   call testutils_assert_eq(out_name_S,'o3','meta w/ name 2 ON')

   name_S = 'MYVARr43d3'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = 'g'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta w/ name 3 and VB')
   call testutils_assert_eq(busname_S,'g','meta w/ name 3 VB')
   call testutils_assert_eq(hgrid_S,'localz','meta w/ name 3 HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','meta w/ name 3 VG')
   call testutils_assert_eq(in_name_S,'i3','meta w/ name 3 IN')
   call testutils_assert_eq(out_name_S,'o3','meta w/ name 3 ON')

   name_S = ' '
   in_name_S = ' '
   out_name_S = 'o3'
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta w/ ON')
   call testutils_assert_eq(name_S,'myvarr43d3','meta w/ ON name')
   call testutils_assert_eq(busname_S,'g','meta w/ ON VB')
   call testutils_assert_eq(hgrid_S,'localz','meta w/ ON HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','meta w/ ON VG')
   call testutils_assert_eq(in_name_S,'i3','meta w/ ON IN')
   call testutils_assert_eq(out_name_S,'o3','meta w/ ON ON')

   name_S = ' '
   in_name_S = 'i2'
   out_name_S = ' '
   busname_S = 'g'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta w/ IN/VB')
   call testutils_assert_eq(name_S,'myvarr43d2','meta w/ IN/VB name')
   call testutils_assert_eq(busname_S,'g','meta w/ IN/VB VB')
   call testutils_assert_eq(hgrid_S,'localz','meta w/ IN/VB HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','meta w/ IN/VB VG')
   call testutils_assert_eq(in_name_S,'i2','meta w/ IN/VB IN')
   call testutils_assert_eq(out_name_S,'o2','meta w/ IN/VB ON')

   name_S = 'MYVARbus1'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta name bus1')

   out_name_S = 'MYVARbus2'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = 'P'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = gmmx_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'meta name bus2')

   !TODO: more complex bus tests

   !- Failed Data Queries
   
   nullify(ptr3d)
   name_S = '__NO_SUCH_VAR__'
   istat = gmmx_data(name_S,ptr3d)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'meta, no such name')

   nullify(ptr2d)
   name_S = 'MYVARr43d'
   istat = gmmx_data(name_S,ptr2d)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'meta, wrong rank')

   nullify(ptr2d)
   name_S = 'MYVARbus2'
   istat = gmmx_data(name_S,ptr2d)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'meta, wrong VB')

   !- Data Queries
   nullify(ptr3d)
   name_S = 'MYVARr43d'
   istat = gmmx_data(name_S,ptr3d)
   call testutils_assert_eq(RMN_IS_OK(istat),.true., 'meta, gmm ptr3d')
   call testutils_assert_eq(associated(ptr3d),.true., 'meta, gmm ptr3d associated')
   if (associated(ptr3d)) then
      call testutils_assert_eq(shape(ptr3d),(/NI0,NJ0,NK0/), 'meta, gmm ptr3d shape')
   endif

   !TODO: gmmx_list tests

   !TODO: bus tests
   ! ---------------------------------------------------------------------
   return
end subroutine test_gmmx_mod
