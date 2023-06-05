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
subroutine test_dict_mod()
   use testutils
   use dict_mod
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
   call testutils_set_name('test_dict_mod')

   !- Fail Creation

   name_S = ''
   meta_S = ''
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'newvar, no name')
   
   name_S = 'MYVAR'
   meta_S = ''
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'newvar, no meta')

   name_S = 'MYVAR2'
   meta_S = 'HGRID=localz'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'newvar, missing meta')

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

   name_S = 'MYVARr43d'
   meta_S = 'HGRID=localz; VGRID=myvgrid'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, r4 3d')

   name_S = 'MYVARr43d2'
   meta_S = 'HGRID=localz; VGRID=myvgrid; in=i2 ; on=o2'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, r4 3d2')

   name_S = 'MYVARr43d3'
   meta_S = 'HGRID=localz; VGRID=myvgrid; in=i3 ; on=o3; VB=G'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, r4 3d3')

   name_S = 'MYVARr43d4'
   meta_S = 'HGRID=localz; VGRID=myvgrid; in=i4 ; on=o4; VB=G'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, r4 3d4')

   name_S = 'MYVARr43d5'
   meta_S = 'VN=myvn; HGRID=localz; VGRID=myvgrid; in=i4 ; on=o4; VB=G'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, r4 3d extra param')

   name_S = 'MYVARbus1'
   meta_S = 'HGRID=localz; VGRID=myvgrid; VB=P'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, bus1')

   name_S = 'MYVARbus2'
   meta_S = 'HGRID=localz; VGRID=myvgrid; VB=P'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'newvar, bus2')

   !- Fail Creation

   name_S = 'MYVARbus2'
   meta_S = 'HGRID=localz; VGRID=myvgrid; VB=P'
   istat = dict_newvar(name_S,meta_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'newvar, duplicate name')

   !TODO: test duplicate in/out name

   !- Fail Meta Queries

   name_S = '__NO_SUCH_VAR__'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'get_meta, no such name')

   name_S = ' '
   in_name_S = 'zz'
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'get_meta, no such IN')

   name_S = 'MYVARr43d'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = 'P'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'get_meta wrong VB')

   name_S = 'MYVARr43d2'
   in_name_S = 'i3'
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'get_meta wrong IN')

   name_S = ' '
   in_name_S = ' '
   out_name_S = 'o3'
   busname_S = 'P'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'get_meta unmatched IN/VB')

  !- Meta Queries

   name_S = 'MYVARr43d'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta w/ name')
   call testutils_assert_eq(busname_S,'g','get_meta w/ name VB')
   call testutils_assert_eq(hgrid_S,'localz','get_meta w/ name HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','get_meta w/ name VG')

   name_S = 'MYVARr43d3'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta w/ name 2')
   call testutils_assert_eq(busname_S,'g','get_meta w/ name 2 VB')
   call testutils_assert_eq(hgrid_S,'localz','get_meta w/ name 2 HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','get_meta w/ name 2 VG')
   call testutils_assert_eq(in_name_S,'i3','get_meta w/ name 2 IN')
   call testutils_assert_eq(out_name_S,'o3','get_meta w/ name 2 ON')

   name_S = 'MYVARr43d3'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = 'g'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta w/ name 3 and VB')
   call testutils_assert_eq(busname_S,'g','get_meta w/ name 3 VB')
   call testutils_assert_eq(hgrid_S,'localz','get_meta w/ name 3 HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','get_meta w/ name 3 VG')
   call testutils_assert_eq(in_name_S,'i3','get_meta w/ name 3 IN')
   call testutils_assert_eq(out_name_S,'o3','get_meta w/ name 3 ON')

   name_S = ' '
   in_name_S = ' '
   out_name_S = 'o3'
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta w/ ON')
   call testutils_assert_eq(name_S,'myvarr43d3','get_meta w/ ON name')
   call testutils_assert_eq(busname_S,'g','get_meta w/ ON VB')
   call testutils_assert_eq(hgrid_S,'localz','get_meta w/ ON HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','get_meta w/ ON VG')
   call testutils_assert_eq(in_name_S,'i3','get_meta w/ ON IN')
   call testutils_assert_eq(out_name_S,'o3','get_meta w/ ON ON')

   name_S = ' '
   in_name_S = 'i2'
   out_name_S = ' '
   busname_S = 'g'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta w/ IN/VB')
   call testutils_assert_eq(name_S,'myvarr43d2','get_meta w/ IN/VB name')
   call testutils_assert_eq(busname_S,'g','get_meta w/ IN/VB VB')
   call testutils_assert_eq(hgrid_S,'localz','get_meta w/ IN/VB HG')
   call testutils_assert_eq(vgrid_S,'myvgrid','get_meta w/ IN/VB VG')
   call testutils_assert_eq(in_name_S,'i2','get_meta w/ IN/VB IN')
   call testutils_assert_eq(out_name_S,'o2','get_meta w/ IN/VB ON')

   name_S = 'MYVARbus1'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = ' '
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta name bus1')

   out_name_S = 'MYVARbus2'
   in_name_S = ' '
   out_name_S = ' '
   busname_S = 'P'
   hgrid_S = ' '
   vgrid_S = ' '
   istat = dict_get_meta(name_S,in_name_S,out_name_S,busname_S,hgrid_S,vgrid_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'get_meta name bus2')

   !TODO: more complex bus tests

   !- Failed Data Queries
   
   nullify(ptr3d)
   name_S = '__NO_SUCH_VAR__'
   istat = dict_get_data(name_S,ptr3d)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'get_meta, no such name')

   nullify(ptr2d)
   name_S = 'MYVARr43d'
   istat = dict_get_data(name_S,ptr2d)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'get_meta, wrong rank')

   nullify(ptr2d)
   name_S = 'MYVARbus2'
   istat = dict_get_data(name_S,ptr2d)
   call testutils_assert_eq(RMN_IS_OK(istat),.false., 'get_meta, wrong VB')

   !- Data Queries
   nullify(ptr3d)
   name_S = 'MYVARr43d'
   istat = dict_get_data(name_S,ptr3d)
   call testutils_assert_eq(RMN_IS_OK(istat),.true., 'get_meta, gmm ptr3d')
   call testutils_assert_eq(associated(ptr3d),.true., 'get_meta, gmm ptr3d associated')
   if (associated(ptr3d)) then
      call testutils_assert_eq(shape(ptr3d),(/NI0,NJ0,NK0/), 'get_meta, gmm ptr3d shape')
   endif

   !TODO: bus tests
   ! ---------------------------------------------------------------------
   return
end subroutine test_dict_mod
