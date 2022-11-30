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
subroutine test_ezgrid_mpi()
use iso_c_binding
   use testutils
   use ezgrid_mod
   use ptopo_utils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@/
#include <rmnlib_basics.hf>
#include <rmn/WhiteBoard.hf>
   include "rpn_comm.inc"

   integer,parameter :: NI0=11,NJ0=13,HALO=2
   integer :: istat,ni,nj,ig1,ig2,ig3,ig4,ig1ref,ig2ref,ig3ref,ig4ref, &
        i0,in,i,j,nij(2),ig14(4),ij0(2),&
        gid_l1,gid_z1,gid_l1b,gid_z1b,myproc,nblocx,nblocy
   real :: ax(NI0,1),ay(1,NJ0),axh(1-HALO:NI0+HALO,1),ayh(1,1-HALO:NJ0+HALO)
   real,pointer :: ax2(:,:),ay2(:,:)
   character(len=2) :: grtyp_S,grref_S
   integer,pointer :: griddata(:)
   logical :: ok_L
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   call ptopo_init_var()
   istat = wb_get('ptopo/nblocx',nblocx)
   istat = wb_get('ptopo/nblocy',nblocy)
   istat = ptopo_bloc_set(nblocx,nblocy)

   do i=1,NI0
      ax(i,1) = float(i)
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)
   enddo

   ig1=1 ; ig2=2 ; ig3=100 ; ig4=200 
   gid_l1 = ezgdef_fmem(NI0,NJ0, 'L', ' ', ig1,ig2,ig3,ig4, ax, ay)

   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
   do i=1,NI0
      ax(i,1) = 10.+float(i)*0.25
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)*0.25
   enddo
   gid_z1 = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)


   gid_l1b = ezgrid_bcast(gid_l1,RPN_COMM_BLOC_COMM)
   gid_z1b = ezgrid_bcast(gid_z1,RPN_COMM_BLOC_COMM)


   !- ezgrid_sameproj: L - L ok
   ok_L = ezgrid_sameproj(gid_l1,gid_l1b)
   call testutils_assert_ok(ok_L,'test_ezgrid','sameproj L')

   !- ezgrid_sameproj: Z/E - Z/E ok
   ok_L = ezgrid_sameproj(gid_z1,gid_z1b)
   call testutils_assert_ok(ok_L,'test_ezgrid','sameproj Z')

   !- ezgrid_params: L
   istat  = ezgrid_params(gid_l1b,nij,grtyp_S,grref_S,ig14,ij0)
   ok_L  = (grtyp_S(1:1)=='L' .and. all(ig14 == (/1,2,100,200/)) .and. all(ij0 == (/1,1/)) .and. all(nij == (/NI0,NJ0/)))
   call testutils_assert_ok(ok_L,'test_ezgrid_bcast','params L')
   if (.not.ok_L) then
      print *,'L : ',grtyp_S(1:1)
      print *,'1,2,100,200 : ',ig14
      print *,'1,1 : ',ij0
      print *,NI0,NJ0,' : ',nij
   endif

   !- ezgrid_params: Z/E
   nullify(ax2,ay2)
   istat  = ezgrid_params(gid_z1b,nij,grtyp_S,grref_S,ig14,ij0,ax2,ay2)
   ok_L  = (grtyp_S(1:1)=='Z' .and. grref_S(1:1)=='E' .and. all(ig14 == (/900,0,43200,43100/)) .and. all(ij0 == (/1,1/)) .and. all(nij == (/NI0,NJ0/)) .and. all(ax2 == ax).and. all(ay2 == ay))
   call testutils_assert_ok(ok_L,'test_ezgrid','params Z/E')
   if (.not.ok_L) then
      print *,'ZE : ',grtyp_S(1:1)//grref_S(1:1)
      print *,'900,0,43200,43100 : ',ig14
      print *,'1,1 : ',ij0
      print *,ax,' : '
      print *,ax2
      print *,ay,' : '
      print *,ay2
   endif

   !- ezgrid_merge
   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43200
   do i=1,NI0
      ax(i,1) = 10.+float(ptopo_bloc_ipex*NI0+i)*0.25
   enddo
   do j=1,NJ0
      ay(1,j) = float(ptopo_bloc_ipey*NJ0+j)*0.25
   enddo
   gid_z1 = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)

   gid_z1b = ezgrid_merge(gid_z1,RPN_COMM_BLOC_COMM)
   call testutils_assert_ok(RMN_IS_OK(gid_z1b),'test_ezgrid','merge bloc')
   ok_L = .true.
   if (ptopo_isblocmaster_L) then
      nullify(ax2,ay2)
      istat  = ezgrid_params(gid_z1b,nij,grtyp_S,grref_S,ig14,ij0,ax2,ay2)
      ok_L  = (grtyp_S(1:1)=='Z' .and. grref_S(1:1)=='E' .and. all(ig14 == (/900,0,43200,43200/)) .and. all(ij0 == (/1,1/)) .and. all(nij == (/ptopo_bloc_npex*NI0,ptopo_bloc_npey*NJ0/)))
      if (ok_L) then
         if (associated(ax2).and.associated(ay2)) then
            do i=1,ptopo_bloc_npex*NI0
               if (ax2(i,1) /= 10.+float(i)*0.25) then
                  print *,'ax2',i,ax2(i,1),10.+float(i)*0.25
                  ok_L = .false.
               endif
            enddo
            do j=1,ptopo_bloc_npey*NJ0
               if (ay2(1,j) /= float(j)*0.25)  then
                  print *,'ay2',j,ay2(1,j),float(j)*0.25
                  ok_L = .false.
               endif
            enddo
         else
            print *,'ax2,ay2: not associated'
            ok_L = .false.
         endif
      endif
   endif
   call testutils_assert_ok(ok_L,'test_ezgrid','merge bloc data')

   !- ezgrid_merge w/ halo
   call ptopo_init_var()
   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43200
   do i=1-HALO,NI0+HALO
      axh(i,1) = 10.+float(ptopo_bloc_ipex*NI0+i)*0.25
   enddo
   do j=1-HALO,NJ0+HALO
      ayh(1,j) = float(ptopo_bloc_ipey*NJ0+j)*0.25
   enddo
   gid_z1 = ezgdef_fmem(NI0+2*HALO,NJ0+2*HALO, 'Z',  'E', ig1,ig2,ig3,ig4, axh, ayh)

   gid_z1b = ezgrid_merge(gid_z1,RPN_COMM_BLOC_COMM,.false.,1+HALO,1+HALO,NI0,NJ0)
   call testutils_assert_ok(RMN_IS_OK(gid_z1b),'test_ezgrid','merge bloc h')
   ok_L = .true.
   if (ptopo_isblocmaster_L) then
      nullify(ax2,ay2)
      istat  = ezgrid_params(gid_z1b,nij,grtyp_S,grref_S,ig14,ij0,ax2,ay2)
      ok_L  = (grtyp_S(1:1)=='Z' .and. grref_S(1:1)=='E' .and. all(ig14 == (/900,0,43200,43200/)) .and. all(ij0 == (/1,1/)) .and. all(nij == (/ptopo_bloc_npex*NI0,ptopo_bloc_npey*NJ0/)))
      if (ok_L) then
         if (associated(ax2).and.associated(ay2)) then
            do i=1,ptopo_bloc_npex*NI0
               if (ax2(i,1) /= 10.+float(i)*0.25) then
                  print *,'ax2',i,ax2(i,1),10.+float(i)*0.25
                  ok_L = .false.
               endif
            enddo
            do j=1,ptopo_bloc_npey*NJ0
               if (ay2(1,j) /= float(j)*0.25)  then
                  print *,'ay2',j,ay2(1,j),float(j)*0.25
                  ok_L = .false.
               endif
            enddo
         else
            print *,'ax2,ay2: not associated'
            ok_L = .false.
         endif
      endif
   endif
   call testutils_assert_ok(ok_L,'test_ezgrid','merge bloc data h')


   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_ezgrid_mpi
