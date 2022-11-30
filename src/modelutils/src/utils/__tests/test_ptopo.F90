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
subroutine test_ptopo()
use iso_c_binding
   use testutils
   use ptopo_utils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-01
   !@/
#include <rmnlib_basics.hf>
#include <rmn/WhiteBoard.hf>
   include "rpn_comm.inc"

   logical,parameter :: PE0_ONLY = .true.
   integer,parameter :: NI0 = 3
   integer,parameter :: NJ0 = 5
   integer,parameter :: NK0 = 2
   integer,parameter :: HALO = 2
   integer :: myproc,istat,i,j,ngrids,igrid,npex,npey,nblocx,nblocy,mycol,myrow,i0,j0,bni,bnj,ib,jb
   integer :: nbloc,iblocx,iblocy,ibloc,bloc_npex,bloc_npey,bloc_npe,bloc_ipex,bloc_ipey,bloc_ipe,ipe_blocmaster
   real :: tmp_r
   real,pointer :: datain(:,:,:),dataout(:,:,:)
   real,pointer :: datain2d(:,:),dataout2d(:,:)
   logical :: ok_L
   ! ---------------------------------------------------------------------
!!$   myproc = testutils_initmpi()
   call testutils_verbosity()
   call testutils_set_name('test_ptopo')
!!$   istat = wb_get('ptopo/ngrids',ngrids)
!!$   istat = wb_get('ptopo/igrid',igrid)
!!$   istat = wb_get('ptopo/npx',npex)
!!$   istat = wb_get('ptopo/npy',npey)
!!$   istat = wb_get('ptopo/nblocx',nblocx)
!!$   istat = wb_get('ptopo/nblocy',nblocy)
!!$   istat = wb_get('ptopo/mycol',mycol)
!!$   istat = wb_get('ptopo/myrow',myrow)
   myproc = 0
   ngrids = 1
   igrid = 0
   npex = 1
   npey = 1
   nblocx = 1
   nblocy = 1
   mycol = 0
   myrow = 0

   bloc_npex = npex / nblocx
   bloc_npey = npey / nblocy
   bloc_npe = bloc_npex * bloc_npey

   nbloc = nblocx*nblocy
   iblocx = mycol/bloc_npex
   iblocy = myrow/bloc_npey
   ibloc = iblocx + iblocy * nblocx

   bloc_ipex = mycol - iblocx * bloc_npex
   bloc_ipey = myrow - iblocy * bloc_npey
   bloc_ipe = bloc_ipex + bloc_ipey * bloc_npex

   ipe_blocmaster = (iblocx * bloc_npex) + (iblocy * bloc_npey) * npex

   istat = RMN_ERR

   call ptopo_init_var()

   call testutils_assert_eq(ptopo_grid_npe,npex*npey,'ptopo_grid_npe')
   call testutils_assert_eq(ptopo_grid_npex,npex,'ptopo_grid_npex')
   call testutils_assert_eq(ptopo_grid_npey,npey,'ptopo_grid_npey')
   call testutils_assert_eq(ptopo_grid_ipe,myproc,'ptopo_grid_ipe')
   call testutils_assert_eq(ptopo_grid_ipex,mycol,'ptopo_grid_ipex')
   call testutils_assert_eq(ptopo_grid_ipey,myrow,'ptopo_grid_ipey')

   call testutils_assert_eq(ptopo_grid_nbloc,nbloc,'ptopo_grid_nbloc')
   call testutils_assert_eq(ptopo_grid_nblocx,nblocx,'ptopo_grid_nblocx')
   call testutils_assert_eq(ptopo_grid_nblocy,nblocy,'ptopo_grid_nblocy')

   call testutils_assert_eq(ptopo_grid_ipe_blocmaster,ipe_blocmaster,'ptopo_grid_ipe_blocmaster')
   call testutils_assert_eq(ptopo_grid_ibloc,ibloc,'ptopo_grid_ibloc')
   call testutils_assert_eq(ptopo_grid_iblocx,iblocx,'ptopo_grid_iblocx')
   call testutils_assert_eq(ptopo_grid_iblocy,iblocy,'ptopo_grid_iblocy')

   call testutils_assert_eq(ptopo_bloc_npe,bloc_npe,'ptopo_bloc_npe')
   call testutils_assert_eq(ptopo_bloc_npex,bloc_npex,'ptopo_bloc_npex')
   call testutils_assert_eq(ptopo_bloc_npey,bloc_npey,'ptopo_bloc_npey')
   call testutils_assert_eq(ptopo_bloc_ipe,bloc_ipe,'ptopo_bloc_ipe')
   call testutils_assert_eq(ptopo_bloc_ipex,bloc_ipex,'ptopo_bloc_ipex')
   call testutils_assert_eq(ptopo_bloc_ipey,bloc_ipey,'ptopo_bloc_ipey')
   call testutils_assert_eq(ptopo_isblocmaster_L,(bloc_ipe==0),'ptopo_isblocmaster_L')

   istat = ptopo_collect_dims(RPN_COMM_WORLD,NI0,NJ0,bni,bnj,i0,j0)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_ptopo','collect_dims not possible')


   !- collect bloc 3d
   istat = ptopo_collect_dims(RPN_COMM_BLOC_COMM,NI0,NJ0,bni,bnj,i0,j0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_ptopo','collect_dims bloc')
   ok_L = (&
        bni == NI0 * ptopo_bloc_npex .and. &
        bnj == NJ0 * ptopo_bloc_npey  .and. &
        i0  == 1 + NI0 * ptopo_bloc_ipex.and. &
        j0  == 1 + NJ0 * ptopo_bloc_ipey)
   call testutils_assert_ok(ok_L,'test_ptopo','collect_dims bloc data')
   if (.not.ok_L) then
      print *,'bni,bnj,i0,j0',bni,bnj,i0,j0
   endif

   allocate(datain(-2:NI0+1,-3:NJ0+4,NK0),stat=istat)
   datain = 0.
   do j=1,NJ0
      do i=1,NI0
         datain(i,j,:) = float(i0-1+i)+float(j0-1+j)/100.
!!         print *,'in ',i,j,ptopo_bloc_ipex,ptopo_bloc_ipey,datain(i,j,1)
      enddo
   enddo

   nullify(dataout)
   if (ptopo_isblocmaster_L) allocate(dataout(bni,bnj,NK0),stat=istat)

   istat = ptopo_collect(dataout,datain,RPN_COMM_BLOC_COMM,i0,j0,NI0,NJ0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_ptopo','bloc_collect')
   ok_L = .true.
   if (ptopo_isblocmaster_L) ok_L = associated(dataout)
   call testutils_assert_ok(ok_L,'test_ptopo','bloc_collect associated')
   if (ptopo_isblocmaster_L .and. ok_L) ok_L = all(ubound(dataout)==(/bni,bnj,NK0/))
   call testutils_assert_ok(ok_L,'test_ptopo','bloc_collect dims')
   if (ptopo_isblocmaster_L .and. ok_L) then
      do j=1,bnj
         do i=1,bni
            tmp_r = float(i)+float(j)/100.
            if (abs(dataout(i,j,1)-tmp_r) >= 0.001) then
               print *,'out',i,j,dataout(i,j,1),tmp_r,abs(dataout(i,j,1)-tmp_r)
               ok_L = .false.
            endif
         enddo
      enddo
   endif
   call testutils_assert_ok(ok_L,'test_ptopo','bloc_collect data')

   if (associated(dataout)) deallocate(dataout,stat=istat)


   !- collect grid 3d
   istat = ptopo_collect_dims(RPN_COMM_GRID,NI0,NJ0,bni,bnj,i0,j0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_ptopo','collect_dims grid')
   ok_L = (&
        bni == NI0 * ptopo_grid_npex .and. &
        bnj == NJ0 * ptopo_grid_npey  .and. &
        i0  == 1 + NI0 * ptopo_grid_ipex.and. &
        j0  == 1 + NJ0 * ptopo_grid_ipey)
   call testutils_assert_ok(ok_L,'test_ptopo','collect_dims grid data')
   if (.not.ok_L) then
      print *,'bni,bnj,i0,j0',bni,bnj,i0,j0
   endif

   datain = 0.
   do j=1,NJ0
      do i=1,NI0
         datain(i,j,:) = float(i0-1+i)+float(j0-1+j)/100.
!!         print *,'in ',i,j,ptopo_bloc_ipex,ptopo_bloc_ipey,datain(i,j,1)
      enddo
   enddo

   if (ptopo_grid_ipe == RPN_COMM_MASTER) allocate(dataout(bni,bnj,NK0),stat=istat)

   istat = ptopo_collect(dataout,datain,RPN_COMM_GRID,i0,j0,NI0,NJ0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_ptopo','grid_collect')
   ok_L = .true.
   if (ptopo_grid_ipe == RPN_COMM_MASTER) ok_L = associated(dataout)
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect associated')

   if (ptopo_grid_ipe == RPN_COMM_MASTER .and. ok_L) ok_L = all(ubound(dataout)==(/bni,bnj,NK0/))
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect dims')
   if (ptopo_grid_ipe == RPN_COMM_MASTER .and. ok_L) then
      do j=1,bnj
         do i=1,bni
            tmp_r = float(i)+float(j)/100.
            if (abs(dataout(i,j,1)-tmp_r) >= 0.001) then
               print *,'out',i,j,dataout(i,j,1),tmp_r,abs(dataout(i,j,1)-tmp_r)
               ok_L = .false.
            endif
         enddo
      enddo
   endif
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect data')


   !- collect grid 2d no halo
   allocate(datain2d(1:NI0,1:NJ0),stat=istat)
   if (ptopo_grid_ipe == RPN_COMM_MASTER) allocate(dataout2d(bni,bnj),stat=istat)
   datain2d(:,:) = datain(1:NI0,1:NJ0,1)
   istat = ptopo_collect(dataout2d,datain2d,RPN_COMM_GRID,i0,j0,NI0,NJ0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_ptopo','grid_collect 2d')
   ok_L = .true.
   if (ptopo_grid_ipe == RPN_COMM_MASTER) ok_L = associated(dataout2d)
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect associated 2d')
   if (ptopo_grid_ipe == RPN_COMM_MASTER .and. ok_L) ok_L = all(ubound(dataout2d)==(/bni,bnj/))
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect dims 2d')
   if (ptopo_grid_ipe == RPN_COMM_MASTER .and. ok_L) then
      do j=1,bnj
         do i=1,bni
            tmp_r = float(i)+float(j)/100.
            if (abs(dataout2d(i,j)-tmp_r) >= 0.001) then
               print *,'out',i,j,dataout2d(i,j),tmp_r,abs(dataout2d(i,j)-tmp_r)
               ok_L = .false.
            endif
         enddo
      enddo
   endif
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect data 2d')


   !- collect grid 2d halo
   if (associated(dataout2d)) deallocate(dataout2d,stat=istat)
   if (associated(datain2d)) deallocate(datain2d,stat=istat)

   allocate(datain2d(1-HALO:NI0+HALO,1-HALO:NJ0+HALO),stat=istat)
   if (ptopo_grid_ipe == RPN_COMM_MASTER) allocate(dataout2d(bni,bnj),stat=istat)
   datain2d = 0.
   datain2d(1:NI0,1:NJ0) = datain(1:NI0,1:NJ0,1)
   istat = ptopo_collect(dataout2d,datain2d,RPN_COMM_GRID,i0,j0,NI0,NJ0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_ptopo','grid_collect 2dh')
   ok_l = .true.
   if (ptopo_grid_ipe == RPN_COMM_MASTER) ok_L = associated(dataout2d)
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect associated 2dh')
   if (ptopo_grid_ipe == RPN_COMM_MASTER .and. ok_L) ok_L = all(ubound(dataout2d)==(/bni,bnj/))
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect dims 2dh')
   if (ptopo_grid_ipe == RPN_COMM_MASTER .and. ok_L) then
      do j=1,bnj
         do i=1,bni
            tmp_r = float(i)+float(j)/100.
            if (abs(dataout2d(i,j)-tmp_r) >= 0.001) then
               print *,'out',i,j,dataout2d(i,j),tmp_r,abs(dataout2d(i,j)-tmp_r)
               ok_L = .false.
            endif
         enddo
      enddo
   endif
   call testutils_assert_ok(ok_L,'test_ptopo','grid_collect data 2dh')

   if (associated(dataout)) deallocate(dataout,stat=istat)
   if (associated(datain)) deallocate(datain,stat=istat)
   if (associated(dataout2d)) deallocate(dataout2d,stat=istat)
   if (associated(datain2d)) deallocate(datain2d,stat=istat)

   !- 
   if (ptopo_grid_ipe == RPN_COMM_MASTER) then
      if (nblocx == npex.or.nblocy == npey) then
         nblocx = 1
         nblocy = 1
      else
         nblocx = npex
         nblocy = npey
      endif
   else
      nblocx = -1
      nblocy = -1
   endif
   istat = ptopo_bloc_set(nblocx,nblocy,.true.)
   call testutils_assert_eq(.true.,RMN_IS_OK(istat),'ptopo_bloc_set status')
   call testutils_assert_eq((/ptopo_grid_nblocx,ptopo_grid_nblocy/),&
        (/nblocx,nblocy/),'ptopo_bloc_set values')

!!$   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_ptopo
