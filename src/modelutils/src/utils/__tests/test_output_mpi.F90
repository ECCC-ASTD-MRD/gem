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
subroutine test_output_mpi()
   use clib_itf_mod
   use testutils
   use output_mod
   use ptopo_utils
   use vGrid_Descriptors
   use vgrid_wb
   use hgrid_wb
   use fstmpi_mod
   use ezgrid_mod
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-09
   !@/
#include <rmnlib_basics.hf>

   real,parameter :: PRES_MB2PA = 100.
   integer, parameter :: MAX_ITEM = 99
   integer, parameter :: N_STEP   = 2
   integer, parameter :: NI   = 9
   integer, parameter :: NJ   = 7
   integer, parameter :: nksoil = 3
   integer :: nk   = 5
   integer, parameter :: HALO = 2
   integer :: my_id(3), istat, istat2,n_items, sum_items, istep,ivar,dateo,dt,myproc,resultdnitems(0:N_STEP),ip1,fileid,gridid,reduc_core(4)
   character(len=99) :: fname_S,msg_S
   character(len=4)  :: mylist_S(MAX_ITEM),resultd_S(MAX_ITEM,0:N_STEP),resultp_S(MAX_ITEM,0:N_STEP)
   logical :: ok_L
   integer,target :: ip1surf(1),ip1soil(nksoil)
   integer,pointer :: p_ip1(:)
   real,pointer :: data(:,:,:),datasoil(:,:,:),data2d(:,:),p0(:,:)
   character(len=32) :: hgrid_S,vgrid_S,dateo_S
   character(len=256) :: nomvar_S,filename_S,dirname_S
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   call testutils_verbosity()
   call testutils_set_name('test_output_mpi')
   call ptopo_init_var()

   resultdnitems = 0
   resultd_S = ' '
   resultd_S(1:11,1) = (/'s0','s3','g0','g3','a0','a1','a2','a3','a4','a5','a6'/)
   resultdnitems(1) = 11

   ip1 = 0
   dateo_S = '20090427.000000'
   call datp2f(dateo,dateo_S)
   dt = 1800
   reduc_core(1:4) = (/2,-2,2,-2/)
   fname_S = 'test_output_mpi.out'
   my_id(:) = 0
   my_id(1) = output_new(fname_S, 'd', dateo, dt,reduc_core)
   call testutils_assert_eq(RMN_IS_OK(minval(my_id(:))),.true.,'new, parse file')

   gridid = priv_defhgrid(hgrid_S,ptopo_grid_npex*NI,ptopo_grid_npey*NJ,NI,NJ,HALO)
   call handle_error(gridid,'test_output_mpi','Problem defining grid')
   call priv_defvgrid(vgrid_S,nk)
   ip1surf(1) = 0
   p_ip1 => ip1surf
   istat = vgrid_wb_put('vgrid_surf',VGRID_SURF_TYPE,p_ip1)
   ip1soil(1:nksoil) = (/1199,1198,1197/)
   p_ip1 => ip1soil
   istat = vgrid_wb_put('vgrid_soil',VGRID_GROUND_TYPE,p_ip1)

   allocate(data(1-HALO:NI+HALO,1-HALO:NJ+HALO,nk),stat=istat)
   allocate(datasoil(1-HALO:NI+HALO,1-HALO:NJ+HALO,nksoil),stat=istat2)
   istat = max(istat,istat2)
   allocate(data2d(1-HALO:NI+HALO,1-HALO:NJ+HALO),stat=istat2)
   istat = max(istat,istat2)
   allocate(p0(1-HALO:NI+HALO,1-HALO:NJ+HALO),stat=istat2)
   istat = -1 * max(istat,istat2)
   call handle_error(istat,'test_output_mpi','allocating mem')

   p0 = 1000.*PRES_MB2PA

   dirname_S = './__test_outout_mpi-to-delete__'
   istat = clib_mkdir(trim(dirname_S))
   istat = clib_chdir(trim(dirname_S))

   ok_L = .true.
   do istep = 0, N_STEP
      write(msg_S,'(a,i3.3)') ' step=',istep
      mylist_S = ' '
      n_items = output_getlist(my_id(1),istep,mylist_S)
      call testutils_assert_eq(n_items,resultdnitems(istep),'getlist'//trim(msg_S)//' nitems')
      call testutils_assert_eq(mylist_S,resultd_S(:,istep),'getlist'//msg_S)
      if (any(mylist_S(:) /= resultd_S(:,istep))) then
         print *,'ERROR: getlist,',istep,'d',n_items,mylist_S(1:n_items)
         ok_L = .false.
      endif

      if (istep == 1) then
         !TODO: change output dir
      endif

      !- controle file
      write(filename_S,'(a,i3.3,a,I3.3)') '__test_output_ctrl__.fst-',istep,'-',ptopo_grid_ipe
      fileid = fstmpi_open(filename_S)
      !-

      do ivar = 1,n_items
         call priv_fill(data,ivar,HALO)
         select case(mylist_S(ivar)(1:1))
            case('s')
               data2d = data(:,:,1)
!!$               print *,'data2d:',shape(data2d),':',lbound(data2d),':',ubound(data2d)

               istat = fstmpi_write(fileid,mylist_S(ivar),data2d,gridid,0,F_lni=NI,F_lnj=NJ)
               istat = output_writevar(my_id(1),istep,mylist_S(ivar),data2d,hgrid_S,'vgrid_surf')
               call testutils_assert_eq(RMN_IS_OK(istat),.true.,'writevar '//trim(mylist_S(ivar))//trim(msg_S))
            case('g')
               datasoil(:,:,1:nksoil) = data(:,:,1:nksoil)
               istat = output_writevar(my_id(1),istep,mylist_S(ivar),datasoil,hgrid_S,'vgrid_soil')
               call testutils_assert_eq(RMN_IS_OK(istat),.true.,'writevar '//trim(mylist_S(ivar))//trim(msg_S))
            case('a')
               istat = output_writevar(my_id(1),istep,mylist_S(ivar),data,hgrid_S,vgrid_S,p0)
               call testutils_assert_eq(RMN_IS_OK(istat),.true.,'writevar '//trim(mylist_S(ivar))//trim(msg_S))
         end select
!!$         if (.not.RMN_IS_OK(istat)) then
!!$            print *,'ERROR: writevar',istep,mylist_S(ivar)
!!$            ok_L = .false.            
!!$         endif
      enddo
      istat = output_close(my_id(1),istep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'close'//trim(msg_S))
      
      !-
      if (fileid > 0) istat = fstmpi_close(fileid)
      !-

      !TODO: call output_post_process_raise_flag

      !TODO: read data and check

      !TODO: clib_unlink created file

      !TODO: output_writestep
   enddo

   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return

contains

   !/@
   function priv_defhgrid(F_name_S,F_gni,F_gnj,F_lni,F_lnj,F_halo) result(F_gridid)
      character(len=*),intent(out) :: F_name_S
      integer,intent(in) :: F_gni,F_gnj,F_lni,F_lnj,F_halo
      integer :: F_gridid
      !@/
      integer :: istat,i,j,ig1,ig2,ig3,ig4,gid,i0,j0
      real :: ax(1-F_Halo:F_gni+F_halo,1),ay(1,1-F_halo:F_gnj+F_halo)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(test_output_mpi) priv_defhgrid')
      F_gridid = RMN_ERR
      F_name_S = 'h1'

      ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
      do i=1-F_halo,F_gni+F_halo
         ax(i,1) = 10.+float(i)*0.25
      enddo
      do j=1-F_halo,F_gnj+F_halo
         ay(1,j) = float(j-F_gnj/2)*0.25
      enddo
      gid = ezgdef_fmem(F_gni+2*F_halo,F_gnj+2*F_halo, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)

      i0 = 1 + ptopo_grid_ipex*F_lni
      j0 = 1 + ptopo_grid_ipey*F_lnj
      F_gridid = ezgrid_sub(gid,i0,j0,i0+F_lni-1+2*F_halo,j0+F_lnj-1+2*F_halo)

      i0 = i0 + F_halo
      j0 = j0 + F_halo
      istat = hgrid_wb_put(F_name_S,gid,i0,j0,F_lni,F_lnj,F_halo,F_halo)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_CRITICAL,'(test_output_mpi) hgrid_wb_put')
         return
      endif

      call msg(MSG_DEBUG,'(test_output_mpi) priv_defhgrid [end]')
      ! ---------------------------------------------------------------------
      return
   end function priv_defhgrid


   !/@
   subroutine priv_defvgrid(F_name_S,F_nk)
      character(len=*),intent(out) :: F_name_S
      integer,intent(inout) :: F_nk
      !@/
      real,parameter :: hyb(5) = (/0.33,0.55,0.73,0.88,0.99/)
      integer :: istat,i0,kk
      type(vgrid_descriptor) :: myvgrid
      integer,pointer :: ip1_m(:),ip1list(:)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(test_output_mpi) priv_defvgrid')
      F_name_S = 'v1'

      istat = vgd_new(myvgrid, &
           kind     = VGRID_HYBS_KIND, &
           version  = VGRID_HYBS_VER, &
           hyb      = hyb, &
           ptop_8   = 9575.0d0, &
           pref_8   = 100000.d0, &
           rcoef1   = 1., &
           rcoef2   = 1.)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_CRITICAL,'(test_output_mpi) vgd_new')
         return
      endif

      nullify(ip1_m)
      istat = vgd_get(myvgrid,key='VIPM',value=ip1_m)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_CRITICAL,'(test_output_mpi) vgd_get vip')
         return
      endif

      F_nk = min(size(ip1_m)-1,F_nk)
      allocate(ip1list(F_nk),stat=istat)
      do kk=1,F_nk
         ip1list(F_nk-kk+1) = ip1_m(kk+1)
      enddo

      istat = vgrid_wb_put(F_name_S,myvgrid,ip1list)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_CRITICAL,'(test_output_mpi) vgrid_wb_put')
         return
      endif
      call msg(MSG_DEBUG,'(test_output_mpi) priv_defvgrid [end]')
      ! ---------------------------------------------------------------------
      return
   end subroutine priv_defvgrid


   !/@
   subroutine priv_fill(F_data,F_ivar,F_halo)
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_ivar,F_halo
      !@/
      integer :: i,j,k,l_ijk(3),u_ijk(3),n,ni2,nj2
      real :: v0
      ! ---------------------------------------------------------------------
      l_ijk = lbound(F_data)
      u_ijk = ubound(F_data)
      ni2 = l_ijk(1) + (u_ijk(1) - l_ijk(1) + 1)/2
      nj2 = l_ijk(2) + (u_ijk(2) - l_ijk(2) + 1)/2
      v0 = float(ptopo_grid_ipe*F_ivar*maxval(u_ijk(1:2))*u_ijk(3))
      F_data = -F_ivar
      do k=l_ijk(3),u_ijk(3)
         do j=l_ijk(2)+F_halo,u_ijk(2)-F_halo
            do i=l_ijk(1)+F_halo,u_ijk(1)-F_halo
               n = max(abs(i-ni2),abs(j-nj2))
               F_data(i,j,k) = v0 + float(n*k*maxval(u_ijk(1:2)))
            enddo
         enddo
      enddo
      ! ---------------------------------------------------------------------
      return
   end subroutine priv_fill


end subroutine test_output_mpi
