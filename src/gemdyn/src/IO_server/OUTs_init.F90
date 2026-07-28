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

!**s/r OUTs_init

      subroutine OUTs_init (F_grid_info)
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use IOs
      use OUTs
      use UTL
      implicit none

      integer, intent(IN) :: F_grid_info(*)
      
#include <rmnlib_basics.hf>
      character(len=4) :: dumc4
      logical :: signal
      integer :: i,j,k,tag,myhost,DISP_UNIT,ierr,ind,ind2,nk,sortie_largest_list
      integer :: size_usr_grid_info,n_hgrid_usr,n_vgrid_usr,index,sizestr,sizei
      integer :: vgd_table_size(3), hgrid_param(9)
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      integer, dimension(100) :: usr_grid_info
      character(len=UTL_string_len_dic), dimension(:), pointer :: list_S
      integer, dimension(:), pointer :: ilist
      real :: dim
      type(C_PTR), save :: basepntr
!     
!--------------------------------------------------------------------
!
      nullify(list_S,ilist,OUT_usr_int_info)
      
! Allocate shared memory among server PEs of the same node
      call MPI_Comm_split_type(MY_WORLD_COMM, MPI_COMM_TYPE_SHARED,0,&
                               MPI_INFO_NULL, myhost, ierr)
      call MPI_COMM_rank (myhost,OUTs_hostmyproc ,ierr)

      disp_unit = 4
      dim = 0.
      ! Use variable sortie_largest_list just to clirifie code
      sortie_largest_list=F_grid_info(29)
      nk=max((G_nk+1)*2,sortie_largest_list)
      if (OUTs_hostmyproc == 0) dim = real(G_ni)*real(G_nj)*real(nk)&
                                     *real(IOS_ncolors)*real(disp_unit)
      WINSIZE = dim
      call MPI_Win_allocate_shared(WINSIZE, disp_unit, MPI_INFO_NULL,&
                                   myhost, basepntr, SHARED_WIN, ierr)
      call MPI_Win_shared_query(SHARED_WIN, MPI_PROC_NULL, WINSIZE  ,&
                                disp_unit, basepntr,ierr)
      call C_F_POINTER ( basepntr, IOs_glbdata, [G_ni,G_nj,nk,IOS_ncolors] )
      call MPI_barrier (MY_WORLD_COMM,ierr)

      Out_path_input_S= IOs_path_input_S
      Out_rewrit_L= TRANSFER(F_grid_info(24), signal)
      Out_deet = F_grid_info(25)
      Out_etik_S(1:12) = TRANSFER(F_grid_info(26), dumc4)//&
                         TRANSFER(F_grid_info(27), dumc4)//&
                         TRANSFER(F_grid_info(28), dumc4)

      if (OUTs_1o1_L) then
         tag=2000
         call MPI_recv (size_usr_grid_info, 1, MPI_INTEGER ,&
              OUTs_gem1o1, tag, OUTs_GEM_COMM, MPI_STATUSES_IGNORE, ierr)
         tag=tag+1
         call MPI_recv (usr_grid_info, size_usr_grid_info, MPI_INTEGER ,&
              OUTs_gem1o1, tag, OUTs_GEM_COMM, MPI_STATUSES_IGNORE, ierr)
      endif
      call MPI_bcast (size_usr_grid_info, 1, MPI_INTEGER, 0,&
           MY_WORLD_COMM, ierr)
      call MPI_bcast (usr_grid_info,size_usr_grid_info, MPI_INTEGER, 0,&
           MY_WORLD_COMM, ierr)
      
      ! Check if there are user vgrids
      n_vgrid_usr=usr_grid_info(1)
      n_hgrid_usr=usr_grid_info(2)
      if(n_vgrid_usr > 0)allocate(vgd_usr(n_vgrid_usr))
      index=3
      do i=1,n_vgrid_usr
         vgd_usr(i)%stag=usr_grid_info(index)
         nullify(vgd_usr(i)%vtbl_8,vgd_usr(i)%levels)
         ierr=vgd_free(vgd_usr(i)%vgd)
         index=index+1
      end do
      if(n_hgrid_usr > 0)allocate(hgd_usr(n_hgrid_usr))
      do i=1,n_hgrid_usr
         hgd_usr(i)%usr_grid_index=usr_grid_info(index)
         index=index+1
      end do

      if (OUTs_1o1_L) then
         if(n_hgrid_usr > 0)then
            if(.not.OUTS_read_hgrid_usr(hgd_usr))then
               write(Lun_out,'("OUTs_server, error in reading user horizontal grid, see above")')
               ! TODO handle error
            endif
         endif
         if(n_vgrid_usr > 0)then
            if(.not.OUTS_read_vgrid_usr(vgd_usr))then
               write(Lun_out,'("OUTs_server, error in reading user vertical grid, see above")')
               ! TODO handle error
            endif
         endif
      endif
      do i=1,n_hgrid_usr
         if (OUTs_1o1_L) then
            hgrid_param=(/hgd_usr(i)%ni,hgd_usr(i)%nj,&
                 hgd_usr(i)%ip1,hgd_usr(i)%ip2,hgd_usr(i)%ip3,&
                 hgd_usr(i)%ig1,hgd_usr(i)%ig2,hgd_usr(i)%ig3,hgd_usr(i)%ig4/)
         endif
         call MPI_bcast (hgrid_param, size(hgrid_param), MPI_INTEGER , 0, MY_WORLD_COMM, ierr)
         if(.not.OUTs_1o1_L)then               
            hgd_usr(i)%ni =hgrid_param(1); hgd_usr(i)%nj =hgrid_param(2)
            hgd_usr(i)%ip1=hgrid_param(3); hgd_usr(i)%ip2=hgrid_param(4); hgd_usr(i)%ip3=hgrid_param(5)
            hgd_usr(i)%ig1=hgrid_param(6); hgd_usr(i)%ig2=hgrid_param(7)
            hgd_usr(i)%ig3=hgrid_param(8); hgd_usr(i)%ig4=hgrid_param(9)
            allocate(hgd_usr(i)%tic(hgd_usr(i)%ni),hgd_usr(i)%tac(hgd_usr(i)%nj))
         end if
         call MPI_bcast (hgd_usr(i)%tic, size(hgd_usr(i)%tic), MPI_REAL , 0, MY_WORLD_COMM, ierr)
         call MPI_bcast (hgd_usr(i)%tac, size(hgd_usr(i)%tac), MPI_REAL , 0, MY_WORLD_COMM, ierr)
      end do
      do i=1,n_vgrid_usr
         if(OUTs_1o1_L)vgd_table_size=(/&
              size(vgd_usr(i)%vtbl_8,1),&
              size(vgd_usr(i)%vtbl_8,2),&
              size(vgd_usr(i)%vtbl_8,3)/)
         call MPI_bcast (vgd_table_size, size(vgd_table_size), MPI_INTEGER , 0, MY_WORLD_COMM, ierr)
         if(.not.OUTs_1o1_L)allocate(vgd_usr(i)%vtbl_8(vgd_table_size(1),vgd_table_size(2),vgd_table_size(3)))
         call MPI_bcast (vgd_usr(i)%vtbl_8, size(vgd_usr(i)%vtbl_8), MPI_DOUBLE_PRECISION , 0, MY_WORLD_COMM, ierr)
         if(.not.OUTs_1o1_L)ierr=vgd_new(vgd_usr(i)%vgd,vgd_usr(i)%vtbl_8)
      end do
      do i=1,n_vgrid_usr
         deallocate (vgd_usr(i)%vtbl_8)
         ierr=vgd_get(vgd_usr(i)%vgd, "VCOD", value = vgd_usr(i)%vcode)
         nullify(vgd_usr(i)%levels)
         ierr=vgd_get(vgd_usr(i)%vgd,'VCDM - vertical coordinate (m)',vgd_usr(i)%levels)
         ierr=vgd_get(vgd_usr(i)%vgd,'KIND',vgd_usr(i)%kind)
      end do
      do i=1,n_hgrid_usr
         hgd_usr(i)%ezgdid = ezgdef_fmem(hgd_usr(i)%ni,hgd_usr(i)%nj,'Z','E',&
              hgd_usr(i)%ig1,hgd_usr(i)%ig2,hgd_usr(i)%ig3,hgd_usr(i)%ig4,&
              hgd_usr(i)%tic,hgd_usr(i)%tac)
         if(Grd_yinyang_L)then
            hgd_usr(i)%same_rotation_L=.false.
         else
            hgd_usr(i)%same_rotation_L= &
                 hgd_usr(i)%ig1 == Rot_ig1 .and. &
                 hgd_usr(i)%ig2 == Rot_ig2 .and. &
                 hgd_usr(i)%ig3 == Rot_ig3 .and. &
                 hgd_usr(i)%ig4 == Rot_ig4
         endif
      end do      

      if(n_hgrid_usr > 0)then
         ! Read interpolation_dictionary     
         sizestr=UTL_string_len_dic/4
         if (OUTs_1o1_L)then
            if(.not.UTL_read_usr_int_dic(list_S,trim(Out_path_input_S)//'/MODEL_INPUT/interpolation_dictionary.txt',F_verbose_L=.false.))then
               print*,'Error in OUTs_init reading user interpolation dirctonary'
               ! TODO handle error gracefully
            endif
            if(mod(UTL_string_len_dic,4) /= 0)then
               print*,'Error in OUTs_init, UTL_string_len_dic must be a multible of 4'
               ! TODO handle error gracefully
               return
            endif
            sizei=sizestr*size(list_S)
            allocate(ilist(sizei))
            do i=1,size(list_S)
               !print*,'from string, list_S(i)',list_S(i)
               do j=1,sizestr
                  ind=(i-1)*sizestr+j
                  ind2=(j-1)*sizestr+1
                  ilist(ind)=TRANSFER(list_S(i)(ind2:ind2+3),k)
               enddo
            end do
         end if         
         call MPI_bcast (sizei, 1, MPI_INTEGER , 0, MY_WORLD_COMM, ierr)
         if (.not.OUTs_1o1_L)then
            allocate(ilist(sizei),list_S(sizei/sizestr))
         endif
         call MPI_bcast (ilist, size(ilist), MPI_INTEGER , 0, MY_WORLD_COMM, ierr)
         do i=1,size(list_S)
            do j=1,sizestr
               ind=(i-1)*sizestr+j
               ind2=(j-1)*sizestr+1
               list_S(i)(ind2:ind2+3)=TRANSFER(ilist(ind),dumc4)
            end do
         end do
         allocate(OUT_usr_int_info(size(list_S)))
         if(.not.UTL_crack_usr_int_dic(OUT_usr_int_info,list_S,'interpolation_dictionary.txt',F_verbose_L=.false.))then
            print*,'Error cracking user interpolation dirctonary'
            error stop 1
         endif
      endif
      Out_previous_usr_grid_S=' '
      Out_previous_src_ezgdid=-1
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_init

