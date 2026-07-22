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

      subroutine OUTs_server
      use iso_c_binding
      use app_mpmd
      use omp_timing
      use IOs
      use OUTs
      use clib_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=  4) :: me_S
      character(len=256) :: component_S
      logical :: process_L,signal
      integer :: myproc,tag,wm(4),ierr,flag,colors(2),COMMs(2),cnt
      integer :: recv(6),tr,nvar,nko,dim_meta,grid_info(100)
      integer :: steps(1000),imin,imax
      real(kind=REAL64) sum,sumd2,moy,var,fijk,npt_8,sume(1000),mind,maxd
!     
!--------------------------------------------------------------------
!
      component_S= 'OUT-server' ; colors= (/3,1/)
      call IOs_mpi_init (component_S, colors, COMMs, 2, OUTs_server_L)

      OUTs_GEM_COMM= COMMs(2)
      call IOs_gem_init (COMMs(2),OUTs_1o1_L,OUTs_me1o1,&
                         OUTs_gem1o1,grid_info)

      if (.not.OUTs_server_L) goto 999

      call OUTs_init (grid_info)
      
      write (me_S,'(I4.4)') myproc_IOS
      process_L= .true.
      if (Lun_out>0) call clock ( Lun_out, 'Init completed', .false. )

      do while (process_L)
 99      if (Lun_out>0) call clock ( Lun_out, 'Waiting for instructions', .false. )

         call gtmg_start ( 10, 'Waiting', 1)
         if (OUTs_1o1_L) then
            tag=601
            call MPI_recv ( recv,6, MPI_INTEGER, gem_1on1, &
                            tag, OUTs_GEM_COMM, MPI_STATUSES_IGNORE, ierr)
            nout_files = recv(1)
            if (nout_files>0) then
               tag=tag+1
               call MPI_recv ( Out_dirname_S,len(Out_dirname_S),&
                  MPI_CHARACTER,gem_1on1, tag, OUTs_GEM_COMM, MPI_STATUSES_IGNORE,ierr )
               tag=tag+1
               call MPI_recv ( Out_filenames_S                   ,&
                            nout_files*len(Out_filenames_S(1)),&
                     MPI_CHARACTER,gem_1on1,tag,OUTs_GEM_COMM,MPI_STATUSES_IGNORE,ierr )
              if (Lun_out>0) call clock ( Lun_out, 'Waiting for data', .false. )
            endif
         endif
         call MPI_bcast (recv, size(recv),MPI_INTEGER,0,MY_WORLD_COMM,ierr)
         call gtmg_stop ( 10 )

         nout_files = recv(1)
         if (nout_files==0) goto 99
         if (nout_files< 0) goto 777
         IOS_events= IOS_events+1
         
         call gtmg_start ( 11, 'FST_OPEN', 1)
         Out_ip2    = recv(2)
         Out_ip3    = recv(3)
         Out_npas   = recv(4)
         Out_typvar_S = TRANSFER(recv(5),chac1)
         Out_dateo  = recv(6)
         cnt= len(Out_dirname_S)
         call MPI_bcast (Out_dirname_S  ,cnt,MPI_CHARACTER,0,&
                         MY_WORLD_COMM,ierr)
         cnt= nout_files*len(Out_filenames_S(1))
         call MPI_bcast (Out_filenames_S,cnt,MPI_CHARACTER,0,&
                         MY_WORLD_COMM,ierr)

         call OUTs_openFST (me_S)
         call MPI_barrier (MY_WORLD_COMM,ierr)
         
         call gtmg_stop ( 11 )
         if (Lun_out>0) call clock ( Lun_out, 'FST files', .false. )
         steps(IOS_events) = Out_npas
         
         if (OUTs_1o1_L) then
            tag=901
            call MPI_recv ( wm, size(wm), MPI_INTEGER, gem_1on1, &
                            tag, OUTs_GEM_COMM, MPI_STATUSES_IGNORE, ierr)
         endif
         call MPI_bcast (wm, size(wm), MPI_INTEGER, 0, MY_WORLD_COMM, ierr)
         dim_meta= wm(2)-wm(1)+1

         if (dim_meta>0) then
            allocate (metaG(dim_meta))
            if (OUTs_1o1_L) then
               tag=tag+1
               call MPI_recv ( metaG, size(metaG), MPI_INTEGER,&
               gem_1on1, tag, OUTs_GEM_COMM, MPI_STATUSES_IGNORE, ierr)
            endif
            call MPI_bcast (metaG, size(metaG), MPI_INTEGER, 0, MY_WORLD_COMM, ierr)
            call gtmg_start ( 13, 'Process', 1)
            call OUTs_process_data (wm(1), wm(2), wm(4))
            call gtmg_stop ( 13 )
            call MPI_barrier (MY_WORLD_COMM,ierr)
            deallocate (metaG)
         endif
         
         call OUTs_closFST (wm(3))
         
         call MPI_barrier (MY_WORLD_COMM,ierr)

         goto 99
      end do

 777  call MPI_Win_free (SHARED_WIN,ierr)
      call gtmg_stop ( 1 )
 999  call gtmg_terminate3( myproc_IOS, MY_WORLD_COMM, trim(component_S))
      if (Lun_out>0) then
         call clock ( Lun_out, 'TERMINATING OUT-server', .true. )
         write(Lun_out,'(80("#")/)')
      endif
      call App_MPMD_Finalize()
      call MPI_FINALIZE(ierr)
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_server
