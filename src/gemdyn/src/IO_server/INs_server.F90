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
      subroutine INs_server
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use omp_timing
      use IOs
      use INs
      implicit none

      logical, external :: INs_open
      
      character(len=256) :: component_S
      character(len=16 ) :: dateV_S
      logical :: process_L
      integer :: colors(2),COMMs(2),ierr,grid_info(100),tag,client
!     
!--------------------------------------------------------------------
!
      component_S= 'IN-server' ; colors= (/4,1/)
      call IOs_mpi_init (component_S, colors, COMMs, 2, INs_server_L)

      INs_GEM_COMM= COMMs(2)
      call IOs_gem_init (COMMs(2),INs_1o1_L,INs_me1o1,&
                         INs_gem1o1,grid_info)
      if (.not.INs_server_L) goto 999
                         
      call INs_init ()

      process_L= .true.
      if (Lun_out>0) call clock ( Lun_out, 'Init completed', .false. )

      do while (process_L)
 99      if (Lun_out>0) call clock ( Lun_out, 'Waiting for instructions', .false. )
         
         call gtmg_start ( 5, 'Ins_waiting', 1)
         call INs_wait (dateV_S, client)
         call gtmg_stop ( 5 )
         if (INs_nreq<0) goto 999
         if (INs_nreq<1) goto 99

         if (client == 1001) then
            call INs_Nes_Data ( dateV_S )
         else if (client == 2001) then
            call INs_Iau_Data  ( dateV_S )
         else if (client == 3001) then
            call INs_Spn_Data  ( dateV_S )
         endif

         goto 99
      end do
      
 999  call gtmg_stop ( 1 )
      call gtmg_terminate3( myproc_IOS, MY_WORLD_COMM, trim(component_S))
      if (Lun_out>0) then
         call clock ( Lun_out, 'TERMINATING IN-server', .true. )
         write(Lun_out,'(80("#")/)')
      endif
      call MPI_FINALIZE(ierr)
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_server
