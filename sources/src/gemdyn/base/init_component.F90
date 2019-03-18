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

!**s/r init_component

      subroutine init_component()
      use clib_itf_mod
      use cpus_options
      use dcst
      use glb_ld
      use HORgrid_options
      use iso_c_binding
      use path
      use ptopo
      use step_options
      use tdpack
      use gem_timing
      implicit none
#include <arch_specific.hf>

      include "rpn_comm.inc"

      external init_ndoms, pe_zero_topo
      integer, external :: model_timeout_alarm

      character(len=256) :: my_dir
      integer :: ierr, mydomain, ngrids
!
!--------------------------------------------------------------------
!
      call gemtime ( 6, '', .false. )
      ierr = model_timeout_alarm (Step_alarm)

!     rpn_comm_mydomain returns mydomain, Grd_ndomains
      call rpn_comm_mydomain (init_ndoms, mydomain)

      write(my_dir,'(a,i4.4)') 'cfg_',mydomain

      ierr = clib_getenv ('TASK_BASEDIR',Path_basedir_S)
      ierr = clib_getenv ('TASK_WORK'   ,Path_work_S   )
      ierr = clib_getenv ('TASK_INPUT'  ,Path_input_S  )
      ierr = clib_getenv ('TASK_OUTPUT' ,Path_output_S )

      Path_input_S  = trim(Path_input_S ) // '/' // trim(my_dir)
      Path_work_S   = trim(Path_work_S  ) // '/' // trim(my_dir)
      Path_output_S = trim(Path_output_S) // '/' // trim(my_dir)

      ierr = clib_chdir (trim(Path_work_S))
      Dcst_rayt_8      = rayt_8
      Dcst_inv_rayt_8  = 1.d0 / rayt_8
      Dcst_omega_8     = omega_8

      ngrids=1
      if (Grd_yinyang_L) ngrids=2

!Special notes for RPN_COMM_init_multi_level
!Must set Ptopo_npex,Ptopo_npey to 0 before calling it!
!It will return myproc, numproc, npex, npey for each PE!
      Ptopo_npex=0
      Ptopo_npey=0
      Ptopo_couleur= RPN_COMM_init_multi_level ( &
                       pe_zero_topo, Ptopo_myproc,Ptopo_numproc, &
                       Ptopo_npex,Ptopo_npey,Grd_ndomains,ngrids )
! print *,'I am proc',Ptopo_myproc, 'npex=',Ptopo_npex, 'numproc=',Ptopo_numproc

      call gemtime_init ( Ptopo_myproc, 'MOD' )
      call gemtime_start ( 1, 'GEMDM', 0)

      if (Grd_yinyang_L) then
         Ptopo_intracomm = RPN_COMM_comm ('GRID')
         Ptopo_intercomm = RPN_COMM_comm ('GRIDPEERS')
         call RPN_COMM_size ('MULTIGRID',Ptopo_world_numproc,ierr)
         call RPN_COMM_rank ('MULTIGRID',Ptopo_world_myproc ,ierr)

         if (Ptopo_couleur == 0) Grd_yinyang_S = 'YIN'
         if (Ptopo_couleur == 1) Grd_yinyang_S = 'YAN'
         ierr= clib_chdir(trim(Grd_yinyang_S))
         Ptopo_ncolors = 2
      else
         Ptopo_couleur = 0
         Ptopo_ncolors = 1
      end if

!     call RPN_COMM_size ( RPN_COMM_GRIDPEERS, Ptopo_nodes, ierr )

      call msg_set_can_write (Ptopo_myproc == 0)

!     PE 0 has these values only from CPU namelist
!     So must broadcast them in PE_ALL_TOPO

      Ptopo_nthreads_phy = Cpus_nthreads_phy
      Ptopo_nthreads_dyn = Cpus_nthreads_dyn

      call pe_all_topo()

      ierr = 0

      G_periodx = .false.
      G_periody = .false.

      l_west  = (0 == Ptopo_mycol)
      l_east  = (Ptopo_npex-1 == Ptopo_mycol)
      l_south = (0 == Ptopo_myrow)
      l_north = (Ptopo_npey-1 == Ptopo_myrow)

      north = 0 ; south = 0 ; east = 0 ; west = 0
      if (l_north) north = 1
      if (l_south) south = 1
      if (l_east ) east  = 1
      if (l_west ) west  = 1
!
!--------------------------------------------------------------------
!
      return
      end

