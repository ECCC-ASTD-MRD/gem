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
      use iso_c_binding
      use clib_itf_mod
      use dcst
      use glb_ld
      use domains
      use HORgrid_options
      use path
      use ptopo
      use step_options
      use tdpack
      use gem_timing
      use version
      implicit none
#include <arch_specific.hf>

      include 'rpn_comm.inc'
      include 'gemdyn_version.inc'
#include <rmnlib_basics.hf>

      integer, external :: model_timeout_alarm, OMP_get_max_threads

      character(len=256) :: my_dir
      logical :: is_official_L
      integer :: ierr, mydomain
      integer, parameter :: nargs=10, npos=0
      character(len=16) listec(nargs)
      character(len=2048) def(nargs), val(nargs)
      character(len=50) :: DSTP,name_S,arch_S,compil_S,user_S
!
!--------------------------------------------------------------------
!
      call gemtime ( 6, '', .false. )
      ierr = model_timeout_alarm (Step_alarm)

      ! List of non-positional optional arguments to the binary maingem
      ! with coresponding first and second default values (val & def)
      ! In MPI the arguments are provided using export CCARD_ARGS=
      listec = [ CHARACTER(LEN=16) :: &
                    'npex.', 'npey.', 'smtdyn.', 'smtphy.','ngrids.',&
                    'dom_start.', 'dom_end.', 'dom_last.',&
                    'input.', 'output.' ]
      val    = [ '1', '1', '0', '0', '1', '0', '0', '0', ' ', ' ' ]
      def    = [ '1', '1', '0', '0', '1', '0', '0', '0', ' ', ' ' ]
      ! Obtain values of calling arguments (listec)
      call ccard (listec,def,val,nargs,npos)
      read (val(1),*) Ptopo_npex
      read (val(2),*) Ptopo_npey
      read (val(3),*) Ptopo_nthreads_dyn
      read (val(4),*) Ptopo_nthreads_phy
      read (val(5),*) Domains_ngrids
      read (val(6),*) Domains_deb
      read (val(7),*) Domains_fin
      read (val(8),*) Domains_last
      Path_input_S = trim(val( 9))
      Path_output_S= trim(val(10))

      ierr = clib_getenv ('PWD', Path_work_S)

      Domains_num= Domains_fin - Domains_deb + 1
      Ptopo_last_domain_L = (Domains_fin == Domains_last)

      Grd_yinyang_L = .false. ; Grd_yinyang_S = ''
      if (Domains_ngrids == 2) Grd_yinyang_L = .true.

      ! Obtain mydomain
      mydomain = RPN_COMM_get_my_domain(Domains_num, Domains_deb)

      write(my_dir,'(a,i4.4)') 'cfg_',mydomain

      Path_input_S  = trim(Path_input_S ) // '/' // trim(my_dir)
      Path_work_S   = trim(Path_work_S  ) // '/' // trim(my_dir)
      Path_output_S = trim(Path_output_S) // '/' // trim(my_dir)

      ! cd to work directory
      ierr = clib_chdir (trim(Path_work_S))

      Dcst_rayt_8      = rayt_8
      Dcst_inv_rayt_8  = 1.d0 / rayt_8
      Dcst_omega_8     = omega_8

      ! Start MPI on a (Ptopo_npex x Ptopo_npey) processors topology and
      ! obtain Ptopo_myproc, Ptopo_numproc, Ptopo_mycol and Ptopo_myrow
      Ptopo_couleur= RPN_COMM_gridinit ( Ptopo_myproc,Ptopo_numproc, &
                     Ptopo_npex,Ptopo_npey,Domains_num,Domains_ngrids )
      ierr = RPN_COMM_mype (Ptopo_myproc, Ptopo_mycol, Ptopo_myrow)

      lun_out     = -1
      Lun_debug_L = .false.

      if (Ptopo_myproc == 0) then
         lun_out = output_unit
         ierr= clib_mkdir (Path_output_S)

         call  open_status_file3 (trim(Path_output_S)//'/status_MOD.dot')
         call write_status_file3 ('_status=ABORT' )

         call atm_model_getversion2 (name_S,Version_number_S,DSTP,&
                              arch_S,compil_S,user_S,is_official_L)
         if (is_official_L) then
            Version_title_S=trim(name_S)//' --- Release of: '//trim(DSTP)
         else
            Version_title_S=trim(name_S)//' --- '//trim(user_S)//&
                              ' Build: '//trim(DSTP)
         end if

         ierr = exdb(trim(Version_title_S),trim(Version_number_S),'NON')
      endif

      call gemtime_init ( Ptopo_myproc, 'MOD' )
      call gemtime_start ( 1, 'GEMDM', 0)

      ! Some MPI cummunicators + init colors
      Ptopo_intracomm = RPN_COMM_comm ('GRID')
      if (Grd_yinyang_L) then
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
         Ptopo_world_numproc = Ptopo_numproc
         Ptopo_world_myproc  = Ptopo_myproc
      end if

      ! Initialize OpenMP
      Ptopo_npeOpenMP = OMP_get_max_threads()
      if (Ptopo_nthreads_dyn < 1) Ptopo_nthreads_dyn=Ptopo_npeOpenMP
      if (Ptopo_nthreads_phy < 1) Ptopo_nthreads_phy=Ptopo_npeOpenMP
      call set_num_threads ( Ptopo_nthreads_dyn, 0 )

      if (Lun_out > 0) then
         write (Lun_out, 8255) Ptopo_npex, Ptopo_npey, &
             Ptopo_npeOpenMP,Ptopo_nthreads_dyn, Ptopo_nthreads_phy
         write (Lun_out, 8256) trim(Path_work_S)
      end if

      call msg_set_can_write (Ptopo_myproc == 0)

      call pe_all_topo()

      ! Initialize local sub domain boundaries flags
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

 8255 format (/," MPI CONFIG (npex x npey): ",i4,' x ',i3,/, &
          " OMP CONFIG (npeOpenMP x nthreads_dyn x nthreads_phy): ",&
            i4,' x ',i3,' x ',i3)
 8256       format (/," WORKING DIRECTORY:"/a/)
!
!--------------------------------------------------------------------
!
      return
      end subroutine init_component
