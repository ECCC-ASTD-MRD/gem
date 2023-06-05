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

!**s/r pe_all_topo - Complete MPI initialization steps

      subroutine pe_all_topo()
      use iso_c_binding
      use clib_itf_mod
      use ctrl
      use HORgrid_options
      use lun
      use omp_timing
      use path
      use ptopo
      use rstr
      use step_options
      use, intrinsic :: iso_fortran_env
      implicit none

      include "rpn_comm.inc"

      integer, external :: fnom, wkoffit!, get_file_unit

      character(len=3)   :: mycol_S, myrow_S
      character(len=1024):: fn!, scratch_dir
      integer :: nc, err, unf
      integer, dimension(0:Ptopo_ncolors-1,-1:Ptopo_npex,-1:Ptopo_npey) :: colrow
!
!-------------------------------------------------------------------
!
      call gemtime ( Lun_out, 'STARTING GEMDM', .false. )
      call gtmg_start ( 2, 'INIT_GEM', 1)

      colrow = 0
      colrow(Ptopo_couleur,Ptopo_mycol,Ptopo_myrow) = Ptopo_myproc
      allocate (Ptopo_colrow(0:Ptopo_ncolors-1,-1:Ptopo_npex,-1:Ptopo_npey))
      nc= Ptopo_ncolors*(Ptopo_npex+2)*(Ptopo_npey+2)
      call rpn_comm_ALLREDUCE (colrow,Ptopo_colrow,nc,&
                               "MPI_INTEGER", "MPI_BOR",'MULTIGRID',err)
      if (Grd_yinyang_L) then
         Ptopo_colrow(:,-1        ,:) = -1
         Ptopo_colrow(:,Ptopo_npex,:) = -1
         Ptopo_colrow(:,:,        -1) = -1
         Ptopo_colrow(:,:,Ptopo_npey) = -1
      else
         Ptopo_colrow(0,-1        ,:) = Ptopo_colrow(0,0           ,:)
         Ptopo_colrow(0,Ptopo_npex,:) = Ptopo_colrow(0,Ptopo_npex-1,:)
         Ptopo_colrow(0,:,        -1) = Ptopo_colrow(0,:,0           )
         Ptopo_colrow(0,:,Ptopo_npey) = Ptopo_colrow(0,:,Ptopo_npey-1)
      endif

      Path_nml_S      = trim(Path_work_S)//'/model_settings.nml'
      Path_outcfg_S   = trim(Path_work_S)//'/output_settings'
      Path_phyincfg_S = trim(Path_input_S)//'/physics_input_table'

! Changing directory to local Ptopo_mycol_Ptopo_myrow

      write(mycol_S,10) Ptopo_mycol
      write(myrow_S,10) Ptopo_myrow
10    format(i3.3)

      if (Ptopo_myproc == 0) call mkdir_gem ('./',Ptopo_npex,Ptopo_npey)
      call rpn_comm_barrier("grid", err)

      err= clib_chdir(mycol_S//'-'//myrow_S)
!
! Determine restart mode and read restart files
!
      Lun_rstrt = 0 ; Rstri_rstn_L= .false. ; Step_kount= 0
      err = wkoffit('gem_restart')
      if (err >= -1) then
         err= fnom( Lun_rstrt,'gem_restart','SEQ+UNF+OLD',0 )
         if (err >= 0) then
           Rstri_rstn_L = .true.
           if (lun_out > 0) write (lun_out,1001)
 1001 format (/' RESTART DETECTED'/)
           call rdrstrt
         end if
      end if
!
! Determine theoretical mode with presence of file ${TASK_WORK}/theoc
!
      unf=0
      Ctrl_theoc_L = .false.
      fn=trim(Path_work_S)//'/theoc'
      if (wkoffit(fn) > -3) then
         if (Ptopo_myproc == 0) write (Lun_out,*) &
                                'Assume Theoretical case'
         Ctrl_theoc_L = .true.
      else
         call fclos (unf)
      end if
!
!-------------------------------------------------------------------
!
      return
      end

      subroutine mkdir_gem (F_path_S,F_npex,F_npey)
      use lun
      use path
      use clib_itf_mod
      use omp_timing
      implicit none

      character(len=*), intent(in) :: F_path_S
      integer      , intent(in) :: F_npex,F_npey

      character(len=2048) ici
      integer :: i,j,ret,err
      integer :: mk_gem_dir

      err= clib_getcwd(ici)
      err= clib_chdir(trim(F_path_S))

      do j=0,F_npey-1
         do i=0,F_npex-1
            ret = mk_gem_dir(i,j)
         end do
      end do

      err= clib_chdir(ici)

      return
      end

      integer function mk_gem_dir(x,y)
      use ISO_C_BINDING
      use lun
      use path
      use omp_timing
      implicit none
      integer, intent(in) :: x, y
      character (len=7) :: filename

      character (len=1), dimension(8), target :: temp

      interface
         integer(C_INT) function mkdir(path,mode) bind(C,name='mkdir')
            use ISO_C_BINDING
            type(C_PTR), value :: path
            integer(C_INT), value :: mode
         end function mkdir
      end interface

      write(filename,100)x,'-',y
100   format(I3.3,A1,I3.3)

      temp = transfer( trim(filename)//achar(0) , temp )
      mk_gem_dir = mkdir(C_LOC(temp),511)

      return
      end function mk_gem_dir
