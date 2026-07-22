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

      subroutine out_outdir()
      use iso_c_binding
      use timestr_mod, only: timestr_prognum,timestr_unitfact
      use step_options
      use HORgrid_options
      use init_options
      use out_options
      use cstv
      use lun
      use svro_mod
      use out_mod
      use out3
      use path
      use clib_itf_mod
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

#include <rmnlib_basics.hf>
      include 'mpif.h'
      include "rpn_comm.inc"

      character(len=1024),save :: dirstep_S=' ', dirbloc_S=' '
      character(len=1024) :: dirname_S
      character(len=10) :: postjob_S
      character(len=7)  :: blocxy_S
      logical :: create_dir_L
      integer :: err,last_step_post,flag_step_post,stepno, &
                 prognum,prognum1,upperlimit
      real :: interval
      real(kind=REAL64) :: fatc_8
!
!----------------------------------------------------------------------
!
      upperlimit = Step_total + Step_initial

      call out_steps ()

      if ( Init_mode_L .and. (Step_kount >= Init_halfspan) ) return

      write (blocxy_S,'(I3.3,"-",I3.3)') Ptopo_mycol, Ptopo_myrow

      if (Out3_postproc_fact <= 0) then
         last_step_post = upperlimit
         flag_step_post = upperlimit
         Out_post_L = .false.
      else
         interval = Out3_close_interval * Out3_postproc_fact
         stepno = max(Step_kount+Step_delay,1)

         err = timestr_prognum(prognum ,Out3_unit_S,interval,Out_dateo,&
                               float(Out_deet),stepno  ,Out_endstepno)
         err = timestr_prognum(prognum1,Out3_unit_S,interval,Out_dateo,&
                               float(Out_deet),stepno+1,Out_endstepno)
         Out_post_L = (prognum1 > prognum .or. stepno == Out_endstepno)

         if (Out3_unit_S(1:3) == 'MON') then
            last_step_post = prognum
         else
            fatc_8 = timestr_unitfact(Out3_unit_S,Cstv_dt_8)
            last_step_post = nint(dble(prognum) * fatc_8)
            last_step_post = min(last_step_post,upperlimit)
         end if
      end if

      if (last_step_post >= 0) then
         write (postjob_S,'(i10.10)') last_step_post
      else
         write (postjob_S,'(i10.9) ') last_step_post
      end if

      Out_laststep_S = 'laststep_'//postjob_S
      dirname_S      = trim(Path_output_S)//'/'//Out_laststep_S
      create_dir_L   = .false.
      if (dirstep_S /= Out_laststep_S) then
         create_dir_L = .true.
         dirstep_S = Out_laststep_S
         if (OUTs_server_L) then
            Out_dirname_S= trim(dirname_S)
         else
            Out_dirname_S= trim(dirname_S)//'/'//blocxy_S
         endif
      endif

      if (create_dir_L) then
! PE0 is responsible for creating shared subdir structure
         if (Ptopo_myproc == 0 .and. Ptopo_couleur == 0) then
            err = clib_mkdir(trim(dirname_S))
            if (Lun_out>0) write(Lun_out,1001) trim(Out_laststep_S),Step_kount
         end if
         call rpn_comm_barrier (RPN_COMM_ALLGRIDS, err)
         if ( .not. OUTs_server_L) then
            err = CLIB_OK
            if (Out3_iome >= 0 .and. dirbloc_S /= Out_dirname_S &
                         .and. Ptopo_couleur == 0) then
               dirbloc_S = Out_dirname_S
               err= clib_mkdir (trim(Out_dirname_S))
               err= clib_isdir (trim(Out_dirname_S))
            endif
            call gem_error (err,'out_outdir','unable to create output directory structure')
         endif
      end if

      call OUTs_start ()

 1001 format (' OUT_OUTDIR: DIRECTORY output/',a,' was created at timestep: ',i9)
!
!----------------------------------------------------------------------
!
      return
      end
