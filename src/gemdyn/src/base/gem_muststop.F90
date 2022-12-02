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
!**s/r gem_muststop

      logical function gem_muststop (F_finalstep)
      use iso_c_binding
      use clib_itf_mod
      use cstv
      use gem_options
      use init_options
      use lun
      use out_mod
      use path
      use ptopo
      use step_options
      use omp_timing
      use timestr_mod, only: timestr_isstep, TIMESTR_MATCH
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_finalstep

      include "rpn_comm.inc"

      character(len=2048) :: filen, filen_link, append
      logical :: flag, pe0_master_L, output_L, finalstep_L, end_of_run_L
      integer :: err
      real(kind=REAL64) :: timeleft, hugetype
      integer, parameter :: unf = 474
!
!     ---------------------------------------------------------------
!
      call gtmg_start ( 70, 'MUSTOP', 1 )

      pe0_master_L = (Ptopo_myproc == 0) .and. (Ptopo_couleur == 0)

      gem_muststop = .false.
      ! Get timeleft to determine if we can continue
      if (Lctl_cktimeleft_L) then
         if (pe0_master_L) then
            filen=trim(Path_work_S)//'/../time_left'

            open (unf,file=trim(filen),access='SEQUENTIAL',&
                  status='OLD',iostat=err,form='FORMATTED')
            if (err == 0) then
               read (unf,*,end=33,err=33) timeleft
 33            close(unf)
            else
               timeleft= huge(hugetype)
            end if
            flag = (timeleft < Step_maxwall)
         end if

         call RPN_COMM_bcast (flag, 1, "MPI_LOGICAL",0,"MULTIGRID",err)
         gem_muststop = flag
      end if

      if ( Fcst_rstrt_S /= 'NIL' ) then
         if ( (Step_kount > 0) .and. (.not. Init_mode_L .or. (Step_kount < Init_halfspan)) )  then
            gem_muststop = gem_muststop .or. &
               timestr_isstep ( Fcst_rstrt_S, Step_CMCdate0, real(Cstv_dt_8), Step_kount ) == TIMESTR_MATCH
         end if
      end if

      filen      = trim(Path_output_S)//'/output_ready_MASTER'
      filen_link = trim(Path_output_S)//'/output_ready'
      finalstep_L = Step_kount == F_finalstep
      end_of_run_L= (finalstep_L .and. (.not.Init_mode_L))

      if ( Out_post_L .or. end_of_run_L .or. gem_muststop) then
         if (pe0_master_L) then
            append=''
            if (finalstep_L .and. Ptopo_last_domain_L) append='^last'
            open  ( unf,file=filen,access='SEQUENTIAL',&
                    form='FORMATTED',position='APPEND' )
            write (unf,'(3(a))') 'NORMAL_OUTPUT ','NA ',trim(Out_laststep_S)//trim(append)
            close (unf)
         end if
      end if

      ! Send a signal to gem_monitor_output
      output_L = Out_post_L .or. end_of_run_L
      if ( output_L .and. (.not.finalstep_L) ) then

         call rpn_comm_barrier (RPN_COMM_ALLGRIDS, err)
         if (pe0_master_L) then
            err = clib_symlink ( trim(filen), trim(filen_link) )
            write (output_unit,1001) trim(Out_laststep_S),lctl_step
         end if

      end if

      gem_muststop = gem_muststop .and. .not.end_of_run_L

      if (gem_muststop) call wrrstrt()

      call gemtime ( Lun_out, 'CURRENT TIMESTEP', .false. )
      call gtmg_stop (70)

 1001 format (' OUT_LAUNCHPOST: DIRECTORY output/',a, &
              ' was released for postprocessing at timestep: ',i9)
!
!     ---------------------------------------------------------------
!
      return
      end
