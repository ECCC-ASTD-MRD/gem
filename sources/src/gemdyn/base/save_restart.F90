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
!**s/r save_restart
!
      subroutine save_restart
      use iso_c_binding
      use timestr_mod, only: timestr_isstep, TIMESTR_MATCH
      use step_options
      use HORgrid_options
      use gem_options
      use init_options
      use cstv
      use lun
      use out_mod
      use path
      use clib_itf_mod
      use ptopo
      implicit none

#include <arch_specific.hf>

      character(len=2048) :: dirname_S,cmd
      character(len=16) :: datev
      character(len=3) :: mycol_S,myrow_S
      logical saverestart, save_additional
      integer err
      real*8 dayfrac
      real*8, parameter :: OV_day = 1.d0/86400.0d0
!
!     ---------------------------------------------------------------
!
      if ( Init_mode_L .and. (Step_kount >= Init_halfspan) ) return
      save_additional= (Step_kount == Step_bkup_additional)
      saverestart= .false.
      if ( Fcst_bkup_S /= 'NIL' ) then
         saverestart = ( timestr_isstep (Fcst_bkup_S, Step_CMCdate0, &
                         real(Cstv_dt_8),Step_kount) == TIMESTR_MATCH )
      end if

      if (saverestart .or. save_additional) then
         call wrrstrt ()
         dayfrac = dble(Step_kount) * Cstv_dt_8 * OV_day
         call incdatsd (datev,Step_runstrt_S,dayfrac)
         dirname_S= trim(Path_output_S)//'/'//Out_laststep_S//'/restart_'//trim(datev)
         if (Grd_yinyang_L) then
            dirname_S=trim(dirname_S)//'/'//trim(Grd_yinyang_S)
         end if
         if (Ptopo_myproc == 0) then
            err= clib_mkdir_r ( trim(dirname_S) )
            call mkdir_gem ( trim(dirname_S), Ptopo_npex, Ptopo_npey )
         end if
         call rpn_comm_Barrier("grid", err)
         write(mycol_S,'(i3.3)') Ptopo_mycol
         write(myrow_S,'(i3.3)') Ptopo_myrow
         cmd='mv *_restart '//trim(dirname_S)//'/'//mycol_S//'-'//myrow_S
         call system (cmd)
!yet another solution...
!err = clib_glob(filelist,nfiles,'_restart',maxnfiles)
!do ifile = 1,nfiles
!    istat = clib_basename(filelist(ifile),myname)
!    pathnew = trim(newdir)//'/'//trim(myname)
!    istat = clib_rename(filelist(ifile),pathnew)
!end do
         if (Lun_out > 0) write (Lun_out,1002) Lctl_step,datev
      end if

 1002 format (' SAVING A RESTART AT TIMESTEP: ',i7,' valid: ',a)
!
!     ---------------------------------------------------------------
!
      return
      end
