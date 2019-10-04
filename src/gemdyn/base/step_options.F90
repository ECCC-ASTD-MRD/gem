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
module step_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Starting date for model run  (yyyymmdd.hhmmss)
   character(len=16) :: Step_runstrt_S = 'NIL'
   namelist /step/ Step_runstrt_S

   !# Starting date for model run slice (yyyymmdd.hhmmss)
   character(len=16) :: Fcst_start_S =  ' '
   namelist /step/ Fcst_start_S

   !# End date for model run slice (yyyymmdd.hhmmss)
   character(len=16) :: Fcst_end_S =  ' '
   namelist /step/ Fcst_end_S

   !# Read nesting data every Fcst_nesdt_S
   character(len=16) :: Fcst_nesdt_S = ' '
   namelist /step/ Fcst_nesdt_S

   !# Output global stat (glbstat) every Fcst_gstat_S
   character(len=16) :: Fcst_gstat_S = ' '
   namelist /step/ Fcst_gstat_S

   !# Save a restart file + stop every Fcst_rstrt_S
   character(len=16) :: Fcst_rstrt_S = 'NIL'
   namelist /step/ Fcst_rstrt_S

   !# Save a restart file + continue every Fcst_bkup_S
   character(len=16) :: Fcst_bkup_S = 'NIL'
   namelist /step/ Fcst_bkup_S

   !#
   character(len=16) :: Fcst_spinphy_S = ' '
   namelist /step/ Fcst_spinphy_S

   !# Save a restart file + continue at that time
   character(len=16) :: Fcst_bkup_additional_S = 'NIL'
   namelist /step/ Fcst_bkup_additional_S

   !# Setting for Fortran alarm time
   integer :: Step_alarm = 600
   namelist /step/ Step_alarm

   !# Account for leap years
   logical :: Step_leapyears_L = .true.
   namelist /step/ Step_leapyears_L

   !# Length of model timestep (sec)
   real(kind=REAL64)  :: Step_dt = -1.
   namelist /step/ Step_dt

   ! Internal variables NOT in step namelist

   integer Step_total, Step_gstat, Step_delay, Step_spinphy, &
           Step_kount, Step_CMCdate0, Step_initial         , &
           Step_bkup_additional, Lctl_step

   real(kind=REAL64)  Step_nesdt, Step_maxwall

contains

!**s/r step_nml - Read namelist step

      integer function step_nml (F_unf)
      use timestr_mod
      use lun
      use rstr
      use clib_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer err
      real(kind=REAL64) nesdt,nsteps
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         step_nml= 0
         if ( Lun_out >= 0) then
            write (Lun_out,nml=step)
            if ( F_unf == -1 ) write (Lun_out,8000) Step_alarm
         end if
         return
      end if

      step_nml= -1 ; nml_must= .true. ; nml_S= 'step'

      print *, '================================================ HAMSTER'
      rewind(F_unf)
      read (F_unf, nml=step, end= 1001, err=1003)
      step_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         step_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 print *, 'Step_runstrt_S=', Step_runstrt_S
      if (step_nml < 0 ) return
      if ((Lun_out>=0).and.(step_nml==0)) write (Lun_out, 6004) trim(nml_S)
      step_nml= -1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
 6009 format (//,' NAMELIST ',A,' IS MANDATORY'//)
! boiler plate - end

      if (Step_dt < 0.) then
         if (Lun_out > 0) write(Lun_out,*)  &
                    ' Step_dt must be specified in namelist &step'
         goto 9999
      end if

      err= 0

      if ( Fcst_start_S  == '' ) Fcst_start_S = '0H'
      if ( Fcst_end_S    == '' ) Fcst_end_S   = Fcst_start_S

      err= min( timestr2step (Step_initial, Fcst_start_S, Step_dt), err)
      err= min( timestr2step (Step_total  , Fcst_end_S  , Step_dt), err)
! transforming the Step_total into actual number of timesteps
      Step_total= Step_total - Step_initial

! Fcst_nesdt_S is transformed into a number of secondes (into Step_nesdt)
      nesdt= 1.d0
      err= min( timestr2step (nsteps, Fcst_nesdt_S, nesdt), err)
      Step_nesdt= dble(nsteps)

      if ( Fcst_rstrt_S /= 'NIL' ) then
         err= timestr_check ( Fcst_rstrt_S )
      end if

      Step_bkup_additional= Step_total+1
      err = clib_toupper ( Fcst_bkup_additional_S )
      if ( Fcst_bkup_additional_S /= 'NIL' ) then
         if (Fcst_bkup_additional_S == 'END' ) then
            Step_bkup_additional= Step_total
         else
            err= min( timestr2step (Step_bkup_additional, &
                      Fcst_bkup_additional_S, Step_dt), err)
         end if
      end if

      err = clib_toupper ( Fcst_bkup_S )
      if ( Fcst_bkup_S == 'END' ) Fcst_bkup_S= Fcst_end_S
      if ( Fcst_bkup_S /= 'NIL' ) then
         err = timestr_check (Fcst_bkup_S)
      end if

      if ( Fcst_gstat_S  == '' ) then
         Step_gstat= Step_total-Step_initial+1
      else
         err= min( timestr2step (Step_gstat, Fcst_gstat_S, Step_dt), err)
      end if
      if ( Fcst_spinphy_S  == '' ) then
         Step_spinphy= Step_total-Step_initial+1
      else
         err= min( timestr2step (Step_spinphy, Fcst_spinphy_S, Step_dt), err)
      end if

      if (err < 0) goto 9999

      Step_delay= Step_initial

      if (.not.Rstri_rstn_L) Lctl_step= Step_initial

      step_nml = 1

 8000 format (/,' MODEL ALARM SET TO: ',i8,' secondes'/)
!
!-------------------------------------------------------------------
!
 9999 return
      end function step_nml

   function step_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function step_options_init

end module step_options
