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

!**s/r gem_run - Performs the integration of the model

      subroutine gem_run (F_rstrt_L)
      use cstv
      use ctrl
      use dcmip_options
      use dynkernel_options
      use glb_ld
      use gmm_itf_mod
      use gmm_vt1
      use HORgrid_options
      use init_options
      use rstr
      use step_options
      use gem_timing
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical F_rstrt_L

!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_rstrt_L     O         Is a restart required
!----------------------------------------------------------------

      logical, external :: gem_muststop
      integer, external :: model_timeout_alarm
      character(len=16) :: datev
      integer :: stepf,last_step, seconds_since, istat
      real(kind=REAL64) :: dayfrac

      real(kind=REAL64), parameter :: sec_in_day = 86400.0d0
!
!     ---------------------------------------------------------------
!
      dayfrac = dble(Step_kount) * Cstv_dt_8 / sec_in_day
      call incdatsd (datev,Step_runstrt_S,dayfrac)

      if (Lun_out > 0) write (Lun_out,900) datev

      call blocstat (.true.)

      call gemtime ( Lun_out, 'STARTING TIME LOOP', .false. )

      stepf= Step_total
      last_step = Step_total + Step_initial
      if (Init_mode_L) then
         stepf     = Init_dfnp-1
         last_step = stepf + Step_initial
      end if

      F_rstrt_L = .false.
      if ( .not. Rstri_rstn_L ) then
         call out_outdir()
         call out_dyn (.true., .true.)
         if (gem_muststop (stepf)) then
            seconds_since = model_timeout_alarm(Step_alarm)
            if (Lun_out > 0) write(Lun_out,4000) Lctl_step
            return
         end if
      end if

      call canonical_cases ("ERR")

      do while (Step_kount < stepf)

         seconds_since = model_timeout_alarm(Step_alarm)

         Lctl_step= Lctl_step + 1  ;  Step_kount= Step_kount + 1

         if (Lun_out > 0) then

            select case (trim(Dynamics_Kernel_S))
               case ('DYNAMICS_EXPO_H')
                  write (Lun_out,1001) Lctl_step, last_step
               case ('DYNAMICS_FISL_H')
                  write (Lun_out,1002) Lctl_step, last_step
               case ('DYNAMICS_FISL_P')
                  write (Lun_out,1003) Lctl_step, last_step
               case default
                  stop 'Something has gone horribly wrong. Back away slowly'
            end select

         end if

         call out_outdir()

         call pw_shuffle()

         call dynstep()

         call out_dyn (.false., .true.) ! casc output

         if ( Ctrl_phyms_L ) call itf_phy_step (Step_kount, Lctl_step)

         call canonical_cases ("PHY")

         call canonical_cases ("ERR")

         call iau_apply (Step_kount)

         if (Grd_yinyang_L .or. Dynamics_Kernel_S(1:13) == 'DYNAMICS_FISL') then
            call yyg_xchng_all()
         end if

         if ( Ctrl_phyms_L.or.Dcmip_physics_L ) then
            istat = gmm_get (gmmk_tt1_s, tt1)
            call tt2virt (tt1, .true., &
            l_minx,l_maxx,l_miny,l_maxy, G_nk)
            call itf_phy_UVupdate()
            call pw_update_GPW()
         end if

         if ( Init_mode_L ) call digflt ! digital filter

         call out_dyn (.true., .false.) ! regular output

         call blocstat (.false.)

         if (Lun_out > 0) write(Lun_out,3000) Lctl_step

         call save_restart()

         F_rstrt_L= gem_muststop (stepf)

         if (F_rstrt_L) exit

      end do

      seconds_since = model_timeout_alarm(Step_alarm)

      if (Lun_out > 0) write(Lun_out,4000) Lctl_step

 900  format (/'STARTING THE INTEGRATION WITH THE FOLLOWING DATA: VALID ',a)
 1001 format(/,'EXPO: PERFORMING TIMESTEP #',I9,' OUT OF ',I9, &
             /,'=========================================================')
 1002 format(/,'FISL-H: PERFORMING TIMESTEP #',I9,' OUT OF ',I9, &
             /,'=========================================================')
 1003 format(/,'FISL-P: PERFORMING TIMESTEP #',I9,' OUT OF ',I9, &
             /,'=========================================================')
 3000 format(/,'THE TIME STEP ',I8,' IS COMPLETED')
 4000 format(/,'GEM_RUN: END OF THE TIME LOOP AT TIMESTEP',I8, &
             /,'===================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
