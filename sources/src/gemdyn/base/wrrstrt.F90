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

!s/r wrrstrt - Write the restart file
!
      subroutine wrrstrt ()
      use iso_c_binding
      use phy_itf, only: phy_restart
      use step_options
      use init_options
      use lun
      use psadjust
      use gmm_itf_mod
      use wb_itf_mod
      use gem_timing
      use adz_BC_deficit
      implicit none
#include <arch_specific.hf>

      include "rpn_comm.inc"

      integer, external :: fnom,fclos
      integer ier,gmmstat,me,howmany,newcomm,i
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,2000) Lctl_step

      call split_on_hostid (RPN_COMM_comm('GRID'),me,howmany,newcomm)

      call gemtime_start ( 33, 'RESTART', 34 )
      do i=0,howmany-1
         if (i == me) then

            Lun_rstrt = 0
            ier = fnom (Lun_rstrt,'gem_restart','SEQ+UNF',0)

            write(Lun_rstrt) Lctl_step,Step_kount,Init_mode_L
            write(Lun_rstrt) PSADJ_g_avg_ps_initial_8,PSADJ_scale_8,PSADJ_fact_8
            write(Lun_rstrt) KEEP_mass_deficit_8(1:MAXTR3D_),tracer_name(1:MAXTR3D_)

            ier = fclos(Lun_rstrt)

            !        Write Gmm-files

            gmmstat = gmm_checkpoint_all(GMM_WRIT_CKPT)
            ier = wb_checkpoint()

            ier = phy_restart ('W', .false.)

         end if

         call mpi_barrier (newcomm,ier)

      end do
      call gemtime_stop (33)

 2000 format(/,'WRITING A RESTART FILE AT TIMESTEP #',I8, &
             /,'=========================================')
!
!     ---------------------------------------------------------------
!
      return
      end
