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

!**s/r fislh_dynstep - Control of the dynamical timestep of the model

      subroutine fislh_dynstep ()
      use dyn_fisl_options
      use gem_options
      use HORgrid_options
      use step_options
      use theo_options
      use omp_timing
      implicit none

      integer icn, keep_itcn
!
!     ---------------------------------------------------------------
!
      keep_itcn = Schm_itcn

      call psadj_init ( Step_kount )

!$omp parallel
      call gtmg_start (10, 'DYNSTEP', 10)
      do icn = 1,Schm_itcn-1

         call fislh_tstpdyn (icn) ! Solver NOT done yet

         call hzd_momentum() ! NOT done yet

      end do

      call fislh_tstpdyn (Schm_itcn)

      if (Ctrl_theoc_L .and. .not.Grd_yinyang_L) call theo_bndry ()

      call adz_tracers_hlt (.true.) ! Mass fixing NOT done yet

      call psadj_hlt ( Step_kount )

      call adz_tracers_hlt (.false.)
      call gtmg_stop (10)

      call t02t1()
!$omp end parallel

      call HOR_bndry ()

      call canonical_cases ("VRD")

      call hzd_main ()

      if (Grd_yinyang_L) call yyg_blend()

      call pw_update_GW()
      call pw_update_UV()
      call pw_update_T()

      if ( Lctl_step-Vtopo_start == Vtopo_ndt) Vtopo_L = .false.

      Schm_itcn = keep_itcn
!
!     ---------------------------------------------------------------
!
      return
      end subroutine fislh_dynstep
