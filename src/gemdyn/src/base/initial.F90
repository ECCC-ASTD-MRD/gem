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

!**s/r initial - Performs initialisation

      subroutine initial (F_rstrt_L)
      use step_options
      use gmm_geof
      use gmm_pw
      use gem_options
      use lam_options
      use HORgrid_options
      use init_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use rstr
      implicit none
#include <arch_specific.hf>

      logical, intent(inout) :: F_rstrt_L

!arguments
!  Name                      Description
!----------------------------------------------------------
! F_rstrt_L         TRUE if a restart is required
!----------------------------------------------------------

      integer n
      real    prn, promegc, prsum, prwin1, prwin2
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,2020) Init_dfnp

      allocate (Init_dfco(0:Init_halfspan))

      promegc = (2.0 * pi_8) / Init_dfpl_8
      prwin1  = pi_8 / real(Init_halfspan + 1)

      Init_dfco(0) = promegc * Cstv_dt_8 / pi_8
      prsum        = Init_dfco(0)

      do n=1,Init_halfspan
         prwin2 = 1.0
         prn    = real(n)

         if ( Init_dfwin_L ) then
            prwin2 = prn * prwin1
            prwin2 = sin(prwin2) / prwin2
         end if

         Init_dfco(n) = prwin2 *dsin(prn * promegc * Cstv_dt_8) /  &
                        (prn * pi_8)
         prsum     = prsum + 2.0 * Init_dfco(n)
      end do

      if (Lun_out > 0) write(Lun_out,2030)
      if (Lun_out > 0) then
         write(Lun_out,*) (Init_dfco(n),n=0,Init_halfspan), prsum
      end if

      do n=0,Init_halfspan
         Init_dfco(n) = Init_dfco(n) / prsum
      end do

      prsum = Init_dfco(0)
      do n=1,Init_halfspan
         prsum = prsum + 2.0 * Init_dfco(n)
      end do

      if (Lun_out > 0) write(Lun_out,2040)
      if (Lun_out > 0) then
         write(Lun_out,*) (Init_dfco(n),n=0,Init_halfspan), prsum
      end if

      if ( .not. Rstri_rstn_L ) call digflt()

      if (Lun_out > 0) then
         write(Lun_out,1000) Lctl_step,Lctl_step+Init_dfnp-1-Step_kount
      end if

      call gem_run (F_rstrt_L)

      if ((Step_kount == Init_dfnp-1).and.(.not.F_rstrt_L)) then
         Init_mode_L = .false.
         call ta2t1tx()
         call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,pw_log_pt, &
                         pw_pm_plus_8,pw_p0_plus_8, &
                         l_minx,l_maxx,l_miny,l_maxy,l_nk,1 )
         call pw_update_GW ()
         call pw_update_UV ()
         call pw_update_T  ()
         Lctl_step = Lctl_step  - Init_halfspan
         Step_kount= Step_kount - Init_halfspan
         if ( .not. Grd_yinyang_L .and. .not. Lam_ctebcs_L) then
            call nest_intt
         endif
         if (Vtopo_start >= 0 .and. Lctl_step-Vtopo_start+1 <= Vtopo_ndt) Vtopo_L = .true.
         call oro_adj ()
!!$         if (Vtopo_L) then
!!$            call var_topo (fis0, F_topo_LS, real(Lctl_step), l_minx,l_maxx,l_miny,l_maxy)
!!$            call rpn_comm_xch_halo (fis0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,1,&
!!$                    G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!!$         end if
         call adz_inittraj ()
         if (Lun_out > 0) write(Lun_out,1050) Lctl_step
      end if

 1000 format(/,' =====> DIGITAL FILTER INITIALIZATION SCHEME: ', &
               'TIMESTEP ',i3,' to TIMESTEP ',i3,/' ',73('='))
 1050 format(/,' =====> DIGITAL FILTER INITIALIZATION SCHEME: ', &
               'COMPLETED',/8x,'LAST TIMESTEP COMPLETED RESET TO:', &
               I5,/' ',55('='))
 2020 format(/,'PREPARATION OF DIGITAL FILTER PARAMETERS', &
      /,'AND COEFFICIENTS (S/R INITIAL)          ', &
      /,'========================================',/,'Init_dfnp  = ',i4)
 2030 format(/,'Digital filter coefficients and sum', &
      /,'before normalization               ')
 2040 format(/,'Digital filter coefficients and sum', &
      /,'after normalization                ')
!
!     ---------------------------------------------------------------
!
      return
      end
