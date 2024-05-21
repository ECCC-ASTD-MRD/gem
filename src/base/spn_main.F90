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
!---------------------------------- LICENCE END --------------------------------
!
!*s/r spn_main - spectral nudging driver

      subroutine spn_main ()
      use step_options
      use spn_options
      use glb_ld
      use lun
      use gmm_vt1
      use mem_nest
      use tdpack
      use HORgrid_options
      use rmn_gmm
      use var_gmm
      use tr3d
      use mem_tracers
      use omp_timing
      implicit none
!
!----------------------------------------------------------------------
!

      real, pointer, dimension(:,:,:) :: nudge_t, nudge_u, nudge_v, nudge_hu
      integer:: istat, ib, ie

      if (.not. Spn_ON_L) return

      call gtmg_start ( 59, 'SPN_main', 1 )

      if (( Lctl_step > 2 ) .and. ( mod(Lctl_step,Spn_interval) == 0 ))  then
      
      
         if (Lun_out > 0) write(Lun_out,1001) Lctl_step
      
         if (Spn_weight_L) then
            Spn_weight = sqrt((cos(pi_8*(float(Lctl_step)/float(Spn_ws))))**2)**Spn_wt_pwr
         end if

         if (.not. Grd_yinyang_L) then
            call spn_apply (tt1, nest_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)
            call spn_apply (ut1, nest_u, l_minx,l_maxx,l_miny,l_maxy, l_nk)
            call spn_apply (vt1, nest_v, l_minx,l_maxx,l_miny,l_maxy, l_nk)
            if (Spn_nudge_HU_L) then
               ib=(Tr3d_hu-1)*G_nk+1
               ie=ib+G_nk-1
               istat=gmm_get('N_hu',nudge_hu)
               call spn_apply (tracers_P(Tr3d_hu)%pntr, nest_tr(:,:,ib:ie), l_minx,l_maxx,l_miny,l_maxy, l_nk)
            end if
         else
            istat=gmm_get('N_tt',nudge_t)
            istat=gmm_get('N_uu',nudge_u)
            istat=gmm_get('N_vv',nudge_v)
            call spn_apply (tt1, nudge_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)
            call spn_apply (ut1, nudge_u, l_minx,l_maxx,l_miny,l_maxy, l_nk)
            call spn_apply (vt1, nudge_v, l_minx,l_maxx,l_miny,l_maxy, l_nk)

            if (Spn_nudge_HU_L) then
               istat=gmm_get('N_hu',nudge_hu)
               call spn_apply (tracers_P(Tr3d_hu)%pntr, nudge_hu, l_minx,l_maxx,l_miny,l_maxy, l_nk)
            end if
         end if
      
      end if

      call gtmg_stop ( 59 )
 1001 format(/' SPN_MAIN: Applying spectral nudging at STEP NO ',I10/)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_main
