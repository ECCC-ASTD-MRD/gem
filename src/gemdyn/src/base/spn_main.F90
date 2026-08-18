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
      use cstv
      implicit none
!
!----------------------------------------------------------------------
!

      real, pointer, dimension(:,:,:) :: nudge_t, nudge_u, nudge_v, nudge_hu
      real(kind=REAL64):: time_now
      integer:: istat, ib, ie
      logical:: spn_end_L

      time_now = float(Lctl_step) * Cstv_dt_8 /3600.
      Spn_end_L=.false.
      if (Spn_end_hour>0.) then
          if (time_now>=Spn_end_hour) Spn_end_L=.true.
      end if

      if ((.not. Spn_ON_L) .or. Spn_end_L) return

      call gtmg_start ( 59, 'SPN_main', 1 )

      if (( Lctl_step > 2 ) .and. ( mod(Lctl_step,Spn_interval) == 0 ))  then
      
      
         time_now = float(Lctl_step) * Cstv_dt_8 /3600.

         if (time_now < Spn_relax_periods(1)) then
             Spn_relax_time =Spn_relax_hours
         elseif (time_now >= Spn_relax_periods(1) .and. time_now<Spn_relax_periods(2)) then
             Spn_relax_time = Spn_relax_hours + (time_now-Spn_relax_periods(1))*&
                             (Spn_relax_hours_end-Spn_relax_hours)/&
                             (Spn_relax_periods(2)-Spn_relax_periods(1))
         else
             Spn_relax_time = Spn_relax_hours_end
         end if

         Spn_relax_time = Cstv_dt_8/(Spn_relax_time*3600.)

         call spn_calfiltre (Lctl_step)
         
         if (Spn_weight_L) then
            Spn_weight = sqrt((cos(pi_8*(float(Lctl_step)/float(Spn_ws))))**2)**Spn_wt_pwr
         end if

         if ((Spn_wt_skip_L) .and. Spn_weight < Spn_skip_threshold) then
            if (Lun_out > 0) write(Lun_out,1002) Lctl_step                
            return
         else
            if (Lun_out > 0) write(Lun_out,1001) Lctl_step

            if (.not. Grd_yinyang_L) then

               if (Spn_nudge_UV_L) then
                  call spn_apply (ut1, nest_u, l_minx,l_maxx,l_miny,l_maxy, l_nk)
                  call spn_apply (vt1, nest_v, l_minx,l_maxx,l_miny,l_maxy, l_nk)
               end if
               if (Spn_nudge_TT_L) then
                  call spn_apply (tt1, nest_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)
               end if              
               if (Spn_nudge_HU_L) then
                  ib=(Tr3d_hu-1)*G_nk+1
                  ie=ib+G_nk-1
                  istat=gmm_get('N_hu',nudge_hu)
                  call spn_apply (tracers_P(Tr3d_hu)%pntr, nest_tr(:,:,ib:ie), l_minx,l_maxx,l_miny,l_maxy, l_nk)
               end if

            else
               if (Spn_nudge_UV_L) then
                  istat=gmm_get('N_uu',nudge_u)
                  istat=gmm_get('N_vv',nudge_v)
                  call spn_apply (ut1, nudge_u, l_minx,l_maxx,l_miny,l_maxy, l_nk)
                  call spn_apply (vt1, nudge_v, l_minx,l_maxx,l_miny,l_maxy, l_nk)
               end if
               if (Spn_nudge_TT_L) then
                  istat=gmm_get('N_tt',nudge_t)
                  call spn_apply (tt1, nudge_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)
               end if
               if (Spn_nudge_HU_L) then
                  istat=gmm_get('N_hu',nudge_hu)
                  call spn_apply (tracers_P(Tr3d_hu)%pntr, nudge_hu, l_minx,l_maxx,l_miny,l_maxy, l_nk)
               end if
            end if
         end if
      
      end if
      

      call gtmg_stop ( 59 )
 1001 format(/' SPN_MAIN: Applying spectral nudging at STEP NO ',I10/)
 1002 format(/' SPN_MAIN: Not applying spectral nudging at STEP NO ',I10/)
      
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_main
