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
      implicit none
!
!----------------------------------------------------------------------
!
      if (.not. Spn_ON_L) return

      if (( Lctl_step > 2 ) .and. ( mod(Lctl_step,Spn_interval) == 0 ))  then
      
      
         if (Lun_out > 0) write(Lun_out,1001) Lctl_step
      
         if (Spn_weight_L) then
            Spn_weight = sqrt((cos(pi_8*(float(Lctl_step)/float(Spn_ws))))**2)**Spn_wt_pwr
         end if
         call spn_apply (tt1, nest_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)
         call spn_apply (ut1, nest_u, l_minx,l_maxx,l_miny,l_maxy, l_nk)
         call spn_apply (vt1, nest_v, l_minx,l_maxx,l_miny,l_maxy, l_nk)
      
      end if

 1001 format(/' SPN_MAIN, Applying spectral nudging at STEP NO ',I10/)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_main
