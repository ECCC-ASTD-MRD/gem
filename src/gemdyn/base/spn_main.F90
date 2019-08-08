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

      subroutine spn_main
      use step_options
      use gem_options
      use spn_options
      use glb_ld
      use cstv
      use lun
      use ldnh
      use sol
      use trp
      use ptopo
      implicit none
#include <arch_specific.hf>

!author
!     Minwei Qian (CCRD) & Bernard Dugas, Syed Husain  (MRB)  - summer 2015
!
!revision
! v4_80 - Qian, Dugas, Hussain            - initial version
! v4_80 - Baek - correction in calls to spn_fld for FFT transpose


      integer offseti, no_steps, tmdt, Nkl
!
!----------------------------------------------------------------------
!
      if (Spn_nudging_S == ' ') return

      tmdt= int(Cstv_dt_8)

      no_steps= Spn_step/tmdt
      no_steps= max(1,no_steps)

      if ( Lctl_step > 2 ) then

         if ( mod(Lctl_step,no_steps) == 0 ) then

            if (Lun_out > 0) write(Lun_out,1001) Lctl_step

            offseti = trp_22n0-1
            Nkl = sol_nk

            if (Lun_debug_L) write(Lun_out,1000)

            if ( index( Spn_nudging_S,'t' ) > 0 )        &
                call spn_fld( ldnh_minx,ldnh_maxx,           &
                ldnh_miny,  ldnh_maxy, ldnh_nj,              &
                trp_12smin, trp_12smax, G_nk,    Nkl, &
                 G_ni, G_nj, trp_22min , trp_22max, trp_22n, &
                 offseti,    Ptopo_npex, Ptopo_npey, 't' )

            if ( index( Spn_nudging_S,'u' ) > 0 )          &
                call spn_fld( ldnh_minx,ldnh_maxx,             &
                  ldnh_miny,  ldnh_maxy, ldnh_nj,              &
                  trp_12smin, trp_12smax, G_nk,    Nkl, &
                  G_ni, G_nj, trp_22min , trp_22max, trp_22n,  &
                  offseti,    Ptopo_npex, Ptopo_npey, 'u' )

            if ( index( Spn_nudging_S,'v' ) > 0 )           &
                 call spn_fld( ldnh_minx,ldnh_maxx,             &
                   ldnh_miny,  ldnh_maxy, ldnh_nj,              &
                   trp_12smin, trp_12smax, G_nk,    Nkl, &
                   G_ni, G_nj, trp_22min , trp_22max, trp_22n,  &
                   offseti,    Ptopo_npex, Ptopo_npey, 'v' )
         end if

      end if

 1000 format( &
           /,'CONTROL OF SPECTRAL NUDGING: (S/R SPN_MAIN)', &
           /,'==========================================='/)
 1001 format(/' In SPN_MAIN, Applying spectral nudging at STEP NO ',I10/)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_main
