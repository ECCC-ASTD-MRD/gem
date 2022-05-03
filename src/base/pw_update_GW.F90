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

!**s/r pw_update_GW - Update physical quantities WZ and GZ

      subroutine pw_update_GW()
      use dynkernel_options
      use gem_options
      use gem_timing
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use metric
      use tdpack
      implicit none
#include <arch_specific.hf>

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk+1) :: fi
!     ________________________________________________________________

      if (Schm_autobar_L) return

      call gemtime_start ( 5, 'PW_UPDATE', 0)

      pw_wz_plus= wt1
      pw_me_plus= fis0

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then
         pw_gz_plus(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1:l_nk)=&
         GVM%zmom_8(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1:l_nk)*grav_8
      else
         call diag_fi (fi, st1, tt1, qt1, l_minx,l_maxx,l_miny,l_maxy,&
                              G_nk, 1-G_halox*west ,l_ni+G_halox*east,&
                              1-G_haloy*south,l_nj+G_haloy*north)
         pw_gz_plus(:,:,1:l_nk)= fi(:,:,1:l_nk)
      end if

      call gemtime_stop (5)
!     ________________________________________________________________
!
      return
      end subroutine pw_update_GW
