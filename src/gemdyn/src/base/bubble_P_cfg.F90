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
!**function bubble_P_cfg - to configure bubble case for Pressure Coord.
!
      integer function bubble_P_cfg()
      use VERgrid_options
      use glb_ld
      use tdpack
      use bubble_options
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer k
      real(kind=REAL64) c1_8,Exner_8,height_8,pres_8,pref_8,ptop_8,htop_8
!
!     ---------------------------------------------------------------
!
      bubble_P_cfg = -1

      ! establish vertical grid configuration
      pref_8 = 1.d5

      G_nk   = bubble_nk
      htop_8 = G_nk*bubble_dz

      if ( hyb(1) < 0 ) then

        !isentropic case
         c1_8=grav_8/(cpd_8*bubble_theta)
         Exner_8=1.d0-c1_8*htop_8
         ptop_8 = Exner_8**(1.d0/cappa_8)*pref_8
!        Uniform distribution of levels in terms of height
         do k=1,G_nk
            height_8=htop_8*(1.d0-(dble(k)-.5d0)/G_nk)
            Exner_8=1.d0-c1_8*height_8
            pres_8=Exner_8**(1.d0/cappa_8)*pref_8
            hyb(k)=(pres_8-ptop_8)/(pref_8-ptop_8)
            hyb(k) = hyb(k) + (1.-hyb(k))*ptop_8/pref_8
         end do

      else

         do k=1024,1,-1
            if(hyb(k) < 0) G_nk=k-1
         end do

      end if

      bubble_P_cfg = 1

      return
      end function bubble_P_cfg

